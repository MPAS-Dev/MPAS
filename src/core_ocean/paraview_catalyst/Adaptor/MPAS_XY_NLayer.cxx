#include "MPASAdaptor.h"
#include "GridUtils.h"

#include "MPASAdaptorAPIMangling.h"

#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkMPI.h"
#include "vtkMPICommunicator.h"
#include "vtkMPIController.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"

#include "vtkTrivialProducer.h"

#include <float.h>
#include <sstream>
#include <string>

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

using namespace std;
using namespace MPAS;

void create_xy3D_mesh(
		 vtkUnstructuredGrid* grid,
                 int meshType,
                 int numberOfCells,
                 int numberOfVertices,
                 int verticesPerCell,
                 int numberOfGhosts,
                 int* ghostCell, int* ghostLevel,
                 double* xVertex, double* yVertex,
                 double* xCenter, double* yCenter,
                 int* nEdgesOnCell,
                 int* vertices,
                 vector<bool> const& makeCell,
                 double* offset,
                 double* poffset,
                 float zFactor);

//////////////////////////////////////////////////////////////////////////
//
// Create the primal and dual cartesian grids
// 2D location with vertex depth making 3D cells
// Input can be locations in x,y or lon,lat
//
//////////////////////////////////////////////////////////////////////////

void create_xy3D_grids(vtkCPDataDescription* data,
		       double* xCell,
                       double* yCell,
                       double* xVertex,
                       double* yVertex,
                       float zFactor)
{
  // Get bounding box for a primal mesh cell used to find vertices out of range
  double pminRange[2], pmaxRange[2], poffset[2];
  pminRange[0] = pmaxRange[0] = xCell[0];
  pminRange[1] = pmaxRange[1] = yCell[0];

  // Assuming first cell is good which may be a bad assumption
  for (int i = 0; i < nEdgesOnCell[0]; i++) {
    float x = xVertex[verticesOnCell[i] - 1];
    float y = yVertex[verticesOnCell[i] - 1];
    if (x < pminRange[0]) pminRange[0] = x;
    if (x > pmaxRange[0]) pmaxRange[0] = x;
    if (y < pminRange[1]) pminRange[1] = y;
    if (y > pmaxRange[1]) pmaxRange[1] = y;
  }
  for (int dim = 0; dim < 2; dim++)
    poffset[dim] = pmaxRange[dim] - pminRange[dim];

  // Insert vertices which delineate primal mesh
  vtkPoints* pPts = vtkPoints::New();
  for (int i = 0; i < nPrimalVerts; i++)
    for (int j = 0; j < (nVertLevels+1); j++)
      pPts->InsertNextPoint(xVertex[i], yVertex[i], zFactor * j);

  // Insert cell centers which delineate dual mesh
  vtkPoints* dPts = vtkPoints::New();
  for (int i = 0; i < nDualVerts; i++)
    for (int j = 0; j < (nVertLevels+1); j++)
      dPts->InsertNextPoint(xCell[i], yCell[i], zFactor * j);

  // Get the overall problem range
  double vertBound[6], cellBound[6];
  double minRange[3], maxRange[3];
  pPts->GetBounds(vertBound);
  dPts->GetBounds(cellBound);
  for (int i = 0; i < 3; i++) {
    minRange[i] = min(cellBound[i*2], vertBound[i*2]);
    maxRange[i] = max(cellBound[i*2 + 1], vertBound[i*2 + 1]);
  }

  // Problem bounding box is the offset for wraparound
  double offset[3];
  for (int dim = 0; dim < 3; dim++)
    offset[dim] = maxRange[dim] - minRange[dim];

  // Create the primal mesh (Voronoi polygons)
  if (data->GetIfGridIsNecessary("X_Y_NLAYER-primal") && 
      data->GetInputDescriptionByName("X_Y_NLAYER-primal")->GetGrid() == NULL) {
    std::cout << "making a X_Y_NLAYER-primal grid\n";
    vtkUnstructuredGrid* pgrid = vtkUnstructuredGrid::New();
    data->GetInputDescriptionByName("X_Y_NLAYER-primal")->SetGrid(pgrid);
    pgrid->Initialize();

    pgrid->SetPoints(pPts);
    pgrid->Allocate(totalPrimalCells, totalPrimalCells);

    // Make arrays to indicate regular cells which will be used to create mesh
    vector<bool>& usePrimalCell = usedCells["X_Y_NLAYER-primal"];
    usePrimalCell.resize(nPrimalCells);
    for (int i = 0; i < nPrimalCells; i++) {
      usePrimalCell[i] = true;
      for (int j = 0; j < nEdgesOnCell[i]; j++) {
        if (verticesOnCell[(i * nPrimalVertsPerCell) + j] >= (nPrimalVerts + 1)) {
          usePrimalCell[i] = false;
        }
      }
    }

    // Allocate ghost level array to be attached to cell data
    vtkUnsignedCharArray* pghost = vtkUnsignedCharArray::New();
    pghost->SetName("vtkGhostLevels");
    pghost->Allocate(totalPrimalCells);
    pgrid->GetCellData()->AddArray(pghost);

    // Create the primal mesh of Voronoi polygons
    create_xy3D_mesh(
		pgrid, PRIMAL, nPrimalCells, nPrimalVerts, nPrimalVertsPerCell,
                nPrimalGhosts, primalGhostCell, primalGhostHalo,
                xVertex, yVertex,
                xCell, yCell,
                nEdgesOnCell, verticesOnCell, usePrimalCell,
                offset, poffset, zFactor);

    int nActualPrimalCells = pgrid->GetNumberOfCells();

    vtkFloatArray* rankarr = vtkFloatArray::New();
    rankarr->SetName(rankName[PRIMAL].c_str());
    rankarr->SetNumberOfComponents(1);
    rankarr->Allocate(totalPrimalCells);
    pgrid->GetCellData()->AddArray(rankarr);

    for (int id = 0; id < totalPrimalCells; id++)
      rankarr->InsertNextValue((float) rank);
    rankarr->Delete();

    pgrid->Delete();
  }

  // Create the dual mesh (Delaunay triangles)
  if (data->GetIfGridIsNecessary("X_Y_NLAYER-dual") && 
      data->GetInputDescriptionByName("X_Y_NLAYER-dual")->GetGrid() == NULL) {
    std::cout << "making a X_Y_NLAYER-dual grid\n";
    vtkUnstructuredGrid* dgrid = vtkUnstructuredGrid::New();
    data->GetInputDescriptionByName("X_Y_NLAYER-dual")->SetGrid(dgrid);
    dgrid->Initialize();

    dgrid->SetPoints(dPts);
    dgrid->Allocate(totalDualCells, totalDualCells);

    // Dual mesh has boundary cells in it with indices = nCells + 1
    // Make arrays to indicate regular cells which will be used to create mesh
    vector<bool>& useDualCell = usedCells["X_Y_NLAYER-dual"];
    useDualCell.resize(nDualCells);
    for (int i = 0; i < nDualCells; i++) {
      useDualCell[i] = true;
      for (int j = 0; j < nDualVertsPerCell; j++) {
        if (cellsOnVertex[(i * nDualVertsPerCell) + j] >= (nDualVerts + 1)) {
          useDualCell[i] = false;
        }
      }
    }

    // Allocate ghost level array to be attached to cell data
    vtkUnsignedCharArray* dghost = vtkUnsignedCharArray::New();
    dghost->SetName("vtkGhostLevels");
    dghost->Allocate(totalDualCells);
    dgrid->GetCellData()->AddArray(dghost);

    // Create the dual mesh of Delaunay triangles
    create_xy3D_mesh(
		dgrid, DUAL, nDualCells, nDualVerts, nDualVertsPerCell,
                nDualGhosts, dualGhostCell, dualGhostHalo,
                xCell, yCell,
                xVertex, yVertex,
                nEdgesOnCell, cellsOnVertex, useDualCell,
                offset, poffset, zFactor);

    // Dual mesh has boundary cells which were not created so get the
    // actual number of cells which is different from nDualCells
    int nActualDualCells = dgrid->GetNumberOfCells();

    vtkFloatArray* drankarr = vtkFloatArray::New();
    drankarr->SetName(rankName[DUAL].c_str());
    drankarr->SetNumberOfComponents(1);
    drankarr->Allocate(nActualDualCells);
    dgrid->GetCellData()->AddArray(drankarr);

    for (int id = 0; id < nActualDualCells; id++)
      drankarr->InsertNextValue((float) rank);
    drankarr->Delete();

    dgrid->Delete();
  }

  pPts->Delete();
  dPts->Delete();
}

//////////////////////////////////////////////////////////////////////////
//
// Create the primal or dual cartesian mesh
//
//////////////////////////////////////////////////////////////////////////

void create_xy3D_mesh(
         	 vtkUnstructuredGrid* grid,
                 int meshType,
                 int numberOfCells,
                 int numberOfVertices,
                 int verticesPerCell,
                 int numberOfGhosts,
                 int* ghostCell, int* ghostLevel,
                 double* xVertex, double* yVertex,
                 double* xCenter, double* yCenter,
                 int* nEdgesOnCell,
                 int* vertices,
                 vector<bool> const& makeCell,
                 double* offset,
                 double* poffset,
                 float zFactor)
{
  std::string suffix = (meshType == PRIMAL) ? "-primal" : "-dual";
  vtkPoints* pts = grid->GetPoints();

  // Allocate ghost level array to be attached to cell data
  vtkUnsignedCharArray* ghosts =  vtkUnsignedCharArray::SafeDownCast(
        grid->GetCellData()->GetArray("vtkGhostLevels"));

  // Add cells for mesh, doing wraparound of vertices when needed
  vtkIdType cell[verticesPerCell];
  vtkIdType cell3d[2*verticesPerCell];
  int numberSkipped = 0;

  // Set type and number of edges for DUAL which is constant
  int cellType;
  int nEdges;
  if (meshType == DUAL) {
    cellType = VTK_WEDGE;
    nEdges = verticesPerCell;
  }

  for (int id = 0; id < numberOfCells; id++) {
    // Set index into vertices which is 2D array
    int indx = id * verticesPerCell;

    // Cell does not have distinct vertex values
    if (makeCell[id] == false) {
      numberSkipped++;
    }
    else {

      // Primal mesh cells are variable size
      if (meshType == PRIMAL) {
        cellType = -1;
        nEdges = nEdgesOnCell[id];
        if (nEdges == 6)
          cellType = VTK_HEXAGONAL_PRISM;
        else if (nEdges == 5)
          cellType = VTK_PENTAGONAL_PRISM;
        else if (nEdges == 4)
          cellType = VTK_HEXAHEDRON;
        else if (nEdges == 3)
          cellType = VTK_WEDGE;
        else
          cout << "Cell on " << meshType << " type " << nEdges << endl;
      }

      if (cellType != -1) {
       for (int vert = 0; vert < nEdges; vert++) {

        // Fortran indexing starts with 1
        int cIndx = vertices[indx] - 1;

        // Verify that point referred to is within bounds of cell
        double xVert = xVertex[cIndx];
        double yVert = yVertex[cIndx];
        bool makeWrapPt = false;
        if ((xVert - poffset[0]) > xCenter[id]) {
          xVert -= offset[0];
          makeWrapPt = true;
        }
        else if ((xVert + poffset[0]) < xCenter[id]) {
          xVert += offset[0];
          makeWrapPt = true;
        }
        if ((yVert - poffset[1]) > yCenter[id]) {
          yVert -= offset[1];
          makeWrapPt = true;
        }
        else if ((yVert + poffset[1]) < yCenter[id]) {
          yVert += offset[1];
          makeWrapPt = true;
        }

        // Vertex is out of the bounding box so add wraparound points all levels
        if (makeWrapPt == true) {
          cIndx = pts->InsertNextPoint(xVert, yVert, 0.0);
          cell[vert] = cIndx;
          for (int j = 1; j < (nVertLevels+1); j++)
            pts->InsertNextPoint(xVert, yVert, zFactor * j);
        } else {
          cell[vert] = cIndx * (nVertLevels + 1);
        }
        indx++;
       }

       // At this point I have the top surface of the cell
       // which must be expressed as one 3D cell per vertex level
       // So this is the model for creating the rest of the column
       for (int level = 0; level < nVertLevels; level++) {
         int indx3d = 0;
         // Top plane
         for (int j = 0; j < nEdges; j++) {
           cell3d[indx3d++] = cell[j] + level;
         }
         // Bottom plane
         for (int j = 0; j < nEdges; j++) {
           cell3d[indx3d++] = cell[j] + level + 1;
         }
         vtkIdType newCellId = grid->InsertNextCell(cellType, 2*nEdges, cell3d);
         if (meshType == DUAL) {
           orient_wedge_cell(grid, newCellId, cell3d);
         }
         ghosts->InsertNextValue(0);
       }
 
       // Check list of ghost cell indices to see if this id is in it
       // Bad cells must be skipped so alter the id index
       for (int g = 0; g < numberOfGhosts; g++) {
         int ghostIndx = ghostCell[g] - 1;
         if (id == ghostIndx) {
           int thisIndx = (id - numberSkipped) * nVertLevels;
           for (int level = 0; level < nVertLevels; level++)
             ghosts->SetValue(thisIndx + level, ghostLevel[g]);
         }
       }
      }
     }
   }
}
