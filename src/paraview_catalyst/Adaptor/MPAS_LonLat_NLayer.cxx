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

void create_lonlat3D_mesh(
		 vtkUnstructuredGrid* grid,
                 int meshType,
                 int numberOfCells,
                 int numberOfVertices,
                 int verticesPerCell,
                 int numberOfGhosts,
                 int* ghostCell, int* ghostLevel,
                 double* xVertex, double* yVertex,
                 int* nEdgesOnCell,
                 int* vertices,
                 vector<bool> const& makeCell,
                 float zFactor);

//////////////////////////////////////////////////////////////////////////
//
// Create the primal and dual cartesian grids
// 2D location with vertex depth making 3D cells
// Input can be locations in x,y or lon,lat
//
//////////////////////////////////////////////////////////////////////////

void create_lonlat3D_grids(vtkCPDataDescription* data,
			   double* xCell,
                           double* yCell,
                           double* xVertex,
                           double* yVertex,
                           float zFactor)
{
  // Create the primal mesh (Voronoi polygons)
  if (data->GetIfGridIsNecessary("LON_LAT_NLAYER-primal") &&
      data->GetInputDescriptionByName("LON_LAT_NLAYER-primal")->GetGrid() == NULL) {
    vtkUnstructuredGrid* pgrid = vtkUnstructuredGrid::New();
    data->GetInputDescriptionByName("LON_LAT_NLAYER-primal")->SetGrid(pgrid);
    pgrid->Initialize();

    // Insert vertices which delineate primal mesh
    vtkPoints* pPts = vtkPoints::New();
    pgrid->SetPoints(pPts);
    pgrid->Allocate(totalPrimalCells, totalPrimalCells);
    for (int i = 0; i < nPrimalVerts; i++)
      for (int j = 0; j < (nVertLevels + 1); j++)
        pPts->InsertNextPoint(xVertex[i], yVertex[i], zFactor * j);

    // Make arrays to indicate regular cells which will be used to create mesh
    vector<bool>& usePrimalCell = usedCells["LON_LAT_NLAYER-primal"];
    usePrimalCell.resize(nPrimalCells);
    for (int i = 0; i < nPrimalCells; i++) {
      usePrimalCell[i] = true;
      for (int j = 0; j < nPrimalVertsPerCell; j++) {
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
    create_lonlat3D_mesh(
		pgrid, PRIMAL, nPrimalCells, nPrimalVerts, nPrimalVertsPerCell,
                nPrimalGhosts, primalGhostCell, primalGhostHalo,
                xVertex, yVertex,
                nEdgesOnCell, verticesOnCell, usePrimalCell, zFactor);

    vtkFloatArray* rankarr = vtkFloatArray::New();
    rankarr->SetName(rankName[PRIMAL].c_str());
    rankarr->SetNumberOfComponents(1);
    rankarr->Allocate(totalPrimalCells);
    pgrid->GetCellData()->AddArray(rankarr);

    vtkFloatArray* maskarr = vtkFloatArray::New();
    maskarr->SetName(maskName[PRIMAL].c_str());
    maskarr->SetNumberOfComponents(1);
    maskarr->Allocate(totalPrimalCells);
    pgrid->GetCellData()->AddArray(maskarr);

    for (int id = 0; id < totalPrimalCells; id++)
      rankarr->InsertNextValue((float) rank);
    rankarr->Delete();

    for (int j = 0; j < nPrimalCells; j++) {
      for (int lev = 0; lev < nVertLevels; lev++) {
        if (usePrimalCell[j] == 1) {
          int findx = j * nVertLevels + lev;
          maskarr->InsertNextValue((float) vertexMask[findx]);
        }
      }
    }

    pPts->Delete();
    pgrid->Delete();
  }

  // Create the dual mesh (Delaunay triangles)
  if (data->GetIfGridIsNecessary("LON_LAT_NLAYER-dual") &&
      data->GetInputDescriptionByName("LON_LAT_NLAYER-dual")->GetGrid() == NULL) {
    vtkUnstructuredGrid* dgrid = vtkUnstructuredGrid::New();
    data->GetInputDescriptionByName("LON_LAT_NLAYER-dual")->SetGrid(dgrid);
    dgrid->Initialize();

    // Insert cell centers which delineate dual mesh
    vtkPoints* dPts = vtkPoints::New();
    dgrid->SetPoints(dPts);
    dgrid->Allocate(totalDualCells, totalDualCells);
    for (int i = 0; i < nDualVerts; i++)
      for (int j = 0; j < (nVertLevels + 1); j++)
        dPts->InsertNextPoint(xCell[i], yCell[i], zFactor * j);

    // Dual mesh has boundary cells in it with indices = nCells + 1
    // Make arrays to indicate regular cells which will be used to create mesh
    vector<bool>& useDualCell = usedCells["LON_LAT_NLAYER-dual"];
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
    create_lonlat3D_mesh(
		dgrid, DUAL, nDualCells, nDualVerts, nDualVertsPerCell,
                nDualGhosts, dualGhostCell, dualGhostHalo,
                xCell, yCell,
                nEdgesOnCell, cellsOnVertex, useDualCell, zFactor);

    // Dual mesh has boundary cells which were not created so get the
    // actual number of cells which is different from nDualCells
    int nActualDualCells = dgrid->GetNumberOfCells();

    vtkFloatArray* drankarr = vtkFloatArray::New();
    drankarr->SetName(rankName[DUAL].c_str());
    drankarr->SetNumberOfComponents(1);
    drankarr->Allocate(nActualDualCells);
    dgrid->GetCellData()->AddArray(drankarr);

    vtkFloatArray* dmaskarr = vtkFloatArray::New();
    dmaskarr->SetName(maskName[DUAL].c_str());
    dmaskarr->SetNumberOfComponents(1);
    dmaskarr->Allocate(totalDualCells);
    dgrid->GetCellData()->AddArray(dmaskarr);

    for (int id = 0; id < nActualDualCells; id++)
      drankarr->InsertNextValue((float) rank);
    drankarr->Delete();

    for (int j = 0; j < nDualCells; j++) {
      for (int lev = 0; lev < nVertLevels; lev++) {
        if (useDualCell[j] == 1) {
          int findx = j * nVertLevels + lev;
          dmaskarr->InsertNextValue((float) cellMask[findx]);
        }
      }
    }

    dPts->Delete();
    dgrid->Delete();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Create the primal or dual cartesian mesh
//
//////////////////////////////////////////////////////////////////////////

void create_lonlat3D_mesh(
		 vtkUnstructuredGrid* grid,
                 int meshType,
                 int numberOfCells,
                 int numberOfVertices,
                 int verticesPerCell,
                 int numberOfGhosts,
                 int* ghostCell, int* ghostLevel,
                 double* xVertex, double* yVertex,
                 int* nEdgesOnCell,
                 int* vertices,
                 vector<bool> const& makeCell,
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
  int cellType = VTK_WEDGE;
  int nEdges = verticesPerCell;

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
        nEdges = nEdgesOnCell[id];
        if (nEdges == 3)
          cellType = VTK_WEDGE;
        else if (nEdges == 4)
          cellType = VTK_HEXAHEDRON;
        else if (nEdges == 5)
          cellType = VTK_PENTAGONAL_PRISM;
        else {
          nEdges = 6;
          cellType = VTK_HEXAGONAL_PRISM;
        }
      }

      // First vertex will be correct and others will wrap if necessary
      double xRef = xVertex[vertices[indx] - 1];
      
      for (int vert = 0; vert < nEdges; vert++) {
        // Fortran indexing starts with 1
        int cIndx = vertices[indx] - 1;

        // If all vertices are not in bounds create wraparound point
        double xVert = xVertex[cIndx];
        double yVert = yVertex[cIndx];
        bool makeWrapPt = false;

        if ((xRef - xVert) < -M_PI) {
          xVert -= 2.0 * M_PI;
          makeWrapPt = true;
        }
        else if ((xRef - xVert) > M_PI) {
          xVert += 2.0 * M_PI;
          makeWrapPt = true;
        }

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
      for (int g = 0; g < numberOfGhosts; g++) {
        int ghostIndx = ghostCell[g] - 1;
        if (id == ghostIndx) {
          int thisIndx = (id - numberSkipped) * nVertLevels;
          for (int lev = 0; lev < nVertLevels; lev++)
            ghosts->SetValue(thisIndx + lev, ghostLevel[g]);
        }
      }
    }
  }
}
