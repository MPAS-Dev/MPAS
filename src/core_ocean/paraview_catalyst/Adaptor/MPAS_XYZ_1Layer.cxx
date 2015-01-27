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

void create_xyz2D_mesh(
		 vtkUnstructuredGrid* grid,
                 int meshType,
                 int numberOfCells,
                 int numberOfVertices,
                 int verticesPerCell,
                 int numberOfGhosts,
                 int* ghostCell, int* ghostLevel,
                 double* xVertex, double* yVertex, double* zVertex,
                 double* xCenter, double* yCenter, double* zCenter,
                 int* nEdgesOnCell,
                 int* vertices,
                 vector<bool> const& makeCell);

//////////////////////////////////////////////////////////////////////////
//
// Create the primal and dual spherical grids
// 2D cells with vertex depth as components
//
//////////////////////////////////////////////////////////////////////////

void create_xyz2D_grids(
			vtkCPDataDescription* data,
                        double* xCell, double* yCell, double* zCell,
                        double* xVertex, double* yVertex, double* zVertex)
{
  // Create the primal mesh (Voronoi polygons)
  if (data->GetIfGridIsNecessary("X_Y_Z_1LAYER-primal") &&
      data->GetInputDescriptionByName("X_Y_Z_1LAYER-primal")->GetGrid() == NULL) {
    vtkUnstructuredGrid* pgrid = vtkUnstructuredGrid::New();
    data->GetInputDescriptionByName("X_Y_Z_1LAYER-primal")->SetGrid(pgrid);
    pgrid->Initialize();

    // Insert vertices which delineate primal mesh
    vtkPoints* pPts = vtkPoints::New();
    pgrid->SetPoints(pPts);
    pgrid->Allocate(nPrimalCells, nPrimalCells);
    for (int i = 0; i < nPrimalVerts; i++)
      pPts->InsertNextPoint(xVertex[i], yVertex[i], zVertex[i]);

    // Allocate ghost level array to be attached to cell data
    vtkUnsignedCharArray* pghost = vtkUnsignedCharArray::New();
    pghost->SetName("vtkGhostLevels");
    pghost->Allocate(nPrimalCells);
    pgrid->GetCellData()->AddArray(pghost);

    // Make arrays to indicate good cells which will be used to create mesh
    vector<bool>& usePrimalCell = usedCells["X_Y_Z_1LAYER-primal"];
    usePrimalCell.resize(nPrimalCells);
    for (int i = 0; i < nPrimalCells; i++) {
      usePrimalCell[i] = true;
      for (int j = 0; j < nEdgesOnCell[i]; j++) {
        if (verticesOnCell[(i * nPrimalVertsPerCell) + j] >= (nPrimalVerts + 1)) {
          usePrimalCell[i] = false;
        }
      }
    }

    // Create the primal mesh of Voronoi polygons
    create_xyz2D_mesh(
		pgrid, PRIMAL, nPrimalCells, nPrimalVerts, nPrimalVertsPerCell,
                nPrimalGhosts, primalGhostCell, primalGhostHalo,
                xVertex, yVertex, zVertex,
                xCell, yCell, zCell,
                nEdgesOnCell, verticesOnCell, usePrimalCell);

    vtkFloatArray* rankarr = vtkFloatArray::New();
    rankarr->SetName(rankName[PRIMAL].c_str());
    rankarr->SetNumberOfComponents(1);
    rankarr->Allocate(nPrimalCells);
    pgrid->GetCellData()->AddArray(rankarr);

    for (int id = 0; id < nPrimalCells; id++)
      rankarr->InsertNextValue((float) rank);
    rankarr->Delete();

    pPts->Delete();
    pgrid->Delete();
  }

  // Create the dual mesh (Delaunay triangles)
  if (data->GetIfGridIsNecessary("X_Y_Z_1LAYER-dual") &&
      data->GetInputDescriptionByName("X_Y_Z_1LAYER-dual")->GetGrid() == NULL ) {
    vtkUnstructuredGrid* dgrid = vtkUnstructuredGrid::New();
    data->GetInputDescriptionByName("X_Y_Z_1LAYER-dual")->SetGrid(dgrid);
    dgrid->Initialize();

    // Insert cell centers which delineate dual mesh
    vtkPoints* dPts = vtkPoints::New();
    dgrid->SetPoints(dPts);
    dgrid->Allocate(nDualCells, nDualCells);
    for (int i = 0; i < nDualVerts; i++)
      dPts->InsertNextPoint(xCell[i], yCell[i], zCell[i]);

    // Allocate ghost level array to be attached to cell data
    vtkUnsignedCharArray* dghost = vtkUnsignedCharArray::New();
    dghost->SetName("vtkGhostLevels");
    dghost->Allocate(nDualCells);
    dgrid->GetCellData()->AddArray(dghost);

    // Dual mesh has boundary cells in it with indices = nCells + 1
    // Make arrays to indicate good cells which will be used to create mesh
    vector<bool>& useDualCell = usedCells["X_Y_Z_1LAYER-dual"];
    useDualCell.resize(nDualCells);
    for (int i = 0; i < nDualCells; i++) {
      useDualCell[i] = true;
      for (int j = 0; j < nDualVertsPerCell; j++) {
        if (cellsOnVertex[(i * nDualVertsPerCell) + j] >= (nDualVerts + 1)) {
          useDualCell[i] = false;
        }
      }
    }

    // Create the dual mesh of Delaunay triangles
    create_xyz2D_mesh(
		dgrid, DUAL, nDualCells, nDualVerts, nDualVertsPerCell,
                nDualGhosts, dualGhostCell, dualGhostHalo,
                xCell, yCell, zCell,
                xVertex, yVertex, zVertex,
                nEdgesOnCell, cellsOnVertex, useDualCell);

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

    dPts->Delete();
    dgrid->Delete();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Create the primal or dual spherical mesh
//
//////////////////////////////////////////////////////////////////////////

void create_xyz2D_mesh(
		 vtkUnstructuredGrid* grid,
                 int meshType,
                 int numberOfCells,
                 int numberOfVertices,
                 int verticesPerCell,
                 int numberOfGhosts,
                 int* ghostCell, int* ghostLevel,
                 double* xVertex, double* yVertex, double* zVertex,
                 double* xCenter, double* yCenter, double* zCenter,
                 int* nEdgesOnCell,
                 int* vertices,
                 vector<bool> const& makeCell)
{
  std::string suffix = (meshType == PRIMAL) ? "-primal" : "-dual";
  vtkPoints* pts = grid->GetPoints();

  // Allocate ghost level array to be attached to cell data
  vtkUnsignedCharArray* ghosts =  vtkUnsignedCharArray::SafeDownCast(
        grid->GetCellData()->GetArray("vtkGhostLevels"));

  vtkIdType cell[verticesPerCell];
  int nEdges = verticesPerCell;
  int cellType = VTK_TRIANGLE;
  int numberSkipped = 0;

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
        cellType = VTK_POLYGON;
      }

      for (int vert = 0; vert < nEdges; vert++) {
          // Fortran indexing starts with 1
          int cIndx = vertices[indx + vert] - 1;
          cell[vert] = cIndx;
       }
      vtkIdType newCellId = grid->InsertNextCell(cellType, nEdges, cell);
      if (meshType == DUAL) {
        orient_sphere_triangle_cell(grid, newCellId, cell);
      }
       ghosts->InsertNextValue(0);
 
       // Check list of ghost cell indices to see if this id is in it
       // Boundary cells must be skipped so alter the id index
       for (int g = 0; g < numberOfGhosts; g++) {
         int ghostIndx = ghostCell[g] - 1;
         if (id == ghostIndx) {
           int thisIndx = id - numberSkipped;
           ghosts->SetValue(thisIndx, ghostLevel[g]);
         }
       }
    }
  }
}
