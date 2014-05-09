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
#include "vtkMath.h"
#include "vtkMPI.h"
#include "vtkMPICommunicator.h"
#include "vtkMPIController.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"

#include <float.h>
#include <sstream>
#include <string>

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

using namespace std;
using namespace MPAS;

void create_xyz3D_mesh(
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
// To get vertical levels starting from spherical x,y,z
// Convert to spherical coordinates, increment rho, convert back to cartesian
//
//////////////////////////////////////////////////////////////////////////

int CartesianToSpherical(double x, double y, double z,
                         double* rho, double* phi, double* theta)
{
  double trho, ttheta, tphi;

  trho = sqrt((x*x) + (y*y) + (z*z));
  ttheta = atan2(y, x);
  tphi = acos(z/(trho));
  if (vtkMath::IsNan(trho) || vtkMath::IsNan(ttheta) || vtkMath::IsNan(tphi)) {
    return -1;
  }
  *rho = trho;
  *theta = ttheta;
  *phi = tphi;
  return 0;
}

int SphericalToCartesian(double rho, double phi, double theta,
                         double* x, double* y, double* z)
{
  double tx, ty, tz;

  tx = rho* sin(phi) * cos(theta);
  ty = rho* sin(phi) * sin(theta);
  tz = rho* cos(phi);
  if (vtkMath::IsNan(tx) || vtkMath::IsNan(ty) || vtkMath::IsNan(tz)) {
    return -1;
  }

  *x = tx;
  *y = ty;
  *z = tz;
  return 0;
}

//////////////////////////////////////////////////////////////////////////
//
// Create the primal and dual spherical grids
// 2D cells with vertex depth as components
//
//////////////////////////////////////////////////////////////////////////

void create_xyz3D_grids(
			vtkCPDataDescription* data,
                        double* xCell,double* yCell, double* zCell,
                        double* xVertex, double* yVertex, double* zVertex,
                        float zFactor)
{
  // Create the primal mesh (Voronoi polygons)
  if (data->GetIfGridIsNecessary("X_Y_Z_NLAYER-primal") &&
      data->GetInputDescriptionByName("X_Y_Z_NLAYER-primal")->GetGrid() == NULL) {
    vtkUnstructuredGrid* pgrid = vtkUnstructuredGrid::New();
    data->GetInputDescriptionByName("X_Y_Z_NLAYER-primal")->SetGrid(pgrid);
    pgrid->Initialize();

    // Insert vertices which delineate primal mesh
    vtkPoints* pPts = vtkPoints::New();
    pgrid->SetPoints(pPts);
    pgrid->Allocate(totalPrimalCells, totalPrimalCells);

    double rho, phi, theta, rholevel;
    double x, y, z;
    for (int i = 0; i < nPrimalVerts; i++) {
      CartesianToSpherical(xVertex[i], yVertex[i], zVertex[i],
                           &rho, &phi, &theta);
      for (int j = 0; j < (nVertLevels + 1); j++) {
        rholevel = rho + (zFactor * j);
        SphericalToCartesian(rholevel, phi, theta, &x, &y, &z);
        pPts->InsertNextPoint(x, y, z);
      }
    }

    // Allocate ghost level array to be attached to cell data
    vtkUnsignedCharArray* pghost = vtkUnsignedCharArray::New();
    pghost->SetName("vtkGhostLevels");
    pghost->Allocate(totalPrimalCells);
    pgrid->GetCellData()->AddArray(pghost);

    // Make arrays to indicate good cells which will be used to create mesh
    vector<bool>& usePrimalCell = usedCells["X_Y_Z_NLAYER-primal"];
    usePrimalCell.resize(nPrimalCells);
    for (int i = 0; i < nPrimalCells; i++) {
      usePrimalCell[i] = true;
      for (int j = 0; j < nPrimalVertsPerCell; j++) {
        if (verticesOnCell[(i * nPrimalVertsPerCell) + j] >= (nPrimalVerts + 1)) {
          usePrimalCell[i] = false;
        }
      }
    }

    // Create the primal mesh of Voronoi polygons
    create_xyz3D_mesh(
                pgrid, PRIMAL, nPrimalCells, nPrimalVerts, nPrimalVertsPerCell,
                nPrimalGhosts, primalGhostCell, primalGhostHalo,
                xVertex, yVertex, zVertex,
                xCell, yCell, zCell,
                nEdgesOnCell, verticesOnCell, usePrimalCell);

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
  if (data->GetIfGridIsNecessary("X_Y_Z_NLAYER-dual") &&
      data->GetInputDescriptionByName("X_Y_Z_NLAYER-dual")->GetGrid() == NULL) {
    vtkUnstructuredGrid* dgrid = vtkUnstructuredGrid::New();
    data->GetInputDescriptionByName("X_Y_Z_NLAYER-dual")->SetGrid(dgrid);
    dgrid->Initialize();

    // Insert cell centers which delineate dual mesh
    vtkPoints* dPts = vtkPoints::New();
    dgrid->SetPoints(dPts);
    dgrid->Allocate(nDualCells, nDualCells);

    double rho, phi, theta, rholevel;
    double x, y, z;
    for (int i = 0; i < nDualVerts; i++) {
      CartesianToSpherical(xCell[i], yCell[i], zCell[i],
                           &rho, &phi, &theta);
      for (int j = 0; j < (nVertLevels + 1); j++) {
        rholevel = rho + (zFactor * j);
        SphericalToCartesian(rholevel, phi, theta, &x, &y, &z);
        dPts->InsertNextPoint(x, y, z);
      }
    }

    // Allocate ghost level array to be attached to cell data
    vtkUnsignedCharArray* dghost = vtkUnsignedCharArray::New();
    dghost->SetName("vtkGhostLevels");
    dghost->Allocate(totalDualCells);
    dgrid->GetCellData()->AddArray(dghost);

    // Dual mesh has boundary cells in it with indices = nCells + 1
    // Make arrays to indicate good cells which will be used to create mesh
    vector<bool>& useDualCell = usedCells["X_Y_Z_NLAYER-dual"];
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
    create_xyz3D_mesh(
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
// Create the primal or dual spherical mesh
//
//////////////////////////////////////////////////////////////////////////

void create_xyz3D_mesh(
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
  vtkIdType cell3d[2*verticesPerCell];
  int numberSkipped = 0;

  // Set type and number of edges for DUAL which is constant
  int nEdges = verticesPerCell;
  int cellType = VTK_WEDGE;

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

      for (int vert = 0; vert < nEdges; vert++) {
        // Fortran indexing starts with 1
        int cIndx = vertices[indx + vert] - 1;
        cell[vert] = cIndx * (nVertLevels + 1);
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
      // Boundary cells must be skipped so alter the id index
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
