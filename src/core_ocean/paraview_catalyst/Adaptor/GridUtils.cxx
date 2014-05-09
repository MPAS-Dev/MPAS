#include "GridUtils.h"

#include "vtkMath.h"
#include "vtkPoints.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkWedge.h"

void orient_sphere_triangle_cell(
                 vtkUnstructuredGrid *grid,
                 vtkIdType id,
                 vtkIdType cellPts[3])
{
  vtkCell *newCell = grid->GetCell(id);

  // XXX: Assuming the center of the sphere is the origin.

  double norm[3];
  vtkPoints *pts = newCell->GetPoints();
  vtkIdType ids[] = {0, 1, 2};
  vtkTriangle::ComputeNormal(pts, 3, ids, norm);
  double pt[3];
  pts->GetPoint(0, pt);
  vtkMath::Normalize(pt);
  if (abs(acos(vtkMath::Dot(pt, norm))) > M_PI_2) {
    vtkIdType newCell2d[3];
    newCell2d[0] = cellPts[1];
    newCell2d[1] = cellPts[0];
    newCell2d[2] = cellPts[2];
    grid->ReplaceCell(id, 3, newCell2d);
  }
}

void orient_plane_triangle_cell(
                 vtkUnstructuredGrid *grid,
                 vtkIdType id,
                 vtkIdType cellPts[3])
{
  vtkCell *newCell = grid->GetCell(id);

  double norm[3];
  vtkPoints *pts = newCell->GetPoints();
  vtkIdType ids[] = {0, 1, 2};
  vtkTriangle::ComputeNormal(pts, 3, ids, norm);
  if (norm[2] < 0) {
    vtkIdType newCell2d[3];
    newCell2d[0] = cellPts[1];
    newCell2d[1] = cellPts[0];
    newCell2d[2] = cellPts[2];
    grid->ReplaceCell(id, 3, newCell2d);
  }
}

void orient_wedge_cell(
                 vtkUnstructuredGrid *grid,
                 vtkIdType id,
                 vtkIdType cellPts[6])
{
  vtkCell *newCell = grid->GetCell(id);
  vtkWedge *wedge = vtkWedge::SafeDownCast(newCell);
  double center[3];
  double *jacobian[3];
  double jacobian0[3];
  double jacobian1[3];
  double jacobian2[3];
  jacobian[0] = jacobian0;
  jacobian[1] = jacobian1;
  jacobian[2] = jacobian2;
  double derivs[18];
  wedge->GetParametricCenter(center);
  wedge->JacobianInverse(center, jacobian, derivs);
  double determinant = vtkMath::Determinant3x3(jacobian0, jacobian1, jacobian2);
  if (determinant < 0) {
    vtkIdType newCell3d[6];
    newCell3d[0] = cellPts[1];
    newCell3d[1] = cellPts[0];
    newCell3d[2] = cellPts[2];
    newCell3d[3] = cellPts[4];
    newCell3d[4] = cellPts[3];
    newCell3d[5] = cellPts[5];
    grid->ReplaceCell(id, 6, newCell3d);
  }
}
