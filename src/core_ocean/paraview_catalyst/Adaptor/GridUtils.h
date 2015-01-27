#ifndef MPASGridUtils_h
#define MPASGridUtils_h

#include "vtkType.h"

class vtkUnstructuredGrid;

void orient_sphere_triangle_cell(
                 vtkUnstructuredGrid *grid,
                 vtkIdType id,
                 vtkIdType cellPts[3]);

void orient_plane_triangle_cell(
                 vtkUnstructuredGrid *grid,
                 vtkIdType id,
                 vtkIdType cellPts[3]);

void orient_wedge_cell(
                 vtkUnstructuredGrid *grid,
                 vtkIdType id,
                 vtkIdType cellPts[6]);

#endif
