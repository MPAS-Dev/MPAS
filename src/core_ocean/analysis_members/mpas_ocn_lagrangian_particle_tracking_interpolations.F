! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  ocn_lagrangian_particle_tracking_interpolations
!
!> \brief   LIGHT Vector reconstruction and filtering module
!> \author  Phillip J. Wolfram
!> \date    07/21/2015
!> \details
!> This module provides routines for performing vector interpolations
!> and spatial filtering.
!
!-----------------------------------------------------------------------
module ocn_lagrangian_particle_tracking_interpolations

  use mpas_derived_types
  use mpas_constants
  use mpas_rbf_interpolation
  use mpas_geometry_utils
  use mpas_vector_reconstruction
  use mpas_dmpar

  implicit none

  contains

!***********************************************************************
!
!  routine ocn_vector_vertex_to_cell_center
!
!> \brief   Interpolate vertex values to cell center
!> \author  Phillip Wolfram
!> \date    08/01/2014
!> \details
!>  Purpose: interpolate vector field at cell center locations from vertex values
!>           using Wachspress interpolation
!>  Input: vertex vector data and mesh information
!>  Output: interpolated cell values
!-----------------------------------------------------------------------
  subroutine ocn_vector_vertex_to_cell_center(meshPool, & !{{{
      uvReconstructX, uvReconstructY, uvReconstructZ, &
      ucReconstructX, ucReconstructY, ucReconstructZ)

    implicit none

    ! input variables
    type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh information
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: uvReconstructX, & !< Input: Vertex x values
                                                              uvReconstructY, & !< Input: Vertex y values
                                                              uvReconstructZ    !< Input: Vertex z values

    ! output variables
    real (kind=RKIND), dimension(:,:), pointer, intent(out) :: ucReconstructX !< Output: X Component of velocity
                                                                              !<         reconstructed to cells
    real (kind=RKIND), dimension(:,:), pointer, intent(out) :: ucReconstructY !< Output: Y Component of velocity
                                                                              !<         reconstructed to cells
    real (kind=RKIND), dimension(:,:), pointer, intent(out) :: ucReconstructZ !< Output: Z Component of velocity
                                                                              !<         reconstructed to cells

    ! local variables
    integer, pointer :: nCellsSolve, nVertLevels
    integer, dimension(:), pointer :: nEdgesOnCell
    integer :: aVertex, aCell, aLevel, nLocalVertices
    integer, dimension(:,:), pointer :: verticesOnCell
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex
    ! temporary arrays needed in the (to be constructed) init procedure
    ! note that lambda is going to be constant for this and could be cached
    real (kind=RKIND), dimension(:), allocatable :: lambda
    real (kind=RKIND), dimension(3) :: pointInterp
    real (kind=RKIND), dimension(:,:), allocatable :: pointVertex
    real (kind=RKIND) :: xp,yp,zp

    ucReconstructX = 0.0_RKIND
    ucReconstructY = 0.0_RKIND
    ucReconstructZ = 0.0_RKIND

    ! setup pointers
    call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)

    call mpas_pool_get_array(meshPool, 'xCell', xCell)
    call mpas_pool_get_array(meshPool, 'yCell', yCell)
    call mpas_pool_get_array(meshPool, 'zCell', zCell)

    call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
    call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
    call mpas_pool_get_array(meshPool, 'zVertex', zVertex)


    ! loop over all vertices
    do aCell = 1, nCellsSolve
      ! could precompute the list as an optimization to
      ! remove the following lines !{{{
      nLocalVertices = nEdgesOnCell(aCell)
      ! really, condition is any boundaryVertex in column greater than 0
      allocate(lambda(nLocalVertices), pointVertex(3,nLocalVertices))
      ! get vertex location and cell center locations
      do aVertex = 1, nLocalVertices
        pointVertex(1,aVertex) = xVertex(verticesOnCell(aVertex, aCell))
        pointVertex(2,aVertex) = yVertex(verticesOnCell(aVertex, aCell))
        pointVertex(3,aVertex) = zVertex(verticesOnCell(aVertex, aCell))
      end do
      ! vertex point for reconstruction
      pointInterp(1) = xCell(aCell)
      pointInterp(2) = yCell(aCell)
      pointInterp(3) = zCell(aCell)
      ! get interpolation constants (should be cached as an optimization!)
      lambda = mpas_wachspress_coordinates(nLocalVertices, pointVertex , pointInterp, meshPool)
      !}}}

      do aLevel = 1, nVertLevels
          ! perform interpolation
          ucReconstructX(aLevel,aCell) = sum(uvReconstructX(aLevel,verticesOnCell(1:nLocalVertices,aCell)) * lambda)
          ucReconstructY(aLevel,aCell) = sum(uvReconstructY(aLevel,verticesOnCell(1:nLocalVertices,aCell)) * lambda)
          ucReconstructZ(aLevel,aCell) = sum(uvReconstructZ(aLevel,verticesOnCell(1:nLocalVertices,aCell)) * lambda)
      end do

      deallocate(lambda, pointVertex)
    end do

  end subroutine ocn_vector_vertex_to_cell_center !}}}

!***********************************************************************
!
!  routine ocn_simple_vector_shapiro_filter
!
!> \brief   Do 1 pass of simple shapiro filter
!> \author  Phillip Wolfram
!> \date    08/01/2014
!> \details
!>  Purpose: one pass of digital shapiro filter to vertexes, back to cells
!>  Input: cell centered data and mesh information
!>  Output: filtered cell values
!-----------------------------------------------------------------------
    subroutine ocn_simple_vector_shapiro_filter(meshPool, scratchPool, boundaryVertex, boundaryCell, &
        ucReconstructX, ucReconstructY, ucReconstructZ) !{{{
      implicit none

      type (mpas_pool_type), pointer, intent(in) :: meshPool, scratchPool
      type (field2DInteger), pointer, intent(in) :: boundaryVertex, boundaryCell
      type (field2DReal), pointer, intent(inout) :: ucReconstructX, ucReconstructY, ucReconstructZ ! cell center values

      ! local variables
      type (field2DReal), pointer :: uvX , uvY, uvZ ! cell center values

      ! allocate scratch memory
      call mpas_pool_get_field(scratchPool, 'uvX', uvX)
      call mpas_pool_get_field(scratchPool, 'uvY', uvY)
      call mpas_pool_get_field(scratchPool, 'uvZ', uvZ)
      call mpas_allocate_scratch_field(uvX,.True.)
      call mpas_allocate_scratch_field(uvY,.True.)
      call mpas_allocate_scratch_field(uvZ,.True.)

      uvX % array = 0.0_RKIND
      uvY % array = 0.0_RKIND
      uvZ % array = 0.0_RKIND

      ! perform filtering

      ! CC -> vertices
      call ocn_vector_cell_center_to_vertex(meshPool, boundaryVertex % array, boundaryCell % array, &
        ucReconstructX % array, ucReconstructY % array, ucReconstructZ % array, &
        uvX % array, uvY % array, uvZ % array)
      ! do halo exchanges
      call mpas_dmpar_exch_halo_field(uvX)
      call mpas_dmpar_exch_halo_field(uvY)
      call mpas_dmpar_exch_halo_field(uvZ)
      ! vertices -> CC
      call ocn_vector_vertex_to_cell_center(meshPool, &
        uvX % array, uvY % array, uvZ % array, &
        ucReconstructX % array, ucReconstructY % array, ucReconstructZ % array)
      ! do halo exchanges
      call mpas_dmpar_exch_halo_field(ucReconstructX)
      call mpas_dmpar_exch_halo_field(ucReconstructY)
      call mpas_dmpar_exch_halo_field(ucReconstructZ)

      ! N.B., effect of forgetting halo exchange may be subtle for a single pass

      ! deallocate scratch memory
      call mpas_deallocate_scratch_field(uvX , .True.)
      call mpas_deallocate_scratch_field(uvY , .True.)
      call mpas_deallocate_scratch_field(uvZ , .True.)

    end subroutine ocn_simple_vector_shapiro_filter !}}}

end module ocn_lagrangian_particle_tracking_interpolations

