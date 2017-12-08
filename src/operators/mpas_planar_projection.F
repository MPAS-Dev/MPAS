! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_mpas_planar_projection
!
!> \brief  incremental remapping transport scheme for cice
!> \author William Lipscomb
!> \date   March 2015
!> \details
!>  This module contains routines for transport of sea-ice mass and tracers
!>  using incremental remapping
!
!-----------------------------------------------------------------------

module mpas_planar_projection

   use mpas_kind_types
   use mpas_dmpar
   use mpas_io_units

   implicit none
   private
   save 
   
   public :: define_local_to_global_transformations, &
             unit_vector_3d,                         &
             cross_product_3d

   real(kind=RKIND), parameter :: &
        eps11 = 1.0e-11_RKIND           ! small number (= puny in CICE)

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine define_local_to_global_transformations
!
!> \brief  compute matrices for transforming between local and global coordinates
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  This routine computes 3x3 matrices for transforming between global xyz coordinates
!>  and local coordinates based at a cell center, vertex or edge.
!
!-----------------------------------------------------------------------

  subroutine define_local_to_global_transformations(&
       xGlobal, yGlobal, zGlobal,  &
       nElements,                  &
       transLocalToGlobal,         &
       transGlobalToLocal)

    real(kind=RKIND), dimension(nElements), intent(in) ::  &
         xGlobal, yGlobal, zGlobal  !< Input: global coordinates of points on the surface of the sphere

    integer, intent(in) :: &
         nElements                  !< Input: number of points (typically nCells, nVertices or nEdges)

    real(kind=RKIND), dimension(3, 3, nElements), intent(out) ::  &
         transLocalToGlobal         !< Output: 3 x 3 matrix for transforming from local to global coordinates

    real(kind=RKIND), dimension(3, 3, nElements), intent(out) ::  &
         transGlobalToLocal         !< Output: 3 x 3 matrix for transforming from global to local coordinates

    ! local variables

    integer :: iElement

    real(kind=RKIND), dimension(3,3) :: &
         unitVectorGlobal,       &  ! unit vectors for global coordinate system
         unitVectorLocal            ! unit vectors for local coordinate system
                                    ! (expressed in global coordinates)

    transLocalToGlobal(:,:,:) = 0.0_RKIND
    transGlobalToLocal(:,:,:) = 0.0_RKIND

    ! set unit vectors for the global coordinate system

    ! global x
    unitVectorGlobal(1,1) = 1.0_RKIND
    unitVectorGlobal(2,1) = 0.0_RKIND
    unitVectorGlobal(3,1) = 0.0_RKIND

    ! global y
    unitVectorGlobal(1,2) = 0.0_RKIND
    unitVectorGlobal(2,2) = 1.0_RKIND
    unitVectorGlobal(3,2) = 0.0_RKIND

    ! global z
    unitVectorGlobal(1,3) = 0.0_RKIND
    unitVectorGlobal(2,3) = 0.0_RKIND
    unitVectorGlobal(3,3) = 1.0_RKIND

    do iElement = 1, nElements

       ! compute unit vectors for the local coordinate system
       ! Note: If using a rotating Cartesian grid, then the input global vectors must already have been rotated.

       ! local vertical
       ! this is the same as the vector from the center of the sphere to the point on the surface
       unitVectorLocal(1,3) = xGlobal(iElement)
       unitVectorLocal(2,3) = yGlobal(iElement)
       unitVectorLocal(3,3) = zGlobal(iElement)
       call unit_vector_3d(unitVectorLocal(:,3))

       if (abs(xGlobal(iElement)**2 + yGlobal(iElement)**2) > eps11) then  ! not at a pole

          ! local east
          ! this is obtained by a 90-degree rotation of (xGlobal, yGlobal) in a plane of constant latitude

          unitVectorLocal(1,1) = -yGlobal(iElement)
          unitVectorLocal(2,1) =  xGlobal(iElement)
          unitVectorLocal(3,1) =  0.0_RKIND
          call unit_vector_3d(unitVectorLocal(:,1))

          ! local north
          ! this is found by taking the cross product of the local vertical and local east unit vectors
          call cross_product_3d(unitVectorLocal(:,3), unitVectorLocal(:,1), unitVectorLocal(:,2))

       else  ! xGlobal = yGlobal = 0; at a pole

          ! take local east and local north to lie in a plane parallel to the equator, as in unitVectorGlobal
          if (zGlobal(iElement) > 0.0_RKIND) then  ! North Pole
             unitVectorLocal(:,1) = unitVectorGlobal(:,1)
             unitVectorLocal(:,2) = unitVectorGlobal(:,2)
          else   ! South Pole (reverse the order of east and north to obtain right-handed coordinates)
             unitVectorLocal(:,1) = unitVectorGlobal(:,2)
             unitVectorLocal(:,2) = unitVectorGlobal(:,1)
          endif

       endif   ! at the pole or not

       ! compute matrix for transforming from local to global coordinates
       transLocalToGlobal(:,:,iElement) = matmul(transpose(unitVectorGlobal), unitVectorLocal)

       ! compute matrix for transforming from local to global coordinates
       transGlobalToLocal(:,:,iElement) = matmul(transpose(unitVectorLocal), unitVectorGlobal)

    enddo   ! iElement

  end subroutine define_local_to_global_transformations

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine unit_vector_3d
!
!> \brief  normalize a 3D vector to have unit length
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  Given a 3D vector, return a vector with unit length
!
!-----------------------------------------------------------------------

  subroutine unit_vector_3d(vector)

    real(kind=RKIND), dimension(3), intent(inout) ::  &
         vector        !< Input/output: on input, any 3D Cartesian vector
                       !                on output, a vector of unit length

    ! local variables

    real(kind=RKIND) :: magnitude

    magnitude = sqrt(vector(1)**2 + vector(2)**2 + vector(3)**2)

    if (magnitude > 0.0_RKIND) then
       vector(:) = vector(:) / magnitude
    else
       vector(:) = 0.0_RKIND
    endif

  end subroutine unit_vector_3d

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine cross_product_3d
!
!> \brief  compute the cross product of two 3D vectors
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  Given two 3D vectors, return the cross product
!
!-----------------------------------------------------------------------

  subroutine cross_product_3d(&
       vector1, vector2,      &
       vectorOut)

    real(kind=RKIND), dimension(3), intent(in) ::  &
         vector1, vector2        !< Input: x/y/z coordinates of two 3D Cartesian vectorx

    real(kind=RKIND), dimension(3), intent(out) ::  &
         vectorOut               !< Output: x/y/z coordinates of the cross product

    vectorOut(1) = vector1(2)*vector2(3) - vector1(3)*vector2(2)
    vectorOut(2) = vector1(3)*vector2(1) - vector1(1)*vector2(3)
    vectorOut(3) = vector1(1)*vector2(2) - vector1(2)*vector2(1)

  end subroutine cross_product_3d

end module mpas_planar_projection
