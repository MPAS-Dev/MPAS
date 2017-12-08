! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_cdg_types
!
!> \brief MPAS ocean derived data types
!> \author Dave Lee
!> \date   December 2015
!> \details
!>  This module contains derived data types for the quadrature and basis
!>	functions. These contain local data and pointers to global fields 
!>	which are assigned from the pools.
!
!-----------------------------------------------------------------------

module ocn_cdg_types

   use mpas_kind_types
   use mpas_pool_routines

   implicit none
   private

   type :: cdgQuadrature
      character (len=10) :: fieldName
      integer, pointer :: order
      integer, pointer :: nq
      real (kind=RKIND), dimension(:), allocatable :: qw
      real (kind=RKIND), dimension(:), allocatable :: qx
      real (kind=RKIND), dimension(:), allocatable :: qy
      real (kind=RKIND), dimension(:), allocatable :: qz
      real (kind=RKIND), dimension(:,:,:), pointer :: xQuadL
      real (kind=RKIND), dimension(:,:,:), pointer :: yQuadL
      real (kind=RKIND), dimension(:,:,:), pointer :: xQuadG
      real (kind=RKIND), dimension(:,:,:), pointer :: yQuadG
      real (kind=RKIND), dimension(:,:,:), pointer :: zQuadG
   end type cdgQuadrature

   type, extends(cdgQuadrature) :: cdgQuad3D
      integer, pointer :: order_horiz
      integer, pointer :: order_vert
      integer, pointer :: nh
      integer, pointer :: nv
      real (kind=RKIND), dimension(:), allocatable :: qwv
      real (kind=RKIND), dimension(:), allocatable :: qzv
   end type cdgQuad3D

   type :: cdgBasis
      character (len=10) :: fieldName
      integer, pointer :: order
      integer, pointer :: nf
      integer, dimension(:), allocatable :: power_x
      integer, dimension(:), allocatable :: power_y
      real (kind=RKIND), dimension(:,:), pointer :: factor
      real (kind=RKIND), dimension(:,:), pointer :: mean
      real (kind=RKIND), dimension(:,:,:), pointer :: matInv
   end type cdgBasis

   type :: cdgBasis3D
      integer, pointer :: order_horiz
      integer, pointer :: order_vert
      integer, pointer :: nf
      !integer, dimension(:), allocatable :: power_x
      !integer, dimension(:), allocatable :: power_y
      !integer, dimension(:), allocatable :: power_z
      integer, dimension(:), pointer :: power_x
      integer, dimension(:), pointer :: power_y
      integer, dimension(:), pointer :: power_z
      real (kind=RKIND), dimension(:,:,:), pointer :: factor
      real (kind=RKIND), dimension(:,:,:), pointer :: mean_cur
      real (kind=RKIND), dimension(:,:,:), pointer :: mean_new
   end type cdgBasis3D

   public :: cdgQuadrature, cdgQuad3D, cdgBasis, cdgBasis3D
   public :: cdg_quadrature_populate, cdg_basis_populate, cdg_basis_populate_3d

   contains

   subroutine cdg_quadrature_populate(cdgPool, quad)!{{{

      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuadrature), intent(inout) :: quad

      call mpas_pool_get_array(cdgPool, 'xQuad_'//quad % fieldName, quad % xQuadL)
      call mpas_pool_get_array(cdgPool, 'yQuad_'//quad % fieldName, quad % yQuadL)
      call mpas_pool_get_array(cdgPool, 'xQuadGlobal_'//quad % fieldName, quad % xQuadG)
      call mpas_pool_get_array(cdgPool, 'yQuadGlobal_'//quad % fieldName, quad % yQuadG)
      call mpas_pool_get_array(cdgPool, 'zQuadGlobal_'//quad % fieldName, quad % zQuadG)

   end subroutine cdg_quadrature_populate!}}}

   subroutine cdg_quad3d_populate(cdgPool, quad)!{{{

      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuad3D), intent(inout) :: quad

      call mpas_pool_get_array(cdgPool, 'xQuad_hq', quad % xQuadL)
      call mpas_pool_get_array(cdgPool, 'yQuad_hq', quad % yQuadL)
      call mpas_pool_get_array(cdgPool, 'xQuadGlobal_hq', quad % xQuadG)
      call mpas_pool_get_array(cdgPool, 'yQuadGlobal_hq', quad % yQuadG)
      call mpas_pool_get_array(cdgPool, 'zQuadGlobal_hq', quad % zQuadG)

   end subroutine cdg_quad3d_populate!}}}

   subroutine cdg_basis_populate(cdgPool, basis)!{{{

      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgBasis), intent(inout) :: basis

      call mpas_pool_get_array(cdgPool, 'factor_'//basis % fieldName, basis % factor)
      call mpas_pool_get_array(cdgPool, 'mean_'//basis % fieldName, basis % mean)
      call mpas_pool_get_array(cdgPool, 'matInv_'//basis % fieldName, basis % matInv)

   end subroutine cdg_basis_populate!}}}

   subroutine cdg_basis_populate_3d(cdgPool, basis)!{{{

      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgBasis3D), intent(inout) :: basis

      call mpas_pool_get_dimension(cdgPool, 'nBasisFuncs3D', basis % nf)  
      call mpas_pool_get_array(cdgPool, 'factor3D', basis % factor)
      call mpas_pool_get_array(cdgPool, 'meanCur3D', basis % mean_cur)
      call mpas_pool_get_array(cdgPool, 'meanNew3D', basis % mean_new)

   end subroutine cdg_basis_populate_3d!}}}

end module ocn_cdg_types

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
