! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

module ocn_cdg_basis

   use mpas_timer
   use mpas_constants
   
   use ocn_constants

   use mpas_planar_projection
   use mpas_matrix_operations
   use mpas_geometry_utils

   use ocn_cdg_types
   use ocn_cdg_triangle

   implicit none
   private 
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_cdg_basis_init,             &
             ocn_cdg_basis_init_3d,          &
             ocn_cdg_basis_init_mean,        &
             ocn_cdg_basis_init_mean_3d,     &
             ocn_cdg_basis_init_factor_3d,   &
             ocn_cdg_basis_eval,             &
             ocn_cdg_basis_eval_i,           &
             ocn_cdg_basis_eval_3d,          &
             ocn_cdg_basis_eval_3d_i,        &
             ocn_cdg_basis_test

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_cdg_basis_init
!
!> \brief   Initialize the basis for a given cell
!> \author  Dave Lee
!> \date    June 2015
!> \details 
!>  This routine initializes the basis for a given cell.
!>  For its input arguments it takes the basis order, and the cell width.
!>  It then generates the powers, scaling factor and mean value for each term in
!>  the basis.
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_basis_init(domain, basis, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain
      type (cdgBasis), intent(inout) :: basis
      integer, intent(out) :: err

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

      call mpas_pool_get_config(domain % configs, 'config_cdg_basis_order_horiz', basis % order)
      call mpas_pool_get_config(domain % configs, 'config_cdg_num_basis_funcs_horiz', basis % nf)

      ! Initalize the private module variables
      allocate(basis % power_x(basis % nf), basis % power_y(basis % nf))

      select case (basis % order)
         case (1)
            basis % power_x = (/0/)
            basis % power_y = (/0/)
         case (2)
            basis % power_x = (/0, 1, 0/)
            basis % power_y = (/0, 0, 1/)
         case (3)
            basis % power_x = (/0, 1, 0, 2, 1, 0/)
            basis % power_y = (/0, 0, 1, 0, 1, 2/)
         case (4)
            basis % power_x = (/0, 1, 0, 2, 1, 0, 3, 2, 1, 0/)
            basis % power_y = (/0, 0, 1, 0, 1, 2, 0, 1, 2, 3/)
         case (5)
            basis % power_x = (/0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0/)
            basis % power_y = (/0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4/)
         end select

      call mpas_log_write( 'cdg: basis init $i $i', intArgs=(/ basis % order, basis % nf /) )
      call mpas_log_write( 'cdg: basis power x $i', intArgs=basis % power_x )
      call mpas_log_write( 'cdg: basis power y $i', intArgs=basis % power_y )

   end subroutine ocn_cdg_basis_init!}}}

   subroutine ocn_cdg_basis_init_3d(domain, basis, err)!{{{

      type (domain_type), intent(in) :: domain
      type (cdgBasis3D), intent(inout) :: basis
      integer, intent(out) :: err

      err = 0

      call mpas_pool_get_config(domain % configs, 'config_cdg_basis_order_vert', basis % order_vert)
      call mpas_pool_get_config(domain % configs, 'config_cdg_basis_order_horiz', basis % order_horiz)
      call mpas_pool_get_config(domain % configs, 'config_cdg_num_basis_funcs_3d', basis % nf)

      allocate(basis % power_x(basis % nf), basis % power_y(basis % nf), basis % power_z(basis % nf))

      if (basis % order_horiz .eq. 1 .and. basis % order_vert .eq. 1) then
         basis % power_x = (/0/)
         basis % power_y = (/0/)
         basis % power_z = (/0/)
      else if (basis % order_horiz .eq. 1 .and. basis % order_vert .eq. 2) then
         basis % power_x = (/0, 0/)
         basis % power_y = (/0, 0/)
         basis % power_z = (/0, 1/)
      else if (basis % order_horiz .eq. 1 .and. basis % order_vert .eq. 3) then
         basis % power_x = (/0, 0, 0/)
         basis % power_y = (/0, 0, 0/)
         basis % power_z = (/0, 1, 2/)
      else if (basis % order_horiz .eq. 2 .and. basis % order_vert .eq. 1) then
         basis % power_x = (/0, 1, 0/)
         basis % power_y = (/0, 0, 1/)
         basis % power_z = (/0, 0, 0/)
      else if (basis % order_horiz .eq. 2 .and. basis % order_vert .eq. 2) then
         basis % power_x = (/0, 1, 0, 0/)
         basis % power_y = (/0, 0, 1, 0/)
         basis % power_z = (/0, 0, 0, 1/)
      else if (basis % order_horiz .eq. 2 .and. basis % order_vert .eq. 3) then
         basis % power_x = (/0, 1, 0, 0, 0/)
         basis % power_y = (/0, 0, 1, 0, 0/)
         basis % power_z = (/0, 0, 0, 1, 2/)
      else if (basis % order_horiz .eq. 3 .and. basis % order_vert .eq. 1) then
         basis % power_x = (/0, 1, 0, 2, 1, 0/)
         basis % power_y = (/0, 0, 1, 0, 1, 2/)
         basis % power_z = (/0, 0, 0, 0, 0, 0/)
      else if (basis % order_horiz .eq. 3 .and. basis % order_vert .eq. 2) then
         basis % power_x = (/0, 1, 0, 2, 1, 0, 0/)
         basis % power_y = (/0, 0, 1, 0, 1, 2, 0/)
         basis % power_z = (/0, 0, 0, 0, 0, 0, 1/)
      else if (basis % order_horiz .eq. 3 .and. basis % order_vert .eq. 3) then
         basis % power_x = (/0, 1, 0, 2, 1, 0, 0, 0/)
         basis % power_y = (/0, 0, 1, 0, 1, 2, 0, 0/)
         basis % power_z = (/0, 0, 0, 0, 0, 0, 1, 2/)
      end if
      
      call mpas_log_write( 'cdg: basis init $i $i $i', intArgs=(/ basis % order_horiz, basis % order_vert, basis % nf /) )
      call mpas_log_write( 'cdg: basis power x $i', intArgs=basis % power_x )
      call mpas_log_write( 'cdg: basis power y $i', intArgs=basis % power_y )
      call mpas_log_write( 'cdg: basis power z $i', intArgs=basis % power_z )

   end subroutine ocn_cdg_basis_init_3d!}}}

   subroutine ocn_cdg_basis_init_factor_3d(meshPool, statePool, cdgPool, basis, err)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgBasis3D), intent(inout) :: basis
      integer, intent(out) :: err

      integer :: iCell, iLevel, iBasis
      real (kind=RKIND) :: dxInv, cellArea, dzInv

      integer, pointer :: nCells
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: triAreas, layerThickness

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

      ! Calculate the normalized triangle areas
      do iCell = 1, nCells
         cellArea = sum(triAreas(1:nEdgesOnCell(iCell),iCell))

         ! Calculate the basis scaling factors
         dxInv = sqrt(1.0_RKIND/cellArea)
         do iLevel = 1, maxLevelCell(iCell)
            dzInv = 1.0_RKIND/layerThickness(iLevel,iCell)
            basis % factor(1,iLevel,iCell) = 1.0_RKIND
            do iBasis = 2, basis % nf
               basis % factor(iBasis,iLevel,iCell) = ( dxInv**basis % power_x(iBasis) )*( dxInv**basis % power_y(iBasis) )*( dzInv**basis % power_z(iBasis) )
            end do
         end do
      end do
basis % factor(:,:,:) = 1.0_RKIND

   end subroutine ocn_cdg_basis_init_factor_3d!}}}

   subroutine ocn_cdg_basis_init_mean(meshPool, cdgPool, quad, basis, on_a_sphere, err)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuadrature), intent(in) :: quad
      type (cdgBasis), intent(inout) :: basis
      logical, intent(in) :: on_a_sphere
      integer, intent(out) :: err

      integer, pointer :: nCells, nEdges, nVertices
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell, edgesOnCell
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:,:), pointer :: triAreas

      integer :: iCell, iEdge, jEdge, iVertex, jVertex, iQuad, iBasis
      real (kind=RKIND) :: xl, yl
      real (kind=RKIND), dimension(:), allocatable :: weight

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)

      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

      err = 0

      allocate(weight(quad % nq))

      ! Initialize the mean for the basis funcs of this cell
      do iCell = 1, nCells
         basis % mean(:,iCell) = 0.0_RKIND

         do iBasis = 2, basis % nf
            do iEdge = 1, nEdgesOnCell(iCell)
               jEdge = mod(iEdge,nEdgesOnCell(iCell))+1
               iVertex = verticesOnCell(iEdge,iCell)
               jVertex = verticesOnCell(jEdge,iCell)

               weight(:) = triAreas(iEdge,iCell)*quad % qw(:)/sum(triAreas(1:nEdgesOnCell(iCell),iCell))

               do iQuad = 1, quad % nq
                  xl = quad % xQuadL(iQuad,iEdge,iCell)
                  yl = quad % yQuadL(iQuad,iEdge,iCell)

                  basis % mean(iBasis,iCell) = basis % mean(iBasis,iCell) + &
                                               weight(iQuad)*( xl**basis % power_x(iBasis) )*( yl**basis % power_y(iBasis) )
               end do
            end do ! iEdge
         end do ! iBasis
      end do ! iCell

      deallocate(weight)

   end subroutine ocn_cdg_basis_init_mean!}}}

   subroutine ocn_cdg_basis_init_mean_3d(meshPool, statePool, cdgPool, quad, basis, timeLev)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuad3D), intent(in) :: quad
      type (cdgBasis3D), intent(inout) :: basis
      integer, intent(in) :: timeLev

      integer, pointer :: nCells
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: triAreas, layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: mean

      integer :: iCell, iEdge, iQuad, jQuad, iBasis, iLevel
      real (kind=RKIND) :: xl, yl, zl, wt
      real (kind=RKIND), dimension(:), allocatable :: weight_h

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

      allocate(weight_h(quad % nh))

      if (timeLev .eq. 1) then
         mean => basis % mean_cur
         call mpas_pool_get_array(cdgPool, 'layerThicknessOld', layerThickness)
      else if (timeLev .eq. 2) then
         mean => basis % mean_new
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      end if
      mean(:,:,:) = 0.0_RKIND

      ! Initialize the mean for the basis funcs of this cell
      do iCell = 1, nCells
         do iLevel = 1, maxLevelCell(iCell)
            do iBasis = 2, basis % nf
               do iEdge = 1, nEdgesOnCell(iCell)
                  weight_h(:) = triAreas(iEdge,iCell)*quad % qw(:)/sum(triAreas(1:nEdgesOnCell(iCell),iCell))

                  do iQuad = 1, quad % nh
                     xl = quad % xQuadL(iQuad,iEdge,iCell)
                     yl = quad % yQuadL(iQuad,iEdge,iCell)

                     do jQuad = 1, quad % nv
                        wt = weight_h(iQuad)*quad % qwv(jQuad)
                        zl = layerThickness(iLevel,iCell)*(quad % qzv(jQuad) - 0.5_RKIND)

                        mean(iBasis,iLevel,iCell) = mean(iBasis,iLevel,iCell) + &
                           wt*( xl**basis % power_x(iBasis) )*( yl**basis % power_y(iBasis) )*( zl**basis % power_z(iBasis) )
                     end do
                  end do
               end do ! iEdge
            end do ! iBasis
         end do ! iLevel
      end do ! iCell

      deallocate(weight_h)

   end subroutine ocn_cdg_basis_init_mean_3d!}}}

!***********************************************************************
!
!  routine ocn_cdg_basis_eval_i
!
!> \brief   Evaluate the ith basis function at a given x,y coordinate
!> \author  Dave Lee
!> \date    June 2015
!> \details 
!
!-----------------------------------------------------------------------

   function ocn_cdg_basis_eval_i( basis, x, y, scale_i, mean_i, i ) result(phi_i) !{{{

      type (cdgBasis), intent(in) :: basis
      real (kind=RKIND), intent(in) :: x, y, scale_i, mean_i
      integer, intent(in) :: i

      real (kind=RKIND) :: phi_i

      phi_i = scale_i*( (x**basis % power_x(i))*(y**basis % power_y(i)) - mean_i)

   end function ocn_cdg_basis_eval_i!}}}

   function ocn_cdg_basis_eval_3d_i( basis, x, y, z, scale_i, mean_i, i ) result(phi_i) !{{{

      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), intent(in) :: x, y, z, scale_i, mean_i
      integer, intent(in) :: i

      real (kind=RKIND) :: phi_i

      phi_i = scale_i*( (x**basis % power_x(i))*(y**basis % power_y(i))*(z**basis % power_z(i)) - mean_i)

   end function ocn_cdg_basis_eval_3d_i!}}}

!***********************************************************************
!
!  routine ocn_cdg_basis_eval
!
!> \brief   Evaluate the full basis at a given x,y coord
!> \author  Dave Lee
!> \date    June 2015
!> \details 
!
!-----------------------------------------------------------------------

   function ocn_cdg_basis_eval( basis, x, y, factor, mean, ci ) result(phi) !{{{

      type (cdgBasis), intent(in) :: basis
      real (kind=RKIND), intent(in) :: x, y
      real (kind=RKIND), dimension(:), intent(in) :: factor, mean, ci

      integer :: iBasis
      real (kind=RKIND) :: phi

      phi = 0.0_RKIND

      do iBasis = 1, basis % nf
         phi = phi + ci(iBasis)*factor(iBasis)*( (x**basis % power_x(iBasis))*(y**basis % power_y(iBasis)) - mean(iBasis))
      end do

   end function ocn_cdg_basis_eval!}}}

   function ocn_cdg_basis_eval_3d( basis, x, y, z, factor, mean, ci ) result(phi) !{{{

      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), intent(in) :: x, y, z
      real (kind=RKIND), dimension(:), intent(in) :: factor, mean, ci

      integer :: iBasis
      real (kind=RKIND) :: phi

      phi = 0.0_RKIND

      do iBasis = 1, basis % nf
         phi = phi + ci(iBasis)*factor(iBasis)*( (x**basis % power_x(iBasis))*(y**basis % power_y(iBasis))*(z**basis % power_z(iBasis)) - mean(iBasis))
      end do

   end function ocn_cdg_basis_eval_3d!}}}

!***********************************************************************
!
!  routine ocn_cdg_basis_test
!
!> \brief   Unit test for the cdg basis functions
!> \author  Dave Lee
!> \date    July 2015
!> \details Ensure that the higher order (>1) terms are massless
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_basis_test(domain, quad, basis, err)!{{{

      type (domain_type), intent(in) :: domain
      type (cdgQuadrature), intent(in) :: quad
      type (cdgBasis), intent(in) :: basis
      integer, intent(out) :: err

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: cdgPool

      integer, pointer :: nCellsSolve
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell
      real (kind=RKIND), pointer :: radius
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:,:), pointer :: triAreas

      integer :: iCell, iEdge, jEdge, iBasis, iQuad, iVertex, jVertex
      real (kind=RKIND) :: detJ, xl, yl
      real (kind=RKIND), dimension(3) :: p1, p2, p3
      real (kind=RKIND), dimension(64) :: iBasisArea, basisAreaMinLocal, basisAreaMinGlobal, basisAreaMaxLocal, basisAreaMaxGlobal
      real (kind=RKIND), dimension(:), allocatable :: weight

      err = 0

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'cdg', cdgPool)

         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_config(meshPool, 'sphere_radius', radius)

         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)

         call mpas_pool_get_array(meshPool, 'xCell', xCell)
         call mpas_pool_get_array(meshPool, 'yCell', yCell)
         call mpas_pool_get_array(meshPool, 'zCell', zCell)
         call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
         call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
         call mpas_pool_get_array(meshPool, 'zVertex', zVertex)

         call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

         allocate(weight(quad % nq))

         basisAreaMinLocal(:) = +1.0E+9_RKIND
         basisAreaMaxLocal(:) = -1.0E+9_RKIND

         do iCell = 1, nCellsSolve
            iBasisArea(:) = 0.0_RKIND

            do iEdge = 1, nEdgesOnCell(iCell)
               jEdge = mod(iEdge,nEdgesOnCell(iCell))+1
               iVertex = verticesOnCell(iEdge,iCell)
               jVertex = verticesOnCell(jEdge,iCell)

               if (.not. domain % on_a_sphere) then
                  p1(1) = xVertex(iVertex) - xCell(iCell)
                  p1(2) = yVertex(iVertex) - yCell(iCell)
                  p2(1) = xVertex(jVertex) - xCell(iCell)
                  p2(2) = yVertex(jVertex) - yCell(iCell)
                  p3(:) = 0.0_RKIND
                  call ocn_cdg_tri_quad_pts(quad, p1, p2, p3, quad % xQuadL(:,iEdge,iCell), quad % yQuadL(:,iEdge,iCell))
               end if

               weight(:) = triAreas(iEdge,iCell)*quad % qw(:)

               do iQuad = 1, quad % nq
                  xl = quad % xQuadL(iQuad,iEdge,iCell)
                  yl = quad % yQuadL(iQuad,iEdge,iCell)

                  do iBasis = 2, basis % nf
                     iBasisArea(iBasis) = iBasisArea(iBasis) + &
                        weight(iQuad)*ocn_cdg_basis_eval_i(basis, xl, yl, basis % factor(iBasis,iCell), basis % mean(iBasis,iCell), iBasis)
                  end do
               end do
            end do ! iEdge

            do iBasis = 2, basis % nf
               if (iBasisArea(iBasis) .lt. basisAreaMinLocal(iBasis)) basisAreaMinLocal(iBasis) = iBasisArea(iBasis)
               if (iBasisArea(iBasis) .gt. basisAreaMaxLocal(iBasis)) basisAreaMaxLocal(iBasis) = iBasisArea(iBasis)
            end do
         end do ! iCell

         deallocate(weight)

         block => block % next
      end do

      do iBasis = 2, basis % nf
        call mpas_dmpar_min_real(domain % dminfo, basisAreaMinLocal(iBasis), basisAreaMinGlobal(iBasis))
        call mpas_dmpar_max_real(domain % dminfo, basisAreaMaxLocal(iBasis), basisAreaMaxGlobal(iBasis))
      end do
      write(0,*) 'cdg: basis min', basisAreaMinGlobal(2:basis % nf)
      write(0,*) 'cdg: basis max', basisAreaMaxGlobal(2:basis % nf)

   end subroutine ocn_cdg_basis_test!}}}

   subroutine ocn_cdg_basis_gs_proj( basis, fu, mu, fv, mv, qx, qy, qw, gs )!{{{

      type(cdgBasis), intent(in) :: basis
      real (kind=RKIND), dimension(:), intent(in) :: fu, mu, fv, mv, qx, qy, qw
      real (kind=RKIND), dimension(:), intent(out) :: gs

      integer :: i
      real (kind=RKIND) :: vu, uu, vi, ui, fac

      vu = 0.0_RKIND
      uu = 0.0_RKIND

      do i = 1, basis % nf
         vi = ocn_cdg_basis_eval_i(basis, qx(i), qy(i), fv(i), mv(i), i)
         ui = ocn_cdg_basis_eval_i(basis, qx(i), qy(i), fu(i), mu(i), i)

         vu = vu + qw(i)*vi*ui
         uu = uu + qw(i)*ui*ui
      end do

      fac = vu/uu

      do i = 1, basis % nf
         ui = ocn_cdg_basis_eval_i(basis, qx(i), qy(i), fu(i), mu(i), i)
         gs(i) = fac*ui
      end do

   end subroutine ocn_cdg_basis_gs_proj!}}}

!***********************************************************************

end module ocn_cdg_basis

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
