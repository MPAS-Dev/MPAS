! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!#define NEW

module ocn_cdg_limiter

   use mpas_timer
   use mpas_constants
   
   use ocn_constants

   use mpas_planar_projection
   use mpas_geometry_utils

   use ocn_cdg_types
   use ocn_cdg_triangle
   use ocn_cdg_basis

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

   public :: ocn_cdg_limiter_bj3d,        &
             ocn_cdg_limiter_bj3d_planar, &
             ocn_cdg_limiter_bj3d_split_planar, &
             ocn_cdg_limiter_weno

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

   subroutine ocn_cdg_limiter_bj3d(meshPool, statePool, cdgPool, basis, coeffs)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: coeffs

      integer, pointer :: nCells, nCellsSolve, nVertices, nVertLevels, nTracers
      integer, dimension(:), pointer :: maxLevelCell, maxLevelVertexBot, nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnVertex
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: transGtoC

      integer :: iCell, jCell, kCell, iLevel, iVertex, jVertex, iTracer, iBasis
      integer, dimension(3) :: firstOrd
      real (kind=RKIND) :: zl, phiAtVertex, phiDiff, alpha_i, alpha_j
      real (kind=RKIND), dimension(2) :: xl

      real (kind=RKIND), dimension(:), allocatable :: alpha
      real (kind=RKIND), dimension(:,:,:), allocatable :: phiMax, phiMin

      if (basis % nf .eq. 1) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
      call mpas_pool_get_array(cdgPool, 'transGlobalToCell', transGtoC)

      firstOrd(1) = 2
      firstOrd(2) = 3
      firstOrd(3) = 4
!      if (basis % nf .gt. 5) then ! 3rd order in vertical and horizontal
!         firstOrd(1) = 2
!         firstOrd(2) = 3
!         firstOrd(3) = 7
!         call ocn_cdg_limiter_bj3d_planar_gradients(meshPool, statePool, cdgPool, basis, coeffs)
!      else if (basis % nf .gt. 4) then ! 3rd order in vertical only
      if (basis % nf .eq. 5) then ! 3rd order in vertical only
         call ocn_cdg_limiter_bj3d_vert_grad(meshPool, statePool, cdgPool, basis, coeffs)
      end if

      allocate(alpha(nTracers))
      allocate(phiMax(nTracers,nVertLevels+1,nVertices), phiMin(nTracers,nVertLevels+1,nVertices))

      ! Compute that max and min mean element values for the elements incident on each vertex
      phiMax(:,:,:) = -1.0E+9_RKIND
      phiMin(:,:,:) = +1.0E+9_RKIND
      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            do iVertex = 1, nEdgesOnCell(iCell)
               jVertex = verticesOnCell(iVertex,iCell)
               do iTracer = 1, nTracers
                  if (coeffs(1,iTracer,iLevel,iCell) .gt. phiMax(iTracer,iLevel,  jVertex)) phiMax(iTracer,iLevel,  jVertex) = coeffs(1,iTracer,iLevel,iCell)
                  if (coeffs(1,iTracer,iLevel,iCell) .lt. phiMin(iTracer,iLevel,  jVertex)) phiMin(iTracer,iLevel,  jVertex) = coeffs(1,iTracer,iLevel,iCell)
                  if (coeffs(1,iTracer,iLevel,iCell) .gt. phiMax(iTracer,iLevel+1,jVertex)) phiMax(iTracer,iLevel+1,jVertex) = coeffs(1,iTracer,iLevel,iCell)
                  if (coeffs(1,iTracer,iLevel,iCell) .lt. phiMin(iTracer,iLevel+1,jVertex)) phiMin(iTracer,iLevel+1,jVertex) = coeffs(1,iTracer,iLevel,iCell)
               end do
            end do
         end do
      end do

      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            alpha(:) = 1.0_RKIND

            ! Determine limiter coefficients for the top face vertices
            zl = +0.5_RKIND*layerThickness(iLevel,iCell)
            do iVertex = 1, nEdgesOnCell(iCell)
               jVertex = verticesOnCell(iVertex,iCell)

               call ocn_cdg_tri_trans_g2l(transGtoC(:,:,iCell), xCell(iCell), yCell(iCell), zCell(iCell), &
                  xVertex(jVertex), yVertex(jVertex), zVertex(jVertex), xl(:))

               do iTracer = 1, nTracers
                  phiAtVertex = ocn_cdg_basis_eval_3d(basis, xl(1), xl(2), zl, basis % factor(:,iLevel,iCell), &
                     basis % mean_new(:,iLevel,iCell), coeffs(:,iTracer,iLevel,iCell))
                  phiDiff = phiAtVertex - coeffs(1,iTracer,iLevel,iCell)

                  if (abs(phiDiff) .lt. 1.0E-10_RKIND) then
                     alpha_j = 1.0_RKIND
                  else if (phiDiff .gt. 0.0_RKIND) then
                     alpha_i = (phiMax(iTracer,iLevel,jVertex) - coeffs(1,iTracer,iLevel,iCell))/phiDiff
                     alpha_j = min(1.0_RKIND,alpha_i)
                  else if (phiDiff .lt. 0.0_RKIND) then
                     alpha_i = (phiMin(iTracer,iLevel,jVertex) - coeffs(1,iTracer,iLevel,iCell))/phiDiff
                     alpha_j = min(1.0_RKIND,alpha_i)
                  end if
                  alpha(iTracer) = min(alpha(iTracer),alpha_j)
               end do
            end do

            ! Determine limiter coefficients for the bottom face vertices
            zl = -0.5_RKIND*layerThickness(iLevel,iCell)
            do iVertex = 1, nEdgesOnCell(iCell)
               jVertex = verticesOnCell(iVertex,iCell)

               call ocn_cdg_tri_trans_g2l(transGtoC(:,:,iCell), xCell(iCell), yCell(iCell), zCell(iCell), &
                  xVertex(jVertex), yVertex(jVertex), zVertex(jVertex), xl(:))

               do iTracer = 1, nTracers
                  phiAtVertex = ocn_cdg_basis_eval_3d(basis, xl(1), xl(2), zl, basis % factor(:,iLevel,iCell), &
                     basis % mean_new(:,iLevel,iCell), coeffs(:,iTracer,iLevel,iCell))
                  phiDiff = phiAtVertex - coeffs(1,iTracer,iLevel,iCell)

                  if (abs(phiDiff) .lt. 1.0E-10_RKIND) then
                     alpha_j = 1.0_RKIND
                  else if (phiDiff .gt. 0.0_RKIND) then
                     alpha_i = (phiMax(iTracer,iLevel+1,jVertex) - coeffs(1,iTracer,iLevel,iCell))/phiDiff
                     alpha_j = min(1.0_RKIND,alpha_i)
                  else if (phiDiff .lt. 0.0_RKIND) then
                     alpha_i = (phiMin(iTracer,iLevel+1,jVertex) - coeffs(1,iTracer,iLevel,iCell))/phiDiff
                     alpha_j = min(1.0_RKIND,alpha_i)
                  end if
                  alpha(iTracer) = min(alpha(iTracer),alpha_j)
               end do
            end do

            ! Apply the scaling factor to all higher order coefficients
            do iTracer = 1, nTracers
               !do iBasis = 1, 3
               !do iBasis = 1, basis % nf - 1
               !   coeffs(firstOrd(iBasis),iTracer,iLevel,iCell) = alpha(iTracer)*coeffs(firstOrd(iBasis),iTracer,iLevel,iCell)
               !end do
               coeffs(2,iTracer,iLevel,iCell) = alpha(iTracer)*coeffs(2,iTracer,iLevel,iCell)
               coeffs(3,iTracer,iLevel,iCell) = alpha(iTracer)*coeffs(3,iTracer,iLevel,iCell)
            end do
         end do ! iLevel
      end do ! iCell

      deallocate(alpha)
      deallocate(phiMax, phiMin)

   end subroutine ocn_cdg_limiter_bj3d!}}}

#ifndef NEW
   subroutine ocn_cdg_limiter_bj3d_planar(meshPool, statePool, cdgPool, basis, coeffs)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: coeffs

      integer, pointer :: nCells, nCellsSolve, nVertices, nVertLevels, vertexDegree, nTracers
      integer, dimension(:), pointer :: maxLevelCell, maxLevelVertexBot, nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnVertex, cellsOnCell
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, xVertex, yVertex, bottomDepth
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness

      logical, dimension(10) :: troubled_x, troubled_y, troubled_z
      integer :: iCell, jCell, kCell, iLevel, iVertex, jVertex, iTracer, iBasis
      integer, dimension(3) :: firstOrd
      real (kind=RKIND) :: xl, yl, zl, phiAtVertex, phiDiff, alpha_i, alpha_j
      real (kind=RKIND), dimension(4,4) :: transGtoC

      real (kind=RKIND), dimension(:), allocatable :: alpha
      real (kind=RKIND), dimension(:,:), allocatable :: zMid
      real (kind=RKIND), dimension(:,:,:), allocatable :: phiMax, phiMin

      if (basis % nf .eq. 1) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)

      allocate(zMid(nVertLevels,nCells))
      do iCell = 1, nCells
         zMid(maxLevelCell(iCell),iCell) = 0.5_RKIND*layerThickness(maxLevelCell(iCell),iCell) - bottomDepth(iCell)
         do iLevel = maxLevelCell(iCell)-1, 1, -1
            zMid(iLevel,iCell)  = zMid(iLevel+1,iCell) + 0.5_RKIND*(layerThickness(iLevel+1,iCell) + layerThickness(iLevel,iCell))
         end do
      end do

      firstOrd(1) = 2
      firstOrd(2) = 3
      firstOrd(3) = 4
      if (basis % nf .gt. 5) then ! 3rd order in vertical and horizontal
         firstOrd(1) = 2
         firstOrd(2) = 3
         firstOrd(3) = 7
         call ocn_cdg_limiter_bj3d_planar_gradients(meshPool, statePool, cdgPool, basis, coeffs)
      else if (basis % nf .gt. 4) then ! 3rd order in vertical only
         call ocn_cdg_limiter_bj3d_vert_grad(meshPool, statePool, cdgPool, basis, coeffs)
      end if

      allocate(alpha(nTracers))
      allocate(phiMax(nTracers,nVertLevels+1,nVertices), phiMin(nTracers,nVertLevels+1,nVertices))

      ! Compute that max and min mean element values for the elements incident on each vertex
      phiMax(:,:,:) = -1.0E+9_RKIND
      phiMin(:,:,:) = +1.0E+9_RKIND
      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            do iVertex = 1, nEdgesOnCell(iCell)
               jVertex = verticesOnCell(iVertex,iCell)
               do iTracer = 1, nTracers
                  if (coeffs(1,iTracer,iLevel,iCell) .gt. phiMax(iTracer,iLevel,  jVertex)) phiMax(iTracer,iLevel,  jVertex) = coeffs(1,iTracer,iLevel,iCell)
                  if (coeffs(1,iTracer,iLevel,iCell) .lt. phiMin(iTracer,iLevel,  jVertex)) phiMin(iTracer,iLevel,  jVertex) = coeffs(1,iTracer,iLevel,iCell)
                  if (coeffs(1,iTracer,iLevel,iCell) .gt. phiMax(iTracer,iLevel+1,jVertex)) phiMax(iTracer,iLevel+1,jVertex) = coeffs(1,iTracer,iLevel,iCell)
                  if (coeffs(1,iTracer,iLevel,iCell) .lt. phiMin(iTracer,iLevel+1,jVertex)) phiMin(iTracer,iLevel+1,jVertex) = coeffs(1,iTracer,iLevel,iCell)
               end do
            end do
         end do
      end do

      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            call ocn_cdg_limiter_troubled_cell_horiz(iCell, iLevel, .false., transGtoC, xCell, yCell, yCell, &
               cellsOnCell(:,iCell), nEdgesOnCell(iCell), nTracers, basis % factor(:,iLevel,iCell), coeffs, troubled_x, troubled_y)
            call ocn_cdg_limiter_troubled_cell_vert(iCell, iLevel, zMid(:,iCell), maxLevelCell(iCell), nTracers, 4, &
               basis % factor(:,iLevel,iCell), coeffs, troubled_z)

            alpha(:) = 1.0_RKIND

            ! Determine limiter coefficients for the top face vertices
            zl = +0.5_RKIND*layerThickness(iLevel,iCell)
            do iVertex = 1, nEdgesOnCell(iCell)
               jVertex = verticesOnCell(iVertex,iCell)
               xl = xVertex(jVertex) - xCell(iCell)
               yl = yVertex(jVertex) - yCell(iCell)

               do iTracer = 1, nTracers
                  phiAtVertex = ocn_cdg_basis_eval_3d(basis, xl, yl, zl, basis % factor(:,iLevel,iCell), &
                     basis % mean_new(:,iLevel,iCell), coeffs(:,iTracer,iLevel,iCell))
                  phiDiff = phiAtVertex - coeffs(1,iTracer,iLevel,iCell)

                  if (abs(phiDiff) .lt. 1.0E-10_RKIND) then
                     alpha_j = 1.0_RKIND
                  else if (phiDiff .gt. 0.0_RKIND) then
                     alpha_i = (phiMax(iTracer,iLevel,jVertex) - coeffs(1,iTracer,iLevel,iCell))/phiDiff
                     alpha_j = min(1.0_RKIND,alpha_i)
                  else if (phiDiff .lt. 0.0_RKIND) then
                     alpha_i = (phiMin(iTracer,iLevel,jVertex) - coeffs(1,iTracer,iLevel,iCell))/phiDiff
                     alpha_j = min(1.0_RKIND,alpha_i)
                  end if
                  alpha(iTracer) = min(alpha(iTracer),alpha_j)
               end do
            end do

            ! Determine limiter coefficients for the bottom face vertices
            zl = -0.5_RKIND*layerThickness(iLevel,iCell)
            do iVertex = 1, nEdgesOnCell(iCell)
               jVertex = verticesOnCell(iVertex,iCell)
               xl = xVertex(jVertex) - xCell(iCell)
               yl = yVertex(jVertex) - yCell(iCell)

               do iTracer = 1, nTracers
                  phiAtVertex = ocn_cdg_basis_eval_3d(basis, xl, yl, zl, basis % factor(:,iLevel,iCell), &
                     basis % mean_new(:,iLevel,iCell), coeffs(:,iTracer,iLevel,iCell))
                  phiDiff = phiAtVertex - coeffs(1,iTracer,iLevel,iCell)

                  if (abs(phiDiff) .lt. 1.0E-10_RKIND) then
                     alpha_j = 1.0_RKIND
                  else if (phiDiff .gt. 0.0_RKIND) then
                     alpha_i = (phiMax(iTracer,iLevel+1,jVertex) - coeffs(1,iTracer,iLevel,iCell))/phiDiff
                     alpha_j = min(1.0_RKIND,alpha_i)
                  else if (phiDiff .lt. 0.0_RKIND) then
                     alpha_i = (phiMin(iTracer,iLevel+1,jVertex) - coeffs(1,iTracer,iLevel,iCell))/phiDiff
                     alpha_j = min(1.0_RKIND,alpha_i)
                  end if
                  alpha(iTracer) = min(alpha(iTracer),alpha_j)
               end do
            end do

            ! Apply the scaling factor to all higher order coefficients
            do iTracer = 1, nTracers
               !do iBasis = 1, 3
               !   coeffs(firstOrd(iBasis),iTracer,iLevel,iCell) = alpha(iTracer)*coeffs(firstOrd(iBasis),iTracer,iLevel,iCell)
               !end do
               if (troubled_x(iTracer)) then
                  coeffs(2,iTracer,iLevel,iCell) = alpha(iTracer)*coeffs(2,iTracer,iLevel,iCell)
               end if
               if (troubled_y(iTracer)) then
                  coeffs(3,iTracer,iLevel,iCell) = alpha(iTracer)*coeffs(3,iTracer,iLevel,iCell)
               end if
               if (troubled_z(iTracer)) then
                  coeffs(4,iTracer,iLevel,iCell) = alpha(iTracer)*coeffs(4,iTracer,iLevel,iCell)
               end if
            end do
         end do ! iLevel
      end do ! iCell

      deallocate(alpha)
      deallocate(phiMax, phiMin)
      deallocate(zMid)

   end subroutine ocn_cdg_limiter_bj3d_planar!}}}
#else
   subroutine ocn_cdg_limiter_bj3d_planar(meshPool, statePool, cdgPool, basis, coeffs)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: coeffs

      integer, pointer :: nCells, nCellsSolve, nVertices, nVertLevels, vertexDegree, nTracers
      integer, dimension(:), pointer :: maxLevelCell, maxLevelVertexBot, nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnVertex, cellsOnCell, cellsOnEdge, edgesOnCell
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, xVertex, yVertex, bottomDepth, xEdge, yEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness

      integer :: iCell, jCell, kCell, iLevel, iVertex, jVertex, iTracer, iBasis, iEdge, jEdge
      real (kind=RKIND) :: xl, yl, zl, phiDiff, alpha_i, alpha_j, phiMaxZ, phiMinZ

      real (kind=RKIND), dimension(:), allocatable :: alpha_x, alpha_y, alpha_z, mean_x, mean_z
      real (kind=RKIND), dimension(:,:,:), allocatable :: phiMax, phiMin

      if (basis % nf .eq. 1) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)

      allocate(alpha_x(nTracers), alpha_y(nTracers), alpha_z(nTracers))
      allocate(phiMax(nTracers,nVertLevels+1,nVertices), phiMin(nTracers,nVertLevels+1,nVertices))
      allocate(mean_x(basis % nf), mean_z(basis % nf))

      ! Compute that max and min mean element values for the elements incident on each vertex
      phiMax(:,:,:) = -1.0E+9_RKIND
      phiMin(:,:,:) = +1.0E+9_RKIND
      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            do iVertex = 1, nEdgesOnCell(iCell)
               jVertex = verticesOnCell(iVertex,iCell)
               do iTracer = 1, nTracers
                  if (coeffs(1,iTracer,iLevel,iCell) .gt. phiMax(iTracer,iLevel,jVertex)) phiMax(iTracer,iLevel,jVertex) = coeffs(1,iTracer,iLevel,iCell)
                  if (coeffs(1,iTracer,iLevel,iCell) .lt. phiMin(iTracer,iLevel,jVertex)) phiMin(iTracer,iLevel,jVertex) = coeffs(1,iTracer,iLevel,iCell)
               end do
            end do
         end do
      end do

      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            alpha_x(:) = basis % factor(2,iLevel,iCell)*coeffs(2,:,iLevel,iCell)
            alpha_y(:) = basis % factor(3,iLevel,iCell)*coeffs(3,:,iLevel,iCell)
            alpha_z(:) = 1.0_RKIND

            ! Determine limiter coefficients for the top of element
            do iTracer = 1, nTracers
               mean_x(:) = 0.0_RKIND
               mean_z(:) = 0.0_RKIND
               do iBasis = 1, basis % nf
                  if (basis % power_x(iBasis) + basis % power_y(iBasis) .gt. 0) mean_x(iBasis) = basis % mean_new(iBasis,iLevel,iCell)
                  if (basis % power_z(iBasis) .gt. 0) mean_z(iBasis) = basis % mean_new(iBasis,iLevel,iCell)
               end do

               if (iLevel .gt. 1 .and. iLevel .lt. maxLevelCell(iCell)) then
                  phiMaxZ = max(coeffs(1,iTracer,iLevel+1,iCell), coeffs(1,iTracer,iLevel,iCell), coeffs(1,iTracer,iLevel-1,iCell))
                  phiMinZ = min(coeffs(1,iTracer,iLevel+1,iCell), coeffs(1,iTracer,iLevel,iCell), coeffs(1,iTracer,iLevel-1,iCell))
               else if (iLevel .eq. 1) then
                  phiMaxZ = max(coeffs(1,iTracer,iLevel+1,iCell), coeffs(1,iTracer,iLevel,iCell))
                  phiMinZ = min(coeffs(1,iTracer,iLevel+1,iCell), coeffs(1,iTracer,iLevel,iCell))
               else
                  phiMaxZ = max(coeffs(1,iTracer,iLevel,iCell), coeffs(1,iTracer,iLevel-1,iCell))
                  phiMinZ = min(coeffs(1,iTracer,iLevel,iCell), coeffs(1,iTracer,iLevel-1,iCell))
               end if

               if (iLevel .gt. 1) then
                  zl = +0.5_RKIND*layerThickness(iLevel,iCell)
                  phiDiff = ocn_cdg_basis_eval_3d(basis, 0.0_RKIND, 0.0_RKIND, zl, basis % factor(:,iLevel,iCell), &
                     mean_z(:), coeffs(:,iTracer,iLevel,iCell)) - coeffs(1,iTracer,iLevel,iCell)

                  if (abs(phiDiff) .lt. 1.0E-10_RKIND) then
                     alpha_i = 1.0_RKIND
                  else if (phiDiff .gt. 0.0_RKIND) then
                     alpha_i = min(1.0_RKIND, (phiMaxZ - coeffs(1,iTracer,iLevel,iCell))/phiDiff)
                  else if (phiDiff .lt. 0.0_RKIND) then
                     alpha_i = min(1.0_RKIND, (phiMinZ - coeffs(1,iTracer,iLevel,iCell))/phiDiff)
                  end if
                  alpha_z(iTracer) = min(alpha_z(iTracer),alpha_i)
               end if

               if (iLevel .lt. maxLevelCell(iCell)) then
                  zl = -0.5_RKIND*layerThickness(iLevel,iCell)
                  phiDiff = ocn_cdg_basis_eval_3d(basis, 0.0_RKIND, 0.0_RKIND, zl, basis % factor(:,iLevel,iCell), &
                     mean_z(:), coeffs(:,iTracer,iLevel,iCell)) - coeffs(1,iTracer,iLevel,iCell)

                  if (abs(phiDiff) .lt. 1.0E-10_RKIND) then
                     alpha_i = 1.0_RKIND
                  else if (phiDiff .gt. 0.0_RKIND) then
                     alpha_i = min(1.0_RKIND, (phiMaxZ - coeffs(1,iTracer,iLevel,iCell))/phiDiff)
                  else if (phiDiff .lt. 0.0_RKIND) then
                     alpha_i = min(1.0_RKIND, (phiMinZ - coeffs(1,iTracer,iLevel,iCell))/phiDiff)
                  end if
                  alpha_z(iTracer) = min(alpha_z(iTracer),alpha_i)
               end if

               ! Apply the scaling factor to all higher order coefficients
               coeffs(4,iTracer,iLevel,iCell) = alpha_z(iTracer)*coeffs(4,iTracer,iLevel,iCell)
            end do !iTracer

            alpha_x(:) = basis % factor(2,iLevel,iCell)*coeffs(2,:,iLevel,iCell)
            alpha_y(:) = basis % factor(3,iLevel,iCell)*coeffs(3,:,iLevel,iCell)
            do iVertex = 1, nEdgesOnCell(iCell)
               jVertex = verticesOnCell(iVertex,iCell)
               xl = xVertex(jVertex) - xCell(iCell)
               yl = yVertex(jVertex) - yCell(iCell)

               do iTracer = 1, nTracers
                  call limit_anisotropic_horiz(xl, yl, coeffs(:,iTracer,iLevel,iCell)*basis % factor(:,iLevel,iCell), basis % mean_new(:,iLevel,iCell), &
                                               phiMin(iTracer,iLevel,jVertex), phiMax(iTracer,iLevel,jVertex), alpha_i, alpha_j)

                  if (abs(alpha_i) .lt. abs(alpha_x(iTracer))) alpha_x(iTracer) = alpha_i
                  if (abs(alpha_j) .lt. abs(alpha_y(iTracer))) alpha_y(iTracer) = alpha_j
               end do
            end do

            do iTracer = 1, nTracers
               !coeffs(2,iTracer,iLevel,iCell) = alpha_x(iTracer)*coeffs(2,iTracer,iLevel,iCell)
               !coeffs(3,iTracer,iLevel,iCell) = alpha_y(iTracer)*coeffs(3,iTracer,iLevel,iCell)
               coeffs(2,iTracer,iLevel,iCell) = alpha_x(iTracer)
               coeffs(3,iTracer,iLevel,iCell) = alpha_y(iTracer)
            end do
         end do ! iLevel
      end do ! iCell

      deallocate(alpha_x, alpha_y, alpha_z)
      deallocate(phiMax, phiMin)
      deallocate(mean_x, mean_z)

   end subroutine ocn_cdg_limiter_bj3d_planar!}}}
#endif
   subroutine ocn_cdg_limiter_bj3d_planar_gradients(meshPool, statePool, cdgPool, basis, coeffs)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: coeffs

      integer :: iVertex, jVertex, iCell, jCell, iLevel, iTracer, iBasis, jBasis
      integer, dimension(3) :: firstOrd
      real (kind=RKIND) :: alpha_i, xl, yl, zl, grad
      real (kind=RKIND), dimension(3) :: dx, dxDiff

      integer, pointer :: nCells, nCellsSolve, nVertices, nVertLevels, vertexDegree, nTracers
      integer, dimension(:), pointer :: maxLevelCell, maxLevelVertexBot, nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnVertex
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, xVertex, yVertex
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness

      real (kind=RKIND), dimension(:,:,:), pointer :: fac
      real (kind=RKIND), dimension(:,:,:,:), allocatable :: gradMax, gradMin, alpha_x

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)

      allocate(alpha_x(3,nTracers,nVertLevels,nCellsSolve))
      allocate(gradMax(3,nTracers,nVertLevels+1,nVertices), gradMin(3,nTracers,nVertLevels+1,nVertices))

      fac => basis % factor
      firstOrd(1) = 2
      firstOrd(2) = 3
      firstOrd(3) = 7

      ! Note: assumes 3rd order in all dimensions
      do iVertex = 1, nVertices
         do iLevel = 1, maxLevelVertexBot(iVertex) + 1
            gradMax(:,:,iLevel,iVertex) = -1.0E+9_RKIND
            gradMin(:,:,iLevel,iVertex) = +1.0E+9_RKIND

            do iCell = 1, vertexDegree
               jCell = cellsOnVertex(iCell,iVertex)
               if (jCell .le. nCells) then

                  if (iLevel .gt. 1 .and. iLevel-1 .le. maxLevelCell(jCell)) then
                     do iTracer = 1, nTracers
                        do iBasis = 1, 3
                           jBasis = firstOrd(iBasis)
                           grad = fac(jBasis,iLevel-1,jCell)*coeffs(jBasis,iTracer,iLevel-1,jCell)
                           if (grad .gt. gradMax(iBasis,iTracer,iLevel,iVertex)) gradMax(iBasis,iTracer,iLevel,iVertex) = grad
                           if (grad .lt. gradMin(iBasis,iTracer,iLevel,iVertex)) gradMin(iBasis,iTracer,iLevel,iVertex) = grad
                        end do
                     end do
                  end if

                  if (iLevel .lt. maxLevelVertexBot(iVertex)+1 .and. iLevel .le. maxLevelCell(jCell)) then
                     do iTracer = 1, nTracers
                        do iBasis = 1, 3
                           jBasis = firstOrd(iBasis)
                           grad = fac(jBasis,iLevel,jCell)*coeffs(jBasis,iTracer,iLevel,jCell)
                           if (grad .gt. gradMax(iBasis,iTracer,iLevel,iVertex)) gradMax(iBasis,iTracer,iLevel,iVertex) = grad
                           if (grad .lt. gradMin(iBasis,iTracer,iLevel,iVertex)) gradMin(iBasis,iTracer,iLevel,iVertex) = grad
                        end do
                     end do
                  end if
               end if
            end do
         end do
      end do

      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            alpha_x(:,:,iLevel,iCell) = 1.0_RKIND

            zl = +0.5_RKIND*layerThickness(iLevel,iCell)
            do iVertex = 1, nEdgesOnCell(iCell)
               jVertex = verticesOnCell(iVertex,iCell)
               xl = xVertex(jVertex) - xCell(iCell)
               yl = yVertex(jVertex) - yCell(iCell)

               do iTracer = 1, nTracers
                  do iBasis = 1, 3   
                     dx(iBasis) = fac(firstOrd(iBasis),iLevel,iCell)*coeffs(firstOrd(iBasis),iTracer,iLevel,iCell)
                  end do

                  dxDiff(1) = 2.0_RKIND*xl*fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell) + yl*fac(5,iLevel,iCell)*coeffs(5,iTracer,iLevel,iCell)
                  dxDiff(2) = 2.0_RKIND*yl*fac(6,iLevel,iCell)*coeffs(6,iTracer,iLevel,iCell) + xl*fac(5,iLevel,iCell)*coeffs(5,iTracer,iLevel,iCell)
                  dxDiff(3) = 2.0_RKIND*zl*fac(8,iLevel,iCell)*coeffs(8,iTracer,iLevel,iCell)

                  do iBasis = 1, 3
                     if (abs(dxDiff(iBasis)) .lt. 1.0E-10_RKIND) then
                        alpha_i = 1.0_RKIND
                     else if (dxDiff(iBasis) .gt. 0.0_RKIND) then
                        alpha_i = (gradMax(iBasis,iTracer,iLevel,jVertex) - dx(iBasis))/dxDiff(iBasis)
                     else if (dxDiff(iBasis) .lt. 0.0_RKIND) then
                        alpha_i = (gradMin(iBasis,iTracer,iLevel,jVertex) - dx(iBasis))/dxDiff(iBasis)
                     end if
                     alpha_x(iBasis,iTracer,iLevel,iCell) = min(alpha_x(iBasis,iTracer,iLevel,iCell),alpha_i)
                  end do
               end do
            end do

            zl = -0.5_RKIND*layerThickness(iLevel,iCell)
            do iVertex = 1, nEdgesOnCell(iCell)
               jVertex = verticesOnCell(iVertex,iCell)
               xl = xVertex(jVertex) - xCell(iCell)
               yl = yVertex(jVertex) - yCell(iCell)

               do iTracer = 1, nTracers
                  do iBasis = 1, 3   
                     dx(iBasis) = fac(firstOrd(iBasis),iLevel,iCell)*coeffs(firstOrd(iBasis),iTracer,iLevel,iCell)
                  end do

                  dxDiff(1) = 2.0_RKIND*xl*fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell) + yl*fac(5,iLevel,iCell)*coeffs(5,iTracer,iLevel,iCell)
                  dxDiff(2) = 2.0_RKIND*yl*fac(6,iLevel,iCell)*coeffs(6,iTracer,iLevel,iCell) + xl*fac(5,iLevel,iCell)*coeffs(5,iTracer,iLevel,iCell)
                  dxDiff(3) = 2.0_RKIND*zl*fac(8,iLevel,iCell)*coeffs(8,iTracer,iLevel,iCell)

                  do iBasis = 1, 3
                     if (abs(dxDiff(iBasis)) .lt. 1.0E-10_RKIND) then
                        alpha_i = 1.0_RKIND
                     else if (dxDiff(iBasis) .gt. 0.0_RKIND) then
                        alpha_i = (gradMax(iBasis,iTracer,iLevel+1,jVertex) - dx(iBasis))/dxDiff(iBasis)
                     else if (dxDiff(iBasis) .lt. 0.0_RKIND) then
                        alpha_i = (gradMin(iBasis,iTracer,iLevel+1,jVertex) - dx(iBasis))/dxDiff(iBasis)
                     end if
                     alpha_x(iBasis,iTracer,iLevel,iCell) = min(alpha_x(iBasis,iTracer,iLevel,iCell),alpha_i)
                  end do
               end do
            end do ! iTracer
         end do ! iLevel
      end do ! iCell

      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            do iTracer = 1, nTracers
               ! Limit the 3rd order horizontal components
               alpha_i = min(alpha_x(1,iTracer,iLevel,iCell),alpha_x(2,iTracer,iLevel,iCell))
               coeffs(4:6,iTracer,iLevel,iCell) = alpha_i*coeffs(4:6,iTracer,iLevel,iCell)

               ! Limit the 3rd order vertical component
               coeffs(8,iTracer,iLevel,iCell) = alpha_x(3,iTracer,iLevel,iCell)*coeffs(8,iTracer,iLevel,iCell)

!               alpha_i = minval(alpha_x(:,iTracer,iLevel,iCell))
!               coeffs(4,iTracer,iLevel,iCell) = alpha_i*coeffs(4,iTracer,iLevel,iCell)
!               coeffs(5,iTracer,iLevel,iCell) = alpha_i*coeffs(5,iTracer,iLevel,iCell)
!               coeffs(6,iTracer,iLevel,iCell) = alpha_i*coeffs(6,iTracer,iLevel,iCell)
!               coeffs(8,iTracer,iLevel,iCell) = alpha_i*coeffs(8,iTracer,iLevel,iCell)
            end do
         end do
      end do

      deallocate(alpha_x)
      deallocate(gradMin, gradMax)

   end subroutine ocn_cdg_limiter_bj3d_planar_gradients!}}}

   subroutine ocn_cdg_limiter_bj3d_vert_grad(meshPool, statePool, cdgPool, basis, coeffs)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: coeffs

      integer :: iCell, jCell, iLevel, iTracer, iBasis, jBasis
      real (kind=RKIND) :: alpha_i, alpha_j, zl, dz, dzDiff, dzVert, gradMax, gradMin

      integer, pointer :: nCells, nCellsSolve, nVertices, nVertLevels, vertexDegree, nTracers
      integer, dimension(:), pointer :: maxLevelCell, maxLevelVertexBot, nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnVertex
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, xVertex, yVertex
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness

      real (kind=RKIND), dimension(:,:,:), pointer :: fac
      real (kind=RKIND), dimension(:,:,:), allocatable :: alpha_z

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)

      allocate(alpha_z(nTracers,nVertLevels,nCellsSolve))

      fac => basis % factor

      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            alpha_z(:,iLevel,iCell) = 1.0_RKIND

            zl = +0.5_RKIND*layerThickness(iLevel,iCell)
            do iTracer = 1, nTracers
               dz = fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell)
               dzVert = dz + 2.0_RKIND*zl*fac(5,iLevel,iCell)*coeffs(5,iTracer,iLevel,iCell)
               dzDiff = dzVert - dz

               if (iLevel .gt. 1) then
                  gradMax = max(fac(4,iLevel-1,iCell)*coeffs(4,iTracer,iLevel-1,iCell),fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell))
                  gradMin = min(fac(4,iLevel-1,iCell)*coeffs(4,iTracer,iLevel-1,iCell),fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell))
               else
                  gradMax = fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell)
                  gradMin = fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell)
               end if

               if (abs(dzDiff) .lt. 1.0E-10_RKIND) then
                  alpha_j = 1.0_RKIND
               else if (dzDiff .gt. 0.0_RKIND) then
                  alpha_i = (gradMax - dz)/dzDiff
                  alpha_j = min(1.0_RKIND,alpha_i)
               else if (dzDiff .lt. 0.0_RKIND) then
                  alpha_i = (gradMin - dz)/dzDiff
                  alpha_j = min(1.0_RKIND,alpha_i)
               end if
               alpha_z(iTracer,iLevel,iCell) = min(alpha_z(iTracer,iLevel,iCell),alpha_j)
            end do

            zl = -0.5_RKIND*layerThickness(iLevel,iCell)
            do iTracer = 1, nTracers
               dz = fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell)
               dzVert = dz + 2.0_RKIND*zl*fac(5,iLevel,iCell)*coeffs(5,iTracer,iLevel,iCell)
               dzDiff = dzVert - dz

               if (iLevel .lt. maxLevelCell(iCell)) then
                  gradMax = max(fac(4,iLevel+1,iCell)*coeffs(4,iTracer,iLevel+1,iCell),fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell))
                  gradMin = min(fac(4,iLevel+1,iCell)*coeffs(4,iTracer,iLevel+1,iCell),fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell))
               else
                  gradMax = fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell)
                  gradMin = fac(4,iLevel,iCell)*coeffs(4,iTracer,iLevel,iCell)
               end if

               if (abs(dzDiff) .lt. 1.0E-10_RKIND) then
                  alpha_j = 1.0_RKIND
               else if (dzDiff .gt. 0.0_RKIND) then
                  alpha_i = (gradMax - dz)/dzDiff
                  alpha_j = min(1.0_RKIND,alpha_i)
               else if (dzDiff .lt. 0.0_RKIND) then
                  alpha_i = (gradMin - dz)/dzDiff
                  alpha_j = min(1.0_RKIND,alpha_i)
               end if
               alpha_z(iTracer,iLevel,iCell) = min(alpha_z(iTracer,iLevel,iCell),alpha_j)
            end do
         end do ! iLevel
      end do ! iCell

      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            do iTracer = 1, nTracers
               coeffs(5,iTracer,iLevel,iCell) = alpha_z(iTracer,iLevel,iCell)*coeffs(5,iTracer,iLevel,iCell)
            end do
         end do
      end do

      deallocate(alpha_z)

   end subroutine ocn_cdg_limiter_bj3d_vert_grad!}}}

   subroutine ocn_cdg_limiter_bj3d_split_planar(meshPool, statePool, cdgPool, basis, coeffs)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: coeffs

      integer, pointer :: nVertices, nVertLevels, nCellsSolve, nCells, nTracers
      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell, px, py, pz
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnVertex, edgesOnCell, cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, xCell, yCell, xEdge, yEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: mean, fac

      integer :: iVertex, jVertex, iEdge, jEdge, iCell, jCell, kCell, iLevel, iBasis, iTracer
      real (kind=RKIND) :: phiMinZ, phiMaxZ, phiDiff, phiZ, alpha_i, xl, yl, zl, alpha_x, alpha_z
      real (kind=RKIND), dimension(:,:,:), allocatable :: phiMin, phiMax

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)

      allocate(phiMax(nTracers,nVertLevels+1,nVertices), phiMin(nTracers,nVertLevels+1,nVertices))

      px => basis % power_x
      py => basis % power_y
      pz => basis % power_z
      fac => basis % factor
      mean => basis % mean_new

      ! Compute that max and min mean element values for the elements incident on each vertex
      phiMax(:,:,:) = -1.0E+9_RKIND
      phiMin(:,:,:) = +1.0E+9_RKIND
      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            do iVertex = 1, nEdgesOnCell(iCell)
               jVertex = verticesOnCell(iVertex,iCell)
               do iTracer = 1, nTracers
                  if (coeffs(1,iTracer,iLevel,iCell) .gt. phiMax(iTracer,iLevel,  jVertex)) phiMax(iTracer,iLevel,  jVertex) = coeffs(1,iTracer,iLevel,iCell)
                  if (coeffs(1,iTracer,iLevel,iCell) .lt. phiMin(iTracer,iLevel,  jVertex)) phiMin(iTracer,iLevel,  jVertex) = coeffs(1,iTracer,iLevel,iCell)
                  if (coeffs(1,iTracer,iLevel,iCell) .gt. phiMax(iTracer,iLevel+1,jVertex)) phiMax(iTracer,iLevel+1,jVertex) = coeffs(1,iTracer,iLevel,iCell)
                  if (coeffs(1,iTracer,iLevel,iCell) .lt. phiMin(iTracer,iLevel+1,jVertex)) phiMin(iTracer,iLevel+1,jVertex) = coeffs(1,iTracer,iLevel,iCell)
               end do
            end do
         end do
      end do

      ! Limit the vertical first
      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            do iTracer = 1, nTracers
               alpha_z = 1.0_RKIND

               if (iLevel .gt. 1 .and. iLevel .lt. maxLevelCell(iCell)) then
                  phiMaxZ = max(coeffs(1,iTracer,iLevel+1,iCell), coeffs(1,iTracer,iLevel,iCell), coeffs(1,iTracer,iLevel-1,iCell))
                  phiMinZ = min(coeffs(1,iTracer,iLevel+1,iCell), coeffs(1,iTracer,iLevel,iCell), coeffs(1,iTracer,iLevel-1,iCell))
               else if (iLevel .eq. 1) then
                  phiMaxZ = max(coeffs(1,iTracer,iLevel+1,iCell), coeffs(1,iTracer,iLevel,iCell))
                  phiMinZ = min(coeffs(1,iTracer,iLevel+1,iCell), coeffs(1,iTracer,iLevel,iCell))
               else
                  phiMaxZ = max(coeffs(1,iTracer,iLevel,iCell), coeffs(1,iTracer,iLevel-1,iCell))
                  phiMinZ = min(coeffs(1,iTracer,iLevel,iCell), coeffs(1,iTracer,iLevel-1,iCell))
               end if

               if (iLevel .gt. 1) then
                  zl = +0.5_RKIND*layerThickness(iLevel,iCell)
                  phiDiff = 0.0_RKIND
                  do iBasis = 2, basis % nf
                     if (pz(iBasis) .gt. 0) then
                        phiDiff = phiDiff + ocn_cdg_basis_eval_3d_i(basis, 0.0_RKIND, 0.0_RKIND, zl, fac(iBasis,iLevel,iCell), mean(iBasis,iLevel,iCell), iBasis)
                     end if
                  end do

                  if (abs(phiDiff) .lt. 1.0E-10_RKIND) then
                     alpha_i = 1.0_RKIND
                  else if (phiDiff .gt. 0.0_RKIND) then
                     alpha_i = min(1.0_RKIND, (phiMaxZ - coeffs(1,iTracer,iLevel,iCell))/phiDiff)
                  else if (phiDiff .lt. 0.0_RKIND) then
                     alpha_i = min(1.0_RKIND, (phiMinZ - coeffs(1,iTracer,iLevel,iCell))/phiDiff)
                  end if
                  alpha_z = min(alpha_z,alpha_i)
               end if

               if (iLevel .lt. maxLevelCell(iCell)) then
                  zl = -0.5_RKIND*layerThickness(iLevel,iCell)
                  phiDiff = 0.0_RKIND
                  do iBasis = 2, basis % nf
                     if (pz(iBasis) .gt. 0) then
                        phiDiff = phiDiff + ocn_cdg_basis_eval_3d_i(basis, 0.0_RKIND, 0.0_RKIND, zl, fac(iBasis,iLevel,iCell), mean(iBasis,iLevel,iCell), iBasis)
                     end if
                  end do

                  if (abs(phiDiff) .lt. 1.0E-10_RKIND) then
                     alpha_i = 1.0_RKIND
                  else if (phiDiff .gt. 0.0_RKIND) then
                     alpha_i = min(1.0_RKIND, (phiMaxZ - coeffs(1,iTracer,iLevel,iCell))/phiDiff)
                  else if (phiDiff .lt. 0.0_RKIND) then
                     alpha_i = min(1.0_RKIND, (phiMinZ - coeffs(1,iTracer,iLevel,iCell))/phiDiff)
                  end if
                  alpha_z = min(alpha_z,alpha_i)
               end if

               do iBasis = 1, basis % nf
                  if (basis % power_z(iBasis) .gt. 0) then
                     coeffs(iBasis,iTracer,iLevel,iCell) = alpha_z*coeffs(iBasis,iTracer,iLevel,iCell)
                  end if
               end do ! iBasis
            end do ! iTracer
         end do ! iLevel
      end do ! iCell

      ! Limit in the horizontal second
      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            do iTracer = 1, nTracers
               alpha_x = 1.0_RKIND

!              do iVertex = 1, nEdgesOnCell(iCell)
!                 jVertex = verticesOnCell(iVertex,iCell)

!                 xl = xVertex(jVertex) - xCell(iCell)
!                 yl = yVertex(jVertex) - yCell(iCell)

!                 phiDiff = 0.0_RKIND
!                 do iBasis = 2, basis % nf
!                    if (px(iBasis) + py(iBasis) .gt. 0) then
!                       phiDiff = phiDiff + ocn_cdg_basis_eval_3d_i(basis, xl, yl, 0.0_RKIND, fac(iBasis,iLevel,iCell), mean(iBasis,iLevel,iCell), iBasis)
!                    end if
!                 end do

!                 ! Top of cell
!                 zl = +0.5_RKIND*layerThickness(iLevel,iCell)
!                 phiZ = coeffs(1,iTracer,iLevel,iCell)
!                 do iBasis = 2, basis % nf
!                    if (pz(iBasis) .gt. 0) then
!                       phiZ = phiZ + ocn_cdg_basis_eval_3d_i(basis, 0.0_RKIND, 0.0_RKIND, zl, fac(iBasis,iLevel,iCell), mean(iBasis,iLevel,iCell), iBasis)
!                    end if
!                 end do

!                 if (abs(phiDiff) .lt. 1.0E-10_RKIND) then
!                    alpha_i = 1.0_RKIND
!                 else if (phiDiff .gt. 0.0_RKIND) then
!                    alpha_i = min(1.0_RKIND, (phiMax(iTracer,iLevel,jVertex) - phiZ)/phiDiff)
!                 else if (phiDiff .lt. 0.0_RKIND) then
!                    alpha_i = min(1.0_RKIND, (phiMin(iTracer,iLevel,jVertex) - phiZ)/phiDiff)
!                 end if
!                 alpha_x = min(alpha_x,alpha_i)

!                 ! Bottom of cell
!                 zl = -0.5_RKIND*layerThickness(iLevel,iCell)
!                 phiZ = coeffs(1,iTracer,iLevel,iCell)
!                 do iBasis = 2, basis % nf
!                    if (pz(iBasis) .gt. 0) then
!                       phiZ = phiZ + ocn_cdg_basis_eval_3d_i(basis, 0.0_RKIND, 0.0_RKIND, zl, fac(iBasis,iLevel,iCell), mean(iBasis,iLevel,iCell), iBasis)
!                    end if
!                 end do

!                 if (abs(phiDiff) .lt. 1.0E-10_RKIND) then
!                    alpha_i = 1.0_RKIND
!                 else if (phiDiff .gt. 0.0_RKIND) then
!                    alpha_i = min(1.0_RKIND, (phiMax(iTracer,iLevel+1,jVertex) - phiZ)/phiDiff)
!                 else if (phiDiff .lt. 0.0_RKIND) then
!                    alpha_i = min(1.0_RKIND, (phiMin(iTracer,iLevel+1,jVertex) - phiZ)/phiDiff)
!                 end if
!                 alpha_x = min(alpha_x,alpha_i)
!              end do ! iVertex

               do iEdge = 1, nEdgesOnCell(iCell)
                  jEdge = edgesOnCell(iEdge,iCell)

                  xl = xEdge(jEdge) - xCell(iCell)
                  yl = yEdge(jEdge) - yCell(iCell)

                  phiZ = coeffs(1,iTracer,iLevel,iCell)
                  do iBasis = 2, basis % nf
                     if (pz(iBasis) .gt. 0) then
                        phiZ = phiZ + ocn_cdg_basis_eval_3d_i(basis, 0.0_RKIND, 0.0_RKIND, 0.0_RKIND, fac(iBasis,iLevel,iCell), mean(iBasis,iLevel,iCell), iBasis)
                     end if
                  end do
                  phiDiff = 0.0_RKIND
                  do iBasis = 2, basis % nf
                     if (px(iBasis) + py(iBasis) .gt. 0) then
                        phiDiff = phiDiff + ocn_cdg_basis_eval_3d_i(basis, xl, yl, 0.0_RKIND, fac(iBasis,iLevel,iCell), mean(iBasis,iLevel,iCell), iBasis)
                     end if
                  end do
                  phiMaxZ = -1.0E+9
                  phiMinZ = +1.0E+9
                  !do jCell = 1, 3
                  !   kCell = cellsOnVertex(jCell,jVertex)
                  !   if (kCell .lt. 1 .or. kCell .gt. nCells) cycle
                  !   if (coeffs(1,iTracer,iLevel,kCell) .lt. phiMinZ) phiMinZ = coeffs(1,iTracer,iLevel,kCell)
                  !   if (coeffs(1,iTracer,iLevel,kCell) .gt. phiMaxZ) phiMaxZ = coeffs(1,iTracer,iLevel,kCell)
                  !end do
                  if (cellsOnEdge(1,jEdge) .lt. 1 .or. cellsOnEdge(1,jEdge) .gt. nCells .or. &
                      cellsOnEdge(2,jEdge) .lt. 1 .or. cellsOnEdge(2,jEdge) .gt. nCells) then
                     phiMaxZ = coeffs(1,iTracer,iLevel,iCell)
                     phiMinZ = coeffs(1,iTracer,iLevel,iCell)
                  else
                     phiMaxZ = max(coeffs(1,iTracer,iLevel,cellsOnEdge(1,jEdge)), coeffs(1,iTracer,iLevel,cellsOnEdge(2,jEdge)))
                     phiMinZ = min(coeffs(1,iTracer,iLevel,cellsOnEdge(1,jEdge)), coeffs(1,iTracer,iLevel,cellsOnEdge(2,jEdge)))
                  end if

                  if (abs(phiDiff) .lt. 1.0E-10_RKIND) then
                     alpha_i = 1.0_RKIND
                  else if (phiDiff .gt. 0.0_RKIND) then
                     alpha_i = min(1.0_RKIND, (phiMaxZ - phiZ)/phiDiff)
                  else if (phiDiff .lt. 0.0_RKIND) then
                     alpha_i = min(1.0_RKIND, (phiMinZ - phiZ)/phiDiff)
                  end if
                  alpha_x = min(alpha_x,alpha_i)
               end do


               do iBasis = 1, basis % nf
                  if (px(iBasis) + py(iBasis) .gt. 0) then
                     coeffs(iBasis,iTracer,iLevel,iCell) = alpha_x*coeffs(iBasis,iTracer,iLevel,iCell)
                  end if
               end do ! iBasis
            end do ! iTracer
         end do ! iLevel
      end do ! iCell

      deallocate(phiMax, phiMin)

   end subroutine ocn_cdg_limiter_bj3d_split_planar!}}}

   subroutine ocn_cdg_limiter_weno(meshPool, statePool, cdgPool, basis, coeffs)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: coeffs

      logical, dimension(100) :: troubled_x, troubled_y, troubled_z
      integer :: iCell, jCell, iLevel, iOmega, jOmega, nOmega, iBasis, iTracer
      real (kind=RKIND) :: gamma_neighbour, gamma_cell, eps, r, cell_vol, beta_l
      real (kind=RKIND), dimension(12) :: omega_bar
      real (kind=RKIND), dimension(12,12) :: polys

      real (kind=RKIND), dimension(:,:), allocatable :: zMid
      real (kind=RKIND), dimension(:,:,:,:), allocatable :: coeffsTemp

      logical, pointer :: onSphere
      integer, pointer :: nCells, nCellsSolve, nVertLevels, nTracers
      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell
      real (kind=RKIND), dimension(:), pointer :: areaCell, bottomDepth, xCell, yCell, zCell
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: transGtoC

      r = 2.0_RKIND
      eps = 1.0E-6_RKIND
      if (basis % nf .lt. 6) then
         gamma_neighbour = 1.0E-2_RKIND ! 2nd order polynomial 
      else
         gamma_neighbour = 1.0E-3_RKIND ! 3rd order polynomial
      end if

      call mpas_pool_get_config(meshPool, 'on_a_sphere', onSphere)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_array(cdgPool, 'transGlobalToCell', transGtoC)

      troubled_z(:) = .false.

      allocate(zMid(nVertLevels,nCells))
      allocate(coeffsTemp(basis % nf,nTracers,nVertLevels,nCells))
      coeffsTemp(1:basis % nf,1:nTracers,1:nVertLevels,1:nCells) = coeffs(1:basis % nf,1:nTracers,1:nVertLevels,1:nCells)
      do iCell = 1, nCells
         zMid(maxLevelCell(iCell),iCell) = 0.5_RKIND*layerThickness(maxLevelCell(iCell),iCell) - bottomDepth(iCell)
         do iLevel = maxLevelCell(iCell)-1, 1, -1
            zMid(iLevel,iCell)  = zMid(iLevel+1,iCell) + 0.5_RKIND*(layerThickness(iLevel+1,iCell) + layerThickness(iLevel,iCell))
         end do
      end do

      ! TODO: assumes second order only for now...
      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            cell_vol = areaCell(iCell)*layerThickness(iLevel,iCell)

            ! Determine if cell is troubled
            call ocn_cdg_limiter_troubled_cell_horiz(iCell, iLevel, onSphere, transGtoC(:,:,iCell), xCell, yCell, zCell, &
               cellsOnCell(:,iCell), nEdgesOnCell(iCell), nTracers, basis % factor(:,iLevel,iCell), coeffs, troubled_x, troubled_y)

            call ocn_cdg_limiter_troubled_cell_vert(iCell, iLevel, zMid(:,iCell), maxLevelCell(iCell), nTracers, 4, &
               basis % factor(:,iLevel,iCell), coeffs, troubled_z)

            do iTracer = 1, nTracers
               if (.not.troubled_x(iTracer) .and. .not.troubled_y(iTracer) .and. .not.troubled_z(iTracer)) cycle

               nOmega = 0
               do iOmega = 1, nEdgesOnCell(iCell)
                  jCell = cellsOnCell(iOmega,iCell)
                  if (jCell .gt. nCells) cycle
                  nOmega = nOmega + 1

                  beta_l = 0.0_RKIND
                  do iBasis = 2, basis % nf
                     beta_l = beta_l + coeffsTemp(iBasis,iTracer,iLevel,jCell)*coeffsTemp(iBasis,iTracer,iLevel,jCell)
                  end do
                  beta_l = cell_vol*cell_vol*beta_l
                  omega_bar(nOmega) = gamma_neighbour/((eps + beta_l)**r)

                  polys(1,nOmega) = coeffsTemp(1,iTracer,iLevel,iCell)
                  polys(2:basis % nf,nOmega) = coeffsTemp(2:basis % nf,iTracer,iLevel,jCell)
               end do

               if (iLevel .gt. 1) then
                  nOmega = nOmega + 1

                  beta_l = 0.0_RKIND
                  do iBasis = 2, basis % nf
                     beta_l = beta_l + coeffsTemp(iBasis,iTracer,iLevel-1,iCell)*coeffsTemp(iBasis,iTracer,iLevel-1,iCell)
                  end do
                  beta_l = cell_vol*cell_vol*beta_l
                  omega_bar(nOmega) = gamma_neighbour/((eps + beta_l)**r)

                  polys(1,nOmega) = coeffsTemp(1,iTracer,iLevel,iCell)
                  polys(2:basis % nf,nOmega) = coeffsTemp(2:basis % nf,iTracer,iLevel-1,iCell)
               end if

               if (iLevel .lt. maxLevelCell(iCell)) then
                  nOmega = nOmega + 1

                  beta_l = 0.0_RKIND
                  do iBasis = 2, basis % nf
                     beta_l = beta_l + coeffsTemp(iBasis,iTracer,iLevel+1,iCell)*coeffsTemp(iBasis,iTracer,iLevel+1,iCell)
                  end do
                  beta_l = cell_vol*cell_vol*beta_l
                  omega_bar(nOmega) = gamma_neighbour/((eps + beta_l)**r)

                  polys(1,nOmega) = coeffsTemp(1,iTracer,iLevel,iCell)
                  polys(2:basis % nf,nOmega) = coeffsTemp(2:basis % nf,iTracer,iLevel+1,iCell)
               end if

               gamma_cell = 1.0_RKIND - nOmega*gamma_neighbour

               ! compute weight for local cell
               nOmega = nOmega + 1
               beta_l = 0.0_RKIND
               do iBasis = 2, basis % nf
                  beta_l = beta_l + coeffsTemp(iBasis,iTracer,iLevel,iCell)*coeffsTemp(iBasis,iTracer,iLevel,iCell)
               end do
               beta_l = cell_vol*cell_vol*beta_l
               omega_bar(nOmega) = gamma_cell/((eps + beta_l)**r)
               polys(1:basis % nf,nOmega) = coeffsTemp(1:basis % nf,iTracer,iLevel,iCell)

               ! normallize the weights
               omega_bar(1:nOmega) = omega_bar(1:nOmega)/sum(omega_bar(1:nOmega))

               ! compute new weno polynomial
               coeffs(:,iTracer,iLevel,iCell) = 0.0_RKIND
               do jOmega = 1, nOmega
                  coeffs(:,iTracer,iLevel,iCell) = coeffs(:,iTracer,iLevel,iCell) + omega_bar(jOmega)*polys(1:basis % nf,jOmega)
               end do
            end do ! iTracer
         end do ! iLevel
      end do ! iCell

      deallocate(zMid)
      deallocate(coeffsTemp)

   end subroutine ocn_cdg_limiter_weno!}}}

   function is_troubled(m,n) result(troubled)!{{{

      integer, intent(in) :: n
      real (kind=RKIND), dimension(:), intent(in) :: m
      logical :: troubled

      integer :: i, j
      real (kind=RKIND) :: a

      troubled = .false.

      j = 1
      a = m(1)
      do i = 2, n
         if (abs(m(i)) .lt. abs(m(j))) then
            j = i
            troubled = .true.
         end if
         if (m(i)*a .lt. 0.0_RKIND) then
            troubled = .true.
         end if
         a = m(i)*a
      end do

   end function is_troubled!}}}

   subroutine ocn_cdg_limiter_bp_filter_planar(iCell, nEdgesOnCell, verticesOnCell, theta)!{{{

      integer, intent(in) :: iCell, nEdgesOnCell
      integer, dimension(:), intent(in) :: verticesOnCell
      real (kind=RKIND), dimension(:), intent(out) :: theta

   end subroutine ocn_cdg_limiter_bp_filter_planar!}}}

! Determine if cell is troubled for all tracers based on horiztonal neighbours
   subroutine ocn_cdg_limiter_troubled_cell_horiz(iCell, iLevel, onSphere, transGtoC, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, &!{{{
      nTracers, factor, coeffs, troubled_x, troubled_y)

      integer, intent(in) :: iCell, iLevel, nTracers, nEdgesOnCell
      integer, dimension(:), intent(in) :: cellsOnCell
      logical, intent(in) :: onSphere
      real (kind=RKIND), dimension(:), intent(in) :: xCell, yCell, zCell, factor
      real (kind=RKIND), dimension(:,:), intent(in) :: transGtoC
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: coeffs
      logical, dimension(:), intent(out) :: troubled_x, troubled_y

      integer :: jCell, kCell, iTracer
      real (kind=RKIND), dimension(2) :: dx
      real (kind=RKIND), dimension(12,100) :: ax, ay

      do jCell = 1, nEdgesOnCell
         kCell = cellsOnCell(jCell)

         if (onSphere) then
            call ocn_cdg_tri_trans_g2l(transGtoC(:,:), xCell(iCell), yCell(iCell), zCell(iCell), xCell(kCell), yCell(kCell), zCell(kCell), dx)
         else
            dx(1) = xCell(kCell) - xCell(iCell)
            dx(2) = yCell(kCell) - yCell(iCell)
         end if

         do iTracer = 1, nTracers
            if (abs(dx(1)) .gt. 1.0E-8_RKIND) then
               ax(jCell+1,iTracer) = (coeffs(1,iTracer,iLevel,kCell) - coeffs(1,iTracer,iLevel,iCell))/dx(1)
            else
               ax(jCell+1,iTracer) = 1000.0_RKIND*factor(2)*coeffs(2,iTracer,iLevel,iCell)
            end if
            if (abs(dx(2)) .gt. 1.0E-8_RKIND) then
               ay(jCell+1,iTracer) = (coeffs(1,iTracer,iLevel,kCell) - coeffs(1,iTracer,iLevel,iCell))/dx(2)
            else
               ay(jCell+1,iTracer) = 1000.0_RKIND*factor(3)*coeffs(3,iTracer,iLevel,iCell)
            end if
         end do
      end do

      do iTracer = 1, nTracers
         ax(1,iTracer) = factor(2)*coeffs(2,iTracer,iLevel,iCell)
         ay(1,iTracer) = factor(3)*coeffs(3,iTracer,iLevel,iCell)
         troubled_x(iTracer) = is_troubled(ax(:,iTracer),nEdgesOnCell+1)
         troubled_y(iTracer) = is_troubled(ay(:,iTracer),nEdgesOnCell+1)
      end do

   end subroutine ocn_cdg_limiter_troubled_cell_horiz!}}}

   subroutine ocn_cdg_limiter_troubled_cell_vert(iCell, iLevel, zMid, nVertLevels, nTracers, zInd, factor, coeffs, troubled_z)!{{{

      integer, intent(in) :: iCell, iLevel, nVertLevels, nTracers, zInd
      real (kind=RKIND), dimension(:), intent(in) :: factor, zMid
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: coeffs
      logical, dimension(:), intent(out) :: troubled_z

      integer :: iTracer, iGrad
      real (kind=RKIND) :: dzInv
      real (kind=RKIND), dimension(3,100) :: az

      iGrad = 1
      if (iLevel .gt. 1) then
         iGrad = iGrad + 1
         dzInv = 1.0_RKIND/(zMid(iLevel-1) - zMid(iLevel))
         do iTracer = 1, nTracers
            az(iGrad,iTracer) = (coeffs(1,iTracer,iLevel-1,iCell) - coeffs(1,iTracer,iLevel,iCell))*dzInv
         end do
      end if
      if (iLevel .lt. nVertLevels) then
         iGrad = iGrad + 1
         dzInv = 1.0_RKIND/(zMid(iLevel+1) - zMid(iLevel))
         do iTracer = 1, nTracers
            az(iGrad,iTracer) = (coeffs(1,iTracer,iLevel+1,iCell) - coeffs(1,iTracer,iLevel,iCell))*dzInv
         end do
      end if

      do iTracer = 1, nTracers
         az(1,iTracer) = factor(zInd)*coeffs(zInd,iTracer,iLevel,iCell)
         troubled_z(iTracer) = is_troubled(az(:,iTracer),iGrad)
      end do

   end subroutine ocn_cdg_limiter_troubled_cell_vert!}}}

! Determine the slope limiter coefficients at a single vertex from a rotate
! coordinate system
!
!   [J] = [ cos(theta) -sin(theta) ]
!         [ sin(theta)  cos(theta) ], counter clockwise rotation into [a,b]
!         coordinate system
!
!   [ a ] = [J].[ x ]
!   [ b ]       [ y ]
!
!   u(x,y) = U + u_{x}x + u_{y}y
!          = U + [ u_{x}, u_{y} ].[J^{-1}].[ alpha_a    0    ].[J].[ x ]
!                                          [    0       1    ]     [ y ]
!          = U + u_{a}a + u_{b}b
!
!          => [ u_{x}, u_{y} ] = [ u_{a}, u_{b} ].[J]
!
! reference: Aizinger, Kosik, Kuzmin and Reuter (2017) IJNMF 84, 543--565

   subroutine limit_anisotropic_horiz(xi, yi, coef, mean, min_i, max_i, alpha_x, alpha_y)!{{{

      real (kind=RKIND), intent(in) :: xi, yi, min_i, max_i
      real (kind=RKIND), dimension(:), intent(in) :: coef, mean
      real (kind=RKIND), intent(out) :: alpha_x, alpha_y

      real (kind=RKIND) :: theta, st, ct, ai, bi, ca, alpha_a, alpha_b, detInv
      real (kind=RKIND), dimension(3) :: coef_a
      real (kind=RKIND), dimension(2,2) :: J, Jinv

      theta = atan2(yi,xi)
      st = sin(theta)
      ct = cos(theta)
    
      J(1,1) = +ct
      J(1,2) = -st
      J(2,1) = +st
      J(2,2) = +ct

      detInv = 1.0_RKIND/(J(1,1)*J(2,2) - J(1,2)*J(2,1))

      Jinv(1,1) = +ct*detInv
      Jinv(1,2) = +st*detInv
      Jinv(2,1) = -st*detInv
      Jinv(2,2) = +ct*detInv

      ! rotated coordinates
      ai = J(1,1)*xi + J(1,2)*yi
      bi = J(2,1)*xi + J(2,2)*yi

      ! rotated coefficients
      coef_a(1) = coef(1)
      coef_a(2) = coef(2)*Jinv(1,1) + coef(3)*Jinv(2,1)
      coef_a(3) = coef(2)*Jinv(1,2) + coef(3)*Jinv(2,2)

      ca = coef_a(2)*ai
      if (abs(ca) .lt. 1.0E-8_RKIND) then
         alpha_a = 1.0_RKIND
      else if (ca .gt. 0.0_RKIND) then
         alpha_a = min(1.0_RKIND, (max_i - coef(1))/ca)
      else if (ca .lt. 0.0_RKIND) then
         alpha_a = min(1.0_RKIND, (min_i - coef(1))/ca)
      end if
      alpha_b = 1.0_RKIND

      ! original coefficients
      alpha_x = coef_a(2)*J(1,1) + coef_a(3)*J(2,1)
      alpha_y = coef_a(2)*J(1,2) + coef_a(3)*J(2,2)

   end subroutine limit_anisotropic_horiz!}}}

   subroutine mean_extrema_on_vertex(nCellsOnVert, cellsOnVert, nCells, iTracer, iLevel, coeffs, cMin, cMax)!{{{

      integer, intent(in) :: nCellsOnVert, nCells, iTracer, iLevel
      integer, dimension(:), intent(in) :: cellsOnVert
      real(kind=RKIND), dimension(:,:,:,:), intent(in) :: coeffs
      real(kind=RKIND), intent(out) :: cMin, cMax

      integer :: iCell, jCell

      cMin = +1.0E+9_RKIND
      cMax = -1.0E+9_RKIND

      do iCell = 1, nCellsOnVert
         jCell = cellsOnVert(iCell)

         if(jCell .lt. 1 .or. jCell .gt. nCells) cycle

         cMin = min(cMin, coeffs(1,iTracer,iLevel,jCell))
         cMax = max(cMax, coeffs(1,iTracer,iLevel,jCell))
      end do

   end subroutine mean_extrema_on_vertex!}}}

!***********************************************************************

end module ocn_cdg_limiter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
