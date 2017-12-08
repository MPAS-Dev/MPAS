! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the
! LICENSE file
! distributed with this code, or at
! http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_cdg_advection_vert
!
!> \brief MPAS ocean characteristic discontinuous Galerkin tracer
!>        advection assembly routines
!> \author Dave Lee
!> \date   June 15
!> \details
!>  This module contains the routines to construct swept region fluxes
!>  and cell pre-volumes for the CDG scheme for fields of given 
!>  quadrature and basis orders
!
!-----------------------------------------------------------------------

#define CDG_EDGE_SCALE
#define CDG_QUAD_PT_3D_ADV

module ocn_cdg_advection_vert

   use mpas_kind_types
   use mpas_sort
   use mpas_hash
   use mpas_log

   use ocn_constants
   
   use mpas_matrix_operations  
   use mpas_spline_interpolation

   use mpas_geometry_utils

   use ocn_cdg_types
   use ocn_cdg_basis
   use ocn_cdg_characteristic

   implicit none
   private
   save

   public :: ocn_cdg_advection_vert_init,              &
             ocn_cdg_advection_vert_fluxes,            &
             ocn_cdg_advection_vert_swept_region

   contains

   subroutine ocn_cdg_advection_vert_init(domain, quad_z, quad_w)!{{{

      type (domain_type), intent (inout) :: domain !< Input/Output
      real (kind=RKIND), dimension(:), intent(out) :: quad_z, quad_w

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: cdgPool

      integer, pointer :: num_quad_pts

      ! Assign the quadrature points and weights
      call mpas_pool_get_config(domain % configs, 'config_cdg_quad_order_vert', num_quad_pts)

      ! Gauss Legendre quadrature points between 0 and 1 with weights summing to 1
      select case(num_quad_pts)
         case(1)
            quad_w = (/ 1.0_RKIND /)
            quad_z = (/ 0.5_RKIND /)
         case(2)
            quad_w = (/ 0.5_RKIND, 0.5_RKIND /)
            quad_z = (/ 0.5_RKIND*(1.0_RKIND - sqrt(1.0_RKIND/3.0_RKIND)), 0.5_RKIND*(1.0_RKIND + sqrt(1.0_RKIND/3.0_RKIND)) /)
         case(3)
            quad_w = (/ 5.0_RKIND/18.0_RKIND, 8.0_RKIND/18.0_RKIND, 5.0_RKIND/18.0_RKIND /)
            quad_z = (/ 0.5_RKIND*(1.0_RKIND - sqrt(0.6_RKIND)), 0.5_RKIND, 0.5_RKIND*(1.0_RKIND + sqrt(0.6_RKIND)) /)
         case(4)
            quad_w = (/ 18.0_RKIND - sqrt(30.0_RKIND), 18.0_RKIND + sqrt(30.0_RKIND), 18.0_RKIND + sqrt(30.0_RKIND), 18.0_RKIND - sqrt(30.0_RKIND) /)/72.0_RKIND
            quad_z = (/ 0.5_RKIND*(1.0_RKIND - sqrt(3.0_RKIND/7.0_RKIND + 2.0_RKIND/7.0_RKIND*sqrt(1.2_RKIND))), &
               0.5_RKIND*(1.0_RKIND - sqrt(3.0_RKIND/7.0_RKIND - 2.0_RKIND/7.0_RKIND*sqrt(1.2_RKIND))), &
               0.5_RKIND*(1.0_RKIND + sqrt(3.0_RKIND/7.0_RKIND - 2.0_RKIND/7.0_RKIND*sqrt(1.2_RKIND))), &
               0.5_RKIND*(1.0_RKIND + sqrt(3.0_RKIND/7.0_RKIND + 2.0_RKIND/7.0_RKIND*sqrt(1.2_RKIND))) /)
         case(6)
            quad_w = (/ 0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704 /)/2.0_RKIND
            quad_z = (/ 1.0_RKIND - 0.9324695142031521, 1.0_RKIND - 0.6612093864662645, 1.0_RKIND - 0.2386191860831969, &
               1.0_RKIND + 0.2386191860831969, 1.0_RKIND + 0.6612093864662645, 1.0_RKIND + 0.9324695142031521 /)/2.0_RKIND
         case default
            write(0,*) 'cdg: ERROR! invalid number of quadrature points for the vertical thickness field, must be 1, 2, 3, 4 or 6, not $d', num_quad_pts
            stop
      end select

      call mpas_log_write( 'cdg: initializing vertical quadrature points' )
      call mpas_log_write( '     weights: $r', realArgs=quad_w )
      call mpas_log_write( '     points:  $r', realArgs=quad_z )

   end subroutine ocn_cdg_advection_vert_init!}}}

   subroutine ocn_cdg_advection_vert_fluxes(meshPool, iCell, nVertLevels, nEdges, nTracers, zMid, wTop, depPts, &!{{{
      quad, basis, triArea, dt, fluxes, coeffs, uVertexOld, vVertexOld, wVertexOld, uVertexMid, vVertexMid, wVertexMid, &
      nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, &
      verticesOnCell, edgesOnCell, cellsOnCell, cellsOnVertex, nEdgesOnCell, areaBArray, mean, contFluxCell, &
      zFullOld, zFullMid, wTopNew, wTopOld, wTopMid, onSphere)

      type(mpas_pool_type), intent(in) :: meshPool
      logical, intent(in) :: onSphere
      integer, intent(in) :: iCell, nVertLevels, nEdges, nTracers
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:), intent(in) :: depPts, zMid, wTop, triArea, zFullOld, zFullMid, wTopNew, wTopOld, wTopMid
      real (kind=RKIND), dimension(:,:,:), intent(in) :: coeffs
      type (cdgQuad3D), intent(in) :: quad
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:), intent(out) :: fluxes
      integer, intent(in) :: nCells
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell, edgesOnCell, cellsOnCell, cellsOnVertex
      real (kind=RKIND), dimension(:), intent(in) :: xCell, yCell, zCell, xVertex, yVertex, zVertex, contFluxCell
      real (kind=RKIND), dimension(:,:), intent(in) :: uVertexOld, vVertexOld, wVertexOld, uVertexMid, vVertexMid, wVertexMid, areaBArray, mean

      logical :: done
      integer :: iLevel, iEdge, jEdge, jLevel, iTracer, iBasis, iQuad, jQuad, norm, lInto, lFrom
      real (kind=RKIND) :: dz, pt, wt, bf, bi, zi, pti, ptf
      real (kind=RKIND), dimension(6) :: qzi, qzf, qzh, qw, wi, wh, tracer
      real (kind=RKIND), dimension(3,12) :: qFinal
      real (kind=RKIND), dimension(:), allocatable :: fluxScale

      allocate(fluxScale(nVertLevels))
      !call ocn_cdg_advection_vert_flux_scaling(nVertLevels, zFullNew, wTopNew, depPts, contFluxCell, fluxScale)
      call ocn_cdg_advection_vert_flux_scaling(nVertLevels, zFullOld, wTopNew, depPts, contFluxCell, fluxScale)

      ! Assume z is -ve with depth (and layer index)
      do iLevel = 2, nVertLevels
         !if (wTopNew(iLevel) .gt. 0.0_RKIND) then
         if (zFullOld(iLevel) - depPts(iLevel) .gt. 0.0_RKIND) then
            norm = +1           ! Direction of the vertical velocity (oppsite direction in which swept regions are integrated)
            lInto = iLevel - 1  ! Cell flux is going into
            lFrom = iLevel      ! Cell flux is going from
            jLevel = iLevel     ! Cell to integrate flux from
         else
            norm = -1
            lInto = iLevel
            lFrom = iLevel - 1
            jLevel = iLevel - 1
         end if

         done = .false.
         do while (.not. done)
            ! Departure point is located in jLevel
            if (norm .eq. +1 .and. jLevel .eq. nVertLevels) then
               done = .true.
            else if (norm .eq. -1 .and. jLevel .eq. 1) then
               done = .true.
            else if (norm .eq. +1 .and. zFullOld(jLevel+1) .lt. depPts(iLevel)) then
               done = .true.
            else if (norm .eq. -1 .and. zFullOld(jLevel+0) .gt. depPts(iLevel)) then
               done = .true.
            end if

            if (done .and. norm .eq. +1) then
               dz = zFullOld(jLevel) - depPts(iLevel)
               zi = depPts(iLevel)
            else if (done .and. norm .eq. -1) then
               dz = depPts(iLevel) - zFullOld(jLevel+1)
               zi = zFullOld(jLevel+1)
            else if (.not. done) then
               dz = zFullOld(jLevel) - zFullOld(jLevel+1)
               zi = zFullOld(jLevel+1)
            end if

            ! Initialize and integrate the quadrature points (in ascending order)
            qzi(1:quad % nv) = zi + dz*quad % qzv(1:quad % nv)
            call ocn_cdg_char_integrate_quad_pts_vertical(quad, nVertLevels, jLevel, zFullOld, zFullMid, wTopOld, wTopMid, dt, zi, dz, qzf)

            do iEdge = 1, nEdges
               jEdge = edgesOnCell(iEdge,iCell)
#ifdef CDG_QUAD_PT_3D_ADV
               if (onSphere) then
                  call ocn_cdg_char_integrate(meshPool, quad % nh, quad % xQuadG(:,iEdge,iCell), quad % yQuadG(:,iEdge,iCell), quad % zQuadG(:,iEdge,iCell), &
                     uVertexOld, vVertexOld, wVertexOld, uVertexMid, vVertexMid, wVertexMid, nCells, xCell, yCell, zCell, &
                     xVertex, yVertex, zVertex, verticesOnCell, cellsOnCell, nEdgesOnCell, iLevel, iCell, areaBArray, dt, qFinal)
               else
                  call ocn_cdg_char_integrate_planar(meshPool, quad % nh, quad % xQuadL(:,iEdge,iCell) + xCell(iCell), quad % yQuadL(:,iEdge,iCell) + yCell(iCell), &
                     uVertexOld, vVertexOld, uVertexMid, vVertexMid, nCells, xCell, yCell, zCell, xVertex, yVertex, &
                     verticesOnCell, cellsOnCell, cellsOnVertex, nEdgesOnCell, iLevel, iCell, areaBArray, dt, qFinal)
               end if
#endif

               do iQuad = 1, quad % nv
                  do jQuad = 1, quad % nh
                     pt = qzi(iQuad) - zMid(jLevel)
                     ! Evaluate the tracer values at the initial quadrature points
                     do iTracer = 1, nTracers
                        tracer(iTracer) = ocn_cdg_basis_eval_3d(basis, quad % xQuadL(jQuad,iEdge,iCell), quad % yQuadL(jQuad,iEdge,iCell), pt, &
                           basis % factor(:,jLevel,iCell), mean(:,jLevel), coeffs(:,iTracer,jLevel))
                     end do

#ifdef CDG_EDGE_SCALE
                     ! Rescaling by volume flux per unit area
                     wt = quad % qwv(iQuad)*triArea(iEdge)*quad % qw(jQuad)*dz*fluxScale(iLevel)
#else
                     wt = quad % qwv(iQuad)*triArea(iEdge)*quad % qw(jQuad)*dz
#endif

                     ! Evaluate the basis values
                     do iBasis = 1, basis % nf
                        pti = qzf(iQuad) - zMid(lInto)
                        ptf = qzf(iQuad) - zMid(lFrom)
#ifdef CDG_QUAD_PT_3D_ADV
                        bi = ocn_cdg_basis_eval_3d_i(basis, qFinal(1,jQuad) - xCell(iCell), qFinal(2,jQuad) - yCell(iCell), pti, &
                                basis % factor(iBasis,lInto,iCell), mean(iBasis,lInto), iBasis)
                        bf = ocn_cdg_basis_eval_3d_i(basis, qFinal(1,jQuad) - xCell(iCell), qFinal(2,jQuad) - yCell(iCell), ptf, &
                                basis % factor(iBasis,lFrom,iCell), mean(iBasis,lFrom), iBasis)
#else
                        bi = ocn_cdg_basis_eval_3d_i(basis, quad % xQuadL(jQuad,iEdge,iCell), quad % yQuadL(jQuad,iEdge,iCell), pti, &
                                basis % factor(iBasis,lInto,iCell), mean(iBasis,lInto), iBasis)
                        bf = ocn_cdg_basis_eval_3d_i(basis, quad % xQuadL(jQuad,iEdge,iCell), quad % yQuadL(jQuad,iEdge,iCell), ptf, &
                                basis % factor(iBasis,lFrom,iCell), mean(iBasis,lFrom), iBasis)
#endif
                        do iTracer = 1, nTracers
                           fluxes(iBasis,iTracer,lFrom) = fluxes(iBasis,iTracer,lFrom) - wt*bf*tracer(iTracer)
                           fluxes(iBasis,iTracer,lInto) = fluxes(iBasis,iTracer,lInto) + wt*bi*tracer(iTracer)
                        end do ! iTracer
                     end do ! iBasis
                  end do ! jQuad
               end do ! iQuad
            end do ! iEdge

            jLevel = jLevel + norm
         end do
      end do

      deallocate(fluxScale)

   end subroutine ocn_cdg_advection_vert_fluxes!}}}

   subroutine ocn_cdg_advection_vert_swept_region(meshPool, particlesPool, cdgPool, quad, basis, dt, flux, coeffs)!{{{

      type(mpas_pool_type), intent(in) :: meshPool
      type(mpas_pool_type), intent(in) :: particlesPool
      type(mpas_pool_type), intent(in) :: cdgPool
      type(cdgQuad3d), intent(in) :: quad
      type(cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: flux
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: coeffs

      integer :: iCell, iLevel

      logical, pointer :: onSphere
      integer, pointer :: nCells, nVertLevels, nTracers
      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, edgesOnCell, cellsOnVertex
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, bottomDepth
      real (kind=RKIND), dimension(:,:), pointer :: tend_layerThickness, triAreas
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessOld, layerThicknessMid, wTopOld, wTopNew, wTopMid, &
                                                    uVertexOld, uVertexMid, wVertexMid, vVertexOld, vVertexMid, wVertexOld

      integer, dimension(:,:), pointer :: verticesOnCell
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex, xCell, yCell, zCell
      real (kind=RKIND), dimension(:,:), pointer :: areaBArray, contFluxCell
      real (kind=RKIND), dimension(:,:,:), pointer :: mean

      real (kind=RKIND), dimension(:), allocatable :: zHalf, wHalf, depPts, zMidOld, zFullOld, zFullMid

      call mpas_pool_get_config(meshPool, 'on_a_sphere', onSphere)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(cdgPool, 'layerThicknessOld', layerThicknessOld)
      call mpas_pool_get_array(cdgPool, 'layerThicknessMid', layerThicknessMid)
      call mpas_pool_get_array(cdgPool, 'vertAleVelocitySum', wTopNew)
      call mpas_pool_get_array(cdgPool, 'vertAleVelocityOld', wTopOld)
      call mpas_pool_get_array(cdgPool, 'vertAleVelocityMid', wTopMid)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityOld', uVertexOld)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityOld', vVertexOld)
      call mpas_pool_get_array(cdgPool, 'wVertexVelocityOld', wVertexOld)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityMid', uVertexMid)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityMid', vVertexMid)
      call mpas_pool_get_array(cdgPool, 'wVertexVelocityMid', wVertexMid)
      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)
      call mpas_pool_get_array(particlesPool, 'wachspressAreaB', areaBArray)
      call mpas_pool_get_array(cdgPool, 'contFluxCell', contFluxCell)

      mean => basis % mean_cur

      allocate(zHalf(nVertLevels), wHalf(nVertLevels), depPts(nVertLevels), zMidOld(nVertLevels))
      allocate(zFullOld(nVertLevels+1), zFullMid(nVertLevels+1))

      do iCell = 1, nCells
         zFullOld(maxLevelCell(iCell)+1) = -bottomDepth(iCell)
         zFullMid(maxLevelCell(iCell)+1) = -bottomDepth(iCell)
         do iLevel = maxLevelCell(iCell), 1, -1
            zFullOld(iLevel) = zFullOld(iLevel+1) + layerThicknessOld(iLevel,iCell)
            zFullMid(iLevel) = zFullMid(iLevel+1) + layerThicknessMid(iLevel,iCell)
            zMidOld(iLevel)  = zFullOld(iLevel+1) + 0.5_RKIND*layerThicknessOld(iLevel,iCell)
         end do

         call ocn_cdg_char_integrate_dep_pts_vertical(maxLevelCell(iCell), layerThicknessOld(:,iCell), layerThicknessMid(:,iCell), &
            wTopNew(:,iCell), wTopMid(:,iCell), bottomDepth(iCell), dt, depPts)

         ! Integrate the fluxes for the tracers
         call ocn_cdg_advection_vert_fluxes(meshPool, iCell, maxLevelCell(iCell), nEdgesOnCell(iCell), nTracers, zMidOld, &
            wTopOld(:,iCell), depPts, quad, basis, triAreas(:,iCell), dt, flux(:,:,:,iCell), coeffs(:,:,:,iCell), &
            uVertexOld, vVertexOld, wVertexOld, uVertexMid, vVertexMid, wVertexMid, nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, verticesOnCell, &
            edgesOnCell, cellsOnCell, cellsOnVertex, nEdgesOnCell, areaBArray, mean(:,:,iCell), contFluxCell(:,iCell), &
            zFullOld, zFullMid, wTopNew(:,iCell), wTopOld(:,iCell), wTopMid(:,iCell), onSphere)
      end do

      deallocate(zHalf, wHalf, depPts, zMidOld)
      deallocate(zFullOld, zFullMid)

   end subroutine ocn_cdg_advection_vert_swept_region!}}}

   subroutine ocn_cdg_advection_vert_flux_scaling(nVertLevels, zFull, wTop, depPts, contFluxCell, fluxScale)!{{{

      integer, intent(in) :: nVertLevels
      real (kind=RKIND), dimension(:), intent(in) :: depPts, zFull, wTop, contFluxCell
      real (kind=RKIND), dimension(:), intent(out) :: fluxScale

      integer :: iLevel
      real (kind=RKIND) :: dz

      fluxScale(:) = 0.0_RKIND

      ! Assume z is -ve with depth (and layer index)
      do iLevel = 2, nVertLevels
         if (depPts(iLevel) .gt. zFull(1)) then
            dz = zFull(1) - zFull(iLevel)
         else if (depPts(iLevel) .lt. zFull(nVertLevels+1)) then
            dz = zFull(iLevel) - zFull(nVertLevels+1)
         else
            dz = zFull(iLevel) - depPts(iLevel)
         end if

         if (abs(dz) .lt. 1.0E-12) then
            fluxScale(iLevel) = 0.0_RKIND
         ! continuity and cdg fluxes are of opposite direction across the cell top
         !else if (contFluxCell(iLevel)*wTop(iLevel) .lt. 0.0_RKIND) then
         !   fluxScale(iLevel) = 0.0_RKIND
         else
            fluxScale(iLevel) = abs(contFluxCell(iLevel)/dz)
         end if
      end do

   end subroutine ocn_cdg_advection_vert_flux_scaling!}}}

end module ocn_cdg_advection_vert

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
