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
!  ocn_cdg_manager
!
!> \brief MPAS ocean characteristic discontinuous Galerkin tracer
!>        advection manager
!> \author Dave Lee
!> \date   December 15
!> \details
!>  This module manages the initialization and execution of the CDG
!>  advection routine for the h and tracer field. Specifically
!>  it executes the advection of the h field (h) and h
!>  weighted tracer fields (hq), and then backs out the tracer fields
!>  (q).
!
!-----------------------------------------------------------------------

#define CDG_TESTING
!#define CDG_PRESCRIBED_VELOCITY
!#define CDG_ADVECTION_TEST
!#define CDG_INIT_ANALYTIC

! Lock exchange
#define NX_LOCAL 8
#define XLEN 8000.0_RKIND
#define DELX 500.0_RKIND
! Overflow
!#define NX_LOCAL 8
!#define XLEN 16000.0_RKIND
!#define DELX 1000.0_RKIND
!#define XLEN 128000.0_RKIND
!#define DELX 8000.0_RKIND
! Baroclinic channel
!#define NX_LOCAL 80
!#define XLEN 400000.0_RKIND
!#define DELX 2500.0_RKIND

#define EPS 1.0E-3_RKIND

module ocn_cdg_manager

   use mpas_kind_types
   use mpas_sort
   use mpas_hash
   use mpas_timer
   use mpas_log

   use mpas_pool_routines
   use mpas_rbf_interpolation
   use mpas_vector_reconstruction

   use ocn_constants

   use mpas_matrix_operations
   use mpas_spline_interpolation

   use ocn_particle_helpers

   use ocn_cdg_types   
   use ocn_cdg_triangle
   use ocn_cdg_basis
   use ocn_cdg_characteristic
   use ocn_cdg_limiter
   use ocn_cdg_advection
   use ocn_cdg_advection_vert

   implicit none
   private
   save

   public :: ocn_cdg_manager_init,            &
             ocn_cdg_manager_finalize,        &
             ocn_cdg_manager_advect,          &
             ocn_cdg_manager_test,            &
             ocn_cdg_manager_run,             &
             ocn_cdg_manager_pre_solve,       &
             ocn_cdg_manager_post_solve,      &
             ocn_cdg_manager_init_periodic_geom

   logical :: firstStep
   integer :: step
   real (kind=RKIND) :: hqInitVol
   type (cdgQuadrature) :: hqQuad
   type (cdgQuad3D) :: hqQuad3D
   type (cdgBasis) :: hqBasis
   type (cdgBasis3D) :: hqBasis3D

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_manager_init
!
!> \brief MPAS ocean CDG advection initialization routine manager
!> \author Dave Lee
!> \date   December 2015
!> \details
!>  This routine calls the initialization routines for the advection of
!>  the thickenss (h) and theickness weighted tracer fields (hq).
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_manager_run(domain, dt)!{{{

      type (domain_type), intent (inout) :: domain !< Input/Output
      real (kind=RKIND), intent (in) :: dt

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: cdgPool

      logical, pointer :: config_cdg_do_passive
      integer, pointer :: config_cdg_lim, config_cdg_do_every

      real (kind=RKIND), dimension(:,:,:,:), pointer :: hqCur, hqNew
      type (field4DReal), pointer :: hqCoeffsField

      call mpas_pool_get_config(domain % configs, 'config_cdg_do_passive', config_cdg_do_passive)
      call mpas_pool_get_config(domain % configs, 'config_cdg_limiter_type', config_cdg_lim)
      call mpas_pool_get_config(domain % configs, 'config_cdg_do_every', config_cdg_do_every)

      if (.not. config_cdg_do_passive) return

      ! 1. Update cdg tracer coefficients at previous time level
      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'cdg', cdgPool)
         call mpas_pool_get_array(cdgPool, 'hqCoeffsCur', hqCur)
         call mpas_pool_get_array(cdgPool, 'hqCoeffsNew', hqNew)
         hqCur(:,:,:,:) = hqNew(:,:,:,:)

         block => block % next
      end do

      ! 2. Increment cdg transport velocities and departure points
      call mpas_timer_start("cdg pre solve")
      call ocn_cdg_manager_pre_solve(domain, dt, config_cdg_do_every)
      call mpas_timer_stop("cdg pre solve")

      ! 3. Apply cdg transport routines
      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'cdg', cdgPool)

         call mpas_pool_get_array(cdgPool, 'hqCoeffsCur', hqCur)
         call mpas_pool_get_array(cdgPool, 'hqCoeffsNew', hqNew)

         call mpas_timer_start("cdg advection")
         call ocn_cdg_manager_advect(meshPool, statePool, diagnosticsPool, tendPool, scratchPool, cdgPool, dt, &
            config_cdg_lim, config_cdg_do_every, hqCur, hqNew)
         call mpas_timer_stop("cdg advection")

         block => block % next
      end do

      ! 4. Halo update and cleanup
      call mpas_pool_get_subpool(domain % blocklist % structs, 'cdg', cdgPool)
      call mpas_pool_get_field(cdgPool, 'hqCoeffsNew', hqCoeffsField)
      call mpas_timer_start("cdg halo")
      call mpas_dmpar_exch_halo_field(hqCoeffsField)
      call mpas_timer_stop("cdg halo")

      call mpas_timer_start("cdg post solve")
      call ocn_cdg_manager_post_solve(domain, dt)
      call mpas_timer_stop("cdg post solve")

   end subroutine ocn_cdg_manager_run!}}}

   subroutine ocn_cdg_manager_init(domain, err)!{{{

      type (domain_type), intent (inout) :: domain !< Input/Output
      integer, intent(out) :: err !< Output: the error flag

      type (block_type), pointer :: block
      type (dm_info), pointer :: dminfo

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: cdgPool

      logical, pointer :: onSphere, config_do_restart
      integer, pointer :: nCells, nVertLevels, nTracers, index_temp
      integer, dimension(:), pointer :: nIncCells, maxLevelCell
      integer, dimension(:,:), pointer :: incCells
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, layerThicknessOld, qField, contFluxEdge, contFluxCell, &
                                                    normVelocitySum, vertAleVelocitySum, normVelocityOld, vertAleVelocityOld, &
                                                    uVertexOld, vVertexOld, wVertexOld
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracersCur, activeTracersNew, debugTracersCur, debugTracersNew
      real (kind=RKIND), dimension(:,:,:,:), pointer :: hqCoeffsNew, hqCoeffsCur
      type (field4DReal), pointer :: hqCoeffsField

      integer :: iCell, iLevel, iTracer, jTracer
      real (kind=RKIND) :: hqInitVolLocal

      err = 0

      step = 0
      firstStep = .true.

      !call ocn_cdg_manager_init_periodic_geom( domain )

      hqQuad % fieldName = 'hq'
      hqQuad3D % fieldName = 'hq'
      call ocn_cdg_tri_init(domain, hqQuad, err)
      call ocn_cdg_tri_init_3d(domain, hqQuad3D, err)

      hqBasis % fieldName = 'hq'
      call ocn_cdg_basis_init_3d(domain, hqBasis3D, err)

      ! Initialize the triangle areas and the global to cell transformations
      call ocn_cdg_init_areas(domain, err)

      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)

      block => domain % blockList
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'cdg', cdgPool)

         call mpas_pool_get_config(meshPool, 'on_a_sphere', onSphere)

         ! Find the incident cells for each edge
         call mpas_pool_get_array(cdgPool, 'nIncCells', nIncCells)
         call mpas_pool_get_array(cdgPool, 'incCells', incCells)
         call ocn_cdg_manager_init_inc_cells(meshPool, cdgPool, nIncCells, incCells)

         ! Populate the quadrature and basis fields
         call cdg_quadrature_populate(cdgPool, hqQuad)
         call cdg_basis_populate_3d(cdgPool, hqBasis3D)

         ! Calculate the 2D local and 3D global quadrature points for each triangle of each cell
         call ocn_cdg_tri_init_quad_pts(meshPool, cdgPool, onSphere, hqQuad)
         hqQuad3D % xQuadL => hqQuad % xQuadL
         hqQuad3D % yQuadL => hqQuad % yQuadL
         hqQuad3D % xQuadG => hqQuad % xQuadG
         hqQuad3D % yQuadG => hqQuad % yQuadG
         hqQuad3D % zQuadG => hqQuad % zQuadG

         allocate(hqQuad3D % qwv(hqQuad3D % nv), hqQuad3D % qzv(hqQuad3D % nv))
         call ocn_cdg_advection_vert_init(domain, hqQuad3d % qzv, hqQuad3D % qwv)

         ! Calculate the basis scaling factors
         call ocn_cdg_basis_init_factor_3d(meshPool, statePool, cdgPool, hqBasis3D, err)

         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temp)
         call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersCur, 1)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersNew, 2)
         call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracersCur, 1)
         call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracersNew, 2)
         call mpas_pool_get_array(cdgPool, 'hqCoeffsNew', hqCoeffsNew)
         call mpas_pool_get_array(cdgPool, 'hqCoeffsCur', hqCoeffsCur)

         ! Initialize the h weighted tracer field
         if (.not. config_do_restart) then
            hqCoeffsNew = 0.0_RKIND
            hqCoeffsCur = 0.0_RKIND
#ifdef CDG_INIT_ANALYTIC
            call ocn_cdg_advection_init_c_kj(meshPool, statePool, cdgPool, hqQuad3D, hqBasis3D, hqCoeffsNew)
            hqCoeffsCur(:,:,:,:) = hqCoeffsNew(:,:,:,:)
#endif
 
            debugTracersCur = 0.0_RKIND
            debugTracersNew = 0.0_RKIND
            do iCell = 1, nCells
               do iLevel = 1, nVertLevels
#ifdef CDG_INIT_ANALYTIC
                  activeTracersCur(index_temp,iLevel,iCell) = hqCoeffsCur(1,1,iLevel,iCell)
                  activeTracersNew(index_temp,iLevel,iCell) = hqCoeffsNew(1,1,iLevel,iCell)
#else
                  hqCoeffsNew(1,1,iLevel,iCell) = activeTracersNew(index_temp,iLevel,iCell)
                  hqCoeffsCur(1,1,iLevel,iCell) = activeTracersCur(index_temp,iLevel,iCell)
#endif
                  debugTracersCur(1,iLevel,iCell) = hqCoeffsCur(1,1,iLevel,iCell)
                  debugTracersNew(1,iLevel,iCell) = hqCoeffsNew(1,1,iLevel,iCell)
               end do
            end do
         end if
!for the vertical advection test
#ifdef CDG_ADVECTION_TEST
call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
call ocn_cdg_vert_adv(meshPool, statePool, diagnosticsPool, cdgPool, .true.)
#endif
         ! Integrate the tracer field to test mass conservation at each time step
         call ocn_cdg_advection_integrate(meshPool, statePool, cdgPool, 1, hqQuad3D, hqBasis3D, hqCoeffsNew, hqInitVolLocal)

         call mpas_pool_get_array(cdgPool, 'qField', qField, 1)
         do iCell = 1, nCells
            do iLevel = 1, maxLevelCell(iCell)
               qField(iLevel,iCell) = hqCoeffsNew(1,1,iLevel,iCell)
            end do
         end do

         call mpas_pool_get_array(cdgPool, 'contFluxEdge', contFluxEdge)
         call mpas_pool_get_array(cdgPool, 'contFluxCell', contFluxCell)
         contFluxEdge(:,:) = 0.0_RKIND
         contFluxCell(:,:) = 0.0_RKIND

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(cdgPool, 'layerThicknessOld', layerThicknessOld)
         layerThicknessOld(:,:) = layerThickness(:,:)

         call mpas_pool_get_array(cdgPool, 'normVelocitySum', normVelocitySum)
         call mpas_pool_get_array(cdgPool, 'vertAleVelocitySum', vertAleVelocitySum)
         normVelocitySum(:,:) = 0.0_RKIND
         vertAleVelocitySum(:,:) = 0.0_RKIND
         call mpas_pool_get_array(cdgPool, 'normVelocityOld', normVelocityOld)
         call mpas_pool_get_array(cdgPool, 'vertAleVelocityOld', vertAleVelocityOld)
         normVelocityOld(:,:) = 0.0_RKIND
         vertAleVelocityOld(:,:) = 0.0_RKIND

         call mpas_pool_get_array(cdgPool, 'uVertexVelocityOld', uVertexOld)
         call mpas_pool_get_array(cdgPool, 'vVertexVelocityOld', vVertexOld)
         call mpas_pool_get_array(cdgPool, 'wVertexVelocityOld', wVertexOld)
         uVertexOld(:,:) = 0.0_RKIND
         vVertexOld(:,:) = 0.0_RKIND
         wVertexOld(:,:) = 0.0_RKIND

         block => block % next
      end do

      ! Update the halo regions for the layer h basis coefficients
      call mpas_pool_get_subpool(domain % blocklist % structs, 'cdg', cdgPool)
      call mpas_pool_get_field(cdgPool, 'hqCoeffsNew', hqCoeffsField)
      call mpas_dmpar_exch_halo_field(hqCoeffsField)
      call mpas_pool_get_field(cdgPool, 'hqCoeffsCur', hqCoeffsField)
      call mpas_dmpar_exch_halo_field(hqCoeffsField)

      dminfo => domain % dminfo
      call mpas_dmpar_sum_real(dminfo, hqInitVolLocal, hqInitVol)
      call mpas_log_write( 'cdg: inititial hq-mass $r', realArgs=(/ hqInitVol /) )

      ! Initialize the normal velocity field
      block => domain % blockList
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

#ifdef CDG_PRESCRIBED_VELOCITY
         call ocn_cdg_init_normal_velocity(meshPool, diagnosticsPool, cdgPool, time)
#endif

         ! Reconstruct the vertex based velocity vectors
         call mpas_rbf_interp_initialize(meshPool)
         call mpas_init_reconstruct(meshPool)

         block => block % next
      end do

      ! Initialize the cell triangle barycentric areas for the Wachspress interpolants
      !call initialize_wachspress_coefficients(domain, err)
      call mpas_initialize_wachspress_coefficients(domain, err)

   end subroutine ocn_cdg_manager_init!}}}

   subroutine ocn_cdg_manager_init_inc_cells(meshPool, cdgPool, nIncCells, incCells)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: cdgPool
      integer, dimension(:), intent(out) :: nIncCells
      integer, dimension(:,:), intent(out) :: incCells

      integer, pointer :: nEdges, nCells, vertexDegree, maxIncCells
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnEdge, cellsOnVertex, edgesOnCell, cellsOnCell

      logical :: addCell
      integer :: iEdge, jEdge, iVertex, jVertex, iCell, jCell, kCell, lCell, testCell, nIncCells1

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

      call mpas_pool_get_dimension(cdgPool, 'maxIncCells', maxIncCells)

      do iEdge = 1, nEdges
         iVertex = verticesOnEdge(1,iEdge)
         jVertex = verticesOnEdge(2,iEdge)

         nIncCells(iEdge) = 0
         do iCell = 1, nCells
            addCell = .false.
            do jCell = 1, vertexDegree
               if (iCell .eq. cellsOnVertex(jCell,iVertex)) addCell = .true.
               if (iCell .eq. cellsOnVertex(jCell,jVertex)) addCell = .true.
            end do
            if (addCell) then
               nIncCells(iEdge) = nIncCells(iEdge) + 1
               incCells(nIncCells(iEdge),iEdge) = iCell
            end if
         end do
      end do ! iEdge

      if (maxIncCells .le. 6) return

      ! Construct neighbour cell lists for CFL > 1
      do iEdge = 1, nEdges
         nIncCells1 = nIncCells(iEdge)

         do iCell = 1, nIncCells1
            jCell = incCells(iCell,iEdge)

            if (jCell .lt. 1 .or. jCell .gt. nCells) cycle

            do jEdge = 1, nEdgesOnCell(jCell)
               testCell = cellsOnCell(jEdge,jCell)

               addCell = .true.

               if (testCell .lt. 1) addCell = .false.
               if (testCell .gt. nCells) addCell = .false.

               do kCell = 1, nIncCells(iEdge)
                  if (incCells(kCell,iEdge) .eq. testCell) addCell = .false.
               end do

               if (addCell) then
                  nIncCells(iEdge) = nIncCells(iEdge) + 1
                  incCells(nIncCells(iEdge),iEdge) = testCell
               end if
            end do
         end do
      end do

      if (maxIncCells .le. 20) return 

      ! Construct neighbour cell lists for CFL > 2
      do iEdge = 1, nEdges
         nIncCells1 = nIncCells(iEdge)

         do iCell = 1, nIncCells1
            jCell = incCells(iCell,iEdge)

            if (jCell .lt. 1 .or. jCell .gt. nCells) cycle

            do jEdge = 1, nEdgesOnCell(jCell)
               testCell = cellsOnCell(jEdge,jCell)

               addCell = .true.

               if (testCell .lt. 1) addCell = .false.
               if (testCell .gt. nCells) addCell = .false.

               do kCell = 1, nIncCells(iEdge)
                  if (incCells(kCell,iEdge) .eq. testCell) addCell = .false.
               end do

               if (addCell) then
                  nIncCells(iEdge) = nIncCells(iEdge) + 1
                  incCells(nIncCells(iEdge),iEdge) = testCell
               end if
            end do
         end do
      end do

   end subroutine ocn_cdg_manager_init_inc_cells!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_manager_advect
!
!> \brief MPAS ocean CDG advection
!> \author Dave Lee
!> \date   December 2015
!> \details
!>  This routine executes the advection of the h field (h) and
!>  the h weighter tracer field (hq), and then backs out the
!>  tracer field (q)
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_manager_pre_solve(domain, dt, do_every)!{{{

      type (domain_type), intent(inout) :: domain !< Input/Output
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: do_every

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: cdgPool

      type (block_type), pointer :: block
      logical, pointer :: onSphere
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, contFluxEdge, contFluxCell, normVelocitySum, normVelocityMid, normVelocityOld, vertAleVelocitySum
      real (kind=RKIND), dimension(:,:), pointer :: uVertexNew, vVertexNew, wVertexNew, uVertexOld, vVertexOld, wVertexOld, uVertexMid, vVertexMid, wVertexMid
      real (kind=RKIND), dimension(:,:), pointer :: vertAleVelocityOld, vertAleVelocityMid, layerThicknessOld, layerThicknessMid
      real (kind=RKIND), dimension(:,:,:), pointer :: departurePoints, departureEdges
      type (field2DReal), pointer :: uVertex, vVertex, wVertex

      integer :: iEdge, iLevel

      block => domain % blockList
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'cdg', cdgPool)

         ! Reconstruct the vertex based velocity vectors
#ifdef CDG_PRESCRIBED_VELOCITY
         call ocn_cdg_init_normal_velocity(meshPool, diagnosticsPool, cdgPool, time)
#endif
!for the vertical advection test
#ifdef CDG_ADVECTION_TEST
call ocn_cdg_thick(meshPool, diagnosticsPool, statePool, dt, time)
call ocn_cdg_vert_adv(meshPool, statePool, diagnosticsPool, cdgPool, .false.)
#endif

         call mpas_pool_get_array(cdgPool, 'contFluxEdge', contFluxEdge)
         call mpas_pool_get_array(cdgPool, 'contFluxCell', contFluxCell)
         call mpas_pool_get_array(cdgPool, 'normVelocitySum', normVelocitySum)
         call mpas_pool_get_array(cdgPool, 'vertAleVelocitySum', vertAleVelocitySum)

         call mpas_timer_start("cdg cont flux")
         call ocn_cdg_manager_increment_cont_flux(meshPool, diagnosticsPool, cdgPool, dt, contFluxEdge, contFluxCell, normVelocitySum, vertAleVelocitySum)
         call mpas_timer_stop("cdg cont flux")

         if (mod(step,do_every) .eq. do_every - 1) then
            call mpas_rbf_interp_initialize(meshPool)
            call mpas_init_reconstruct(meshPool)

            call mpas_pool_get_config(meshPool, 'on_a_sphere', onSphere)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
            call mpas_pool_get_array(cdgPool, 'layerThicknessOld', layerThicknessOld)
            call mpas_pool_get_array(cdgPool, 'layerThicknessMid', layerThicknessMid)
            call mpas_pool_get_array(cdgPool, 'vertAleVelocitySum', vertAleVelocitySum)
            call mpas_pool_get_array(cdgPool, 'vertAleVelocityOld', vertAleVelocityOld)
            call mpas_pool_get_array(cdgPool, 'vertAleVelocityMid', vertAleVelocityMid)
            call mpas_pool_get_array(cdgPool, 'normVelocitySum', normVelocitySum)
            call mpas_pool_get_array(cdgPool, 'normVelocityMid', normVelocityMid)
            call mpas_pool_get_array(cdgPool, 'normVelocityOld', normVelocityOld)
            call mpas_pool_get_array(cdgPool, 'uVertexVelocityOld', uVertexOld)
            call mpas_pool_get_array(cdgPool, 'vVertexVelocityOld', vVertexOld)
            call mpas_pool_get_array(cdgPool, 'wVertexVelocityOld', wVertexOld)
            call mpas_pool_get_array(cdgPool, 'uVertexVelocityMid', uVertexMid)
            call mpas_pool_get_array(cdgPool, 'vVertexVelocityMid', vVertexMid)
            call mpas_pool_get_array(cdgPool, 'wVertexVelocityMid', wVertexMid)
            call mpas_pool_get_array(cdgPool, 'uVertexVelocityNew', uVertexNew)
            call mpas_pool_get_array(cdgPool, 'vVertexVelocityNew', vVertexNew)
            call mpas_pool_get_array(cdgPool, 'wVertexVelocityNew', wVertexNew)

            ! Integrate the vertices backwards along velocity characteristics to their departure points
            call mpas_pool_get_field(cdgPool, 'uVertexVelocityNew', uVertex)
            call mpas_pool_get_field(cdgPool, 'vVertexVelocityNew', vVertex)
            call mpas_pool_get_field(cdgPool, 'wVertexVelocityNew', wVertex)
            call ocn_vertex_reconstruction(0, meshPool, scratchPool, scratchPool, layerThickness, normVelocitySum, uVertex, vVertex, wVertex)

            if (firstStep) then
               vertAleVelocityOld(:,:) = vertAleVelocityMid(:,:)
               normVelocityOld(:,:) = normVelocitySum(:,:)
               uVertexOld(:,:) = uVertexNew(:,:)
               vVertexOld(:,:) = vVertexNew(:,:)
               wVertexOld(:,:) = wVertexNew(:,:)
               firstStep = .false.
            end if
            layerThicknessMid(:,:) = 0.5_RKIND*(layerThickness(:,:) + layerThicknessOld(:,:))
            vertAleVelocityMid(:,:) = 0.5_RKIND*(vertAleVelocitySum(:,:) + vertAleVelocityOld(:,:))
            normVelocityMid(:,:) = 0.5_RKIND*(normVelocitySum(:,:) + normVelocityOld(:,:))
            uVertexMid(:,:) = 0.5_RKIND*(uVertexNew(:,:) + uVertexOld(:,:))
            vVertexMid(:,:) = 0.5_RKIND*(vVertexNew(:,:) + vVertexOld(:,:))
            wVertexMid(:,:) = 0.5_RKIND*(wVertexNew(:,:) + wVertexOld(:,:))
         end if
         block => block % next
      end do

      block => domain % blockList
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'cdg', cdgPool)

         if (mod(step,do_every) .eq. do_every - 1) then
            call mpas_pool_get_array(cdgPool, 'departurePoints', departurePoints)
            call mpas_pool_get_array(cdgPool, 'departureEdges', departureEdges)

            if (onSphere) then
               call ocn_cdg_char_departure_points(meshPool, cdgPool, dt, departurePoints)
            else
               call ocn_cdg_char_departure_points_planar(meshPool, cdgPool, dt, departurePoints)
               call ocn_cdg_char_integrate_edge_planar(meshPool, cdgPool, dt, departureEdges)
            end if

            call ocn_cdg_basis_init_mean_3d(meshPool, statePool, cdgPool, hqQuad3D, hqBasis3D, 1)
            call ocn_cdg_basis_init_mean_3d(meshPool, statePool, cdgPool, hqQuad3D, hqBasis3D, 2)
         end if

         block => block % next
      end do

   end subroutine ocn_cdg_manager_pre_solve!}}}

   subroutine ocn_cdg_manager_post_solve(domain, dt)!{{{

      type (domain_type), intent(inout) :: domain !< Input/Output
      real (kind=RKIND), intent(in) :: dt

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: cdgPool

      type (block_type), pointer :: block

      integer, pointer :: nCells, nCellsSolve, nVertLevels, config_cdg_do_every
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: qField, layerThicknessOld, layerThickness, normVelocitySum, normVelocityOld, vertAleVelocitySum, vertAleVelocityOld, &
                                                    uVertex, vVertex, wVertex, uVertexOld, vVertexOld, wVertexOld, contFluxEdge, contFluxCell
      real (kind=RKIND), dimension(:,:,:,:), pointer :: hqCoeffs

      integer :: iCell, iLevel
      real (kind=RKIND) :: localMin, localMax, globalMin, globalMax, hqVolLocal, hqVolGlobal, q

      call mpas_pool_get_config(domain % configs, 'config_cdg_do_every', config_cdg_do_every)

#ifdef CDG_TESTING
      if (mod(step,config_cdg_do_every) .eq. config_cdg_do_every - 1) then
         block => domain % blockList
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'cdg', cdgPool)
            call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
            call mpas_pool_get_array(cdgPool, 'hqCoeffsNew', hqCoeffs)

            ! Check the extrema
            localMin = +1.0E+9; localMax = -1.0E+9
            do iCell = 1, nCellsSolve
               do iLevel = 1, maxLevelCell(iCell)
                  q = hqCoeffs(1,1,iLevel,iCell)
                  if (q .lt. localMin) localMin = q
                  if (q .gt. localMax) localMax = q
               end do
            end do

            ! Test mass conservation
            call ocn_cdg_advection_integrate(meshPool, statePool, cdgPool, 2, hqQuad3D, hqBasis3D, hqCoeffs, hqVolLocal)

            block => block % next
         end do
         call mpas_dmpar_min_real(domain % dminfo, localMin, globalMin)
         call mpas_dmpar_max_real(domain % dminfo, localMax, globalMax)
         call mpas_log_write( 'cdg: global extrema (q).. $r $r', realArgs=(/ globalMin, globalMax /) )

         if (hqInitVol .gt. 0.0_RKIND) then
            call mpas_dmpar_sum_real(domain % dminfo, hqVolLocal, hqVolGlobal)
            call mpas_log_write( 'cdg: mass conservation .. $r', realArgs=(/ 1.0_RKIND - hqVolGlobal/hqInitVol /) )
         end if
      end if
#endif

      block => domain % blockList
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'cdg', cdgPool)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(cdgPool, 'qField', qField, 1)
         call mpas_pool_get_array(cdgPool, 'hqCoeffsNew', hqCoeffs)

         do iCell = 1, nCells
            do iLevel = 1, maxLevelCell(iCell)
               qField(iLevel,iCell) = hqCoeffs(1,1,iLevel,iCell)
            end do
         end do

         if (mod(step,config_cdg_do_every) .eq. config_cdg_do_every - 1) then
            call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
            call mpas_pool_get_array(cdgPool, 'layerThicknessOld', layerThicknessOld)
            call mpas_pool_get_array(cdgPool, 'normVelocitySum', normVelocitySum)
            call mpas_pool_get_array(cdgPool, 'vertAleVelocitySum', vertAleVelocitySum)
            call mpas_pool_get_array(cdgPool, 'normVelocityOld', normVelocityOld)
            call mpas_pool_get_array(cdgPool, 'vertAleVelocityOld', vertAleVelocityOld)
            call mpas_pool_get_array(cdgPool, 'contFluxEdge', contFluxEdge)
            call mpas_pool_get_array(cdgPool, 'contFluxCell', contFluxCell)
            call mpas_pool_get_array(cdgPool, 'uVertexVelocityNew', uVertex)
            call mpas_pool_get_array(cdgPool, 'vVertexVelocityNew', vVertex)
            call mpas_pool_get_array(cdgPool, 'wVertexVelocityNew', wVertex)
            call mpas_pool_get_array(cdgPool, 'uVertexVelocityOld', uVertexOld)
            call mpas_pool_get_array(cdgPool, 'vVertexVelocityOld', vVertexOld)
            call mpas_pool_get_array(cdgPool, 'wVertexVelocityOld', wVertexOld)

            layerThicknessOld(:,:) = layerThickness(:,:)
            normVelocityOld(:,:) = normVelocitySum(:,:)
            vertAleVelocityOld(:,:) = vertAleVelocitySum(:,:)
            normVelocitySum(:,:) = 0.0_RKIND
            vertAleVelocitySum(:,:) = 0.0_RKIND
            contFluxEdge(:,:) = 0.0_RKIND
            contFluxCell(:,:) = 0.0_RKIND

            uVertexOld(:,:) = uVertex(:,:)
            vVertexOld(:,:) = vVertex(:,:)
            wVertexOld(:,:) = wVertex(:,:)
         end if

         block => block % next
      end do

      step = step + 1

   end subroutine ocn_cdg_manager_post_solve!}}}

! hq_prev: thickness weighter tracer concentration field at the previous time level
! hq_data: thickness weighter tracer concentration field used to assemble the flux vector
! hq_curr: thickness weighter tracer concentration field at the new time level

   subroutine ocn_cdg_manager_advect(meshPool, statePool, diagnosticsPool, tendPool, scratchPool, cdgPool, dt, limType, doEvery, hqCur, hqNew)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: diagnosticsPool
      type (mpas_pool_type), intent(in) :: tendPool
      type (mpas_pool_type), intent(in) :: scratchPool
      type (mpas_pool_type), intent(inout) :: cdgPool
      integer, intent(in) :: limType, doEvery
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: hqCur
      real (kind=RKIND), dimension(:,:,:,:), intent(out) :: hqNew

      type(mpas_pool_type), pointer :: tracersPool
      logical, pointer :: onSphere
      integer, pointer :: nVertLevels, nEdges, nCells, nCellsSolve, nTracers
      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell, maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, edgesOnCell, highOrderAdvectionMask, edgeSignOnCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: wTop, layerThickness, normalVelocity, tend_layerThickness, contFluxEdge, contFluxCell

      integer :: iCell, iLevel, iTracer, iBasis
      real (kind=RKIND), dimension(:,:,:,:), allocatable :: flux

      call mpas_pool_get_config(meshPool, 'on_a_sphere', onSphere)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_array(cdgPool, 'contFluxEdge', contFluxEdge)
      call mpas_pool_get_array(cdgPool, 'contFluxCell', contFluxCell)

      if (mod(step,doEvery) .ne. doEvery - 1) return

      allocate(flux(hqBasis3D % nf, nTracers, nVertLevels, nCells))
      flux = 0.0_RKIND

      ! Construct the vertical fluxes
      call mpas_timer_start("cdg vert flux")
      call ocn_cdg_advection_vert_swept_region(meshPool, cdgPool, hqQuad3D, hqBasis3D, dt, flux, hqNew)
      call mpas_timer_stop("cdg vert flux")

      ! Construct the horizontal fluxes
      call mpas_timer_start("cdg edge scale")
      call ocn_cdg_advection_edge_scaling(meshPool, diagnosticsPool, cdgPool, dt)
      call mpas_timer_stop("cdg edge scale")
      if (onSphere) then
         call ocn_cdg_advection_swept_regions(meshPool, cdgPool, hqQuad3D, hqBasis3D, dt, flux, hqNew)
      else
         call mpas_timer_start("cdg horiz flux", .false.)
         call ocn_cdg_advection_swept_regions_planar(meshPool, cdgPool, hqQuad3D, hqBasis3D, dt, flux, hqNew)
         call mpas_timer_stop("cdg horiz flux")
      end if

      ! Construct the integral over the elements at the previous time level
      call mpas_timer_start("cdg lhs assemble")
      if (onSphere) then
         call ocn_cdg_advection_prev_flux(meshPool, diagnosticsPool, cdgPool, hqQuad3D, hqBasis3D, dt, flux, hqCur)
      else
         call ocn_cdg_advection_prev_flux_planar(meshPool, cdgPool, hqQuad3D, hqBasis3D, dt, flux, hqCur)
      end if
      call mpas_timer_stop("cdg lhs assemble")

      ! Construct the mass matrices and update the coefficients at the new time level
      call mpas_timer_start("cdg mass matrix")
      !call ocn_cdg_advection_mass_matrix(meshPool, statePool, cdgPool, hqQuad3D, hqBasis3D, flux, hqNew)
      call ocn_cdg_advection_mass_matrix_lumped(meshPool, statePool, cdgPool, hqQuad3D, hqBasis3D, flux, hqNew)
      call mpas_timer_stop("cdg mass matrix")

      ! Apply the slope limiter
      call mpas_timer_start("cdg limiter")
      if (limType .eq. 1) then
         if (onSphere) then
            call ocn_cdg_limiter_bj3d(meshPool, statePool, cdgPool, hqBasis3D, hqNew)
         else
            call ocn_cdg_limiter_bj3d_planar(meshPool, statePool, cdgPool, hqBasis3D, hqNew)
         end if
      else if (limType .eq. 2) then
         call ocn_cdg_limiter_weno(meshPool, statePool, cdgPool, hqBasis3D, hqNew)
      end if
      call mpas_timer_stop("cdg limiter")

      deallocate(flux)

   end subroutine ocn_cdg_manager_advect!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_manager_test
!
!> \brief MPAS ocean CDG advection tests
!> \author Dave Lee
!> \date   December 2015
!> \details
!>  This routine launches the tests for the different CDG modules
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_manager_test(domain, iErr)!{{{

      type (domain_type), intent(in) :: domain
      integer, intent(out) :: iErr

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: cdgPool

      real (kind=RKIND), dimension(:,:,:,:), pointer :: hqCoeffs

      real (kind=RKIND) :: volLocal, volGlobal, volA, pi

      iErr = 0

      pi = 3.1415926535897932384626433832795028841971

      ! Test the different quadrature objects
      call ocn_cdg_tri_test(domain, hqQuad, iErr)

      ! Test the basis function means
      call ocn_cdg_basis_test(domain, hqQuad, hqBasis, iErr)

      ! Test the volume of the tracer field determined from the initial
      ! condition
      block => domain % blockList
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'cdg', cdgPool)

         call mpas_pool_get_array(cdgPool, 'hqCoeffsNew', hqCoeffs)

         call ocn_cdg_advection_integrate(meshPool, statePool, cdgPool, 2, hqQuad3D, hqBasis3D, hqCoeffs, volLocal)

         block => block % next
      end do

      call mpas_dmpar_sum_real(domain % dminfo, volLocal, volGlobal)

      volA = 2.0_RKIND*sqrt(2.0_RKIND)*pi/15.0_RKIND

      ! Write result to the log file
      call mpas_log_write( 'cdg: advection unit test volume $r $r $r', realArgs=(/ volGlobal, volA, 1.0_RKIND - volGlobal/volA /) )

   end subroutine ocn_cdg_manager_test!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_manager_finalize
!
!> \brief MPAS ocean CDG advection finalize routine
!> \author Dave Lee
!> \date   December 2015
!> \details
!>  This routine calculates the normalized L_2 error for the tracer
!>  field after a single rotation and cleans up the local memory 
!>  associated with each of the quadrature and basis objects
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_manager_finalize(domain, iErr)!{{{

      type (domain_type), intent(in) :: domain
      integer, intent(out) :: iErr

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: cdgPool
      type (block_type), pointer :: block

      real (kind=RKIND) :: errLocal, normLocal, errGlobal, normGlobal
      real (kind=RKIND), dimension(:,:,:,:), pointer :: hqCoeffs

      iErr = 0

#ifdef CDG_PRESCRIBED_VELOCITY
      ! Calculate the errors against the analytic function
      block => domain % blockList
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', statePool)
         call mpas_pool_get_subpool(block % structs, 'cdg', cdgPool)
         call mpas_pool_get_array(cdgPool, 'hqCoeffsNew', hqCoeffs)
         call ocn_cdg_advection_error(meshPool, statePool, cdgPool, hqQuad3D, hqBasis3D, hqCoeffs, errLocal, normLocal)
         block => block % next
      end do

      call mpas_log_write( 'cdg: thickness weighted tracer field.....' )
      call mpas_dmpar_sum_real(domain % dminfo, errLocal, errGlobal)
      call mpas_dmpar_sum_real(domain % dminfo, normLocal, normGlobal)
      call mpas_log_write( 'cdg: normalized L_2 error norm for basis $r $r', realArgs=(/ sqrt(errGlobal)/sqrt(normGlobal) /) )
#endif

      deallocate(hqQuad3D % qw, hqQuad3D % qx, hqQuad3D % qy, hqQuad3D % qz, hqQuad3D % qwv, hqQuad3D % qzv)
      deallocate(hqBasis3D % power_x, hqBasis3D % power_y, hqBasis3D % power_z)

   end subroutine ocn_cdg_manager_finalize!}}}

   subroutine ocn_cdg_manager_increment_cont_flux(meshPool, diagnosticsPool, cdgPool, dt, contFluxEdge, contFluxCell, normVelocitySum, vertAleVelocitySum)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: diagnosticsPool
      type (mpas_pool_type), intent(in) :: cdgPool
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:), intent(inout) :: contFluxEdge, contFluxCell, normVelocitySum, vertAleVelocitySum

      integer, pointer :: nEdges, nCells
      integer, dimension(:), pointer :: maxLevelEdgeTop, maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, normalVelocity, vertAleTransportTop, layerThicknessEdge

      integer :: iEdge, iLevel, iCell, jCell
      real (kind=RKIND) :: hEdge

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)

      ! Increment the horizontal edge volume fluxes
      do iEdge = 1, nEdges
         do iLevel = 1, maxLevelEdgeTop(iEdge)
            contFluxEdge(iLevel,iEdge) = contFluxEdge(iLevel,iEdge) + &
               dt*normalVelocity(iLevel,iEdge)*dvEdge(iEdge)*layerThicknessEdge(iLevel,iEdge)
         end do
      end do

      ! Increment the vertical cell volume fluxes
      do iCell = 1, nCells
         do iLevel = 2, maxLevelCell(iCell)
            ! Note: this is a volume flux per unit cell area!
            contFluxCell(iLevel,iCell) = contFluxCell(iLevel,iCell) + dt*vertAleTransportTop(iLevel,iCell)
         end do
      end do

      ! Increment the CDG transport velocity fields
      normVelocitySum(:,:) = normVelocitySum(:,:) + normalVelocity(:,:)
      vertAleVelocitySum(:,:) = vertAleVelocitySum(:,:) + vertAleTransportTop(:,:)

   end subroutine ocn_cdg_manager_increment_cont_flux!}}}

   subroutine ocn_cdg_manager_init_periodic_geom(domain)!{{{

      type (domain_type), intent (inout) :: domain !< Input/Output

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool

      logical :: check
      integer :: sgn, iVertex, iEdge, iCell, jCell, kCell, cell1, cell2

      integer, pointer :: nVertices, nVerticesSolve, nEdges, nEdgesSolve, nCells, nCellsSolve, vertexDegree
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, cellsOnVertex
      real (kind=RKIND), dimension(:), pointer :: xVertex, xEdge, xCell

      integer, dimension(:), allocatable :: side, sideOld

!#ifdef CDG_RESTART
!return
!#endif

      if (domain % on_a_sphere) return

      ! Even processors are on the left side of the domain, odd on the right
      if (mod(domain % dminfo % my_proc_id,2) .eq. 0) then
         sgn = -1
      else
         sgn = +1
      end if

      ! Loop through the halo cells and update the geometry
      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

         call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         
         call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
         call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
         call mpas_pool_get_array(meshPool, 'xCell', xCell)

         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)

         allocate(side(nCells-nCellsSolve), sideOld(nCells-nCellsSolve))
         side(:) = 0
         sideOld(:) = 0

         ! Determine which halo cells are on the left (-1) and right
         ! (+1) sides of the domain. Assumes a halo depth of 3.

         ! First layer of halos
        do iCell = nCellsSolve + 1, nCells
           do jCell = 1, nEdgesOnCell(iCell)    
              kCell = cellsOnCell(jCell,iCell)

              ! Left side
              if (mod(kCell,NX_LOCAL) .eq. 1 .and. kCell .le. nCellsSolve) then
                 side(iCell-nCellsSolve) = -1
              ! Right side
              else if (mod(kCell,NX_LOCAL) .eq. 0 .and. kCell .le. nCellsSolve) then
                 side(iCell-nCellsSolve) = +1
              end if
           end do
        end do
        sideOld(:) = side(:)

        ! Second layer of halos
        do iCell = nCellsSolve + 1, nCells
           ! Skip existing halos
           if (sideOld(iCell-nCellsSolve) .ne. 0) cycle

           do jCell = 1, nEdgesOnCell(iCell)
              kCell = cellsOnCell(jCell,iCell)
              if (kCell .le. nCellsSolve .or. kCell .gt. nCells) cycle 

              if (sideOld(kCell-nCellsSolve) .eq. -1) then
                 side(iCell-nCellsSolve) = -1
              else if (sideOld(kCell-nCellsSolve) .eq. +1) then
                 side(iCell-nCellsSolve) = +1
              end if 
           end do
        end do
        sideOld(:) = side(:)

        ! Third layer of halos
        do iCell = nCellsSolve + 1, nCells
           ! Skip existing halos
           if (sideOld(iCell-nCellsSolve) .ne. 0) cycle

           do jCell = 1, nEdgesOnCell(iCell)
              kCell = cellsOnCell(jCell,iCell)
              if (kCell .le. nCellsSolve .or. kCell .gt. nCells) cycle 

              if (sideOld(kCell-nCellsSolve) .eq. -1) then
                 side(iCell-nCellsSolve) = -1
              else if (sideOld(kCell-nCellsSolve) .eq. +1) then
                 side(iCell-nCellsSolve) = +1
              end if 
           end do
        end do

         ! Move the left halo cells left on the left decomp and 
         ! move the right halo cells right on the right decomp
         do iCell = nCellsSolve + 1, nCells
            if (side(iCell - nCellsSolve) .eq. sgn .and. &
               (xCell(iCell) .lt. xCell(1) - 0.5_RKIND*DELX - EPS .or. &
                xCell(iCell) .gt. xCell(NX_LOCAL) + 0.5_RKIND*DELX + EPS)) then
               xCell(iCell) = xCell(iCell) + sgn*XLEN
            end if
         end do

         !do iEdge = nEdgesSolve + 1, nEdges
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            if (cell1 .le. nCells .and. abs(xEdge(iEdge) - xCell(cell1)) .gt. 0.5_RKIND*XLEN) then
               xEdge(iEdge) = xEdge(iEdge) + sgn*XLEN
            else if (cell2 .le. nCells .and. abs(xEdge(iEdge) - xCell(cell2)) .gt. 0.5_RKIND*XLEN) then
               xEdge(iEdge) = xEdge(iEdge) + sgn*XLEN
            end if
         end do

         !do iVertex = nVerticesSolve + 1, nVertices
         do iVertex = 1, nVertices
            check = .true.
            do iCell = 1, vertexDegree
               jCell = cellsOnVertex(iCell,iVertex)

               if (.not. check) cycle
               if (jCell .gt. nCells) cycle

               if (abs(xVertex(iVertex) - xCell(jCell)) .gt. 0.5_RKIND*XLEN) then
                  xVertex(iVertex) = xVertex(iVertex) + sgn*XLEN
                  check = .false.
               end if
            end do
         end do

         deallocate(side, sideOld)

         block => block % next
      end do

   end subroutine ocn_cdg_manager_init_periodic_geom!}}}

end module ocn_cdg_manager

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
