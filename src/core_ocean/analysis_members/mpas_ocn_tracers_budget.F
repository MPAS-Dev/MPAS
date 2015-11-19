! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracers_budget
!
!> \brief MPAS ocean analysis mode member: tracers_budget
!> \author Milena Veneziani
!> \date   Sept 2015
!> \details
!>  MPAS ocean analysis mode member: tracers_budget
!>  
!>  5. In src/core_ocean/analysis_members/Makefile, add your
!>     new analysis member to the list of members. See another analysis member
!>     in that file for an example. 
!>     NOTE: If your analysis member depends on other files, add a dependency
!>           line for the member and list them there. See okubo weiss for an example.
!>
!-----------------------------------------------------------------------

module ocn_tracers_budget

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_routines
   use ocn_tracer_hmix
   use ocn_high_freq_thickness_hmix_del2
   use ocn_tracer_advection
   use ocn_tracer_short_wave_absorption
   use ocn_tracer_nonlocalflux
   use ocn_tracer_surface_restoring
   use ocn_tracer_interior_restoring
   use ocn_tracer_exponential_decay
   use ocn_tracer_ideal_age
   use ocn_tracer_TTD
   use ocn_tracer_surface_flux_to_tend
   use ocn_surface_bulk_forcing
   use ocn_surface_land_ice_fluxes

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

   public :: ocn_init_tracers_budget, &
             ocn_compute_tracers_budget, &
             ocn_restart_tracers_budget, &
             ocn_finalize_tracers_budget

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_tracers_budget
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Milena Veneziani
!> \date    Sept 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_tracers_budget(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_init_tracers_budget!}}}

!***********************************************************************
!
!  routine ocn_compute_tracers_budget
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Milena Veneziani
!> \date    Sept 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_tracers_budget(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      ! Note: variables that are needed at 2 time-levels (n and n+1) have been
      ! named 'Cur' if time-level=n and 'New' if time-level=n+1. Now, the
      ! forward model time integration step has already been done when this 
      ! AM is computed, BUT time-levels have also been switched before the AM
      ! is called, therefore: timeLevel=1 corresponds to the New time (n+1) and
      ! timeLevel=2 corresponds to the Old (or Cur) time (n). Finally, any other 
      ! 2-time-level variable that does not have a 'Cur' or a 'New' suffix 
      ! is associated with time-level=n+1.
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: tracersBudgetAMPool
      type (mpas_pool_type), pointer :: tracersBudgetScratchPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: forcingPool

      !
      ! additional pools
      !
      type (mpas_pool_type), pointer :: tracersPool, tracersTendPool            ! tracers and their tendency
      type (mpas_pool_type), pointer :: tracersSurfaceFluxPool                  ! surface fluxes
      type (mpas_pool_type), pointer :: tracersSurfaceRestoringFieldsPool       ! surface restoring
      type (mpas_pool_type), pointer :: tracersInteriorRestoringFieldsPool      ! interior restoring
      type (mpas_pool_type), pointer :: tracersExponentialDecayFieldsPool       ! exponential decay
      type (mpas_pool_type), pointer :: tracersIdealAgeFieldsPool               ! ideal age
      type (mpas_pool_type), pointer :: tracersTTDFieldsPool                    ! transit time distribution

      ! scalar pointers
      integer :: nTracerGroup, iCell, iEdge, k, iTracer
      integer, pointer :: nVertLevels, nEdges, nCellsSolve, indexTemperature
! nEdgesSolve, nVerticesSolve
      logical, pointer :: config_disable_tr_all_tend, config_use_cvmix_kpp
      logical, pointer :: config_use_tracerGroup, config_use_tracerGroup_surface_bulk_forcing, config_use_tracerGroup_surface_restoring, &
                          config_use_tracerGroup_interior_restoring, config_use_tracerGroup_exponential_decay, config_use_tracerGroup_idealAge_forcing, &
                          config_use_tracerGroup_ttd_forcing
      integer :: timeLevelCur, timeLevelNew

      ! time step
      type (MPAS_timeInterval_type) :: timeStep
      real (kind=RKIND) :: dt
      character(len=StrKIND), pointer :: config_dt

      ! iterator for tracer categories
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: modifiedConfigName

      !
      ! one dimensional pointers
      !
      integer, dimension(:), pointer :: maxLevelCell
!, maxLevelEdgeTop, maxLevelVertexBot
      real (kind=RKIND), dimension(:), pointer :: penetrativeTemperatureFlux
      real (kind=RKIND), dimension(:), pointer :: tracerGroupExponentialDecayRate
!      real (kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge

      !
      ! two dimensional pointers
      !
      real (kind=RKIND), dimension(:,:), pointer :: tracerGroupPistonVelocity, tracerGroupSurfaceRestoringValue, tracerGroupIdealAgeMask, tracerGroupTTDMask

      real (kind=RKIND), dimension(:,:), pointer :: &
        normalTransportVelocity, layerThicknessCur, layerThicknessNew, vertAleTransportTop, layerThicknessEdge, vertDiffTopOfCell, &
        tend_layerThickness, normalThicknessFlux, tracerGroupSurfaceFlux, fractionAbsorbed, zMid, relativeSlopeTopOfEdge, &
        relativeSlopeTapering, relativeSlopeTaperingCell, tracerProvisional

      !
      ! three dimensional pointers
      !
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        tracerGroupCur, tracerGroupNew, tracerGroupTend, vertNonLocalFlux, &
        tracerGroupBudgetHadv, tracerGroupBudgetVadv, tracerGroupBudgetHdif, &
        tracerGroupBudgetVdif, tracerGroupBudgetForc, tracerGroupBudgetRate

      real (kind=RKIND), dimension(:,:,:), pointer :: tracerGroupInteriorRestoringRate, tracerGroupInteriorRestoringValue

      !
      ! Field pointers
      !
      type (field2DReal), pointer :: normalThicknessFluxField
      type (field2DReal), pointer :: tracerProvisionalField

      err = 0

      ! get dt
      call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)
      call mpas_set_timeInterval(timeStep, timeString=config_dt, ierr=err)
      call mpas_get_timeInterval(timeStep, dt=dt)

      ! Define time levels:
      timeLevelNew = 1 ! n+1
      timeLevelCur = 2 ! n

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         !
         ! get pools
         !
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
         call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)
         call mpas_pool_get_subpool(block % structs, 'tracersBudgetAM', tracersBudgetAMPool)
         call mpas_pool_get_subpool(block % structs, 'tracersBudgetScratch', tracersBudgetScratchPool)

         !
         ! get dimensions
         !
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
!         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
!         call mpas_pool_get_dimension(block % dimensions, 'nVerticesSolve', nVerticesSolve)
         call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)

         !
         ! get arrays
         !
!         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
!         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
!         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
!         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
!         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, timeLevelCur)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, timeLevelNew)
         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
         call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
         call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
         call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
         call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTopOfEdge', relativeSlopeTopOfEdge)
         call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTapering', relativeSlopeTapering)
         call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTaperingCell', relativeSlopeTaperingCell)
         call mpas_pool_get_array(diagnosticsPool, 'vertNonLocalFlux', vertNonLocalFlux)
         call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)
         call mpas_pool_get_array(forcingPool, 'fractionAbsorbed', fractionAbsorbed)
         call mpas_pool_get_array(tendPool, 'layerThickness', tend_layerThickness)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         !
         ! get configure options
         !
         call mpas_pool_get_config(ocnConfigs, 'config_disable_tr_all_tend', config_disable_tr_all_tend)
         call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)


         call mpas_pool_get_field(scratchPool, 'normalThicknessFlux', normalThicknessFluxField)
         call mpas_allocate_scratch_field(normalThicknessFluxField, .true.)
         normalThicknessFlux => normalThicknessFluxField % array

         call mpas_pool_get_field(tracersBudgetScratchPool, 'tracerProvisional', tracerProvisionalField)
         call mpas_allocate_scratch_field(tracerProvisionalField, .true.)
         tracerProvisional => tracerProvisionalField % array

         if(config_disable_tr_all_tend) return

         !
         ! transport velocity for the tracer.
         !
         do iEdge = 1, nEdges
            do k = 1, nVertLevels
               normalThicknessFlux(k, iEdge) = normalTransportVelocity(k, iEdge) * layerThicknessEdge(k, iEdge)
            end do
         end do

         !
         ! begin iterate over tracer categories
         !
         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
           if ( groupItr % memberType == MPAS_POOL_FIELD ) then

              ! load configure setting for this category
              !
              modifiedConfigName = 'config_use_' // trim(groupItr % memberName)
              call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup)

              if ( config_use_tracerGroup ) then
                 modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_surface_bulk_forcing' 
                 call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_surface_bulk_forcing)
                 modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_surface_restoring'
                 call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_surface_restoring)
                 modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_interior_restoring'
                 call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_interior_restoring)
                 modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_exponential_decay'
                 call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_exponential_decay)
                 modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_idealAge_forcing'
                 call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_idealAge_forcing)
                 modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_ttd_forcing'
                 call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_ttd_forcing)


                 ! Get tracer group, and other groups (tendencies, etc.)
                 call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracerGroupCur, timeLevelCur)
                 call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracerGroupNew, timeLevelNew)
                 nTracerGroup = size(tracerGroupCur, dim=1)

                 ! Get Tendency array
                 modifiedGroupName = trim(groupItr % memberName) // "Tend"
                 call mpas_pool_get_array(tracersTendPool, trim(modifiedGroupName), tracerGroupTend)

                 ! Get surface flux array
                 modifiedGroupName = trim(groupItr % memberName) // "SurfaceFlux"
                 call mpas_pool_get_array(tracersSurfaceFluxPool, trim(modifiedGroupName), tracerGroupSurfaceFlux)

                 ! Get tracer budget arrays
                 modifiedGroupName = trim(groupItr % memberName) // "BudgetHadv"
                 call mpas_pool_get_array(tracersBudgetAMPool, trim(modifiedGroupName), tracerGroupBudgetHadv)
!                 modifiedGroupName = trim(groupItr % memberName) // "BudgetVadv"
!                 call mpas_pool_get_array(tracersBudgetAMPool, trim(modifiedGroupName), tracerGroupBudgetVadv)
                 modifiedGroupName = trim(groupItr % memberName) // "BudgetHdif"
                 call mpas_pool_get_array(tracersBudgetAMPool, trim(modifiedGroupName), tracerGroupBudgetHdif)
                 modifiedGroupName = trim(groupItr % memberName) // "BudgetVdif"
                 call mpas_pool_get_array(tracersBudgetAMPool, trim(modifiedGroupName), tracerGroupBudgetVdif)
                 modifiedGroupName = trim(groupItr % memberName) // "BudgetForc"
                 call mpas_pool_get_array(tracersBudgetAMPool, trim(modifiedGroupName), tracerGroupBudgetForc)
                 modifiedGroupName = trim(groupItr % memberName) // "BudgetRate"
                 call mpas_pool_get_array(tracersBudgetAMPool, trim(modifiedGroupName), tracerGroupBudgetRate)

                 !
                 ! initialize tracer surface flux and tendencies to zero.
                 !
                 tracerGroupTend(:,:,:) = 0.0
                 tracerGroupSurfaceFlux(:,:) = 0.0
                 tracerGroupBudgetHadv(:,:,:) = 0.0
!                 tracerGroupBudgetVadv(:,:,:) = 0.0
                 tracerGroupBudgetHdif(:,:,:) = 0.0
                 tracerGroupBudgetVdif(:,:,:) = 0.0
                 tracerGroupBudgetForc(:,:,:) = 0.0
                 tracerGroupBudgetRate(:,:,:) = 0.0

                 !
                 ! tracer tendency: external forcing term
                 !
                 !
                 ! First, fill components of surface tracer flux
                 !
                 if (config_use_tracerGroup_surface_bulk_forcing) then
                    call ocn_surface_bulk_forcing_tracers(meshPool, groupItr % memberName, forcingPool, tracerGroupSurfaceFlux, err)
                 end if

                 !
                 ! ocean surface restoring
                 !
                 if (config_use_tracerGroup_surface_restoring) then 
                     call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "PistonVelocity"
                     call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, trim(modifiedGroupName), tracerGroupPistonVelocity)
                     modifiedGroupName = trim(groupItr % memberName) // "SurfaceRestoringValue"
                     call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, trim(modifiedGroupName), tracerGroupSurfaceRestoringValue)
                     call ocn_tracer_surface_restoring_compute(nTracerGroup, nCellsSolve, tracerGroupNew, tracerGroupPistonVelocity, tracerGroupSurfaceRestoringValue, tracerGroupSurfaceFlux, err)
                 endif

                 ! land-ice / ocean interface flux
                 ! this is a flux at the top ocean surface -- so these fluxes should be added into tracerGroupSurfaceFlux
                 !if (put correct logic here, only 'active' when coupling is turned on)
                 !   call ocn_tracer_landIce_ocean_coupling(tracerGroupNew, tracerGroupSurfaceFlux)
                 !endif

                 !
                 ! other additions to tracerGroupSurfaceFlux should be added here
                 !

                 !
                 ! Second, convert the surface tracer flux into a tracer tendency by distributing the flux across some number of surface layers
                 !
                 call ocn_tracer_surface_flux_tend(meshPool, fractionAbsorbed, layerThicknessNew, tracerGroupSurfaceFlux, tracerGroupBudgetForc, err)

                 !
                 ! Third, add interior forcing terms
                 !
                 !
                 ! interior restoring forcing tendency
                 !
                 if (config_use_tracerGroup_interior_restoring) then
                     call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "InteriorRestoringRate"
                     call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, trim(modifiedGroupName), tracerGroupInteriorRestoringRate)
                     modifiedGroupName = trim(groupItr % memberName) // "InteriorRestoringValue"
                     call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, trim(modifiedGroupName),tracerGroupInteriorRestoringValue)
                     call ocn_tracer_interior_restoring_compute(nTracerGroup, nCellsSolve, maxLevelCell, layerThicknessNew, &
                        tracerGroupNew, tracerGroupInteriorRestoringRate, tracerGroupInteriorRestoringValue, tracerGroupBudgetForc, err)
                 endif

                 !
                 ! exponential decay tendency
                 !
                 if (config_use_tracerGroup_exponential_decay) then
                     write (stderrUnit,'(a)') 'exponential decay not fully tested'
                     call mpas_pool_get_subpool(forcingPool, 'tracersExponentialDecayFields', tracersExponentialDecayFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "ExponentialDecayRate"
                     call mpas_pool_get_array(tracersExponentialDecayFieldsPool, trim(modifiedGroupName), tracerGroupExponentialDecayRate)
                     call ocn_tracer_exponential_decay_compute(nTracerGroup, nCellsSolve, maxLevelCell, layerThicknessNew, &
                        tracerGroupNew, tracerGroupExponentialDecayRate, tracerGroupBudgetForc, err)
                 endif

                 !
                 ! ideal age forcing tendency
                 !   note: ocn_tracer_ideal_age_compute resets tracers in top layer to zero
                 !
                 if (config_use_tracerGroup_idealAge_forcing) then
                     write (stderrUnit,'(a)') 'ideal age not fully tested'
                     call mpas_pool_get_subpool(forcingPool, 'tracersIdealAgeFields', tracersIdealAgeFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "IdealAgeMask"
                     call mpas_pool_get_array(tracersIdealAgeFieldsPool, trim(modifiedGroupName), tracerGroupIdealAgeMask)
                     call ocn_tracer_ideal_age_compute(nTracerGroup, nCellsSolve, maxLevelCell, layerThicknessNew, &
                        tracerGroupIdealAgeMask, tracerGroupNew, tracerGroupBudgetForc, err)
                 endif

                 !
                 ! transit-time distribution (TTD) forcing tendency
                 !   note: no tendency is actually computed in ocn_tracer_TTD_compute
                 !   note: rather, tracerGroupCur is reset to tracerGroupTTDMask in top-most layer
                 !   Milena Veneziani's note: check what this would entail for
                 !   tracers budget AM
                 !if (config_use_tracerGroup_ttd_forcing) then
                 !    write (stderrUnit,'(a)') 'ideal age not fully tested'
                 !    call mpas_pool_get_subpool(forcingPool, 'tracersTTDFields', tracersTTDFieldsPool)
                 !    modifiedGroupName = trim(groupItr % memberName) // "TTDMask"
                 !    call mpas_pool_get_array(tracersTTDFieldsPool, trim(modifiedGroupName), tracerGroupTTDMask)
                 !    call ocn_tracer_TTD_compute(nTracerGroup, nCellsSolve, maxLevelCell, layerThicknessNew, &
                 !       tracerGroupTTDMask, tracerGroupNew, err)
                 !endif

                 !
                 ! Performing shortwave absorption
                 !
                 if ( trim(groupItr % memberName) == 'activeTracers' ) then
                    call ocn_tracer_short_wave_absorption_tend(meshPool, indexTemperature, layerThicknessNew, penetrativeTemperatureFlux, tracerGroupBudgetForc, err)
                 endif

                 !
                 ! tracer tendency: horizontal advection term -div( layerThicknessCur \phi u)
                 !
                 ! Monotonic Advection, or standard advection
                 call ocn_tracer_advection_tend(tracerGroupNew, normalThicknessFlux, vertAleTransportTop, &
                    layerThicknessNew, layerThicknessNew, dt, meshPool, scratchPool, tend_layerThickness, tracerGroupBudgetHadv)
!                  tracerGroupBudgetHadv, tracerGroupBudgetVadv)

                 !
                 ! tracer tendency: del2 horizontal tracer diffusion, div(h \kappa_2 \nabla \phi)
                 !
                 call ocn_tracer_hmix_tend(meshPool, scratchPool, layerThicknessEdge, zMid, tracerGroupNew, &
                    relativeSlopeTopOfEdge, relativeSlopeTapering, relativeSlopeTaperingCell, tracerGroupBudgetHdif, err)

                 !
                 ! tracer tendency: vertical diffusion term
                 !
                 !
                 ! Compute tracer tendency due to non-local flux computed in KPP
                 ! (explicit vertical diffusion term)
                 !
                 if (config_use_cvmix_kpp) then
                   call ocn_tracer_nonlocalflux_tend(meshPool, vertNonLocalFlux, tracerGroupSurfaceFlux, tracerGroupBudgetVdif, err)
                 end if

                 !
                 ! Add implicit vertical diffusion term
                 ! (the final result will be already in tracer-budget-term units)
                 !
                 ! Within the tracer loop, also do the following:
                 ! 1) compute time rate-of-change term
                 ! 2) convert all tendencies to tracer budget terms
                 !
                 do iTracer = 1, nTracerGroup
                   do iCell = 1, nCellsSolve
                     do k = 1, nVertLevels
                   write(*,*) '***********'
                   write(*,*) 'k, iCell, iTracer, layerThicknessCur, layerThicknessNew = ', &
                    k,iCell,iTracer,layerThicknessCur(k, iCell),layerThicknessNew(k, iCell)
                   write(*,*) '***********'
                       tracerProvisional(k, iCell) = ( tracerGroupCur(iTracer, k, iCell) * layerThicknessCur(k, iCell) + &
                                                       dt * ( tracerGroupBudgetForc(iTracer, k, iCell) +                 &
                                                              tracerGroupBudgetHadv(iTracer, k, iCell) +                 &
                                                              tracerGroupBudgetHdif(iTracer, k, iCell) +                 &
                                                              tracerGroupBudgetVdif(iTracer, k, iCell) ) ) / layerThicknessNew(k, iCell)

                       tracerGroupBudgetVdif(iTracer, k, iCell) = tracerGroupBudgetVdif(iTracer, k, iCell) / layerThicknessNew(k, iCell) + &
                                                                  ( tracerGroupNew(iTracer, k, iCell) - tracerProvisional(k, iCell)) / dt

                       tracerGroupBudgetRate(iTracer, k, iCell) = ( tracerGroupNew(iTracer, k, iCell) - tracerGroupCur(iTracer, k, iCell) ) / dt

                       tracerGroupBudgetForc(iTracer, k, iCell) = tracerGroupBudgetForc(iTracer, k, iCell) / layerThicknessNew(k, iCell)
                       tracerGroupBudgetHadv(iTracer, k, iCell) = tracerGroupBudgetHadv(iTracer, k, iCell) / layerThicknessNew(k, iCell)
                       tracerGroupBudgetHdif(iTracer, k, iCell) = tracerGroupBudgetHdif(iTracer, k, iCell) / layerThicknessNew(k, iCell)
                     end do
                   end do
                 end do

              end if
           end if
         end do
         !
         ! end iterate over tracer categories
         !

         !
         ! deallocate workspace
         !
         call mpas_deallocate_scratch_field(normalThicknessFluxField, .true.)
         call mpas_deallocate_scratch_field(tracerProvisionalField, .true.)

         block => block % next
      end do

   end subroutine ocn_compute_tracers_budget!}}}

!***********************************************************************
!
!  routine ocn_restart_tracers_budget
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Milena Veneziani
!> \date    Sept 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_tracers_budget(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_tracers_budget!}}}

!***********************************************************************
!
!  routine ocn_finalize_tracers_budget
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Milena Veneziani
!> \date    Sept 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_tracers_budget(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_finalize_tracers_budget!}}}

end module ocn_tracers_budget

! vim: foldmethod=marker
