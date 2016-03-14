! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_mep
!
!> \brief Ocean meso-scale eddy parameterization module
!> \author Juan A. Saenz, Todd Ringler
!> \date   14 September 2015
!> \details
!>  This module contains interface routines for computing the meso-scale eddy
!>  parameterization.
!
!-----------------------------------------------------------------------

module ocn_mep

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_vector_operations
   use ocn_constants
   use ocn_mep_gm
   use ocn_diagnostics
   use ocn_mep_vmix_mom

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

   public :: ocn_mep_init, &
             ocn_mep_climate_time_filter_update, &
             ocn_mep_slow_time_filter_update, &
             ocn_mep_slow_time_filter_accum, &
             ocn_mep_update, &
             ocn_mep_transportVelocity_update, &
             ocn_mep_vel_tend, &
             ocn_mep_reconstruct, &
             ocn_mep_vel_hmix_tend, &
             ocn_mep_vmix_coefs_build

   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: mepStdGMOn, mepForceTWAdragGMOn, mepForceVertViscOn, mepImplicitVertVisc
   logical :: mepRediOn, mepUseGMmodule
   integer :: NstepsSlowTimeFilter


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_mep_init
!
!> \brief  Initialize the mep module
!> \author Juan A. Saenz, Todd Ringler
!> \date   14 September 2015
!> \details
!> This routine initializes control flags, parameters and quantities related
!> to the meso-scale eddy parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_init(domain, err)!{{{

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
      type (mpas_pool_type), pointer :: mepPool
      type (block_type), pointer :: block
      integer :: counter
      logical, pointer :: config_mep_Redi_enable
      character (len=StrKIND), pointer :: config_mep_formulation
      real (kind=RKIND), pointer :: config_mep_slowtime_Tc, config_dt
      real (kind=RKIND), pointer :: uSlowTimeFilter, uSlowTimeFilterAccum, eddyViscosityForce

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_mep_formulation', config_mep_formulation)
      call mpas_pool_get_config(ocnConfigs, 'config_mep_Redi_enable', config_mep_Redi_enable)

      mepStdGMOn = .false.
      mepForceTWAdragGMOn = .false.
      mepForceVertViscOn = .false.
      mepRediOn = .false.
      mepUseGMmodule = .false.
      mepImplicitVertVisc = .false.

      if ( trim(config_mep_formulation) == 'none' ) then
         return
      end if
      counter = 0
      if ( trim(config_mep_formulation) == 'standardGM' ) then
         counter = counter + 1
         mepStdGMOn = .true.
         mepUseGMmodule = .true.
         call ocn_gm_init(err)
      end if
      if ( trim(config_mep_formulation) == 'twaDragGM' ) then
         counter = counter + 1
         mepForceTWAdragGMOn = .true.
         mepUseGMmodule= .true.
         call ocn_gm_init(err)
      end if
      if ( trim(config_mep_formulation) == 'forceVertVisc' ) then
         counter = counter + 1
         mepForceVertViscOn = .true.
         call ocn_mep_vmix_mom_init(err)
         call mpas_pool_get_config(ocnConfigs, 'config_mep_slowtime_Tc', config_mep_slowtime_Tc)
         call mpas_pool_get_config(ocnConfigs, 'config_dt', config_dt)
         NstepsSlowTimeFilter = floor( config_mep_slowtime_Tc/config_dt )
         
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mep', mepPool)
            call mpas_pool_get_array(mepPool, 'uSlowTimeFilterAccum', uSlowTimeFilterAccum)
            call mpas_pool_get_array(mepPool, 'uSlowTimeFilter', uSlowTimeFilter)
            call mpas_pool_get_array(mepPool, 'eddyViscosityForce', eddyViscosityForce)
            uSlowTimeFilterAccum = 0.0_RKIND
            uSlowTimeFilter = 0.0_RKIND
            eddyViscosityForce = 0.0_RKIND
            block => block % next
         end do

      end if
      if ( trim(config_mep_formulation) == 'implicitVertVisc' ) then
         counter = counter + 1
         mepImplicitVertVisc = .true.
         call ocn_mep_vmix_mom_init(err)
      end if

      ! make sure incompatible combinations are not requested by user
      if (counter>1) then
         write (stderrUnit,*) 'Invalid value in config_mep_formulation:'
         write (stderrUnit,*) config_mep_formulation
         err = 1
         !juan issue: crash out
      end if

      ! Tracer eddy parameterization flags
      if (config_mep_Redi_enable) then
         mepRediOn = .true.
         mepUseGMmodule = .true.
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_mep_init!}}}


!***********************************************************************
!
!  routine ocn_mep_slow_time_filter_update
!
!> \brief  Update the climate of the flow
!> \author Juan A. Saenz, Todd Ringler
!> \date   March, 2016
!> \details
!> This routine does the following
!> 1) update slow time filtered quantities by retrieving 
!> the work variables that carry the accumulated sums, and normalizing
!> them so to produce the slow time filtered variables.
!> 2) reset the working variables to zero for the next cycle.
!> 3) use the updated slow time filter variables to calculate the 
!> eddyViscosityForce
!> There are two sets of variables: 
!>   one that contains work variables with accumulated quantities
!>   one that contains the time filtered quantities
!> It is assumed that work variables are updated elseshere.
!-----------------------------------------------------------------------

   subroutine ocn_mep_slow_time_filter_update(meshPool, diagnosticsPool, mepPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: diagnostics information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: mepPool structure


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
      real (kind=RKIND), dimension(:,:), pointer :: uSlowTimeFilter, uSlowTimeFilterAccum
      real (kind=RKIND), dimension(:,:), pointer :: divEPFT31SlowTimeFilter, divEPFT31SlowTimeFilterAccum


      if (mepForceVertViscOn) then
         call mpas_pool_get_array(mepPool, 'uSlowTimeFilter', uSlowTimeFilter)
         call mpas_pool_get_array(mepPool, 'uSlowTimeFilterAccum', uSlowTimeFilterAccum)
         call mpas_pool_get_array(mepPool, 'divEPFT31SlowTimeFilter', divEPFT31SlowTimeFilter)
         call mpas_pool_get_array(mepPool, 'divEPFT31SlowTimeFilterAccum', divEPFT31SlowTimeFilterAccum)


         ! Update the slow time filters by normalizing the accumulated quantities
         uSlowTimeFilter = uSlowTimeFilterAccum / NstepsSlowTimeFilter
         divEPFT31SlowTimeFilter = divEPFT31SlowTimeFilterAccum / NstepsSlowTimeFilter

         ! Reset work variables with accumulated quantities
         uSlowTimeFilterAccum = 0.0_RKIND
         divEPFT31SlowTimeFilterAccum = 0.0_RKIND

         ! update the eddy force
         call ocn_mep_vmix_mom_calculate_eddy_force(meshPool, diagnosticsPool, mepPool, err)

      end if


   !--------------------------------------------------------------------

   end subroutine ocn_mep_slow_time_filter_update!}}}


!***********************************************************************
!
!  routine ocn_mep_slow_time_filter_accum
!
!> \brief  Update the climate of the flow
!> \author Juan A. Saenz, Todd Ringler
!> \date   March, 2016
!> \details
!> This routine will 
!> 1) update slow time filtered quantities by retrieving 
!> the work variables that carry the accumulated sums, and normalizing
!> then so to produce the slow time filtered variables.
!> 2) reset the working variables to zero for the next cycle.
!> There are two set of variables: 
!>   one that contains work variables with accumulated quantities
!>   one that contains the time filtered quantities
!> It is assumed that work variables are updated elseshere.
!-----------------------------------------------------------------------

   subroutine ocn_mep_slow_time_filter_accum(statePool, mepPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: statePool !< Input/Output: statePool structure

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: mepPool structure


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
      real (kind=RKIND), dimension(:,:), pointer :: uSlowTimeFilterAccum, normalVelocity
      real (kind=RKIND), dimension(:,:), pointer :: divEPFT31SlowTimeFilterAccum, dEPFT31dz

      err = 0

      if (mepForceVertViscOn) then

         call mpas_pool_get_array(mepPool, 'uSlowTimeFilterAccum', uSlowTimeFilterAccum)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity) ! make sure using right one
         call mpas_pool_get_array(mepPool, 'divEPFT31SlowTimeFilterAccum', divEPFT31SlowTimeFilterAccum)
         call mpas_pool_get_array(mepPool, 'dEPFT31dz', dEPFT31dz)

         ! Update the slow time filters by normalizing the accumulated quantities
         uSlowTimeFilterAccum = uSlowTimeFilterAccum + normalVelocity
         divEPFT31SlowTimeFilterAccum = divEPFT31SlowTimeFilterAccum + dEPFT31dz

      end if


   !--------------------------------------------------------------------

   end subroutine ocn_mep_slow_time_filter_accum!}}}


!***********************************************************************
!
!  routine ocn_mep_climate_time_filter_update
!
!> \brief  Update the climate of the flow
!> \author Juan A. Saenz, Todd Ringler
!> \date   December, 2015
!> \details
!> This routine updates the climate of the flow, variables used by the
!> prognostic mep.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_climate_time_filter_update(statePool, meshPool, diagnosticsPool, mepPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: statePool !< Input: state information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: diagnostics information


      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: mepPool !< Input: MEP information

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


      ! Update the filtered, low frequency variables that are
      ! used to in the EVP and in vertIntEddyEnergy tendency calculations
      call ocn_mep_vmix_mom_update_climate( statePool, meshPool, diagnosticsPool, mepPool, err )

      ! update the solution of the eigenvalue problem
      call ocn_mep_vmix_solve_evp( statePool, meshPool, diagnosticsPool, mepPool, err )

      ! Update the vertically integrated eddy energy by integrating in time
      call ocn_mep_vmix_time_integrate_vertInt_eddy_energy( meshPool, diagnosticsPool, mepPool, err )

      ! Get eddy energy as a function of depth, using evp solution and vertIntEddyEnergy
      call ocn_mep_vmix_eddyEnergy( statePool, meshPool, mepPool, err )

      ! update viscocity from mep
      call ocn_mep_vmix_vert_viscocity( meshPool, mepPool, err )


   !--------------------------------------------------------------------

   end subroutine ocn_mep_climate_time_filter_update!}}}


!***********************************************************************
!
!  routine ocn_mep_update
!
!> \brief  This routine updates the mep module given an ocean state
!> \author Juan A. Saenz, Todd Ringler
!> \date   17 September 2015
!> \details
!>  This routine updates the mep given a current ocean state.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_update(diagnosticsPool, meshPool, mepPool, scratchPool)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: diagnosticsPool !< Input/Output: Diagnostics information
      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: mepPool structure
      type (mpas_pool_type), intent(inout) :: scratchPool !< Input/Output: Scratch structure

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------


      if (mepUseGMmodule) then
         ! might want to give the routine called here a more general name
         call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, mepPool, scratchPool)
      end if


   end subroutine ocn_mep_update!}}}


!***********************************************************************
!
!  routine ocn_mep_transportVelocity_update
!
!> \brief  This routine updates the transport velocity
!> \author Juan A. Saenz, Todd Ringler
!> \date   17 September 2015
!> \details
!>  This routine receives a velocity (eulerian velocity) and adds
!>  the GM bolus velocity, resulting in the residual velocity
!>  It is assumed that the mep module has already been updated.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_transportVelocity_update(normalTransportVelocity, mepPool)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: mepPool !< Input/Output: mepPool structure

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: normalTransportVelocity

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:,:), pointer :: normalGMBolusVelocity


      call mpas_pool_get_array(mepPool, 'normalGMBolusVelocity', normalGMBolusVelocity)


      ! NOTE: This routine is called within block loops

      if (mepStdGMOn) then
         normalTransportVelocity(:,:) = normalTransportVelocity(:,:) + normalGMBolusVelocity(:,:)
      end if


   end subroutine ocn_mep_transportVelocity_update!}}}


!***********************************************************************
!
!  routine ocn_mep_vel_tend
!
!> \brief  This routine increments normalVelocity tend by the tendency due to eddy forces.
!> \author Juan A. Saenz, Todd Ringler
!> \date   October, 2015
!> \details
!>  This routine increments normalVelocity tend by the tendency due to eddy forces.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vel_tend(mepPool, tendPool, meshPool, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool        !< Input: Mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: Tendency structure
      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency structure
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend !< Input/Output: normal velocity tendency

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


      ! calculate the eddy force from parameterized div EPFT and apply it to momentum tendency

      ! force from vertical viscocity
      if (mepForceVertViscOn) then
         call ocn_mep_force_vert_visc(mepPool, tendPool, tend, err)
   end if

      ! TWA GM drag
      if (mepForceTWAdragGMOn) then
         call ocn_mep_force_twa_drag_GM(mepPool, tendPool, meshPool, tend, err)
      end if

   end subroutine ocn_mep_vel_tend!}}}


!***********************************************************************
!
!  routine ocn_mep_force_vert_visc
!
!> \brief  This routine calculates normalVelocity tend due to vertical viscocity force.
!> \author Juan A. Saenz, Todd Ringler
!> \date   Feb, 2016
!> \details
!>  This routine calculates normalVelocity tend due to vertical viscosity force.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_force_vert_visc(mepPool, tendPool, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: mepPool  !< Input: Diagnostic information

      !jas issue: it might be that only need uBolus or normalBolusGM,
      !not entire diagnostics?

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency structure
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend          !< Input/Output: normal velocity tendency

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

      real (kind=RKIND), dimension(:,:), pointer :: tend_vel_mep, eddyViscosityForce

      err = 0

      ! retrieve velocity tendencry by MEP
      call mpas_pool_get_array(tendPool, 'normalVelocityMEP', tend_vel_mep)
      call mpas_pool_get_array(mepPool, 'eddyViscosityForce', eddyViscosityForce)

      tend_vel_mep = 0.0

      ! calculate normal velocity time tendency do to vertical eddy viscosity force
      ! right now, only including form drag terms, i.e. vertical transfer of horizontal momentum
      ! tendency = -dE3i/dz dot ni = - d/dz( vertVisc * dubar/dz )
      tend_vel_mep = - eddyViscosityForce

      ! update tendency with div EPFT parameterized eddy viscosity force
      tend = tend + tend_vel_mep

   end subroutine ocn_mep_force_vert_visc!}}}


!***********************************************************************
!
!  routine ocn_mep_force_twa_drag_GM
!
!> \brief  This routine calculates normalVelocity tend due GM bolus coriolis.
!> \author Juan A. Saenz, Todd Ringler
!> \date   November, 2015
!> \details
!>  This routine calculates normalVelocity tend due coriolis force caused by
!>  the GM bolus velocity.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_force_twa_drag_GM(mepPool, tendPool, meshPool, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: mepPool  !< Input: Diagnostic information
      type (mpas_pool_type), intent(in) :: meshPool            !< Input: Mesh information

      !jas issue: it might be that only need uBolus or normalBolusGM,
      !not entire diagnostics?

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency structure
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend          !< Input/Output: normal velocity tendency

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

      logical  :: includeHalo
      integer :: iEdge, kLevel
      integer, pointer :: nVertLevels, nCells, nEdges

      integer, dimension(:,:), pointer :: edgeMask
      real (kIND=RKIND), dimension(:),   pointer :: fEdge
      real (kind=RKIND), dimension(:,:), pointer :: edgeTangentVectors
      real (kind=RKIND), dimension(:,:), pointer :: tend_vel_mep
      real (kind=RKIND), dimension(:,:), pointer :: GMBolusVelocityX
      real (kind=RKIND), dimension(:,:), pointer :: GMBolusVelocityY
      real (kind=RKIND), dimension(:,:), pointer :: GMBolusVelocityZ

      real (kind=RKIND), dimension(:,:),   allocatable :: normalBolusEdge
      real (kind=RKIND), dimension(:,:),   allocatable :: tangentBolusEdge
      real (kind=RKIND), dimension(:,:,:), allocatable :: bolusCell


      err = 0


      call mpas_pool_get_dimension(meshPool, 'nCells', nCells) ! includes halos
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges) ! includes halos
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'edgeTangentVectors', edgeTangentVectors)


      allocate(normalBolusEdge(nVertLevels,nEdges))
      allocate(tangentBolusEdge(nVertLevels,nEdges))
      allocate(bolusCell(3,nVertLevels,nCells+1))


      includeHalo = .true.

      ! retrieve velocity tendencry by MEP
      call mpas_pool_get_array(tendPool, 'normalVelocityMEP', tend_vel_mep)
      tend_vel_mep = 0.0


      !retrieve x,y,z bolus velocity at cell centers
      call mpas_pool_get_array(mepPool, 'GMBolusVelocityX', GMBolusVelocityX)
      call mpas_pool_get_array(mepPool, 'GMBolusVelocityY', GMBolusVelocityY)
      call mpas_pool_get_array(mepPool, 'GMBolusVelocityZ', GMBolusVelocityZ)

      bolusCell(1,:,:) = GMBolusVelocityX
      bolusCell(2,:,:) = GMBolusVelocityY
      bolusCell(3,:,:) = GMBolusVelocityZ

      !interpolate bolus velocity at cell centers, mid layers, to Edges
      !make sure to average the x,y,z components, not lat/lon components.
      call mpas_vector_R3Cell_to_2DEdge(bolusCell, meshPool, edgeTangentVectors, &
        includeHalo, normalBolusEdge, tangentBolusEdge)

      !calculate tendency = -dE3i/dz dot ni = - f*uBolus_tangent
      do iEdge = 1, nEdges
        do kLevel = 1,nVertLevels
          tend_vel_mep(kLevel,iEdge) = &
              -fEdge(iEdge) * tangentBolusEdge(kLevel,iEdge) * edgeMask(kLevel, iEdge)
        end do
      end do

      !update tendency with div EPFT parameterized with GM
      tend = tend + tend_vel_mep


      deallocate(normalBolusEdge)
      deallocate(tangentBolusEdge)
      deallocate(bolusCell)

    end subroutine ocn_mep_force_twa_drag_GM!}}}


!***********************************************************************
!
!  routine ocn_mep_vmix_coefs_build
!
!> \brief  This routine calculates the vertical viscocity from the mep.
!> \author Juan A. Saenz, Todd Ringler
!> \date   November, 2015
!> \details
!>  This routine calculates the vertical viscocity from the mep.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vmix_coefs_build(meshPool, statePool, diagnosticsPool, mepPool, err, timeLevel)!{{{
      ! This routine receives more arguments than it needs, because this is the standard list of
      ! arguments sent to other ocn_vmix_coefs_build routines in ocn_vmix_coefs(...)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool        !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: statePool        !< Input: state
      type (mpas_pool_type), intent(in) :: mepPool        !< Input/Output: mep stuff
      integer, intent(in), optional :: timeLevel         !< Input: Time level for state pool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: diagnosticsPool        !< Input/Output: diagnostics

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
      real (kind=RKIND), dimension(:,:), pointer :: vertViscTopOfEdge


      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfEdge', vertViscTopOfEdge)


      ! MEP vertical viscocity
      call ocn_mep_vmix_mom_build(vertViscTopOfEdge, mepPool, err)

      return

   end subroutine ocn_mep_vmix_coefs_build!}}}


!***********************************************************************
!
!  routine ocn_mep_vel_hmix_tend
!
!> \brief  This routine increments the horizontal mixing tendency for momentum due to mep.
!> \author Juan A. Saenz, Todd Ringler
!> \date   November, 2015
!> \details
!>  Placeholder routine, this doesn't do anything right now.
!>  This routine increments the horizontal mixing tendency for momentum due to
!>  meso scale eddies. 
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vel_hmix_tend(meshPool, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool        !< Input: Mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: tend  !< Input/Output: normal velocity tendency

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

      return

   end subroutine ocn_mep_vel_hmix_tend!}}}


!***********************************************************************
!
!  routine ocn_mep_reconstruct
!
!> \brief  This routine interfaces the cell center reconstruction
!> \author Juan A. Saenz, Todd Ringler
!> \date   17 September 2015
!> \details
!>  This routine interfaces the cell center reconstruction for diagnostics
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_reconstruct(diagnosticsPool, meshPool, mepPool)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostic information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: mepPool !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      ! NOTE: This routine is called within block loops

      if (mepUseGMmodule) then
         call ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool, mepPool)
      end if


   end subroutine ocn_mep_reconstruct!}}}


!***********************************************************************
!
!  routine ocn_mep_twa_drag_VertVisc_init
!
!> \brief  This routine interfaces the cell center reconstruction
!> \author Juan A. Saenz, Todd Ringler
!> \date   17 September 2015
!> \details
!>  This routine interfaces the cell center reconstruction for diagnostics
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_twa_drag_VertVisc_init(err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_mep_twa_drag_VertVisc_init!}}}


!***********************************************************************

end module ocn_mep

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
