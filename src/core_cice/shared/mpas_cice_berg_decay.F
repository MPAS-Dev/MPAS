!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_berg_decay
!
!> \brief   Calculate iceberg velocities
!> \author  Darin Comeau, LANL
!> \date    September 5, 2017
!> \details This module calculates iceberg decay by parameterizations
!> for basal melting, lateral melting, and wave erosion.
!> See Adcroft & Martin 2010 Ocean Modeling for details.
!
!-----------------------------------------------------------------------

module cice_berg_decay

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_dmpar
  use mpas_timer
  use mpas_log, only: mpas_log_write

  implicit none

  private
  save

  public :: &
       cice_run_berg_decay

  ! berg decay mask thresholds
  real(kind=RKIND), parameter, private :: &
       bergMassMinimum = 0.01_RKIND, &
       bergLengthMinimum = 0.01_RKIND, &
       bergHeightMinimum = 0.01_RKIND

contains

!-----------------------------------------------------------------------
! Time stepping
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_run_berg_decay
!
!> \brief   Calculate iceberg decay
!> \author  Darin Comeau, LANL
!> \date    September 5, 2017
!> \details This subroutine is a wrapper for the iceberg decay timestep.
!
!-----------------------------------------------------------------------

  subroutine cice_run_berg_decay(domain, clock)!{{{

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type(block_type), pointer :: &
         block

    logical, pointer :: &
         config_use_berg_decay, & ! flag to turn on / off berg decay parameterization
         config_berg_allow_rollover ! flag to turn on / off berg rollover due to instability

    ! determine if decay is switched on
    call MPAS_pool_get_config(domain % configs, "config_use_berg_decay", config_use_berg_decay)
    call MPAS_pool_get_config(domain % configs, "config_berg_allow_rollover", config_berg_allow_rollover)

    if (config_use_berg_decay) then

       block => domain % blocklist
       do while (associated(block))

          ! reset iceberg decay variables
          call init_berg_decay(block)

          ! get mask for iceberg decay calculations
          call berg_decay_calculation_mask(block)

          ! wrapper for calculating each of three iceberg decay terms
          call calculate_berg_decay_terms(block)

          ! apply iceberg decay terms
          call berg_calculate_melt(block)

          ! convert lost iceberg mass to freshwater and latent heat fluxes
          call berg_melt_fluxes(block)

          ! determine if iceberg is stable, and rollover if needed
          if (config_berg_allow_rollover) then
             call mpas_timer_start("Berg rollover check")
             call berg_rollover_check(block)
             call mpas_timer_stop("Berg rollover check")
          endif

          ! update iceberg area state variable
          call berg_update_area_after_decay(block)

          block => block % next
       enddo

    endif ! config_use_berg_decay

  end subroutine cice_run_berg_decay!}}}

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_decay
!
!> \brief   Initialize berg decay
!> \author  Darin Comeau, LANL
!> \date    September 5, 2017
!> \details This subroutine resets iceberg decay terms and fluxes from 
!> last timestep as 0.
!
!-----------------------------------------------------------------------

  subroutine init_berg_decay(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         bergDecayPool, &
         bergFluxesPool

    real(kind=RKIND), dimension (:,:), pointer :: &
         bergBasalMeltRate, &   ! berg decay rate due to basal melting (m/s)
         bergLateralMeltRate, & ! berg decay rate due to lateral melting (m/s)
         bergWaveErosionRate, & ! berg decay rate due to wave erosion (m/s)
         bergMeltMass, &        ! berg mass lost to decay in one timestep (kg)
         bergFreshwaterFluxCategory, & ! freshwater flux due to berg decay (kg/m^2/s)
         bergLatentHeatFluxCategory    ! latent heat flux due to berg decay (J/m^2/s)

    real(kind=RKIND), dimension (:), pointer :: &
         bergFreshwaterFluxCell, & ! aggregate freshwater flux (kg/m^2/s)
         bergLatentHeatFluxCell    ! aggregate latent heat flux (J/m^2/s)

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
    call MPAS_pool_get_subpool(block % structs, "berg_decay", bergDecayPool)
    call MPAS_pool_get_subpool(block % structs, "berg_fluxes", bergFluxesPool)

    call MPAS_pool_get_dimension(meshPool, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(meshPool, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(bergDecayPool, "bergBasalMeltRate", bergBasalMeltRate)
    call MPAS_pool_get_array(bergDecayPool, "bergLateralMeltRate", bergLateralMeltRate)
    call MPAS_pool_get_array(bergDecayPool, "bergWaveErosionRate", bergwaveErosionRate)
    call MPAS_pool_get_array(bergDecayPool, "bergMeltMass", bergMeltMass)

    call MPAS_pool_get_array(bergFluxesPool, "bergFreshwaterFluxCategory", bergFreshwaterFluxCategory)
    call MPAS_pool_get_array(bergFluxesPool, "bergLatentHeatFluxCategory", bergLatentHeatFluxCategory)
    call MPAS_pool_get_array(bergFluxesPool, "bergFreshwaterFluxCell", bergFreshwaterFluxCell)
    call MPAS_pool_get_array(bergFluxesPool, "bergLatentHeatFluxCell", bergLatentHeatFluxCell)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          ! reset berg decay rates & fluxes from previous timestep
          bergBasalMeltRate(iCategory,iCell) = 0.0_RKIND
          bergLateralMeltRate(iCategory,iCell) = 0.0_RKIND
          bergWaveErosionRate(iCategory,iCell) = 0.0_RKIND
          bergMeltMass(iCategory,iCell) = 0.0_RKIND

          bergFreshwaterFluxCategory(iCategory, iCell) = 0.0_RKIND
          bergLatentHeatFluxCategory(iCategory, iCell) = 0.0_RKIND

       enddo

       bergFreshwaterFluxCell(iCell) = 0.0_RKIND
       bergLatentHeatFluxCell(iCell) = 0.0_RKIND

    enddo

  end subroutine init_berg_decay

!-----------------------------------------------------------------------
! Calculation mask
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_decay_calculation_mask
!
!> \brief   Create mask for berg decay calculation
!> \author  Darin Comeau, LANL
!> \date    September 5, 2017
!> \details This subroutine creates a mask for calculating iceberg decay.
!
!-----------------------------------------------------------------------

  subroutine berg_decay_calculation_mask(block)

    type(block_type), intent(inout)  :: &
         block

    type(MPAS_pool_type), pointer :: &
         bergTracersPool, &
         bergStatePool, &
         bergDecayPool

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         bergLength, & ! berg state variable
         bergHeight    ! berg state variable

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory ! berg state variable

    integer, dimension(:,:), pointer :: &
         bergDecayMask ! mask = 1 indicates berg decay to occur, = 0 elsewhere

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer :: &
         iCategory, &
         iCell

    call MPAS_pool_get_subpool(block % structs, "berg_tracers", bergTracersPool)
    call MPAS_pool_get_subpool(block % structs, "berg_state", bergStatePool)
    call MPAS_pool_get_subpool(block % structs, "berg_decay", bergDecayPool)

    call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(bergTracersPool, "bergLength", bergLength, 1)
    call MPAS_pool_get_array(bergTracersPool, "bergHeight", bergHeight, 1)
    call MPAS_pool_get_array(bergStatePool, "bergMassCategory", bergMassCategory)
    call MPAS_pool_get_array(bergDecayPool, "bergDecayMask", bergDecayMask)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          bergDecayMask(iCategory,iCell) = 0

          if (bergMassCategory(iCategory,iCell) > bergMassMinimum .and. &
             bergLength(1,iCategory,iCell) > bergLengthMinimum .and. &
             bergHeight(1,iCategory,iCell) > bergHeightMinimum) then

             ! this cell has sufficient bergs for decay
             bergDecayMask(iCategory,iCell) = 1

          endif

       enddo

    enddo ! iCell

  end subroutine berg_decay_calculation_mask

!-----------------------------------------------------------------------
! Calculate decay terms
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calculate_berg_decay_terms
!
!> \brief   Calculate decay terms in berg decay equation
!> \author  Darin Comeau, LANL
!> \date    Sep 11, 2017
!> \details This subroutine is a wrapper for calculating basal melting,
!> lateral melting, and wave erosion terms in the iceberg decay equation.
!> See Adcroft & Martin 2010 Ocean Modeling for details.
!
!-----------------------------------------------------------------------

  subroutine calculate_berg_decay_terms(block)

    type(block_type), intent(inout) :: &
         block

    ! calculate basal melting
    call mpas_timer_start("berg basal melting")
    call berg_basal_melting(block)
    call mpas_timer_stop("berg basal melting")

    ! calculate lateral melting
    call mpas_timer_start("berg lateral melting")
    call berg_lateral_melting(block)
    call mpas_timer_stop("berg lateral melting")

    ! calculate decay due to wave erosion
    call mpas_timer_start("berg wave erosion")
    call berg_wave_erosion(block)
    call mpas_timer_stop("berg wave erosion")

  end subroutine calculate_berg_decay_terms!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_relative_speed
!
!> \brief   Get relative speed between two fields
!> \author  Darin Comeau, LANL
!> \date    Sep 12, 2017
!> \details This subroutine calculates the relative speed between two
!> fields at a point.
!
!-----------------------------------------------------------------------

  subroutine get_relative_speed(&
       uVelocity1, &
       vVelocity1, &
       uVelocity2, &
       vVelocity2, &
       relSpeed)

    real(kind=RKIND), intent(in) :: &
         uVelocity1, &
         vVelocity1, &
         uVelocity2, &
         vVelocity2

    real(kind=RKIND), intent(out) :: &
         relSpeed

    relSpeed = sqrt((uVelocity1 - uVelocity2)**2 + (vVelocity1 - vVelocity2)**2)

  end subroutine get_relative_speed

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_basal_melting
!
!> \brief   Calculate berg basal melting
!> \author  Darin Comeau, LANL
!> \date    Sep 11, 2017
!> \details This subroutine calculates berg decay due to basal melting,
!> in units of m/s.
!
!-----------------------------------------------------------------------

  subroutine berg_basal_melting(block)

    use cice_constants, only: &
         ciceSecondsPerDay

    use cice_mesh, only: &
         cice_interpolate_vertex_to_cell

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         boundaryPool, &
         bergDecayPool, &
         bergVelocitySolverPool, &
         bergTracersPool, &
         oceanCouplingPool

    integer, dimension(:,:), pointer :: &
         bergDecayMask ! mask = 1 indicates berg decay to occur, = 0 elsewhere

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergBasalMeltRate, & ! decay rate due to basal melting (m/s)
         uBergVelocityCell, & ! berg velocity interpolated to cell centers
         vBergVelocityCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         bergLength           ! berg length dimension (m)

    real(kind=RKIND), dimension(:,:), pointer :: &
         uBergVelocity, & ! berg velocity at vertices
         vBergVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         uOceanVelocity, &
         vOceanVelocity, &
         seaSurfaceTemperature

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    real, pointer :: &
         bergTemperature ! berg temperature, assumed constant

    integer :: &
         iCell, &
         iCategory

    real(kind=RKIND) :: &
         relSpeed

    real(kind=RKIND), parameter :: &
         basalMeltCoef = 0.58_RKIND, & ! empirical coefficients in decay parameterization
         basalMeltExponent1 = 0.8_RKIND, &
         basalMeltExponent2 = 0.2_RKIND

    call MPAS_pool_get_config(block % configs, "config_berg_temperature", bergTemperature)

    call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

    call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundaryPool)
    call MPAS_pool_get_subpool(block % structs, "berg_decay", bergDecayPool)
    call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)
    call MPAS_pool_get_subpool(block % structs, "berg_tracers", bergTracersPool)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)

    call MPAS_pool_get_array(bergDecayPool, "bergDecayMask", bergDecayMask)
    call MPAS_pool_get_array(bergDecayPool, "bergBasalMeltRate", bergBasalMeltRate)

    call MPAS_pool_get_array(bergVelocitySolverPool, "uBergVelocity", uBergVelocity)
    call MPAS_pool_get_array(bergVelocitySolverPool, "vBergVelocity", vBergVelocity)
    call MPAS_pool_get_array(bergVelocitySolverPool, "uBergVelocityCell", uBergVelocityCell)
    call MPAS_pool_get_array(bergVelocitySolverPool, "vBergVelocityCell", vBergVelocityCell)

    call MPAS_pool_get_array(bergTracersPool, "bergLength", bergLength, 1)

    call MPAS_pool_get_array(oceanCouplingPool, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(oceanCouplingPool, "vOceanVelocity", vOceanVelocity)
    call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTemperature", seaSurfaceTemperature)

    ! interpolate berg velocities to cell centers
    do iCategory = 1, nBergCategories

       call cice_interpolate_vertex_to_cell(&
            meshPool, &
            boundaryPool, &
            uBergVelocityCell(iCategory,:), &
            uBergVelocity(iCategory,:))

       call cice_interpolate_vertex_to_cell(&
            meshPool, &
            boundaryPool, &
            vBergVelocityCell(iCategory,:), &
            vBergVelocity(iCategory,:))

    enddo

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          if (bergDecayMask(iCategory,iCell) == 1) then

             ! get relative speed between ocean and bergs
             call get_relative_speed(uOceanVelocity(iCell), &
                                     vOceanVelocity(iCell), &
                                     uBergVelocityCell(iCategory,iCell), &
                                     vBergVelocityCell(iCategory,iCell), &
                                     relSpeed)

             ! calculate berg decay (m/day)
             bergBasalMeltRate(iCategory,iCell) = basalMeltCoef*relSpeed**basalMeltExponent1 &
                                                * (seaSurfaceTemperature(iCell) - bergTemperature) &
                                                / (bergLength(1,iCategory,iCell)**basalMeltExponent2)

             ! convert from m/day to m/sec
             bergBasalMeltRate(iCategory,iCell) = bergBasalMeltRate(iCategory,iCell)/real(ciceSecondsPerDay,RKIND)

             ! melt term must be positive
             bergBasalMeltRate(iCategory,iCell) = max(bergBasalMeltRate(iCategory,iCell), 0.0_RKIND)

          endif

       enddo

    enddo

  end subroutine berg_basal_melting

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_lateral_melting
!
!> \brief   Calculate berg lateral melting
!> \author  Darin Comeau, LANL
!> \date    Sep 11, 2017
!> \details This subroutine calculates berg decay due to lateral melting,
!> in units of m/s.
!
!-----------------------------------------------------------------------

  subroutine berg_lateral_melting(block)

    use cice_constants, only: &
         ciceSecondsPerDay

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         bergDecayPool, &
         oceanCouplingPool

    integer, dimension(:,:), pointer :: &
         bergDecayMask ! mask = 1 indicates berg decay to occur, = 0 elsewhere

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergLateralMeltRate ! decay rate due to lateral melting (m/s)

    real(kind=RKIND), dimension(:), pointer :: &
         seaSurfaceTemperature

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer :: &
         iCell, &
         iCategory

    real(kind=RKIND) :: &
         relSpeed

    real(kind=RKIND), parameter :: &
         lateralMeltCoef1 = 0.00762_RKIND, & ! empirical coefficients in expression
         lateralMeltCoef2 = 0.00129_RKIND    ! for lateral melting

    call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)
      
    call MPAS_pool_get_subpool(block % structs, "berg_decay", bergDecayPool)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)

    call MPAS_pool_get_array(bergDecayPool, "bergDecayMask", bergDecayMask)
    call MPAS_pool_get_array(bergDecayPool, "bergLateralMeltRate", bergLateralMeltRate)
     
    call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTemperature", seaSurfaceTemperature)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          if (bergDecayMask(iCategory,iCell) == 1) then

             ! calculate decay rate due to lateral melting (m/day)
             bergLateralMeltRate(iCategory,iCell) = lateralMeltCoef1*seaSurfaceTemperature(iCell) + &
                                                    lateralMeltCoef2*seaSurfaceTemperature(iCell)**2

             ! convert from m/day to m/sec
             bergLateralMeltRate(iCategory,iCell) = bergLateralMeltRate(iCategory,iCell)/real(ciceSecondsPerDay,RKIND)

             ! melt term must be positive
             bergLateralMeltRate(iCategory,iCell) = max(bergLateralMeltRate(iCategory,iCell), 0.0_RKIND)

          endif

       enddo

    enddo

  end subroutine berg_lateral_melting

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_wave_erosion
!
!> \brief   Calculate berg wave erosion
!> \author  Darin Comeau, LANL
!> \date    Sep 11, 2017
!> \details This subroutine calculates berg decay due to wave erosion,
!> in units of m/s.
!
!-----------------------------------------------------------------------

  subroutine berg_wave_erosion(block)

    use cice_constants, only: &
         ciceSecondsPerDay, &
         pii

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         bergDecayPool, &
         tracersAggregatePool, & 
         oceanCouplingPool, &
         atmosCouplingPool

    integer, dimension(:,:), pointer :: &
         bergDecayMask ! mask = 1 indicates berg decay to occur, = 0 elsewhere

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergWaveErosionRate ! decay rate due to wave erosion (m/s)

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         seaSurfaceTemperature, &
         uOceanVelocity, &
         vOceanVelocity, &
         uAirVelocity, &
         vAirVelocity

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer :: &
         iCell, &
         iCategory

    real(kind=RKIND) :: &
         relSpeed, &
         seaState ! estimate sea state, fit to Beaufort scale

    real(kind=RKIND), parameter :: &
         waveCoef1 = 1.0_RKIND/12.0_RKIND, & ! empirical coefficients in expression
         waveCoef2 = 3.0_RKIND/2.0_RKIND, &  ! for decay due to wave erosion
         waveCoef3 = 1.0_RKIND/10.0_RKIND

    call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

    call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
    call MPAS_pool_get_subpool(block % structs, "berg_decay", bergDecayPool)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCouplingPool)

    call MPAS_pool_get_array(bergDecayPool, "bergDecayMask", bergDecayMask)
    call MPAS_pool_get_array(bergDecayPool, "bergWaveErosionRate", bergWaveErosionRate)

    call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)

    call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTemperature", seaSurfaceTemperature)
    call MPAS_pool_get_array(oceanCouplingPool, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(oceanCouplingPool, "vOceanVelocity", vOceanVelocity)

    call MPAS_pool_get_array(atmosCouplingPool, "uAirVelocity", uAirVelocity)
    call MPAS_pool_get_array(atmosCouplingPool, "vAirVelocity", vAirVelocity)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          if (bergDecayMask(iCategory,iCell) == 1) then

             ! get relative speed of air and bergs
             call get_relative_speed(uAirVelocity(iCell), &
                                     vAirVelocity(iCell), &
                                     uOceanVelocity(iCell), &
                                     vOceanVelocity(iCell), &
                                     relSpeed)

             ! estimate sea state fit to Beaufort scale
             seaState = waveCoef2*sqrt(relSpeed) + waveCoef3*relSpeed

             ! calculate wave erosion rate (m/day)
             bergWaveErosionRate(iCategory,iCell) = waveCoef1*seaState &
                                                  *(1.0_RKIND + cos(iceAreaCell(iCell)**3*pii)) &
                                                  *(seaSurfaceTemperature(iCell) + 2.0_RKIND)

             ! convert from m/day to m/sec
             bergWaveErosionRate(iCategory,iCell) = bergWaveErosionRate(iCategory,iCell)/real(ciceSecondsPerDay,RKIND)

             ! melt term must be positive
             bergWaveErosionRate(iCategory,iCell) = max(bergWaveErosionRate(iCategory,iCell), 0.0_RKIND)

          endif

       enddo

    enddo

  end subroutine berg_wave_erosion

!-----------------------------------------------------------------------
! Apply decay terms
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_calculate_melt
!
!> \brief   Calculate berg melt.
!> \author  Darin Comeau, LANL
!> \date    Sep 11, 2017
!> \details This subroutine calculates change in berg mass due to
!> decay in a single timestep.
!
!-----------------------------------------------------------------------

  subroutine berg_calculate_melt(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         bergDecayPool, &
         bergStatePool, &
         bergTracersPool

    integer, dimension(:,:), pointer :: &
         bergDecayMask ! mask = 1 indicates berg decay to occur, = 0 elsewhere

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory, &    ! berg mass state variable
         bergBasalMeltRate, &   ! decay due to basal melting (m/s)
         bergLateralMeltRate, & ! decay due to lateral melting (m/s)
         bergWaveErosionRate, & ! decay due to wave erosion (m/s)
         bergMeltMass           ! berg mass lost to decay in one timestep (kg)

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         bergLength, &       ! berg category length dimension
         bergHeight          ! berg category height dimension

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    real(kind=RKIND), pointer :: &
         config_dt

    integer :: &
         iCell, &
         iCategory

    real(kind=RKIND) :: &
         dBergLength, & ! change in berg length due to decay (m)
         dBergHeight, & ! change in berg height due to decay (m)
         dBergMass      ! fraction of berg mass lost

    call MPAS_pool_get_config(block % configs, "config_dt", config_dt)

    call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

    call MPAS_pool_get_subpool(block % structs, "berg_decay", bergDecayPool)
    call MPAS_pool_get_subpool(block % structs, "berg_state", bergStatePool)
    call MPAS_pool_get_subpool(block % structs, "berg_tracers", bergTracersPool)

    call MPAS_pool_get_array(bergDecayPool, "bergDecayMask", bergDecayMask)
    call MPAS_pool_get_array(bergDecayPool, "bergBasalMeltRate", bergBasalMeltRate)
    call MPAS_pool_get_array(bergDecayPool, "bergLateralMeltRate", bergLateralMeltRate)
    call MPAS_pool_get_array(bergDecayPool, "bergWaveErosionRate", bergWaveErosionRate)
    call MPAS_pool_get_array(bergDecayPool, "bergMeltMass", bergMeltMass)
    call MPAS_pool_get_array(bergStatePool, "bergMassCategory", bergMassCategory)
    call MPAS_pool_get_array(bergTracersPool, "bergLength", bergLength, 1)
    call MPAS_pool_get_array(bergTracersPool, "bergHeight", bergHeight, 1)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          if (bergDecayMask(iCategory,iCell) == 1) then

             ! calculate change in dimensions
             dBergLength = config_dt*1.25_RKIND*(bergLateralMeltRate(iCategory,iCell) + bergWaveErosionRate(iCategory,iCell))
             dBergHeight = config_dt*bergBasalMeltRate(iCategory,iCell)

             ! terms must be positive
             if (dBergLength < 0.0_RKIND) call mpas_log_write(&
                                               "berg_decay: negative berg length $r:", realArgs=(/dBergLength/), &
                                               messageType=MPAS_LOG_CRIT)
             if (dBergHeight < 0.0_RKIND) call mpas_log_write(&
                                               "berg_decay: negative berg height $r:", realArgs=(/dBergHeight/), &
                                               messageType=MPAS_LOG_CRIT)

             ! fraction of mass lost
             dBergMass = (dBergLength**2*dBergHeight) / (bergLength(1,iCategory,iCell)**2*bergHeight(1,iCategory,iCell))

             ! update berg dimensions
             bergLength(1,iCategory,iCell) = bergLength(1,iCategory,iCell) - dBergLength
             bergHeight(1,iCategory,iCell) = bergHeight(1,iCategory,iCell) - dBergHeight
  
             ! calculate berg mass lost (kg)
             bergMeltMass(iCategory,iCell) = dBergMass*bergMassCategory(iCategory,iCell)

             ! update berg mass state variable
             bergMassCategory(iCategory,iCell) = bergMassCategory(iCategory,iCell) - bergMeltMass(iCategory,iCell)

          endif

       enddo

    enddo

  end subroutine berg_calculate_melt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_melt_fluxes
!
!> \brief   Calculate berg melt fluxes
!> \author  Darin Comeau, LANL
!> \date    Sep 14, 2017
!> \details Convert berg melt to freshwater and latent heat fluxes
!> intended for the ocean.
!
!-----------------------------------------------------------------------

  subroutine berg_melt_fluxes(block)

    use cice_constants, only: &
         ciceLatentHeatMelting

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         bergDecayPool, &
         bergFluxesPool

    integer, dimension(:,:), pointer :: &
         bergDecayMask ! mask = 1 indicates berg decay to occur, = 0 elsewhere

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMeltMass, & ! berg mass lost to decay in one timestep (kg)
         bergFreshwaterFluxCategory, & ! freshwater flux due to berg decay (kg/m^2/s)
         bergLatentHeatFluxCategory    ! latent heat flux due to berg decay (J/m^2/s)

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell, &
         bergFreshwaterFluxCell, & ! aggregate freshwater flux (kg/m^2/s)
         bergLatentHeatFluxCell    ! aggregate latent heat flux (J/m^2/s)

    real(kind=RKIND), pointer :: &
         config_dt, &
         bergTemperature, &        ! berg temperature, assumed constant
         specificHeatFreshIce

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
    call MPAS_pool_get_config(block % configs, "config_berg_temperature", bergTemperature)
    call MPAS_pool_get_config(block % configs, "config_specific_heat_fresh_ice", specificHeatFreshIce)

    call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

    call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
    call MPAS_pool_get_subpool(block % structs, "berg_decay", bergDecayPool)
    call MPAS_pool_get_subpool(block % structs, "berg_fluxes", bergFluxesPool)

    call MPAS_pool_get_array(meshPool, "areaCell", areaCell)

    call MPAS_pool_get_array(bergDecayPool, "bergDecayMask", bergDecayMask)
    call MPAS_pool_get_array(bergDecayPool, "bergMeltMass", bergMeltMass)

    call MPAS_pool_get_array(bergFluxesPool, "bergFreshwaterFluxCategory", bergFreshwaterFluxCategory)
    call MPAS_pool_get_array(bergFluxesPool, "bergLatentHeatFluxCategory", bergLatentHeatFluxCategory)
    call MPAS_pool_get_array(bergFluxesPool, "bergFreshwaterFluxCell", bergFreshwaterFluxCell)
    call MPAS_pool_get_array(bergFluxesPool, "bergLatentHeatFluxCell", bergLatentHeatFluxCell)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          if (bergDecayMask(iCategory,iCell) == 1) then

             ! calculate freshwater flux (kg/m^2/s)
             bergFreshwaterFluxCategory(iCategory,iCell) = bergMeltMass(iCategory,iCell) / (areaCell(iCell)*config_dt)

             ! calculate latent heat flux (J/m^2/s)
             bergLatentHeatFluxCategory(iCategory,iCell) = bergMeltMass(iCategory,iCell) &
                                                         * (ciceLatentHeatMelting - specificHeatFreshIce*bergTemperature) &
                                                         / (areaCell(iCell)*config_dt)

             ! aggregate quantities
             bergFreshwaterFluxCell(iCell) = bergFreshwaterFluxCell(iCell) + bergFreshwaterFluxCategory(iCategory,iCell)
             bergLatentHeatFluxCell(iCell) = bergLatentHeatFluxCell(iCell) + bergLatentHeatFluxCategory(iCategory,iCell)

          endif

       enddo

    enddo

  end subroutine berg_melt_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_stability_check
!
!> \brief   Check iceberg stability
!> \author  Darin Comeau, LANL
!> \date    Sep 18, 2017
!> \details This subroutine checks the Weeks-Mellor stability criterion
!> to see if the iceberg is unstable, indicating rollover.
!
!-----------------------------------------------------------------------

  subroutine berg_stability_check(&
       bergLength, &
       bergHeight, &
       bergStabilityCheck)

    real(kind=RKIND), intent(in) :: &
         bergLength, &
         bergHeight

    integer, intent(out) :: &
         bergStabilityCheck ! 1 indicates rollover should occur

    real(kind=RKIND), parameter :: &
         rolloverCoef1 = 0.92_RKIND, & ! empirical coefficients in the
         rolloverCoef2 = 58.32_RKIND   ! Weeks-Mellor stability criterion

    bergStabilityCheck = 0

    ! check Weeks-Mellor stability criterion
    if (bergLength < sqrt(rolloverCoef1*bergHeight**2 + rolloverCoef2*bergHeight)) &
       bergStabilityCheck = 1

  end subroutine berg_stability_check

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_rollover_check
!
!> \brief   Check if berg rollover needed
!> \author  Darin Comeau, LANL
!> \date    Sep 18, 2017
!> \details Check if berg dimensions violate Weeks-Mellor stability
!> criterion, and need to rollover.
!
!-----------------------------------------------------------------------

  subroutine berg_rollover_check(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         bergDecayPool, &
         bergTracersPool

    integer, dimension(:,:), pointer :: &
         bergDecayMask ! mask = 1 indicates berg decay to occur, = 0 elsewhere

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         bergLength, & ! berg category length
         bergHeight    ! berg category height

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer :: &
         iCell, &
         iCategory

    integer :: &
         bergStabilityCheck ! 1 indicates berg is unstable and rollover should occur

    real(kind=RKIND) :: &
         bergLengthOld, & ! placeholders during rollover when length and height
         bergHeightOld    ! are switched

    call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

    call MPAS_pool_get_subpool(block % structs, "berg_decay", bergDecayPool)
    call MPAS_pool_get_subpool(block % structs, "berg_tracers", bergTracersPool)

    call MPAS_pool_get_array(bergDecayPool, "bergDecayMask", bergDecayMask)

    call MPAS_pool_get_array(bergTracersPool, "bergLength", bergLength, 1)
    call MPAS_pool_get_array(bergTracersPool, "bergHeight", bergHeight, 1)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          if (bergDecayMask(iCategory,iCell) == 1) then

             ! check if berg is unstable
             call berg_stability_check(bergLength(1,iCategory,iCell), &
                                       bergHeight(1,iCategory,iCell), &
                                       bergStabilityCheck)

            ! berg is unstable, needs to rollover
             if (bergStabilityCheck == 1) then

                bergLengthOld = bergLength(1,iCategory,iCell)
                bergHeightOld = bergHeight(1,iCategory,iCell)

                ! switch height and width in such a way as to conserve volume
                bergHeight(1,iCategory,iCell) = (2.0_RKIND/3.0_RKIND)*bergLengthOld
                bergLength(1,iCategory,iCell) = sqrt((3.0_RKIND/2.0_RKIND)*bergLengthOld*bergHeightOld)

             endif

          endif

       enddo

    enddo

  end subroutine berg_rollover_check

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_update_area_after_decay
!
!> \brief   Update berg area after decay
!> \author  Darin Comeau, LANL
!> \date    September 26, 2017
!> \details This subroutine updates the iceberg area state variable
!> after the iceberg mass state variable has changed due to decay.
!
!-----------------------------------------------------------------------

  subroutine berg_update_area_after_decay(block)!{{{

    use cice_berg_state, only: cice_berg_area_update

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         bergTracersPool, &
         bergTracersAggregatePool, &
         bergStatePool, &
         bergDecayPool

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory    ! berg mass state variable

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         bergAreaCategory, & ! berg area state variable
         bergHeight          ! berg category height

    real(kind=RKIND), dimension(:), pointer :: &
         bergAreaCell, &     ! aggregate berg area
         areaCell

    integer, dimension(:,:), pointer :: &
         bergDecayMask ! mask = 1 indicates berg decay occured, = 0 elsewhere

    real(kind=RKIND), pointer :: &
         bergDensity   ! berg density (kg/m^3)

    integer :: &
         iCell, &
         iCategory

    logical :: &
         updateBergArea ! flag to keep track of when bergArea variables need updating

    call MPAS_pool_get_config(block % configs, "config_berg_density", bergDensity)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "berg_tracers", bergTracersPool)
    call MPAS_pool_get_subpool(block % structs, "berg_tracers_aggregate", bergTracersAggregatePool)
    call MPAS_pool_get_subpool(block % structs, "berg_state", bergStatePool)
    call MPAS_pool_get_subpool(block % structs, "berg_decay", bergDecayPool)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)
    call MPAS_pool_get_array(bergStatePool, "bergMassCategory", bergMassCategory)
    call MPAS_pool_get_array(bergTracersPool, "bergAreaCategory", bergAreaCategory,1)
    call MPAS_pool_get_array(bergTracersPool, "bergHeight", bergHeight,1)
    call MPAS_pool_get_array(bergTracersAggregatePool, "bergAreaCell", bergAreaCell)
    call MPAS_pool_get_array(bergDecayPool, "bergDecayMask", bergDecayMask)

    do iCell = 1, nCellsSolve

       updateBergArea = .false.

       do iCategory = 1, nBergCategories

          ! if decay occured, berg area state variables need to be updated
          if (bergDecayMask(iCategory,iCell) == 1) updateBergArea = .true.

       enddo

       ! update iceberg area if any category changed due to decay
       if (updateBergArea) then
          call cice_berg_area_update(nBergCategories, &
                                     bergMassCategory(:,iCell), &
                                     bergHeight(1,:,iCell), &
                                     bergDensity, &
                                     areaCell(iCell), &
                                     bergAreaCategory(1,:,iCell), &
                                     bergAreaCell(iCell))

       endif

    enddo

  end subroutine berg_update_area_after_decay!}}}

!-----------------------------------------------------------------------

end module cice_berg_decay
