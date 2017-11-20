!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_berg_state
!
!> \brief   Update iceberg model
!> \author  Darin Comeau, LANL
!> \date    26 Sept 2017
!> \details This module is intended to update the iceberg state variables.
!
!-----------------------------------------------------------------------

module cice_berg_state

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_log, only: mpas_log_write

   use mpas_stream_manager
   use mpas_io_units

   implicit none

   private
   save

   public :: &
      cice_berg_accumulate_calving_mass, &
      cice_berg_release_calving_mass, &
      cice_berg_dimensions_update, &
      cice_berg_area_update

contains

!-----------------------------------------------------------------------
! Calving subroutines
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_berg_accumulate_calving_mass
!
!> \brief   Accumulate calving mass
!> \author  Darin Comeau, LANL
!> \date    3 Oct 2017
!> \details This routine is intended to accumulate iceberg calving mass
!> from calving fluxes over a single timestep and add them to a stored
!> calving mass variable.
!
!-----------------------------------------------------------------------

  subroutine cice_berg_accumulate_calving_mass(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         berg_calving

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer, dimension(:,:), pointer :: &
         bergCalvingMask ! mask = 1 indicates calving cell, = 0 elsewhere

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergCalvingMass ! storage variable to collect mass from calving fluxes (kg)

    real(kind=RKIND), dimension(:), pointer :: &
         bergCalvingDistribution ! berg category calving distribution (dimensionless)

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell, &
         bergCalvingFlux ! calving flux (kg/m^2/s)

    real(kind=RKIND), pointer :: &
         config_dt, &
         bergDensity ! berg density (kg/m^3)

    integer :: &
         iCell, &
         iCategory

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_berg_density", bergDensity)

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

       call MPAS_pool_get_array(mesh, "areaCell", areaCell)
       call MPAS_pool_get_array(berg_calving, "bergCalvingMask", bergCalvingMask)
       call MPAS_pool_get_array(berg_calving, "bergCalvingMass", bergCalvingMass)
       call MPAS_pool_get_array(berg_calving, "bergCalvingFlux", bergCalvingFlux)
       call MPAS_pool_get_array(berg_calving, "bergCalvingDistribution", bergCalvingDistribution)

       do iCell = 1, nCellsSolve

          do iCategory = 1, nBergCategories

             if (bergCalvingMask(iCategory,iCell) == 1) then

                ! add mass based on calving flux to each category's storage variable
                bergCalvingMass(iCategory,iCell) = bergCalvingMass(iCategory,iCell) + &
                                                   bergCalvingDistribution(iCategory) * &
                                                   bergCalvingFlux(iCell)*areaCell(iCell)*config_dt

             endif

          enddo

       enddo

       block => block % next
    enddo

  end subroutine cice_berg_accumulate_calving_mass

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_berg_release_calving_mass
!
!> \brief   Release calving fluxes
!> \author  Darin Comeau, LANL
!> \date    3 Oct 2017
!> \details This routine is intended to release iceberg calving masses
!> accumulated from calving fluxes.
!
!-----------------------------------------------------------------------

  subroutine cice_berg_release_calving_mass(domain, itimestep)

    use cice_constants, only: &
         ciceSecondsPerYear

    type(domain_type), intent(inout) :: &
         domain

    integer, intent(in) :: &
         itimestep

    type(block_type), pointer :: &
         block

    character(len=strKIND), pointer :: &
         config_berg_calving_type ! indicates how calving mass is released to model

    real(kind=RKIND), pointer :: &
         dt, &
         annualCalvingFrequency ! number of calving events per year

    integer :: &
         calvingFlag ! flag to indicate if calving is to be done

    call MPAS_pool_get_config(domain % configs, "config_berg_calving_type", config_berg_calving_type)

    ! for continuous calving, release every time step
    if (trim(config_berg_calving_type) == "continuous") then

       block => domain % blocklist
       do while (associated(block))

          call berg_calving_event(block)

       block => block % next
       enddo

    ! for episodic calving, check if this timestep is a calving event
    elseif (trim(config_berg_calving_type) == "episodic") then

       call MPAS_pool_get_config(domain % configs, "config_dt", dt)
       call MPAS_pool_get_config(domain % configs, "config_berg_annual_calving_frequency", annualCalvingFrequency)

       calvingFlag = NINT(ciceSecondsPerYear/(dt*annualCalvingFrequency))

       ! check if calving even is to occur
       if (mod(itimestep,calvingFlag)==0) then

          call mpas_log_write("episodic berg calving event")

          block => domain % blocklist
          do while (associated(block))

             call berg_calving_event(block)

             block => block % next
          enddo

       endif

    ! for calving based on mass threshold, check if calving mass exceeds a single berg's mass
    elseif (trim(config_berg_calving_type) == "mass_threshold") then

       block => domain % blocklist
       do while (associated(block))

          call berg_calving_mass_threshold_check(block)
          call berg_calving_event(block)

       block => block % next
       enddo

    endif

  end subroutine cice_berg_release_calving_mass

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_calving_event
!
!> \brief   Iceberg calving event
!> \author  Darin Comeau, LANL
!> \date    26 Sept 2017
!> \details This routine is intended to add iceberg mass through calving,
!> and reset the variable bergCalvingMass to 0.
!
!-----------------------------------------------------------------------

  subroutine berg_calving_event(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         berg_tracers, &
         berg_tracers_aggregate, &
         berg_state, &
         berg_calving

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory, & ! berg mass state variable (kg)
         bergCalvingMass     ! storage variable for calving mass (kg)

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         bergAreaCategory, & ! berg area state variable (dimensionless)
         bergLength, &       ! berg category length (m)
         bergHeight          ! berg category height (m)

    integer, dimension(:,:), pointer :: &
         bergCalvingMask ! mask = 1 indicates calving cell, = 0 elsewhere

    real(kind=RKIND), dimension(:), pointer :: &
         bergMassCell, & ! aggregate berg mass state variable (kg)
         bergAreaCell, & ! aggregate berg area state variable (dimensionless)
         areaCell

    real(kind=RKIND), dimension(:), pointer :: &
         bergCalvingLength, & ! category length of berg when it calves (m)
         bergCalvingHeight    ! category height of berg when it calves (m)

    real(kind=RKIND), pointer :: &
         bergDensity ! berg density (kg/m^3)

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_config(block % configs, "config_berg_density", bergDensity)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "berg_tracers", berg_tracers)
    call MPAS_pool_get_subpool(block % structs, "berg_tracers_aggregate", berg_tracers_aggregate)
    call MPAS_pool_get_subpool(block % structs, "berg_state", berg_state)
    call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)
    call MPAS_pool_get_array(berg_state, "bergMassCategory", bergMassCategory)
    call MPAS_pool_get_array(berg_state, "bergMassCell", bergMassCell)
    call MPAS_pool_get_array(berg_tracers, "bergAreaCategory", bergAreaCategory,1)
    call MPAS_pool_get_array(berg_tracers, "bergLength", bergLength,1)
    call MPAS_pool_get_array(berg_tracers, "bergHeight", bergHeight,1)
    call MPAS_pool_get_array(berg_tracers_aggregate, "bergAreaCell", bergAreaCell)
    call MPAS_pool_get_array(berg_calving, "bergCalvingMask", bergCalvingMask)
    call MPAS_pool_get_array(berg_calving, "bergCalvingMass", bergCalvingMass)
    call MPAS_pool_get_array(berg_calving, "bergCalvingLength", bergCalvingLength)
    call MPAS_pool_get_array(berg_calving, "bergCalvingHeight", bergCalvingHeight)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          ! calving occurs
          if (bergCalvingMask(iCategory,iCell) == 1) then

             ! average dimensions with existing bergs
             call cice_berg_dimensions_update(bergMassCategory(iCategory,iCell), &
                                              bergCalvingMass(iCategory,iCell), &
                                              bergLength(1,iCategory,iCell), &
                                              bergCalvingLength(iCategory), &
                                              bergHeight(1,iCategory,iCell), &
                                              bergCalvingHeight(iCategory), &
                                              bergLength(1,iCategory,iCell), &
                                              bergHeight(1,iCategory,iCell))

             ! add calving mass
             bergMassCategory(iCategory,iCell) = bergMassCategory(iCategory,iCell) &
                                               + bergCalvingMass(iCategory,iCell)

             ! reset calving mass
             bergCalvingMass(iCategory,iCell) = 0.0_RKIND

          endif

       enddo

       ! integrated quantities
       bergMassCell(iCell) = sum(bergMassCategory(:,iCell))

       ! update iceberg area if any new bergs calved
       if (maxval(bergCalvingMask(:,iCell)) > 0) then

          call cice_berg_area_update(nBergCategories, &
                                     bergMassCategory(:,iCell), &
                                     bergHeight(1,:,iCell), &
                                     bergDensity, &
                                     areaCell(iCell), &
                                     bergAreaCategory(1,:,iCell), &
                                     bergAreaCell(iCell))

       endif

    enddo

  end subroutine berg_calving_event

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_calving_mass_threshold_check
!
!> \brief   Check calving mass threshold
!> \author  Darin Comeau, LANL
!> \date    3 Sept 2017
!> \details This routine is to check if a calving mass exceeds the
!> threshold of a single berg and therefore should be calved.
!
!-----------------------------------------------------------------------

  subroutine berg_calving_mass_threshold_check(block)

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type(MPAS_pool_type), pointer :: &
         mesh, &
         berg_calving

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer, dimension(:,:), pointer :: &
         bergCalvingMask ! mask = 1 indicates calving cell, = 0 elsewhere

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergCalvingMass ! storage variable for calving mass (kg)

    real(kind=RKIND), dimension(:), pointer :: &
         oneCalvedBergMass ! mass of one newly calved berg (kg)

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(berg_calving, "bergCalvingMask", bergCalvingMask)
    call MPAS_pool_get_array(berg_calving, "oneCalvedBergMass", oneCalvedBergMass)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          ! check if mass in storage variable exceeds that of a single berg
          if (bergCalvingMass(iCategory,iCell) > oneCalvedBergMass(iCategory)) then

             bergCalvingMask(iCategory,iCell) = 1

          else

             bergCalvingMask(iCategory,iCell) = 0

          endif

       enddo

    enddo

  end subroutine berg_calving_mass_threshold_check

!-----------------------------------------------------------------------
! Iceberg state variable updates
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_berg_dimensions_update
!
!> \brief   Update iceberg dimensions
!> \author  Darin Comeau, LANL
!> \date    Sep 26, 2017
!> \details This subroutine updates iceberg dimensions as the weighted
!> average of two iceberg masses in the same cell.
!
!-----------------------------------------------------------------------

  subroutine cice_berg_dimensions_update(&
       bergMass1, &
       bergMass2, &
       bergLength1, &
       bergLength2, &
       bergHeight1, &
       bergHeight2, &
       newBergLength, &
       newBergHeight)

    real(kind=RKIND), intent(in) :: &
         bergMass1, &
         bergMass2, &
         bergLength1, &
         bergLength2, &
         bergHeight1, &
         bergHeight2

    real(kind=RKIND), intent(out) :: &
         newBergLength, &
         newBergHeight

    real(kind=RKIND) :: &
         weight

    weight = bergMass1 / (bergMass1 + bergMass2)

    ! update dimensions as weighted sum of existing two dimensions
    newBergLength = weight*bergLength1 + (1.0_RKIND-weight)*bergLength2
    newBergHeight = weight*bergHeight1 + (1.0_RKIND-weight)*bergHeight2

  end subroutine cice_berg_dimensions_update

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_berg_area_update
!
!> \brief   Update iceberg area
!> \author  Darin Comeau, LANL
!> \date    Sep 26, 2017
!> \details This subroutine updates iceberg area state variable from the
!> iceberg mass state variable.
!
!-----------------------------------------------------------------------

  subroutine cice_berg_area_update(&
       nBergCategories, &
       bergMassCategory, &
       bergHeight, &
       bergDensity, &
       areaCell, &
       bergAreaCategory, &
       bergAreaCell)

    use cice_constants, only: cicePuny

    integer, intent(in) :: &
         nBergCategories

    real(kind=RKIND), dimension(:), intent(in) :: &
         bergMassCategory, & ! berg mass state variable (kg)
         bergHeight          ! berg category height (m)

    real(kind=RKIND), intent(in) :: &
         bergDensity, & ! berg density (kg/m^3)
         areaCell

    real(kind=RKIND), dimension(:), intent(out) :: &
         bergAreaCategory ! berg area state variable (dimensionless)

    real(kind=RKIND), intent(out) :: &
         bergAreaCell ! aggregate berg area (dimensionless)

    integer :: &
         iCategory

    do iCategory = 1, nBergCategories

       if (bergHeight(iCategory) > cicePuny) then

          ! calculate fractional berg area
          bergAreaCategory(iCategory) = bergMassCategory(iCategory) / (bergHeight(iCategory)*bergDensity)
          bergAreaCategory(iCategory) = bergAreaCategory(iCategory) / areaCell

       else

          bergAreaCategory(iCategory) = 0.0_RKIND

       endif

    enddo

    ! integrate quantity
    bergAreaCell = sum(bergAreaCategory(:))

    ! check that berg area does not exceed cell area
    if (bergAreaCell > 1.0_RKIND) then
       call mpas_log_write(&
            "berg area greater than one: $r ", realArgs=(/bergAreaCell/), &
            messageType = MPAS_LOG_CRIT)

    endif

  end subroutine cice_berg_area_update

!-----------------------------------------------------------------------

end module cice_berg_state
