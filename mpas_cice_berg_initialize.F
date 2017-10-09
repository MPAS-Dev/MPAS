!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_bergs_initialize
!
!> \brief   Initialize iceberg model
!> \author  Darin Comeau, LANL
!> \date    19 May 2017
!> \details This module is intended to initialize the iceberg model.
!
!-----------------------------------------------------------------------

module cice_bergs_initialize

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
      cice_bergs_init

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_bergs_init
!
!> \brief   Initialize iceberg model         
!> \author  Darin Comeau, LANL
!> \date    19 May 2017
!> \details This routine is intended to initialize the iceberg model if
!> icebergs are turned on based on the namelist option config_use_bergs.
!
!-----------------------------------------------------------------------

  subroutine cice_bergs_init(domain)

   use cice_berg_velocity_solver, only: cice_init_berg_velocity_solver

    type(domain_type), intent(inout) :: &
         domain

    logical, pointer :: &
         config_use_bergs

    call MPAS_pool_get_config(domain % configs, "config_use_bergs", config_use_bergs)

    if (config_use_bergs) then

       call mpas_log_write(" Initialize icebergs...")

       call cice_init_berg_velocity_solver(domain)

       call init_berg_state(domain)

    endif

  end subroutine cice_bergs_init

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_state
!
!> \brief   Initialize iceberg model state         
!> \author  Darin Comeau, LANL
!> \date    19 May 2017
!> \details This routine is intended to initialize the iceberg model
!> physical state variables based on the namelist option
!> config_berg_initial_condition.
!
!-----------------------------------------------------------------------

  subroutine init_berg_state(domain)

    use cice_berg_state, only: berg_calving_event

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    character(len=strKIND), pointer :: &
         config_berg_initial_condition

    call MPAS_pool_get_config(domain % configs, "config_berg_initial_condition", config_berg_initial_condition)

    ! create uniform calving mask if needed
    call mpas_timer_start("Berg uniform calving mask")
    if (trim(config_berg_initial_condition) == "coastal_uniform") call init_berg_calving_mask(domain)
    if (trim(config_berg_initial_condition) == "regional_uniform") call init_berg_calving_mask(domain)
    call mpas_timer_stop("Berg uniform calving mask")

    block => domain % blocklist
    do while(associated(block))

       ! set berg calving masses
       if (trim(config_berg_initial_condition) == "coastal_uniform") then

          call init_berg_state_calving_uniform(block)

       elseif (trim(config_berg_initial_condition) == "regional_uniform") then

          call init_berg_state_calving_uniform(block)

       else

          call mpas_log_write(&
               "init_berg_state: config_berg_initial_condition unknown:"//trim(config_berg_initial_condition), &
               MPAS_LOG_CRIT)

       endif

       ! initial calving event
       call berg_calving_event(block)

       block => block % next
    enddo

    ! halo exchanges (needed?)
    call MPAS_dmpar_field_halo_exch(domain, 'bergMassCategory')
    call MPAS_dmpar_field_halo_exch(domain, 'bergAreaCategory')
    call MPAS_dmpar_field_halo_exch(domain, 'bergMassCell')
    call MPAS_dmpar_field_halo_exch(domain, 'bergAreaCell')

   end subroutine init_berg_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_calving_mask
!
!> \brief   Initialize iceberg calving mask
!> \author  Darin Comeau, LANL
!> \date    25 Sept 2017
!> \details This routine is intended to initialize the iceberg calving
!> mask, either coastal or regional (includes interior cells) for
!> spatially uniform calving.
!
!-----------------------------------------------------------------------

  subroutine init_berg_calving_mask(domain)

    use cice_constants, only: &
         ciceRadiansToDegrees

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         boundary, &
         berg_calving, &
         ocean_coupling

    integer, pointer :: &
         nCellsSolve

    integer, pointer :: &
         nBergCalvingCells, &
         nBergCalvingCellsTotal

    real(kind=RKIND), pointer :: &
         bergMaskLatMin, &
         bergMaskLatMax, &
         bergMaskLonMin, &
         bergMaskLonMax

    real(kind=RKIND), dimension(:), pointer :: &
         latCell, &
         lonCell

    integer, dimension(:), pointer :: &
         interiorCell, &
         bergCalvingMask, &
         landIceMask

    type (dm_info), pointer :: dminfo

    character(len=strKIND), pointer :: &
         config_berg_initial_condition

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_config(domain % configs, "config_berg_initial_condition", config_berg_initial_condition)
    call MPAS_pool_get_config(domain % configs, "config_berg_calving_mask_latitude_minimum", bergMaskLatMin)
    call MPAS_pool_get_config(domain % configs, "config_berg_calving_mask_latitude_maximum", bergMaskLatMax)
    call MPAS_pool_get_config(domain % configs, "config_berg_calving_mask_longitude_minimum", bergMaskLonMin)
    call MPAS_pool_get_config(domain % configs, "config_berg_calving_mask_longitude_maximum", bergMaskLonMax)

    if (bergMaskLonMin < 0.0_RKIND) bergMaskLonMin = bergMaskLonMin + 360.0_RKIND
    if (bergMaskLonMax < 0.0_RKIND) bergMaskLonMax = bergMaskLonMax + 360.0_RKIND

    block => domain % blocklist
    do while(associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
       call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_array(mesh, "latCell", latCell)
       call MPAS_pool_get_array(mesh, "lonCell", lonCell)
       call MPAS_pool_get_array(boundary, "interiorCell", interiorCell)
       call MPAS_pool_get_array(berg_calving, "bergCalvingMask", bergCalvingMask)
       call MPAS_pool_get_array(berg_calving, "nBergCalvingCells", nBergCalvingCells)
       call MPAS_pool_get_array(berg_calving, "nBergCalvingCellsTotal", nBergCalvingCellsTotal)
       call MPAS_pool_get_array(ocean_coupling, "landIceMask", landIceMask)

       if (trim(config_berg_initial_condition) == "coastal_uniform") then

          do iCell = 1, nCellsSolve

             ! coastal cells only
             if (interiorCell(iCell) == 0 .and. &
                landIceMask(iCell) == 0 .and. &
                latCell(iCell)*ciceRadiansToDegrees >= bergMaskLatMin .and. &
                latCell(iCell)*ciceRadiansToDegrees <= bergMaskLatMax .and. &
                lonCell(iCell)*ciceRadiansToDegrees >= bergMaskLonMin .and. &
                lonCell(iCell)*ciceRadiansToDegrees <= bergMaskLonMax) then

                bergCalvingMask(iCell) = 1
                nBergCalvingCells = nBergCalvingCells + 1

             endif

          enddo

       elseif (trim(config_berg_initial_condition) == "regional_uniform") then

          do iCell = 1, nCellsSolve

             ! include interior cells
             if (landIceMask(iCell) == 0 .and. &
                latCell(iCell)*ciceRadiansToDegrees >= bergMaskLatMin .and. &
                latCell(iCell)*ciceRadiansToDegrees <= bergMaskLatMax .and. &
                lonCell(iCell)*ciceRadiansToDegrees >= bergMaskLonMin .and. &
                lonCell(iCell)*ciceRadiansToDegrees <= bergMaskLonMax) then

                bergCalvingMask(iCell) = 1
                nBergCalvingCells = nBergCalvingCells + 1

             endif

          enddo

       else

          call mpas_log_write(&
               "config_berg_initial_condition unknown:"//trim(config_berg_initial_condition), &
                MPAS_LOG_CRIT)

       endif

       block => block % next
    enddo

    ! sum up calving cells over all blocks
    call mpas_dmpar_sum_int(domain % dminfo, nBergCalvingCells, nBergCalvingCellsTotal)
    call mpas_log_write(' # Berg calving cells $i', intArgs=(/nBergCalvingCellsTotal/))


  end subroutine init_berg_calving_mask

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_state_calving_uniform
!
!> \brief   Initialize with icebergs calving mass
!> \author  Darin Comeau, LANL
!> \date    25 Sept 2017
!> \details This routine is intended to initialize the iceberg calving
!> mass variables uniformly distributed on the calving mask.
!
!-----------------------------------------------------------------------

  subroutine init_berg_state_calving_uniform(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         berg_calving

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer, dimension(:), pointer :: &
         bergCalvingMask

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergCalvingMass

    real(kind=RKIND), dimension(:), pointer :: &
         bergCalvingDistribution

    integer, pointer :: &
         nBergCalvingCellsTotal

    real(kind=RKIND), pointer :: &
         annualCalvingRate, &
         annualCalvingFrequency

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_config(block % configs, "config_berg_annual_calving_rate", annualCalvingRate)
    call MPAS_pool_get_config(block % configs, "config_berg_annual_calving_frequency", annualCalvingFrequency)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(berg_calving, "bergCalvingMask", bergCalvingMask)
    call MPAS_pool_get_array(berg_calving, "nBergCalvingCellsTotal", nBergCalvingCellsTotal)
    call MPAS_pool_get_array(berg_calving, "bergCalvingDistribution", bergCalvingDistribution)
    call MPAS_pool_get_array(berg_calving, "bergCalvingMass", bergCalvingMass)

    do iCell = 1, nCellsSolve

       if (bergCalvingMask(iCell) == 1) then

          do iCategory = 1, nBergCategories

             bergCalvingMass(iCategory,iCell) = bergCalvingDistribution(iCategory)*annualCalvingRate &
                                              / (annualCalvingFrequency*nBergCalvingCellsTotal)

          enddo

       endif

    enddo

  end subroutine init_berg_state_calving_uniform

!-----------------------------------------------------------------------

end module cice_bergs_initialize
