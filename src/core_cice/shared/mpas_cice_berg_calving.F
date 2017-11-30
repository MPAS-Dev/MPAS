!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_berg_calving
!
!> \brief   Set up calving fluxes for the iceberg model
!> \author  Darin Comeau, LANL
!> \date    19 May 2017
!> \details This module is intended to set up and handle calving fluxes
!> for the iceberg model.
!
!-----------------------------------------------------------------------

module cice_berg_calving

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
      cice_init_berg_calving_fluxes

contains

!-----------------------------------------------------------------------
! Set up calving fluxes
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_berg_calving_fluxes
!
!> \brief   Initialize iceberg calving fluxes
!> \author  Darin Comeau, LANL
!> \date    19 May 2017
!> \details This routine is intended to initialize the iceberg model
!> calving fluxes based on the namelist option
!> config_berg_calving_flux_type.
!
!-----------------------------------------------------------------------

  subroutine cice_init_berg_calving_fluxes(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    character(len=strKIND), pointer :: &
         config_berg_calving_flux_type, & ! indicates how calving fluxes will be set
         config_berg_calving_type         ! indicates how calved mass will be released to model

    call MPAS_pool_get_config(domain % configs, "config_berg_calving_flux_type", config_berg_calving_flux_type)
    call MPAS_pool_get_config(domain % configs, "config_berg_calving_type", config_berg_calving_type)

    ! create uniform calving mask if needed
    call mpas_timer_start("Berg uniform calving mask")
    if ((trim(config_berg_calving_flux_type) == "coastal_uniform") .or. &
        (trim(config_berg_calving_flux_type) == "regional_uniform")) then

       call init_berg_calving_mask(domain)

    endif
    call mpas_timer_stop("Berg uniform calving mask")

    block => domain % blocklist
    do while(associated(block))

       ! set berg calving fluxes
       if ((trim(config_berg_calving_flux_type) == "coastal_uniform") .or. &
           (trim(config_berg_calving_flux_type) == "regional_uniform")) then

          call init_berg_calving_flux_uniform(block)

       elseif (trim(config_berg_calving_flux_type) == "obs_Rignot2013") then

          call init_berg_calving_flux_obs(block)

       else

          call mpas_log_write(&
               "init_berg_state: config_berg_calving_flux_type unknown:"//trim(config_berg_calving_flux_type), &
               MPAS_LOG_CRIT)

       endif

       ! calculate berg mass threshold for calving if needed
       if (trim(config_berg_calving_type) == "mass_threshold") &
          call init_berg_calving_mass_threshold(block)

       block => block % next
    enddo

   end subroutine cice_init_berg_calving_fluxes

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
         nBergCalvingCells, &   ! number of calving cells on this processor
         nBergCalvingCellsTotal ! number of calving cells on all processors

    real(kind=RKIND), pointer :: &
         bergMaskLatMin, & ! latitude range used to define calving cells
         bergMaskLatMax, &
         bergMaskLonMin, & ! longitude range used to define calving cells
         bergMaskLonMax

    real(kind=RKIND), dimension(:), pointer :: &
         latCell, &
         lonCell

    integer, dimension(:), pointer :: &
         interiorCell, &
         landIceMask

    integer, dimension(:,:), pointer :: &
         bergCalvingMask ! mask = 1 indicates a calving cell, = 0 elsewhere

    type (dm_info), pointer :: dminfo

    character(len=strKIND), pointer :: &
         config_berg_calving_flux_type ! indicates how calved mass will be released to model

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_config(domain % configs, "config_berg_calving_flux_type", config_berg_calving_flux_type)
    call MPAS_pool_get_config(domain % configs, "config_berg_calving_mask_latitude_minimum", bergMaskLatMin)
    call MPAS_pool_get_config(domain % configs, "config_berg_calving_mask_latitude_maximum", bergMaskLatMax)
    call MPAS_pool_get_config(domain % configs, "config_berg_calving_mask_longitude_minimum", bergMaskLonMin)
    call MPAS_pool_get_config(domain % configs, "config_berg_calving_mask_longitude_maximum", bergMaskLonMax)

    ! convert longitudes to 0 - 360 indexing
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

       if (trim(config_berg_calving_flux_type) == "coastal_uniform") then

          do iCell = 1, nCellsSolve

             ! identify coastal cells only in designated range
             if (interiorCell(iCell) == 0 .and. &
                landIceMask(iCell) == 0 .and. &
                latCell(iCell)*ciceRadiansToDegrees >= bergMaskLatMin .and. &
                latCell(iCell)*ciceRadiansToDegrees <= bergMaskLatMax .and. &
                lonCell(iCell)*ciceRadiansToDegrees >= bergMaskLonMin .and. &
                lonCell(iCell)*ciceRadiansToDegrees <= bergMaskLonMax) then

                bergCalvingMask(:,iCell) = 1
                nBergCalvingCells = nBergCalvingCells + 1

             endif

          enddo

       elseif (trim(config_berg_calving_flux_type) == "regional_uniform") then

          do iCell = 1, nCellsSolve

             ! identify coastal and interior cells in designated range
             if (landIceMask(iCell) == 0 .and. &
                latCell(iCell)*ciceRadiansToDegrees >= bergMaskLatMin .and. &
                latCell(iCell)*ciceRadiansToDegrees <= bergMaskLatMax .and. &
                lonCell(iCell)*ciceRadiansToDegrees >= bergMaskLonMin .and. &
                lonCell(iCell)*ciceRadiansToDegrees <= bergMaskLonMax) then

                bergCalvingMask(:,iCell) = 1
                nBergCalvingCells = nBergCalvingCells + 1

             endif

          enddo

       else

          call mpas_log_write(&
               "config_berg_calving_flux_type unknown:"//trim(config_berg_calving_flux_type), &
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
!  init_berg_calving_flux_uniform
!
!> \brief   Initialize with icebergs calving mass
!> \author  Darin Comeau, LANL
!> \date    25 Sept 2017
!> \details This routine is intended to initialize the iceberg calving
!> flux uniformly distributed on the calving mask.
!
!-----------------------------------------------------------------------

  subroutine init_berg_calving_flux_uniform(block)

    use cice_constants, only: ciceSecondsPerYear

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         berg_calving

    integer, pointer :: &
         nCellsSolve

    integer, dimension(:,:), pointer :: &
         bergCalvingMask ! mask = 1 indicates a calving cell, = 0 elsewhere

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell, &
         bergCalvingFlux ! calving flux forcing (kg/m^2/s)

    integer, pointer :: &
         nBergCalvingCellsTotal ! total number of cells that receive calving flux forcing

    real(kind=RKIND), pointer :: &
         annualCalvingRate ! total annual calving rate (kg/yr) to be divided among calving cells

    integer :: &
         iCell

    call MPAS_pool_get_config(block % configs, "config_berg_annual_calving_rate", annualCalvingRate)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)
    call MPAS_pool_get_array(berg_calving, "bergCalvingMask", bergCalvingMask)
    call MPAS_pool_get_array(berg_calving, "bergCalvingFlux", bergCalvingFlux)
    call MPAS_pool_get_array(berg_calving, "nBergCalvingCellsTotal", nBergCalvingCellsTotal)

    do iCell = 1, nCellsSolve

       if (maxval(bergCalvingMask(:,iCell)) > 0) then

          ! distribute calving flux evenly among all calving cells
          bergCalvingFlux(iCell) = annualCalvingRate / nBergCalvingCellsTotal ! kg/yr
          bergCalvingFlux(iCell) = bergCalvingFlux(iCell) / (areaCell(iCell)*ciceSecondsPerYear) ! kg/m^2/s

       endif

    enddo

  end subroutine init_berg_calving_flux_uniform

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_calving_flux_obs
!
!> \brief   Initialize with icebergs calving mass
!> \author  Darin Comeau, LANL
!> \date    25 Sept 2017
!> \details This routine is intended to initialize the iceberg calving
!> flux from observations (Rignot et al 2013).
!
!-----------------------------------------------------------------------

  subroutine init_berg_calving_flux_obs(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         berg_calving

    integer, pointer :: &
         nCellsSolve

    integer, dimension(:,:), pointer :: &
         bergCalvingMask ! mask = 1 indicates a calving cell, = 0 elsewhere

    real(kind=RKIND), dimension(:), pointer :: &
         bergCalvingFlux, & ! ! calving flux forcing (kg/m^2/s)
         calvingFluxObs ! ! calving flux forcing read in from file (kg/m^2/s)

    integer :: &
         iCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

    call MPAS_pool_get_array(berg_calving, "bergCalvingMask", bergCalvingMask)
    call MPAS_pool_get_array(berg_calving, "bergCalvingFlux", bergCalvingFlux)
    call MPAS_pool_get_array(berg_calving, "calvingFluxObs", calvingFluxObs)

    do iCell = 1, nCellsSolve

       ! if input calving flux is positive, indicate as a calving cell
       if (calvingFluxObs(iCell) > 0.0_RKIND) then

         bergCalvingMask(:,iCell) = 1

         bergCalvingFlux(iCell) = calvingFluxObs(iCell)

       endif

    enddo

  end subroutine init_berg_calving_flux_obs

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_calving_mass_threshold
!
!> \brief   Calculate mass of one berg in each category
!> \author  Darin Comeau, LANL
!> \date    3 Oct 2017
!> \details This routine calculates the mass of one berg in each category
!> based on the calving distribution. This mass is intended to be used
!> as a threshold for releasing calving mass that has been accumulated by
!> the calving flux.
!
!-----------------------------------------------------------------------

  subroutine init_berg_calving_mass_threshold(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         berg_calving

    integer, pointer :: &
         nBergCategories ! number of berg categories

    real(kind=RKIND), dimension(:), pointer :: &
         bergCalvingLength, & ! berg category's length from calving distribution
         bergCalvingHeight, & ! berg category's height from calving distribution
         oneCalvedBergMass ! mass of one categoruy's berg from calving distribution

    real(kind=RKIND), pointer :: &
         bergDensity ! berg density (kg/m^3)

    integer :: &
         iCell, &
         iCategory

    ! this parameter assume a horizontal aspect ratio of 1:1.5
    real(kind=RKIND), parameter :: &
         bergWidthFactor = 2.0_RKIND/3.0_RKIND

    call MPAS_pool_get_config(block % configs, "config_berg_density", bergDensity)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)

    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(berg_calving, "bergCalvingLength", bergCalvingLength)
    call MPAS_pool_get_array(berg_calving, "bergCalvingHeight", bergCalvingHeight)
    call MPAS_pool_get_array(berg_calving, "oneCalvedBergMass", oneCalvedBergMass)

    do iCategory = 1, nBergCategories

       ! calculate mass of a single berg in each category from the calving distribution
       oneCalvedBergMass(iCategory) = bergWidthFactor*bergCalvingLength(iCategory)**2 &
                                    * bergCalvingHeight(iCategory)*bergDensity

    enddo

  end subroutine init_berg_calving_mass_threshold

!-----------------------------------------------------------------------

end module cice_berg_calving
