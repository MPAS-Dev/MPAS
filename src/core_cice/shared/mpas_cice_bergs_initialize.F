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

   use mpas_dmpar
   use mpas_derived_types
   use mpas_pool_routines
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

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    character(len=strKIND), pointer :: &
         config_berg_initial_condition

    call MPAS_pool_get_config(domain % configs, "config_berg_initial_condition", config_berg_initial_condition)

    if (trim(config_berg_initial_condition) == "bergs_coastal_uniform") call init_berg_coastal_mask(domain)

    block => domain % blocklist
    do while(associated(block))

       ! set berg masses
       if (trim(config_berg_initial_condition) == "no_bergs") then

          call init_berg_state_no_bergs(block)

       elseif (trim(config_berg_initial_condition) == "bergs_everywhere") then

          call init_berg_state_bergs_everywhere(block)

       elseif (trim(config_berg_initial_condition) == "bergs_coastal_uniform") then

          call init_berg_state_bergs_coastal_uniform(block)

       else

          call mpas_log_write(&
               "init_berg_state: config_berg_initial_condition unknown:"//trim(config_berg_initial_condition), &
               MPAS_LOG_CRIT)

       endif

       ! set berg size dimensions
       call init_berg_size(block)

       block => block % next
    enddo

    ! halo exchanges
    call MPAS_dmpar_field_halo_exch(domain, 'bergMassCategory')
    call MPAS_dmpar_field_halo_exch(domain, 'bergMassCell')

   end subroutine init_berg_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_coastal_mask
!
!> \brief   Initialize iceberg coastal mask
!> \author  Darin Comeau, LANL
!> \date    25 Sept 2017
!> \details This routine is intended to initialize the iceberg coastal
!> mask for spatially uniform calving.
!
!-----------------------------------------------------------------------

  subroutine init_berg_coastal_mask(domain)

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

!     integer, dimension(:), pointer :: &
!          nBergCalvingCells, &
!          nBergCalvingCellsTotal

    real(kind=RKIND), pointer :: &
         bergMaskLatMax

    real(kind=RKIND), dimension(:), pointer :: &
         latCell

    integer, dimension(:), pointer :: &
         interiorCell, &
         bergCalvingMask, &
         landIceMask

    type (dm_info), pointer :: dminfo

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_config(domain % configs, "config_berg_antarctic_mask_latitude_maximum", bergMaskLatMax)

    block => domain % blocklist
    do while(associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
       call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_array(mesh, "latCell", latCell)
       call MPAS_pool_get_array(boundary, "interiorCell", interiorCell)
       call MPAS_pool_get_array(berg_calving, "bergCalvingMask", bergCalvingMask)
       call MPAS_pool_get_array(berg_calving, "nBergCalvingCells", nBergCalvingCells)
       call MPAS_pool_get_array(berg_calving, "nBergCalvingCellsTotal", nBergCalvingCellsTotal)
       call MPAS_pool_get_array(ocean_coupling, "landIceMask", landIceMask)

! only allow if neighboring a land ice cell
!        do iCell = 1, nCellsSolve

!           if (interiorCell(iCell) == 0 .and. &
!              (landIceMask(iCell) == 0 .and. &
!              latCell(iCell) < bergMaskLatMax) then

!              ! see if neighboring a land ice cell
!              do iCellonCell = 1, nEdgesonCell(iCell)

!                 iCellNeighbor = cellsOnCells(iCellonCell,iCell)

!                 if (landIceMask(iCellNeighbor) == 1) bergCalvingMask(iCell) = 1

!              enddo

!           endif

!           ! counter for calving cells
!           if (bergCalvingMask(iCell) == 1) nBergCalvingCells = nBergCalvingCells + 1

!        enddo

! count all coastal cells
       do iCell = 1, nCellsSolve

          if (interiorCell(iCell) == 0 .and. &
             landIceMask(iCell) == 0 .and. &
             latCell(iCell) < bergMaskLatMax) then

             bergCalvingMask(iCell) = 1
             nBergCalvingCells = nBergCalvingCells + 1

          endif

       enddo

       nBergCalvingCellsTotal = nBergCalvingCells ! for now until do global sum

       block => block % next
    enddo

    ! sum up calving cells over all blocks
    call mpas_dmpar_sum_int(domain % dminfo, nBergCalvingCells, nBergCalvingCellsTotal)
    call mpas_log_write('nBergCalvingCellsTotal $i', intArgs=(/nBergCalvingCellsTotal/))


  end subroutine init_berg_coastal_mask

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_state_bergs_coastal_uniform
!
!> \brief   Initialize with icebergs along coastal cells
!> \author  Darin Comeau, LANL
!> \date    25 Sept 2017
!> \details This routine is intended to initialize the iceberg mass
!> variables uniformly distributed around the coast.
!
!-----------------------------------------------------------------------

  subroutine init_berg_state_bergs_coastal_uniform(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate, &
         berg_calving

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer, dimension(:), pointer :: &
         bergCalvingMask

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory

    real(kind=RKIND), dimension(:), pointer :: &
         bergMassCell

    real(kind=RKIND), dimension(:), pointer :: &
         bergCalvingDistribution

    integer, pointer :: &
         nBergCalvingCellsTotal

    real(kind=RKIND), pointer :: &
         annualAntarcticCalvingRate, &
         annualAntarcticCalvingFrequency

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
    call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(tracers, "bergMassCategory", bergMassCategory,1)
    call MPAS_pool_get_array(tracers_aggregate, "bergMassCell", bergMassCell)
    call MPAS_pool_get_array(berg_calving, "bergCalvingMask", bergCalvingMask)
    call MPAS_pool_get_array(berg_calving, "nBergCalvingCellsTotal", nBergCalvingCellsTotal)
    call MPAS_pool_get_array(berg_calving, "bergCalvingDistribution", bergCalvingDistribution)

    call MPAS_pool_get_config(block % configs, "config_berg_annual_antarctic_calving_rate", annualAntarcticCalvingRate)
    call MPAS_pool_get_config(block % configs, "config_berg_annual_antarctic_calving_frequency", annualAntarcticCalvingFrequency)

    do iCell = 1, nCellsSolve

       if (bergCalvingMask(iCell) == 1) then

          do iCategory = 1, nBergCategories

             bergMassCategory(iCategory,iCell) = bergCalvingDistribution(iCategory)*annualAntarcticCalvingRate &
                                               / (annualAntarcticCalvingFrequency*nBergCalvingCellsTotal)

          enddo

       endif

       ! integrated quantities
       bergMassCell(iCell) = sum(bergMassCategory(:,iCell))

    enddo


  end subroutine init_berg_state_bergs_coastal_uniform

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_state_no_bergs
!
!> \brief   Initialize with no icebergs everywhere
!> \author  Darin Comeau, LANL
!> \date    19 May 2017
!> \details This routine is intended to initialize the iceberg mass
!> variables as 0 at all grid cells.
!
!-----------------------------------------------------------------------

  subroutine init_berg_state_no_bergs(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory

    real(kind=RKIND), dimension(:), pointer :: &
         bergMassCell

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(tracers, "bergMassCategory", bergMassCategory,1)
    call MPAS_pool_get_array(tracers_aggregate, "bergMassCell", bergMassCell)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          bergMassCategory(iCategory,iCell) = 0.0_RKIND

       enddo

       ! integrated quantities
       bergMassCell(iCell) = sum(bergMassCategory(:,iCell))

    enddo


  end subroutine init_berg_state_no_bergs

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_state_bergs_everywhere
!
!> \brief   Initialize with icebergs everywhere (for testing)
!> \author  Darin Comeau, LANL
!> \date    19 May 2017
!> \details This routine is intended to initialize the iceberg mass
!> variables as 100 kg at all grid cells.
!
!-----------------------------------------------------------------------

  subroutine init_berg_state_bergs_everywhere(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory

    real(kind=RKIND), dimension(:), pointer :: &
         bergMassCell

    real(kind=RKIND), pointer :: &
         config_berg_initial_mass

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(tracers, "bergMassCategory", bergMassCategory,1)
    call MPAS_pool_get_array(tracers_aggregate, "bergMassCell", bergMassCell)

    call MPAS_pool_get_config(block % configs, "config_berg_initial_mass", config_berg_initial_mass)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          bergMassCategory(iCategory,iCell) = config_berg_initial_mass

       enddo

       ! integrated quantities
       bergMassCell(iCell) = sum(bergMassCategory(:,iCell))

    enddo


  end subroutine init_berg_state_bergs_everywhere

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_size
!
!> \brief   Initialize icebergs size dimensions
!> \author  Darin Comeau, LANL
!> \date    6 June 2017
!> \details This routine is intended to initialize the iceberg size
!> dimensions as newly calved bergs.
!
!-----------------------------------------------------------------------

  subroutine init_berg_size(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate, &
         berg_calving

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory, &
         bergAreaCategory, &
         bergLength, &
         bergHeight

    real(kind=RKIND), dimension(:), pointer :: &
         bergCalvingLength, &
         bergCalvingHeight

    real(kind=RKIND), dimension(:), pointer :: &
         bergAreaCell

    real(kind=RKIND), pointer :: &
         bergDensity

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_config(block % configs, "config_berg_density", bergDensity)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
    call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(tracers, "bergMassCategory", bergMassCategory,1)
    call MPAS_pool_get_array(tracers, "bergAreaCategory", bergAreaCategory,1)
    call MPAS_pool_get_array(tracers, "bergLength", bergLength,1)
    call MPAS_pool_get_array(tracers, "bergHeight", bergHeight,1)
    call MPAS_pool_get_array(tracers_aggregate, "bergAreaCell", bergAreaCell)
    call MPAS_pool_get_array(berg_calving, "bergCalvingLength", bergCalvingLength)
    call MPAS_pool_get_array(berg_calving, "bergCalvingHeight", bergCalvingHeight)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          if (bergMassCategory(iCategory,iCell) > 0.0_RKIND) then

             bergLength(iCategory,iCell) = bergCalvingLength(iCategory)
             bergHeight(iCategory,iCell) = bergCalvingHeight(iCategory)
             bergAreaCategory(iCategory,iCell) = bergMassCategory(iCategory,iCell)/(bergHeight(iCategory,iCell)*bergDensity)

          else

             bergLength(iCategory,iCell) = 0.0_RKIND
             bergHeight(iCategory,iCell) = 0.0_RKIND
             bergAreaCategory(iCategory,iCell) = 0.0_RKIND

          endif

       enddo

       ! integrated quantities
       bergAreaCell(iCell) = sum(bergAreaCategory(:,iCell))

    enddo

  end subroutine init_berg_size

!-----------------------------------------------------------------------

end module cice_bergs_initialize
