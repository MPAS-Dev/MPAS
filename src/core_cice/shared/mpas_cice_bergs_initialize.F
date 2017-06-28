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

    type(domain_type), intent(inout) :: &
         domain

    logical, pointer :: &
         config_use_bergs

    call MPAS_pool_get_config(domain % configs, "config_use_bergs", config_use_bergs)

    if (config_use_bergs) then

       call mpas_log_write("Initialize iceberg model")

       ! initialize berg state
       call mpas_log_write(" Initialize iceberg state...")
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

    block => domain % blocklist
    do while(associated(block))

       ! set berg masses
       if (trim(config_berg_initial_condition) == "no_bergs") then

          call init_berg_state_no_bergs(block)

       else

          call mpas_log_write(&
               "init_berg_state: config_berg_initial_condition unknown:"//trim(config_berg_initial_condition), &
               MPAS_LOG_CRIT)

       endif

       ! set berg size dimensions
       call init_berg_size(block)

       ! set berg velocities
       call init_berg_velocity(block)

       ! set berg fluxes
       call init_berg_fluxes(block)

       block => block % next
    enddo

    ! halo exchanges
    call MPAS_dmpar_field_halo_exch(domain, 'bergMassCategory')
    call MPAS_dmpar_field_halo_exch(domain, 'bergMassCell')

   end subroutine init_berg_state

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

    call MPAS_pool_get_array(tracers, "bergMassCategory", bergMassCategory)
    call MPAS_pool_get_array(tracers_aggregate, "bergMassCell", bergMassCell)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          bergMassCategory(iCategory,iCell) = 0.0_RKIND

       enddo

    enddo

    ! integrated quantities
    do iCell = 1, nCellsSolve

       bergMassCell(iCell) = sum(bergMassCategory(:,iCell))

    enddo

  end subroutine init_berg_state_no_bergs

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
         bergs

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory, &
         bergLength, &
         bergHeight

    real(kind=RKIND), dimension(:), pointer :: &
         bergCalvingLength, &
         bergCalvingHeight

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "bergs", bergs)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(tracers, "bergMassCategory", bergMassCategory)
    call MPAS_pool_get_array(tracers, "bergLength", bergLength)
    call MPAS_pool_get_array(tracers, "bergHeight", bergHeight)
    call MPAS_pool_get_array(bergs, "bergCalvingLength", bergCalvingLength)
    call MPAS_pool_get_array(bergs, "bergCalvingHeight", bergCalvingHeight)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          if (bergMassCategory(iCategory,iCell) > 0.0_RKIND) then

             bergLength(iCategory,iCell) = bergCalvingLength(iCategory)
             bergHeight(iCategory,iCell) = bergCalvingHeight(iCategory)

          else

             bergLength(iCategory,iCell) = 0.0_RKIND
             bergHeight(iCategory,iCell) = 0.0_RKIND

          endif

       enddo

    enddo

  end subroutine init_berg_size

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_velocity
!
!> \brief   Initialize icebergs velocity
!> \author  Darin Comeau, LANL
!> \date    6 June 2017
!> \details This routine is intended to initialize iceberg velocities.
!
!-----------------------------------------------------------------------

  subroutine init_berg_velocity(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         bergs

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    real(kind=RKIND), dimension(:,:), pointer :: &
         uBergVelocity, &
         vBergVelocity

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "bergs", bergs)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(tracers, "uBergVelocity", uBergVelocity)
    call MPAS_pool_get_array(tracers, "vBergVelocity", vBergVelocity)

    do iCell = 1, nCellsSolve

       do iCategory = 1, nBergCategories

          uBergVelocity(iCategory,iCell) = 0.0_RKIND
          vBergVelocity(iCategory,iCell) = 0.0_RKIND

       enddo

    enddo

    ! initialize velocity solver type (later)

  end subroutine init_berg_velocity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_fluxes
!
!> \brief   Initialize icebergs fluxes
!> \author  Darin Comeau, LANL
!> \date    6 June 2017
!> \details This routine is intended to initialize freshwater and latent
!> heat fluxes due to icebergs.
!
!-----------------------------------------------------------------------

  subroutine init_berg_fluxes(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate, &
         bergs

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergFreshwaterFluxCategory, &
         bergLatentHeatFluxCategory

    real(kind=RKIND), dimension(:), pointer :: &
         bergFreshwaterFluxCell, &
         bergLatentHeatFluxCell

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
    call MPAS_pool_get_subpool(block % structs, "bergs", bergs)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(tracers, "bergFreshwaterFluxCategory", bergFreshwaterFluxCategory)
    call MPAS_pool_get_array(tracers, "bergLatentHeatFluxCategory", bergLatentHeatFluxCategory)
    call MPAS_pool_get_array(tracers_aggregate, "bergFreshwaterFluxCell", bergFreshwaterFluxCell)
    call MPAS_pool_get_array(tracers_aggregate, "bergLatentHeatFluxCell", bergLatentHeatFluxCell)

    do iCell = 1, nCellsSolve

       bergFreshwaterFluxCell(iCell) = 0.0_RKIND
       bergLatentHeatFluxCell(iCell) = 0.0_RKIND

       do iCategory = 1, nBergCategories

          bergLatentHeatFluxCategory(iCategory,iCell) = 0.0_RKIND
          bergLatentHeatFluxCategory(iCategory,iCell) = 0.0_RKIND

       enddo

    enddo

  end subroutine init_berg_fluxes

!-----------------------------------------------------------------------

end module cice_bergs_initialize
