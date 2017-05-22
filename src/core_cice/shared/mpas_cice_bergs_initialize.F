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
!> \details This routine is intended to initialize the iceberg model.
!
!-----------------------------------------------------------------------

   subroutine cice_bergs_init(domain)

      type(domain_type), intent(inout) :: &
         domain

      call mpas_log_write("cice_bergs_init: Initialize iceberg model")

      ! initialize berg state
      call mpas_log_write(" Initialize iceberg state...")
      call init_berg_state(domain)

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
!> berg_initial_condition_type.
!
!-----------------------------------------------------------------------

   subroutine init_berg_state(domain)

      type(domain_type), intent(inout) :: &
         domain

      type(block_type), pointer :: &
         block

      logical, pointer :: &
         config_use_bergs, &
         config_do_restart

      character(len=strKIND), pointer :: &
         config_berg_initial_condition_type

      call MPAS_pool_get_config(domain % configs, "config_use_bergs", config_use_bergs)
      call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

      if (.not. config_do_restart .and. config_use_bergs) then

         call MPAS_pool_get_config(domain % configs, "config_berg_initial_condition_type", config_berg_initial_condition_type)

         ! set berg masses
         block => domain % blocklist
         do while(associated(block))

            if (trim(config_berg_initial_condition_type) == "test_everywhere") then

               call init_berg_state_everywhere(&
                  block, &
                  domain % configs)

            else

               call mpas_log_write(&
                  "init_berg_state: config_berg_initial_condition_type unknown:"//trim(config_berg_initial_condition_type), &
                  MPAS_LOG_CRIT)

            endif

            block => block % next
         enddo

         ! set berg velocities (later)

         ! halo exchanges
         call MPAS_dmpar_field_halo_exch(domain, 'bergMassCategory', timeLevel=1)
         call MPAS_dmpar_field_halo_exch(domain, 'bergMassCell', timeLevel=1)

      endif ! restart / use_bergs

   end subroutine init_berg_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_state_everywhere
!
!> \brief   Initialize icebergs everywhere            
!> \author  Darin Comeau, LANL
!> \date    19 May 2017
!> \details This routine is intended to initialize the iceberg mass 
!> variables as 0 at all grid cells.
!
!-----------------------------------------------------------------------

   subroutine init_berg_state_everywhere(&
      block, &
      configs)

      type(block_type), intent(inout) :: &
         block

      type(MPAS_pool_type), pointer, intent(in) :: &
         configs

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

      call MPAS_pool_get_array(tracers, "bergMassCategory", bergMassCategory, 1)
      call MPAS_pool_get_array(tracers_aggregate, "bergMassCell", bergMassCell, 1)


      do iCell = 1, nCellsSolve

         do iCategory = 1, nBergCategories

            bergMassCategory(iCategory,iCell) = 0.0_RKIND

         enddo

      enddo

      ! integrated quantities
      do iCell = 1, nCellsSolve

         bergMassCell(iCell) = sum(bergMassCategory(:,iCell))

      enddo

   end subroutine init_berg_state_everywhere

!-----------------------------------------------------------------------

end module cice_bergs_initialize
