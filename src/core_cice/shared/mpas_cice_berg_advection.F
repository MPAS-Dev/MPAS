!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_berg_advection
!
!> \brief   Iceberg model advection
!> \author  Darin Comeau, LANL
!> \date    29 Nov 2017
!> \details This module is intended to advect the iceberg model.
!
!-----------------------------------------------------------------------

module cice_berg_advection

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
      cice_init_berg_tracers, &
      cice_run_berg_advection

contains

!-----------------------------------------------------------------------
! Set up iceberg tracers
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_berg_tracers
!
!> \brief   Initialize iceberg tracers
!> \author  Darin Comeau, LANL
!> \date    29 Nov 2017
!> \details This routine is intended to construct iceberg tracers.
!
!-----------------------------------------------------------------------

  subroutine cice_init_berg_tracers(domain)

    use cice_advection_incremental_remap_tracers, only: &
         cice_add_berg_tracers_to_linked_list, &
         bergTracersHead1, &
         bergTracersHead2, &
         bergTracersHead3, &
         bergTracersHead4, &
         bergTracersHead5, &
         bergTracersHead6, &
         bergTracersHead7, &
         bergTracersHead8, &
         bergTracersHead9, &
         bergTracersHead10

    type(domain_type), intent(inout) :: &
         domain

    ! ! dcnote - possible this way?
    ! integer :: iCategory

    ! call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    ! do iCategory = 1, nBergCategories

    !    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead//trim(iCategory), iCategory)

    ! enddo

    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead1, 1)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead2, 2)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead3, 3)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead4, 4)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead5, 5)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead6, 6)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead7, 7)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead8, 8)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead9, 9)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead10, 10)

  end subroutine cice_init_berg_tracers

!-----------------------------------------------------------------------
! Iceberg advection
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_run_berg_advection
!
!> \brief   Advect icebergs
!> \author  Darin Comeau, LANL
!> \date    29 Nov 2017
!> \details This routine is intended to advect icebergs by incremental
!> remapping (currently only option).
!
!-----------------------------------------------------------------------

  subroutine cice_run_berg_advection(domain, clock, ierr)

    use cice_advection_incremental_remap, only: &
         cice_run_advection_incremental_remap

    use cice_advection_incremental_remap_tracers, only: &
         bergTracersHead1, &
         bergTracersHead2, &
         bergTracersHead3, &
         bergTracersHead4, &
         bergTracersHead5, &
         bergTracersHead6, &
         bergTracersHead7, &
         bergTracersHead8, &
         bergTracersHead9, &
         bergTracersHead10

    type(domain_type), intent(inout) :: &
         domain

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(inout) :: &
         ierr  !< Input/Output:

    call cice_run_advection_incremental_remap(domain, bergTracersHead1, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead2, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead3, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead4, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead5, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead6, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead7, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead8, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead9, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead10, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)

  end subroutine cice_run_berg_advection

!-----------------------------------------------------------------------

end module cice_berg_advection
