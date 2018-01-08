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
         bergTracersHead01, &
         bergTracersHead02, &
         bergTracersHead03, &
         bergTracersHead04, &
         bergTracersHead05, &
         bergTracersHead06, &
         bergTracersHead07, &
         bergTracersHead08, &
         bergTracersHead09, &
         bergTracersHead10

    type(domain_type), intent(inout) :: &
         domain

    ! dcnote - possible this way? open to improvement
    ! integer :: iCategory

    ! call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    ! do iCategory = 1, nBergCategories

    !    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead//trim(iCategory), iCategory)

    ! enddo

    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead01, 1)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead02, 2)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead03, 3)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead04, 4)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead05, 5)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead06, 6)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead07, 7)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead08, 8)
    call cice_add_berg_tracers_to_linked_list(domain, bergTracersHead09, 9)
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
         bergTracersHead01, &
         bergTracersHead02, &
         bergTracersHead03, &
         bergTracersHead04, &
         bergTracersHead05, &
         bergTracersHead06, &
         bergTracersHead07, &
         bergTracersHead08, &
         bergTracersHead09, &
         bergTracersHead10

    type(domain_type), intent(inout) :: &
         domain

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(inout) :: &
         ierr  !< Input/Output:


    call cice_run_advection_incremental_remap(domain, bergTracersHead01, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead02, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead03, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead04, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead05, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead06, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead07, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead08, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead09, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)
    call cice_run_advection_incremental_remap(domain, bergTracersHead10, clock, ierr, &
                                              updateHaloInitIn = .true., updateHaloFinalIn = .true.)

  end subroutine cice_run_berg_advection

!-----------------------------------------------------------------------

end module cice_berg_advection
