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
         bergTracersHeadArray

    type(domain_type), intent(inout) :: &
         domain

    integer, pointer :: &
         nBergCategories

    integer :: iCategory

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nBergCategories", nBergCategories)

    allocate(bergTracersHeadArray(nBergCategories))

    do iCategory = 1, nBergCategories

       call cice_add_berg_tracers_to_linked_list(domain, &
                                                 bergTracersHeadArray(iCategory) % tracerTypePtr, &
                                                 iCategory)

    enddo

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
         bergTracersHeadArray

    type(domain_type), intent(inout) :: &
         domain

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(inout) :: &
         ierr  !< Input/Output:

    integer, pointer :: &
         nBergCategories

    integer :: iCategory

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nBergCategories", nBergCategories)

    do iCategory = 1, nBergCategories

       call cice_run_advection_incremental_remap(domain, &
                                                 bergTracersHeadArray(iCategory) % tracerTypePtr, &
                                                 clock, &
                                                 ierr, &
                                                 updateHaloInitIn = .true., &
                                                 updateHaloFinalIn = .true.)

    enddo

  end subroutine cice_run_berg_advection

!-----------------------------------------------------------------------

end module cice_berg_advection
