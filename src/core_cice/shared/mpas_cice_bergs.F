!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_bergs
!
!> \brief   Iceberg model
!> \author  Darin Comeau, LANL
!> \date    29 Sept 2017
!> \details This module is intended to interface the iceberg model with 
!> MPAS-seaice.
!
!-----------------------------------------------------------------------

module cice_bergs

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
      cice_bergs_init, &
      cice_run_berg_dynamics, &
      cice_run_berg_postdynamics

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

    use cice_berg_initialize, only: init_berg_state
    use cice_berg_velocity_solver, only: init_berg_velocity_solver

    type(domain_type), intent(inout) :: &
         domain

    logical, pointer :: &
         config_use_bergs

    call MPAS_pool_get_config(domain % configs, "config_use_bergs", config_use_bergs)

    if (config_use_bergs) then

       call mpas_timer_start("Berg initialization")

       call init_berg_state(domain)

       call init_berg_velocity_solver(domain)

       call mpas_timer_stop("Berg initialization")

    endif

  end subroutine cice_bergs_init

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_run_berg_dynamics
!
!> \brief   Wrapper for iceberg dynamics
!> \author  Darin Comeau, LANL
!> \date    29 Sept 2017
!> \details This subroutine is a wrapper for the iceberg velocity solver
!> over a sea ice time step.
!
!-----------------------------------------------------------------------

  subroutine cice_run_berg_dynamics(domain, clock)!{{{

    use cice_berg_velocity_solver, only: run_berg_velocity_solver

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    ! berg velocity solver
    call mpas_timer_start("Berg velocity solver")
    call run_berg_velocity_solver(domain, clock)
    call mpas_timer_stop("Berg velocity solver")

    ! advection (to add)

  end subroutine cice_run_berg_dynamics!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_run_berg_postdynamics
!
!> \brief   Wrapper for post-dynamics routines
!> \author  Darin Comeau, LANL
!> \date    29 Sept 2017
!> \details This subroutine is a wrapper for the iceberg velocity solver
!> over a sea ice time step.
!
!-----------------------------------------------------------------------

  subroutine cice_run_berg_postdynamics(domain, clock, itimestep)!{{{

    use cice_berg_decay, only: run_berg_decay
    use cice_berg_state, only: berg_calving_event_check 

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(in) :: &
         itimestep !< Input:

    ! berg decay solver
    call mpas_timer_start("Berg decay")
    call run_berg_decay(domain, clock)
    call mpas_timer_stop("Berg decay")

    ! berg calving check
    call mpas_timer_start("Berg calving")
    call berg_calving_event_check(domain, itimestep)
    call mpas_timer_stop("Berg calving")

  end subroutine cice_run_berg_postdynamics!}}}  

!-----------------------------------------------------------------------

end module cice_bergs
