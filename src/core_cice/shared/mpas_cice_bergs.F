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

    use cice_berg_calving, only: &
         cice_init_berg_calving_fluxes

    use cice_berg_velocity_solver, only: &
         cice_init_berg_velocity_solver

    use cice_berg_advection, only: &
         cice_init_berg_tracers

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    logical, pointer :: &
         config_use_bergs ! flag to turn on / off iceberg model

    call MPAS_pool_get_config(domain % configs, "config_use_bergs", config_use_bergs)

    if (config_use_bergs) then

       call mpas_timer_start("Berg initialization")

       ! initialize calving fluxes
       call cice_init_berg_calving_fluxes(domain)

       ! initialize velocity solver
       call cice_init_berg_velocity_solver(domain)

       ! initialize tracers
       call cice_init_berg_tracers(domain)

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

  subroutine cice_run_berg_dynamics(domain, clock, ierr)

    use cice_berg_velocity_solver, only: &
         cice_run_berg_velocity_solver

    use cice_berg_advection, only: &
         cice_run_berg_advection

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(inout) :: &
         ierr  !< Input/Output:

    logical, pointer :: &
         config_use_berg_velocity_solver ! indicates if berg velocity solver is used

    logical, pointer :: &
         config_use_berg_advection ! flag for advection

    ! determine if velocity solver switched on
    call MPAS_pool_get_config(domain % configs, "config_use_berg_velocity_solver", config_use_berg_velocity_solver)
    call MPAS_pool_get_config(domain % configs, "config_use_berg_advection", config_use_berg_advection)

    ! berg velocity solver
    call mpas_timer_start("Berg velocity solver")
    if (config_use_berg_velocity_solver) call cice_run_berg_velocity_solver(domain, clock)
    call mpas_timer_stop("Berg velocity solver")

    ! advection by incremental remapping
    call mpas_timer_start("Berg advection")
    if (config_use_berg_advection) call cice_run_berg_advection(domain, clock, ierr)
    call mpas_timer_stop("Berg advection")

  end subroutine cice_run_berg_dynamics

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

  subroutine cice_run_berg_postdynamics(domain, clock, itimestep)

    use cice_berg_decay, only: &
         cice_run_berg_decay

    use cice_berg_state, only: &
         cice_berg_update_mass_and_area, &
         cice_berg_accumulate_calving_mass, &
         cice_berg_release_calving_mass

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(in) :: &
         itimestep !< Input:

    ! berg mass and area state variables update
    call mpas_timer_start("Berg update mass")
    call cice_berg_update_mass_and_area(domain)
    call mpas_timer_stop("Berg update mass")

    ! berg decay solver
    call mpas_timer_start("Berg decay")
    call cice_run_berg_decay(domain, clock)
    call mpas_timer_stop("Berg decay")

    ! berg calving
    call mpas_timer_start("Berg calving")
    call cice_berg_accumulate_calving_mass(domain)
    call cice_berg_release_calving_mass(domain, itimestep)
    call mpas_timer_stop("Berg calving")

  end subroutine cice_run_berg_postdynamics

!-----------------------------------------------------------------------

end module cice_bergs
