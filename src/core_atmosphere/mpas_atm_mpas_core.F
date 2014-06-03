! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_core

   use mpas_framework

   type (io_output_object), save :: restart_obj
   type (io_input_object), save :: sfc_update_obj
   integer :: current_outfile_frames
   
   type (MPAS_Clock_type) :: clock

   integer, parameter :: outputAlarmID = 1 
   integer, parameter :: restartAlarmID = 2 
   integer, parameter :: sfcAlarmID = 3 
   integer, parameter :: hifreqAlarmID = 4

   contains


   subroutine mpas_core_init(domain, startTimeStamp)

      use mpas_configure
      use mpas_kind_types
      use mpas_grid_types

      implicit none

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp

      real (kind=RKIND) :: dt
      type (block_type), pointer :: block

      character(len=StrKIND) :: timeStamp
      integer :: i
      integer :: ierr

      if (.not. config_do_restart) then

         ! Code that was previously in atm_setup_test_case()

         block => domain % blocklist
         do while (associated(block))
            do i=2,nTimeLevs
               call mpas_copy_state(block % state % time_levs(i) % state, block % state % time_levs(1) % state)
            end do
            block => block % next
         end do

      end if


      !
      ! Initialize core
      !
      dt = config_dt

      call atm_simulation_clock_init(domain, dt, startTimeStamp)

      call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(1) % state % u)

      block => domain % blocklist
      do while (associated(block))
         call atm_mpas_init_block(domain % dminfo, block, block % mesh, dt)
         block % state % time_levs(1) % state % xtime % scalar = startTimeStamp
         block => block % next
      end do

      call mpas_dmpar_exch_halo_field(domain % blocklist % diag % pv_edge)
      call mpas_dmpar_exch_halo_field(domain % blocklist % diag % ru)
      call mpas_dmpar_exch_halo_field(domain % blocklist % diag % rw)

      current_outfile_frames = 0

      if (config_sfc_update_interval /= "none") then

         sfc_update_obj % filename = trim(config_sfc_update_name)
         sfc_update_obj % stream = STREAM_SFC

         call mpas_io_input_init(sfc_update_obj, domain % blocklist, domain % dminfo)

         !     
         ! We need to decide which time slice to read from the surface file - read the most recent time slice that falls before or on the start time
         !
         sfc_update_obj % time = MPAS_seekStream(sfc_update_obj % io_stream, trim(startTimeStamp), MPAS_STREAM_LATEST_BEFORE, timeStamp, ierr)
         if (ierr == MPAS_IO_ERR) then
            write(0,*) 'Error: surface update file '//trim(sfc_update_obj % filename)//' did not contain any times at or before '//trim(startTimeStamp)
            call mpas_dmpar_abort(domain % dminfo)
         end if

         write(0,*) 'Starting model with surface time ', trim(timeStamp)

      end if

   end subroutine mpas_core_init


   subroutine atm_simulation_clock_init(domain, dt, startTimeStamp)

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(out) :: startTimeStamp

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      integer :: ierr

      if(trim(config_start_time) == 'file') then
         open(22,file=trim(config_restart_timestamp_name),form='formatted',status='old')
         read(22,*) startTimeStamp
         close(22)
      else
        startTimeStamp = config_start_time
      end if
      call mpas_set_time(curr_time=startTime, dateTimeString=startTimeStamp, ierr=ierr)
      call mpas_set_timeInterval(timeStep, dt=dt, ierr=ierr)

      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=ierr)
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=ierr)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
            if(startTime + runduration /= stopTime) then
               write(0,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=ierr)
      else
          write(0,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          call mpas_dmpar_abort(domain % dminfo)
      end if

      ! set output alarm
      call mpas_set_timeInterval(alarmTimeStep, timeString=config_output_interval, ierr=ierr)
      alarmStartTime = startTime + alarmTimeStep
      call mpas_add_clock_alarm(clock, outputAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)

      ! set restart alarm, if necessary
      if (trim(config_restart_interval) /= "none") then
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_restart_interval, ierr=ierr)
         alarmStartTime = startTime + alarmTimeStep
         call mpas_add_clock_alarm(clock, restartAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      end if

      ! set high-frequency history alarm, if necessary
      if (trim(config_hifreq_output_interval) /= "none") then
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_hifreq_output_interval, ierr=ierr)
         alarmStartTime = startTime + alarmTimeStep
         call mpas_add_clock_alarm(clock, hifreqAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      end if

      ! set sfc alarm, if necessary
      if (trim(config_sfc_update_interval) /= "none") then      
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_sfc_update_interval, ierr=ierr)
         alarmStartTime = startTime
         call mpas_add_clock_alarm(clock, sfcAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      end if
      
      !TODO: set phyics alarms here...
      !....
      !....

      call mpas_get_time(curr_time=startTime, dateTimeString=startTimeStamp, ierr=ierr)

   end subroutine atm_simulation_clock_init


   subroutine atm_mpas_init_block(dminfo, block, mesh, dt)
   
      use mpas_grid_types
   !   use atm_advection
      use atm_time_integration
      use mpas_configure
      use mpas_rbf_interpolation
      use mpas_vector_reconstruction
#ifdef DO_PHYSICS
!     use mpas_atmphys_aquaplanet
      use mpas_atmphys_control
      use mpas_atmphys_init
      use mpas_atmphys_manager
#endif
   
      implicit none
   
      type (dm_info), intent(in) :: dminfo
      type (block_type), intent(inout) :: block
      type (mesh_type), intent(inout) :: mesh
      real (kind=RKIND), intent(in) :: dt
   
      if (.not. config_do_restart .or. (config_do_restart .and. config_do_DAcycling)) then
         call atm_init_coupled_diagnostics( block % state % time_levs(1) % state, block % diag, mesh)
      end if
      call atm_compute_solve_diagnostics(dt, block % state % time_levs(1) % state, block % diag, mesh)

      call mpas_rbf_interp_initialize(mesh)
      call mpas_init_reconstruct(mesh)
      call mpas_reconstruct(mesh, block % state % time_levs(1) % state % u % array, &
                            block % diag % uReconstructX % array,                   &
                            block % diag % uReconstructY % array,                   &
                            block % diag % uReconstructZ % array,                   &
                            block % diag % uReconstructZonal % array,               &
                            block % diag % uReconstructMeridional % array           &
                           )
   
   !
   ! Note: The following initialization calls have been moved to mpas_setup_test_case()
   !       since values computed by these routines are needed to produce initial fields
   !
   !   call atm_initialize_advection_rk(mesh)
   !   call atm_initialize_deformation_weights(mesh)

#ifdef DO_PHYSICS
      !check that all the physics options are correctly defined and that at least one physics
      !parameterization is called (using the logical moist_physics):
      call physics_namelist_check

      !proceed with initialization of physics parameterization if moist_physics is set to true:
      if(moist_physics) then
         !initialization of seom input variables in registry:
         call physics_registry_init(config_do_restart, mesh, block % sfc_input)
         call physics_run_init(mesh,block % state % time_levs(1) % state,clock)

         !initialization of all physics:
         call physics_init(dminfo, clock, config_do_restart, mesh, block % diag, &
                           block % state % time_levs(1) % state, &
                           block % state % time_levs(1) % state, &
                           block % diag_physics, block % atm_input, &
                           block % sfc_input)

      endif
#endif
   
      call atm_compute_mesh_scaling(mesh)

      call atm_compute_damping_coefs(mesh)

      call atm_compute_pgf_coefs(mesh)

      write(0,*) 'min/max of meshScalingDel2 = ', minval(mesh % meshScalingDel2 % array(1:mesh%nEdges)), &
                                                  maxval(mesh % meshScalingDel2 % array(1:mesh%nEdges))
      write(0,*) 'min/max of meshScalingDel4 = ', minval(mesh % meshScalingDel4 % array(1:mesh%nEdges)), &
                                                  maxval(mesh % meshScalingDel4 % array(1:mesh%nEdges))

      call atm_adv_coef_compression(mesh)
   
   end subroutine atm_mpas_init_block
   
   
   subroutine mpas_core_run(domain, output_obj, output_frame)
   
      use mpas_grid_types
      use mpas_kind_types
      use mpas_io_output
      use mpas_io_input, only: insert_string_suffix => mpas_insert_string_suffix
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
      integer, intent(inout) :: output_frame
   
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: itimestep
      integer :: ierr

      ! For high-frequency diagnostics output
      character (len=StrKIND) :: tempfilename

      ! Eventually, dt should be domain specific
      dt = config_dt

      call atm_write_output_frame(output_obj, output_frame, domain)

      if (trim(config_hifreq_output_interval) /= 'none') then
         block_ptr => domain % blocklist
         call atm_compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, &
                                        block_ptr % diag_physics, block_ptr % mesh, &
                                        block_ptr % tend_physics, block_ptr % tend)
         !call atm_compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, &
         !                                    block_ptr % diag_physics, block_ptr % mesh)

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)

         call insert_string_suffix('diagnostics.nc', trim(timeStamp), tempfilename)
         call write_hifreq_output(block_ptr, tempfilename)
      end if

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 1
      do while (.not. mpas_is_clock_stop_time(clock))

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
         write(0,*) 'Begin timestep ', trim(timeStamp)
         
         !nick 23-02-2014: track vortex
         block_ptr => domain % blocklist
         !call atm_compute_vortexLoc(block_ptr % state % time_levs(1) % state, block_ptr % diag, &
         !                               block_ptr % mesh, domain)

         ! Input external updates (i.e. surface)
         if (mpas_is_alarm_ringing(clock, sfcAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, sfcAlarmID, ierr=ierr)

            call mpas_read_and_distribute_fields(sfc_update_obj)
            sfc_update_obj % time = sfc_update_obj % time + 1
         end if

         call mpas_timer_start("time integration")
         call atm_do_timestep(domain, dt, itimestep)
         call mpas_timer_stop("time integration")   

         ! Move time level 2 fields back into time level 1 for next time step
         call mpas_shift_time_levels_state(domain % blocklist % state)
         

         ! Advance clock before writing output
         itimestep = itimestep + 1
         call mpas_advance_clock(clock)
         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         

         !TODO: MPAS_getClockRingingAlarms is probably faster than multiple MPAS_isAlarmRinging...
         if (mpas_is_alarm_ringing(clock, outputAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, outputAlarmID, ierr=ierr)
            ! output_frame will always be > 1 here unless it was reset after the maximum number of frames per outfile was reached
            if(output_frame == 1) then
               call mpas_output_state_finalize(output_obj, domain % dminfo)
               call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(timeStamp))
            end if
            call atm_write_output_frame(output_obj, output_frame, domain)
         end if

         if (mpas_is_alarm_ringing(clock, hifreqAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, hifreqAlarmID, ierr=ierr)

            block_ptr => domain % blocklist
            call atm_compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, &
                                        block_ptr % diag_physics, block_ptr % mesh, &
                                        block_ptr % tend_physics, block_ptr % tend)
            !call atm_compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, &
            !                                    block_ptr % diag_physics, block_ptr % mesh)

            call insert_string_suffix('diagnostics.nc', trim(timeStamp), tempfilename)
            call write_hifreq_output(block_ptr, tempfilename)
         end if

         if (mpas_is_alarm_ringing(clock, restartAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, restartAlarmID, ierr=ierr)

            block_ptr => domain % blocklist
            do while (associated(block_ptr))
               call atm_compute_restart_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, block_ptr % mesh)
               block_ptr => block_ptr % next
            end do

            ! Write one restart time per file
            call mpas_output_state_init(restart_obj, domain, "RESTART", trim(timeStamp))
            call mpas_output_state_for_domain(restart_obj, domain, 1)
            call mpas_output_state_finalize(restart_obj, domain % dminfo)
         end if

      end do
   
   end subroutine mpas_core_run
   
   
   subroutine atm_write_output_frame(output_obj, output_frame, domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain and write model state to output file
   !
   ! Input/Output: domain - contains model state; diagnostic field are computed
   !                        before returning
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
      use mpas_io_output
   
      implicit none
   
      integer, intent(inout) :: output_frame
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
   
      type (block_type), pointer :: block_ptr
   
      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         call atm_compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, &
                                        block_ptr % diag_physics, block_ptr % mesh, &
                                        block_ptr % tend_physics, block_ptr % tend)
         !call atm_compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, &
         !                                    block_ptr % diag_physics, block_ptr % mesh)
         block_ptr => block_ptr % next
      end do
   
      call mpas_output_state_for_domain(output_obj, domain, output_frame)
      output_frame = output_frame + 1
   
      ! reset frame if the maximum number of frames per outfile has been reached
      if (config_frames_per_outfile > 0) then
         current_outfile_frames = current_outfile_frames + 1            
         if(current_outfile_frames >= config_frames_per_outfile) then
            current_outfile_frames = 0
            output_frame = 1
         end if
      end if

   end subroutine atm_write_output_frame
   
   
   !change atm_compute_output_diagnostics() function argument list in core.F
   !subroutine atm_compute_output_diagnostics(state, diag, diag_physics, grid)
   subroutine atm_compute_output_diagnostics(state, diag, diag_physics, grid, tend_physics, tend)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain to be written to history files
   !
   ! Input: state - contains model prognostic fields
   !        grid  - contains grid metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
      use mpas_constants
      use mpas_atm_interp_diagnostics
   
      implicit none
   
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      type (diag_physics_type), intent(inout) :: diag_physics
      type (mesh_type), intent(in) :: grid
      type (tend_physics_type), intent(inout) :: tend_physics
      type (tend_type), intent(in) :: tend
   
      integer :: iCell, k
      real(kind=RKIND) :: pvuVal, missingVal

      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1._RKIND + rvord * state % scalars % array(state % index_qv,k,iCell))
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * grid % zz % array(k,iCell)
         end do
      end do
      
      call interp_diagnostics(grid,state,diag,diag_physics)
      
      write(0,*) "Calculating ertel pv"
      !need halo cells for everything w/ horizontal derivative
      call mpas_dmpar_exch_halo_field(diag % theta)
      call mpas_dmpar_exch_halo_field(diag % uReconstructX)
      call mpas_dmpar_exch_halo_field(diag % uReconstructY)
      call mpas_dmpar_exch_halo_field(diag % uReconstructZ)
      call mpas_dmpar_exch_halo_field(state % w)
      
      call calc_ertel_pv(grid, state, diag)
      
      pvuVal = 2.0_RKIND
      missingVal = -99999.0_RKIND
      call interp_pv_diagnostics(grid, diag, pvuVal, missingVal)
      !call interp_DT_floodFill(grid, diag, pvuVal, missingVal)
      !call interp_pv(grid%nCells, grid%nVertLevels, pvuVal, grid % latCell % array, &
      !               diag % ertel_pv % array, diag % theta % array, diag % theta_pv % array, &
      !               missingVal)
                     
      write(0,*) "Calculating dpv_dt"
      
      !Update tendencies in halo cells
      !call mpas_timer_start("halo updates for depv_dt")
      call mpas_dmpar_exch_halo_field(tend_physics % rthratenlw)
      call mpas_dmpar_exch_halo_field(tend_physics % rthratensw)
      call mpas_dmpar_exch_halo_field(tend_physics % rthcuten)
      call mpas_dmpar_exch_halo_field(tend_physics % rthblten)
      call mpas_dmpar_exch_halo_field(diag % dtheta_dt_mp)
      call mpas_dmpar_exch_halo_field(tend % theta_euler)
      
      !call mpas_timer_stop("halo updates for depv_dt")
      !add to registry: diag % depv_dt_lw
      call calc_ertel_dpv_dt(grid, state, diag, tend_physics, tend, diag_physics)
      
   end subroutine atm_compute_output_diagnostics
   
   
   subroutine atm_compute_restart_diagnostics(state, diag, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain to be written to restart files
   !
   ! Input: state - contains model prognostic fields
   !        grid  - contains grid metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
      use mpas_constants
   
      implicit none
   
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      type (mesh_type), intent(in) :: grid
   
      integer :: iCell, k

      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1._RKIND + rvord * state % scalars % array(state % index_qv,k,iCell))
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * grid % zz % array(k,iCell)
         end do
      end do
   
   end subroutine atm_compute_restart_diagnostics


   subroutine write_hifreq_output(block_ptr, fname)

      implicit none

      type (block_type), pointer :: block_ptr
      character (len=*), intent(in) :: fname

      integer :: ierr
      type (MPAS_Stream_type) :: hifreq_stream

      call MPAS_createStream(hifreq_stream, trim(fname), MPAS_IO_PNETCDF, MPAS_IO_WRITE, 1, ierr)

      if (block_ptr % mesh % on_a_sphere) then
         call MPAS_writeStreamAtt(hifreq_stream, 'on_a_sphere', 'YES             ', ierr)
      else
         call MPAS_writeStreamAtt(hifreq_stream, 'on_a_sphere', 'NO              ', ierr)
      end if
      call MPAS_writeStreamAtt(hifreq_stream, 'sphere_radius', block_ptr % mesh % sphere_radius, ierr)
      call MPAS_writeStreamAtt(hifreq_stream, 'model_name', block_ptr % domain % modelName, ierr)
      call MPAS_writeStreamAtt(hifreq_stream, 'core_name', block_ptr % domain % coreName, ierr)
      call MPAS_writeStreamAtt(hifreq_stream, 'model_version', block_ptr % domain % modelVersion, ierr)
      call MPAS_writeStreamAtt(hifreq_stream, 'source', 'MPAS', ierr)
      call MPAS_writeStreamAtt(hifreq_stream, 'Conventions', 'MPAS', ierr)

      call MPAS_streamAddField(hifreq_stream, block_ptr % state % time_levs(1) % state % xtime, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag_physics % olrtoa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag_physics % rainc, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag_physics % rainnc, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag_physics % refl10cm_max, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag_physics % precipw, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag_physics % t2m, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % temperature_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % temperature_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % temperature_700hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % temperature_850hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % height_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % height_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % height_850hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % uzonal_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % uzonal_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % uzonal_700hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % uzonal_850hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % umeridional_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % umeridional_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % umeridional_700hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % umeridional_850hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % w_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % w_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % w_850hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % vorticity_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % vorticity_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % vorticity_850hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % theta_pv, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % u_pv, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % v_pv, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % vort_pv, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % depv_dt_lw, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % depv_dt_sw, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % depv_dt_bl, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % depv_dt_cu, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % depv_dt_mp, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % depv_dt_mix, ierr)       
      call MPAS_writeStream(hifreq_stream, 1, ierr)
      call MPAS_closeStream(hifreq_stream, ierr)

   end subroutine write_hifreq_output
   
   
   subroutine atm_do_timestep(domain, dt, itimestep)
   
      use mpas_grid_types
      use mpas_kind_types
      use atm_time_integration
#ifdef DO_PHYSICS
      use mpas_atmphys_control
      use mpas_atmphys_driver
      use mpas_atmphys_manager
      use mpas_atmphys_update
#endif
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: itimestep
      
      type (MPAS_Time_Type) :: startTime, currTime
      type (MPAS_TimeInterval_Type) :: xtimeTime
      character(len=StrKIND) :: timeStamp
      integer :: s, s_n, s_d
      real (kind=RKIND) :: xtime_s
      integer :: ierr

      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         
      xtimeTime = currTime - startTime
      call mpas_get_timeInterval(interval=xtimeTime, S=s, S_n=s_n, S_d=s_d, ierr=ierr)         
      xtime_s = (s + s_n / s_d)

      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         


#ifdef DO_PHYSICS
      !proceed with physics if moist_physics is set to true:
      if(moist_physics) then
         call physics_timetracker(domain,dt,clock,itimestep,xtime_s)
         call physics_driver(domain,itimestep,xtime_s)
      endif
#endif

      call atm_timestep(domain, dt, timeStamp, itimestep)

   end subroutine atm_do_timestep
   
   
   subroutine mpas_core_finalize(domain)
   
      use mpas_grid_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      integer :: ierr

      if (config_sfc_update_interval /= "none") call mpas_io_input_finalize(sfc_update_obj, domain % dminfo)

      call mpas_destroy_clock(clock, ierr)
   
   end subroutine mpas_core_finalize


   subroutine atm_compute_mesh_scaling(mesh)

      use mpas_grid_types

      implicit none

      type (mesh_type), intent(inout) :: mesh

      integer :: iEdge, cell1, cell2
      real (kind=RKIND), dimension(:), pointer :: meshDensity, meshScalingDel2, meshScalingDel4

      meshDensity => mesh % meshDensity % array
      meshScalingDel2 => mesh % meshScalingDel2 % array
      meshScalingDel4 => mesh % meshScalingDel4 % array

      !
      ! Compute the scaling factors to be used in the del2 and del4 dissipation
      !
      meshScalingDel2(:) = 1.0
      meshScalingDel4(:) = 1.0
      if (config_h_ScaleWithMesh) then
         do iEdge=1,mesh%nEdges
            cell1 = mesh % cellsOnEdge % array(1,iEdge)
            cell2 = mesh % cellsOnEdge % array(2,iEdge)
            meshScalingDel2(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)**0.5
            meshScalingDel4(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)
         end do
      end if

   end subroutine atm_compute_mesh_scaling


   subroutine atm_compute_damping_coefs(mesh)

      use mpas_grid_types
      use mpas_configure

      implicit none

      type (mesh_type), intent(inout) :: mesh

      integer :: iCell, k
      real (kind=RKIND) :: z, zt, m1, pii
      real (kind=RKIND), dimension(:,:), pointer :: dss, zgrid

      m1 = -1.0
      pii = acos(m1)

      dss => mesh % dss % array
      zgrid => mesh % zgrid % array

      dss(:,:) = 0.0
      do iCell=1,mesh%nCells
         zt = zgrid(mesh%nVertLevels+1,iCell)
         do k=1,mesh%nVertLevels
            z = 0.5*(zgrid(k,iCell) + zgrid(k+1,iCell))
            if (z > config_zd) then
               dss(k,iCell) = config_xnutr*sin(0.5*pii*(z-config_zd)/(zt-config_zd))**2.0
            end if
         end do
      end do

   end subroutine atm_compute_damping_coefs


   subroutine atm_compute_pgf_coefs(mesh)

      use mpas_grid_types
      use mpas_configure

      implicit none

      type (mesh_type), intent(inout) :: mesh

      integer :: iEdge, iCell1, iCell2, k, iCell, nz, nz1
      real (kind=RKIND) :: d1, d2, d3
      real (kind=RKIND), dimension(:,:), pointer :: cpr, cpl, zgrid, pzp, pzm

      cpr   => mesh % cpr % array
      cpl   => mesh % cpl % array
      pzp   => mesh % pzp % array
      pzm   => mesh % pzm % array
      zgrid => mesh % zgrid % array

!**** coefficient arrays for new pressure gradient calculation

      cpr(:,:) = 0.0
      cpl(:,:) = 0.0

      if (config_newpx) then
         do iEdge=1,mesh%nEdges

            iCell1 = mesh % cellsOnEdge % array(1,iEdge)
            iCell2 = mesh % cellsOnEdge % array(2,iEdge)

            d1       = .25*(zgrid(1,iCell2)+zgrid(2,iCell2)-zgrid(1,iCell1)-zgrid(2,iCell1))
            d2       = d1+.5*(zgrid(3,iCell2)-zgrid(1,iCell2))
            d3       = d2+.5*(zgrid(4,iCell2)-zgrid(2,iCell2))
!            cpr(1,iEdge) = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!            cpr(2,iEdge) = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!            cpr(3,iEdge) = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

            cpr(1,iEdge) =  d2/(d2-d1)
            cpr(2,iEdge) = -d1/(d2-d1)
            cpr(3,iEdge) =  0.

            d1       = .25*(zgrid(1,iCell1)+zgrid(2,iCell1)-zgrid(1,iCell2)-zgrid(2,iCell2))
            d2       = d1+.5*(zgrid(3,iCell1)-zgrid(1,iCell1))
            d3       = d2+.5*(zgrid(4,iCell1)-zgrid(2,iCell1))
!            cpl(1,iEdge) = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!            cpl(2,iEdge) = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!            cpl(3,iEdge) = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

            cpl(1,iEdge) =  d2/(d2-d1)
            cpl(2,iEdge) = -d1/(d2-d1)
            cpl(3,iEdge) =  0.

         end do

!         write(6,*) 'cpr1 = ',cpr(1,1),'  cpl1 = ',cpl(1,1)
!         write(6,*) 'cpr2 = ',cpr(2,1),'  cpl2 = ',cpl(2,1)
!         write(6,*) 'cpr3 = ',cpr(3,1),'  cpl3 = ',cpl(3,1)

      else

!        Coefficients for computing vertical pressure gradient dp/dz
!        dp/dz (k,iCell) = pzp(k,iCell) * (p(k+1,iCell) - p(k,iCell)) +pzm(k,iCell) * (p(k,iCell) - p(k-1,iCell))

         nz1 = mesh % nVertLevels
         nz = nz1 + 1

         do iCell=1, mesh % nCells

            d1 = zgrid(3,iCell)-zgrid(1,iCell)
            d2 = zgrid(4,iCell)-zgrid(2,iCell)
            d3 = d1+d2
            pzm(1,iCell) =  2.*d3/(d1*d2)
            pzp(1,iCell) = -2.*d1/(d2*d3)

            do k=2,nz1-1
               pzp(k,iCell) = 2.*(zgrid(k+1,iCell)-zgrid(k-1,iCell))/     &
     &                      ((zgrid(k+2,iCell)-zgrid(k  ,iCell))*     &
     &                       (zgrid(k+2,iCell)-zgrid(k  ,iCell)       &
     &                       +zgrid(k+1,iCell)-zgrid(k-1,iCell)))
               pzm(k,iCell) = 2.*(zgrid(k+2,iCell)-zgrid(k  ,iCell))/     &
     &                      ((zgrid(k+1,iCell)-zgrid(k-1,iCell))*     &
     &                       (zgrid(k+2,iCell)-zgrid(k  ,iCell)       &
     &                       +zgrid(k+1,iCell)-zgrid(k-1,iCell)))
            end do

            pzp(nz1,iCell) = 0.
            pzm(nz1,iCell) = 2./(zgrid(nz,iCell)-zgrid(nz1-1,iCell))

         end do

      end if

   end subroutine atm_compute_pgf_coefs


   subroutine atm_adv_coef_compression( grid )

      implicit none

      type (mesh_type), intent(inout) :: grid


      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, advCellsForEdge
      integer, dimension(:), pointer :: nEdgesOnCell, nAdvCellsForEdge

      integer :: cell1, cell2, iEdge, n, i, j, j_in, iCell
      integer :: cell_list(20), ordered_cell_list(20)
      logical :: addcell

      deriv_two => grid % deriv_two % array
      adv_coefs => grid % adv_coefs % array
      adv_coefs_3rd => grid % adv_coefs_3rd % array
      cellsOnCell => grid % cellsOnCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      advCellsForEdge => grid % advCellsForEdge % array
      nEdgesOnCell => grid % nEdgesOnCell % array
      nAdvCellsForEdge => grid % nAdvCellsForEdge % array

      do iEdge = 1, grid % nEdges
         nAdvCellsForEdge(iEdge) = 0
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         !
         ! do only if this edge flux is needed to update owned cells
         !
         if (cell1 <= grid%nCells .or. cell2 <= grid%nCells) then
 
            cell_list(1) = cell1
            cell_list(2) = cell2
            n = 2 
  
          !  add cells surrounding cell 1.  n is number of cells currently in list
            do i = 1, nEdgesOnCell(cell1)
               if (cellsOnCell(i,cell1) /= cell2) then
                  n = n + 1
                  cell_list(n) = cellsOnCell(i,cell1)
               end if
            end do
  
          !  add cells surrounding cell 2 (brute force approach)
            do iCell = 1, nEdgesOnCell(cell2)
               addcell = .true.
               do i=1,n
                  if (cell_list(i) == cellsOnCell(iCell,cell2)) addcell = .false.
               end do
               if (addcell) then
                  n = n+1
                  cell_list(n) = cellsOnCell(iCell,cell2)
               end if
            end do
  
          ! order the list by increasing cell number (brute force approach)
  
            do i=1,n
               ordered_cell_list(i) = grid % nCells + 2
               j_in = 1
               do j=1,n
                  if (ordered_cell_list(i) > cell_list(j) ) then
                     j_in = j
                     ordered_cell_list(i) = cell_list(j)
                  end if
               end do
!               ordered_cell_list(i) = cell_list(j_in)
               cell_list(j_in) = grid % nCells + 3
            end do
  
            nAdvCellsForEdge(iEdge) = n
            do iCell = 1, nAdvCellsForEdge(iEdge)
               advCellsForEdge(iCell,iEdge) = ordered_cell_list(iCell)
            end do
  
          ! we have the ordered list, now construct coefficients
  
            adv_coefs(:,iEdge) = 0.
            adv_coefs_3rd(:,iEdge) = 0.
          
          ! pull together third and fourth order contributions to the flux
          ! first from cell1
  
            j_in = 0
            do j=1, n
               if( ordered_cell_list(j) == cell1 ) j_in = j
            end do
            adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(1,1,iEdge)
            adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) + deriv_two(1,1,iEdge)
  
            do iCell = 1, nEdgesOnCell(cell1)
               j_in = 0
               do j=1, n
                 if( ordered_cell_list(j) == cellsOnCell(iCell,cell1) ) j_in = j
               end do
               adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(iCell+1,1,iEdge)
               adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) + deriv_two(iCell+1,1,iEdge)
            end do
  
          ! pull together third and fourth order contributions to the flux
          ! now from cell2
  
            j_in = 0
            do j=1, n
               if( ordered_cell_list(j) == cell2 ) j_in = j
            end do
            adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(1,2,iEdge)
            adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) - deriv_two(1,2,iEdge)
  
            do iCell = 1, nEdgesOnCell(cell2)
               j_in = 0
               do j=1, n
                  if( ordered_cell_list(j) == cellsOnCell(iCell,cell2) ) j_in = j
               end do
               adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(iCell+1,2,iEdge)
               adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) - deriv_two(iCell+1,2,iEdge)
            end do
  
            do j = 1,n
               adv_coefs    (j,iEdge) = - (grid % dcEdge % array (iEdge) **2) * adv_coefs    (j,iEdge) / 12.
               adv_coefs_3rd(j,iEdge) = - (grid % dcEdge % array (iEdge) **2) * adv_coefs_3rd(j,iEdge) / 12.
            end do
  
          ! 2nd order centered contribution - place this in the main flux weights
  
            j_in = 0
            do j=1, n
               if( ordered_cell_list(j) == cell1 ) j_in = j
            end do
            adv_coefs(j_in,iEdge) = adv_coefs(j_in,iEdge) + 0.5
  
            j_in = 0
            do j=1, n
               if( ordered_cell_list(j) == cell2 ) j_in = j
            end do
            adv_coefs(j_in,iEdge) = adv_coefs(j_in,iEdge) + 0.5
  
          !  multiply by edge length - thus the flux is just dt*ru times the results of the vector-vector multiply
  
            do j=1,n
               adv_coefs    (j,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs    (j,iEdge)
               adv_coefs_3rd(j,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs_3rd(j,iEdge)
            end do
 
         end if  ! only do for edges of owned-cells
         
      end do ! end loop over edges

   end subroutine atm_adv_coef_compression
   
   subroutine interp_DT_floodFill(grid, diag, pvuVal, missingVal)
      !searching down columns to find 2pvu surface is buggy with stratospheric wave breaking,
      !since it won't find the dynamic tropopause we'd like to find.
      !a quick fix to check if a few levels below is also < 2pvu ends up finding the surface
      !if there's high pv at the surface that gets close to DT.
      !to alleviate the (hopefully) pockets of wave breaking, we can flood fill from a known
      !region (ie model top) and hopefully cascade around any trouble regions.
      !we can also add in communication between blocks in the future if need be.
     
      implicit none
      
      type (mesh_type), intent(in) :: grid
      type (diag_type), intent(inout) :: diag
      real(kind=RKIND) ::  pvuVal, missingVal
      
      integer :: iCell, k, nCells, nVertLevels, nChanged, iNbr, iCellNbr
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex!, inStrato
                                          
      real(kind=RKIND) :: valInterpCell, sgnHemi, sgn, valh, vall, dv_dl, levFrac
      real(kind=RKIND),dimension(:),pointer:: areaCell, latCell, u_pv, v_pv, theta_pv, vort_pv
      real(kind=RKIND),dimension(:,:),pointer:: uZonal, uMeridional, vorticity, theta, ertel_pv, &
                                                kiteAreasOnVertex
                          
      integer, dimension(:,:), allocatable :: candInStrato, inStrato
      real(kind=RKIND), dimension(:,:), allocatable :: vVort
      
      nCells = grid % nCellsSolve 
      nVertLevels = grid % nVertLevels
      !maxEdgesCell = grid % maxEdges
      
      nEdgesOnCell => grid % nEdgesOnCell % array
      cellsOnCell => grid % cellsOnCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      verticesOnCell => grid % verticesOnCell % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnVertex => grid % cellsOnVertex % array
      areaCell => grid % areaCell % array
      
      latCell => grid % latCell % array
      
      uZonal      => diag % uReconstructZonal % array
      uMeridional => diag % uReconstructMeridional % array
      vorticity   => diag % vorticity % array
      theta       => diag % theta % array
      ertel_pv    => diag % ertel_pv % array
      
      u_pv => diag % u_pv % array
      v_pv => diag % v_pv % array
      theta_pv => diag % theta_pv % array
      vort_pv => diag % vort_pv % array
      !inStrato => diag % inStrato % array
      
      !flood fill above dynamic tropopause.
      !candInStrato>0 if above 2pvu. inStrato(k,iCell)>0 if above DT
      allocate(candInStrato(nVertLevels, nCells))
      allocate(inStrato(nVertLevels, nCells))
      candInStrato(:,:) = 0
      inStrato(:,:) = 0
      !store whether candidate above DT to avoid repeating logic.
      do iCell=1,nCells
         sgnHemi = sign(1.0_RKIND, latCell(iCell)) !at the equator, sign(0)=0
         if (sgnHemi .EQ. 0.0) sgnHemi = 1.0
         valInterpCell = pvuVal*sgnHemi
         do k=1,nVertLevels
            sgn = (ertel_pv(k,iCell)-valInterpCell)*sgnHemi
            if (sgn .GE. 0) candInStrato(k,iCell) = 1
         end do
      end do
      
      !seed flood fill with model top that's above DT.
      !can have model top below 2pvu (eg, tropics)
      nChanged = 0
      do iCell=1,nCells
         do k=nVertLevels-5,nVertLevels
            if (candInStrato(k,iCell) .GT. 0) then
               inStrato(k,iCell) = 1
               candInStrato(k,iCell) = 0
               nChanged = nChanged+1
            end if
         end do
      end do
      
      !flood fill from the given seeds. since I don't know enough fortran,
      !we'll just brute force a continuing loop.
      do while(nChanged .GT. 0)
        nChanged = 0
        do iCell=1,nCells
          do k=nVertLevels,1,-1
             !update if candidate and neighbor in strato
             if (candInStrato(k,iCell) .GT. 0) then
                !nbr above
                if (k .LT. nVertLevels) then
                  if (inStrato(k+1,iCell) .GT. 0) then
                    inStrato(k,iCell) = 1
                    candInStrato(k,iCell) = 0
                    nChanged = nChanged+1
                    cycle
                  end if
                end if
                
                !side nbrs
                do iNbr = 1, nEdgesOnCell(iCell)
                  iCellNbr = cellsOnCell(iNbr,iCell)
                  if (inStrato(k,iCellNbr) .GT. 0) then
                    inStrato(k,iCell) = 1
                    candInStrato(k,iCell) = 0
                    nChanged = nChanged+1
                    cycle
                  end if
                end do
                
                !nbr below
                if (k .GT. 1) then
                  if (inStrato(k-1,iCell) .GT. 0) then
                    inStrato(k,iCell) = 1
                    candInStrato(k,iCell) = 0
                    nChanged = nChanged+1
                    cycle
                  end if
                end if
             end if !candInStrato
          end do !levels
        end do !cells
      end do !while
      
      deallocate(candInStrato)
      
      !interpolate fields to the bottom of this strato
      allocate(vVort(nVertLevels, nCells))
      do iCell=1,nCells
         do k=1,nVertLevels
            vVort(k,iCell) = calc_verticalVorticity_cell(iCell, k, nEdgesOnCell(iCell), verticesOnCell, cellsOnVertex, &
                                                         kiteAreasOnVertex, areaCell(iCell), vorticity)
         end do
      end do
      
      u_pv(:)= missingVal
      v_pv(:)= missingVal
      theta_pv(:)= missingVal
      vort_pv(:)= missingVal
      do iCell=1,nCells
         do k=1,nVertLevels
            !find lowest strato level
            if (inStrato(k,iCell) .GT. 0) exit
         end do
         if (k .EQ. nVertLevels) cycle
         
         !for interpolation to 2 pvu,
         !take value if at lowest value. linear if above.
         if (k .EQ. 1) then
           u_pv(iCell) = uZonal(k,iCell)
           v_pv(iCell) = uMeridional(k,iCell)
           theta_pv(iCell) = theta(k,iCell)
           vort_pv(iCell) = vVort(k,iCell)
         else
           !interpolate between this level and beneath
           k = k-1
           valh = ertel_pv(k+1,iCell)
           vall = ertel_pv(k,iCell)
           dv_dl = valh-vall
           if (abs(dv_dl)<1.e-6) then
             levFrac = .5
           else
             sgnHemi = sign(1.0_RKIND, latCell(iCell)) !at the equator, sign(0)=0
             if (sgnHemi .EQ. 0.0) sgnHemi = 1.0
             valInterpCell = pvuVal*sgnHemi
             levFrac = (valInterpCell-vall)/dv_dl
           end if
           
           valh = uZonal(k+1,iCell)
           vall = uZonal(k,iCell)
           dv_dl = valh-vall
           u_pv(iCell) = vall+dv_dl*levFrac
           
           valh = uMeridional(k+1,iCell)
           vall = uMeridional(k,iCell)
           dv_dl = valh-vall
           v_pv(iCell) = vall+dv_dl*levFrac
           
           valh = theta(k+1,iCell)
           vall = theta(k,iCell)
           dv_dl = valh-vall
           theta_pv(iCell) = vall+dv_dl*levFrac
           
           valh = vVort(k+1,iCell)
           vall = vVort(k,iCell)
           dv_dl = valh-vall
           vort_pv(iCell) = vall+dv_dl*levFrac
           
         end if !column interp
      
      end do !nCells
      
      deallocate(vVort)
      deallocate(inStrato)
      
   end subroutine interp_DT_floodFill
   
   subroutine interp_pv_diagnostics(grid, diag, pvuVal, missingVal)
      !compute various fields on 2pvu surface using calculated PVU field
      !potential temperature, uZonal, uMeridional, vertical vorticity
      
      implicit none
      
      type (mesh_type), intent(in) :: grid
      type (diag_type), intent(inout) :: diag
      real(kind=RKIND) ::  pvuVal, missingVal
      
      integer :: iCell, k, nCells, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex
                                          
      real(kind=RKIND),dimension(:),pointer:: areaCell, latCell, u_pv, v_pv, theta_pv, vort_pv
      real(kind=RKIND),dimension(:,:),pointer:: uZonal, uMeridional, vorticity, theta, ertel_pv, &
                                                kiteAreasOnVertex
      
      real(kind=RKIND), dimension(:,:), allocatable :: vVort
      
      nCells = grid % nCellsSolve 
      nVertLevels = grid % nVertLevels
      !maxEdgesCell = grid % maxEdges
      
      nEdgesOnCell => grid % nEdgesOnCell % array
      cellsOnCell => grid % cellsOnCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      verticesOnCell => grid % verticesOnCell % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnVertex => grid % cellsOnVertex % array
      areaCell => grid % areaCell % array
      
      latCell => grid % latCell % array
      
      uZonal      => diag % uReconstructZonal % array
      uMeridional => diag % uReconstructMeridional % array
      vorticity   => diag % vorticity % array
      theta       => diag % theta % array
      ertel_pv    => diag % ertel_pv % array
      
      u_pv => diag % u_pv % array
      v_pv => diag % v_pv % array
      theta_pv => diag % theta_pv % array
      vort_pv => diag % vort_pv % array
      
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, uZonal, u_pv, missingVal)
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, uMeridional, v_pv, missingVal)
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, theta, theta_pv, missingVal)
                     
      allocate(vVort(nVertLevels, nCells))
      do iCell=1,nCells
         do k=1,nVertLevels
            vVort(k,iCell) = calc_verticalVorticity_cell(iCell, k, nEdgesOnCell(iCell), verticesOnCell, cellsOnVertex, &
                                                         kiteAreasOnVertex, areaCell(iCell), vorticity)
         end do
      end do
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, vVort, vort_pv, missingVal)
      deallocate(vVort)
      
   end subroutine interp_pv_diagnostics
   
   subroutine calc_ertel_dpv_dt(grid, state, diag, tend_physics, tend, diag_physics)
      !Driver for Cavallo's way of doing PV budgets.
      !eg Pedlosky Geophysical Fluid Dynamics, Depv_Dt = gradxu.gradDtheta_Dt/density + FrictionTerm
      !is linear in the Dtheta_Dt physics forcings so we can piece out the contributions of individual physics params.
   
      !this inefficient version is straightforward in that we'll recalculate vorticity and such
      !for every component.
      !This function depends on adding depv_dt* terms to registry.
      !Since the default microphysics tendency rt_diabatic_tend is for theta_m, 
      !track dtheta_dt_microphysics in src/core_atmosphere/physics/mpas_atmphys_interface_nhyd.F:microphysics_to_MPAS() 
      !through th_p.
      !Could also add dqv_micro to registry, and storing the change in qv from
      !the microphysics scheme over a timestep (do in core_atmosphere/dynamics/mpas/mpas_atm_time_integration.f90: call microphysics_driver()).
      
      use mpas_constants !for rvord
      
      implicit none

      type (mesh_type), intent(in) :: grid
      type (state_type), intent(in) :: state
      type (diag_type), intent(inout) :: diag
      type(tend_physics_type),intent(in):: tend_physics
      type (tend_type), intent(in) :: tend
      type(diag_physics_type),intent(in):: diag_physics

      integer :: iCell, k
      integer :: nCellsSolve, nVertLevels, nCells
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex
      !real(kind=RKIND) :: rvord
      real(kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real(kind=RKIND), dimension(:,:), pointer :: u, w, rho, vorticity, zgrid, &
                                                   localVerticalUnitVectors, edgeNormalVectors, kiteAreasOnVertex, &
                                                   theta, uReconstructX, uReconstructY, uReconstructZ
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_lw, depv_dt_sw, depv_dt_bl, depv_dt_cu, depv_dt_mp, &
                                                   depv_dt_mix
      real(kind=RKIND), dimension(:,:), pointer :: rthblten, rthcuten, rthratenlw, rthratensw, &
                                                   dtheta_dt_mp, tend_theta_euler, qv
      real(kind=RKIND), dimension(:,:,:), pointer :: cellTangentPlane
      
      real(kind=RKIND), dimension(:,:), allocatable :: thetaTendency

      nCellsSolve = grid % nCellsSolve 
      nVertLevels = grid % nVertLevels
      nCells = grid % nCells
      !maxEdgesCell = grid % maxEdges
      
      nEdgesOnCell => grid % nEdgesOnCell % array
      cellsOnCell => grid % cellsOnCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      edgesOnCell => grid % edgesOnCell % array
      verticesOnCell => grid % verticesOnCell % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnVertex => grid % cellsOnVertex % array
      dvEdge => grid % dvEdge % array
      areaCell => grid % areaCell % array
      !ertel_pv => diag % ertel_pv % array
      u => state % u % array
      w => state % w % array
      theta => diag % theta % array
      rho => diag % rho % array
      vorticity    => diag % vorticity % array
      uReconstructX => diag % uReconstructX % array
      uReconstructY => diag % uReconstructY % array
      uReconstructZ => diag % uReconstructZ % array
      cellTangentPlane => grid % cellTangentPlane % array
      localVerticalUnitVectors => grid % localVerticalUnitVectors % array
      edgeNormalVectors => grid % edgeNormalVectors % array
      zgrid => grid % zgrid % array
      
      depv_dt_lw => diag % depv_dt_lw % array
      depv_dt_sw => diag % depv_dt_sw % array
      depv_dt_bl => diag % depv_dt_bl % array
      depv_dt_cu => diag % depv_dt_cu % array
      depv_dt_mp => diag % depv_dt_mp % array
      depv_dt_mix => diag % depv_dt_mix % array

      rthblten   => tend_physics % rthblten   % array
      rthcuten   => tend_physics % rthcuten   % array
      rthratenlw => tend_physics % rthratenlw % array
      rthratensw => tend_physics % rthratensw % array
      dtheta_dt_mp => diag % dtheta_dt_mp % array
      tend_theta_euler => tend % theta_euler % array
      qv => state % scalars % array(state%index_qv,:,:)
      
     !Physics schemes ----------------------------- 
     !radiation ------
     !lw
     do iCell=1,nCellsSolve
         do k=1,nVertLevels
            depv_dt_lw(k,iCell) = driver_ertelPV_cell_reconstruct(iCell, k, nVertLevels, &
                                                 nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                                 cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                                 cellsOnVertex, &
                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                                 uReconstructX, uReconstructY, uReconstructZ, w,rthratenlw,rho(k,iCell),vorticity)

        end do 
     end do
     
     !sw
     do iCell=1,nCellsSolve
         do k=1,nVertLevels
            depv_dt_sw(k,iCell) = driver_ertelPV_cell_reconstruct(iCell, k, nVertLevels, &
                                                 nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                                 cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                                 cellsOnVertex, &
                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                                 uReconstructX, uReconstructY, uReconstructZ, w,rthratensw,rho(k,iCell),vorticity)

        end do 
     end do
     
     !boundary layer -----------
     do iCell=1,nCellsSolve
         do k=1,nVertLevels
            depv_dt_bl(k,iCell) = driver_ertelPV_cell_reconstruct(iCell, k, nVertLevels, &
                                                 nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                                 cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                                 cellsOnVertex, &
                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                                 uReconstructX, uReconstructY, uReconstructZ, w,rthblten,rho(k,iCell),vorticity)

        end do 
     end do
     
     !convection scheme --------
     do iCell=1,nCellsSolve
         do k=1,nVertLevels
            depv_dt_cu(k,iCell) = driver_ertelPV_cell_reconstruct(iCell, k, nVertLevels, &
                                                 nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                                 cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                                 cellsOnVertex, &
                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                                 uReconstructX, uReconstructY, uReconstructZ, w,rthcuten,rho(k,iCell),vorticity)

        end do 
     end do
     
     !microphysics -------
     !requires adding to /src/core_atmosphere/physics/mpas_atmphys_driver_microphysics.F:microphysics_driver()
     !to get dtheta/dt instead of dtheta_m/dt
     do iCell=1,nCellsSolve
         do k=1,nVertLevels
            depv_dt_mp(k,iCell) = driver_ertelPV_cell_reconstruct(iCell, k, nVertLevels, &
                                                 nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                                 cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                                 cellsOnVertex, &
                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                                 uReconstructX, uReconstructY, uReconstructZ, w,dtheta_dt_mp,rho(k,iCell),vorticity)

        end do 
     end do
     
     !mixing -----------------------
     
     allocate(thetaTendency(nVertLevels,nCells))
     
     do iCell=1,nCells
        do k=1,nVertLevels
          !if qv not mixed
          thetaTendency(k,iCell) = tend_theta_euler(k,iCell)/( 1._RKIND+rvord*qv(k,iCell) )
        end do
     end do
     
     do iCell=1,nCellsSolve
         do k=1,nVertLevels
            depv_dt_mix(k,iCell) = driver_ertelPV_cell_reconstruct(iCell, k, nVertLevels, &
                                                 nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                                 cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                                 cellsOnVertex, &
                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                                 uReconstructX, uReconstructY, uReconstructZ, w, thetaTendency,rho(k,iCell),vorticity)

        end do 
     end do
     
     deallocate(thetaTendency)
     
   end subroutine calc_ertel_dpv_dt
   
   subroutine calc_ertel_pv(grid, state, diag)

      implicit none

      type (mesh_type), intent(in) :: grid
      type (state_type), intent(in) :: state
      type (diag_type), intent(inout) :: diag

      integer :: iCell, k
      integer :: nCellsSolve, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex
      real(kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real(kind=RKIND), dimension(:,:), pointer :: ertel_pv, u, w, theta, rho, vorticity, zgrid, &
                                                   localVerticalUnitVectors, edgeNormalVectors, kiteAreasOnVertex, &
                                                   uReconstructX, uReconstructY, uReconstructZ
      real(kind=RKIND), dimension(:,:,:), pointer :: cellTangentPlane


      nCellsSolve = grid % nCellsSolve 
      nVertLevels = grid % nVertLevels
      !maxEdgesCell = grid % maxEdges
      
      nEdgesOnCell => grid % nEdgesOnCell % array
      cellsOnCell => grid % cellsOnCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      edgesOnCell => grid % edgesOnCell % array
      verticesOnCell => grid % verticesOnCell % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnVertex => grid % cellsOnVertex % array
      dvEdge => grid % dvEdge % array
      areaCell => grid % areaCell % array
      ertel_pv => diag % ertel_pv % array
      u => state % u % array
      w => state % w % array
      theta => diag % theta % array
      rho => diag % rho % array
      vorticity    => diag % vorticity % array
      uReconstructX => diag % uReconstructX % array
      uReconstructY => diag % uReconstructY % array
      uReconstructZ => diag % uReconstructZ % array
      cellTangentPlane => grid % cellTangentPlane % array
      localVerticalUnitVectors => grid % localVerticalUnitVectors % array
      edgeNormalVectors => grid % edgeNormalVectors % array
      zgrid => grid % zgrid % array

      do iCell=1,nCellsSolve
         do k=1,nVertLevels
!            ertel_pv(k,iCell) = driver_ertelPV_cell_uv(iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, &
!                                                 kiteAreasOnVertex, &
!                                                 cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
!                                                 cellsOnVertex, &
!                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
!                                                 u, w,theta,rho(k,iCell),vorticity)
            !ertel_pv(k,iCell) = 3.1415925
!            ertel_pv(k,iCell) = driver_ertelPV_cell_approx(iCell, k,nVertLevels, nEdgesOnCell(iCell), verticesOnCell, &
!                                                 kiteAreasOnVertex, &
!                                                 cellsOnCell, edgesOnCell,cellsOnEdge, edgeNormalVectors, &
!                                                 cellsOnVertex, &
!                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
!                                                 theta,rho(k,iCell),vorticity)
            ertel_pv(k,iCell) = driver_ertelPV_cell_reconstruct(iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, &
                                                  kiteAreasOnVertex, &
                                                 cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                                 cellsOnVertex, &
                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                                 uReconstructX, uReconstructY, uReconstructZ, w,theta,rho(k,iCell),vorticity)

         end do 
      end do 

   end subroutine calc_ertel_pv
 
   real(kind=RKIND) function driver_ertelPV_cell_approx(iCell, level, nVertLevels, nEdgesCell0, verticesOnCell, kiteAreasOnVertex, &
                                                 cellsOnCell, edgesOnCell, cellsOnEdge, edgeNormalVectors, &
                                                 cellsOnVertex, &
                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell0, &
                                                 theta,rho0,vorticity)

      !We'll use the local xyz system for each cell since it's well aligned with physical gradients.
      !x,y are the tangent plane and z is the radially out.

      !Instead of dot(gradxu,gradtheta)/density, we'll do
      !verticalVorticity*dtheta_dz/density as a "large scale" approximation
                                                 
      implicit none

      real(kind=RKIND), intent(in) :: rho0, areaCell0
      real(kind=RKIND), dimension(3,2,*), intent(in) :: cellTangentPlane
      real(kind=RKIND), dimension(3,*), intent(in) :: localVerticalUnitVectors, edgeNormalVectors
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid, theta, vorticity, kiteAreasOnVertex
      integer, intent(in) :: iCell, level, nVertLevels, nEdgesCell0
      integer, dimension(:,:), intent(in) :: cellsOnCell, edgesOnCell, cellsOnEdge, verticesOnCell, cellsOnVertex
      
      integer :: i, iNbr, iEdge
      real(kind=RKIND) :: val0, valNbr, volumeCell, areaFactor, z0, zp, zm, valp, valm, dw_dx, dw_dy, du_dz, dv_dz
      real(kind=RKIND), dimension(3) :: gradxu, gradtheta, unitDeriv
      real(kind=RKIND), dimension(3,3) :: xyzLocal
      real(kind=RKIND), dimension(nEdgesCell0) :: valEdges, dvEdgeCell, dhEdge
      real(kind=RKIND), dimension(3,nEdgesCell0) :: normalEdgeCell

      !local coordinate system
      call coordinateSystem_cell(cellTangentPlane, localVerticalUnitVectors, iCell, xyzLocal)
      
      !vertical vorticity to cell center
      gradxu(1) = 0.0_RKIND
      gradxu(2) = 0.0_RKIND
      gradxu(3) = calc_verticalVorticity_cell(iCell, level, nEdgesCell0, verticesOnCell, cellsOnVertex, &
                                              kiteAreasOnVertex, areaCell0, vorticity)
      
      !dtheta_dz
      gradtheta(1) = 0.0_RKIND
      gradtheta(2) = 0.0_RKIND
      gradtheta(3) = 0.0_RKIND
      z0 = calc_heightCellCenter(iCell, level, zgrid)
      val0 = theta(level, iCell)
      if (level>1) then
         !have cell beneath
         valm = theta(level-1, iCell)
         zm = calc_heightCellCenter(iCell, level-1, zgrid)
      end if
      if (level<nVertLevels) then
         !have cell above
         valp = theta(level+1, iCell)
         zp = calc_heightCellCenter(iCell, level+1, zgrid)
      end if

      if (level==1) then
         !calc_vertDeriv_one(valp, valm, dz)
         gradtheta(3) = calc_vertDeriv_one(valp, val0, zp-z0)
      else if (level==nVertLevels) then
         gradtheta(3) = calc_vertDeriv_one(val0, valm, z0-zm)
      else
         gradtheta(3) = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      end if
      
      driver_ertelPV_cell_approx = formErtelPV(gradxu, gradtheta, rho0, xyzLocal(:,1), xyzLocal(:,2), xyzLocal(:,3))
   end function driver_ertelPV_cell_approx

   
   real(kind=RKIND) function driver_ertelPV_cell_reconstruct(iCell, level, nVertLevels, nEdgesCell0, verticesOnCell, kiteAreasOnVertex, &
                                                 cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                                 cellsOnVertex, &
                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell0, &
                                                 uReconstructX, uReconstructY, uReconstructZ, w,theta,rho0,vorticity)

      !We'll use the local xyz system for each cell since it's well aligned with physical gradients.
      !x,y are the tangent plane and z is the radially out.
      !For velocity du/dz and dv/dz, we'll use the reconstructed horizontal
      !velocity at cell centers

      !Full curl is (dw/dy-dv/dz)i + (du/dz-dw/dx)j + (dv/dx-du/dy)k
      !k component is solved for already through circulation to vertices. We'll use
      !the kite_areas weighting as in pv_cell to get it at the cell centers.

      !The other components of vorticity are:
      !(dw/dy-dv/dz)i + (du/dz-dw/dx)j where velocities are in the local coordinate system.
      !We'll also need dTheta/dX, dTheta/dY, dTheta/dZ.
      !For all these, it's finite volume gradients in the
      !horizontal and finite differencing in the vertical.
                                                 
      implicit none

      real(kind=RKIND), intent(in) :: rho0, areaCell0
      real(kind=RKIND), dimension(:), intent(in) :: dvEdge
      real(kind=RKIND), dimension(3,2,*), intent(in) :: cellTangentPlane
      real(kind=RKIND), dimension(3,*), intent(in) :: localVerticalUnitVectors, edgeNormalVectors
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid,uReconstructX, uReconstructY, uReconstructZ, &
                                                      w, theta, vorticity, kiteAreasOnVertex
      integer, intent(in) :: iCell, level, nVertLevels, nEdgesCell0
      integer, dimension(:,:), intent(in) :: cellsOnCell, edgesOnCell, cellsOnEdge, verticesOnCell, cellsOnVertex
      
      integer :: i, iNbr, iEdge
      real(kind=RKIND) :: val0, valNbr, volumeCell, areaFactor, z0, zp, zm, valp, valm, dw_dx, dw_dy, du_dz, dv_dz
      real(kind=RKIND), dimension(3) :: gradxu, gradtheta, unitDeriv, velCell0, velCellp, velCellm
      real(kind=RKIND), dimension(3,3) :: xyzLocal
      real(kind=RKIND), dimension(nEdgesCell0) :: valEdges, dvEdgeCell, dhEdge
      real(kind=RKIND), dimension(3,nEdgesCell0) :: normalEdgeCell

      !local coordinate system
      call coordinateSystem_cell(cellTangentPlane, localVerticalUnitVectors, iCell, xyzLocal)
      !normal vectors at voronoi polygon edges pointing out of cell
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         !dhEdge(i) = calc_heightVerticalEdge(iCell, iNbr, level, zgrid) !vertical thickness of that face
         !if don't want to consider 3d cell since we haven't calculated the cell
         !volume well, set all thicknesses to be the same
         dhEdge(i) = 100.0_RKIND

         iEdge = edgesOnCell(i,iCell)
         dvEdgeCell(i) = dvEdge(iEdge)
         val0 = fluxSign(iCell, iEdge, cellsOnEdge)
         normalEdgeCell(:,i) = edgeNormalVectors(:,iEdge)
         call normalizeVector(normalEdgeCell(:,i),3)
         normalEdgeCell(:,i) = normalEdgeCell(:,i)*val0
      end do

      volumeCell = calcVolumeCell(areaCell0, nEdgesCell0, dhEdge)

      !Need to get 3d curl and grad theta
      !horizontal derivatives
      !calc_horizDeriv_fv(valEdges, nNbrs, dvEdge, dhEdge, &
      !                                         normalEdge, unitDeriv, volumeCell)
      !theta
      val0 = theta(level, iCell)
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         valNbr = theta(level,iNbr)
         valEdges(i) = 0.5*(valNbr+val0)
      end do
      unitDeriv(:) = xyzLocal(:,1)
      gradtheta(1) = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)
      unitDeriv(:) = xyzLocal(:,2)
      gradtheta(2) = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)

      !w
      val0 = .5*(w(level+1, iCell)+w(level, iCell))
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         valNbr = .5*(w(level+1, iNbr)+w(level, iNbr))
         valEdges(i) = 0.5*(valNbr+val0)
      end do
      unitDeriv(:) = xyzLocal(:,1)
      dw_dx = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)
      unitDeriv(:) = xyzLocal(:,2)
      dw_dy = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)

      !vertical derivatives
      !calc_heightCellCenter(c0, level, zgrid) calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      !theta
      gradtheta(3) = 0.0_RKIND
      z0 = calc_heightCellCenter(iCell, level, zgrid)
      val0 = theta(level, iCell)
      if (level>1) then
         !have cell beneath
         valm = theta(level-1, iCell)
         zm = calc_heightCellCenter(iCell, level-1, zgrid)
      end if
      if (level<nVertLevels) then
         !have cell above
         valp = theta(level+1, iCell)
         zp = calc_heightCellCenter(iCell, level+1, zgrid)
      end if

      if (level==1) then
         !calc_vertDeriv_one(valp, valm, dz)
         gradtheta(3) = calc_vertDeriv_one(valp, val0, zp-z0)
      else if (level==nVertLevels) then
         gradtheta(3) = calc_vertDeriv_one(val0, valm, z0-zm)
      else
         gradtheta(3) = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      end if

      !du/dz and dv/dz
      velCell0(1) = uReconstructX(level,iCell)
      velCell0(2) = uReconstructY(level,iCell)
      velCell0(3) = uReconstructZ(level,iCell)
      if (level>1) then
         !have cell beneath
         velCellm(1) = uReconstructX(level-1,iCell)
         velCellm(2) = uReconstructY(level-1,iCell)
         velCellm(3) = uReconstructZ(level-1,iCell)
      end if
      if (level<nVertLevels) then
         !have cell above
         velCellp(1) = uReconstructX(level+1,iCell)
         velCellp(2) = uReconstructY(level+1,iCell)
         velCellp(3) = uReconstructZ(level+1,iCell)
      end if

      if (level==1) then
         !calc_vertDeriv_one(valp, valm, dz)
         !u
         val0 = dotProduct(velCell0, xyzLocal(:,1),3)
         valp = dotProduct(velCellp, xyzLocal(:,1),3)
         du_dz = calc_vertDeriv_one(valp, val0, zp-z0)
         !v
         val0 = dotProduct(velCell0, xyzLocal(:,2),3)
         valp = dotProduct(velCellp, xyzLocal(:,2),3)
         dv_dz = calc_vertDeriv_one(valp, val0, zp-z0)
      else if (level==nVertLevels) then
         !u
         val0 = dotProduct(velCell0, xyzLocal(:,1),3)
         valm = dotProduct(velCellm, xyzLocal(:,1),3)
         du_dz = calc_vertDeriv_one(val0, valm, z0-zm)
         !v
         val0 = dotProduct(velCell0, xyzLocal(:,2),3)
         valm = dotProduct(velCellp, xyzLocal(:,2),3)
         dv_dz = calc_vertDeriv_one(val0, valm, z0-zm)

      else
         !u
         val0 = dotProduct(velCell0, xyzLocal(:,1),3)
         valp = dotProduct(velCellp, xyzLocal(:,1),3)
         valm = dotProduct(velCellm, xyzLocal(:,1),3)
         du_dz = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
         !v
         val0 = dotProduct(velCell0, xyzLocal(:,2),3)
         valp = dotProduct(velCellp, xyzLocal(:,2),3)
         valm = dotProduct(velCellm, xyzLocal(:,2),3)
         dv_dz = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      end if

      gradxu(3) = calc_verticalVorticity_cell(iCell, level, nEdgesCell0, verticesOnCell, cellsOnVertex, &
                                              kiteAreasOnVertex, areaCell0, vorticity)

      gradxu(1) = dw_dy-dv_dz
      gradxu(2) = du_dz-dw_dx
      
      driver_ertelPV_cell_reconstruct = formErtelPV(gradxu, gradtheta, rho0, xyzLocal(:,1), xyzLocal(:,2), xyzLocal(:,3))
   end function driver_ertelPV_cell_reconstruct

   real(kind=RKIND) function driver_ertelPV_cell_uv(iCell, level, nVertLevels, nEdgesCell0, verticesOnCell, kiteAreasOnVertex, &
                                                 cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                                 cellsOnVertex, &
                                                 cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell0, &
                                                 u, w,theta,rho0,vorticity)

      !We'll use the local xyz system for each cell since it's well aligned with physical gradients.
      !x,y are the tangent plane and z is the radially out.

      !Full curl is (dw/dy-dv/dz)i + (du/dz-dw/dx)j + (dv/dx-du/dy)k
      !k component is solved for already through circulation to vertices. We'll use
      !the kite_areas weighting as in pv_cell to get it at the cell centers.

      !The other components of vorticity are:
      !(dw/dy-dv/dz)i + (du/dz-dw/dx)j where velocities are in the local coordinate system.
      !We'll also need dTheta/dX, dTheta/dY, dTheta/dZ.
      !For all these, it's finite volume gradients in the
      !horizontal and finite differencing in the vertical.
                                                 
      implicit none

      real(kind=RKIND), intent(in) :: rho0, areaCell0
      real(kind=RKIND), dimension(:), intent(in) :: dvEdge
      real(kind=RKIND), dimension(3,2,*), intent(in) :: cellTangentPlane
      real(kind=RKIND), dimension(3,*), intent(in) :: localVerticalUnitVectors, edgeNormalVectors
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid,u,w, theta, vorticity, kiteAreasOnVertex
      integer, intent(in) :: iCell, level, nVertLevels, nEdgesCell0
      integer, dimension(:,:), intent(in) :: cellsOnCell, edgesOnCell, cellsOnEdge, verticesOnCell, cellsOnVertex
      
      integer :: i, iNbr, iEdge
      real(kind=RKIND) :: val0, valNbr, volumeCell, areaFactor, z0, zp, zm, valp, valm, dw_dx, dw_dy, du_dz, dv_dz
      real(kind=RKIND), dimension(3) :: gradxu, gradtheta, unitDeriv
      real(kind=RKIND), dimension(3,3) :: xyzLocal
      real(kind=RKIND), dimension(nEdgesCell0) :: valEdges, dvEdgeCell, dhEdge
      real(kind=RKIND), dimension(3,nEdgesCell0) :: normalEdgeCell

      write(0,*) ' Dont use unless make full velocity at cell edges through u and v '

      !local coordinate system
      call coordinateSystem_cell(cellTangentPlane, localVerticalUnitVectors, iCell, xyzLocal)

      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         dhEdge(i) = calc_heightVerticalEdge(iCell, iNbr, level, zgrid)

         iEdge = edgesOnCell(i,iCell)
         dvEdgeCell(i) = dvEdge(iEdge)
         val0 = fluxSign(iCell, iEdge, cellsOnEdge)
         normalEdgeCell(:,i) = edgeNormalVectors(:,iEdge)
         call normalizeVector(normalEdgeCell(:,i),3)
         normalEdgeCell(:,i) = normalEdgeCell(:,i)*val0
      end do

      volumeCell = calcVolumeCell(areaCell0, nEdgesCell0, dhEdge)

      !Need to get 3d curl and grad theta
      !horizontal derivatives
      !calc_horizDeriv_fv(valEdges, nNbrs, dvEdge, dhEdge, &
      !                                         normalEdge, unitDeriv, volumeCell)
      !theta
      val0 = theta(level, iCell)
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         valNbr = theta(level,iNbr)
         valEdges(i) = 0.5*(valNbr+val0)
      end do
      unitDeriv(:) = xyzLocal(:,1)
      gradtheta(1) = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)
      unitDeriv(:) = xyzLocal(:,2)
      gradtheta(2) = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)

      !w
      val0 = .5*(w(level+1, iCell)+w(level, iCell))
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         valNbr = .5*(w(level+1, iNbr)+w(level, iNbr))
         valEdges(i) = 0.5*(valNbr+val0)
      end do
      unitDeriv(:) = xyzLocal(:,1)
      dw_dx = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)
      unitDeriv(:) = xyzLocal(:,2)
      dw_dy = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)

      !vertical derivatives
      !calc_heightCellCenter(c0, level, zgrid) calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      !theta
      gradtheta(3) = 0.0_RKIND
      z0 = calc_heightCellCenter(iCell, level, zgrid)
      val0 = theta(level, iCell)
      if (level>1) then
         !have cell beneath
         valm = theta(level-1, iCell)
         zm = calc_heightCellCenter(iCell, level-1, zgrid)
      end if
      if (level<nVertLevels) then
         !have cell above
         valp = theta(level+1, iCell)
         zp = calc_heightCellCenter(iCell, level+1, zgrid)
      end if

      if (level==1) then
         !calc_vertDeriv_one(valp, valm, dz)
         gradtheta(3) = calc_vertDeriv_one(valp, val0, zp-z0)
      else if (level==nVertLevels) then
         gradtheta(3) = calc_vertDeriv_one(val0, valm, z0-zm)
      else
         gradtheta(3) = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      end if

      !du/dz and dv/dz
      !get dux/dz at faces and then weighted average by area in x-direction to cell center
      du_dz = calc_vertDeriv_velocity(nEdgesCell0, xyzLocal(:,1), dhEdge, dvEdgeCell, cellsOnCell, &
                                      edgesOnCell, edgeNormalVectors, u, zgrid, iCell, level, nVertLevels)
      dv_dz = calc_vertDeriv_velocity(nEdgesCell0, xyzLocal(:,2), dhEdge, dvEdgeCell, cellsOnCell, &
                                      edgesOnCell, edgeNormalVectors, u, zgrid, iCell, level, nVertLevels)

      !(dw/dy-dv/dz)i + (du/dz-dw/dx)j + (dv/dx-du/dy)k
      gradxu(3) = calc_verticalVorticity_cell(iCell, level, nEdgesCell0, verticesOnCell, cellsOnVertex, &
                                              kiteAreasOnVertex, areaCell0, vorticity)

      gradxu(1) = dw_dy-dv_dz
      gradxu(2) = du_dz-dw_dx
      
      driver_ertelPV_cell_uv = formErtelPV(gradxu, gradtheta, rho0, xyzLocal(:,1), xyzLocal(:,2), xyzLocal(:,3))
   end function driver_ertelPV_cell_uv
   
   real(kind=RKIND) function dotProduct(a, b, sz)

      implicit none

      real(kind=RKIND), dimension(:), intent(in) :: a, b
      integer, intent(in) :: sz

      integer :: i
      real(kind=RKIND) :: rsum

      rsum = 0.0_RKIND

      do i=1,sz
         rsum = rsum + a(i)*b(i)
      end do

      dotProduct = rsum
   end function dotProduct

   integer function elementIndexInArray(val, array, sz)

      implicit none

      integer, intent(in) :: val
      integer, dimension(:), intent(in) :: array
      integer, intent(in) :: sz

      integer :: i, ind
      ind = -1
      do i=1,sz
         if (array(i)==val) then
            ind = i
            elementIndexInArray = ind !This returns, right?
         end if
      end do
      elementIndexInArray = ind
   end function elementIndexInArray

   real(kind=RKIND) function formErtelPV(gradxu, gradtheta, density, unitX, unitY, unitZ)

      use mpas_constants, only : omega_e => omega

      implicit none

      real(kind=RKIND), dimension(3), intent(inout) :: gradxu
      real(kind=RKIND), dimension(3), intent(in) :: gradtheta
      real(kind=RKIND), intent(in) :: density
      real(kind=RKIND), dimension(3), intent(in) :: unitX, unitY, unitZ

      real(kind=RKIND) :: epv, eVort
      real(kind=RKIND), dimension(3) :: eVortDir, eVortComponents

      !earth vorticity is in z-direction in global Cartesian space
      eVort = 2.0 * omega_e
      eVortDir(1) = 0.0_RKIND
      eVortDir(2) = 0.0_RKIND
      eVortDir(3) = eVort

      eVortComponents(1) = dotProduct(eVortDir, unitX,3)
      eVortComponents(2) = dotProduct(eVortDir, unitY,3)
      eVortComponents(3) = dotProduct(eVortDir, unitZ,3)

      gradxu(:) = gradxu(:) + eVortComponents(:)

      epv = dotProduct(gradxu, gradtheta,3) / density

      epv = epv * 1.0e6 !SI to PVUs
    
      formErtelPV = epv
   end function formErtelPV

   real(kind=RKIND) function calc_verticalVorticity_cell(c0, level, nVerticesOnCell, verticesOnCell, cellsOnVertex, &
                                                         kiteAreasOnVertex, areaCell, vVortVertex)
      !area weighted average of vorticity at vertices (really midpts of faces) to cell center for the specified cell
      !
      implicit none

      real(kind=RKIND), intent(in) :: areaCell
      integer, intent(in) :: c0, level, nVerticesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell, cellsOnVertex
      real(kind=RKIND), dimension(:,:), intent(in) :: kiteAreasOnVertex, vVortVertex

      real(kind=RKIND) :: vVortCell
      integer :: i, iVertex, cellIndOnVertex

      vVortCell = 0.0_RKIND
      do i = 1,nVerticesOnCell
         iVertex = verticesOnCell(i,c0)
         cellIndOnVertex = elementIndexInArray(c0, cellsOnVertex(:,iVertex), 3)
         vVortCell = vVortCell + kiteAreasOnVertex(cellIndOnVertex, iVertex)*vVortVertex(level, iVertex)/areaCell
      end do

      calc_verticalVorticity_cell = vVortCell
   end function calc_verticalVorticity_cell

   subroutine coordinateSystem_cell(cellTangentPlane, localVerticalUnitVectors, c0, xyz)

      implicit none

      real(kind=RKIND), dimension(3,2,*), intent(in) :: cellTangentPlane
      real(kind=RKIND), dimension(3,*), intent(in) :: localVerticalUnitVectors
      integer, intent(in) :: c0
      real(kind=RKIND), dimension(3,3), intent(out) :: xyz

      integer :: i

      xyz(:,1) = cellTangentPlane(:,1,c0) !are these guaranteed unit vectors?
      xyz(:,2) = cellTangentPlane(:,2,c0)
      xyz(:,3) = localVerticalUnitVectors(:,c0)
      do i=1,2
         call normalizeVector(xyz(:,i), 3)
      end do
   end subroutine coordinateSystem_cell

   real(kind=RKIND) function fluxSign(c0, iEdge, cellsOnEdge)
      
      !For finite volume computations, we'll use a normal pointing out of the cell
      implicit none

      integer, intent(in) :: c0
      integer, intent(in) :: iEdge
      integer, dimension(:,:), intent(in) :: cellsOnEdge

      if (c0 == cellsOnEdge(1,iEdge)) then
         fluxSign = 1.0_RKIND
      else
         fluxSign = -1.0_RKIND
      end if
   end function fluxSign

   real(kind=RKIND) function calc_heightCellCenter(c0, level, zgrid)

      implicit none

      integer, intent(in) :: c0, level
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid

      calc_heightCellCenter = 0.5*(zgrid(level,c0)+zgrid(level+1,c0))
   end function calc_heightCellCenter

   real(kind=RKIND) function calc_heightVerticalEdge(c0, c1, level, zgrid)

      implicit none

      integer, intent(in) :: c0, c1, level
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid

      real(kind=RKIND) :: hTop, hBottom

      hTop = .5*(zgrid(level+1,c0)+zgrid(level+1,c1))
      hBottom = .5*(zgrid(level,c0)+zgrid(level,c1))

      calc_heightVerticalEdge = hTop-hBottom
   end function calc_heightVerticalEdge

   subroutine normalizeVector(vals, sz)
      !normalize a vector to unit magnitude
      implicit none

      real (kind=RKIND), dimension(:), intent(inout) :: vals
      integer, intent(in) :: sz

      integer :: i
      real (kind=RKIND) :: mag

      mag = 0.0_RKIND !sqrt(sum(squares))
      do i=1,sz
         mag = mag+vals(i)*vals(i)
      end do
      mag = sqrt(mag)
      vals(:) = vals(:)/mag
   end subroutine normalizeVector

   real(kind=RKIND) function calcVolumeCell(areaCell, nEdges, hEdge)

      implicit none

      integer, intent(in) :: nEdges
      real(kind=RKIND), intent(in) :: areaCell
      real(kind=RKIND), dimension(nEdges), intent(in) :: hEdge

      integer :: i
      real(kind=RKIND) :: avgHt, vol

      avgHt = 0.0_RKIND
      do i=1,nEdges
         avgHt = avgHt + hEdge(i)
      end do
      avgHt = avgHt/nEdges

      vol = areaCell*avgHt
      calcVolumeCell = vol
   end function calcVolumeCell

   real(kind=RKIND) function calc_horizDeriv_fv(valEdges, nNbrs, dvEdge, dhEdge, &
                                                normalEdge, unitDeriv, volumeCell)
      !normals to edges point out of cell
      implicit none

      integer, intent(in) :: nNbrs
      real(kind=RKIND), dimension(:), intent(in) :: valEdges, dvEdge, dhEdge
      real(kind=RKIND), dimension(3,nNbrs), intent(in) :: normalEdge
      real(kind=RKIND), dimension(3), intent(in) :: unitDeriv
      real(kind=RKIND), intent(in) :: volumeCell

      integer :: i
      real(kind=RKIND) :: vale, rsum, areaFace
      real(kind=RKIND), dimension(3) :: unitNormalEdge

      rsum = 0.0_RKIND
      do i=1,nNbrs
         vale = valEdges(i) !0.5 * (val0 + valNbrs(i))
         areaFace = dvEdge(i) * dhEdge(i)
         unitNormalEdge(:) = normalEdge(:,i)
         call normalizeVector(unitNormalEdge,3)
         areaFace = areaFace*dotProduct(unitNormalEdge, unitDeriv,3)  !* abs(dotProduct(unitNormalEdge, unitDeriv,3))
         rsum = rsum + vale * areaFace
      end do
      rsum = rsum / volumeCell

      calc_horizDeriv_fv = rsum
   end function calc_horizDeriv_fv

   real(kind=RKIND) function calc_vertDeriv_velocity(nEdgesCell0, unitDeriv, hEdge, dvEdge, cellsOnCell, &
                                                     edgesOnCell, &
                                                     edgeNormalVectors, u, zgrid, iCell, level, nVertLevels)
      implicit none

      integer, intent(in) :: nEdgesCell0, iCell, level, nVertLevels
      integer, dimension(:,:), intent(in) :: cellsOnCell, edgesOnCell
      real(kind=RKIND), dimension(3), intent(in) :: unitDeriv
      real(kind=RKIND), dimension(nEdgesCell0), intent(in) :: hEdge, dvEdge
      real(kind=RKIND), dimension(3,*), intent(in) :: edgeNormalVectors 
      real(kind=RKIND), dimension(:,:), intent(in) :: u, zgrid

      integer :: i, iNbr, iEdge
      real(kind=RKIND) :: areaFace, val0, valp, valm, z0, zp, zm, dval_dz, netArea
      real(kind=RKIND), dimension(3) :: edgeUnitNormal
      real(kind=RKIND), dimension(nEdgesCell0) :: areaFactor, valEdges
      
      if (level==1) then
         !just cell and above.
         !get derivatives at faces
         do i = 1, nEdgesCell0
            iNbr = cellsOnCell(i, iCell)
            z0 = .5*(calc_heightCellCenter(iCell, level, zgrid)+calc_heightCellCenter(iNbr, level, zgrid))
            zp = .5*(calc_heightCellCenter(iCell, level+1, zgrid)+calc_heightCellCenter(iNbr, level+1, zgrid))
            iEdge = edgesOnCell(i,iCell)
            val0 = u(level, iEdge)
            valp = u(level+1, iEdge)

            !velocity component in direction of derivative
            edgeUnitNormal(:) = edgeNormalVectors(:,iEdge)
            call normalizeVector(edgeUnitNormal,3)
            areaFactor(i) = dotProduct(edgeUnitNormal, unitDeriv,3)
            val0 = val0*areaFactor(i)
            valp = valp*areaFactor(i)
            valEdges(i) = calc_vertDeriv_one(valp, valm, zp-z0)
         end do

      else if (level==nVertLevels) then
         !cell and below
         do i = 1, nEdgesCell0
            iNbr = cellsOnCell(i, iCell)
            z0 = .5*(calc_heightCellCenter(iCell, level, zgrid)+calc_heightCellCenter(iNbr, level, zgrid))
            zm = .5*(calc_heightCellCenter(iCell, level-1, zgrid)+calc_heightCellCenter(iNbr, level-1, zgrid))
            iEdge = edgesOnCell(i,iCell)
            val0 = u(level, iEdge)
            valm = u(level-1, iEdge)

            !velocity component in direction of derivative
            edgeUnitNormal(:) = edgeNormalVectors(:,iEdge)
            call normalizeVector(edgeUnitNormal,3)
            areaFactor(i) = dotProduct(edgeUnitNormal, unitDeriv,3)
            val0 = val0*areaFactor(i)
            valm = valm*areaFactor(i)
            valEdges(i) = calc_vertDeriv_one(val0, valm, z0-zm)
         end do

      else
         !cells above and below
         do i = 1, nEdgesCell0
            iNbr = cellsOnCell(i, iCell)
            z0 = .5*(calc_heightCellCenter(iCell, level, zgrid)+calc_heightCellCenter(iNbr, level, zgrid))
            zp = .5*(calc_heightCellCenter(iCell, level+1, zgrid)+calc_heightCellCenter(iNbr, level+1, zgrid))
            zm = .5*(calc_heightCellCenter(iCell, level-1, zgrid)+calc_heightCellCenter(iNbr, level-1, zgrid))
            iEdge = edgesOnCell(i,iCell)
            val0 = u(level, iEdge)
            valp = u(level+1, iEdge)
            valm = u(level-1, iEdge)

            !velocity component in direction of derivative
            edgeUnitNormal(:) = edgeNormalVectors(:,iEdge)
            call normalizeVector(edgeUnitNormal,3)
            areaFactor(i) = dotProduct(edgeUnitNormal, unitDeriv,3)
            val0 = val0*areaFactor(i)
            valp = valp*areaFactor(i)
            valm = valm*areaFactor(i)
            valEdges(i) = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
         end do
      end if

      !weight derivatives to cell centers
      netArea = 0.0_RKIND
      dval_dz = 0.0_RKIND
      do i = 1, nEdgesCell0
         areaFace = hEdge(i)*dvEdge(i)*abs(areaFactor(i))
         dval_dz = dval_dz + valEdges(i)*areaFace
         netArea = netArea + areaFace
      end do
      dval_dz = dval_dz / netArea

      calc_vertDeriv_velocity = dval_dz
   end function calc_vertDeriv_velocity

   !cell centers are halfway between w faces
   real(kind=RKIND) function calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)

      implicit none

      real(kind=RKIND), intent(in) :: val0, valp, valm, z0,zp,zm !center, plus, minus
      
      real(kind=RKIND) :: dval_dzp, dval_dzm

      !Average 1 sided differences to below and above since not equally spaced pts
      dval_dzp = calc_vertDeriv_one(valp, val0, zp-z0)
      dval_dzm = calc_vertDeriv_one(val0, valm, z0-zm)
      calc_vertDeriv_center = 0.5*(dval_dzp+dval_dzm)

   end function calc_vertDeriv_center


   real(kind=RKIND) function calc_vertDeriv_one(valp, valm, dz)
      !1 sided finite difference

      implicit none

      real(kind=RKIND), intent(in) :: valp, valm, dz

      calc_vertDeriv_one = (valp - valm) / dz

   end function calc_vertDeriv_one
   
   subroutine interp_pv( nCells, nLevels, interpVal, &
                         latCell, field0, field1,field, &
                         missingVal)

      implicit none
      ! interpolate columns of field1 to where field0 is interpVal*sign(lat) such that
      ! vals in cells above index are above interpVal and cells below index are below interpVal.
      ! to limit higher level spikes (eg from mixing ~17km), require change in sign to persist >1 levels.
      !if whole column is above value, use surface value.
      
      ! input

      integer :: nCells, nLevels
      real(kind=RKIND) ::  interpVal, missingVal
      real(kind=RKIND), intent(in) ::latCell(nCells)
      real(kind=RKIND), intent(in) :: field0(nLevels,nCells), field1(nLevels,nCells)
      real(kind=RKIND), intent(out) :: field(nCells)

      !  local
      
      integer :: iCell, iLev, levInd, indlNbr
      real(kind=RKIND) :: valh, vall, vallNbr, sgnh, sgnl, sgnlNbr
      real(kind=RKIND) :: dv_dl, levFrac, valInterpCell, sgnHemi

      do iCell = 1, nCells
        !starting from top, trap val if values on opposite side
        levInd = -1 !what should happen with missing values?
        levFrac = 0.0
        sgnHemi = sign(1.0_RKIND, latCell(iCell)) !problem at the equator...is sign(0)=0?
        if (sgnHemi .EQ. 0.0) sgnHemi = 1.0
        valInterpCell = interpVal*sgnHemi
        do iLev = nLevels,2,-1
        !do iLev = 2,nLevels-1,1
          valh = field0(iLev,iCell)
          vall = field0(iLev-1,iCell)
          indlNbr = max(iLev-2,1)
          vallNbr = field0(indlNbr, iCell)
          sgnh = (valh-valInterpCell)*sgnHemi
          sgnl = (vall-valInterpCell)*sgnHemi
          sgnlNbr = (vallNbr-valInterpCell)*sgnHemi
          if (sgnh>= 0.0 .AND. sgnl<=0.0 .AND. sgnlNbr<=0.0) then
            !sandwiched value. equal in case val0 is a vals[l].
            !get linear interpolation: val0 = vals[l]+dvals/dl * dl
            !Avoid divide by 0 by just assuming value is 
            !halfway between...
     
            dv_dl = valh-vall;
            if (abs(dv_dl)<1.e-6) then
              levFrac = 0.5;
            else
              levFrac = (valInterpCell-vall)/dv_dl
            end if
            
            levInd = iLev-1
            
            exit !once found, stop searching column
          end if
        end do ! done searching column

        !find value of field using index we just found
        if (levInd<0) then
          if (sgnl>0.0) then
            field(iCell) = field1(1,iCell)
          else
            field(iCell) = missingVal
          end if
        else
          valh = field1(levInd+1,iCell)
          vall = field1(levInd,iCell)
        
          dv_dl = valh-vall
          field(iCell) = vall+dv_dl*levFrac
        end if
      end do
      
   end subroutine interp_pv
   
   real (kind=RKIND) function sphere_distance(lat1, lon1, lat2, lon2, radius)

     !Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
     !sphere with given radius.
     !==================================================================================================
     implicit none

     real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius
     real (kind=RKIND) :: arg1

     arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                  cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
     sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance


  subroutine atm_compute_vortexLoc(state, diag, grid, domain)
     
      use mpas_grid_types
      use mpas_constants
   
      implicit none
      
      !Track the min of the vortex on 2pvu from previous time and
      !update inVortex[iCell]>0 if cell w/in radius
   
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      type (mesh_type), intent(in) :: grid
      type (domain_type), intent(in) :: domain
   
      integer :: iCell, k
      real(kind=RKIND) :: pvuVal, missingVal, valCell, theta_min, global_theta_min
      real(kind=RKIND) :: vortexLat, vortexLon, lat0, lon0, bigVal, dist, maxDist
      maxDist = 2.0e5 !200km

      do iCell=1,grid%nCellsSolve
         do k=1,grid%nVertLevels
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1._RKIND + rvord * state % scalars % array(state % index_qv,k,iCell))
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * grid % zz % array(k,iCell)
         end do
      end do
      
      write(0,*) "Calculating ertel pv"
      !need halo cells for everything w/ horizontal derivative
      call mpas_dmpar_exch_halo_field(diag % theta)
      call mpas_dmpar_exch_halo_field(diag % uReconstructX)
      call mpas_dmpar_exch_halo_field(diag % uReconstructY)
      call mpas_dmpar_exch_halo_field(diag % uReconstructZ)
      call mpas_dmpar_exch_halo_field(state % w)
      
      call calc_ertel_pv(grid, state, diag)
      
      pvuVal = 2.0_RKIND
      missingVal = -99999.0_RKIND
      call interp_pv(grid%nCells, grid%nVertLevels, pvuVal, grid % latCell % array, &
                     diag % ertel_pv % array, diag % theta % array, diag % theta_pv % array, &
                     missingVal)
      
      !now we track min within a spatial disk
      !since i don't know how to make this work for >1 block/process,
      !only use this for 1 block for process :)
      vortexLat = diag % vortexLat % scalar
      vortexLon = diag % vortexLon % scalar
      
      bigVal = 99999.0_RKIND
      theta_min = bigVal
      do iCell = 1, grid % nCellsSolve
        !include cell in min only if
        !(1) cell within search disk by distance and (2) value not missing value
        if (diag % theta_pv % array(iCell)==missingVal) cycle
        
        dist = sphere_distance(grid%latCell%array(iCell), grid%lonCell%array(iCell), vortexLat, vortexLon, a)
        if (dist>maxDist) cycle
        
        valCell = diag % theta_pv % array(iCell)
        if (valCell<theta_min) then
          theta_min = valCell
          lat0 = grid%latCell%array(iCell)
          lon0 = grid%lonCell%array(iCell)
        end if
      end do

      call mpas_dmpar_min_real(domain%dminfo, theta_min, global_theta_min)
      if (theta_min>global_theta_min) then
       !then block doesn't have the minimum cell.
       !rig this so the block w/ the minimum will update all others
       lat0 = bigVal
       lon0 = bigVal
       !n.b. as is, this can mix lat/lon of >1 cells with equal min values
       !but what are the odds of 2 floats being bit equal?
       !They're close in space as well
      end if

      call mpas_dmpar_min_real(domain%dminfo, lat0, vortexLat)
      call mpas_dmpar_min_real(domain%dminfo, lon0, vortexLon)
      
      diag % vortexLat % scalar = vortexLat
      diag % vortexLon % scalar = vortexLon
      
      write(0,*) 'vortex lat/lon: ', vortexLat, vortexLon
      
      !update inVortex array
      maxDist = 8.0e5 !800km
      diag % inVortex % array(:) = 0
      do iCell = 1, grid % nCellsSolve
        dist = sphere_distance(grid%latCell%array(iCell), grid%lonCell%array(iCell), vortexLat, vortexLon, a)
        if (dist .le. maxDist) then
          diag % inVortex % array(iCell) = 1
        end if
      end do
        
  end subroutine atm_compute_vortexLoc

!***********************************************************************
!
!  routine mpas_core_setup_packages
!
!> \brief   Pacakge setup routine
!> \author  Doug Jacobsen
!> \date    September 2011
!> \details 
!>  This routine is intended to correctly configure the packages for this MPAS
!>   core. It can use any Fortran logic to properly configure packages, and it
!>   can also make use of any namelist options. All variables in the model are
!>   *not* allocated until after this routine is called.
!
!-----------------------------------------------------------------------
   subroutine mpas_core_setup_packages(ierr)!{{{

      use mpas_packages

      implicit none

      integer, intent(out) :: ierr

      ierr = 0

   end subroutine mpas_core_setup_packages!}}}

end module mpas_core
