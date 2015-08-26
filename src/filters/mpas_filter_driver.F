! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!***********************************************************************
!
!  mpas_filter_driver
!
!> \brief   MPAS Filters Module
!> \author  Doug Jacobsen
!> \date    07/12/2015
!> \details
!> This module defines the routines that perform actions using a filter list
!> created from the mpas_filter_list module.
!
!-----------------------------------------------------------------------
module mpas_filter_driver
   use mpas_kind_types
   use mpas_derived_types
   use mpas_filter_list
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_io_units
   use mpas_stream_manager
   use mpas_timekeeping
#ifdef MPAS_TIME_FILTER_DRIVER
   use mpas_timer
#endif

   implicit none
   private

   public :: mpas_filter_driver_setup_packages, mpas_filter_driver_init
   public :: mpas_filter_driver_compute, mpas_filter_driver_restart, mpas_filter_driver_finalize

   contains


!***********************************************************************
!
!  routine mpas_filter_driver_setup_packages
!
!> \brief   MPAS Filters Setup Packages routine
!> \author  Doug Jacobsen
!> \date    07/12/2015
!> \details
!> This routine iterates over a list of filers, and sets up packages for each
!> depending on config option values. The optional argument 'filterList' can be
!> used to setup the packages of a subset of filters defined in the inList
!> mpas_filter_list_type.
!> Each filter defined is not required to have a package. If the packge is not
!> found in the packagePool, it is skipped.
!> If a package is found, and the _enable namelist option is true, the package
!> is activated.
!
!-----------------------------------------------------------------------
   recursive subroutine mpas_filter_driver_setup_packages(inList, packagePool, configPool, filterList, iErr)!{{{
      type (mpas_filter_list_type), pointer :: inList
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), pointer, optional :: filterList
      integer, intent(out) :: iErr

      type (mpas_pool_type), pointer :: filterListLocal

      type (mpas_pool_type), pointer :: subPool
      integer :: nameLength
      character (len=StrKIND) :: configName, packageName, filterName
      type (mpas_pool_iterator_type) :: poolItr
      logical, pointer :: filterPackage
      logical, pointer :: config_filterEnable

      iErr = MPAS_FILTER_NOERR

      if ( present(filterList) ) then
         filterListLocal => filterList
      else
         call mpas_pool_create_pool(filterListLocal)
         call mpas_filter_list_get_filter_names(inList, filterListLocal)
      end if

      call mpas_pool_begin_iteration(filterListLocal)
      do while ( mpas_pool_get_next_member(filterListLocal, poolItr) )
         if ( poolItr % memberType == MPAS_POOL_CONFIG ) then
            nameLength = len_trim(poolItr % memberName)
            configName = 'config_' // poolItr % memberName(1:nameLength) // '_enable'
            call mpas_pool_get_config(configPool, configName, config_filterEnable)
            if ( associated(config_filterEnable) ) then
               if ( config_filterEnable ) then
                  packageName = poolItr % memberName(1:nameLength) // 'PKG'
                  call mpas_pool_get_package(packagePool, packageName, filterPackage)
                  if ( associated(filterPackage) ) then
                     filterPackage = .true.
                  end if
               end if
            end if
         else if ( poolItr % memberType == MPAS_POOL_SUBPOOL ) then
            call mpas_pool_get_subpool(filterListLocal, poolItr % memberName, subPool)
            call mpas_filter_driver_setup_packages(inList, packagePool, configPool, subPool, iErr)
         end if
      end do

      if ( .not. present(filterList) ) then
         call mpas_pool_destroy_pool(filterListLocal)
      end if

   end subroutine mpas_filter_driver_setup_packages!}}}


!***********************************************************************
!
!  routine mpas_filter_driver_init
!
!> \brief   MPAS Filters init routine
!> \author  Doug Jacobsen
!> \date    07/12/2015
!> \details
!> This routine iterates over a list of filers, and initializes each.
!> The optional argument 'filterList' can be used to act on a subset of filters
!> defined in the inList mpas_filter_list_type.
!>
!> This routine also verifies that the stream specified in the filter's
!> stream_name exists in the domain's streamManager.
!>
!> This routine also sets up the compute interval for the filter.
!>
!> This routine will also compute the filter, if it's compute_on_startup
!> namelist option is true, and write the filter if it's write_on_startup
!> namelist option is true.
!
!-----------------------------------------------------------------------
   recursive subroutine mpas_filter_driver_init(inList, domain, filterList, iErr)!{{{
      type (mpas_filter_list_type), pointer :: inList
      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer, optional :: filterList
      integer, intent(out) :: iErr

      type (mpas_pool_type), pointer :: filterListLocal
      type (mpas_pool_type), pointer :: subPool

      character (len=StrKIND) :: configName, alarmName, streamName, timerName, filterName
      character (len=StrKIND), pointer :: filterComputeInterval, filterStreamName
      logical, pointer :: filterEnable, filterComputeStartup, filterWriteStartup
      type (mpas_pool_iterator_type) :: poolItr
      integer :: nameLength

      logical :: streamFound
      character  (len=StrKIND) :: referenceTimeString, outputIntervalString
      type (MPAS_Time_Type) :: referenceTime 
      type (MPAS_TimeInterval_type) :: alarmTimeStep

      procedure (mpas_filter_list_init_function), pointer :: initFunc
      procedure (mpas_filter_list_compute_function), pointer :: computeFunc

      integer :: err_local

      iErr = MPAS_FILTER_NOERR

      if ( present(filterList) ) then
         filterListLocal => filterList
      else
         call mpas_pool_create_pool(filterListLocal)
         call mpas_filter_list_get_filter_names(inList, filterListLocal)
      end if

      call mpas_pool_begin_iteration(filterListLocal)
      do while ( mpas_pool_get_next_member(filterListLocal, poolItr) )
         if ( poolItr % memberType == MPAS_POOL_CONFIG ) then
            nameLength = len_trim(poolItr % memberName)
   
            configName = 'config_' // poolItr % memberName(1:nameLength) // '_enable'
            call mpas_pool_get_config(domain % configs, configName, filterEnable)
   
            if ( .not. associated(filterEnable) ) then
               call mpas_dmpar_global_abort('ERROR: Namelist option ''' // trim(configName) // ''' not found for filter ''' &
                                            // poolItr % memberName(1:nameLength) // '''')
            end if
   
            if ( filterEnable ) then
#ifdef MPAS_TIME_FILTER_DRIVER
               timerName = trim(initTimerPrefix) // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName, .false.)
#endif

               ! Initialize filter
               call mpas_filter_list_get_filter_init(inList, poolItr % memberName(1:nameLength), initFunc, err_local)
               err_local = initFunc(domain, poolItr % memberName(1:nameLength))
   
               ! Setup compute interval
               configName = 'config_' // poolItr % memberName(1:nameLength) // '_compute_interval'
               call mpas_pool_get_config(domain % configs, configName, filterComputeInterval)
   
               if ( .not. associated(filterComputeInterval) ) then
                  call mpas_dmpar_global_abort('ERROR: Namelist option ''' // trim(configName) // ''' not found for filter ''' &
                                               // poolItr % memberName(1:nameLength) // '''')
               end if
   
               configName = 'config_' // poolItr % memberName(1:nameLength) // '_stream_name'
               call mpas_pool_get_config(domain % configs, configName, filterStreamName)
   
               if ( .not. associated(filterStreamName) ) then
                  call mpas_dmpar_global_abort('ERROR: Namelist option ''' // trim(configName) // ''' not found for filter ''' &
                                               // poolItr % memberName(1:nameLength) // '''')
               end if
   
               ! Ensure stream exists in stream manager
               if ( filterStreamName /= 'none' ) then
                  streamFound = .false.
   
                  call mpas_stream_mgr_begin_iteration(domain % streamManager)
                  do while (mpas_stream_mgr_get_next_stream(domain % streamManager, streamName) )
                     if ( trim(streamName) == trim(filterStreamName) ) then
                        streamFound = .true.
                     end if
                  end do
   
                  if ( .not. streamFound ) then
                     call mpas_dmpar_global_abort('ERROR: Stream ''' // trim(filterStreamName) // ''' not found for filter ''' &
                                                  // poolItr % memberName(1:nameLength) // '''')
                  end if
               end if
   
               ! Define alarm for compute interval
               if ( filterComputeInterval /= 'output_interval' .and. filterComputeInterval /= 'dt' &
                                                               .and. filterStreamName /= 'none' ) then
                  alarmName = poolItr % memberName(1:nameLength) // computeAlarmSuffix
                  call mpas_set_timeInterval(alarmTimeStep, timeString=filterComputeInterval, ierr=err_local)
                  call MPAS_stream_mgr_get_property(domain % streamManager, filterStreamName, MPAS_STREAM_PROPERTY_REF_TIME, &
                                                    referenceTimeString, err_local)
                  call mpas_set_time(referenceTime, dateTimeString=referenceTimeString, ierr=err_local)
                  call mpas_add_clock_alarm(domain % clock, alarmName, referenceTime, alarmTimeStep, ierr=err_local)
                  call mpas_reset_clock_alarm(domain % clock, alarmName, ierr=err_local)
               end if
   
               ! Compute and write if requested on startup.
               configName = 'config_' // poolItr % memberName(1:nameLength) // '_compute_on_startup'
               call mpas_pool_get_config(domain % configs, configName, filterComputeStartup)
   
               if ( .not. associated(filterComputeStartup) ) then
                  call mpas_dmpar_global_abort('ERROR: Namelist option ''' // trim(configName) // ''' not found for filter ''' &
                                               // poolItr % memberName(1:nameLength) // '''')
               end if
   
               configName = 'config_' // poolItr % memberName(1:nameLength) // '_write_on_startup'
               call mpas_pool_get_config(domain % configs, configName, filterWriteStartup)
   
               if ( .not. associated(filterWriteStartup) ) then
                  call mpas_dmpar_global_abort('ERROR: Namelist option ''' // trim(configName) // ''' not found for filter ''' &
                                               // poolItr % memberName(1:nameLength) // '''')
               end if
   
               if ( filterComputeStartup ) then
                  ! Call compute function with a time level of 1. On init, a time level of 2 or greater doesn't make a lot of sense.
                  call mpas_filter_list_get_filter_compute(inList, poolItr % memberName(1:nameLength), computeFunc, iErr=err_local)
                  err_local = computeFunc(domain, poolItr % memberName(1:nameLength), 1)
   
                  if ( filterWriteStartup ) then
                     if ( filterStreamName /= 'none' ) then
                        call mpas_stream_mgr_write(domain % streamManager, streamID=filterStreamName, forceWriteNow=.true., &
                                                   iErr=err_local)
                     end if
                  end if
               else
                  if ( filterWriteStartup ) then
                     write(stderrUnit, *) 'WARNING: write_on_startup called for filter ''' // poolItr % memberName(1:nameLength) &
                                          // ''' without compute_on_startup. Skipping write...'
                  end if
               end if
   
   
               ! Ensure stream exists in stream manager
               if ( filterStreamName /= 'none' ) then
                  call mpas_stream_mgr_reset_alarms(domain % streamManager, streamName, direction=MPAS_STREAM_OUTPUT, iErr=err_local)
               end if

#ifdef MPAS_TIME_FILTER_DRIVER
               call mpas_timer_stop(timerName)
#endif
            end if
         else if ( poolItr % memberType == MPAS_POOL_SUBPOOL ) then
            call mpas_pool_get_subpool(filterListLocal, poolItr % memberName, subPool)
            call mpas_filter_driver_init(inList, domain, subPool, iErr=err_local)
         end if
      end do

      if ( .not. present(filterList) ) then
         call mpas_pool_destroy_pool(filterListLocal)
      end if

   end subroutine mpas_filter_driver_init!}}}


!***********************************************************************
!
!  routine mpas_filter_driver_compute
!
!> \brief   MPAS Filters compute routine
!> \author  Doug Jacobsen
!> \date    07/12/2015
!> \details
!> This routine iterates over a list of filers, and computes each based on it's
!> compute_interval.
!> The optional argument 'filterList' can be used to act on a subset of filters
!> defined in the inList mpas_filter_list_type.
!
!-----------------------------------------------------------------------
   recursive subroutine mpas_filter_driver_compute(inList, domain, timeLevel, filterList, iErr)!{{{
      type (mpas_filter_list_type), pointer :: inList
      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer, optional :: filterList
      integer, intent(in), optional :: timeLevel
      integer, intent(out) :: iErr

      type (mpas_pool_type), pointer :: filterListLocal
      integer :: timeLevelLocal

      type (mpas_pool_type), pointer :: subPool

      character (len=StrKIND) :: configName, alarmName, timerName
      character (len=StrKIND), pointer :: filterStreamName, filterComputeInterval
      logical, pointer :: filterEnable
      type (mpas_pool_iterator_type) :: poolItr
      integer :: nameLength

      procedure (mpas_filter_list_compute_function), pointer :: computeFunc

      integer :: err_local

      iErr = MPAS_FILTER_NOERR

      if ( present(filterList) ) then
         filterListLocal => filterList
      else
         call mpas_pool_create_pool(filterListLocal)
         call mpas_filter_list_get_filter_names(inList, filterListLocal)
      end if

      if ( present(timeLevel) ) then
         timeLevelLocal = timeLevel
      else
         timeLevelLocal = 1
      end if

      call mpas_pool_begin_iteration(filterListLocal)
      do while ( mpas_pool_get_next_member(filterListLocal, poolItr) )

         if ( poolItr % memberType == MPAS_POOL_CONFIG ) then
            nameLength = len_trim(poolItr % memberName)
            configName = 'config_' // poolItr % memberName(1:nameLength) // '_enable'
            call mpas_pool_get_config(domain % configs, configName, filterEnable)

            if ( .not. associated(filterEnable) ) then
               call mpas_dmpar_global_abort('ERROR: Namelist option ''' // trim(configName) // ''' not found for filter ''' &
                                            // poolItr % memberName(1:nameLength) // '''')
            end if

#ifdef MPAS_TIME_FILTER_DRIVER
            timerName = computeTimerPrefix // poolItr % memberName(1:nameLength)
#endif
            alarmName = poolItr % memberName(1:nameLength) // computeAlarmSuffix
   
            if ( filterEnable ) then
               configName = 'config_' // poolItr % memberName(1:nameLength) // '_compute_interval'
               call mpas_pool_get_config(domain % configs, configName, filterComputeInterval)
               if ( .not. associated(filterComputeInterval) ) then
                  call mpas_dmpar_global_abort('ERROR: Namelist option ''' // trim(configName) // ''' not found for filter ''' &
                                               // poolItr % memberName(1:nameLength) // '''')
               end if
   
               configName = 'config_' // poolItr % memberName(1:nameLength) // '_stream_name'
               call mpas_pool_get_config(domain % configs, configName, filterStreamName)
               if ( .not. associated(filterStreamName) ) then
                  call mpas_dmpar_global_abort('ERROR: Namelist option ''' // trim(configName) // ''' not found for filter ''' &
                                               // poolItr % memberName(1:nameLength) // '''')
               end if
   
               call mpas_filter_list_get_filter_compute(inList, poolItr % memberName(1:nameLength), computeFunc, iErr=err_local)
   
               if ( filterComputeInterval == 'output_interval' .and. filterStreamName /= 'none' ) then
                  if ( mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID=filterStreamName, &
                                                      direction=MPAS_STREAM_OUTPUT, iErr=err_local) ) then
#ifdef MPAS_TIME_FILTER_DRIVER
                     call mpas_timer_start(timerName, .false.)
#endif
                     err_local = computeFunc(domain, poolItr % memberName(1:nameLength), timeLevelLocal)
#ifdef MPAS_TIME_FILTER_DRIVER
                     call mpas_timer_stop(timerName, .false.)
#endif
                  end if
               else if ( filterComputeInterval == 'dt' ) then
#ifdef MPAS_TIME_FILTER_DRIVER
                  call mpas_timer_start(timerName, .false.)
#endif
                  err_local = computeFunc(domain, poolItr % memberName(1:nameLength), timeLevelLocal)
#ifdef MPAS_TIME_FILTER_DRIVER
                  call mpas_timer_stop(timerName, .false.)
#endif

               else if ( mpas_is_alarm_ringing(domain % clock, alarmName, iErr=err_local) ) then
                  call mpas_reset_clock_alarm(domain % clock, alarmName, iErr=err_local)
#ifdef MPAS_TIME_FILTER_DRIVER
                  call mpas_timer_start(timerName, .false.)
#endif
                  err_local = computeFunc(domain, poolItr % memberName(1:nameLength), timeLevelLocal)
#ifdef MPAS_TIME_FILTER_DRIVER
                  call mpas_timer_stop(timerName, .false.)
#endif
               end if
            end if
         else if ( poolItr % memberType == MPAS_POOL_SUBPOOL ) then
            call mpas_pool_get_subpool(filterListLocal, poolItr % memberName, subPool)
            call mpas_filter_driver_compute(inList, domain, timeLevelLocal, subPool, iErr=err_local)
         end if
      end do

      if ( .not. present(filterList) ) then
         call mpas_pool_destroy_pool(filterListLocal)
      end if
   end subroutine mpas_filter_driver_compute!}}}


!***********************************************************************
!
!  routine mpas_filter_driver_restart
!
!> \brief   MPAS Filters restart routine
!> \author  Doug Jacobsen
!> \date    07/12/2015
!> \details
!> This routine iterates over a list of filers, and calls each filter's restart function.
!> The optional argument 'filterList' can be used to act on a subset of filters
!> defined in the inList mpas_filter_list_type.
!
!-----------------------------------------------------------------------
   recursive subroutine mpas_filter_driver_restart(inList, domain, filterList, iErr)!{{{
      type (mpas_filter_list_type), pointer :: inList
      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer, optional :: filterList
      integer, intent(out) :: iErr

      type (mpas_pool_type), pointer :: filterListLocal
      type (mpas_pool_type), pointer :: subPool

      character (len=StrKIND) :: configName
      logical, pointer :: filterEnable
      integer :: nameLength
      type (mpas_pool_iterator_type) :: poolItr

      procedure (mpas_filter_list_restart_function), pointer :: restartFunc

      integer :: err_local

      iErr = MPAS_FILTER_NOERR

      if ( present(filterList) ) then
         filterListLocal => filterList
      else
         call mpas_pool_create_pool(filterListLocal)
         call mpas_filter_list_get_filter_names(inList, filterListLocal)
      end if

      call mpas_pool_begin_iteration(filterListLocal)
      do while ( mpas_pool_get_next_member(filterListLocal, poolItr) )
         if ( poolItr % memberType == MPAS_POOL_CONFIG ) then
            nameLength = len_trim(poolItr % memberName)
            configName = 'config_' // poolItr % memberName(1:nameLength) // '_enable'
            call mpas_pool_get_config(domain % configs, configname, filterEnable)
   
            if ( .not. associated(filterEnable) ) then
               call mpas_dmpar_global_abort('ERROR: Namelist option ''' // trim(configName) // ''' not found for filter ''' &
                                            // poolItr % memberName(1:nameLength) // '''')
            end if
   
            if ( filterEnable ) then
#ifdef MPAS_TIME_FILTER_DRIVER
               timerName = restartTimerPrefix // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName, .false.)
#endif

               call mpas_filter_list_get_filter_restart(inList, poolItr % memberName(1:nameLength), restartFunc, iErr=err_local)
               err_local = restartFunc(domain, poolItr % memberName(1:nameLength))

#ifdef MPAS_TIME_FILTER_DRIVER
               call mpas_timer_stop(timerName, .false.)
#endif
            end if
         else if ( poolItr % memberType == MPAS_POOL_SUBPOOL ) then
            call mpas_pool_get_subpool(filterListLocal, poolItr % memberName, subPool)
            call mpas_filter_driver_restart(inList, domain, subPool, iErr=err_local)
         end if
      end do

      if ( .not. present(filterList) ) then
         call mpas_pool_destroy_pool(filterListLocal)
      end if
   end subroutine mpas_filter_driver_restart!}}}


!***********************************************************************
!
!  routine mpas_filter_driver_finalize
!
!> \brief   MPAS Filters finalize routine
!> \author  Doug Jacobsen
!> \date    07/12/2015
!> \details
!> This routine iterates over a list of filers, and calls each filter's finalize function.
!> The optional argument 'filterList' can be used to act on a subset of filters
!> defined in the inList mpas_filter_list_type.
!
!-----------------------------------------------------------------------
   recursive subroutine mpas_filter_driver_finalize(inList, domain, filterList, iErr)!{{{
      type (mpas_filter_list_type), pointer :: inList
      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer, optional :: filterList
      integer, intent(out) :: iErr

      type (mpas_pool_type), pointer :: filterListLocal
      type (mpas_pool_type), pointer :: subPool
      type (mpas_pool_iterator_type) :: poolItr
      character (len=StrKIND) :: configName, timerName
      logical, pointer :: filterEnable
      integer :: nameLength

      procedure (mpas_filter_list_finalize_function), pointer :: finalizeFunc

      integer :: err_local

      iErr = MPAS_FILTER_NOERR

      if ( present(filterList) ) then
         filterListLocal => filterList
      else
         call mpas_pool_create_pool(filterListLocal)
         call mpas_filter_list_get_filter_names(inList, filterListLocal)
      end if

      call mpas_pool_begin_iteration(filterListLocal)
      do while ( mpas_pool_get_next_member(filterListLocal, poolItr) )
         if ( poolItr % memberType == MPAS_POOL_CONFIG ) then
            nameLength = len_trim(poolItr % memberName)
            configName = 'config_' // poolItr % memberName(1:nameLength) // '_enable'
            call mpas_pool_get_config(domain % configs, configName, filterEnable)
   
            if ( .not. associated(filterEnable) ) then
               call mpas_dmpar_global_abort('ERROR: Namelist option ''' // trim(configName) // ''' not found for filter ''' &
                                            // poolItr % memberName(1:nameLength) // '''')
            end if
   
            if ( filterEnable ) then
#ifdef MPAS_TIME_FILTER_DRIVER
               timerName = finalizeTimerPrefix // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName, .false.)
#endif

               call mpas_filter_list_get_filter_finalize(inList, poolItr % memberName(1:nameLength), finalizeFunc, iErr=err_local)
               err_local = finalizeFunc(domain, poolItr % memberName(1:nameLength))

#ifdef MPAS_TIME_FILTER_DRIVER
               call mpas_timer_stop(timerName)
#endif
            end if
         else if ( poolItr % memberType == MPAS_POOL_SUBPOOL ) then
            call mpas_pool_get_subpool(filterListLocal, poolItr % memberName, subPool)
            call mpas_filter_driver_finalize(inList, domain, subPool, iErr=err_local)
         end if
      end do

      if ( .not. present(filterList) ) then
         call mpas_pool_destroy_pool(filterListLocal)
      end if

   end subroutine mpas_filter_driver_finalize!}}}


end module mpas_filter_driver
