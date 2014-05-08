! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_paraview_catalyst
!
!> \brief MPAS ocean analysis core member: paraview_catalyst
!> \author Mark Petersen
!> \date   May 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module ocn_paraview_catalyst

   use mpas_grid_types
   use mpas_timer
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_paraview_catalyst

   use ocn_constants
   use ocn_diagnostics_routines

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_setup_packages_paraview_catalyst, &
             ocn_init_paraview_catalyst, &
             ocn_init_alarms_paraview_catalyst, &
             ocn_compute_paraview_catalyst, &
             ocn_write_paraview_catalyst, &
             ocn_restart_paraview_catalyst, &
             ocn_finalize_paraview_catalyst

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: amParaviewCatalystTimer
   integer :: imageFrame

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_setup_packages_paraview_catalyst
!
!> \brief   Set up packages for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    May 2014
!> \details 
!>  This routine is intended to configure the packages for this MPAS
!>   ocean analysis member
!
!-----------------------------------------------------------------------

   subroutine ocn_setup_packages_paraview_catalyst(configPool, packagePool, err)!{{{

      use mpas_packages

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      logical, pointer :: amParaviewCatalystActive

      err = 0

      call mpas_pool_get_package(packagePool, 'amParaviewCatalystActive', amParaviewCatalystActive)

      ! turn on package for this analysis member
      amParaviewCatalystActive = .true.

   end subroutine ocn_setup_packages_paraview_catalyst!}}}


!***********************************************************************
!
!  routine ocn_init_paraview_catalyst
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    May 2014
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_paraview_catalyst(domain, err)!{{{

      use mpas_packages

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

      call mpas_init_paraview_catalyst
      call mpas_init_mesh_paraview_catalyst(domain)
      call mpas_init_data_paraview_catalyst(domain)

      imageFrame = 1

   end subroutine ocn_init_paraview_catalyst!}}}

!***********************************************************************
!
!  routine ocn_init_alarms_paraview_catalyst
!
!> \brief   Initialize alarms for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    May 2014
!> \details 
!>  This routine conducts all alarm initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_alarms_paraview_catalyst(clock, startTime, computeAlarmID, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: computeAlarmID
      type (MPAS_Time_Type), intent(in) :: startTime

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (MPAS_Clock_type), intent(inout) :: clock

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: ierr
      character(len=StrKIND) :: compute_interval
      type (MPAS_Time_Type) :: alarmStartTime
      type (MPAS_TimeInterval_type) :: alarmTimeStep

      character (len=StrKIND), pointer :: config_paraview_catalyst_compute_interval, config_output_interval

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_paraview_catalyst_compute_interval', config_paraview_catalyst_compute_interval)
      call mpas_pool_get_config(ocnConfigs, 'config_output_interval', config_output_interval)

      ! set compute alarm for this analysis member
      if (config_paraview_catalyst_compute_interval=='same_as_output') then
         compute_interval = config_output_interval
      else
         compute_interval = config_paraview_catalyst_compute_interval
      endif

      call mpas_set_timeInterval(alarmTimeStep, timeString=compute_interval, ierr=ierr)
      alarmStartTime = startTime + alarmTimeStep
      call mpas_add_clock_alarm(clock, computeAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)

   end subroutine ocn_init_alarms_paraview_catalyst!}}}

!***********************************************************************
!
!  routine ocn_compute_paraview_catalyst
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    May 2014
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_paraview_catalyst(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: doWork, fileID

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: diagnosticsPool
      character (len=StrKIND), pointer :: xtime

      real (kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge

      err = 0

      dminfo = domain % dminfo

      call mpas_timer_start("compute_paraview_catalyst", .false., amParaviewCatalystTimer)

      doWork = 0
      call mpas_check_coprocess(imageFrame, doWork)
      if(doWork .eq. 0) then
         return
      endif

      call mpas_load_data_paraview_catalyst(domain, imageFrame)
      call mpas_coprocess

      if (dminfo % my_proc_id == IO_NODE) then
         call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
         fileID = getFreeUnit()
         open(fileID,file='paraview_catalyst_time.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(i10,10x,a,100es24.14)') imageFrame, trim(xtime)
         close (fileID)
      end if

      imageFrame = imageFrame + 1

      call mpas_timer_stop("paraview_catalyst", amParaviewCatalystTimer)

   end subroutine ocn_compute_paraview_catalyst!}}}

!***********************************************************************
!
!  routine ocn_restart_paraview_catalyst
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    May 2014
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_paraview_catalyst(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_paraview_catalyst!}}}

!***********************************************************************
!
!  routine ocn_write_paraview_catalyst
!
!> \brief   MPAS-Ocean analysis output
!> \author  Mark Petersen
!> \date    May 2014
!> \details 
!>  This routine writes all output for this MPAS-Ocean analysis member.
!>  At this time this is just a stub, and all analysis output is written
!>  to the output file specified by config_output_name.
!
!-----------------------------------------------------------------------

   subroutine ocn_write_paraview_catalyst(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_write_paraview_catalyst!}}}

!***********************************************************************
!
!  routine ocn_finalize_paraview_catalyst
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    May 2014
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_paraview_catalyst(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

      call mpas_finalize_paraview_catalyst

   end subroutine ocn_finalize_paraview_catalyst!}}}

   integer function getFreeUnit()!{{{
      implicit none

      integer :: index
      logical :: isOpened

      getFreeUnit = 0
      do index = 1,99
         if((index /= 5) .and. (index /= 6)) then
            inquire(unit = index, opened = isOpened)
            if( .not. isOpened) then
               getFreeUnit = index
               return
            end if
         end if
      end do
   end function getFreeUnit!}}}

end module ocn_paraview_catalyst

! vim: foldmethod=marker
