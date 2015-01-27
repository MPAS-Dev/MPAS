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
             ocn_compute_paraview_catalyst, &
             ocn_write_paraview_catalyst, &
             ocn_restart_paraview_catalyst, &
             ocn_finalize_paraview_catalyst

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: paraview_catalystTimer
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
      logical, pointer :: paraview_catalystActive

      err = 0

      call mpas_pool_get_package(packagePool, 'paraview_catalystActive', paraview_catalystActive)

      ! turn on package for this analysis member
      paraview_catalystActive = .true.
!mrp
print *, 'ocn_setup_packages_paraview_catalyst'

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

!mrp
print *, 'ocn_init_paraview_catalyst'

      imageFrame = 1

   end subroutine ocn_init_paraview_catalyst!}}}

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

      call mpas_timer_start("compute_paraview_catalyst", .false., paraview_catalystTimer)

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

      call mpas_timer_stop("paraview_catalyst", paraview_catalystTimer)
!mrp
print *, 'ocn_compute_paraview_catalyst'
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

!***********************************************************************
!
!  routine mpas_init_data_paraview_catalyst
!
!> \brief   Create mesh variables for paraview catalyst
!> \author  Patricia Fasel, Mark Petersen
!> \date    April 2014
!> \details 
!>  Create all the mesh information needed by the paraview catalyst
!>  for in-situ visualization and analysis.
!
!-----------------------------------------------------------------------
   subroutine mpas_init_data_paraview_catalyst(domain)!{{{
   
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, pointer :: indexTemperature, indexSalinity, num_tracers
      integer :: timeLevel

      type (block_type), pointer :: block
      type (field2DReal), pointer :: densityField, kineticEnergyCellField, relativeVorticityCellField
      type (field3DReal), pointer :: tracersField

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: diagnosticsPool

      timeLevel=1  ! mrp change to input argument later

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)
         call mpas_pool_get_field(statePool, 'tracers', tracersField, timeLevel)
         call mpas_pool_get_dimension(statePool, 'index_temperature', indexTemperature)

         ! subtract one from dimsSizes to use nCells rather than nCells+1
         call coprocessor_register_tracer_data( \
                       indexTemperature, \
                       tracersField % constituentNames(indexTemperature), \
                       num_tracers, \
                       tracersField % dimSizes(2), \
                       tracersField % dimSizes(3)-1, & 
                       tracersField % array)
         call mpas_pool_get_dimension(statePool, 'index_salinity', indexSalinity)
         call coprocessor_register_tracer_data( \
                       indexSalinity, \
                       tracersField % constituentNames(indexSalinity), \
                       num_tracers, \
                       tracersField % dimSizes(2), \
                       tracersField % dimSizes(3)-1, & 
                       tracersField % array)

         call mpas_pool_get_field(diagnosticsPool, 'density', densityField)
         call coprocessor_register_data( \
                       densityField % fieldName, \
                       densityField % dimSizes(1), \
                       densityField % dimSizes(2)-1, & 
                       densityField % array)

         call mpas_pool_get_field(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCellField)
         call coprocessor_register_data( \
                       kineticEnergyCellField % fieldName, \
                       kineticEnergyCellField % dimSizes(1), \
                       kineticEnergyCellField % dimSizes(2)-1, & 
                       kineticEnergyCellField % array)

         call mpas_pool_get_field(diagnosticsPool, 'relativeVorticityCell', relativeVorticityCellField)
         call coprocessor_register_data( \
                       relativeVorticityCellField % fieldName, \
                       relativeVorticityCellField % dimSizes(1), \
                       relativeVorticityCellField % dimSizes(2)-1, & 
                       relativeVorticityCellField % array)

       block => block % next
    end do

!mrp
print *, 'mpas_init_data_paraview_catalyst'
   end subroutine mpas_init_data_paraview_catalyst!}}}


!***********************************************************************
!
!  routine mpas_load_data_paraview_catalyst
!
!> \brief   Pass data to paraview catalyst
!> \author  Patricia Fasel, Mark Petersen
!> \date    April 2014
!> \details 
!>  Pass data to paraview catalyst
!
!-----------------------------------------------------------------------
   subroutine mpas_load_data_paraview_catalyst(domain, itime)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (domain_type), intent(in) :: domain
      integer, intent(in) :: itime
!      integer, intent(in) :: timeLevel  ! mrp add later

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      ! Dimensions
      integer, pointer :: indexTemperature, indexSalinity, num_tracers

      type (block_type), pointer :: block
      integer :: timeLevel
      type (field2DReal), pointer :: densityField, kineticEnergyCellField, relativeVorticityCellField
      type (field3DReal), pointer :: tracersField

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: diagnosticsPool

timeLevel=1  ! mrp change to input argument later

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)
         call mpas_pool_get_field(statePool, 'tracers', tracersField, timeLevel)
         call mpas_pool_get_dimension(statePool, 'index_temperature', indexTemperature)

         ! subtract one from dimsSizes to use nCells rather than nCells+1
         call coprocessor_add_tracer_data( \
                       itime, \
                       indexTemperature, \
                       tracersField % constituentNames(indexTemperature), \
                       num_tracers, \
                       tracersField % dimSizes(2), \
                       tracersField % dimSizes(3)-1, & 
                       tracersField % array)
         call mpas_pool_get_dimension(statePool, 'index_salinity', indexSalinity)
         call coprocessor_add_tracer_data( \
                       itime, \
                       indexSalinity, \
                       tracersField % constituentNames(indexSalinity), \
                       num_tracers, \
                       tracersField % dimSizes(2), \
                       tracersField % dimSizes(3)-1, & 
                       tracersField % array)

         call mpas_pool_get_field(diagnosticsPool, 'density', densityField)
         call coprocessor_add_data( \
                       itime, \
                       densityField % fieldName, \
                       densityField % dimSizes(1), \
                       densityField % dimSizes(2)-1, & 
                       densityField % array)

         call mpas_pool_get_field(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCellField)
         call coprocessor_add_data( \
                       itime, \
                       kineticEnergyCellField % fieldName, \
                       kineticEnergyCellField % dimSizes(1), \
                       kineticEnergyCellField % dimSizes(2)-1, & 
                       kineticEnergyCellField % array)

         call mpas_pool_get_field(diagnosticsPool, 'relativeVorticityCell', relativeVorticityCellField)
         call coprocessor_add_data( \
                       itime, \
                       relativeVorticityCellField % fieldName, \
                       relativeVorticityCellField % dimSizes(1), \
                       relativeVorticityCellField % dimSizes(2)-1, & 
                       relativeVorticityCellField % array)

     block => block % next
   end do

!mrp
print *, 'mpas_load_data_paraview_catalyst'

   end subroutine mpas_load_data_paraview_catalyst!}}}

end module ocn_paraview_catalyst

! vim: foldmethod=marker
