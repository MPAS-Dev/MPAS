! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_TEM_PLATE
!
!> \brief MPAS ocean analysis mode member: TEM_PLATE
!> \author FILL_IN_AUTHOR
!> \date   FILL_IN_DATE
!> \details
!>  MPAS ocean analysis mode member: TEM_PLATE
!>  In order to add a new analysis member, do the following:
!>  1. Copy these to your new analysis member name:
!>     cp mpas_ocn_TEMPLATE.F mpas_ocn_your_new_name.F
!>     cp Registry_ocn_TEMPLATE.xml Registry_ocn_your_new_name.xml
!>
!>  2. In those two new files, replace the following text:
!>     temPlate, TEM_PLATE, FILL_IN_AUTHOR, FILL_IN_DATE
!>     Typically temPlate uses camel case (variable names), like yourNewName,
!>     while TEM_PLATE uses underscores (subroutine names), like your_new_name.
!>     note: do not replace 'filename_template' in Registry_ocn_yourNewName.xml
!>
!>  3. Add a #include line for your registry to
!>     Registry_analysis_members.xml
!>
!>  4. In mpas_ocn_analysis_driver.F, add a use statement for your new analysis member.
!>     In addition, add lines for your analysis member, and replace TEM_PLATE
!>     and temPlate as described in step 2. There should be 5 places that need additions:
!>      - Adding the analysis member name to the analysis member list
!>      - Adding an init if test can subroutine call
!>      - Adding a compute if test can subroutine call
!>      - Adding a restart if test can subroutine call
!>      - Adding a finalize if test can subroutine call
!>
!>  5. In src/core_ocean/analysis_members/Makefile, add your
!>     new analysis member to the list of members. See another analysis member
!>     in that file for an example.
!>     NOTE: If your analysis member depends on other files, add a dependency
!>           line for the member and list them there. See okubo weiss for an example.
!>
!-----------------------------------------------------------------------

module ocn_moc_streamfunction

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

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

   public :: ocn_init_moc_streamfunction, &
             ocn_compute_moc_streamfunction, &
             ocn_restart_moc_streamfunction, &
             ocn_finalize_moc_streamfunction

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   
   integer :: nMocStreamfunctionBinsUsed

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_moc_streamfunction
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Nils H. Feige
!> \date    2016-04-08
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_moc_streamfunction(domain, err)!{{{

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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mocStreamfunctionAMPool
      type (mpas_pool_type), pointer :: meshPool

      integer ::  iBin
     
      real (kind=RKIND) :: binWidth
      ! These are array size 1 because mpas_dmpar_min_real_array calls require arrays.
      real (kind=RKIND), dimension(1) :: minBin, maxBin, minBinDomain, maxBinDomain
      ! the variable used to discriminate cells into Bins (either the y-value or the latitude)
      real (kind=RKIND), dimension(:), pointer :: binVariable, binBoundaryMocStreamfunction

      !number of latitude bins specified in the config
      integer, pointer :: config_AM_mocStreamfunction_num_bins
      !smallest and highest latitude specified in the config
      real (kind=RKIND), pointer :: config_AM_mocStreamfunction_min_bin, config_AM_mocStreamfunction_max_bin

      !determines if the simulation was run on a sphere or on a plane
      logical, pointer :: on_a_sphere

      dminfo = domain % dminfo

      err = 0

      minBin =  1.0e34_RKIND
      maxBin = -1.0e34_RKIND

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)

      call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_num_bins', &
                                config_AM_mocStreamfunction_num_bins)
      call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_min_bin', &
                                config_AM_mocStreamfunction_min_bin)
      call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_max_bin', &
                                config_AM_mocStreamfunction_max_bin)

      nMocStreamfunctionBinsUsed = config_AM_mocStreamfunction_num_bins

      call mpas_pool_get_array(mocStreamfunctionAMPool, 'binBoundaryMocStreamfunction', binBoundaryMocStreamfunction)

      ! Find min and max values of binning variable.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

         ! Bin by latitude on a sphere, by yCell otherwise.
         if (on_a_sphere) then
            call mpas_pool_get_array(meshPool, 'latCell', binVariable)
         else
            call mpas_pool_get_array(meshPool, 'yCell', binVariable)
         end if

         minBin = min(minBin, minval(binVariable) )
         maxBin = max(maxBin, maxval(binVariable) )

         block => block % next
      end do

      call mpas_dmpar_min_real_array(dminfo, 1, minBin, minBinDomain)
      call mpas_dmpar_max_real_array(dminfo, 1, maxBin, maxBinDomain)

      ! Set up bins.
      binBoundaryMocStreamfunction = -1.0e34_RKIND

      ! Change min and max bin bounds to configuration settings, if applicable.
      if (config_AM_mocStreamfunction_min_bin > -1.0e33_RKIND) then
         minBinDomain(1) = config_AM_mocStreamfunction_min_bin
      else
         ! use measured min value, but decrease slightly to include least value.
         minBinDomain(1) = minBinDomain(1) - 1.0e-10_RKIND * abs(minBinDomain(1))
      end if

      if (config_AM_mocStreamfunction_max_bin > -1.0e33_RKIND) then
         maxBinDomain(1) = config_AM_mocStreamfunction_max_bin
      else
         ! use measured max value, but increase slightly to include max value.
         maxBinDomain(1) = maxBinDomain(1) + 1.0e-10_RKIND * abs(maxBinDomain(1))
      end if

      binBoundaryMocStreamfunction(1) = minBinDomain(1)
      binWidth = (maxBinDomain(1) - minBinDomain(1)) / nMocStreamfunctionBinsUsed

      do iBin = 2, nMocStreamfunctionBinsUsed
         binBoundaryMocStreamfunction(iBin) = binBoundaryMocStreamfunction(iBin-1) + binWidth
      end do
      binBoundaryMocStreamfunction(nMocStreamfunctionBinsUsed+1) = binBoundaryMocStreamfunction(nMocStreamfunctionBinsUsed) + binWidth

   end subroutine ocn_init_moc_streamfunction!}}}

!***********************************************************************
!
!  routine ocn_compute_moc_streamfunction
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Nils H. Feige
!> \date    2016-04-08
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_moc_streamfunction(domain, timeLevel, err)!{{{
	  implicit none
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

      type (mpas_pool_type), pointer :: mocStreamfunctionAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: meshPool
      
      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nCellsSolve
      integer :: iCell, iBin, genericCounter, k
      real (kind=RKIND) :: binWidth
      real (kind=RKIND), dimension(:,:), pointer :: mocStreamValLatAndDepthLocal
	  real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer ::  areaCell, binBoundaryMocStreamfunction
      real (kind=RKIND), dimension(:,:), pointer :: mocStreamvalLatAndDepth, mocStreamValLatAndDepthTotal
      real (kind=RKIND), dimension(:,:), pointer :: vertVelocityTop
      real (kind=RKIND), dimension(:,:), pointer :: sumVertBinVelocity
      
      character (len=30) :: verticalVelocityArrayName
      integer, pointer :: verticalVelocityArrayNumber

      err = 0
      
      dminfo = domain % dminfo
      
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)
      
      allocate(mocStreamValLatAndDepthLocal(nMocStreamfunctionBinsUsed, nVertLevels + 1))
      allocate(sumVertBinVelocity(nMocStreamfunctionBinsUsed, nVertLevels))
      allocate(mocStreamValLatAndDepthTotal(nMocStreamfunctionBinsUsed, nVertLevels + 1))
      
      do iBin = 1, nMocStreamfunctionBinsUsed
         mocStreamValLatAndDepthLocal(iBin,:) = 0.0_RKIND
         mocStreamValLatAndDepthTotal(iBin,:) = 0.0_RKIND
         sumVertBinVelocity(iBin, :) = 0.0_RKIND
      end do

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         call mpas_pool_get_array(mocStreamfunctionAMPool, 'binBoundaryMocStreamfunction', binBoundaryMocStreamfunction)

      	 call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_velocity_value', &
                                verticalVelocityArrayNumber)
                                
      	 !1 for 'vertVelocityTop', 2 for 'vertTransportVelocityTop' or 3 for 'vertGMBolusVelocityTop'                          
      	 if (verticalVelocityArrayNumber == 1) then
      	    verticalVelocityArrayName = 'vertVelocityTop'
      	 else if (verticalVelocityArrayNumber == 2) then
      	    verticalVelocityArrayName = 'vertTransportVelocityTop'
     	 else if (verticalVelocityArrayNumber == 3) then
     	    verticalVelocityArrayName = 'vertGMBolusVelocityTop'
     	 end if   

         binWidth = (binBoundaryMocStreamfunction(nMocStreamfunctionBinsUsed + 1) - binBoundaryMocStreamfunction(1)) & 
         / nMocStreamfunctionBinsUsed
         
         !call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere) might come in handy

         ! Here are some example variables which may be needed for your analysis member

         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(diagnosticsPool, verticalVelocityArrayName, vertVelocityTop)
		 call mpas_pool_get_array(meshPool, 'latCell', latCell)
		 !not needed for global moc, but later on for zonal moc
		 call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

         do iCell = 1,nCellsSolve
			do iBin = 1, nMocStreamfunctionBinsUsed
               if (latCell(iCell) .gt. binBoundaryMocStreamfunction(iBin) .and. latCell(iCell) .lt. binBoundaryMocStreamfunction(iBin + 1)) then
                  sumVertBinVelocity(iBin, :) = sumVertBinVelocity(iBin, :) + (vertVelocityTop(:, iCell) * areaCell(iCell))
               end if
            end do
         end do
         
         do iBin = 2, nMocStreamfunctionBinsUsed
            do k = 1, nVertLevels
               mocStreamValLatAndDepthLocal(iBin, k) = mocStreamValLatAndDepthLocal(iBin-1, k) &
                  + sumVertBinVelocity(iBin, k)
            end do
         end do  

         block => block % next
     end do
      
     call mpas_dmpar_sum_real_array(dminfo, (nVertLevels+1) * nMocStreamfunctionBinsUsed, mocStreamValLatAndDepthLocal, &
         mocStreamvalLatAndDepthTotal)
         
!     print *, 'sumVertBinVelocitry', sumVertBinVelocity
!     print *, 'mocStreamvalLatAndDepthLocal', mocStreamvalLatAndDepthLocal
!     print *, 'mocStreamvalLatAndDepthTotal', mocStreamvalLatAndDepthTotal
     
     !block => domain % blocklist
     !do while (associated(block))
     !   call mpas_pool_get_subpool(block % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)
     !   call mpas_pool_get_array(mocStreamfunctionAMPool, 'mocStreamvalLatAndDepth', mocStreamvalLatAndDepth)
     !    
     !   do iBin = 1, nMocStreamfunctionBinsUsed
     !      do k = 1, nVertLevels + 1
     !         mocStreamvalLatAndDepth(iBin, k) = mocStreamvalLatAndDepthTotal(iBin, k)
     !      end do
     !   end do
     !    
     !   block => block % next
     !end do
     
     
     
     call mpas_pool_get_subpool(domain % blocklist % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)
     call mpas_pool_get_array(mocStreamfunctionAMPool, 'mocStreamvalLatAndDepth', mocStreamvalLatAndDepth)
     mocStreamvalLatAndDepth = mocStreamvalLatAndDepthTotal / 1000000_RKIND
!     print *, 'mocStreamvalLatAndDepth', mocStreamvalLatAndDepth
     deallocate(mocStreamvalLatAndDepthTotal)
     deallocate(mocStreamvalLatAndDepthLocal)
     deallocate(sumVertBinVelocity)

   end subroutine ocn_compute_moc_streamfunction!}}}

!***********************************************************************
!
!  routine ocn_restart_moc_streamfunction
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Nils H. Feige
!> \date    2016-04-08
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_moc_streamfunction(domain, err)!{{{

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

   end subroutine ocn_restart_moc_streamfunction!}}}

!***********************************************************************
!
!  routine ocn_finalize_moc_streamfunction
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Nils H. Feige
!> \date    2016-04-08
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_moc_streamfunction(domain, err)!{{{

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

   end subroutine ocn_finalize_moc_streamfunction!}}}

end module ocn_moc_streamfunction