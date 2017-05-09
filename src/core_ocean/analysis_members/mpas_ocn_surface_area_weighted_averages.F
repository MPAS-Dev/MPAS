! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_surface_area_weighted_averages
!
!> \brief MPAS ocean analysis member: areal-{min,max,avg} of 2D fields
!> \author Todd Ringler, Anne Berres
!> \date   Mar 2, 2017
!> \details
!>  MPAS ocean analysis member: surface-area-weighted averages
!
!-----------------------------------------------------------------------

module ocn_surface_area_weighted_averages

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

   public :: ocn_init_surface_area_weighted_averages, &
             ocn_compute_surface_area_weighted_averages, &
             ocn_restart_surface_area_weighted_averages, &
             ocn_finalize_surface_area_weighted_averages

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains


!***********************************************************************
!
!  routine ocn_init_surface_area_weighted_averages
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_surface_area_weighted_averages(domain, err)!{{{

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
      type (mpas_pool_type), pointer :: surfaceAreaWeightedAveragesAMPool
      type (mpas_pool_type), pointer :: meshPool
      ! pointers to data in mesh pool
      integer, pointer :: nCells, nCellsSolve, nSfcAreaWeightedAvgFields, nOceanRegions
      real (kind=RKIND), dimension(:), pointer :: lonCell, latCell
      logical, pointer :: predefRegions
      real (kind=RKIND), dimension(:,:), pointer :: surfaceMask

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
      logical, pointer :: predefReg

      err = 0

      ! get pointers to analysis member arrays
      call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceAreaWeightedAveragesAM', surfaceAreaWeightedAveragesAMPool)

      ! region config for SAWA
      call mpas_pool_get_config(domain % configs, 'config_AM_surfaceAreaWeightedAverages_compute_predef_regions', &
                                   predefReg)

      call mpas_pool_get_array(surfaceAreaWeightedAveragesAMPool, 'surfaceMask', surfaceMask)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nOceanRegions', nOceanRegions)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)

      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)

      ! compute mask
      if(predefReg == .true.) then
         call compute_mask(nCells, nCellsSolve, nOceanRegions, lonCell, latCell, surfaceMask)
      endif


   contains

   subroutine compute_mask(nCells, nCellsSolve, nOceanRegions, lonCell, latCell, surfaceMask)!{{{
   ! this subroutines produces a 0/1 mask that is multiplied with workArray to
   ! allow for min/max/avg to represent specific regions of the ocean domain
   !
   ! NOTE: computes_mask is temporary. workMask should be intent(in) to this entire module !
   !
   integer, intent(in) :: nCells, nCellsSolve
   real(kind=RKIND), dimension(:), intent(in) :: lonCell, latCell
   real(kind=RKIND), dimension(:,:), intent(out) :: surfaceMask
   integer, pointer :: nOceanRegions
   integer :: iCell
   real(kind=RKIND) :: dtr

   dtr = 4.0_RKIND*atan(1.0_RKIND) / 180.0_RKIND
   surfaceMask(:,:) = 0.0_RKIND
   do iCell=1,nCellsSolve
      surfaceMask(:,iCell) = 1.0_RKIND
   enddo

   ! Arctic
   do iCell=1,nCellsSolve
      if(latCell(iCell).lt. 60.0_RKIND*dtr) surfaceMask(1,iCell) = 0.0_RKIND
   enddo
   ! Equatorial
   do iCell=1,nCellsSolve
      if(latCell(iCell).gt. 15.0_RKIND*dtr) surfaceMask(2,iCell) = 0.0_RKIND
      if(latCell(iCell).lt.-15.0_RKIND*dtr) surfaceMask(2,iCell) = 0.0_RKIND
   enddo
   ! Southern Ocean
   do iCell=1,nCellsSolve
      if(latCell(iCell).gt.-50.0_RKIND*dtr) surfaceMask(3,iCell) = 0.0_RKIND
   enddo
   ! Nino 3
   do iCell=1,nCellsSolve
      if(latCell(iCell).gt.  5.0_RKIND*dtr) surfaceMask(4,iCell) = 0.0_RKIND
      if(latCell(iCell).lt. -5.0_RKIND*dtr) surfaceMask(4,iCell) = 0.0_RKIND
      if(lonCell(iCell).lt.210.0_RKIND*dtr) surfaceMask(4,iCell) = 0.0_RKIND
      if(lonCell(iCell).gt.270.0_RKIND*dtr) surfaceMask(4,iCell) = 0.0_RKIND
   enddo
   ! Nino 4
   do iCell=1,nCellsSolve
      if(latCell(iCell).gt.  5.0_RKIND*dtr) surfaceMask(5,iCell) = 0.0_RKIND
      if(latCell(iCell).lt. -5.0_RKIND*dtr) surfaceMask(5,iCell) = 0.0_RKIND
      if(lonCell(iCell).lt.160.0_RKIND*dtr) surfaceMask(5,iCell) = 0.0_RKIND
      if(lonCell(iCell).gt.210.0_RKIND*dtr) surfaceMask(5,iCell) = 0.0_RKIND
   enddo
      ! Nino 3.4
   do iCell=1,nCellsSolve
      if(latCell(iCell).gt.  5.0_RKIND*dtr) surfaceMask(6,iCell) = 0.0_RKIND
      if(latCell(iCell).lt. -5.0_RKIND*dtr) surfaceMask(6,iCell) = 0.0_RKIND
      if(lonCell(iCell).lt.190.0_RKIND*dtr) surfaceMask(6,iCell) = 0.0_RKIND
      if(lonCell(iCell).gt.240.0_RKIND*dtr) surfaceMask(6,iCell) = 0.0_RKIND
   enddo

   end subroutine compute_mask!}}}

   end subroutine ocn_init_surface_area_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_compute_surface_area_weighted_averages
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Todd Ringler, Anne Berres
!> \date    Mar 2, 2017
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_surface_area_weighted_averages(domain, timeLevel, err)!{{{

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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: surfaceAreaWeightedAveragesAMPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: tracersSurfaceFluxPool
      type (mpas_pool_type), pointer :: regionPool

      !! Hard-wired Region variables
      real (kind=RKIND), dimension(:,:), pointer :: minValueWithinOceanRegion
      real (kind=RKIND), dimension(:,:), pointer :: maxValueWithinOceanRegion
      real (kind=RKIND), dimension(:,:), pointer :: avgValueWithinOceanRegion

      !! Region File variables
      real (kind=RKIND), dimension(:,:), pointer :: minValueWithinOceanRegionMask
      real (kind=RKIND), dimension(:,:), pointer :: maxValueWithinOceanRegionMask
      real (kind=RKIND), dimension(:,:), pointer :: avgValueWithinOceanRegionMask

      ! pointers to data in pools to be analyzed
      ! bulkForcing Pkg
      real (kind=RKIND), dimension(:), pointer :: latentHeatFlux
      real (kind=RKIND), dimension(:), pointer :: sensibleHeatFlux
      real (kind=RKIND), dimension(:), pointer :: longWaveHeatFluxUp
      real (kind=RKIND), dimension(:), pointer :: longWaveHeatFluxDown
      real (kind=RKIND), dimension(:), pointer :: seaIceHeatFlux
      real (kind=RKIND), dimension(:), pointer :: shortWaveHeatFlux
      real (kind=RKIND), dimension(:), pointer :: evaporationFlux
      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux
      real (kind=RKIND), dimension(:), pointer :: riverRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: iceRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: rainFlux
      real (kind=RKIND), dimension(:), pointer :: snowFlux

      ! frazilIce Pkg
      real (kind=RKIND), dimension(:), pointer :: seaIceEnergy

      real (kind=RKIND), dimension(:), pointer :: surfaceThicknessFlux
      real (kind=RKIND), dimension(:,:), pointer :: activeTracersSurfaceFlux
      real (kind=RKIND), dimension(:), pointer ::  penetrativeTemperatureFlux

      real (kind=RKIND), dimension(:), pointer :: seaIceSalinityFlux
      real (kind=RKIND), dimension(:), pointer :: surfaceStressMagnitude
      real (kind=RKIND), dimension(:), pointer :: windStressZonal
      real (kind=RKIND), dimension(:), pointer :: windStressMeridional
      real (kind=RKIND), dimension(:), pointer :: atmosphericPressure
      real (kind=RKIND), dimension(:), pointer :: ssh
      real (kind=RKIND), dimension(:), pointer :: boundaryLayerDepth
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      ! pointers to data in mesh pool
      integer, pointer :: nCells, nCellsSolve, nSfcAreaWeightedAvgFields, nOceanRegions
      integer, pointer :: indexTemperature, indexSalinity
      real (kind=RKIND), dimension(:), pointer ::  areaCell, lonCell, latCell
      logical, pointer :: predefRegions

      ! pointers to data in AM pool
      real (kind=RKIND), dimension(:,:), pointer :: surfaceMask

      ! scratch space
      type(field2DReal), pointer :: workArrayField, workArrayFieldRegionMask
      real (kind=RKIND), dimension(:,:), pointer :: workArray, workArrayRegionMask
      type(field2DReal), pointer :: workMinField, workMaxField, workSumField
      type(field2DReal), pointer :: workMinFieldRegionMask, workMaxFieldRegionMask, workSumFieldRegionMask
      real (kind=RKIND), dimension(:,:), pointer :: workMin, workMax, workSum
      real (kind=RKIND), dimension(:,:), pointer :: workMinRegionMask, workMaxRegionMask, workSumRegionMask

      ! local variables
      integer :: iDataField, nDefinedDataFields
      integer :: iCell, iRegion, iTracer, err_tmp

      ! package flag
      logical, pointer :: activeTracersBulkRestoringPKG
      logical, pointer :: frazilIcePkgActive

      ! buffers data for message passaging
      integer :: kBuffer, kBufferLength
      real (kind=RKIND), dimension(:), allocatable :: workBufferSum, workBufferSumReduced
      real (kind=RKIND), dimension(:), allocatable :: workBufferSumRegionMask, workBufferSumReducedRegionMask
      real (kind=RKIND), dimension(:), allocatable :: workBufferMin, workBufferMinReduced
      real (kind=RKIND), dimension(:), allocatable :: workBufferMinRegionMask, workBufferMinReducedRegionMask
      real (kind=RKIND), dimension(:), allocatable :: workBufferMax, workBufferMaxReduced
      real (kind=RKIND), dimension(:), allocatable :: workBufferMaxRegionMask, workBufferMaxReducedRegionMask

      !!! region file variables
      integer, pointer :: regionsInAddGroup
      integer :: regionGroupOffset, regionGroupNumber, i, j
      character (len=STRKIND), dimension(:), pointer :: regionNames, regionGroupNames
      integer, dimension(:,:), pointer :: regionCellMasks, regionsInGroup
      integer, dimension(:), pointer ::  nRegionsInGroup
      integer, pointer :: nRegions, nRegionGroups, maxRegionsInGroup
      character (len=STRKIND), pointer :: additionalRegion

      ! assume no error
      err = 0

      ! get status of other packages
      call mpas_pool_get_package(ocnPackages, 'activeTracersBulkRestoringPKGActive', activeTracersBulkRestoringPKG)
      call mpas_pool_get_package(ocnPackages, 'frazilIceActive', frazilIcePkgActive)

      ! set highest level pointer
      dminfo = domain % dminfo

      ! get pointers to analysis member arrays
      call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceAreaWeightedAveragesAM', surfaceAreaWeightedAveragesAMPool)

      ! region config for SAWA
      call mpas_pool_get_config(domain % configs, 'config_AM_surfaceAreaWeightedAverages_compute_predef_regions', &
                                   predefRegions)
      call mpas_pool_get_config(domain % configs, 'config_AM_surfaceAreaWeightedAverages_region_group', &
                                   additionalRegion)


      ! find the number of data fields
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nSfcAreaWeightedAvgFields', nSfcAreaWeightedAvgFields)


      if (predefRegions == .true.) then
         !!! hard-wired regions
         ! find the number of regions
         call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nOceanRegions', nOceanRegions)
         call mpas_pool_get_array(surfaceAreaWeightedAveragesAMPool, 'surfaceMask', surfaceMask)
         ! allocate buffer for message passing
         kBuffer=0
         kBufferLength=nOceanRegions*nSfcAreaWeightedAvgFields
         allocate(workBufferSum(kBufferLength))
         allocate(workBufferMin(kBufferLength))
         allocate(workBufferMax(kBufferLength))
         allocate(workBufferSumReduced(kBufferLength))
         allocate(workBufferMinReduced(kBufferLength))
         allocate(workBufferMaxReduced(kBufferLength))
         workBufferSum=0.0_RKIND
         workBufferMin=0.0_RKIND
         workBufferMax=0.0_RKIND
         workBufferSumReduced=0.0_RKIND
         workBufferMinReduced=0.0_RKIND
         workBufferMaxReduced=0.0_RKIND

         ! get pointers to analysis member arrays
         call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceAreaWeightedAveragesAM', &
               surfaceAreaWeightedAveragesAMPool)
         call mpas_pool_get_array(surfaceAreaWeightedAveragesAMPool, 'minValueWithinOceanRegion', minValueWithinOceanRegion)
         call mpas_pool_get_array(surfaceAreaWeightedAveragesAMPool, 'maxValueWithinOceanRegion', maxValueWithinOceanRegion)
         call mpas_pool_get_array(surfaceAreaWeightedAveragesAMPool, 'avgValueWithinOceanRegion', avgValueWithinOceanRegion)

         ! get pointers to scratch variables
         call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceAreaWeightedAveragesAMScratch', scratchPool)
         call mpas_pool_get_field(scratchPool, 'workArray', workArrayField)
         call mpas_pool_get_field(scratchPool, 'workMin', workMinField)
         call mpas_pool_get_field(scratchPool, 'workMax', workMaxField)
         call mpas_pool_get_field(scratchPool, 'workSum', workSumField)
         call mpas_allocate_scratch_field(workArrayField, .true.)
         call mpas_allocate_scratch_field(workMinField, .true.)
         call mpas_allocate_scratch_field(workMaxField, .true.)
         call mpas_allocate_scratch_field(workSumField, .true.)
         workArray => workArrayField % array
         workMin => workMinField % array
         workMax => workMaxField % array
         workSum => workSumField % array

         ! loop over blocks
         ! NOTE: code is not valid for multiple blocks !
         block => domain % blocklist
         do while (associated(block))
            ! get pointers to pools
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)

            ! get pointers to mesh
            call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
            call mpas_pool_get_dimension(block % dimensions, 'nSfcAreaWeightedAvgFields', nSfcAreaWeightedAvgFields)
            call mpas_pool_get_dimension(block % dimensions, 'nOceanRegions', nOceanRegions)
            call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
            call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)
            call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
            call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
            call mpas_pool_get_array(meshPool, 'latCell', latCell)

            ! test to make sure the arrays are big enough
            nDefinedDataFields = size(avgValueWithinOceanRegion,dim=1)
            if (nDefinedDataFields > nSfcAreaWeightedAvgFields) then
                write (stderrUnit,*) 'Abort: nDefinedDataFields > nSfcAreaWeightedAvgFields'
                write (stderrUnit,*) '     :  increase size of ocn_surface_area_weighted_averages scratch space'
                call mpas_dmpar_global_abort('MPAS-ocean: Abort: nDefinedDataFields > nSfcAreaWeightedAvgFields')
            endif

            ! get pointers to data that will be analyzed
            ! listed in the order in which the fields appear in {avg,min,max}SurfaceStatistics
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUp)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', &
                  longWaveHeatFluxDown)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', &
                  seaIceFreshWaterFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
            if (frazilIcePkgActive)   call mpas_pool_get_array(forcingPool, 'seaIceEnergy', seaIceEnergy)
            call mpas_pool_get_array(forcingPool, 'surfaceThicknessFlux', surfaceThicknessFlux)
            call mpas_pool_get_array(tracersSurfaceFluxPool, 'activeTracersSurfaceFlux', activeTracersSurfaceFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFlux)
            call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', surfaceStressMagnitude)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'windStressMeridional', &
                  windStressMeridional)
            call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
            call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
            call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
            call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)

            ! initialize buffers
            workBufferSum(:) = 0.0_RKIND
            workBufferMin(:) = +1.0e20_RKIND
            workBufferMax(:) = -1.0e20_RKIND

            ! copy data into work array
            workArray( :,:) = 0.0_RKIND
            ! workArray( 1,:) used to contain workMask, which is now written out as variable surfaceMask
            ! kept blank to keep old analysis scripts working
            workArray( 2,:) = areaCell(:)
            if (activeTracersBulkRestoringPKG) workArray( 3,:) = latentHeatFlux(:)
            if (activeTracersBulkRestoringPKG) workArray( 4,:) = sensibleHeatFlux(:)
            if (activeTracersBulkRestoringPKG) workArray( 5,:) = longWaveHeatFluxUp(:)
            if (activeTracersBulkRestoringPKG) workArray( 6,:) = longWaveHeatFluxDown(:)
            if (activeTracersBulkRestoringPKG) workArray( 7,:) = seaIceHeatFlux(:)
            if (activeTracersBulkRestoringPKG) workArray( 8,:) = shortWaveHeatFlux(:)
            if (activeTracersBulkRestoringPKG) workArray( 9,:) = evaporationFlux(:)
            if (activeTracersBulkRestoringPKG) workArray(10,:) = seaIceFreshWaterFlux(:)
            if (activeTracersBulkRestoringPKG) workArray(11,:) = riverRunoffFlux(:)
            if (activeTracersBulkRestoringPKG) workArray(12,:) = iceRunoffFlux(:)
            if (activeTracersBulkRestoringPKG) workArray(13,:) = rainFlux(:)
            if (activeTracersBulkRestoringPKG) workArray(14,:) = snowFlux(:)
            if (frazilIcePkgActive) workArray(15,:) = seaIceEnergy(:)
            workArray(16,:) = surfaceThicknessFlux(:)
            workArray(17,:) = activeTracersSurfaceFlux(indexTemperature,:)
            workArray(18,:) = activeTracersSurfaceFlux(indexSalinity,:)
            if (activeTracersBulkRestoringPKG) workArray(19,:) = seaIceSalinityFlux(:)
            workArray(20,:) = surfaceStressMagnitude(:)
            if (activeTracersBulkRestoringPKG) workArray(21,:) = windStressZonal(:)
            if (activeTracersBulkRestoringPKG) workArray(22,:) = windStressMeridional(:)
            workArray(23,:) = atmosphericPressure(:)
            workArray(24,:) = ssh(:)
            workArray(25,:) = activeTracers(indexTemperature,1,:)
            workArray(26,:) = activeTracers(indexSalinity,1,:)
            workArray(27,:) = boundaryLayerDepth(:)

            ! build net heat, salinity and fresh water budget
            ! net heat into ocean = latentHeatFlux + sensibleHeatFlux + longWaveHeatFluxUp + longWaveHeatFluxDown
            !                     + shortWaveHeatFlux + seaIceHeatFlux + (?seaIceEnergy?)
            ! net salinity into ocean = seaIceSalinityFlux
            ! net freshwater into ocean = evaporationFlux + seaIceFreshWaterFlux + riverRunoffFlux + iceRunoffFlux
            !                           + rainFlux + snowFlux + (?seaIceEnergy?)
            if (activeTracersBulkRestoringPKG) then
               workArray(28,:) =   latentHeatFlux(:) &
                                 + sensibleHeatFlux(:) &
                                 + longWaveHeatFluxUp(:) &
                                 + longWaveHeatFluxDown(:) &
                                 + shortWaveHeatFlux(:)  &
                                 + seaIceHeatFlux(:)
                        !        + seaIceEnergy
               workArray(29,:) =   seaIceSalinityFlux(:)
               workArray(30,:) =   evaporationFlux(:) &
                                 + seaIceFreshWaterFlux(:) &
                                 + riverRunoffFlux(:) &
                                 + iceRunoffFlux(:) &
                                 + rainFlux(:) &
                                 + snowFlux(:)
                        !        + seaIceEnergy(:)
            end if

            call compute_statistics(nDefinedDataFields, nCellsSolve, workArray, surfaceMask, workMin, workMax, workSum, &
                  nOceanRegions)

            ! store data in buffer in order to allow only three dmpar calls
            ! loop over all ocean regions
            do iRegion=1,nOceanRegions
               do iDataField=1,nDefinedDataFields
                  kBuffer = kBuffer+1
                  workBufferSum(kBuffer) = workBufferSum(iDataField) + workSum(iDataField, iRegion)
                  workBufferMin(kBuffer) = min(workBufferMin(iDataField), workMin(iDataField, iRegion))
                  workBufferMax(kBuffer) = max(workBufferMax(iDataField), workMax(iDataField, iRegion))
               enddo !iDataField
            enddo !iRegion

            block => block % next
         enddo !block
      endif

      if (additionalRegion /= '') then
         !!! region file
         ! if a region file is given and a region is selected, set up region file version
         ! region file dimensions
         call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegions', nRegions)
         call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegionGroups', nRegionGroups)
         call mpas_pool_get_dimension(domain % blocklist % dimensions, 'maxRegionsInGroup', maxRegionsInGroup)

         ! get region file dimensions
         call mpas_pool_get_subpool(domain % blocklist % structs, 'regions', regionPool)
         call mpas_pool_get_array(regionPool, 'regionsInGroup', regionsInGroup)
         call mpas_pool_get_array(regionPool, 'nRegionsInGroup', nRegionsInGroup)
         call mpas_pool_get_array(regionPool, 'regionNames', regionNames)
         call mpas_pool_get_array(regionPool, 'regionGroupNames', regionGroupNames)
         call mpas_pool_get_array(regionPool, 'regionCellMasks', regionCellMasks)

         regionGroupOffset = 1
         regionGroupNumber = 0

         ! region preliminaries
         ! figure out the region group number that matches the configured additional region's name
         do i=1,nRegionGroups
            if (regionGroupNames(i) .eq. additionalRegion) then
               regionGroupNumber = i
               ! determine offset to compensate for several region groups in the
               ! regions file
               do j=1,i-1
                 regionGroupOffset = regionGroupOffset + nRegionsInGroup(j)
               enddo !j
            endif
         enddo !i

         regionsInAddGroup => nRegionsInGroup(regionGroupNumber)

         ! allocate buffer for message passing
         kBuffer=0
         kBufferLength=regionsInAddGroup*nSfcAreaWeightedAvgFields
         allocate(workBufferSumRegionMask(kBufferLength))
         allocate(workBufferMinRegionMask(kBufferLength))
         allocate(workBufferMaxRegionMask(kBufferLength))
         allocate(workBufferSumReducedRegionMask(kBufferLength))
         allocate(workBufferMinReducedRegionMask(kBufferLength))
         allocate(workBufferMaxReducedRegionMask(kBufferLength))
         workBufferSumRegionMask=0.0_RKIND
         workBufferMinRegionMask=0.0_RKIND
         workBufferMaxRegionMask=0.0_RKIND
         workBufferSumReducedRegionMask=0.0_RKIND
         workBufferMinReducedRegionMask=0.0_RKIND
         workBufferMaxReducedRegionMask=0.0_RKIND

         ! get pointers to analysis member arrays
         call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceAreaWeightedAveragesAM', &
               surfaceAreaWeightedAveragesAMPool)
         call mpas_pool_get_array(surfaceAreaWeightedAveragesAMPool, 'minValueWithinOceanRegionMask', minValueWithinOceanRegionMask)
         call mpas_pool_get_array(surfaceAreaWeightedAveragesAMPool, 'maxValueWithinOceanRegionMask', maxValueWithinOceanRegionMask)
         call mpas_pool_get_array(surfaceAreaWeightedAveragesAMPool, 'avgValueWithinOceanRegionMask', avgValueWithinOceanRegionMask)

         ! get pointers to scratch variables
         call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceAreaWeightedAveragesAMScratch', scratchPool)
         call mpas_pool_get_field(scratchPool, 'workArrayRegionMask', workArrayFieldRegionMask)
         call mpas_pool_get_field(scratchPool, 'workMinRegionMask', workMinFieldRegionMask)
         call mpas_pool_get_field(scratchPool, 'workMaxRegionMask', workMaxFieldRegionMask)
         call mpas_pool_get_field(scratchPool, 'workSumRegionMask', workSumFieldRegionMask)
         call mpas_allocate_scratch_field(workArrayFieldRegionMask, .true.)
         call mpas_allocate_scratch_field(workMinFieldRegionMask, .true.)
         call mpas_allocate_scratch_field(workMaxFieldRegionMask, .true.)
         call mpas_allocate_scratch_field(workSumFieldRegionMask, .true.)
         workArrayRegionMask => workArrayFieldRegionMask % array
         workMinRegionMask => workMinFieldRegionMask % array
         workMaxRegionMask => workMaxFieldRegionMask % array
         workSumRegionMask => workSumFieldRegionMask % array
        ! loop over blocks
        ! NOTE: code is not valid for multiple blocks !
         block => domain % blocklist
         do while (associated(block))
            ! get pointers to pools
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)

            ! get pointers to mesh
            call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
            call mpas_pool_get_dimension(block % dimensions, 'nSfcAreaWeightedAvgFields', nSfcAreaWeightedAvgFields)
            call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
            call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)
            call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
            call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
            call mpas_pool_get_array(meshPool, 'latCell', latCell)

            ! test to make sure the arrays are big enough
            nDefinedDataFields = size(avgValueWithinOceanRegionMask,dim=1)
            if (nDefinedDataFields > nSfcAreaWeightedAvgFields) then
                write (stderrUnit,*) 'Abort: nDefinedDataFields > nSfcAreaWeightedAvgFields'
                write (stderrUnit,*) '     :  increase size of ocn_surface_area_weighted_averages scratch space'
                write (stderrUnit,*) nSfcAreaWeightedAvgFields, nDefinedDataFields, nOceanRegions, maxRegionsInGroup
                call mpas_dmpar_global_abort('MPAS-ocean: Abort: nDefinedDataFields > nSfcAreaWeightedAvgFields')
            endif

            ! get pointers to data that will be analyzed
            ! listed in the order in which the fields appear in {avg,min,max}SurfaceStatistics
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUp)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', &
                  longWaveHeatFluxDown)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', &
                  seaIceFreshWaterFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
            if (frazilIcePkgActive)   call mpas_pool_get_array(forcingPool, 'seaIceEnergy', seaIceEnergy)
            call mpas_pool_get_array(forcingPool, 'surfaceThicknessFlux', surfaceThicknessFlux)
            call mpas_pool_get_array(tracersSurfaceFluxPool, 'activeTracersSurfaceFlux', activeTracersSurfaceFlux)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFlux)
            call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', surfaceStressMagnitude)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
            if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'windStressMeridional', &
                  windStressMeridional)
            call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
            call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
            call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
            call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth',boundaryLayerDepth)

            ! initialize buffers
            workBufferSumRegionMask(:) = 0.0_RKIND
            workBufferMinRegionMask(:) = +1.0e20_RKIND
            workBufferMaxRegionMask(:) = -1.0e20_RKIND

            ! copy data into work array
            workArrayRegionMask( :,:) = 0.0_RKIND
            workArrayRegionMask( 1,:) = 0.0_RKIND ! originally used to contain a mask; using regionCellMasks here
            workArrayRegionMask( 2,:) = areaCell(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask( 3,:) = latentHeatFlux(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask( 4,:) = sensibleHeatFlux(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask( 5,:) = longWaveHeatFluxUp(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask( 6,:) = longWaveHeatFluxDown(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask( 7,:) = seaIceHeatFlux(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask( 8,:) = shortWaveHeatFlux(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask( 9,:) = evaporationFlux(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask(10,:) = seaIceFreshWaterFlux(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask(11,:) = riverRunoffFlux(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask(12,:) = iceRunoffFlux(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask(13,:) = rainFlux(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask(14,:) = snowFlux(:)
            if (frazilIcePkgActive) workArrayRegionMask(15,:) = seaIceEnergy(:)
            workArrayRegionMask(16,:) = surfaceThicknessFlux(:)
            workArrayRegionMask(17,:) = activeTracersSurfaceFlux(indexTemperature,:)
            workArrayRegionMask(18,:) = activeTracersSurfaceFlux(indexSalinity,:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask(19,:) = seaIceSalinityFlux(:)
            workArrayRegionMask(20,:) = surfaceStressMagnitude(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask(21,:) = windStressZonal(:)
            if (activeTracersBulkRestoringPKG) workArrayRegionMask(22,:) = windStressMeridional(:)
            workArrayRegionMask(23,:) = atmosphericPressure(:)
            workArrayRegionMask(24,:) = ssh(:)
            workArrayRegionMask(25,:) = activeTracers(indexTemperature,1,:)
            workArrayRegionMask(26,:) = activeTracers(indexSalinity,1,:)
            workArrayRegionMask(27,:) = boundaryLayerDepth(:)

            ! build net heat, salinity and fresh water budget
            ! net heat into ocean = latentHeatFlux + sensibleHeatFlux + longWaveHeatFluxUp + longWaveHeatFluxDown
            !                     + shortWaveHeatFlux + seaIceHeatFlux + (?seaIceEnergy?)
            ! net salinity into ocean = seaIceSalinityFlux
            ! net freshwater into ocean = evaporationFlux + seaIceFreshWaterFlux + riverRunoffFlux + iceRunoffFlux
            !                           + rainFlux + snowFlux + (?seaIceEnergy?)
            if (activeTracersBulkRestoringPKG) then
               workArrayRegionMask(28,:) = latentHeatFlux(:) &
                                 + sensibleHeatFlux(:) &
                                 + longWaveHeatFluxUp(:) &
                                 + longWaveHeatFluxDown(:) &
                                 + shortWaveHeatFlux(:)  &
                                 + seaIceHeatFlux(:)
                        !        + seaIceEnergy
               workArrayRegionMask(29,:) = seaIceSalinityFlux(:)
               workArrayRegionMask(30,:) = evaporationFlux(:) &
                                 + seaIceFreshWaterFlux(:) &
                                 + riverRunoffFlux(:) &
                                 + iceRunoffFlux(:) &
                                 + rainFlux(:) &
                                 + snowFlux(:)
                        !        + seaIceEnergy(:)
            end if

            call compute_statisticsRegionMask(nDefinedDataFields, nCellsSolve, workArrayRegionMask, regionCellMasks, &
                     workMinRegionMask, workMaxRegionMask, workSumRegionMask, regionsInAddGroup)

            ! loop over all ocean regions
            do iRegion=1,regionsInAddGroup
               do iDataField=1,nDefinedDataFields
                  kBuffer = kBuffer+1
                  workBufferSumRegionMask(kBuffer) = workBufferSumRegionMask(iDataField) + workSumRegionMask(iDataField, iRegion)
                  workBufferMinRegionMask(kBuffer) = min(workBufferMinRegionMask(iDataField), workMinRegionMask(iDataField, iRegion))
                  workBufferMaxRegionMask(kBuffer) = max(workBufferMaxRegionMask(iDataField), workMaxRegionMask(iDataField, iRegion))
               enddo !iDataField

            enddo !iRegion

            block => block % next
         enddo !block
      endif


      if (predefRegions == .true.) then
         !!! hard-wired region averaging
         ! communication
         call mpas_dmpar_sum_real_array(dminfo, kBufferLength, workBufferSum, workBufferSumReduced)
         call mpas_dmpar_min_real_array(dminfo, kBufferLength, workBufferMin, workBufferMinReduced)
         call mpas_dmpar_max_real_array(dminfo, kBufferLength, workBufferMax, workBufferMaxReduced)

         ! unpack the buffer into intent(out) of this analysis member
         kBuffer=0
         do iRegion=1,nOceanRegions
           do iDataField=1,nDefinedDataFields
              kBuffer = kBuffer+1
              avgValueWithinOceanRegion(iDataField,iRegion)=workBufferSumReduced(kBuffer)
              minValueWithinOceanRegion(iDataField,iRegion)=workBufferMinReduced(kBuffer)
              maxValueWithinOceanRegion(iDataField,iRegion)=workBufferMaxReduced(kBuffer)
           enddo
         enddo

         ! normalize averages
         do iRegion=1,nOceanRegions
            ! normalize all field by total area
            do iDataField=3,nDefinedDataFields
               avgValueWithinOceanRegion(iDataField,iRegion) = avgValueWithinOceanRegion(iDataField,iRegion) &
                                                             / max(avgValueWithinOceanRegion(2,iRegion),1.0e-8_RKIND)
            enddo
            ! normalize total area by number of cells in region
            avgValueWithinOceanRegion(2,iRegion) = avgValueWithinOceanRegion(2,iRegion) &
                                                 / max(avgValueWithinOceanRegion(1,iRegion),1.0e-8_RKIND)
         enddo
      endif


      if (additionalRegion /= '') then
         !!! region file averaging
         ! communication
         call mpas_dmpar_sum_real_array(dminfo, kBufferLength, workBufferSumRegionMask, workBufferSumReducedRegionMask)
         call mpas_dmpar_min_real_array(dminfo, kBufferLength, workBufferMinRegionMask, workBufferMinReducedRegionMask)
         call mpas_dmpar_max_real_array(dminfo, kBufferLength, workBufferMaxRegionMask, workBufferMaxReducedRegionMask)

         ! unpack the buffer into intent(out) of this analysis member
         kBuffer=0
         do iRegion=1,regionsInAddGroup
           do iDataField=1,nDefinedDataFields
              kBuffer = kBuffer+1
              avgValueWithinOceanRegionMask(iDataField,iRegion) = workBufferSumReducedRegionMask(kBuffer)
              minValueWithinOceanRegionMask(iDataField,iRegion) = workBufferMinReducedRegionMask(kBuffer)
              maxValueWithinOceanRegionMask(iDataField,iRegion) = workBufferMaxReducedRegionMask(kBuffer)
           enddo
         enddo

         ! normalize averages
         do iRegion=1,regionsInAddGroup
            ! normalize all field by total area
            do iDataField=3,nDefinedDataFields
               avgValueWithinOceanRegionMask(iDataField,iRegion) = avgValueWithinOceanRegionMask(iDataField,iRegion) &
                                                             / max(avgValueWithinOceanRegionMask(2,iRegion),1.0e-8_RKIND)
            enddo
            ! normalize total area by number of cells in region
            avgValueWithinOceanRegionMask(2,iRegion) = avgValueWithinOceanRegionMask(2,iRegion) &
                                                 / max(avgValueWithinOceanRegionMask(1,iRegion),1.0e-8_RKIND)
         enddo   ! communication
      endif

      ! deallocate scratch fields
      call mpas_deallocate_scratch_field(workArrayField, .true.)
      call mpas_deallocate_scratch_field(workMinField, .true.)
      call mpas_deallocate_scratch_field(workMaxField, .true.)
      call mpas_deallocate_scratch_field(workSumField, .true.)
      call mpas_deallocate_scratch_field(workArrayFieldRegionMask, .true.)
      call mpas_deallocate_scratch_field(workMinFieldRegionMask, .true.)
      call mpas_deallocate_scratch_field(workMaxFieldRegionMask, .true.)
      call mpas_deallocate_scratch_field(workSumFieldRegionMask, .true.)

      ! deallocate buffers
      deallocate(workBufferSumReduced)
      deallocate(workBufferMinReduced)
      deallocate(workBufferMaxReduced)
      deallocate(workBufferSumReducedRegionMask)
      deallocate(workBufferMinReducedRegionMask)
      deallocate(workBufferMaxReducedRegionMask)


   contains

   subroutine compute_statistics(nDefinedDataFields, nCellsSolve, workArray, surfaceMask, workMin, workMax, workSum, nRegions)
   ! this subroutines does the actual summing, min, max, masking ect
   ! this hides the messy code from the high-level subroutine

   integer, intent(in) :: nDefinedDataFields, nCellsSolve
   real(kind=RKIND), dimension(:,:), intent(in) :: workArray
   real(kind=RKIND), dimension(:,:), intent(in) :: surfaceMask
   real(kind=RKIND), dimension(:,:), intent(out) :: workMin, workMax, workSum
   integer :: iCell, iData
   integer, pointer, intent(in) :: nRegions
   real(kind=RKIND) :: cellMask

   workSum = 0.0_RKIND
   do iCell=1,nCellsSolve
      do iRegion=1,nRegions
         cellMask = surfaceMask(iRegion,iCell)
         ! don't need workSum(2,iRegion) because mask is written out as variable surfaceMask
         workSum(2,iRegion) = workSum(2,iRegion) + workArray(2,iCell)*cellMask
         do iData=3,nDefinedDataFields
            workSum(iData,iRegion) = workSum(iData,iRegion) + workArray(2,iCell)*workArray(iData,iCell)*cellMask
            workMin(iData,iRegion) = minval(workArray(iData,1:nCellsSolve),surfaceMask(iRegion,1:nCellsSolve)>0.5_RKIND)
            workMax(iData,iRegion) = maxval(workArray(iData,1:nCellsSolve),surfaceMask(iRegion,1:nCellsSolve)>0.5_RKIND)
         enddo
      enddo
   enddo
   end subroutine compute_statistics


   subroutine compute_statisticsRegionMask(nDefinedDataFields, nCellsSolve, workArray, surfaceMask, workMin, workMax, workSum, nRegions)
   ! this subroutines does the actual summing, min, max, masking ect
   ! this hides the messy code from the high-level subroutine

   integer, intent(in) :: nDefinedDataFields, nCellsSolve
   real(kind=RKIND), dimension(:,:), intent(in) :: workArray
   integer, dimension(:,:), intent(in) :: surfaceMask
   real(kind=RKIND), dimension(:,:), intent(out) :: workMin, workMax, workSum
   integer :: iCell, iData
   integer, pointer, intent(in) :: nRegions
   real(kind=RKIND) :: cellMask

   workSum = 0.0_RKIND
   do iCell=1,nCellsSolve
      do iRegion=1,nRegions
         cellMask = surfaceMask(iRegion,iCell)
         ! don't need workSum(2,iRegion) because mask is written out as variable surfaceMask
         workSum(2,iRegion) = workSum(2,iRegion) + workArray(2,iCell)*cellMask
         do iData=3,nDefinedDataFields
            workSum(iData,iRegion) = workSum(iData,iRegion) + workArray(2,iCell)*workArray(iData,iCell)*cellMask
            workMin(iData,iRegion) = minval(workArray(iData,1:nCellsSolve),surfaceMask(iRegion,1:nCellsSolve)>0.5_RKIND)
            workMax(iData,iRegion) = maxval(workArray(iData,1:nCellsSolve),surfaceMask(iRegion,1:nCellsSolve)>0.5_RKIND)
         enddo
      enddo
   enddo
   end subroutine compute_statisticsRegionMask

   end subroutine ocn_compute_surface_area_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_restart_surface_area_weighted_averages
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_surface_area_weighted_averages(domain, err)!{{{

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

   end subroutine ocn_restart_surface_area_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_finalize_surface_area_weighted_averages
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_surface_area_weighted_averages(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

   end subroutine ocn_finalize_surface_area_weighted_averages!}}}

end module ocn_surface_area_weighted_averages

! vim: foldmethod=marker
