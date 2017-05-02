! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_layer_volume_weighted_averages
!
!> \brief MPAS ocean analysis member: horizonal layer volume weighted averages at each vertical level
!> \author Todd Ringler
!> \date   April 24, 2015
!> \details
!>  MPAS ocean analysis member: layer_volume_weighted_averages
!
!-----------------------------------------------------------------------

module ocn_layer_volume_weighted_averages

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

   public :: ocn_init_layer_volume_weighted_averages, &
             ocn_compute_layer_volume_weighted_averages, &
             ocn_restart_layer_volume_weighted_averages, &
             ocn_finalize_layer_volume_weighted_averages

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   ! create a multi-block region mask temporarily
   ! this should replaced later by a MPAS ocean-wide region definition

   real (kind=RKIND), dimension(:,:,:), allocatable :: &
      regionMask

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_layer_volume_weighted_averages
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_layer_volume_weighted_averages(domain, err)!{{{

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
      type (mpas_pool_type), pointer :: meshPool

      integer, pointer :: &
         nCells,          &! number of cells in a given block
         nCellsSolve,     &! number of cells owned by a given block
         nRegions          ! number of ocean regions defined

      real (kind=RKIND), dimension(:), pointer ::  &
         lonCell, latCell  ! latitude, longitude of each cell

      integer ::          &
         iBlock,          &! block counter
         numBlocks,       &! track number of blocks and
         maxCells          ! largest cell count to size array

      !-----------------------------------------------------------------

      err = 0

      ! set highest level pointer
      dminfo = domain % dminfo

      ! retrieve number of regions
      call mpas_pool_get_dimension(domain % blocklist % dimensions, &
                                   'nOceanRegionsTmp', nRegions)

      ! count up the number of blocks and the max size of each domain
      numBlocks = 0
      maxCells  = 0
      block => domain % blocklist
      do while (associated(block))

         numBlocks = numBlocks + 1
         call mpas_pool_get_dimension(block % dimensions, 'nCells', &
                                                           nCells)
         maxCells = max(maxCells, nCells)

         block => block % next
      end do

      ! allocate the region mask array
      allocate(regionMask(maxCells,nRegions,numBlocks))

      ! compute the region mask for each block
      iBlock = 0
      block => domain % blocklist
      do while (associated(block))

         iBlock = iBlock + 1

         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', &
                                                           nCellsSolve)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
         call mpas_pool_get_array(meshPool, 'latCell', latCell)

         call ocn_LayerVolWgtAvgRegionMask(iBlock, nCellsSolve, &
                                           lonCell, latCell)

         block => block % next
      end do

      !-----------------------------------------------------------------

   end subroutine ocn_init_layer_volume_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_compute_layer_volume_weighted_averages
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_layer_volume_weighted_averages(domain, timeLevel, err)!{{{

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
      type (mpas_pool_type), pointer :: layerVolumeWeightedAverageAMPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: tracersPool

      ! pointers to result arrays
      real (kind=RKIND), dimension(:,:,:), pointer ::    &
                         avgLayer, minLayer, maxLayer
      real (kind=RKIND), dimension(:,:),   pointer ::     &
                         avgField, minField, maxField

      ! pointers to data in pools to be analyzed
      real (kind=RKIND), dimension(:,:),   pointer :: layerThickness
      real (kind=RKIND), dimension(:,:),   pointer :: density
      real (kind=RKIND), dimension(:,:),   pointer :: potentialDensity
      real (kind=RKIND), dimension(:,:),   pointer :: BruntVaisalaFreqTop
      real (kind=RKIND), dimension(:,:),   pointer :: velocityZonal
      real (kind=RKIND), dimension(:,:),   pointer :: velocityMeridional
      real (kind=RKIND), dimension(:,:),   pointer :: vertVelocityTop
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      real (kind=RKIND), dimension(:,:),   pointer :: kineticEnergyCell
      real (kind=RKIND), dimension(:,:),   pointer :: relativeVorticityCell
      real (kind=RKIND), dimension(:,:),   pointer :: divergence

      ! pointers to data in mesh pool
      integer, pointer :: nVertLevels, nCells, nCellsSolve, &
                          nFields, nOceanRegionsTmp
      integer, pointer :: index_temperature, index_salinity
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer ::  areaCell

      ! local variables
      integer :: iCell, iLevel, iRegion, iTracer, iField, iBlock

      real (kind=RKIND) :: normFactor ! for normalizing averages

      real (kind=RKIND), dimension(:,:), allocatable :: &
         vertMask,                     &!mask for active vertical levels
         fieldTmp                       !field array in local order

      ! buffers data for message passaging
      integer :: bufferAddr, bufferLength
      real (kind=RKIND), dimension(:), allocatable ::  &
                 msgBuffer, msgBufferReduced

      ! end of preamble
      !-----------------------------------------------------------------
      ! begin code

      ! assume no error
      err = 0

      ! set highest level pointer
      dminfo = domain % dminfo

      ! find the number of regions, fields and vertical levels
      call mpas_pool_get_dimension(domain % blocklist % dimensions, &
                                   'nOceanRegionsTmp', nOceanRegionsTmp)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, &
                                   'nLayerVolWeightedAvgFields', nFields)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, &
                                   'nVertLevels', nVertLevels)

      ! get pointers to analysis member arrays to store results
      call mpas_pool_get_subpool(domain % blocklist % structs, &
                               'layerVolumeWeightedAverageAM', &
                                layerVolumeWeightedAverageAMPool)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool, &
                               'minValueWithinOceanLayerRegion', &
                                minLayer)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool, &
                               'maxValueWithinOceanLayerRegion', &
                                maxLayer)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool, &
                               'avgValueWithinOceanLayerRegion', &
                                avgLayer)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool,  &
                               'minValueWithinOceanVolumeRegion', &
                                minField)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool,  &
                               'maxValueWithinOceanVolumeRegion', &
                                maxField)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool,  &
                               'avgValueWithinOceanVolumeRegion', &
                                avgField)

      ! initialize so this routine works for multiple blocks
      avgLayer =  0.0_RKIND
      minLayer =  1.e20_RKIND
      maxLayer = -1.e20_RKIND

      ! loop over blocks
      iBlock = 0
      block => domain % blocklist
      do while (associated(block))

         iBlock = iBlock + 1

         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state',       &
                                                      statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh',        &
                                                      meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', &
                                                      diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing',     &
                                                      forcingPool)
         call mpas_pool_get_subpool(statePool,       'tracers',     &
                                                      tracersPool)

         ! get pointers to mesh
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', &
                                                           nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nCells',      &
                                                           nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', &
                                                           nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nOceanRegionsTmp', &
                                                           nOceanRegionsTmp)
         call mpas_pool_get_dimension(tracersPool, 'index_temperature', &
                                                    index_temperature)
         call mpas_pool_get_dimension(tracersPool, 'index_salinity',    &
                                                    index_salinity)
         call mpas_pool_get_array(meshPool, 'areaCell',     areaCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         ! get pointers to data that will be analyzed
         ! listed in the order in which the fields appear in output
         call mpas_pool_get_array(statePool,       'layerThickness', &
                                                    layerThickness, 1)
         call mpas_pool_get_array(diagnosticsPool, 'density', &
                                                    density)
         call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', &
                                                    potentialDensity)
         call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', &
                                                    BruntVaisalaFreqTop)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', &
                                                    velocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', &
                                                    velocityMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', &
                                                    vertVelocityTop)
         call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', &
                                                    kineticEnergyCell)
         call mpas_pool_get_array(diagnosticsPool, 'relativeVorticityCell', &
                                                    relativeVorticityCell)
         call mpas_pool_get_array(diagnosticsPool, 'divergence', &
                                                    divergence)
         call mpas_pool_get_array(tracersPool,     'activeTracers', &
                                                    activeTracers, 1)


         ! allocate arrays to hold local block masks, and fields
         !$omp single
         allocate(fieldTmp(nVertLevels,nCellsSolve))
         allocate(vertMask(nVertLevels,nCells))
         !$omp end single

         call ocn_LayerVolWgtAvgVertMask(maxLevelCell, nCellsSolve, &
                                         vertMask)

         iField = 1 ! mask
         call ocn_LayerVolWgtAvgFieldMinMaxSum(                         &
          avgLayer, minLayer, maxLayer, vertMask, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 2 ! cell area
         !$omp do collapse(2) private(iCell,iLevel)
         do iCell=1,nCellsSolve
         do iLevel=1,nVertLevels
            fieldTmp(iLevel,iCell) = areaCell(iCell)
         end do
         end do
         !$omp end do
         call ocn_LayerVolWgtAvgFieldMinMaxSum(                         &
          avgLayer, minLayer, maxLayer, fieldTmp, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         ! for field 3, min, max refers to layer thickness, but the
         ! average is for cell volume, so include area in mask
         ! (will be doing something similar for remaining fields)
         iField = 3 ! cell volume
         !$omp do private(iCell,iLevel)
         do iCell=1,nCellsSolve
         do iLevel=1,nVertLevels
            vertMask(iLevel,iCell) = vertMask(iLevel,iCell)*areaCell(iCell)
         end do
         end do
         !$omp end do
         call ocn_LayerVolWgtAvgFieldMinMaxSum(                               &
          avgLayer, minLayer, maxLayer, layerThickness, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         ! for all remaining fields, we want volume-weighted sums
         ! add layer thickness to the vertMask to get a volume mask
         ! use the transposed version of layerThick in fieldTmp
         !$omp do private(iCell,iLevel)
         do iCell=1,nCellsSolve
         do iLevel=1,nVertLevels
            vertMask(iLevel,iCell) = vertMask(iLevel,iCell)* &
                               layerThickness(iLevel,iCell)
         end do
         end do
         !$omp end do

         iField = 4 ! density
         call ocn_LayerVolWgtAvgFieldMinMaxSum(                        &
          avgLayer, minLayer, maxLayer, density, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 5 ! potential density
         call ocn_LayerVolWgtAvgFieldMinMaxSum( avgLayer, minLayer, &
                  maxLayer, potentialDensity, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 6 ! Brunt Vaisala Freq
         call ocn_LayerVolWgtAvgFieldMinMaxSum( avgLayer, minLayer, &
               maxLayer, BruntVaisalaFreqTop, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 7 ! zonal velocity
         call ocn_LayerVolWgtAvgFieldMinMaxSum(avgLayer, minLayer, &
                    maxLayer, velocityZonal, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 8 ! meridional velocity
         call ocn_LayerVolWgtAvgFieldMinMaxSum(avgLayer, minLayer, &
               maxLayer, velocityMeridional, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 9 ! vertical velocity
         call ocn_LayerVolWgtAvgFieldMinMaxSum(avgLayer, minLayer, &
                  maxLayer, vertVelocityTop, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 10 ! temperature
         !$omp do collapse(2) private(iCell,iLevel)
         do iCell=1,nCellsSolve
         do iLevel=1,nVertLevels
            fieldTmp(iLevel,iCell) = &
                 activeTracers(index_temperature,iLevel,iCell)
         end do
         end do
         !$omp end do
         call ocn_LayerVolWgtAvgFieldMinMaxSum(                         &
          avgLayer, minLayer, maxLayer, fieldTmp, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 11 ! salinity
         !$omp do collapse(2) private(iCell,iLevel)
         do iCell=1,nCellsSolve
         do iLevel=1,nVertLevels
            fieldTmp(iLevel,iCell) = &
                 activeTracers(index_salinity,iLevel,iCell)
         end do
         end do
         !$omp end do
         call ocn_LayerVolWgtAvgFieldMinMaxSum(                         &
          avgLayer, minLayer, maxLayer, fieldTmp, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 12 ! kinetic energy
         call ocn_LayerVolWgtAvgFieldMinMaxSum(avgLayer, minLayer, &
                maxLayer, kineticEnergyCell, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 13 ! relative vorticity
         call ocn_LayerVolWgtAvgFieldMinMaxSum(avgLayer, minLayer, &
            maxLayer, relativeVorticityCell, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 14 ! divergence
         call ocn_LayerVolWgtAvgFieldMinMaxSum(                           &
          avgLayer, minLayer, maxLayer, divergence, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         iField = 15 ! enstrophy
         !$omp do private(iCell,iLevel)
         do iCell=1,nCellsSolve
         do iLevel=1,nVertLevels
            fieldTmp(iLevel,iCell) = relativeVorticityCell(iLevel,iCell)* & 
                                     relativeVorticityCell(iLevel,iCell)
         end do
         end do
         !$omp end do
         call ocn_LayerVolWgtAvgFieldMinMaxSum(                         &
          avgLayer, minLayer, maxLayer, fieldTmp, vertMask, regionMask, &
          nCellsSolve, nVertLevels, nOceanRegionsTmp, iField, iBlock)

         ! deallocate fields from this block
         !$omp single
         deallocate(fieldTmp)
         deallocate(vertMask)
         !$omp end single
         block => block % next
      end do

      ! allocate buffer for message passing
      bufferAddr=0
      bufferLength=nOceanRegionsTmp*nFields*nVertLevels
      !$omp single
      allocate(msgBuffer       (bufferLength))
      allocate(msgBufferReduced(bufferLength))
      !$omp end single

      ! communication for sums
      ! pack the buffer 
      !$omp do collapse(3) private(iRegion,iLevel,iField,bufferAddr)
      do iRegion=1,nOceanRegionsTmp
      do iLevel=1,nVertLevels
      do iField=1,nFields
         bufferAddr = (iRegion -1)*nFields*nVertLevels + &
                      (iLevel  -1)*nFields + iField
         msgBuffer(bufferAddr) = avgLayer(iField,iLevel,iRegion)
      enddo
      enddo
      enddo
      !$omp end do
      call mpas_dmpar_sum_real_array(dminfo, bufferLength, msgBuffer, &
                                                           msgBufferReduced )

      ! compute normalized averages

      ! unpack the buffer 
      !$omp do collapse(3) private(iRegion,iLevel,iField,bufferAddr)
      do iRegion=1,nOceanRegionsTmp
      do iLevel=1,nVertLevels
      do iField=1,nFields
         bufferAddr = (iRegion -1)*nFields*nVertLevels + &
                      (iLevel  -1)*nFields + iField
         avgLayer(iField,iLevel,iRegion)=msgBufferReduced(bufferAddr)
      enddo
      enddo
      enddo
      !$omp end do

      ! compute vertical sum before layer-by-layer normalization
      !$omp do private(iLevel,iField,iRegion,normFactor)
      do iRegion=1,nOceanRegionsTmp
        do iField=1,nFields
          avgField(iField,iRegion) = 0.0_RKIND
          do iLevel=1,nVertLevels
             avgField(iField,iRegion) = avgField(iField,iRegion) &
                                      + avgLayer(iField,iLevel,iRegion)
          enddo
        enddo
        normFactor = 1.0_RKIND / max(avgField(3,iRegion),1.0e-8_RKIND)
        do iField=4,nFields
          avgField(iField,iRegion) = avgField(iField,iRegion)*normFactor
        enddo
        ! normalize total region volume by total volume cell area
        avgField(3,iRegion) = avgField(3,iRegion) &
                            / max(avgField(2,iRegion),1.0e-8_RKIND)
        ! normalize total volume cell area by total number of cells
        avgField(2,iRegion) = avgField(2,iRegion) &
                            / max(avgField(1,iRegion),1.0e-8_RKIND)
      enddo
      !$omp end do

      ! normalize averages layer-by-layer
      !$omp do collapse(2) private(iLevel, iRegion, iField, normFactor)
      do iRegion=1,nOceanRegionsTmp
      do iLevel=1,nVertLevels
         ! normalize all field by total volume in each layer
         normFactor = 1.0_RKIND/max(avgLayer(3,iLevel,iRegion),1.0e-8_RKIND)
         do iField=4,nFields
            avgLayer(iField,iLevel,iRegion) = &
            avgLayer(iField,iLevel,iRegion) * normFactor
         enddo
         ! normalize total layer volume by layer area
         avgLayer(3,iLevel,iRegion) = avgLayer(3,iLevel,iRegion) &
                                / max(avgLayer(2,iLevel,iRegion),1.0e-8_RKIND)
         ! normalize total layer area by number of cells in region
         avgLayer(2,iLevel,iRegion) = avgLayer(2,iLevel,iRegion) &
                                / max(avgLayer(1,iLevel,iRegion),1.0e-8_RKIND)
      enddo
      enddo
      !$omp end do

      ! now compute min/max for each field, region
      ! find global min/max

      ! pack the buffer for minval 
      !$omp do collapse(3) private(iRegion,iLevel,iField,bufferAddr)
      do iRegion=1,nOceanRegionsTmp
      do iLevel=1,nVertLevels
      do iField=1,nFields
         bufferAddr = (iRegion -1)*nFields*nVertLevels + &
                      (iLevel  -1)*nFields + iField
         msgBuffer(bufferAddr) = minLayer(iField,iLevel,iRegion)
      enddo
      enddo
      enddo
      !$omp end do
      ! communicate for global min
      call mpas_dmpar_min_real_array(dminfo, bufferLength, msgBuffer, &
                                                           msgBufferReduced )
      ! unpack the buffer
      !$omp do collapse(3) private(iRegion,iLevel,iField,bufferAddr)
      do iRegion=1,nOceanRegionsTmp
      do iLevel=1,nVertLevels
      do iField=1,nFields
         bufferAddr = (iRegion -1)*nFields*nVertLevels + &
                      (iLevel  -1)*nFields + iField
         minLayer(iField,iLevel,iRegion)=msgBufferReduced(bufferAddr)
      enddo
      enddo
      enddo
      !$omp end do

      ! pack the buffer for maxval 
      !$omp do collapse(3) private(iRegion,iLevel,iField,bufferAddr)
      do iRegion=1,nOceanRegionsTmp
      do iLevel=1,nVertLevels
      do iField=1,nFields
         bufferAddr = (iRegion -1)*nFields*nVertLevels + &
                      (iLevel  -1)*nFields + iField
         msgBuffer(bufferAddr) = maxLayer(iField,iLevel,iRegion)
      enddo
      enddo
      enddo
      !$omp end do
      call mpas_dmpar_max_real_array(dminfo, bufferLength, msgBuffer, &
                                                           msgBufferReduced )
      ! unpack the buffer
      !$omp do collapse(3) private(iRegion,iLevel,iField,bufferAddr)
      do iRegion=1,nOceanRegionsTmp
      do iLevel=1,nVertLevels
      do iField=1,nFields
         bufferAddr = (iRegion -1)*nFields*nVertLevels + &
                      (iLevel  -1)*nFields + iField
         maxLayer(iField,iLevel,iRegion)=msgBufferReduced(bufferAddr)
      enddo
      enddo
      enddo
      !$omp end do

      ! deallocate buffer
      !$omp single
      deallocate(msgBuffer)
      deallocate(msgBufferReduced)
      !$omp end single

      ! find min/max with region volume
      !$omp do collapse(2) private(iRegion, iField)
      do iRegion=1,nOceanRegionsTmp
      do iField=1,nFields
         minField(iField,iRegion) = minval(minLayer(iField,:,iRegion))
         maxField(iField,iRegion) = maxval(maxLayer(iField,:,iRegion))
      enddo
      enddo
      !$omp end do

   end subroutine ocn_compute_layer_volume_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_restart_layer_volume_weighted_averages
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_layer_volume_weighted_averages(domain, err)!{{{

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

   end subroutine ocn_restart_layer_volume_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_finalize_layer_volume_weighted_averages
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_layer_volume_weighted_averages(domain, err)!{{{

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

   end subroutine ocn_finalize_layer_volume_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_LayerVolWgtAvgFieldMinMaxSum
!
!> \brief   Computes min, max and sum of an input field
!> \author  Phil Jones
!> \date    April 19, 2017
!> \details
!>  This is a local routine for this analysis member that computes
!>  the min, max, and sum of a given field under the influence of a 
!>  mask that defines active layers and regions. The results are packed
!>  into a buffer for a later reduction call. The calling routine
!>  supplies the field, a precomputed region mask and a base address 
!>  for the buffer. 
!
!-----------------------------------------------------------------------

   subroutine ocn_LayerVolWgtAvgFieldMinMaxSum(                   &
          localSum, localMin, localMax, field, vertMask, regMask, &
          nCellsSolve, numLevels, numRegions, iField, iBlock)!{{{

   !--------------------------------------------------------------------
   !
   ! input variables 
   !
   !--------------------------------------------------------------------

   integer, intent(in) :: &
      numLevels,          &! number of vertical levels
      nCellsSolve,        &! number of cells owned by this field
      numRegions,         &! number of regions for diagnostics
      iField,             &! index of field in local arrays
      iBlock               ! index of current block for region mask

   real (kind=RKIND), dimension(:,:), intent(in) :: &
      vertMask,           &! multiplicative mask defining active levels
      field                ! field for which statistics desired in
                           ! transposed order (nCells, nVertLevels)

   real (kind=RKIND), dimension(:,:,:), intent(in) :: &
      regMask              ! mask defining active regions

   !--------------------------------------------------------------------
   !
   ! input/output variables 
   !
   !--------------------------------------------------------------------

   ! results are accumulated to account for multiple blocks
   ! and are assumed to be initialized elsewhere

   real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
      localSum, localMin, localMax  ! accumulated local results

   !--------------------------------------------------------------------
   !
   ! local variables 
   !
   !--------------------------------------------------------------------

   integer :: iCell, k, iRegion, iBuff   ! loop indices and addresses

   real (kind=RKIND) :: maskTmp          ! temp for holding mask info

   ! end preamble
   !--------------------------------------------------------------------
   ! begin routine

   ! compute min, max, sum over a masked horizontal region

   !$omp do collapse(2) private(iCell,k,iRegion,maskTmp)
   do iRegion=1,numRegions 
   do k=1,numLevels
      
      ! compute sum
      do iCell=1,nCellsSolve
         localSum(iField,k,iRegion) = &
         localSum(iField,k,iRegion) + field(k,iCell)* &
                                   vertMask(k,iCell)* &
                                    regMask(iCell,iRegion,iBlock)
      end do

      ! compute min and max
      do iCell=1,nCellsSolve
         maskTmp = vertMask(k,iCell)* &
                    regMask(iCell,iRegion,iBlock)
         if (maskTmp /= 0.0_RKIND) then
            localMin(iField,k,iRegion) = min( &
            localMin(iField,k,iRegion),field(k,iCell))
            localMax(iField,k,iRegion) = max( &
            localMax(iField,k,iRegion),field(k,iCell))
         endif
      end do

   end do
   end do
   !$omp end do

   !--------------------------------------------------------------------

   end subroutine ocn_LayerVolWgtAvgFieldMinMaxSum !}}}

!***********************************************************************
!
!  routine ocn_LayerVolWgtAvgRegionMask
!
!> \brief   Computes mask to define local regions
!> \author  Todd Ringler, rev. Phil Jones
!> \date    April 19, 2017
!> \details
!>  This is a local routine for this analysis member that computes
!>  a multiplicative mask that is 1 if the point is a member of the
!>  desired region and zero otherwise. 
!>  NOTE: This routine should be temporary and should be replaced
!>        once an MPAS-wide ocean region and mask capability 
!>        are introduced.
!
!-----------------------------------------------------------------------

   subroutine ocn_LayerVolWgtAvgRegionMask(iBlock, nCellsSolve, &
                                           lonCell, latCell)!{{{

     !------------------------------------------------------------------
     !
     ! input variables
     !
     !------------------------------------------------------------------

     integer, intent(in) :: &
        iBlock,             &! current block in block decomposition
        nCellsSolve          ! number of cells owned by this block

     real(kind=RKIND), dimension(:), intent(in) :: &
        lonCell, latCell     ! lat, lon coordinates for each cell
     
     !------------------------------------------------------------------
     !
     ! output variables
     !
     !------------------------------------------------------------------

     ! output in in module variable regionMask

     !------------------------------------------------------------------
     !
     ! local variables
     !
     !------------------------------------------------------------------

     integer ::                       &
        iCell, iRegion,               &! loop indices
        nCells, nRegions, nCellsMax    ! array sizes

     real(kind=RKIND) ::              &
        dtr                            ! degree to radian conversion

     !------------------------------------------------------------------

     nCellsMax = size(regionMask, dim=1)
     nRegions  = size(regionMask, dim=2)

     ! Currently assumes 6 regions, so fail if nRegions not large enough
     if (nRegions < 6) then
        !$omp master
        write (stderrUnit,*) 'Abort: nRegions inconsistent '
        write (stderrUnit,*) '     : in layer volume weighted avg AM'
        !$omp end master
        call mpas_dmpar_global_abort(&
          'MPAS-ocean: Abort: nRegions inconsistent in layer vol wgt avg AM')
     endif

     ! degrees to radians
     dtr = 4.0_RKIND*atan(1.0_RKIND) / 180.0_RKIND

     ! Define masks for each region
     !$omp do 
     do iCell=1,nCellsSolve

        ! Define Arctic region (region 1)
        if (latCell(iCell) < 60.0_RKIND*dtr) then
           regionMask(iCell,1,iBlock) = 0.0_RKIND
        else
           regionMask(iCell,1,iBlock) = 1.0_RKIND
        endif

        ! Equatorial region (region 2)
        if (latCell(iCell) >  15.0_RKIND*dtr .or. &
            latCell(iCell) < -15.0_RKIND*dtr) then
           regionMask(iCell,2,iBlock) = 0.0_RKIND
        else
           regionMask(iCell,2,iBlock) = 1.0_RKIND
        endif

        ! Southern Ocean (region 3)
        if (latCell(iCell) > -50.0_RKIND*dtr) then
           regionMask(iCell,3,iBlock) = 0.0_RKIND
        else
           regionMask(iCell,3,iBlock) = 1.0_RKIND
        endif

        ! Nino 3 (region 4)
        if (latCell(iCell) >   5.0_RKIND*dtr .or. &
            latCell(iCell) <  -5.0_RKIND*dtr .or. &
            lonCell(iCell) < 210.0_RKIND*dtr .or. &
            lonCell(iCell) > 270.0_RKIND*dtr) then
           regionMask(iCell,4,iBlock) = 0.0_RKIND
        else
           regionMask(iCell,4,iBlock) = 1.0_RKIND
        endif

        ! Nino 4 (region 5)
        if (latCell(iCell) >   5.0_RKIND*dtr .or. &
            latCell(iCell) <  -5.0_RKIND*dtr .or. &
            lonCell(iCell) < 160.0_RKIND*dtr .or. &
            lonCell(iCell) > 210.0_RKIND*dtr) then
           regionMask(iCell,5,iBlock) = 0.0_RKIND
        else
           regionMask(iCell,5,iBlock) = 1.0_RKIND
        endif

        ! Nino 3.4 (region 6)
        if (latCell(iCell) >   5.0_RKIND*dtr .or. &
            latCell(iCell) <  -5.0_RKIND*dtr .or. &
            lonCell(iCell) < 190.0_RKIND*dtr .or. &
            lonCell(iCell) > 240.0_RKIND*dtr) then
           regionMask(iCell,6,iBlock) = 0.0_RKIND
        else
           regionMask(iCell,6,iBlock) = 1.0_RKIND
        endif

     enddo
     !$omp end do

     ! since number of cells varies per block, take care of any
     ! extra cells in this block

     !$omp do collapse(2)
     do iRegion=1,nRegions
     do iCell=nCellsSolve+1,nCellsMax
        regionMask(iCell,iRegion,iBlock) = 0.0_RKIND
     end do
     end do
     !$omp end do

   !--------------------------------------------------------------------

   end subroutine ocn_LayerVolWgtAvgRegionMask !}}}

!***********************************************************************
!
!  routine ocn_LayerVolWgtAvgVertMask
!
!> \brief   Computes mask to define active layers
!> \author  Todd Ringler, rev. Phil Jones
!> \date    April 19, 2017
!> \details
!>  This is a local routine for this analysis member that computes
!>  a multiplicative mask that is 1 if a vertical layer is active
!>  and is zero otherwise. 
!>  NOTE: This routine should be temporary. An MPAS ocean-wide
!>        vertical mask will be used eventually.
!
!-----------------------------------------------------------------------

   subroutine ocn_LayerVolWgtAvgVertMask(maxLevelCell, nCellsSolve, &
                                         vertMask)!{{{

     !------------------------------------------------------------------
     !
     ! input variables
     !
     !------------------------------------------------------------------

     integer, intent(in) :: &
        nCellsSolve          ! number of cells owned by this subdomain 

     integer, intent(in), dimension(:) :: &
        maxLevelCell         ! last active cell in each column

     !------------------------------------------------------------------
     !
     ! output variables
     !
     !------------------------------------------------------------------

     real(kind=RKIND), dimension(:,:), intent(out) :: &
        vertMask             ! multiplicative mask for active points

     !------------------------------------------------------------------
     !
     ! local variables
     !
     !------------------------------------------------------------------

     integer :: &
        iLevel, iCell, &! loop indices
        maxLevel,      &! last active level in each column
        nVertLevels,   &! total number of vertical levels
        nCells          ! total number of cells

     !------------------------------------------------------------------

     nVertLevels = size(vertMask,dim=1)

     ! compute mask assuming all levels are active down to maxLevel

     !$omp do private(iLevel,iCell,maxLevel)
     do iCell=1,nCellsSolve
        maxLevel = maxLevelCell(iCell) 
        do iLevel=1,maxLevel
           vertMask(iLevel,iCell) = 1.0_RKIND
        end do 
        do iLevel=maxLevel+1,nVertLevels
           vertMask(iLevel,iCell) = 0.0_RKIND
        end do 
     end do
     !$omp end do

   !--------------------------------------------------------------------

   end subroutine ocn_LayerVolWgtAvgVertMask !}}}

end module ocn_layer_volume_weighted_averages

! vim: foldmethod=marker
