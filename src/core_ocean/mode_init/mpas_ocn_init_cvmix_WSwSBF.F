! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_cvmix_WSwSBF
!
!> \brief MPAS ocean initialize case -- CVMix Unit Test
!> WSwSBF means Wind Stress with Surface Buoyancy Forcing
!> \author Todd Ringler
!> \date   04/23/2015
!> \details
!>  This module contains the routines for initializing the
!>  the cvmix WSwSBF unit test configuration. This in a
!>  single column configuration
!
!-----------------------------------------------------------------------

module ocn_init_cvmix_WSwSBF

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_stream_manager
   use mpas_dmpar

   use ocn_init_cell_markers
   use ocn_init_vertical_grids

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

   public :: ocn_init_setup_cvmix_WSwSBF, &
             ocn_init_validate_cvmix_WSwSBF

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_cvmix_WSwSBF
!
!> \brief   Setup for cvmix WSwSBF unit test configuration
!> \author  Todd Ringler
!> \date    04/23/2015
!> \details
!>  This routine sets up the initial conditions for the cvmix WSwSBF unit test configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_cvmix_WSwSBF(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr
      real (kind=RKIND) :: temperature, salinity

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool, statePool
      type (mpas_pool_type), pointer :: diagnosticsPool, forcingPool

      type (mpas_pool_type), pointer :: tracersPool, &
                                        tracersSurfaceRestoringFieldsPool, &
                                        tracersInteriorRestoringFieldsPool

      integer,parameter :: seed = 2017
      integer :: nCellsSolve, spot
      integer, pointer :: nVertLevels, nVertLevelsP1, nEdgesSolve, nVerticesSolve
      integer, pointer :: index_temperature, index_salinity, index_tracer1

      integer, dimension(:), pointer :: cullCell,maxLevelCell,nCellsArray
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, refZMid, vertCoordMovementWeights
      real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional
      real (kind=RKIND), dimension(:), pointer :: latentHeatFlux, sensibleHeatFlux, shortWaveHeatFlux
      real (kind=RKIND), dimension(:), pointer :: evaporationFlux, rainFlux
      real (kind=RKIND), dimension(:), pointer :: salinityRestore, bottomDepth, angleEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, yEdge, yVertex, fCell, fEdge, fVertex
      real (kind=RKIND), dimension(:, :), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:, :, :), pointer :: activeTracers, debugTracers
      real (kind=RKIND), dimension(:, :), pointer ::    activeTracersPistonVelocity, activeTracersSurfaceRestoringValue
      real (kind=RKIND), dimension(:, :, :), pointer :: activeTracersInteriorRestoringValue, activeTracersInteriorRestoringRate

      real (kind=RKIND) :: randN, xMin, xMax, yMin, yMax, yMinGlobal, yMaxGlobal, dcEdgeMin, dcEdgeMinGlobal

      real (kind=RKIND), dimension(:), pointer :: interfaceLocations
      integer :: directionProperty
      character (len=StrKIND) :: streamID
      integer :: iCell, iEdge, iVertex, k, kML

      real (kind=RKIND) :: BLdepth

      character (len=StrKIND), pointer :: config_init_configuration, &
                                          config_cvmix_WSwSBF_vertical_grid

      integer, pointer ::                 config_cvmix_WSwSBF_vert_levels

      logical, pointer :: config_cvmix_WSwSBF_wind_zonal_variation

      logical, pointer :: config_cvmix_WSwSBF_useSeaMount, config_cvmix_WSwSBF_bottom_noise

      logical, pointer :: config_cvmix_WSwSBF_useIsland, config_write_cull_cell_mask

      real(kind=RKIND), pointer :: y_period

      real (kind=RKIND), pointer ::       config_cvmix_WSwSBF_surface_temperature, &
                                          config_cvmix_WSwSBF_surface_salinity, &
                                          config_cvmix_WSwSBF_surface_restoring_temperature, &
                                          config_cvmix_WSwSBF_surface_restoring_salinity, &
                                          config_cvmix_WSwSBF_temperature_piston_velocity, &
                                          config_cvmix_WSwSBF_salinity_piston_velocity, &
                                          config_cvmix_WSwSBF_sensible_heat_flux, &
                                          config_cvmix_WSwSBF_latent_heat_flux, &
                                          config_cvmix_WSwSBF_shortwave_heat_flux, &
                                          config_cvmix_WSwSBF_rain_flux, &
                                          config_cvmix_WSwSBF_evaporation_flux, &
                                          config_cvmix_WSwSBF_interior_temperature_restoring_rate, &
                                          config_cvmix_WSwSBF_interior_salinity_restoring_rate, &
                                          config_cvmix_WSwSBF_temperature_gradient, &
                                          config_cvmix_WSwSBF_salinity_gradient, &
                                          config_cvmix_WSwSBF_bottom_depth, &
                                          config_cvmix_WSwSBF_max_windstress_zonal, &
                                          config_cvmix_WSwSBF_coriolis_parameter, &
                                          config_cvmix_WSwSBF_temperature_gradient_mixed_layer, &
                                          config_cvmix_WSwSBF_salinity_gradient_mixed_layer, &
                                          config_cvmix_WSwSBF_mixed_layer_depth_temperature, &
                                          config_cvmix_WSwSBF_mixed_layer_depth_salinity,    &
                                          config_cvmix_WSwSBF_mixed_layer_temperature_change,       &
                                          config_cvmix_WSwSBF_mixed_layer_salinity_change,       &
                                          config_cvmix_WSwSBF_domain_y_length, &
                                          config_cvmix_WSwSBF_coriolis_gradient, &
                                          config_cvmix_WSwSBF_max_windstress_meridional, &
                                          config_cvmix_ocean_debug_tracer_xCell_center, &
                                          config_cvmix_ocean_debug_tracer_yCell_center, &
                                          config_cvmix_ocean_debug_tracer_disk_radius, &
                                          config_cvmix_WSwSBF_seaMountHeight, &
                                          config_cvmix_WSwSBF_seaMount_xCenter, &
                                          config_cvmix_WSwSBF_seaMount_yCenter, &
                                          config_cvmix_WSwSBF_seaMountWidth, &
                                          config_cvmix_WSwSBF_smooth_radius, &
                                          config_cvmix_WSwSBF_smooth_xctr, &
                                          config_cvmix_WSwSBF_smooth_yctr, &
                                          config_cvmix_WSwSBF_island_radius, &
                                          config_cvmix_WSwSBF_island_xctr, &
                                          config_cvmix_WSwSBF_island_yctr

      integer, pointer :: config_cvmix_ocean_debug_tracer_min_levels, config_cvmix_ocean_debug_tracer_max_levels	

      integer, pointer :: config_cvmix_WSwSBF_max_dz

      real(kind=RKIND) :: dist, xCenter, yCenter, tracerRadius

      ! assume no error
      iErr = 0

      ! get and test if this is the configuration specified
      call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('cvmix_WSwSBF')) return

      ! build the vertical grid
      ! intent(out) is interfaceLocations. An array ranging from 0 to 1
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_vertical_grid', config_cvmix_WSwSBF_vertical_grid)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevelsP1', nVertLevelsP1)
      allocate(interfaceLocations(nVertLevelsP1))
      call ocn_generate_vertical_grid(config_cvmix_WSwSBF_vertical_grid, interfaceLocations)

      ! load the remaining configuration parameters
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_surface_temperature', &
                                config_cvmix_WSwSBF_surface_temperature)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_surface_salinity', config_cvmix_WSwSBF_surface_salinity)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_surface_restoring_temperature', &
                                config_cvmix_WSwSBF_surface_restoring_temperature)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_surface_restoring_salinity', &
                                config_cvmix_WSwSBF_surface_restoring_salinity)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_temperature_piston_velocity', &
                                config_cvmix_WSwSBF_temperature_piston_velocity)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_salinity_piston_velocity', &
                                config_cvmix_WSwSBF_salinity_piston_velocity)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_sensible_heat_flux', config_cvmix_WSwSBF_sensible_heat_flux)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_latent_heat_flux', config_cvmix_WSwSBF_latent_heat_flux)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_shortwave_heat_flux', &
                                config_cvmix_WSwSBF_shortwave_heat_flux)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_rain_flux', config_cvmix_WSwSBF_rain_flux)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_evaporation_flux', config_cvmix_WSwSBF_evaporation_flux)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_interior_temperature_restoring_rate', &
                                config_cvmix_WSwSBF_interior_temperature_restoring_rate)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_interior_salinity_restoring_rate', &
                                config_cvmix_WSwSBF_interior_salinity_restoring_rate)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_temperature_gradient', &
                                config_cvmix_WSwSBF_temperature_gradient)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_salinity_gradient', config_cvmix_WSwSBF_salinity_gradient)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_bottom_depth', config_cvmix_WSwSBF_bottom_depth)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_max_windstress_zonal', config_cvmix_WSwSBF_max_windstress_zonal)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_max_windstress_meridional', config_cvmix_WSwSBF_max_windstress_meridional)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_coriolis_parameter', config_cvmix_WSwSBF_coriolis_parameter)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_temperature_gradient_mixed_layer', &
                config_cvmix_WSwSBF_temperature_gradient_mixed_layer)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_salinity_gradient_mixed_layer', &
                config_cvmix_WSwSBF_salinity_gradient_mixed_layer)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_mixed_layer_depth_temperature', &
                config_cvmix_WSwSBF_mixed_layer_depth_temperature)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_mixed_layer_depth_salinity', &
                config_cvmix_WSwSBF_mixed_layer_depth_salinity)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_mixed_layer_temperature_change', &
                config_cvmix_WSwSBF_mixed_layer_temperature_change)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_mixed_layer_salinity_change', &
                config_cvmix_WSwSBF_mixed_layer_salinity_change)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_coriolis_gradient', config_cvmix_WSwSBF_coriolis_gradient)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_domain_y_length', config_cvmix_WSwSBF_domain_y_length)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_wind_zonal_variation', config_cvmix_WSwSBF_wind_zonal_variation)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_ocean_debug_tracer_xCell_center', config_cvmix_ocean_debug_tracer_xCell_center)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_ocean_debug_tracer_yCell_center', config_cvmix_ocean_debug_tracer_yCell_center)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_ocean_debug_tracer_disk_radius', config_cvmix_ocean_debug_tracer_disk_radius) 
      call mpas_pool_get_config(domain % configs, 'config_cvmix_ocean_debug_tracer_min_levels', config_cvmix_ocean_debug_tracer_min_levels)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_ocean_debug_tracer_max_levels', config_cvmix_ocean_debug_tracer_max_levels)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_useSeaMount', &
                config_cvmix_WSwSBF_useSeaMount)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_seaMountHeight', &
                config_cvmix_WSwSBF_seaMountHeight)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_seaMountWidth', &
                config_cvmix_WSwSBF_seaMountWidth)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_seaMount_xCenter', &
                config_cvmix_WSwSBF_seaMount_xCenter)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_seaMount_yCenter', &
                config_cvmix_WSwSBF_seaMount_yCenter)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_bottom_noise', &
                config_cvmix_WSwSBF_bottom_noise)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_smooth_radius', &
                config_cvmix_WSwSBF_smooth_radius)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_smooth_xctr', &
                config_cvmix_WSwSBF_smooth_xctr)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_smooth_yctr', &  
                config_cvmix_WSwSBF_smooth_yctr)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_max_dz', &
                config_cvmix_WSwSBF_max_dz)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_useIsland', &
                config_cvmix_WSwSBF_useIsland)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_island_xctr', &
                config_cvmix_WSwSBF_island_xctr)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_island_yctr', &
                config_cvmix_WSwSBF_island_yctr)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_WSwSBF_island_radius', &
                config_cvmix_WSwSBF_island_radius)

      call mpas_pool_get_config(domain % configs, 'config_write_cull_cell_mask',      &
                config_write_cull_cell_mask)
      ! load data that required to initialize the ocean simulation
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)
        call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
        call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)

        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
        call mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)

        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)

        call mpas_pool_get_array(meshPool, 'fCell', fCell)
        call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
        call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

	call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
        call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
        call mpas_pool_get_array(meshPool, 'cullCell', cullCell)

        call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

        call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal, 1)
        call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional, 1)
        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersPistonVelocity', activeTracersPistonVelocity, 1)
        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersSurfaceRestoringValue', &
                                 activeTracersSurfaceRestoringValue, 1)
        call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, 'activeTracersInteriorRestoringRate', &
                                 activeTracersInteriorRestoringRate, 1)
        call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, 'activeTracersInteriorRestoringValue', &
                                 activeTracersInteriorRestoringValue, 1)
        call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
        call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
        call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
        call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
        call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)

      nCellsSolve = nCellsArray(1)
        ! Set refBottomDepth and refBottomDepthTopOfCell
        do k = 1, nVertLevels
           refBottomDepth(k) = config_cvmix_WSwSBF_bottom_depth * interfaceLocations(k+1)
           refZMid(k) = - 0.5_RKIND * config_cvmix_WSwSBF_bottom_depth * (interfaceLocations(k) + interfaceLocations(k+1))
        end do

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND

        yMin = 1.0E10_RKIND
	yMax = -1.0E10_RKIND
	dcEdgeMin = 1.0E10_RKIND
        xMin =1E10_RKIND
	yMax = -1.0E10_RKIND

	
	    ! Determine local min and max values.
 
          call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

          yMin = min( yMin, minval(yCell(1:nCellsSolve)))
          yMax = max( yMax, maxval(yCell(1:nCellsSolve)))
	  xMin = min( xMin, minval(xCell(1:nCellsSolve)))
	  xMax = max( xMax, maxval(xCell(1:nCellsSolve)))

          dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgesSolve)))


!        call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
!        call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
!        call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

        ! mark north / south boundaries
        if(abs(config_cvmix_WSwSBF_coriolis_gradient) .gt. 1.0E-12_RKIND) then
          call ocn_mark_north_boundary(meshPool, yMax, dcEdgeMin, iErr)
          call ocn_mark_south_boundary(meshPool, yMin, dcEdgeMin, iErr)
          call mpas_pool_get_config(meshPool, 'y_period', y_period)
          y_period = 0.0_RKIND
    
        call mpas_stream_mgr_begin_iteration(domain % streamManager)
        do while (mpas_stream_mgr_get_next_stream(domain % streamManager, streamID, directionProperty))
          if ( directionProperty == MPAS_STREAM_OUTPUT .or. directionProperty == MPAS_STREAM_INPUT_OUTPUT ) then
            call mpas_stream_mgr_add_att(domain % streamManager, 'y_period', 0.0_RKIND, streamID)
          end if
        end do
        end if

        if(associated(debugTracers)) then
          do iCell = 1, nCellsSolve
             do k=1,nVertLevels
                debugTracers(index_tracer1,k,iCell) = 0.0_RKIND
             end do
          end do
        do iCell=1,nCellsSolve

        if ( config_cvmix_WSwSBF_useSeaMount ) then
                bottomDepth(iCell) = config_cvmix_WSwSBF_bottom_depth - config_cvmix_WSwSBF_seaMountHeight * &
                 exp(-((yCell(iCell) - config_cvmix_WSwSBF_seaMount_xCenter)**2 /            &
                  config_cvmix_WSwSBF_seaMountWidth**2 + (xCell(iCell) - &
                  config_cvmix_WSwSBF_seaMount_yCenter)**2 / config_cvmix_WSwSBF_seaMountWidth**2))
                bottomDepth(iCell) = max(bottomDepth(iCell),0.0_RKIND)
        elseif ( config_cvmix_WSwSBF_bottom_noise ) then
                dist = sqrt((xCell(iCell) - config_cvmix_WSwSBF_smooth_xctr)**2 +   &
                          (yCell(iCell) - config_cvmix_WSwSBF_smooth_yctr)**2)
                if ( dist < config_cvmix_WSwSBF_smooth_radius ) then
                   bottomDepth(iCell) = refBottomDepth(nVertLevels - config_cvmix_WSwSBF_max_dz)
                else
                   
                   call srand(seed) 
                   randN = rand(seed)

                   if ( randN >= 0.5 ) then
                      spot = nVertLevels - 2*config_cvmix_WSwSBF_max_dz
                      bottomDepth(iCell) = refBottomDepth(spot)
                   else
                      bottomDepth(iCell) = refBottomDepth(nVertLevels)
                   endif
                endif
        else
                bottomDepth(iCell) = config_cvmix_WSwSBF_bottom_depth
        endif

        if (config_write_cull_cell_mask .and. config_cvmix_WSwSBF_useIsland) then
           dist = sqrt((xCell(iCell) - config_cvmix_WSwSBF_island_xctr)**2 +   &
                     (yCell(iCell) - config_cvmix_WSwSBF_island_yctr)**2)
           if (dist < config_cvmix_WSwSBF_island_radius) then
                        cullCell(iCell) = 1
           endif

           bottomDepth(iCell) = config_cvmix_WSwSBF_bottom_depth
        endif

        print *, iCell
        maxLevelCell(iCell) = nVertLevels
        do k=1,nVertLevels
          if (refZmid(k) < -bottomDepth(iCell)) then
            maxLevelCell(iCell) = k-1
            exit
          endif
        enddo

      enddo

          xCenter = (xMax - xMin) * config_cvmix_ocean_debug_tracer_xCell_center + xMin
	  yCenter = (yMax - yMin) * config_cvmix_ocean_debug_tracer_yCell_center + yMin
          tracerRadius = dcEdgeMin * config_cvmix_ocean_debug_tracer_disk_radius
  
   	  do iCell = 1, nCellsSolve
	     dist = sqrt((xCell(iCell) - xCenter)**2.0_RKIND + (yCell(iCell) - yCenter)**2.0_RKIND)  
             if ( dist < tracerRadius ) then
		call mpas_log_write('setting debug tracer for found match')
		do k = 1,nVertLevels
		  if (config_cvmix_ocean_debug_tracer_min_levels > config_cvmix_ocean_debug_tracer_max_levels) then
		  call mpas_log_write("Error, debug_tracer_max_levels of $i is smaller than debug_tracer_min_levels of $i", &
		  	MPAS_LOG_CRIT, intArgs = (/ config_cvmix_ocean_debug_tracer_min_levels, &
				config_cvmix_ocean_debug_tracer_max_levels/))
		  endif
		  if (config_cvmix_ocean_debug_tracer_min_levels <= k .and. &
		     k <= config_cvmix_ocean_debug_tracer_max_levels .and. &
		     k <= maxLevelCell(iCell)) then
		    debugTracers(index_tracer1,k,iCell) = 1.0_RKIND
		  endif
		enddo
	     endif
	  enddo
	endif
 
        do iCell = 1, nCellsSolve
           if(associated(activeTracers) ) then

           ! Loop from surface through surface layer depth
           k=1

           do while (k .le. maxLevelCell(iCell) .and. refZMid(k) > - config_cvmix_WSwSBF_mixed_layer_depth_temperature)
              temperature = config_cvmix_WSwSBF_surface_temperature + refZMid(k) *  &
                            config_cvmix_WSwSBF_temperature_gradient_mixed_layer
              activeTracers(index_temperature, k, iCell) = max(0.0_RKIND, temperature)
              k = k + 1
           enddo

           ! the value of k is now the first layer below the surface layer
           if ( k > 1 ) then
              temperature = activeTracers(index_temperature, k-1, iCell) + config_cvmix_WSwSBF_mixed_layer_temperature_change
              activeTracers(index_temperature, k, iCell) = max(temperature, 0.0_RKIND)
              BLdepth = refZMid(k)
           else
              activeTracers(index_temperature, k, iCell) = config_cvmix_WSwSBF_surface_temperature
              BLdepth = 0.0_RKIND
           endif

           ! find the first level below the mixed layer
           kML = k + 1

           ! now loop from the bottom of the mixed layer thru to the bottom of the domain
           do k = kML, maxLevelCell(iCell) 
              temperature = activeTracers(index_temperature, kML-1, iCell) + (refZMid(k) - BLdepth) * &
              config_cvmix_WSwSBF_temperature_gradient
              activeTracers(index_temperature, k, iCell) = max(0.0_RKIND, temperature)
           enddo

           !
           ! next compute the salinity profile
           !

           ! Loop from surface through surface layer depth
           k=1
           do while (k .le. maxLevelCell(iCell) .and. refZMid(k) > - config_cvmix_WSwSBF_mixed_layer_depth_salinity)
              salinity = config_cvmix_WSwSBF_surface_salinity + refZMid(k) * config_cvmix_WSwSBF_salinity_gradient_mixed_layer
              activeTracers(index_salinity, k, iCell) = salinity
              k = k + 1
           enddo

           ! the value of k is now the first layer below the surface layer
            if ( k > 1 ) then
               salinity = activeTracers(index_salinity, k-1, iCell) + config_cvmix_WSwSBF_mixed_layer_salinity_change
               activeTracers(index_salinity, k, iCell) = salinity
               BLdepth = refZMid(k)
            else
               activeTracers(index_salinity, k, iCell) = config_cvmix_WSwSBF_surface_salinity
               BLdepth = 0.0_RKIND
            endif

            ! find the first level below the mixed layer
            kML = k + 1

            ! now loop from the bottom of the mixed layer thru to the bottom of the domain
            do k = kML, maxLevelCell(iCell) 
               salinity = activeTracers(index_salinity, kML-1, iCell) + (refZMid(k) - BLdepth) *  &
               config_cvmix_WSwSBF_salinity_gradient
               activeTracers(index_salinity, k, iCell) = salinity
            enddo

            activeTracers(:,maxLevelCell(iCell)+1:nVertLevels,iCell) = -1.0E34
        endif ! if (associated(activeTracer))

        ! as a place holder, have some debug tracer in the top few layers and zero below
!        if ( associated(debugTracers) ) then
!           debugTracers(index_tracer1, k, iCell) = 0.0_RKIND
!           do k=1,min(4,nVertLevels)
!              debugTracers(index_tracer1, k, iCell) = 1.0_RKIND
!           enddo
!        endif

        ! set bottom depth
!        if ( config_cvmix_WSwSBF_useSeaMount ) then
!          bottomDepth(iCell) = config_cvmix_WSwSBF_bottom_depth - config_cvmix_WSwSBF_seaMountHeight * &
!                           exp(-((yCell(iCell) - config_cvmix_WSwSBF_seaMount_xCenter)**2 /            &
!                            config_cvmix_WSwSBF_seaMountWidth**2 + (xCell(iCell) - &
!                            config_cvmix_WSwSBF_seaMount_yCenter)**2 / config_cvmix_WSwSBF_seaMountWidth**2))
!          bottomDepth(iCell) = max(bottomDepth(iCell),0.0_RKIND)
!        else
!          bottomDepth(iCell) = config_cvmix_WSwSBF_bottom_depth
!        endif

        layerThickness(:,iCell) = 0.0_RKIND
        ! Set layerThickness
        layerThickness(1,iCell) = refBottomDepth(1)
        do k = 2, maxLevelCell(iCell) 
          
           layerThickness(k, iCell) = abs(refBottomDepth(k-1) - refBottomDepth(k))
           restingThickness(k, iCell) = layerThickness(k, iCell)
        end do

          print *, '******************'
          print *, sum(layerThickness(:,ICell)),bottomDepth(iCell)
          print *, '*****************'
                  bottomDepth(iCell) = sum(layerThickness(:,iCell))
        ! Set surface temperature restoring value and rate
        ! Value in units of C, piston velocity in units of m/s
        if ( associated(activeTracersSurfaceRestoringValue) ) then
           activeTracersSurfaceRestoringValue(index_temperature, iCell) = config_cvmix_WSwSBF_surface_restoring_temperature
        end if
        if ( associated(activeTracersPistonVelocity) ) then
           activeTracersPistonVelocity(index_temperature, iCell) = config_cvmix_WSwSBF_temperature_piston_velocity
        end if

        ! Set surface salinity restoring value and rate
        ! Value in units of PSU, piston velocity in units of m/s
        if ( associated(activeTracersSurfaceRestoringValue) ) then
           activeTracersSurfaceRestoringValue(index_salinity, iCell) = config_cvmix_WSwSBF_surface_restoring_salinity
        end if
        if ( associated(activeTracersPistonVelocity) ) then
           activeTracersPistonVelocity(index_salinity, iCell) = config_cvmix_WSwSBF_salinity_piston_velocity
        end if

        ! Set sensible heat flux
        sensibleHeatFlux(iCell) = config_cvmix_WSwSBF_sensible_heat_flux

        ! Set latent heat flux
        latentHeatFlux(iCell) = config_cvmix_WSwSBF_latent_heat_flux

        ! Set shortwave heat flux
        shortWaveHeatFlux(iCell) = config_cvmix_WSwSBF_shortwave_heat_flux

        ! Set precipation and evaporation
        rainFlux(iCell) = config_cvmix_WSwSBF_rain_flux
        evaporationFlux(iCell) = config_cvmix_WSwSBF_evaporation_flux

        ! Set interior temperature restoring value and rate
        do k = 1, maxLevelCell(iCell) 
           if ( associated(activeTracersInteriorRestoringValue) ) then
              activeTracersInteriorRestoringValue(index_temperature, k, iCell) = activeTracers(index_temperature, k, iCell)
           end if
           if ( associated(activeTracersInteriorRestoringRate) ) then
              activeTracersInteriorRestoringRate(index_temperature, k, iCell) = &
                                config_cvmix_WSwSBF_interior_temperature_restoring_rate
           end if
        enddo

        ! Set interior salinity restoring value and rate
        do k = 1, maxLevelCell(iCell) 
           if ( associated(activeTracersInteriorRestoringValue) ) then
              activeTracersInteriorRestoringValue(index_salinity, k, iCell) = activeTracers(index_salinity, k, iCell)
           end if
           if ( associated(activeTracersInteriorRestoringRate) ) then
              activeTracersInteriorRestoringRate(index_salinity, k, iCell) = config_cvmix_WSwSBF_interior_salinity_restoring_rate
           end if
        enddo

        ! Set Coriolis parameter
        fCell(iCell) = config_cvmix_WSwSBF_coriolis_parameter

        ! Set bottomDepth
!        bottomDepth(iCell) = config_cvmix_WSwSBF_bottom_depth

        ! Set maxLevelCell

     end do  ! do iCell

     do iCell = 1, nCellsSolve
        windStressZonal(iCell) = config_cvmix_WSwSBF_max_windstress_zonal
        windStressMeridional(iCell) = config_cvmix_WSwSBF_max_windstress_meridional
     enddo

     if(config_cvmix_WSwSBF_wind_zonal_variation) then
	config_cvmix_WSwSBF_domain_y_length = yMax - yMin
       do iCell = 1, nCellsSolve
          windStressZonal(iCell) = config_cvmix_WSwSBF_max_windstress_zonal * sin(2.0_RKIND*3.14159265_RKIND &
                                       * yCell(iCell) / config_cvmix_WSwSBF_domain_y_length)
	  windStressMeridional(iCell) = 0.0_RKIND
       enddo
     endif

     do iCell = 1, nCellsSolve
        fCell(iCell) = config_cvmix_WSwSBF_coriolis_parameter + config_cvmix_WSwSBF_coriolis_gradient * yCell(iCell)
     enddo

     do iEdge = 1, nEdgesSolve
        fEdge(iEdge) = config_cvmix_WSwSBF_coriolis_parameter + config_cvmix_WSwSBF_coriolis_gradient * yEdge(iEdge)
     end do

     do iVertex=1, nVerticesSolve
        fVertex(iVertex) = config_cvmix_WSwSBF_coriolis_parameter + config_cvmix_WSwSBF_coriolis_gradient * yVertex(iVertex)
     end do

     block_ptr => block_ptr % next
   end do

   deallocate(interfaceLocations)

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_cvmix_WSwSBF!}}}

!***********************************************************************
!
!  routine ocn_init_validate_cvmix_WSwSBF
!
!> \brief   Validation for CVMix WSwSBF mixing unit test case
!> \author  Doug Jacobsen
!> \date    04/01/2015
!> \details
!>  This routine validates the configuration options for the CVMix WSwSBF mixing unit test configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_cvmix_WSwSBF(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_cvmix_WSwSBF_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('cvmix_WSwSBF')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_cvmix_WSwSBF_vert_levels', config_cvmix_WSwSBF_vert_levels)

      if(config_vert_levels <= 0 .and. config_cvmix_WSwSBF_vert_levels > 0) then
         config_vert_levels = config_cvmix_WSwSBF_vert_levels
      else if(config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for CVMix WSwSBF unit test case. Not given a usable value for ' &
                          // 'vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_cvmix_WSwSBF!}}}

!***********************************************************************

end module ocn_init_cvmix_WSwSBF

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
