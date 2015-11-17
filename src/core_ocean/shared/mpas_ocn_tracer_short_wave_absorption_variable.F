!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_short_wave_absorption_jerlov
!
!> \brief MPAS ocean tracer short wave
!> \author Doug Jacobsen
!> \date   12/17/12
!> \details
!>  This module contains the routine for computing 
!>  short wave tendencies using Jerlov
!
!-----------------------------------------------------------------------

module ocn_tracer_short_wave_absorption_variable

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_forcing
   use mpas_stream_manager
   use ocn_constants

   implicit none
   type(MPAS_forcing_group_type),pointer :: shortwaveForcingGroup
   
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

   public :: ocn_tracer_short_wave_absorption_variable_tend, &
             ocn_tracer_short_wave_absorption_variable_init, &
             ocn_get_variable_sw_fraction,					 &
             ocn_get_m05_coeffs, 							 &
             ocn_get_os00_coeffs,							 &
             ocn_get_o03_coeffs,							 &
             init_shortwave_forcing_m05,					 &
             init_shortwave_forcing_ohlmann,				 &
             get_shortWaveData,								 &
             shortwave_forcing_write_restart

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_short_wave_absorption_jerlov_tend
!
!> \brief   Computes tendency term for surface fluxes
!> \author  Luke Van Roekel
!> \date    11/10/2015
!> \details 
!>  This routine computes the tendency for tracers based on surface fluxes.
!>  This computation is now based on spatially variable chlorophyll, cloud fraction, and zenith angle
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_short_wave_absorption_variable_tend(meshPool, swForcingPool, forcingPool, index_temperature, layerThickness, penetrativeTemperatureFlux, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! LPV for spatially variable, will also need hooks to the forcingPool to get the extinction coefficients for the blue and red
      ! bands (extinction_coeffs_red and extinction_coeffs_bg)
      ! This will have to be updated with times.  It should be okay to just update immediately between months.

      type (mpas_pool_type), intent(in) :: &
         meshPool,						   &          !< Input: mesh information
         swForcingPool,					   &			  !< Input: chlorophyll, cloud, zenith data
         forcingPool
         
      real (kind=RKIND), dimension(:), intent(in) :: &
        penetrativeTemperatureFlux !< Input: penetrative temperature flux through the surface

      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Layer thicknesses

      integer, intent(in) :: index_temperature

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

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

      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND) :: depth
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth
      real (kind=RKIND), dimension(:), allocatable :: weights
	  real (kind=RKIND), dimension(:), pointer :: chlorophyllA, zenithAngle, clearSkyRadiation, shortWaveHeatFlux
	  character (len=StrKIND), pointer :: config_sw_absorption_type
	  real (kind=RKIND), dimension(4) :: Avals, Kvals
	  real (kind=RKIND) :: cloudRatio ! cloud Ratio = 1 - incident_sfc_sw_radiation/clearSkyRadiation
	  
      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

      allocate(weights(nVertLevels+1))
      weights = 0.0_RKIND
      weights(1) = 1.0_RKIND
	  Avals(:)=0.0_RKIND
	  Kvals(:)=0.0_RKIND

	  call mpas_pool_get_config(ocnConfigs, 'config_sw_absorption_type', config_sw_absorption_type)
	  call mpas_pool_get_array(swForcingPool,'Chlorophyll',chlorophyllA)
	  call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
	  
	  select case ( trim(config_sw_absorption_type) ) 
      case('manizza05')
         !$omp do schedule(runtime) private(depth, k)
         do iCell = 1, nCells
           depth = 0.0_RKIND
           call ocn_get_m05_coeffs(chlorophyllA(iCell),Avals,Kvals)
           
           do k = 1, maxLevelCell(iCell)
             depth = depth + layerThickness(k,iCell)
			 
             call ocn_get_variable_sw_fraction(depth, weights(k+1),Avals,Kvals)
             tend(index_temperature, k, iCell) = tend(index_temperature, k, iCell) + penetrativeTemperatureFlux(iCell)*(weights(k) - weights(k+1))
           end do
         end do
         !$omp end do
      case('ohlmann00')
         call mpas_pool_get_array(swForcingPool,'ZenithAngle',zenithAngle)
         call mpas_pool_get_array(swForcingPool,'clearSkyRadiation',clearSkyRadiation)
         !$omp do schedule(runtime) private(depth, k)
         do iCell = 1, nCells
           depth = 0.0_RKIND
           cloudRatio = min(1.0_RKIND, 1.0_RKIND - shortWaveHeatFlux(iCell)/clearSkyRadiation(iCell))
           cloudRatio = max(0.0_RKIND, cloudRatio)
           
           call ocn_get_os00_coeffs(chlorophyllA(iCell),zenithAngle(iCell),cloudRatio,Avals, Kvals)
           
           do k = 1, maxLevelCell(iCell)
             depth = depth + layerThickness(k, iCell)

             call ocn_get_variable_sw_fraction(depth, weights(k+1), Avals, Kvals)
             tend(index_temperature, k, iCell) = tend(index_temperature, k, iCell) + penetrativeTemperatureFlux(iCell)*(weights(k) - weights(k+1))
           end do
         end do
         !$omp end do
      case('ohlmann03')
         call mpas_pool_get_array(swForcingPool,'zenithAngle',zenithAngle)
         call mpas_pool_get_array(swForcingPool,'clearSkyRadiation',clearSkyRadiation)
         !$omp do schedule(runtime) private(depth, k)
         do iCell = 1, nCells
           cloudRatio = min(1.0_RKIND, 1.0_RKIND - shortWaveHeatFlux(iCell)/clearSkyRadiation(iCell))
           cloudRatio = max(0.0_RKIND, cloudRatio)
           
           call ocn_get_o03_coeffs(chlorophyllA(iCell),zenithAngle(iCell),cloudRatio, Avals, Kvals)
           depth = 0.0_RKIND
           do k = 1, maxLevelCell(iCell)
             depth = depth + layerThickness(k, iCell)

             call ocn_get_variable_sw_fraction(depth, weights(k+1), Avals, Kvals)
             tend(index_temperature, k, iCell) = tend(index_temperature, k, iCell) + penetrativeTemperatureFlux(iCell)*(weights(k) - weights(k+1))
           end do
         end do
         !$omp end do
      end select

      deallocate(weights)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_short_wave_absorption_variable_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_short_wave_absorption_variable_init
!
!> \brief   Initializes ocean tracer surface flux quantities
!> \author  Luke Van Roekel
!> \date    11/10/15
!> \details 
!>  This routine initializes quantities related to surface fluxes in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_short_wave_absorption_variable_init(domain,err)!{{{

   !--------------------------------------------------------------------

	  type(domain_type) :: domain
	  
      integer, intent(out) :: err !< Output: error flag

      character (len=StrKIND), pointer :: config_sw_absorption_type

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_sw_absorption_type', config_sw_absorption_type)

	  select case ( trim(config_sw_absorption_type) ) 
      case ('manizza05')
         call init_shortwave_forcing_m05(domain)
      case ('ohlmann00')
         call init_shortwave_forcing_ohlmann(domain)
      case ('ohlmann03')
      	 call init_shortwave_forcing_ohlmann(domain)
      case default
         write(stderrUnit,*) "Shortwave parameterization type unknown: ", trim(config_sw_absorption_type)
         call MPAS_dmpar_global_abort('Unknown config_sw_absortion_type')
      end select
          

   end subroutine ocn_tracer_short_wave_absorption_variable_init!}}}

!***********************************************************************

!***********************************************************************
!
!  routine init_shortwave_forcing_m05
!
!> \brief   Initializes forcing group if parameterization of Manizza et al. (2005)
!> 			This parameterization only requires chlorophyll-a concentrations, so only add that
!> \author  Luke Van Roekel
!> \date    11/10/15
!> \details 
!>  This routine initializes forcing stream for Manizza et al. (2005) parameterization
!
!-----------------------------------------------------------------------

   subroutine init_shortwave_forcing_m05(domain)
   
   type(domain_type) :: domain
   
   character(len=strKind), pointer :: &
   		config_forcing_restart_file
   		
   	logical, pointer :: &
   		config_do_restart
   		
   	character(len=strKIND) :: &
   		forcingIntervalMonthly,  &
   		forcingReferenceTimeMonthly
   		
   	call MPAS_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
   	call MPAS_pool_get_config(domain % configs, 'config_forcing_restart_file', config_forcing_restart_file)
   	
   	forcingIntervalMonthly = "0000-01-00_00:00:00"
   	forcingReferenceTimeMonthly = "0000-01-01_00:00:00"
   	
   	call MPAS_forcing_init_group( shortwaveForcingGroup,  &
   		"shortwave_monthly_observations", &
   		domain, &
   		'0000-01-01_00:00:00', &
   		'0000-01-01_00:00:00', &
   		'0001-00-00_00:00:00', &
   		config_do_restart, &
   		config_forcing_restart_file)
   		
   	call MPAS_forcing_init_field( domain % streamManager, &
   		shortwaveForcingGroup, &
   		'shortwave_monthly_observations', &
   		'Chlorophyll', &
   		'shortwave_forcing_data', &
   		'shortwavePool',  &
   		'Chlorophyll',  &
   		'linear',  &
   		forcingReferenceTimeMonthly,  &
   		forcingIntervalMonthly)
   		
   	call MPAS_forcing_init_field_data( shortwaveForcingGroup, &
   		'shortwave_monthly_observations', &
   		domain % streamManager)
   		
   	end subroutine init_shortwave_forcing_m05
   	
!***********************************************************************

!***********************************************************************
!
!  routine init_shortwave_forcing_ohlmann
!
!> \brief   Initializes forcing group if parameterization of Ohlmann and Siegel (2000) or Ohlmann(2003)
!> 			This parameterization only requires chlorophyll-a concentrations, so only add that
!> \author  Luke Van Roekel
!> \date    11/10/15
!> \details 
!>  This routine initializes forcing stream for Ohlmann and Siegel (2000) or Ohlman (2003) parameterization
!
!-----------------------------------------------------------------------

   subroutine init_shortwave_forcing_ohlmann(domain)
   
   type(domain_type) :: domain
   
   character(len=strKind), pointer :: &
   		config_forcing_restart_file
   		
   	logical, pointer :: &
   		config_do_restart
   		
   	character(len=strKIND) :: &
   		forcingIntervalMonthly,  &
   		forcingReferenceTimeMonthly
   		
   	call MPAS_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
   	call MPAS_pool_get_config(domain % configs, 'config_forcing_restart_file', config_forcing_restart_file)
   	
   	forcingIntervalMonthly = "0000-01-00_00:00:00"
   	forcingReferenceTimeMonthly = "0000-01-01_00:00:00"
   	
   	call MPAS_forcing_init_group( shortwaveForcingGroup,  &
   		"shortwave_monthly_observations", &
   		domain, &
   		'0000-01-01_00:00:00', &
   		'0000-01-01_00:00:00', &
   		'0001-00-00_00:00:00', &
   		config_do_restart, &
   		config_forcing_restart_file)
   		
   	call MPAS_forcing_init_field( domain % streamManager, &
   		shortwaveForcingGroup, &
   		'shortwave_monthly_observations', &
   		'Chlorophyll', &
   		'shortwave_forcing_data', &
   		'shortwavePool',  &
   		'Chlorophyll',  &
   		'linear',  &
   		forcingReferenceTimeMonthly,  &
   		forcingIntervalMonthly)
       	
    call MPAS_forcing_init_field( domain % streamManager, &
   		shortwaveForcingGroup, &
   		'shortwave_monthly_observations', &
   		'CloudFraction', &
   		'shortwave_forcing_data', &
   		'shortwavePool',  &
   		'CloudFraction',  &
   		'linear',  &
   		forcingReferenceTimeMonthly,  &
   		forcingIntervalMonthly)
  
    call MPAS_forcing_init_field( domain % streamManager, &
   		shortwaveForcingGroup, &
   		'shortwave_monthly_observations', &
   		'ZenithAngle', &
   		'shortwave_forcing_data', &
   		'shortwavePool',  &
   		'ZenithAngle',  &
   		'linear',  &
   		forcingReferenceTimeMonthly,  &
   		forcingIntervalMonthly)
   		 		
   	call MPAS_forcing_init_field_data( shortwaveForcingGroup, &
   		'shortwave_monthly_observations', &
   		domain % streamManager)
   		
   	end subroutine init_shortwave_forcing_ohlmann
   	   	
!***********************************************************************

!***********************************************************************
!
!  routine get_shortWaveData
!
!> \brief   retrieve data needed to compute penetration of shortwave radiation
!> \author  Luke Van Roekel
!> \date    11/10/15
!> \details 
!>  This routine calls mpas_forcing routines to acquire needed shortwave data and interpolates
!>    between time levels
!
!-----------------------------------------------------------------------

    subroutine get_shortWaveData( streamManager, &
    	domain, &
    	simulationClock, &
    	firstTimeStep)
    	
    	type (MPAS_streamManager_type), intent(inout) :: streamManager
    	
    	type (domain_type) :: domain
    	
    	type (MPAS_clock_type) :: simulationClock
    	
    	logical, intent(in) :: firstTimeStep
    	character(len=strKind), pointer :: config_sw_absorption_type
    	real(kind=RKIND), pointer :: config_dt
    	
    	call MPAS_pool_get_config(domain%configs, 'config_sw_absorption_type', config_sw_absorption_type)
    	call MPAS_pool_get_config(domain%configs, 'config_dt', config_dt)
    	
    	if(trim(config_sw_absorption_type) == 'manizza05' .or. trim(config_sw_absorption_type) == 'ohlmann00' &
    	     .or. trim(config_sw_absorption_type) == 'ohlmann03') then
    	     call MPAS_forcing_get_forcing(shortwaveForcingGroup, &
    	     		'shortwave_monthly_observations', streamManager, config_dt)
    	endif
    	
    end subroutine get_shortWaveData
    

!***********************************************************************	
    	
    	
!***********************************************************************
!
!  routine ocn_get_variable_fractions
!
!> \brief   Computes short wave absorption fractions
!> \author  Luke Van Roekel
!> \date    11/10/2015
!> \details 
!>  Computes fraction of solar short-wave flux penetrating to
!>  specified depth due to time and space varying chlorophyll, cloud fraction, and zenith angle
!> based on:
!>     Ohlmann and Siegel (2000), Ohlmann (2003), Manizza et al. (2005)

!
!-----------------------------------------------------------------------
   subroutine ocn_get_variable_sw_fraction(depth, weight, Avals, Kvals)!{{{
!  Note: below 200m the solar penetration gets set to zero,
!     otherwise the limit for the exponent ($+/- 5678$) needs to be 
!     taken care of.

      real (kind=RKIND), intent(in) :: depth !< Input: Depth of bottom of cell
      real (kind=RKIND), intent(in),dimension(4) :: Avals  !< Input: spectral partitioning of radiation
      real (kind=RKIND), intent(in),dimension(4) :: Kvals  !< Input: extinction coefficients for different radiation bands
      real (kind=RKIND), intent(out) :: weight !< Output: Weight for Jerlov absorption

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!
      integer :: k

      real (kind=RKIND), parameter :: depth_cutoff = -200.0_RKIND
   
!-----------------------------------------------------------------------
!
!  compute absorption fraction
!
!-----------------------------------------------------------------------

      if (-depth < depth_cutoff) then
         weight = 0.0_RKIND
      else
         weight=0.0_RKIND
         do k=1,4
            weight = weight + Avals(k)*exp(-depth*Kvals(k))
         enddo
      endif
   end subroutine ocn_get_variable_sw_fraction!}}}
   
!***********************************************************************	
    	
    	
!***********************************************************************
!
!  routine ocn_get_m05_coeffs
!
!> \brief   Computes coefficients for spatially varying penetrating shortwave via Manizza et al 2005
!> \author  Luke Van Roekel
!> \date    11/10/2015
!> \details 
!>   This will fill in coefficients for the penetrating shortwave parameterization
!> 		A1*exp(-K1*depth) + A2*exp(-K2*depth) + A3*exp(-K3*depth) + A4*exp(-K4*depth)
!>   For manizza et al (2005) A4 is set to zero
!>   The visible portion of the spectrum is divided into red and blue/green bands

!
!-----------------------------------------------------------------------   

   subroutine ocn_get_m05_coeffs(chlorophyllA,Avals,Kvals)
   
      real(kind=RKIND), intent(in) :: chlorophyllA
      real(kind=RKIND), intent(out), dimension(4) :: Avals, Kvals
      
      
      real(kind=RKIND), parameter ::  &
      	kswr=0.225,		&
      	kswb=0.0232,	&
      	chir=0.037,		&
      	chib=0.074,		&
      	ered=0.629,		&
      	eblue=0.674
      	
      Avals(1)=0.52_RKIND
      Avals(2)=0.24_RKIND
      Avals(3)=0.24_RKIND
      Avals(4)=0.0_RKIND
      Kvals(1)=2.86_RKIND
      Kvals(2)=kswr + chir*chlorophyllA**ered     
      Kvals(3)=kswb + chib*chlorophyllA**eblue
    
    end subroutine ocn_get_m05_coeffs

!***********************************************************************	
    	
    	
!***********************************************************************
!
!  routine ocn_get_os00_coeffs
!
!> \brief   Computes coefficients for spatially varying penetrating shortwave via Ohlmann and Siegel 2000
!> \author  Luke Van Roekel
!> \date    11/10/2015
!> \details 
!>   This will fill in coefficients for the penetrating shortwave parameterization
!> 		A1*exp(-K1*depth) + A2*exp(-K2*depth) + A3*exp(-K3*depth) + A4*exp(-K4*depth)
!>   For Ohlmann and Siegel (2000) A4 no approximations are made
!>   Here the IR portion of the spectrum is also decomposed.  Four exponential terms are used 

!
!-----------------------------------------------------------------------   
  
  subroutine ocn_get_os00_coeffs(chlorophyllA,zenithAngle,cloudFraction,Avals,Kvals)
  
     real(kind=RKIND), intent(in) :: chlorophyllA, zenithAngle, cloudFraction
     real(kind=RKIND), intent(out), dimension(4) :: Avals, Kvals
     
     if(cloudFraction > 0.1) then ! cloudy skies
     	Avals(1) = 0.026_RKIND*chlorophyllA + 0.112_RKIND*cloudFraction + 0.366_RKIND  
     	Avals(2) = -0.009_RKIND*chlorophyllA + 0.034_RKIND*cloudFraction + 0.207_RKIND  
     	Avals(3) = -0.015_RKIND*chlorophyllA -0.006_RKIND*cloudFraction + 0.188_RKIND  
     	Avals(4) = -0.003_RKIND*chlorophyllA -0.131_RKIND*cloudFraction + 0.169_RKIND  
     	Kvals(1) = 0.063_RKIND*chlorophyllA -0.015_RKIND*cloudFraction + 0.082_RKIND  
     	Kvals(2) = 0.278_RKIND*chlorophyllA -0.562_RKIND*cloudFraction + 1.02_RKIND  
     	Kvals(3) = 3.91_RKIND*chlorophyllA -12.91_RKIND*cloudFraction + 16.62_RKIND  
     	Kvals(4) = 16.64_RKIND*chlorophyllA -478.28_RKIND*cloudFraction + 736.56_RKIND  
     else ! clear skies
     	Avals(1) = 0.033_RKIND*chlorophyllA -0.025_RKIND*zenithAngle + 0.419_RKIND  
     	Avals(2) = -0.010_RKIND*chlorophyllA -0.007_RKIND*zenithAngle + 0.231_RKIND  
     	Avals(3) = -0.019_RKIND*chlorophyllA -0.003_RKIND*zenithAngle + 0.195_RKIND  
     	Avals(4) = -0.006_RKIND*chlorophyllA -0.004_RKIND*zenithAngle + 0.154_RKIND  
     	Kvals(1) = 0.066_RKIND*chlorophyllA + 0.006_RKIND*zenithAngle + 0.066_RKIND  
     	Kvals(2) = 0.396_RKIND*chlorophyllA -0.027_RKIND*zenithAngle + 0.866_RKIND  
     	Kvals(3) = 7.68_RKIND*chlorophyllA -2.49_RKIND*zenithAngle + 17.81_RKIND  
     	Kvals(4) = 51.27_RKIND*chlorophyllA + 13.14_RKIND*zenithAngle + 665.19_RKIND  
     endif
   
   end subroutine ocn_get_os00_coeffs
     
!***********************************************************************	
    	
    	
!***********************************************************************
!
!  routine ocn_get_o03_coeffs
!
!> \brief   Computes coefficients for spatially varying penetrating shortwave via Ohlmann (2003)
!> \author  Luke Van Roekel
!> \date    11/10/2015
!> \details 
!>   This will fill in coefficients for the penetrating shortwave parameterization
!> 		A1*exp(-K1*depth) + A2*exp(-K2*depth) + A3*exp(-K3*depth) + A4*exp(-K4*depth)
!>   For Ohlmann (2003) A3 and A4 are set to zero
!>   This parameterization allows chlA to influence the IR, but is an approximation of OS00 (above)
!>   NOTE: Ohlmann (2003) states explicitly that this parameterization is NOT valid for depths
!>             less than two meters

!
!-----------------------------------------------------------------------   

   subroutine ocn_get_o03_coeffs(chlorophyllA,zenithAngle,cloudFraction,Avals,Kvals)
   
     real(kind=RKIND), intent(in) :: chlorophyllA, zenithAngle, cloudFraction
     real(kind=RKIND), intent(inout), dimension(4) :: Avals, Kvals
     

	 Avals(1)=0.571_RKIND + 0.025_RKIND*log(0.149_RKIND*chlorophyllA)
   	 Avals(2)=0.223_RKIND + 0.010_RKIND*log(2.329_RKIND*chlorophyllA)
   		
   	 Kvals(1)=0.015_RKIND + 0.176_RKIND*sqrt(0.462_RKIND*chlorophyllA)
   	 Kvals(2)=0.688_RKIND + 0.060_RKIND*log(0.125_RKIND*chlorophyllA)
   		
   	 if(cloudFraction <= 0.1) then
   	 	Avals(1) = Avals(1)*(0.66_RKIND + 0.42_RKIND*zenithAngle)
   	 	Avals(2) = Avals(2)*(0.66_RKIND + 0.42_RKIND*zenithAngle)
   	 endif
   end subroutine ocn_get_o03_coeffs

!***********************************************************************
!
!  routine shortwave_forcing_write_restart
!
!> \brief   writes restart timestamp for SW data to be read in on future restart
!> \author  Luke Van Roekel
!> \date    11/16/2015

!
!-----------------------------------------------------------------------   
   	    	
   subroutine shortwave_forcing_write_restart(domain)
   
      type(domain_type) :: domain
      
      character(len=strKind), pointer :: config_forcing_restart_file,   &
      		config_sw_absorption_type
      
     
      call MPAS_pool_get_config(domain % configs, "config_sw_absorption_type", config_sw_absorption_type)
      
      if( trim(config_sw_absorption_type) == 'manizza05' .or. trim(config_sw_absorption_type) == 'ohlmann00'  &
      		.or. trim(config_sw_absorption_type) == 'ohlmann03') then
      		call MPAS_pool_get_config(domain % configs, "config_forcing_restart_file", config_forcing_restart_file)
      		
      		call MPAS_forcing_write_restart_times(shortwaveForcingGroup, config_forcing_restart_file)
      		
      endif
      
    end subroutine shortwave_forcing_write_restart
    
    
end module ocn_tracer_short_wave_absorption_variable


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
