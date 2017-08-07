! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_surface_flux_to_tend
!
!> \brief MPAS ocean tracer surface flux
!> \author Doug Jacobsen
!> \date   12/17/12
!> \details
!>  This module contains the routine for computing
!>  surface flux tendencies.
!
!-----------------------------------------------------------------------

module ocn_tracer_surface_flux_to_tend

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines

   use ocn_constants
   use ocn_forcing

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

   public :: ocn_tracer_surface_flux_tend, &
             ocn_correct_surface_heat_flux_tend, &
             ocn_tracer_surface_flux_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: surfaceTracerFluxOn
   logical :: correct_surfaceHeatFlux_RainEvap
   logical :: correct_surfaceHeatFlux_Runoff
   logical :: correct_surfaceHeatFlux_SeaIce

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_surface_flux_tend
!
!> \brief   Computes tendency term for surface fluxes
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \details
!>  This routine computes the tendency for tracers based on surface fluxes.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_surface_flux_tend(meshPool, fractionAbsorbed, fractionAbsorbedRunoff, layerThickness,  &
      surfaceTracerFlux, surfaceTracerFluxRunoff, tend, err)!{{{
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        layerThickness !< Input: Layer thickness

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        surfaceTracerFlux !< Input: surface tracer fluxes

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        fractionAbsorbed !< Input: Coefficients for the application of surface fluxes

      real (kind=RKIND), dimension(:,:), intent(in), pointer :: &
        surfaceTracerFluxRunoff !< Input: surface tracer fluxes from river runoff

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        fractionAbsorbedRunoff !< Input: Coefficients for the application of surface fluxes due to river runoff

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

      integer :: iCell, k, iTracer, nTracers, nCells
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray
      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND) :: remainingFlux

      err = 0

      if (.not. surfaceTracerFluxOn) return

      call mpas_timer_start("surface_tracer_flux")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      nTracers = size(tend, dim=1)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      nCells = nCellsArray( 1 )

      !$omp do schedule(runtime) private(remainingFlux, k, iTracer)
      do iCell = 1, nCells
        remainingFlux = 1.0_RKIND
        do k = 1, maxLevelCell(iCell)
          remainingFlux = remainingFlux - fractionAbsorbed(k, iCell)

          do iTracer = 1, nTracers
            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + surfaceTracerFlux(iTracer, iCell) * fractionAbsorbed(k, iCell)
          end do
        end do

        if(maxLevelCell(iCell) > 0 .and. remainingFlux > 0.0_RKIND) then
          do iTracer = 1, nTracers
            tend(iTracer, maxLevelCell(iCell), iCell) = tend(iTracer, maxLevelCell(iCell), iCell) &
                                                      + surfaceTracerFlux(iTracer, iCell) * remainingFlux
          end do
        end if
      end do
      !$omp end do

      call mpas_timer_stop("surface_tracer_flux")

      ! now do runoff component

      if (associated(surfaceTracerFluxRunoff)) then
        call mpas_timer_start("surface_tracer_runoff_flux")

        !$omp do schedule(runtime) private(remainingFlux, k, iTracer)
        do iCell = 1, nCells
          remainingFlux = 1.0_RKIND
          do k = 1, maxLevelCell(iCell)
            remainingFlux = remainingFlux - fractionAbsorbedRunoff(k, iCell)

            do iTracer = 1, nTracers
              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) +  &
                 surfaceTracerFluxRunoff(iTracer, iCell) * fractionAbsorbedRunoff(k, iCell)
            end do
          end do

          if(maxLevelCell(iCell) > 0 .and. remainingFlux > 0.0_RKIND) then
            do iTracer = 1, nTracers
              tend(iTracer, maxLevelCell(iCell), iCell) = tend(iTracer, maxLevelCell(iCell), iCell) &
                                                        + surfaceTracerFluxRunoff(iTracer, iCell) * remainingFlux
            end do
          end if
        end do
        !$omp end do

        call mpas_timer_stop("surface_tracer_runoff_flux")
      end if


   !--------------------------------------------------------------------

   end subroutine ocn_tracer_surface_flux_tend!}}}

!***********************************************************************
!
!  routine ocn_correct_surface_heat_flux_tend
!
!> \brief   Computes tendency term for surface fluxes
!> \author  Mark Petersen
!> \date    March 16 2017
!> \details
!>  This routine removes the average heat flux from all cells due to rain, evap,
!>  runoff, and sea ice.
!
!-----------------------------------------------------------------------

   subroutine ocn_correct_surface_heat_flux_tend(meshPool, fractionAbsorbed, layerThickness,  &
      surfaceTracerFlux, indexTemperature, surfaceHeatFluxSums, tend, err)!{{{
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      type (mpas_pool_type), pointer :: tracersPool, tracersTendPool
! tracers and their tendency

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        layerThickness !< Input: Layer thickness

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        surfaceTracerFlux !< Input: surface tracer fluxes

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        fractionAbsorbed !< Input: Coefficients for the application of surface fluxes

      integer, intent(in) :: indexTemperature

      real (kind=RKIND), dimension(:), intent(in) :: surfaceHeatFluxSums

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

      integer :: iCell, k, iTracer, nTracers, nCells
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), pointer :: areaCellTotal

      real (kind=RKIND) :: remainingFlux, surfaceHeatFluxCorrection

      err = 0
      if (.not. surfaceTracerFluxOn) return

      call mpas_timer_start("correct_surface_heat_flux_tend")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool,'areaCellTotal',areaCellTotal)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      nTracers = size(tend, dim=1)

      nCells = nCellsArray( 1 )
      iTracer = indexTemperature 
      surfaceHeatFluxCorrection = 0.0_RKIND
      if (correct_surfaceHeatFlux_RainEvap) then
         surfaceHeatFluxCorrection = surfaceHeatFluxCorrection + surfaceHeatFluxSums(1)/areaCellTotal
      endif
      if (correct_surfaceHeatFlux_Runoff) then
         surfaceHeatFluxCorrection = surfaceHeatFluxCorrection + surfaceHeatFluxSums(2)/areaCellTotal
      endif
      if (correct_surfaceHeatFlux_SeaIce) then
         surfaceHeatFluxCorrection = surfaceHeatFluxCorrection + surfaceHeatFluxSums(3)/areaCellTotal
      endif

      !$omp do schedule(runtime) private(remainingFlux, k, iTracer)
      do iCell = 1, nCells

! mrp Right now I am applying correction to the top layer only.  Need to think
! about if the commented lines below are better.  It distributes the heat on the
! top many layers, but is more computation (all the way down to maxLevelCell)
         tend(iTracer, 1, iCell) = tend(iTracer, 1, iCell) - surfaceHeatFluxCorrection


!        remainingFlux = 1.0_RKIND
!        do k = 1, maxLevelCell(iCell)
!          remainingFlux = remainingFlux - fractionAbsorbed(k, iCell)
!
!            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) - surfaceHeatFluxCorrection * fractionAbsorbed(k, iCell)
!        end do
!
!        if(maxLevelCell(iCell) > 0 .and. remainingFlux > 0.0_RKIND) then
!            tend(iTracer, maxLevelCell(iCell), iCell) = tend(iTracer, maxLevelCell(iCell), iCell) &
!                                                      - surfaceHeatFluxCorrection * remainingFlux
!        end if
      end do
      !$omp end do

      call mpas_timer_stop("correct_surface_heat_flux_tend")

   !--------------------------------------------------------------------

   end subroutine ocn_correct_surface_heat_flux_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_surface_flux_init
!
!> \brief   Initializes ocean tracer surface flux quantities
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \details
!>  This routine initializes quantities related to surface fluxes in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_surface_flux_init(err)!{{{

   !--------------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_disable_tr_sflux
      logical, pointer :: config_correct_surfaceHeatFlux_RainEvap
      logical, pointer :: config_correct_surfaceHeatFlux_Runoff
      logical, pointer :: config_correct_surfaceHeatFlux_SeaIce

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_disable_tr_sflux', config_disable_tr_sflux)
      call mpas_pool_get_config(ocnConfigs,'config_correct_surfaceHeatFlux_RainEvap',config_correct_surfaceHeatFlux_RainEvap)
      call mpas_pool_get_config(ocnConfigs,'config_correct_surfaceHeatFlux_Runoff',config_correct_surfaceHeatFlux_Runoff)
      call mpas_pool_get_config(ocnConfigs,'config_correct_surfaceHeatFlux_SeaIce',config_correct_surfaceHeatFlux_SeaIce)

      surfaceTracerFluxOn = .true.
      correct_surfaceHeatFlux_RainEvap = config_correct_surfaceHeatFlux_RainEvap
      correct_surfaceHeatFlux_Runoff = config_correct_surfaceHeatFlux_Runoff
      correct_surfaceHeatFlux_SeaIce = config_correct_surfaceHeatFlux_SeaIce

      if (config_disable_tr_sflux) then
         surfaceTracerFluxOn = .false.
         correct_surfaceHeatFlux_RainEvap = .false.
         correct_surfaceHeatFlux_Runoff = .false.
         correct_surfaceHeatFlux_SeaIce = .false.
      end if

   end subroutine ocn_tracer_surface_flux_init!}}}

!***********************************************************************

end module ocn_tracer_surface_flux_to_tend

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
