! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_mep_vmix_mom
!
!> \brief Ocean meso-scale eddy parameterization of vertical viscocity module
!> \author Juan A. Saenz, Todd Ringler
!> \date   November 2015
!> \details
!>  This module contains routines for computing the meso-scale eddy
!>  parameterization of vertical viscocity, as well as interface routines.
!
!-----------------------------------------------------------------------

module ocn_mep_vmix_mom

   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants


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

   public :: ocn_mep_vmix_mom_update_climate, &
             ocn_mep_vmix_solve_evp, &
             ocn_mep_vmix_time_integrate_vertInt_eddy_energy, &
             ocn_mep_vmix_eddyEnergy, &
             ocn_mep_vmix_vert_viscocity, &
             ocn_mep_vmix_mom_init, &
             ocn_mep_vmix_mom_build, &
             ocn_mep_vmix_mom_calculate_eddy_force


   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: mepVertViscConstOn, mepVertViscProgOn


   character (len=StrKIND), pointer :: config_mep_formulation
   real (kind=RKIND), pointer :: config_mep_climate_Tc, config_mep_climate_dt
   real (kind=RKIND), pointer :: config_mep_climate_eddy_energy_dt
   real (kind=RKIND), pointer :: config_mep_eigenvector_A, config_mep_eigenvector_Le
   real (kind=RKIND), pointer :: config_mep_eigenvector_phi0
   real (kind=RKIND), pointer :: config_mep_vert_visc_const
   real (kind=RKIND), pointer :: config_mep_energy_dissipation_const


!***********************************************************************

contains


!***********************************************************************
!
!  routine ocn_mep_vmix_mom_update_climate
!
!> \brief   Update the climate of the flow.
!> \author  Juan A. Saenz, Todd Ringler
!> \date    December 2015
!> \version
!> \details
!>  This routine updates the climate of the flow, quantities used by the
!>  prognostic mep.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vmix_mom_update_climate(statePool, meshPool, diagnosticsPool, mepPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: diagnosticsPool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: mepPool structure

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
      integer :: iCell, kLevel, iEdge
      real (kind=RKIND) :: uDivEPFT
      real (kind=RKIND), dimension(:,:), allocatable :: dFielddzTopEdges, workVariable

      integer, pointer :: nCells, nEdges, nVertLevels
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeBot, maxLevelEdgeTop
      real(kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge
      real(kind=RKIND), dimension(:,:), pointer :: NsqFiltered, dEPFT31dz
      real(kind=RKIND), dimension(:,:), pointer :: uDivEPFTClimateFilter
      real(kind=RKIND), dimension(:,:), pointer :: vertViscTopOfEdgeMEP
      real(kind=RKIND), dimension(:,:), pointer :: BruntVaisalaFreqTop, normalVelocity

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_array(meshPool, 'maxLevelCell',  maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot',  maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
      call mpas_pool_get_array(mepPool, 'NsqFiltered', NsqFiltered)
      call mpas_pool_get_array(mepPool, 'dEPFT31dz', dEPFT31dz)
      call mpas_pool_get_array(mepPool, 'uDivEPFTClimateFilter', uDivEPFTClimateFilter)
      call mpas_pool_get_array(mepPool, 'vertViscTopOfEdgeMEP', vertViscTopOfEdgeMEP)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity)

      allocate( dFielddzTopEdges(nVertLevels+1, nEdges) )
      allocate( workVariable(nVertLevels+1, nEdges) )

      do iCell = 1, nCells
         do kLevel = 1, maxLevelCell(iCell)
            ! NsqFiltered lives on top of cell centers
            NsqFiltered(kLevel, iCell) = &
                 NsqFiltered(kLevel, iCell) &
               + config_mep_climate_dt/config_mep_climate_Tc &
               * ( BruntVaisalaFreqTop(kLevel, iCell) - NsqFiltered(kLevel, iCell) )
         end do
      end do

      ! EPFT is rotated horizontally onto the normal, tangential unit vectors.
      ! EPFT(:,1) (first column) corresponds to normalVelocity component of momentum.
      ! Since we are not parameterizing the Reynolds stress terms
      ! EPFT(1,1) and EPFT(2,1), the only contribution to Div EPFT is
      ! the vertical derivative of EPFT(3,1)

      call computeVerticalDerivativeEdgeMidToTop( nEdges, maxLevelEdgeBot, &
         layerThicknessEdge, normalVelocity, dFielddzTopEdges )

      workVariable = - vertViscTopOfEdgeMEP * dFielddzTopEdges

      call computeVerticalDerivativeEdgeTopToMid( nEdges, maxLevelEdgeBot, &
         layerThicknessEdge, dFielddzTopEdges, dEPFT31dz)

      do iEdge = 1, nEdges
         do kLevel = 1, maxLevelEdgeTop(iEdge)
            uDivEPFT = normalVelocity(kLevel, iEdge) * dEPFT31dz(kLevel, iEdge)
            uDivEPFTClimateFilter(kLevel, iEdge) = &
               + config_mep_climate_dt/config_mep_climate_Tc &
               * ( uDivEPFT - uDivEPFTClimateFilter(kLevel, iEdge) )
         end do
      end do

      deallocate( dFielddzTopEdges )
      deallocate( workVariable )


   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_mom_update_climate!}}}



!***********************************************************************
!
!  routine ocn_mep_vmix_solve_evp
!
!> \brief   Assign an analytical expression to the eigenvalue problem
!> \author  Juan A. Saenz, Todd Ringler
!> \date    December 2015
!> \version
!> \details
!>  This routine assignes an approximate analytical expression to the
!   eigenvalue problem with the new climate.
!>
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vmix_solve_evp(statePool, meshPool, diagnosticsPool, mepPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: diagnosticsPool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: mepPool structure

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
      integer :: iCell, kLevel, cell1, cell2, iEdge
      integer, pointer :: nCells, nEdges
      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge
      real(kind=RKIND) :: NsqFilteredTopEdge
      real(kind=RKIND), dimension(:), pointer :: ssh, Leddy, fEdge, maxLevelEdgeTop
      real(kind=RKIND), dimension(:,:), pointer :: phi, NsqFiltered, alphaEff
      real(kind=RKIND), dimension(:,:), pointer :: zMid, layerThicknessEdge

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_array(meshPool, 'maxLevelCell',  maxLevelCell)
      call mpas_pool_get_array(meshPool, 'fEdge',  fEdge)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge',  cellsOnEdge)
      call mpas_pool_get_array(mepPool, 'eigenvectorPhi', phi)
      call mpas_pool_get_array(mepPool, 'Leddy', Leddy)
      call mpas_pool_get_array(mepPool, 'alphaEff', alphaEff)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge',  layerThicknessEdge)
      call mpas_pool_get_array(statePool, 'ssh', ssh)

      ! eigenvector
      do iCell = 1, nCells
         do kLevel = 1, maxLevelCell(iCell)
            phi(kLevel, iCell) = &
               - config_mep_eigenvector_A &
               * exp( (zMid(kLevel, iCell)-ssh(iCell) ) / config_mep_eigenvector_Le ) &
               + config_mep_eigenvector_phi0
         end do
      end do

      ! efficiency of eddy mixing is set to 1.0 for now
      alphaEff = 1.0_RKIND

      ! eddy length scale is calculated using the Wentzel-Kramers-Brillouin (WKB)
      ! approximation of the Rossby Radius of deformation
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         Leddy(iEdge) = 0.0_RKIND
         do kLevel = 1, maxLevelEdgeTop(iEdge)
            NsqFilteredTopEdge = 0.5_RKIND * ( NsqFiltered(kLevel,cell1) + NsqFiltered(kLevel,cell2) )
            Leddy(iEdge) = Leddy(iEdge) + sqrt(NsqFilteredTopEdge) * layerThicknessEdge(kLevel, iEdge)
         end do
         Leddy(iEdge) =  Leddy(iEdge) / fEdge(iEdge)
      end do


   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_solve_evp!}}}



!***********************************************************************
!
!  routine ocn_mep_vmix_time_integrate_vertInt_eddy_energy
!
!> \brief   Update the vertically integrated eddy energy to a new timestep
!> \author  Juan A. Saenz, Todd Ringler
!> \date    December 2015
!> \version
!> \details
!>  This routine updates the vertically integrated eddy energy by integrating in time
!>  cell.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vmix_time_integrate_vertInt_eddy_energy(meshPool, diagnosticsPool, mepPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: diagnosticsPool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: mepPool structure

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
      integer :: iCell
      integer, pointer :: nCells
      real(kind=RKIND), dimension(:), pointer :: vertIntEddyEnergy_tend
      real(kind=RKIND), dimension(:), pointer :: vertIntEddyEnergy

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(mepPool, 'vertIntEddyEnergy_tend', vertIntEddyEnergy_tend)
      call mpas_pool_get_array(mepPool, 'vertIntEddyEnergy', vertIntEddyEnergy)

      ! update the time tendency of vertically integrated eddy energy
      call ocn_mep_vmix_mom_vertIntEddyEnergy_tend(mepPool, meshPool, diagnosticsPool, err)

      ! Forward Euler time integration
      do iCell = 1, nCells
         vertIntEddyEnergy(iCell) = vertIntEddyEnergy(iCell) &
            + vertIntEddyEnergy_tend(iCell) * config_mep_climate_eddy_energy_dt
      end do


   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_time_integrate_vertInt_eddy_energy!}}}




!***********************************************************************
!
!  routine ocn_mep_vmix_eddyEnergy
!
!> \brief   Calculate the vertical profile of eddy energy at each cell
!> \author  Juan A. Saenz, Todd Ringler
!> \date    December 2015
!> \version
!> \details
!>  This routine calculates the vertical profile of eddy energy at each
!>  cell.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vmix_eddyEnergy(statePool, meshPool, mepPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: meshPool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: mepPool

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
      integer :: iCell, kLevel
      integer, pointer :: nCells
      integer, dimension(:), pointer :: maxLevelCell
      real(kind=RKIND) :: intPhiSq
      real(kind=RKIND), dimension(:), pointer :: vertIntEddyEnergy
      real(kind=RKIND), dimension(:,:), pointer :: layerThickness, phi
      real(kind=RKIND), dimension(:,:), pointer :: eddyEnergy

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'maxLevelCell',  maxLevelCell)
      call mpas_pool_get_array(statePool, 'layerThickness',  layerThickness)
      call mpas_pool_get_array(mepPool, 'eigenvectorPhi', phi)
      call mpas_pool_get_array(mepPool, 'eddyEnergy', eddyEnergy)
      call mpas_pool_get_array(mepPool, 'vertIntEddyEnergy', vertIntEddyEnergy)

      do iCell = 1, nCells
         intPhiSq = 0.0
         do kLevel = 1, maxLevelCell(iCell)
            intPhiSq = intPhiSq + phi(kLevel,iCell)**2 * layerThickness(kLevel,iCell)
         end do
         do kLevel = 1, maxLevelCell(iCell)
            eddyEnergy(kLevel, iCell) = vertIntEddyEnergy(iCell) * phi(kLevel,iCell) / intPhiSq
         end do
      end do


   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_eddyEnergy!}}}


!***********************************************************************
!
!  routine ocn_mep_vmix_vert_viscocity
!
!> \brief   Calculate the vertical viscocity from mep
!> \author  Juan A. Saenz, Todd Ringler
!> \date    November 2015
!> \version
!> \details
!>  This routine calculates vertical viscocity from the meso-scale eddy
!>  parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vmix_vert_viscocity(meshPool, mepPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: mepPool

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
      integer :: cell1, cell2, iEdge, kLevel

      integer, pointer :: nEdges
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:),   pointer :: fEdge

      real (kind=RKIND) :: NsqFilteredTopEdge
      real (kind=RKIND), dimension(:),   pointer :: Leddy
      real(kind=RKIND), dimension(:,:), pointer :: alphaEff, eddyEnergy
      real(kind=RKIND), dimension(:,:), pointer :: NsqFiltered, vertViscTopOfEdgeMEP

      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge',  cellsOnEdge)

      call mpas_pool_get_array(mepPool, 'Leddy', Leddy)
      call mpas_pool_get_array(mepPool, 'alphaEff', alphaEff)
      call mpas_pool_get_array(mepPool, 'eddyEnergy', eddyEnergy)
      call mpas_pool_get_array(mepPool, 'NsqFiltered', NsqFiltered)
      call mpas_pool_get_array(mepPool, 'vertViscTopOfEdgeMEP', vertViscTopOfEdgeMEP)


      err = 0

      ! Might not need to check if both are on because this is guaranteed in ocn_mep_init
      if (.not. mepVertViscConstOn .and. .not. mepVertViscProgOn) then
         return

      ! Constant vertical viscosity
      else if (mepVertViscConstOn) then
         vertViscTopOfEdgeMEP = config_mep_vert_visc_const

      ! Prognostic eddy closure
      else if (mepVertViscProgOn) then
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do kLevel = 1 , maxLevelEdgeTop(iEdge)
               NsqFilteredTopEdge = 0.5_RKIND * (NsqFiltered(kLevel,cell1) + NsqFiltered(kLevel,cell2))
               NsqFilteredTopEdge = max(NsqFilteredTopEdge, 0.0_RKIND)
               vertViscTopOfEdgeMEP(kLevel, iEdge) = &
                  fEdge(iEdge)**2 / NsqFilteredTopEdge &
                  * sqrt( eddyEnergy(kLevel, iEdge) ) * alphaEff(kLevel, iEdge) * Leddy(iEdge)
            end do
         end do
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_vert_viscocity!}}}


!***********************************************************************
!
!  routine ocn_mep_vmix_mom_init
!
!> \brief   Initialize mep variables related to vertical mixing of momentum
!> \author  Juan A. Saenz, Todd Ringler
!> \date    December 2015
!> \version
!> \details
!>  This routine initializes variables used in this module
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vmix_mom_init(err)!{{{

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
      character (len=StrKIND), pointer :: config_mep_vert_visc_formulation

      call mpas_pool_get_config(ocnConfigs, 'config_mep_vert_visc_formulation', config_mep_vert_visc_formulation)

      err = 0

      mepVertViscConstOn = .false.
      mepVertViscProgOn = .false.

      config_mep_climate_Tc = 0.0_RKIND
      config_mep_climate_dt = 0.0_RKIND
      config_mep_climate_eddy_energy_dt = 0.0_RKIND
      config_mep_eigenvector_A = 0.0_RKIND
      config_mep_eigenvector_Le = 0.0_RKIND
      config_mep_eigenvector_phi0 = 0.0_RKIND
      config_mep_vert_visc_const = 0.0_RKIND
      config_mep_energy_dissipation_const = 0.0_RKIND

      if ( trim(config_mep_formulation) == 'constantVerticalViscocity' ) then
         mepVertViscConstOn = .true.
         call mpas_pool_get_config(ocnConfigs, 'config_mep_vert_visc_const', config_mep_vert_visc_const)
      else if ( trim(config_mep_formulation) == 'prognosticVerticalViscocity' ) then
         mepVertViscProgOn = .true.
         call mpas_pool_get_config(ocnConfigs, 'config_mep_climate_Tc', config_mep_climate_Tc)
         call mpas_pool_get_config(ocnConfigs, 'config_mep_climate_dt', config_mep_climate_dt)
         call mpas_pool_get_config(ocnConfigs, 'config_mep_climate_eddy_energy_dt', config_mep_climate_eddy_energy_dt)
         call mpas_pool_get_config(ocnConfigs, 'config_mep_eigenvector_A', config_mep_eigenvector_A)
         call mpas_pool_get_config(ocnConfigs, 'config_mep_eigenvector_Le', config_mep_eigenvector_Le)
         call mpas_pool_get_config(ocnConfigs, 'config_mep_eigenvector_phi0', config_mep_eigenvector_phi0)
         call mpas_pool_get_config(ocnConfigs, 'config_mep_energy_dissipation_const', config_mep_energy_dissipation_const)
         
      end if


   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_mom_init!}}}


!***********************************************************************
!
!  routine ocn_mep_vmix_mom_build
!
!> \brief   Add mep vertical viscosity coefficient to vertViscTopOfEdge
!> \author  Juan A. Saenz, Todd Ringler
!> \date    November 2015
!> \version
!> \details
!>  This routine adds the vertical viscosity coefficient from meso-scale
!>  eddy parameterization to vertViscTopOfEdge
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vmix_mom_build(vertViscTopOfEdge, mepPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: mepPool !< Input/Output: mepPool structure

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: vertViscTopOfEdge !< Input/Output: Vertical viscosity

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
      real(kind=RKIND), dimension(:,:), pointer :: vertViscTopOfEdgeMEP

      call mpas_pool_get_array(mepPool, 'vertViscTopOfEdgeMEP', vertViscTopOfEdgeMEP)

      err = 0

      ! I don't think I need to check if both are on because this is guaranteed in ocn_mep_init
      if (.not. mepVertViscConstOn .and. .not. mepVertViscProgOn) then
         return
      else
         vertViscTopOfEdge = vertViscTopOfEdge + vertViscTopOfEdgeMEP
      end if


   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_mom_build!}}}


!***********************************************************************
!
!  routine ocn_mep_vmix_mom_vertIntEddyEnergy_tend
!
!> \brief   Compute prognostic vertically integrated eddy energy time tendency
!> \author  Juan A. Saenz, Todd Ringler
!> \date    December 2015
!> \version
!> \details
!>  This routine computes prognostic vertically integrated eddy energy
!>  time tendency.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vmix_mom_vertIntEddyEnergy_tend(mepPool, meshPool, diagnosticsPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: diagnosticsPool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: mepPool structure

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
      integer :: iCell, kLevel, iEdge
      integer, pointer :: nCells
      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
      real(kind=RKIND), dimension(:), pointer :: vertIntEddyEnergy_tend
      real(kind=RKIND), dimension(:), pointer :: vertIntEddyEnergy
      real(kind=RKIND), dimension(:,:), pointer :: uDivEPFTClimateFilter
      real(kind=RKIND), dimension(:,:), pointer :: uSlowTimeFilter
      real(kind=RKIND), dimension(:,:), pointer :: divEPFT31SlowTimeFilter
      real(kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge
      real(kind=RKIND), dimension(:,:), pointer :: edgeMask

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'maxLevelCell',  maxLevelCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell',  nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeMask',  edgeMask)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge',  layerThicknessEdge)

      call mpas_pool_get_array(mepPool, 'vertIntEddyEnergy_tend', vertIntEddyEnergy_tend)
      call mpas_pool_get_array(mepPool, 'vertIntEddyEnergy', vertIntEddyEnergy)
      call mpas_pool_get_array(mepPool, 'uDivEPFTClimateFilter', uDivEPFTClimateFilter)
      call mpas_pool_get_array(mepPool, 'uSlowTimeFilter', uSlowTimeFilter)
      call mpas_pool_get_array(mepPool, 'divEPFT31SlowTimeFilter', divEPFT31SlowTimeFilter)

      err = 0

      vertIntEddyEnergy_tend = 0.0

      ! juan issue : using unfiltered layerThicknessEdge

      ! Time tendency from eddy mean flow interactions
      ! Investigating whether to use ubarDivEPFTFiltered or uSlowTimeFilter*Div(viscosityBar * dubar/dz).
      ! Will implement the latter while the former will be commented out for the time being.
      do iCell = 1, nCells
         do iEdge = 1, nEdgesOnCell(iCell)
            do kLevel = 1, maxLevelCell(iCell)
               vertIntEddyEnergy_tend(iCell) = vertIntEddyEnergy_tend(iCell) &
               !    - uDivEPFTClimateFilter(kLevel,iEdge)*layerThicknessEdge(kLevel,iEdge)*edgeMask(kLevel,iEdge)
                  - uSlowTimeFilter(kLevel,iEdge)*divEPFT31SlowTimeFilter(kLevel, iEdge)*layerThicknessEdge(kLevel,iEdge)*edgeMask(kLevel,iEdge)
            end do
         end do
      end do


      ! vertically integrated dissipation
      do iCell = 1, nCells
         vertIntEddyEnergy_tend(iCell) = vertIntEddyEnergy_tend(iCell) &
            + config_mep_energy_dissipation_const*vertIntEddyEnergy(iCell)
      end do

      ! These terms are zero for now:
      ! vertically integrated advection of eddy energy
      ! vertically integrated transport

   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_mom_vertIntEddyEnergy_tend!}}}


!***********************************************************************
!
!  routine ocn_mep_vmix_mom_calculate_eddy_force
!
!> \brief   Update force from eddy viscisity
!> \author  Juan A. Saenz
!> \date    March, 2016
!> \version
!> \details
!>  This routine computes the eddy force using the mep vertical viscisity
!>  and the updated slow time filtered normal velocity.
!
!-----------------------------------------------------------------------
   subroutine ocn_mep_vmix_mom_calculate_eddy_force(meshPool, diagnosticsPool, mepPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: diagnostics information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: mepPool structure

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
      integer, pointer :: nEdges, nVertLevels
      integer, dimension(:), pointer :: maxLevelEdgeBot
      real(kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge
      real(kind=RKIND), dimension(:,:), pointer :: uSlowTimeFilter
      real(kind=RKIND), dimension(:,:), pointer :: vertViscTopOfEdgeMEP, eddyViscosityForceMEP
      real (kind=RKIND), dimension(:,:), allocatable :: duSlowTimeFilterdz, workVariableTopEdges

      err = 0      

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot',  maxLevelEdgeBot)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(mepPool, 'uSlowTimeFilter', uSlowTimeFilter)
      call mpas_pool_get_array(mepPool, 'vertViscTopOfEdgeMEP', vertViscTopOfEdgeMEP)
      call mpas_pool_get_array(mepPool, 'eddyViscosityForceMEP', eddyViscosityForceMEP)

      allocate( duSlowTimeFilterdz(nVertLevels+1, nEdges) )
      allocate( workVariableTopEdges(nVertLevels+1, nEdges) )

      ! calculate the vertical derivative of uSlowTimeFilter: duSlowTimeFilterdz
      call computeVerticalDerivativeEdgeMidToTop( &
         nEdges, maxLevelEdgeBot, layerThicknessEdge, uSlowTimeFilter, duSlowTimeFilterdz &
         )

      workVariableTopEdges = vertViscTopOfEdgeMEP * duSlowTimeFilterdz
      
      ! calculate eddyViscosityForceMEP = - ddz(viscosity * duSlowTimeFilterdz)
      call computeVerticalDerivativeEdgeTopToMid( &
         nEdges, maxLevelEdgeBot, layerThicknessEdge, eddyViscosityForceMEP, workVariableTopEdges &
         )
      eddyViscosityForceMEP = - eddyViscosityForceMEP

      deallocate( duSlowTimeFilterdz )
      deallocate( workVariableTopEdges )

   end subroutine ocn_mep_vmix_mom_calculate_eddy_force!}}}



!***********************************************************************
!
!  routine computeVerticalDerivativeEdgeMidToTop
!
!> \brief   Compute vertical derivative of field on mid edges
!> \author  Juan A. Saenz
!> \date    January 2016
!> \version
!> \details
!>  This routine computes the vertical derivative of a field that exists
!>  in the middle of edges. The resulting derivative exists at top of Edges.
!
!-----------------------------------------------------------------------

   subroutine computeVerticalDerivativeEdgeMidToTop( &
      nEdges, maxLevelEdgeBot, layerThicknessEdge, fieldMidEdge, dFielddzTopEdges &
      )!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      integer, intent(in) :: nEdges
      integer, dimension(nEdges), intent(in) :: maxLevelEdgeBot
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThicknessEdge
      real (kind=RKIND), dimension(:,:), intent(in) :: fieldMidEdge

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

      real (kind=RKIND), dimension(:,:), intent(out) :: dFielddzTopEdges

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      integer :: iEdge, kLevel
      real (kind=RKIND) :: wrkAbove, wrkBelow, dz

      dFielddzTopEdges = 0.0_RKIND

      do iEdge = 1, nEdges
         ! maxLevelEdgeBot is the maximum (deepest) of the surrounding cells
         do kLevel = 2, maxLevelEdgeBot(iEdge)
            wrkAbove = fieldMidEdge(kLevel-1, iEdge)
            wrkBelow = fieldMidEdge(kLevel, iEdge)
            dz = 0.5_RKIND * ( layerThicknessEdge(kLevel-1, iEdge) + layerThicknessEdge(kLevel, iEdge) )
            dFielddzTopEdges(kLevel,iEdge) = (wrkAbove - wrkBelow) / dz
         end do
         ! juan issue : can make derivative at top of top edge a function of wind stress?
         dFielddzTopEdges(1,iEdge) = dFielddzTopEdges(2,iEdge)
         ! juan issue : can make derivative at bottom of bottom edge a function of bottom drag?
         dFielddzTopEdges(maxLevelEdgeBot(iEdge)+1, iEdge) = dFielddzTopEdges(maxLevelEdgeBot(iEdge), iEdge)
      end do

   !--------------------------------------------------------------------

   end subroutine computeVerticalDerivativeEdgeMidToTop!}}}



!***********************************************************************
!
!  routine computeVerticalDerivativeEdgeTopToMid
!
!> \brief   Compute vertical derivative of field on top edges.
!> \author  Juan A. Saenz
!> \date    January 2016
!> \version
!> \details
!>  This routine computes the vertical derivative of a field that exists on
!>  top of edges. The resulting derivative exists in the middle of Edges.
!
!-----------------------------------------------------------------------

   subroutine computeVerticalDerivativeEdgeTopToMid( &
      nEdges, maxLevelEdgeBot, layerThicknessEdge, fieldTopEdge, dFielddzMidEdges )!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      integer, intent(in) :: nEdges
      integer, dimension(nEdges), intent(in) :: maxLevelEdgeBot
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThicknessEdge
      real (kind=RKIND), dimension(:,:), intent(in) :: fieldTopEdge

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

      real (kind=RKIND), dimension(:,:), intent(out) :: dFielddzMidEdges

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      integer :: iEdge, kLevel
      real (kind=RKIND) :: wrkAbove, wrkBelow

      dFielddzMidEdges = 0.0_RKIND

      do iEdge = 1, nEdges
         ! maxLevelEdgeBot is the maximum (deepest) of the surrounding cells
         do kLevel = 1, maxLevelEdgeBot(iEdge)
            wrkAbove = fieldTopEdge(kLevel, iEdge)
            wrkBelow = fieldTopEdge(kLevel+1, iEdge)
            dFielddzMidEdges(kLevel,iEdge) = (wrkAbove - wrkBelow) / layerThicknessEdge(kLevel, iEdge)
         end do
      end do

   !--------------------------------------------------------------------

   end subroutine computeVerticalDerivativeEdgeTopToMid!}}}


   

!***********************************************************************

end module ocn_mep_vmix_mom

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
