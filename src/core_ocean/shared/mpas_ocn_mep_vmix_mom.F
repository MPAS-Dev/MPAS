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
             ocn_mep_vmix_mom_build

   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: mepVertViscConstOn, mepVertViscProgOn

   character (len=StrKIND), pointer :: config_mep_formulation
   real (kind=RKIND), pointer :: config_mep_climate_Tc, config_mep_climate_dt
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

   subroutine ocn_mep_vmix_mom_update_climate(meshPool, diagnosticsPool, mepPool, err)!{{{

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
      integer :: iCell, kLevel
      integer, pointer :: nCells
      integer, dimension(:), pointer :: maxLevelCell
      real(kind=RKIND), dimension(:,:), pointer :: NsqFiltered
      real(kind=RKIND), dimension(:,:), pointer :: BruntVaisalaFreqTop


      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'maxLevelCell',  maxLevelCell)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
      call mpas_pool_get_array(mepPool, 'NsqFiltered', NsqFiltered)


      do iCell = 1, nCells
         do kLevel = 1, maxLevelCell(iCell)
            ! NsqFiltered lives on top of cell centers
            NsqFiltered(k, iCell) = &
                 NsqFiltered(kLevel, iCell) &
               + config_mep_climate_dt/config_mep_climate_Tc &
               * ( BruntVaisalaFreqTop(kLevel, iCell) - NsqFiltered(kLevel, iCell) )
         end do
         dissipationFiltered(iCell) = 
      end do

      do iEdge = 1, nEdges
         do kLevel = 1, maxLevelCell(iCell)
            uDivEPFTFiltered(kLevel, iEdge) = 
         end do
      end do



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

   subroutine ocn_mep_vmix_solve_evp(diagnosticsPool, mepPool, meshPool, statePool, err)!{{{

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
      integer :: iCell, kLevel
      integer, pointer :: nCells
      integer, dimension(:), pointer :: maxLevelCell
      real(kind=RKIND), dimension(:), pointer :: ssh
      real(kind=RKIND), dimension(:,:), pointer :: phi ! eigenvector
      real(kind=RKIND), dimension(:,:), pointer :: zMid



      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'maxLevelCell',  maxLevelCell)
      call mpas_pool_get_array(mepPool, 'eigenvectorPhi', phi)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
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

         update:
            alphaEff
            Leddy


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

   subroutine ocn_mep_vmix_time_integrate_vertInt_eddy_energy(mepPool, meshPool, diagnosticsPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: diagnosticsPool

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

      need to iterate on blocks

      ! start time integration
      do iter = 1, nIter
         ! substep update vertIntEddyEnergy and variables that depend on it
         vertIntEddyEnergy = vertIntEddyEnergy + dt * subStepCoef(iter) * tendOld
         update viscosity which depends on eddy energy
         call ocn_mep_vmix_mom_vertIntEddyEnergy_tend(mepPool, meshPool, diagnosticsPool, err)
         do iCell = 1, nCells
            vertIntEddyEnergyNew = vertIntEddyEnergy(iCell) + dt * coef(iter)*tend
         end do
      end do
      ! end time integration


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

   subroutine ocn_mep_vmix_eddyEnergy(mepPool, meshPool, statePool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

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

   subroutine ocn_mep_vmix_vert_viscocity(vertViscTopOfEdge, mepPool, meshPool, err)!{{{

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

      call mpas_pool_get_config(ocnConfigs, 'config_mep_formulation', config_mep_formulation)

      err = 0

      mepVertViscConstOn = .false.
      mepVertViscProgOn = .false.

      config_mep_climate_Tc = 0.0
      config_mep_climate_dt = 0.0
      config_mep_eigenvector_A = 0.0
      config_mep_eigenvector_Le = 0.0
      config_mep_eigenvector_phi0 = 0.0
      config_mep_vert_visc_const = 0.0
      config_mep_energy_dissipation_const = 0.0

      if ( trim(config_mep_formulation) == 'constantVerticalViscocity' ) then
         mepVertViscConstOn = .true.
         call mpas_pool_get_config(ocnConfigs, 'config_mep_vert_visc_const', config_mep_vert_visc_const)
      else if ( trim(config_mep_formulation) == 'prognosticVerticalViscocity' ) then
         mepVertViscProgOn = .true.
         call mpas_pool_get_config(ocnConfigs, 'config_mep_climate_Tc', config_mep_climate_Tc)
         call mpas_pool_get_config(ocnConfigs, 'config_mep_climate_dt', config_mep_climate_dt)
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

   subroutine ocn_mep_vmix_mom_build(vertViscTopOfEdge, err)!{{{

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

      err = 0

      ! I don't think I need to check if both are on because this is guaranteed in ocn_mep_init
      if (.not. mepVertViscConstOn .and. .not. mepVertViscProgOn) then
         return
      else then
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
      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: diagnosticsPool

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
      real(kind=RKIND), dimension(:,:), pointer :: uDivEPFTFiltered
      real(kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge
      real(kind=RKIND), dimension(:,:), pointer :: edgeMask

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'maxLevelCell',  maxLevelCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell',  nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeMask',  edgeMask)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge',  layerThicknessEdge)

      call mpas_pool_get_array(mepPool, 'vertIntEddyEnergy_tend', vertIntEddyEnergy_tend)
      call mpas_pool_get_array(mepPool, 'vertIntEddyEnergy', vertIntEddyEnergy)
      call mpas_pool_get_array(mepPool, 'uDivEPFTFiltered', uDivEPFTFiltered)

      err = 0

      vertIntEddyEnergy_tend = 0.0

      ! vertically integrated advection of eddy energy
      ! no contribution for now

      ! vertically integrated eddy-mean flow interactions
      ! juan issue : using unfiltered layerThicknessEdge
      do iCell = 1, nCells
         do iEdge = 1, nEdgesOnCell(iCell)
            do kLevel = 1, maxLevelCell(iCell)
               vertIntEddyEnergy_tend(iCell) = vertIntEddyEnergy_tend(iCell) &
                  + uDivEPFTFiltered(kLevel,iEdge)*layerThicknessEdge(kLevel,iEdge)*edgeMask(kLevel,iEdge)
            end do
         end do
      end do

      ! vertically integrated transport
      ! no contribution for now

      ! vertically integrated dissipation
      do iCell = 1, nCells
         vertIntEddyEnergy_tend(iCell) = vertIntEddyEnergy_tend(iCell) &
            + config_mep_energy_dissipation_const*vertIntEddyEnergy
      end do


   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_mom_vertIntEddyEnergy_tend!}}}



!***********************************************************************

end module ocn_mep_vmix_mom

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
