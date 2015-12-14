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
             ocn_mep_vmix_mom_update_vertIntEddyE, &
             ocn_mep_vert_visc_init, &
             ocn_mep_vmix_mom_build

   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: mepVertViscConstOn, mepVertViscProgOn


   real (kind=RKIND), pointer :: config_mep_climate_Tc, config_mep_climate_dt
   real (kind=RKIND), pointer :: config_mep_vert_visc_const_value

   

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
      integer, dimension(:), pointer   :: maxLevelCell
      real(kind=RKIND), dimension(:,:), pointer :: BruntVaisalaFreqTop


      call mpas_pool_get_array(meshPool, 'maxLevelCell',  maxLevelCell)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
      call mpas_pool_get_array(mepPool, 'NsqClimate', NsqClimate)


      do iCell = 1, nCells
         do kLevel = 1, maxLevelCell(iCell)
            NsqClimate(k, iCell) = NsqClimate(kLevel, iCell) &
               + config_mep_climate_dt/config_mep_climate_Tc &
               * ( BruntVaisalaFreqTop(kLevel, iCell) - NsqClimate(kLevel, iCell) )
         end do
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_mom_update_climate!}}}



!***********************************************************************
!
!  routine ocn_mep_vert_visc_init
!
!> \brief   Initialize mep variables related to vertical mixing of momentum
!> \author  Juan A. Saenz, Todd Ringler
!> \date    December 2015
!> \version 
!> \details 
!>  This routine initializes variables used in this module
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vert_visc_init(err)!{{{

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

      character (len=StrKIND), pointer :: config_mep_formulation

      call mpas_pool_get_config(ocnConfigs, 'config_mep_formulation', config_mep_formulation)

      err = 0

      mepVertViscConstOn = .false.
      mepVertViscProgOn = .false.

      config_mep_climate_Tc = 0.0
      config_mep_climate_dt = 0.0

      call mpas_pool_get_config(ocnConfigs, 'config_mep_formulation', config_mep_formulation)
      if ( trim(config_mep_formulation) == 'constantVerticalViscocity' ) then
         mepVertViscConstOn = .true.
         call mpas_pool_get_config(ocnConfigs, 'config_mep_formulation', config_mep_formulation)
      else if ( trim(config_mep_formulation) == 'prognosticVerticalViscocity' ) then
         mepVertViscProgOn = .true.
         call mpas_pool_get_config(ocnConfigs, 'config_mep_climate_Tc', config_mep_climate_Tc)
         call mpas_pool_get_config(ocnConfigs, 'config_mep_climate_dt', config_mep_climate_dt)
      end if


   !--------------------------------------------------------------------

   end subroutine ocn_mep_vert_visc_init!}}}



!***********************************************************************
!
!  routine ocn_mep_vmix_mom_build
!
!> \brief   Add momentum vertical mixing coefficient to vertDiffTopOfCell
!> \author  Juan A. Saenz, Todd Ringler
!> \date    November 2015
!> \version 
!> \details 
!>  This routine adds constant momentum vertical mixing coefficient to vertDiffTopOfCell
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vmix_mom_build(vertDiffTopOfCell, err)!{{{

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

      real (kind=RKIND), dimension(:,:), intent(inout) :: vertDiffTopOfCell !< Input/Output: Vertical diffusion

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
      elseif (mepVertViscConstOn) then
         vertDiffTopOfCell = vertDiffTopOfCell + config_mep_vert_visc_const_value
      else 

      end if


   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_mom_build!}}}


!***********************************************************************

end module ocn_mep_vmix_mom

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
