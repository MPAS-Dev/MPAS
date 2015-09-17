! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_mep
!
!> \brief Ocean meso-scale eddy parameterization module 
!> \author Juan A. Saenz, Todd Ringler
!> \date   14 September 2015
!> \details
!>  This module contains routines for computing the meso-scale eddy 
!>  parameterization, as well as interface routines.
!
!-----------------------------------------------------------------------

module ocn_mep

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

   public :: ocn_mep_update, &
             ocn_mep_init, &
             ...

   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: mepStdGMOn, mepTWAdragGMOn, mepTWAdragVertViscOn, mepRediOn, mepUseGMmodule

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_mep_init
!
!> \brief  Initialize the mep module
!> \author Juan A. Saenz, Todd Ringler
!> \date   14 September 2015
!> \details 
!> This routine initializes control flags, parameters and quantities related
!> to the meso-scale eddy parameterization.

!
!-----------------------------------------------------------------------

   subroutine ocn_mep_init(err)!{{{

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

      character (len=StrKIND), pointer :: config_eos_type

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_mep_type', config_mep_type)

      mepStdGMOn = .false.
      mepTWAdragGMOn = .false.
      mepTWAdragVertViscOn = .false.
      mepRediOn = .false.
      mepUseGMmodule = .false.

      ! Momentum eddy parameterization flags
      if (config_mep_stdGM_enable) then
         mepStdGMOn = .true.
         mepUseGMmodule = .true.
      end if
      if (config_mep_twaDragGM_enable) then
         mepTWAdragGMOn = .true.
         mepUseGMmodule= .true.
      end if
      if (config_mep_twaDragVertVisc_enable) then
         mepTWAdragVertViscOn = .true.
         call ocn_mep_twa_drag_VertVisc_init(err)
      end if
      
      ! Tracer eddy parameterization flags
      if (config_mep_Redi_enable) then
         mepRediOn = .true.
         mepUseGMmodule = .true.
      end if

      ! Check for invalid parameterizations being enabled
      if (mepStdGMOn .and. (mepTWAdragGMOn .or. mepTWAdragVertViscOn)) then
         write (stderrUnit,*) 'Invalid combination of meso-scale eddy parameterizations:'
         write (stderrUnit,*) '  stdGM is incompatible with twaDragGM and with twaDragVertVisc.'
         err = 1
      end if
      if (mepTWAdragGMOn .and. mepTWAdragVertViscOn)
         write (stderrUnit,*) 'Invalid combination of meso-scale eddy parameterizations:'
         write (stderrUnit,*) '  twaDragGM is incompatible with twaDragVertVisc.'
         err = 1
      end if

      ! Initialize GM module if needed
      if (mepUseGMmodule) then
         call ocn_gm_init(err)
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_mep_init!}}}


!***********************************************************************
!
!  routine ocn_mep_update
!
!> \brief  This routine updates the mep module given an ocean state 
!> \author Juan A. Saenz, Todd Ringler
!> \date   17 September 2015
!> \details 
!>  This routine updates the mep given a current ocean state.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_update(diagnosticsPool, meshPool, scratchPool)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: diagnosticsPool !< Input/Output: Diagnostics information
      type (mpas_pool_type), intent(inout) :: scratchPool !< Input/Output: Scratch structure

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      
      if (mepUseGMmodule) then
         ! might want to give the routine called here a more general name
         call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
      end if
   

   end subroutine ocn_mep_update!}}}


!***********************************************************************
!
!  routine ocn_mep_transportVelocity_update
!
!> \brief  This routine updates the transport velocity 
!> \author Juan A. Saenz, Todd Ringler
!> \date   17 September 2015
!> \details 
!>  This routine receives a velocity (eulerian velocity) and adds 
!>  the GM bolus velocity, resulting in the residual velocity
!>  It is assumed that the mep module has already been updated.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_transportVelocity_update(normalTransportVelocity)!{{{

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

      real (kind=RKIND), dimension(:,:), intent(inout) :: normalTransportVelocity

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      ! NOTE: This routine is called within block loops
      
      if (mepStdGMOn) then
         normalTransportVelocity(:,:) = normalTransportVelocity(:,:) + normalGMBolusVelocity(:,:)
      end if
   

   end subroutine ocn_mep_transportVelocity_update!}}}


!***********************************************************************
!
!  routine ocn_mep_reconstruct
!
!> \brief  This routine interfaces the cell center reconstruction
!> \author Juan A. Saenz, Todd Ringler
!> \date   17 September 2015
!> \details 
!>  This routine interfaces the cell center reconstruction for diagnostics
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_reconstruct(diagnosticsPool, meshPool)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostic information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      ! NOTE: This routine is called within block loops
      
      if (mepUseGMmodule) then
         call ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)
      end if
   

   end subroutine ocn_mep_transportVelocity_update!}}}




!***********************************************************************

end module ocn_mep

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
