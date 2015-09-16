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

   logical :: mepStdGMOn, mepTWAdragGMOn, mepTWAdragVertViscOn 

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

      if (config_mep_stdGM_enable) then
         mepStdGMOn = .true.
         call ocn_mep_stdGM_init(err)
      end if
      if (config_mep_twaDragGM_enable) then
         mepTWAdragGMOn = .true.
         call ocn_mep_twa_drag_GM_init(err)
      end if
      if (config_mep_twaDragVertVisc_enable) then
         mepTWAdragVertViscOn = .true.
         call ocn_mep_twa_drag_VertVisc_init(err)
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


   !--------------------------------------------------------------------

   end subroutine ocn_mep_init!}}}

!***********************************************************************
!
!  routine ocn_mep_update
!
!> \brief  This routine updates the mep module given an ocean state 
!> \author Juan A. Saenz, Todd Ringler
!> \date   14 September 2015
!> \details 
!>  This routine updates the mep given a current ocean state.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_update(err)!{{{


   integer, intent(out) :: err !< Output: error flag

   !--------------------------------------------------------------------
   !
   ! set some local module variables based on input config choices
   !
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------

   end subroutine ocn_mep_update!}}}

!***********************************************************************

end module ocn_mep

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
