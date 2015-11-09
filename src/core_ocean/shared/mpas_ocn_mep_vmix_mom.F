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

   use ...

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

   public :: ocn_mep_vmix_mom_const_build, &
   	ocn_mep_vmix_mom_GL_build, &
   	ocn_mep_vmix_mom_prognostic_closure_build

   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND), pointer :: config_mep_vert_visc_const


!***********************************************************************

contains

ocn_mep_vmix_const_init
initialize vertConstantMEPVisc with value from namelist record

!***********************************************************************
!
!  routine ocn_mep_vert_visc_const_init
!
!> \brief   Add constant momentum vertical mixing coefficient to vertDiffTopOfCell
!> \author  Juan A. Saenz, Todd Ringler
!> \date    November 2015
!> \version 
!> \details 
!>  This routine adds constant momentum vertical mixing coefficient to vertDiffTopOfCell
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vert_visc_const_init(mepPool, err)!{{{

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

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_mep_const_visc', config_mep_vert_visc_const)

   !--------------------------------------------------------------------

   end subroutine ocn_mep_vert_visc_const_init!}}}



!***********************************************************************
!
!  routine ocn_mep_vmix_mom_const_build
!
!> \brief   Add constant momentum vertical mixing coefficient to vertDiffTopOfCell
!> \author  Juan A. Saenz, Todd Ringler
!> \date    November 2015
!> \version 
!> \details 
!>  This routine adds constant momentum vertical mixing coefficient to vertDiffTopOfCell
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vmix_mom_const_build(vertDiffTopOfCell, err)!{{{

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

      vertDiffTopOfCell = vertDiffTopOfCell + config_mep_vert_visc_const

   !--------------------------------------------------------------------

   end subroutine ocn_mep_vmix_mom_const_build!}}}


!***********************************************************************

end module ocn_mep_vmix_mom

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
