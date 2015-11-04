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
   use ocn_mep_gm
   use ocn_diagnostics
   use mpas_vector_operations

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
             ocn_mep_transportVelocity_update, &
             ocn_mep_vel_tend, &
             ocn_mep_reconstruct, &
             ocn_mep_vel_hmix_tend


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

      logical, pointer :: config_mep_stdGM_enable, config_mep_twaDragGM_enable
      logical, pointer :: config_mep_twaDragVertVisc_enable, config_mep_Redi_enable

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_mep_stdGM_enable', config_mep_stdGM_enable)
      call mpas_pool_get_config(ocnConfigs, 'config_mep_twaDragGM_enable', config_mep_twaDragGM_enable)
      call mpas_pool_get_config(ocnConfigs, 'config_mep_twaDragVertVisc_enable', config_mep_twaDragVertVisc_enable)
      call mpas_pool_get_config(ocnConfigs, 'config_mep_Redi_enable', config_mep_Redi_enable)

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
      if (mepTWAdragGMOn .and. mepTWAdragVertViscOn) then
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

   subroutine ocn_mep_update(diagnosticsPool, meshPool, mepPool, scratchPool)!{{{

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
      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: mepPool structure
      type (mpas_pool_type), intent(inout) :: scratchPool !< Input/Output: Scratch structure

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      
      if (mepUseGMmodule) then
         ! might want to give the routine called here a more general name
         call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, mepPool, scratchPool)
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

   subroutine ocn_mep_transportVelocity_update(normalTransportVelocity, mepPool)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      
      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: mepPool structure

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
      real (kind=RKIND), dimension(:,:), pointer :: normalGMBolusVelocity


      call mpas_pool_get_array(mepPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
     
      
      ! NOTE: This routine is called within block loops
      
      if (mepStdGMOn) then
         normalTransportVelocity(:,:) = normalTransportVelocity(:,:) + normalGMBolusVelocity(:,:)
      end if
   

   end subroutine ocn_mep_transportVelocity_update!}}}


!***********************************************************************
!
!  routine ocn_mep_vel_tend
!
!> \brief  This routine increments normalVelocity tend by the tendency due to eddy forces.
!> \author Juan A. Saenz, Todd Ringler
!> \date   October, 2015
!> \details 
!>  This routine increments normalVelocity tend by the tendency due to eddy forces.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vel_tend(mepPool, tendPool, diagnosticsPool, meshPool, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool        !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostic information

      !jas issue: it might be that only need uBolus or normalBolusGM, not entire diagnostics?

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: mepPool !< Input/Output: Tendency structure
      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency structure
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: normal velocity tendency

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


      if (mepTWAdragGMOn) then
         call ocn_mep_twa_drag_GM(mepPool, tendPool, meshPool, tend, err)
      end if

   end subroutine ocn_mep_vel_tend!}}}


!***********************************************************************
!
!  routine ocn_mep_twa_drag_GM
!
!> \brief  This routine calculates normalVelocity tend due GM bolus coriolis.
!> \author Juan A. Saenz, Todd Ringler
!> \date   November, 2015
!> \details 
!>  This routine calculates normalVelocity tend due coriolis force caused by 
!>  the GM bolus velocity.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_twa_drag_GM(mepPool, tendPool, meshPool, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: mepPool  !< Input: Diagnostic information
      type (mpas_pool_type), intent(in) :: meshPool            !< Input: Mesh information

      !jas issue: it might be that only need uBolus or normalBolusGM, 
      !not entire diagnostics?

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency structure
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend          !< Input/Output: normal velocity tendency

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

      logical  :: includeHalo
      integer :: iEdge, kLevel
      integer, pointer :: nVertLevels, nCells, nEdges

      integer, dimension(:,:), pointer :: edgeMask
      real (kIND=RKIND), dimension(:),   pointer :: fEdge
      real (kind=RKIND), dimension(:,:), pointer :: edgeTangentVectors
      real (kind=RKIND), dimension(:,:), pointer :: tend_vel_mep
      real (kind=RKIND), dimension(:,:), pointer :: GMBolusVelocityX
      real (kind=RKIND), dimension(:,:), pointer :: GMBolusVelocityY
      real (kind=RKIND), dimension(:,:), pointer :: GMBolusVelocityZ
      
      real (kind=RKIND), dimension(:,:),   allocatable :: normalBolusEdge
      real (kind=RKIND), dimension(:,:),   allocatable :: tangentBolusEdge
      real (kind=RKIND), dimension(:,:,:), allocatable :: bolusCell


      err = 0


      call mpas_pool_get_dimension(meshPool, 'nCells', nCells) ! includes halos
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges) ! includes halos
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'edgeTangentVectors', edgeTangentVectors)

      
      allocate(normalBolusEdge(nVertLevels,nEdges))
      allocate(tangentBolusEdge(nVertLevels,nEdges))
      allocate(bolusCell(3,nVertLevels,nCells+1))


      includeHalo = .true.

      ! retrieve velocity tendencry by MEP
      call mpas_pool_get_array(tendPool, 'normalVelocityMEP', tend_vel_mep)
      tend_vel_mep = 0.0


      !retrieve x,y,z bolus velocity at cell centers
      call mpas_pool_get_array(mepPool, 'GMBolusVelocityX', GMBolusVelocityX)
      call mpas_pool_get_array(mepPool, 'GMBolusVelocityY', GMBolusVelocityY)
      call mpas_pool_get_array(mepPool, 'GMBolusVelocityZ', GMBolusVelocityZ)

      bolusCell(1,:,:) = GMBolusVelocityX
      bolusCell(2,:,:) = GMBolusVelocityY
      bolusCell(3,:,:) = GMBolusVelocityZ

      !interpolate bolus velocity at cell centers, mid layers, to Edges
      !make sure to average the x,y,z components, not lat/lon components.
      call mpas_vector_R3Cell_to_2DEdge(bolusCell, meshPool, edgeTangentVectors, &
        includeHalo, normalBolusEdge, tangentBolusEdge)

      !calculate tendency = -dE3i/dz dot ni = - f*uBolus_tangent
      do iEdge = 1, nEdges
        do kLevel = 1,nVertLevels
          tend_vel_mep(kLevel,iEdge) = &
              -fEdge(iEdge) * tangentBolusEdge(kLevel,iEdge) * edgeMask(kLevel, iEdge)
        end do
      end do

      !update tendency with div EPFT parameterized with GM
      tend = tend + tend_vel_mep


      deallocate(normalBolusEdge)
      deallocate(tangentBolusEdge)
      deallocate(bolusCell)
    
    end subroutine ocn_mep_twa_drag_GM!}}}


!***********************************************************************
!
!  routine ocn_mep_vel_hmix_tend
!
!> \brief  This routine increments the horizontal mixing tendency for momentum due to mep.
!> \author Juan A. Saenz, Todd Ringler
!> \date   November, 2015
!> \details 
!>  This routine increments the horizontal mixing tendency for momentum due to 
!>  meso scale eddies.
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_vel_hmix_tend(meshPool, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool        !< Input: Mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: tend  !< Input/Output: normal velocity tendency

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


      return

   end subroutine ocn_mep_vel_hmix_tend!}}}


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

   subroutine ocn_mep_reconstruct(diagnosticsPool, meshPool, mepPool)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(in) :: mepPool !< Input: mesh information
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
         call ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool, mepPool)
      end if
   

   end subroutine ocn_mep_reconstruct!}}}


!***********************************************************************
!
!  routine ocn_mep_twa_drag_VertVisc_init
!
!> \brief  This routine interfaces the cell center reconstruction
!> \author Juan A. Saenz, Todd Ringler
!> \date   17 September 2015
!> \details 
!>  This routine interfaces the cell center reconstruction for diagnostics
!
!-----------------------------------------------------------------------

   subroutine ocn_mep_twa_drag_VertVisc_init(err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      
      integer, intent(out) :: err !< Output: error flag

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

      err = 0 

   end subroutine ocn_mep_twa_drag_VertVisc_init!}}}


!***********************************************************************

end module ocn_mep

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
