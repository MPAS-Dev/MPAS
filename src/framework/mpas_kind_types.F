!***********************************************************************
!
!  mpas_kind_types
!
!> \brief   MPAS Kind definition module
!> \author  Michael Duda, Doug Jacobsen
!> \date    03/27/13
!> \details 
!> This module defines the kind types for basic fortran data types within MPAS.
!
!-----------------------------------------------------------------------

module mpas_kind_types

#ifdef SINGLE_PRECISION
   integer, parameter :: RKIND  = selected_real_kind(6)
#else
   integer, parameter :: RKIND  = selected_real_kind(12)
#endif

   integer, parameter :: StrKIND = 512

   contains

!***********************************************************************
!
!  routine dummy_kinds
!
!> \brief   MPAS Dummy kind routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This is a dummy routine that doesn't do anything.
!
!-----------------------------------------------------------------------
   subroutine dummy_kinds()

   end subroutine dummy_kinds

end module mpas_kind_types
