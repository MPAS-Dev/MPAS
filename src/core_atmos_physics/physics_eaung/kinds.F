!==============================================================================
 MODULE kinds

 IMPLICIT NONE
 SAVE

!This module defines variable precision for all common data types.
!Laura D. Fowler (CIRA/Colorado State University).
!fowler@cira.colostate.edu.
!10-08-2007.

!==============================================================================

 INTEGER, PARAMETER:: &
    char_len  = 80,                      &
    int_kind  = kind(1),                 &
    log_kind  = kind(.true.),            &
    real_kind = selected_real_kind(6),   &
    dbl_kind  = selected_real_kind(13)

 END MODULE kinds

!==============================================================================
