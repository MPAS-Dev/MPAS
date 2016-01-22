      end if
    end do
  ! real weights
  else if (regions % weighting_info % fieldType == MPAS_POOL_REAL) then
    call mpas_pool_get_array(block % allFields, &
      regions % weighting_field, real_weights, 1)

    if (regions % weighting_function == MUL_FUNC) then
      do i = 1, solve
        if (mask(m, i) /= 0) then
