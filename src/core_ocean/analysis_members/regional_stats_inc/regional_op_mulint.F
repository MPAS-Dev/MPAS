        end if
      end do
    end if
  ! integer weights
  else 
    call mpas_pool_get_array(block % allFields, &
      regions % weighting_field, integer_weights, 1)

    if (regions % weighting_function == MUL_FUNC) then
      do i = 1, solve
        if (mask(m, i) /= 0) then
