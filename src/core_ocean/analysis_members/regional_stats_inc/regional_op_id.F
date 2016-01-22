  ! no weights
  if (regions % weighting_function == ID_FUNC) then
    do i = 1, solve
      if (mask(m, i) /= 0) then
