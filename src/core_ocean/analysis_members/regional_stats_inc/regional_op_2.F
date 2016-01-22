  ! iterate over blocks
  block => start_block
  do while (associated(block))
    ! get the dimensions for this block
    if (regions % region_element == CELL_REGION) then
      call mpas_pool_get_dimension(block % dimensions, CELL_SOLVE, solve)
    else if (regions % region_element == VERTEX_REGION) then
      call mpas_pool_get_dimension(block % dimensions, VERTEX_SOLVE, solve)
    else
      call mpas_pool_get_dimension(block % dimensions, EDGE_SOLVE, solve)
    end if

    ! get the subpools
    call mpas_pool_get_subpool(block % structs, &
      regions % mask_pool, maskPool)

    ! get the arrays
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)
    call mpas_pool_get_array(maskPool, regions % mask_field, mask, 1)  

