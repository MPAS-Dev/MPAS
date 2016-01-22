  if (mpas_threading_get_thread_num() == 0) then
    ! get the field for this block
    call mpas_pool_get_field(block % allFields, inname, src, 1)

    ! check that elem_name is in the dim list
    if (.not. check_element_dim(src % dimNames, elem_name)) then
      call mpas_dmpar_global_abort('Error: element dimension "' // &
        elem_name // '" was not found in field "' // src % fieldName // &
        ', which is required for the current regional stats AM' // &
        'configuration.')
    end if

    ! allocate the linked list for the field blocks
    allocate(dst) 
    nullify(dst % prev)  
    nullify(dst % next)
  
    ! copy field info 
    dst % fieldName = outname ! name is different than src
    dst % isDecomposed = .false. ! NOT decomposed! reducing across space

    dst % block => src % block
    dst % isVarArray = src % isVarArray
    dst % defaultValue = src % defaultValue
    dst % isActive = src % isActive
    dst % hasTimeDimension = src % hasTimeDimension
    dst % sendList => src % sendList
    dst % recvList => src % recvList
    dst % copyList => src % copyList

    ! no constituent names because dst is a scalar

    ! no isPersistent, dimNames, dimSizes because dst is scalar

    ! no need to allocate memory or loop over blocks, 
    ! as 0D has it already in dst % scalar

    ! add field to this block's pool
    call mpas_pool_get_subpool(block % structs, &
      REGIONAL_STATS_POOL, pool)
    call mpas_pool_add_field(pool, dst % fieldName, dst)  
    call mpas_pool_add_field(block % allFields, dst % fieldName, dst)
  end if
