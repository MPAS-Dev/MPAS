    else
      nullify(dst % array)
    end if

    ! add field to this block's pool
    call mpas_pool_get_subpool(block % structs, &
      REGIONAL_STATS_POOL, pool)
    call mpas_pool_add_field(pool, dst % fieldName, dst)
    call mpas_pool_add_field(block % allFields, dst % fieldName, dst)
  end if
