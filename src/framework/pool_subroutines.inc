!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           POOL SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!-----------------------------------------------------------------------
!  routine mpas_pool_set_error_level
!
!> \brief MPAS Pool Error level set routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!>  This routine sets the internal error level for pools.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_set_error_level(newErrorLevel)!{{{

      implicit none

      integer, intent(in) :: newErrorLevel

      currentErrorLevel = newErrorLevel

   end subroutine!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_create_pool
!
!> \brief MPAS Pool creation routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!>  This routine will create a new empty pool and associate newPool to this new
!>  pool location.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_create_pool(newPool, poolSize)!{{{

      implicit none

      type (mpas_pool_type), pointer :: newPool
      integer, intent(in), optional :: poolSize

      
      allocate(newPool)

      if (present(poolSize)) then
         newPool % size = poolSize
      else
         newPool % size = MPAS_POOL_TABLE_SIZE
      end if
      allocate(newPool % table(newPool % size))

   end subroutine mpas_pool_create_pool!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_destroy_pool
!
!> \brief MPAS Pool deallocation routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!>  This routine will destroy a pool associated with inPool.
!
!-----------------------------------------------------------------------
   recursive subroutine mpas_pool_destroy_pool(inPool)!{{{

      implicit none

      type (mpas_pool_type), pointer :: inPool

      integer :: i, j
      type (mpas_pool_member_type), pointer :: ptr
      type (mpas_pool_data_type), pointer :: dptr


      do i=1,inPool % size

         ptr => inPool % table(i) % head
         do while(associated(inPool % table(i) % head))
            ptr => inPool % table(i) % head
            inPool % table(i) % head => inPool % table(i) % head % next

            if (ptr % contentsType == MPAS_POOL_DIMENSION) then

               if (ptr % data % contentsDims > 0) then
                  deallocate(ptr % data % simple_int_arr)
               else
                  deallocate(ptr % data % simple_int)
               end if

            else if (ptr % contentsType == MPAS_POOL_CONFIG) then

               dptr => ptr % data

               if (dptr % contentsType == MPAS_POOL_REAL) then
                  deallocate(dptr % simple_real)
               else if (dptr % contentsType == MPAS_POOL_INTEGER) then
                  deallocate(dptr % simple_int)
               else if (dptr % contentsType == MPAS_POOL_CHARACTER) then
                  deallocate(dptr % simple_char)
               else if (dptr % contentsType == MPAS_POOL_LOGICAL) then
                  deallocate(dptr % simple_logical)
               end if

            else if (ptr % contentsType == MPAS_POOL_FIELD) then

               dptr => ptr % data

               ! Do this through brute force...
               if (associated(dptr % r0)) then
                  call mpas_deallocate_field(dptr % r0)
               else if (associated(dptr % r1)) then
                  call mpas_deallocate_field(dptr % r1)
               else if (associated(dptr % r2)) then
                  call mpas_deallocate_field(dptr % r2)
               else if (associated(dptr % r3)) then
                  call mpas_deallocate_field(dptr % r3)
               else if (associated(dptr % r4)) then
                  call mpas_deallocate_field(dptr % r4)
               else if (associated(dptr % r5)) then
                  call mpas_deallocate_field(dptr % r5)
               else if (associated(dptr % i0)) then
                  call mpas_deallocate_field(dptr % i0)
               else if (associated(dptr % i1)) then
                  call mpas_deallocate_field(dptr % i1)
               else if (associated(dptr % i2)) then
                  call mpas_deallocate_field(dptr % i2)
               else if (associated(dptr % i3)) then
                  call mpas_deallocate_field(dptr % i3)
               else if (associated(dptr % c0)) then
                  call mpas_deallocate_field(dptr % c0)
               else if (associated(dptr % c1)) then
                  call mpas_deallocate_field(dptr % c1)
               else if (associated(dptr % l0)) then
                  call mpas_deallocate_field(dptr % l0)
               else if (associated(dptr % r0a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % r0 => dptr % r0a(j)
                     call mpas_deallocate_field(dptr % r0)
                  end do
                  deallocate(dptr % r0a)
               else if (associated(dptr % r1a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % r1 => dptr % r1a(j)
                     call mpas_deallocate_field(dptr % r1)
                  end do
                  deallocate(dptr % r1a)
               else if (associated(dptr % r2a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % r2 => dptr % r2a(j)
                     call mpas_deallocate_field(dptr % r2)
                  end do
                  deallocate(dptr % r2a)
               else if (associated(dptr % r3a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % r3 => dptr % r3a(j)
                     call mpas_deallocate_field(dptr % r3)
                  end do
                  deallocate(dptr % r3a)
               else if (associated(dptr % r4a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % r4 => dptr % r4a(j)
                     call mpas_deallocate_field(dptr % r4)
                  end do
                  deallocate(dptr % r4a)
               else if (associated(dptr % r5a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % r5 => dptr % r5a(j)
                     call mpas_deallocate_field(dptr % r5)
                  end do
                  deallocate(dptr % r5a)
               else if (associated(dptr % i0a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % i0 => dptr % i0a(j)
                     call mpas_deallocate_field(dptr % i0)
                  end do
                  deallocate(dptr % i0a)
               else if (associated(dptr % i1a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % i1 => dptr % i1a(j)
                     call mpas_deallocate_field(dptr % i1)
                  end do
                  deallocate(dptr % i1a)
               else if (associated(dptr % i2a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % i2 => dptr % i2a(j)
                     call mpas_deallocate_field(dptr % i2)
                  end do
                  deallocate(dptr % i2a)
               else if (associated(dptr % i3a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % i3 => dptr % i3a(j)
                     call mpas_deallocate_field(dptr % i3)
                  end do
                  deallocate(dptr % i3a)
               else if (associated(dptr % c0a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % c0 => dptr % c0a(j)
                     call mpas_deallocate_field(dptr % c0)
                  end do
                  deallocate(dptr % c0a)
               else if (associated(dptr % c1a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % c1 => dptr % c1a(j)
                     call mpas_deallocate_field(dptr % c1)
                  end do
                  deallocate(dptr % c1a)
               else if (associated(dptr % l0a)) then
                  do j=1,dptr % contentsTimeLevs
                     dptr % l0 => dptr % l0a(j)
                     call mpas_deallocate_field(dptr % l0)
                  end do
                  deallocate(dptr % l0a)
               else
                  call pool_mesg('While destroying pool, member '//trim(ptr % key)//' has no valid field pointers.')
               end if

            else if (ptr % contentsType == MPAS_POOL_SUBPOOL) then

               call mpas_pool_destroy_pool(ptr % data % p)
               deallocate(ptr % data % p)

            end if
            deallocate(ptr)
         end do

      end do

      deallocate(inPool % table)
      deallocate(inPool)

   end subroutine mpas_pool_destroy_pool!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_empty_pool
!
!> \brief MPAS Pool empty routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!>  This routine will remove all memebers from within a pool associated with inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_empty_pool(inPool)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool

      integer :: i
      type (mpas_pool_member_type), pointer :: ptr


      do i=1,inPool % size

         ptr => inPool % table(i) % head
         do while(associated(inPool % table(i) % head))
            ptr => inPool % table(i) % head
            inPool % table(i) % head => inPool % table(i) % head % next
!TODO: now we need to deallocate all simple_* variables?
            if (ptr % contentsType == MPAS_POOL_DIMENSION) then
               if (ptr % data % contentsDims > 0) then
                  deallocate(ptr % data % simple_int_arr)
               end if
!TODO: should we empty sub-pools?
            end if
            deallocate(ptr)
         end do

      end do

      inPool % iteratorIndex = 1
      nullify(inPool % iterator)

   end subroutine mpas_pool_empty_pool!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_clone_pool
!
!> \brief MPAS Pool clone routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine assumes destPool is an empty pool. It will clone all of the members
!> from srcPool into destPool.
!
!-----------------------------------------------------------------------
   recursive subroutine mpas_pool_clone_pool(srcPool, destPool)!{{{

      implicit none

      type (mpas_pool_type), pointer :: srcPool
      type (mpas_pool_type), pointer :: destPool


      integer :: i, j
      type (mpas_pool_member_type), pointer :: ptr
      type (mpas_pool_data_type), pointer :: dptr
      type (mpas_pool_member_type), pointer :: newmem

!TODO: should we force destPool to have the same table size as srcPool?

      do i=1,srcPool % size

         ptr => srcPool % table(i) % head
         do while(associated(ptr))

            allocate(newmem)
            newmem % key = ptr % key
            newmem % keyLen = ptr % keyLen
            newmem % contentsType = ptr % contentsType
            allocate(newmem % data)
            newmem % data % contentsType = ptr % data % contentsType
            newmem % data % contentsDims = ptr % data % contentsDims
            newmem % data % contentsTimeLevs = ptr % data % contentsTimeLevs

            if (ptr % contentsType == MPAS_POOL_DIMENSION) then

               if (ptr % data % contentsDims > 0) then
                  allocate(newmem % data % simple_int_arr(size(ptr % data % simple_int_arr)))
                  newmem % data % simple_int_arr(:) = ptr % data % simple_int_arr(:)
               else
                  allocate(newmem % data % simple_int)
                  newmem % data % simple_int = ptr % data % simple_int
               end if

            else if (ptr % contentsType == MPAS_POOL_CONFIG) then

               dptr => ptr % data

               if (dptr % contentsType == MPAS_POOL_REAL) then
                  allocate(newmem % data % simple_real)
                  newmem % data % simple_real = dptr % simple_real
               else if (dptr % contentsType == MPAS_POOL_INTEGER) then
                  allocate(newmem % data % simple_int)
                  newmem % data % simple_int = dptr % simple_int
               else if (dptr % contentsType == MPAS_POOL_CHARACTER) then
                  allocate(newmem % data % simple_char)
                  newmem % data % simple_char = dptr % simple_char
               else if (dptr % contentsType == MPAS_POOL_LOGICAL) then
                  allocate(newmem % data % simple_logical)
                  newmem % data % simple_logical = dptr % simple_logical
               end if

            else if (ptr % contentsType == MPAS_POOL_FIELD) then

               dptr => ptr % data

               ! Do this through brute force...
               if (associated(dptr % r0)) then
                  call mpas_duplicate_field(dptr % r0, newmem % data % r0)
               else if (associated(dptr % r1)) then
                  call mpas_duplicate_field(dptr % r1, newmem % data % r1)
               else if (associated(dptr % r2)) then
                  call mpas_duplicate_field(dptr % r2, newmem % data % r2)
               else if (associated(dptr % r3)) then
                  call mpas_duplicate_field(dptr % r3, newmem % data % r3)
               else if (associated(dptr % r4)) then
                  call mpas_duplicate_field(dptr % r4, newmem % data % r4)
               else if (associated(dptr % r5)) then
                  call mpas_duplicate_field(dptr % r5, newmem % data % r5)
               else if (associated(dptr % i0)) then
                  call mpas_duplicate_field(dptr % i0, newmem % data % i0)
               else if (associated(dptr % i1)) then
                  call mpas_duplicate_field(dptr % i1, newmem % data % i1)
               else if (associated(dptr % i2)) then
                  call mpas_duplicate_field(dptr % i2, newmem % data % i2)
               else if (associated(dptr % i3)) then
                  call mpas_duplicate_field(dptr % i3, newmem % data % i3)
               else if (associated(dptr % c0)) then
                  call mpas_duplicate_field(dptr % c0, newmem % data % c0)
               else if (associated(dptr % c1)) then
                  call mpas_duplicate_field(dptr % c1, newmem % data % c1)
               else if (associated(dptr % l0)) then
                  call mpas_duplicate_field(dptr % l0, newmem % data % l0)
               else if (associated(dptr % r0a)) then
                  allocate(newmem % data % r0a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % r0a(j), newmem % data % r0)
                     newmem % data % r0a(j) = newmem % data % r0
                     deallocate(newmem % data % r0)
                  end do
               else if (associated(dptr % r1a)) then
                  allocate(newmem % data % r1a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % r1a(j), newmem % data % r1)
                     newmem % data % r1a(j) = newmem % data % r1
                     deallocate(newmem % data % r1)
                  end do
               else if (associated(dptr % r2a)) then
                  allocate(newmem % data % r2a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % r2a(j), newmem % data % r2)
                     newmem % data % r2a(j) = newmem % data % r2
                     deallocate(newmem % data % r2)
                  end do
               else if (associated(dptr % r3a)) then
                  allocate(newmem % data % r3a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % r3a(j), newmem % data % r3)
                     newmem % data % r3a(j) = newmem % data % r3
                     deallocate(newmem % data % r3)
                  end do
               else if (associated(dptr % r4a)) then
                  allocate(newmem % data % r4a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % r4a(j), newmem % data % r4)
                     newmem % data % r4a(j) = newmem % data % r4
                     deallocate(newmem % data % r4)
                  end do
               else if (associated(dptr % r5a)) then
                  allocate(newmem % data % r5a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % r5a(j), newmem % data % r5)
                     newmem % data % r5a(j) = newmem % data % r5
                     deallocate(newmem % data % r5)
                  end do
               else if (associated(dptr % i0a)) then
                  allocate(newmem % data % i0a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % i0a(j), newmem % data % i0)
                     newmem % data % i0a(j) = newmem % data % i0
                     deallocate(newmem % data % i0)
                  end do
               else if (associated(dptr % i1a)) then
                  allocate(newmem % data % i1a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % i1a(j), newmem % data % i1)
                     newmem % data % i1a(j) = newmem % data % i1
                     deallocate(newmem % data % i1)
                  end do
               else if (associated(dptr % i2a)) then
                  allocate(newmem % data % i2a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % i2a(j), newmem % data % i2)
                     newmem % data % i2a(j) = newmem % data % i2
                     deallocate(newmem % data % i2)
                  end do
               else if (associated(dptr % i3a)) then
                  allocate(newmem % data % i3a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % i3a(j), newmem % data % i3)
                     newmem % data % i3a(j) = newmem % data % i3
                     deallocate(newmem % data % i3)
                  end do
               else if (associated(dptr % c0a)) then
                  allocate(newmem % data % c0a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % c0a(j), newmem % data % c0)
                     newmem % data % c0a(j) = newmem % data % c0
                     deallocate(newmem % data % c0)
                  end do
               else if (associated(dptr % c1a)) then
                  allocate(newmem % data % c1a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % c1a(j), newmem % data % c1)
                     newmem % data % c1a(j) = newmem % data % c1
                     deallocate(newmem % data % c1)
                  end do
               else if (associated(dptr % l0a)) then
                  allocate(newmem % data % l0a(dptr % contentsTimeLevs))
                  do j=1,dptr % contentsTimeLevs
                     call mpas_duplicate_field(dptr % l0a(j), newmem % data % l0)
                     newmem % data % l0a(j) = newmem % data % l0
                     deallocate(newmem % data % c0)
                  end do
               else
                  call pool_mesg('While cloning pool, member '//trim(ptr % key)//' has no valid field pointers.')
               end if

            else if (ptr % contentsType == MPAS_POOL_SUBPOOL) then

                call mpas_pool_create_pool(newmem % data % p, poolSize = ptr % data % p % size)
                call mpas_pool_clone_pool(ptr % data % p, newmem % data % p)

            end if

            if (.not. pool_add_member(destPool, newmem % key, newmem)) then
               call pool_mesg('Error: Had problems adding '//trim(newmem % key)//' to clone of pool.')
            end if

            ptr => ptr % next
         end do

      end do

   end subroutine mpas_pool_clone_pool!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_copy_pool
!
!> \brief MPAS Pool copy routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine assumes srcPool and destPool have identical members. It will
!> copy the data from the members of srcPool into the members of destPool.
!
!-----------------------------------------------------------------------
   recursive subroutine mpas_pool_copy_pool(srcPool, destPool)!{{{

      implicit none

      type (mpas_pool_type), pointer :: srcPool
      type (mpas_pool_type), pointer :: destPool


      integer :: i, j
      type (mpas_pool_member_type), pointer :: ptr
      type (mpas_pool_data_type), pointer :: dptr
      type (mpas_pool_data_type), pointer :: mem


      do i=1,srcPool % size

         ptr => srcPool % table(i) % head
         do while(associated(ptr))

            if (ptr % contentsType == MPAS_POOL_DIMENSION) then

               mem => pool_get_member(destPool, ptr % key, MPAS_POOL_DIMENSION)
!TODO: ensure that mem is associated
               if (ptr % data % contentsDims > 0) then
                  mem % simple_int_arr(:) = ptr % data % simple_int_arr(:)
               else
                  mem % simple_int = ptr % data % simple_int
               end if

            else if (ptr % contentsType == MPAS_POOL_CONFIG) then

               dptr => ptr % data

               mem => pool_get_member(destPool, ptr % key, MPAS_POOL_CONFIG)
               if (dptr % contentsType == MPAS_POOL_REAL) then
                  mem % simple_real = dptr % simple_real
               else if (dptr % contentsType == MPAS_POOL_INTEGER) then
                  mem % simple_int = dptr % simple_int
               else if (dptr % contentsType == MPAS_POOL_CHARACTER) then
                  mem % simple_char = dptr % simple_char
               else if (dptr % contentsType == MPAS_POOL_LOGICAL) then
                  mem % simple_logical = dptr % simple_logical
               end if

            else if (ptr % contentsType == MPAS_POOL_FIELD) then

               dptr => ptr % data

               ! Do this through brute force...
               mem => pool_get_member(destPool, ptr % key, MPAS_POOL_FIELD)
               if (associated(dptr % r0)) then
                  call mpas_duplicate_field(dptr % r0, mem % r0, copy_array_only=.true.)
               else if (associated(dptr % r1)) then
                  call mpas_duplicate_field(dptr % r1, mem % r1, copy_array_only=.true.)
               else if (associated(dptr % r2)) then
                  call mpas_duplicate_field(dptr % r2, mem % r2, copy_array_only=.true.)
               else if (associated(dptr % r3)) then
                  call mpas_duplicate_field(dptr % r3, mem % r3, copy_array_only=.true.)
               else if (associated(dptr % r4)) then
                  call mpas_duplicate_field(dptr % r4, mem % r4, copy_array_only=.true.)
               else if (associated(dptr % r5)) then
                  call mpas_duplicate_field(dptr % r5, mem % r5, copy_array_only=.true.)
               else if (associated(dptr % i0)) then
                  call mpas_duplicate_field(dptr % i0, mem % i0, copy_array_only=.true.)
               else if (associated(dptr % i1)) then
                  call mpas_duplicate_field(dptr % i1, mem % i1, copy_array_only=.true.)
               else if (associated(dptr % i2)) then
                  call mpas_duplicate_field(dptr % i2, mem % i2, copy_array_only=.true.)
               else if (associated(dptr % i3)) then
                  call mpas_duplicate_field(dptr % i3, mem % i3, copy_array_only=.true.)
               else if (associated(dptr % c0)) then
                  call mpas_duplicate_field(dptr % c0, mem % c0, copy_array_only=.true.)
               else if (associated(dptr % c1)) then
                  call mpas_duplicate_field(dptr % c1, mem % c1, copy_array_only=.true.)
               else if (associated(dptr % l0)) then
                  call mpas_duplicate_field(dptr % l0, mem % l0, copy_array_only=.true.)
               else if (associated(dptr % r0a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % r0 => mem % r0a(j)
                     call mpas_duplicate_field(dptr % r0a(j), mem % r0, copy_array_only=.true.)
                     nullify(mem % r0)
                  end do
               else if (associated(dptr % r1a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % r1 => mem % r1a(j)
                     call mpas_duplicate_field(dptr % r1a(j), mem % r1, copy_array_only=.true.)
                     nullify(mem % r1)
                  end do
               else if (associated(dptr % r2a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % r2 => mem % r2a(j)
                     call mpas_duplicate_field(dptr % r2a(j), mem % r2, copy_array_only=.true.)
                     nullify(mem % r2)
                  end do
               else if (associated(dptr % r3a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % r3 => mem % r3a(j)
                     call mpas_duplicate_field(dptr % r3a(j), mem % r3, copy_array_only=.true.)
                     nullify(mem % r3)
                  end do
               else if (associated(dptr % r4a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % r4 => mem % r4a(j)
                     call mpas_duplicate_field(dptr % r4a(j), mem % r4, copy_array_only=.true.)
                     nullify(mem % r4)
                  end do
               else if (associated(dptr % r5a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % r5 => mem % r5a(j)
                     call mpas_duplicate_field(dptr % r5a(j), mem % r5, copy_array_only=.true.)
                     nullify(mem % r5)
                  end do
               else if (associated(dptr % i0a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % i0 => mem % i0a(j)
                     call mpas_duplicate_field(dptr % i0a(j), mem % i0, copy_array_only=.true.)
                     nullify(mem % i0)
                  end do
               else if (associated(dptr % i1a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % i1 => mem % i1a(j)
                     call mpas_duplicate_field(dptr % i1a(j), mem % i1, copy_array_only=.true.)
                     nullify(mem % i1)
                  end do
               else if (associated(dptr % i2a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % i2 => mem % i2a(j)
                     call mpas_duplicate_field(dptr % i2a(j), mem % i2, copy_array_only=.true.)
                     nullify(mem % i2)
                  end do
               else if (associated(dptr % i3a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % i3 => mem % i3a(j)
                     call mpas_duplicate_field(dptr % i3a(j), mem % i3, copy_array_only=.true.)
                     nullify(mem % i3)
                  end do
               else if (associated(dptr % c0a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % c0 => mem % c0a(j)
                     call mpas_duplicate_field(dptr % c0a(j), mem % c0, copy_array_only=.true.)
                     nullify(mem % c0)
                  end do
               else if (associated(dptr % c1a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % c1 => mem % c1a(j)
                     call mpas_duplicate_field(dptr % c1a(j), mem % c1, copy_array_only=.true.)
                     nullify(mem % c1)
                  end do
               else if (associated(dptr % l0a)) then
                  do j=1,dptr % contentsTimeLevs
                     mem % l0 => mem % l0a(j)
                     call mpas_duplicate_field(dptr % l0a(j), mem % l0, copy_array_only=.true.)
                     nullify(mem % l0)
                  end do
               else
                  call pool_mesg('While copying pool, member '//trim(ptr % key)//' has no valid field pointers.')
               end if

            else if (ptr % contentsType == MPAS_POOL_SUBPOOL) then

                mem => pool_get_member(destPool, ptr % key, MPAS_POOL_SUBPOOL)
                call mpas_pool_copy_pool(ptr % data % p, mem % p)

            end if

            ptr => ptr % next
         end do

      end do

   end subroutine mpas_pool_copy_pool!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_0d_real
!
!> \brief MPAS Pool 0D Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 0D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_0d_real(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field0DReal), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = 1
      newmem % data % r0 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_0d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_1d_real
!
!> \brief MPAS Pool 1D Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 1D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_1d_real(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field1DReal), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 1
      newmem % data % contentsTimeLevs = 1
      newmem % data % r1 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_1d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_2d_real
!
!> \brief MPAS Pool 2D Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 2D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_2d_real(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field2DReal), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 2
      newmem % data % contentsTimeLevs = 1
      newmem % data % r2 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_2d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_3d_real
!
!> \brief MPAS Pool 3D Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 3D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_3d_real(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field3DReal), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 3
      newmem % data % contentsTimeLevs = 1
      newmem % data % r3 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_3d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_4d_real
!
!> \brief MPAS Pool 4D Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 4D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_4d_real(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field4DReal), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 4
      newmem % data % contentsTimeLevs = 1
      newmem % data % r4 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_4d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_5d_real
!
!> \brief MPAS Pool 5D Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 5D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_5d_real(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field5DReal), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 5
      newmem % data % contentsTimeLevs = 1
      newmem % data % r5 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_5d_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_0d_int
!
!> \brief MPAS Pool 0D Integer field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 0D integer field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_0d_int(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field0DInteger), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_INTEGER
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = 1
      newmem % data % i0 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_0d_int!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_1d_int
!
!> \brief MPAS Pool 1D Integer field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 1D integer field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_1d_int(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field1DInteger), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_INTEGER
      newmem % data % contentsDims = 1
      newmem % data % contentsTimeLevs = 1
      newmem % data % i1 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_1d_int!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_2d_int
!
!> \brief MPAS Pool 2D Integer field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 2D integer field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_2d_int(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field2DInteger), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_INTEGER
      newmem % data % contentsDims = 2
      newmem % data % contentsTimeLevs = 1
      newmem % data % i2 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_2d_int!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_3d_int
!
!> \brief MPAS Pool 3D Integer field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 3D integer field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_3d_int(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field3DInteger), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_INTEGER
      newmem % data % contentsDims = 3
      newmem % data % contentsTimeLevs = 1
      newmem % data % i3 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_3d_int!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_0d_char
!
!> \brief MPAS Pool 0D Character field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 0D character field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_0d_char(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field0DChar), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_CHARACTER
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = 1
      newmem % data % c0 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_0d_char!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_1d_char
!
!> \brief MPAS Pool 1D Character field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts field into inPool when field is a 1D character field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_1d_char(inPool, key, field)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field1DChar), pointer :: field

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_CHARACTER
      newmem % data % contentsDims = 1
      newmem % data % contentsTimeLevs = 1
      newmem % data % c1 => field
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_1d_char!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_0d_reals
!
!> \brief MPAS Pool 0D Multi-level Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 0D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_0d_reals(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field0DReal), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % r0 => fields(1)
      else
         newmem % data % r0a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_0d_reals!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_1d_reals
!
!> \brief MPAS Pool 1D Multi-level Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 1D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_1d_reals(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field1DReal), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 1
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % r1 => fields(1)
      else
         newmem % data % r1a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_1d_reals!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_2d_reals
!
!> \brief MPAS Pool 2D Multi-level Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 2D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_2d_reals(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field2DReal), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 2
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % r2 => fields(1)
      else
         newmem % data % r2a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_2d_reals!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_3d_reals
!
!> \brief MPAS Pool 3D Multi-level Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 3D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_3d_reals(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field3DReal), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 3
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % r3 => fields(1)
      else
         newmem % data % r3a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_3d_reals!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_4d_reals
!
!> \brief MPAS Pool 4D Multi-level Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 4D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_4d_reals(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field4DReal), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 4
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % r4 => fields(1)
      else
         newmem % data % r4a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_4d_reals!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_5d_reals
!
!> \brief MPAS Pool 5D Multi-level Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 5D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_5d_reals(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field5DReal), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 5
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % r5 => fields(1)
      else
         newmem % data % r5a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_5d_reals!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_0d_ints
!
!> \brief MPAS Pool 0D Multi-level Integer field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 0D integer field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_0d_ints(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field0DInteger), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_INTEGER
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % i0 => fields(1)
      else
         newmem % data % i0a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_0d_ints!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_1d_ints
!
!> \brief MPAS Pool 1D Multi-level Integer field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 1D integer field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_1d_ints(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field1DInteger), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_INTEGER
      newmem % data % contentsDims = 1
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % i1 => fields(1)
      else
         newmem % data % i1a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_1d_ints!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_2d_ints
!
!> \brief MPAS Pool 2D Multi-level integer field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 2D integer field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_2d_ints(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field2DInteger), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_INTEGER
      newmem % data % contentsDims = 2
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % i2 => fields(1)
      else
         newmem % data % i2a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_2d_ints!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_3d_ints
!
!> \brief MPAS Pool 3D Multi-level Integer field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 3D integer field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_3d_ints(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field3DInteger), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_INTEGER
      newmem % data % contentsDims = 3
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % i3 => fields(1)
      else
         newmem % data % i3a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_3d_ints!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_0d_chars
!
!> \brief MPAS Pool 0D Multi-level Character field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 0D character field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_0d_chars(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field0DChar), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_CHARACTER
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % c0 => fields(1)
      else
         newmem % data % c0a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_0d_chars!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_1d_chars
!
!> \brief MPAS Pool 1D Multi-level Character field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level 1D character field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_1d_chars(inPool, key, fields)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (field1DChar), dimension(:), pointer :: fields

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_FIELD

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_CHARACTER
      newmem % data % contentsDims = 1
      newmem % data % contentsTimeLevs = size(fields)
      if (newmem % data % contentsTimeLevs == 1) then
         newmem % data % c1 => fields(1)
      else
         newmem % data % c1a => fields
      end if
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_field_1d_chars!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_info
!
!> \brief MPAS Pool Field Information Query subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a data structure containing information related to the
!>  field in inPool with the name key
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_info(inPool, key, info)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (mpas_pool_field_info_type), intent(out) :: info

      integer :: hash, endl
      type (mpas_pool_member_type), pointer :: ptr


      endl = len_trim(key)
      call pool_hash(hash, key, endl)

      hash = mod(hash, inPool % size) + 1

      ptr => inPool % table(hash) % head
      do while (associated(ptr))
         if (ptr % contentsType == MPAS_POOL_FIELD) then
            if (endl == ptr % keyLen) then
               if (key(1:endl) == ptr % key(1:endl)) then

                  info % fieldType = ptr % data % contentsType
                  info % nDims = ptr % data % contentsDims
                  info % nTimeLevels = ptr % data % contentsTimeLevs

                  exit
               end if
            end if
         end if
         ptr => ptr % next
      end do

      if (.not. associated(ptr)) then
         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_get_field_info!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_0d_real
!
!> \brief MPAS Pool 0D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_0d_real(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field0DReal), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_REAL) then
            call pool_mesg('Error: Field '//trim(key)//' is not type real.')
         end if
         if (mem % contentsDims /= 0) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 0-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % r0
         else
            field => mem % r0a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_0d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_1d_real
!
!> \brief MPAS Pool 1D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_1d_real(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field1DReal), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_REAL) then
            call pool_mesg('Error: Field '//trim(key)//' is not type real.')
         end if
         if (mem % contentsDims /= 1) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 1-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % r1
         else
            field => mem % r1a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_1d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_2d_real
!
!> \brief MPAS Pool 2D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_2d_real(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field2DReal), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_REAL) then
            call pool_mesg('Error: Field '//trim(key)//' is not type real.')
         end if
         if (mem % contentsDims /= 2) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 2-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % r2
         else
            field => mem % r2a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_2d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_3d_real
!
!> \brief MPAS Pool 3D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_3d_real(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field3DReal), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_REAL) then
            call pool_mesg('Error: Field '//trim(key)//' is not type real.')
         end if
         if (mem % contentsDims /= 3) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 3-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % r3
         else
            field => mem % r3a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_3d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_4d_real
!
!> \brief MPAS Pool 4D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_4d_real(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field4DReal), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_REAL) then
            call pool_mesg('Error: Field '//trim(key)//' is not type real.')
         end if
         if (mem % contentsDims /= 4) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 4-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % r4
         else
            field => mem % r4a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_4d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_5d_real
!
!> \brief MPAS Pool 5D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_5d_real(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field5DReal), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_REAL) then
            call pool_mesg('Error: Field '//trim(key)//' is not type real.')
         end if
         if (mem % contentsDims /= 5) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 5-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % r5
         else
            field => mem % r5a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_5d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_0d_int
!
!> \brief MPAS Pool 0D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_0d_int(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field0DInteger), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_INTEGER) then
            call pool_mesg('Error: Field '//trim(key)//' is not type integer.')
         end if
         if (mem % contentsDims /= 0) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 0-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % i0
         else
            field => mem % i0a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_0d_int!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_1d_int
!
!> \brief MPAS Pool 1D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_1d_int(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field1DInteger), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_INTEGER) then
            call pool_mesg('Error: Field '//trim(key)//' is not type integer.')
         end if
         if (mem % contentsDims /= 1) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 1-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % i1
         else
            field => mem % i1a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_1d_int!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_2d_int
!
!> \brief MPAS Pool 2D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_2d_int(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field2DInteger), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_INTEGER) then
            call pool_mesg('Error: Field '//trim(key)//' is not type integer.')
         end if
         if (mem % contentsDims /= 2) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 2-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % i2
         else
            field => mem % i2a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_2d_int!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_3d_int
!
!> \brief MPAS Pool 3D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_3d_int(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field3DInteger), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_INTEGER) then
            call pool_mesg('Error: Field '//trim(key)//' is not type integer.')
         end if
         if (mem % contentsDims /= 3) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 3-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % i3
         else
            field => mem % i3a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_3d_int!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_0d_char
!
!> \brief MPAS Pool 0D Character field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_0d_char(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field0DChar), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_CHARACTER) then
            call pool_mesg('Error: Field '//trim(key)//' is not type character.')
         end if
         if (mem % contentsDims /= 0) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 0-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % c0
         else
            field => mem % c0a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_0d_char!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_1d_char
!
!> \brief MPAS Pool 1D Character field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_1d_char(inPool, key, field, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field1DChar), pointer :: field
      integer, intent(in), optional :: timeLevel

      type (mpas_pool_data_type), pointer :: mem
      integer :: local_timeLevel


      if (present(timeLevel)) then
         local_timeLevel = timeLevel
      else
         local_timeLevel = 1
      end if

      mem => pool_get_member(inPool, key, MPAS_POOL_FIELD)

      nullify(field)
      if (associated(mem)) then

         if (mem % contentsType /= MPAS_POOL_CHARACTER) then
            call pool_mesg('Error: Field '//trim(key)//' is not type character.')
         end if
         if (mem % contentsDims /= 1) then
            call pool_mesg('Error: Field '//trim(key)//' is not a 1-d field.')
         end if
         if ((mem % contentsTimeLevs > 1) .and. (.not. present(timeLevel))) then
            call pool_mesg('Error: Field '//trim(key)//' has more than one time level, but no timeLevel argument given.')
         end if
         if (mem % contentsTimeLevs < local_timeLevel) then
            call pool_mesg('Error: Field '//trim(key)//' has too few time levels.')
         end if
         
         if (mem % contentsTimeLevs == 1) then
            field => mem % c1
         else
            field => mem % c1a(local_timeLevel)
         end if

      else

         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')

      end if

   end subroutine mpas_pool_get_field_1d_char!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_0d_real
!
!> \brief MPAS Pool 0D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_0d_real(inPool, key, scalar, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      real (kind=RKIND), pointer :: scalar
      integer, intent(in), optional :: timeLevel

      type (field0DReal), pointer :: field


      call mpas_pool_get_field_0d_real(inPool, key, field, timeLevel)

!TODO: the lack of a pointer here leads to different behavior from array fields...
      if (associated(field)) scalar = field % scalar

   end subroutine mpas_pool_get_array_0d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_1d_real
!
!> \brief MPAS Pool 1D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_1d_real(inPool, key, array, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      real (kind=RKIND), dimension(:), pointer :: array
      integer, intent(in), optional :: timeLevel

      type (field1DReal), pointer :: field


      call mpas_pool_get_field_1d_real(inPool, key, field, timeLevel)

      if (associated(field)) array => field % array

   end subroutine mpas_pool_get_array_1d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_2d_real
!
!> \brief MPAS Pool 2D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_2d_real(inPool, key, array, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      real (kind=RKIND), dimension(:,:), pointer :: array
      integer, intent(in), optional :: timeLevel

      type (field2DReal), pointer :: field


      call mpas_pool_get_field_2d_real(inPool, key, field, timeLevel)

      if (associated(field)) array => field % array

   end subroutine mpas_pool_get_array_2d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_3d_real
!
!> \brief MPAS Pool 3D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_3d_real(inPool, key, array, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      real (kind=RKIND), dimension(:,:,:), pointer :: array
      integer, intent(in), optional :: timeLevel

      type (field3DReal), pointer :: field


      call mpas_pool_get_field_3d_real(inPool, key, field, timeLevel)

      if (associated(field)) array => field % array

   end subroutine mpas_pool_get_array_3d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_4d_real
!
!> \brief MPAS Pool 4D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_4d_real(inPool, key, array, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      real (kind=RKIND), dimension(:,:,:,:), pointer :: array
      integer, intent(in), optional :: timeLevel

      type (field4DReal), pointer :: field


      call mpas_pool_get_field_4d_real(inPool, key, field, timeLevel)

      if (associated(field)) array => field % array

   end subroutine mpas_pool_get_array_4d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_5d_real
!
!> \brief MPAS Pool 5D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_5d_real(inPool, key, array, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: array
      integer, intent(in), optional :: timeLevel

      type (field5DReal), pointer :: field


      call mpas_pool_get_field_5d_real(inPool, key, field, timeLevel)

      if (associated(field)) array => field % array

   end subroutine mpas_pool_get_array_5d_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_0d_int
!
!> \brief MPAS Pool 0D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_0d_int(inPool, key, scalar, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, pointer :: scalar
      integer, intent(in), optional :: timeLevel

      type (field0DInteger), pointer :: field


      call mpas_pool_get_field_0d_int(inPool, key, field, timeLevel)

!TODO: the lack of a pointer here leads to different behavior from array fields...
      if (associated(field)) scalar = field % scalar

   end subroutine mpas_pool_get_array_0d_int!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_1d_int
!
!> \brief MPAS Pool 1D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_1d_int(inPool, key, array, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, dimension(:), pointer :: array
      integer, intent(in), optional :: timeLevel

      type (field1DInteger), pointer :: field


      call mpas_pool_get_field_1d_int(inPool, key, field, timeLevel)

      if (associated(field)) array => field % array

   end subroutine mpas_pool_get_array_1d_int!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_2d_int
!
!> \brief MPAS Pool 2D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_2d_int(inPool, key, array, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, dimension(:,:), pointer :: array
      integer, intent(in), optional :: timeLevel

      type (field2DInteger), pointer :: field


      call mpas_pool_get_field_2d_int(inPool, key, field, timeLevel)

      if (associated(field)) array => field % array

   end subroutine mpas_pool_get_array_2d_int!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_3d_int
!
!> \brief MPAS Pool 3D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_3d_int(inPool, key, array, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, dimension(:,:,:), pointer :: array
      integer, intent(in), optional :: timeLevel

      type (field3DInteger), pointer :: field


      call mpas_pool_get_field_3d_int(inPool, key, field, timeLevel)

      if (associated(field)) array => field % array

   end subroutine mpas_pool_get_array_3d_int!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_0d_char
!
!> \brief MPAS Pool 0D Character field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_0d_char(inPool, key, string, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      character (len=StrKIND), pointer :: string
      integer, intent(in), optional :: timeLevel

      type (field0DChar), pointer :: field


      call mpas_pool_get_field_0d_char(inPool, key, field, timeLevel)

!TODO: the lack of a pointer here leads to different behavior from array fields...
      if (associated(field)) string = field % scalar

   end subroutine mpas_pool_get_array_0d_char!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_1d_char
!
!> \brief MPAS Pool 1D Character field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associated with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_1d_char(inPool, key, array, timeLevel)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      character (len=StrKIND), dimension(:), pointer :: array
      integer, intent(in), optional :: timeLevel

      type (field1DChar), pointer :: field


      call mpas_pool_get_field_1d_char(inPool, key, field, timeLevel)

      if (associated(field)) array => field % array

   end subroutine mpas_pool_get_array_1d_char!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_config_real
!
!> \brief MPAS Pool Real Config Insertion Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts a real value as a config option into inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_config_real(inPool, key, value)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      real (kind=RKIND), intent(in) :: value

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_CONFIG

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_REAL
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = 0
      allocate(newmem % data % simple_real)
      newmem % data % simple_real = value
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_config_real!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_config_int
!
!> \brief MPAS Pool Integer Config Insertion Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts a integer value as a config option into inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_config_int(inPool, key, value)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in) :: value

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_CONFIG

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_INTEGER
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = 0
      allocate(newmem % data % simple_int)
      newmem % data % simple_int = value
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_config_int!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_config_char
!
!> \brief MPAS Pool Character Config Insertion Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts a character string as a config option into inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_config_char(inPool, key, value)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      character (len=*), intent(in) :: value

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_CONFIG

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_CHARACTER
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = 0
      allocate(newmem % data % simple_char)
      newmem % data % simple_char = value    !TODO: should we check that value is not more than StrKIND characters?
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_config_char!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_config_logical
!
!> \brief MPAS Pool Logical Config Insertion Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts a logical flag as a config option into inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_config_logical(inPool, key, value)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      logical, intent(in) :: value

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_CONFIG

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_LOGICAL
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = 0
      allocate(newmem % data % simple_logical)
      newmem % data % simple_logical = value
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_config_logical!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_config_real
!
!> \brief MPAS Pool Real Config Access Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns the value associated with a config option with the
!> name key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_config_real(inPool, key, value)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      real (kind=RKIND), pointer :: value

      type (mpas_pool_data_type), pointer :: mem


      mem => pool_get_member(inPool, key, MPAS_POOL_CONFIG)

      if (associated(mem)) then
         if (mem % contentsType /= MPAS_POOL_REAL) then
            call pool_mesg('Error: Config '//trim(key)//' is not type real.')
         end if
         value => mem % simple_real
      else
         call pool_mesg('Error: Config '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_get_config_real!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_config_int
!
!> \brief MPAS Pool Integer Config Access Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns the value associated with a config option with the
!> name key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_config_int(inPool, key, value)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, pointer :: value

      type (mpas_pool_data_type), pointer :: mem


      mem => pool_get_member(inPool, key, MPAS_POOL_CONFIG)

      if (associated(mem)) then
         if (mem % contentsType /= MPAS_POOL_INTEGER) then
            call pool_mesg('Error: Config '//trim(key)//' is not type integer.')
         end if
         value => mem % simple_int
      else
         call pool_mesg('Error: Config '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_get_config_int!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_config_char
!
!> \brief MPAS Pool Character Config Access Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns the value associated with a config option with the
!> name key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_config_char(inPool, key, value)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      character (len=StrKIND), pointer :: value

      type (mpas_pool_data_type), pointer :: mem


      mem => pool_get_member(inPool, key, MPAS_POOL_CONFIG)

      if (associated(mem)) then
         if (mem % contentsType /= MPAS_POOL_CHARACTER) then
            call pool_mesg('Error: Config '//trim(key)//' is not type character.')
         end if
         value => mem % simple_char  ! TODO: Should we check that value is large enough?
      else
         call pool_mesg('Error: Config '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_get_config_char!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_config_logical
!
!> \brief MPAS Pool Logical Config Access Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns the value associated with a config option with the
!> name key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_config_logical(inPool, key, value)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      logical, pointer :: value

      type (mpas_pool_data_type), pointer :: mem


      mem => pool_get_member(inPool, key, MPAS_POOL_CONFIG)

      if (associated(mem)) then
         if (mem % contentsType /= MPAS_POOL_LOGICAL) then
            call pool_mesg('Error: Config '//trim(key)//' is not type logical.')
         end if
         value => mem % simple_logical
      else
         call pool_mesg('Error: Config '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_get_config_logical!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_dimension_0d
!
!> \brief MPAS Pool 0D Dimension Insertion routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts a 0D dimension into inPool, and associated it with key.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_dimension_0d(inPool, key, dim)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in) :: dim

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_DIMENSION

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_INTEGER
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = 0
      allocate(newmem % data % simple_int)
      newmem % data % simple_int = dim
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_dimension_0d!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_dimension_1d
!
!> \brief MPAS Pool 1D Dimension Insertion routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts a 1D dimension into inPool, and associated it with key.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_dimension_1d(inPool, key, dims)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      integer, dimension(:), intent(in) :: dims

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_DIMENSION

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_INTEGER
      newmem % data % contentsDims = size(dims)
      newmem % data % contentsTimeLevs = 0
      allocate(newmem % data % simple_int_arr(newmem % data % contentsDims))
      newmem % data % simple_int_arr(:) = dims(:)
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data % simple_int_arr)
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_dimension_1d!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_dimension_0d
!
!> \brief MPAS Pool 0D Dimension Access subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns the value of the 0D dimension associated with key in
!>  inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_dimension_0d(inPool, key, dim)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, intent(out) :: dim

      type (mpas_pool_data_type), pointer :: mem


      mem => pool_get_member(inPool, key, MPAS_POOL_DIMENSION)

      if (associated(mem)) then
         if (mem % contentsDims /= 0) then
            call pool_mesg('Error: Dimension '//trim(key)//' is not a scalar.')
         else
            dim = mem % simple_int
         end if
      else
         call pool_mesg('Error: Dimension '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_get_dimension_0d!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_dimension_1d
!
!> \brief MPAS Pool 1D Dimension Access subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns the value of the 1D dimension associated with key in
!>  inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_dimension_1d(inPool, key, dims)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, intent(out), dimension(:) :: dims

      type (mpas_pool_data_type), pointer :: mem


      mem => pool_get_member(inPool, key, MPAS_POOL_DIMENSION)

      if (associated(mem)) then
         if (mem % contentsDims /= size(dims)) then
            call pool_mesg('Error: Dimension '//trim(key)//' is not the correct size.')
         else
            dims(:) = mem % simple_int_arr(:)
         end if
      else
         call pool_mesg('Error: Dimension '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_get_dimension_1d!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_add_subpool
!
!> \brief MPAS Pool Subpool insertion routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts a subpool (subPool) into inPool and associated it with
!>  the name key.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_subpool(inPool, key, subPool)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (mpas_pool_type), intent(in), target :: subPool


      type (mpas_pool_member_type), pointer :: newmem

      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_SUBPOOL

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_SUBPOOL
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = 0
      newmem % data % p => subPool
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_subpool!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_subpool
!
!> \brief MPAS Pool Subpool access subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the subpool named key within inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_subpool(inPool, key, subPool)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (mpas_pool_type), pointer :: subPool

      type (mpas_pool_data_type), pointer :: mem


      mem => pool_get_member(inPool, key, MPAS_POOL_SUBPOOL)

      if (associated(mem)) then
         subPool => mem % p
      else
         call pool_mesg('Error: Sub-pool '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_get_subpool!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_add_package
!
!> \brief MPAS Pool Package insertion subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine inserts a package into a inPool and associates it with the
!>  name key.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_package(inPool, key, value)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      logical, intent(in) :: value

      type (mpas_pool_member_type), pointer :: newmem


      allocate(newmem)
      newmem % key = key
      newmem % keyLen = len_trim(key)
      newmem % contentsType = MPAS_POOL_PACKAGE

      allocate(newmem % data)
      newmem % data % contentsType = MPAS_POOL_LOGICAL
      newmem % data % contentsDims = 0
      newmem % data % contentsTimeLevs = 0
      allocate(newmem % data % simple_logical)
      newmem % data % simple_logical = value
   
      if (.not. pool_add_member(inPool, key, newmem)) then
         deallocate(newmem % data)
         deallocate(newmem)
      end if

   end subroutine mpas_pool_add_package!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_package
!
!> \brief MPAS Pool Package access subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine sets the package pointer to point to the logical associated
!>  with the package in inPool with name key.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_package(inPool, key, package)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      logical, pointer :: package

      type (mpas_pool_data_type), pointer :: mem


      mem => pool_get_member(inPool, key, MPAS_POOL_PACKAGE)

      if (associated(mem)) then
         package => mem % simple_logical
      else
         call pool_mesg('Error: Package '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_get_package!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_remove_field
!
!> \brief MPAS Pool Field Removal Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine removes a field with the name key from inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_remove_field(inPool, key)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key


      if (.not. pool_remove_member(inPool, key, MPAS_POOL_FIELD)) then
         call pool_mesg('Error: Field '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_remove_field!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_remove_config
!
!> \brief MPAS Pool Config Removal Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine removes a config with the name key from inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_remove_config(inPool, key)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key


      if (.not. pool_remove_member(inPool, key, MPAS_POOL_CONFIG)) then
         call pool_mesg('Error: Config '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_remove_config!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_remove_dimension
!
!> \brief MPAS Pool Dimension Removal Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine removes a dimension with the name key from inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_remove_dimension(inPool, key)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key


!TODO: special case: if we have a 1-d dimension, need to deallocate simple_int_arr
      if (.not. pool_remove_member(inPool, key, MPAS_POOL_DIMENSION)) then
         call pool_mesg('Error: Dimension '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_remove_dimension!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_remove_subpool
!
!> \brief MPAS Pool Subpool Removal Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine removes a subpool with the name key from inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_remove_subpool(inPool, key)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key


      if (.not. pool_remove_member(inPool, key, MPAS_POOL_SUBPOOL)) then
         call pool_mesg('Error: Sub-pool '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_remove_subpool!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_remove_package
!
!> \brief MPAS Pool Package Removal Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine removes a package with the name key from inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_remove_package(inPool, key)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key


      if (.not. pool_remove_member(inPool, key, MPAS_POOL_PACKAGE)) then
         call pool_mesg('Error: Package '//trim(key)//' not found in pool.')
      end if

   end subroutine mpas_pool_remove_package!}}}


!-----------------------------------------------------------------------
!  routine mpas_pool_begin_iteration
!
!> \brief MPAS Pool Begin Iteration Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine sets up the pool's internal iterator to iterate over fields.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_begin_iteration(inPool)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool

      integer :: i


      do i=1,inPool % size
         if (associated(inPool % table(i) % head)) exit
      end do
      inPool % iteratorIndex = i
      if (i <= inPool % size) then
         inPool % iterator => inPool % table(i) % head
      else
         nullify(inPool % iterator)
      end if

   end subroutine mpas_pool_begin_iteration!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_next_member
!
!> \brief MPAS Pool Iterate To Next Member subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This function advances the internal iterator to the next member in the pool,
!>  and returns an iterator type for the current member, if one exists. The function
!>  returns .true. if a valid member was returned, and .false. if there are no members
!>  left to be iterated over.
!
!-----------------------------------------------------------------------
   logical function mpas_pool_get_next_member(inPool, iterator)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      type (mpas_pool_iterator_type), intent(inout) :: iterator

      integer :: i
 

      !
      ! As long as there are members left to be iterated over, the inPool%iterator
      !   should always be pointing to the next member to be returned
      !
      if (associated(inPool % iterator)) then
         iterator % memberName = inPool % iterator % key
         iterator % memberType = inPool % iterator % contentsType
         iterator % dataType = inPool % iterator % data % contentsType
         if (iterator % memberType == MPAS_POOL_FIELD) then
            iterator % nDims = inPool % iterator % data % contentsDims
            iterator % nTimeLevels = inPool % iterator % data % contentsTimeLevs
         else
            iterator % nDims = 0
            iterator % nTimeLevels = 0
         end if
         mpas_pool_get_next_member = .true.

         ! Advance iterator to next item
         inPool % iterator => inPool % iterator % next

         ! We may have reached the end of list for current head pointer
         if (.not. associated(inPool % iterator)) then
            do i=inPool % iteratorIndex+1, inPool % size
               if (associated(inPool % table(i) % head)) exit
            end do
            inPool % iteratorIndex = i
            if (i <= inPool % size) then
               inPool % iterator => inPool % table(i) % head
            else
               nullify(inPool % iterator)
            end if
         end if

      else
         mpas_pool_get_next_member = .false.
      end if

   end function mpas_pool_get_next_member!}}}


!-----------------------------------------------------------------------
!  subroutine mpas_pool_shift_time_levels
!
!> \brief MPAS Pool Time level shift routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine shifts the time levels of all multi-level fields contained within.
!>   When shifting, time level 1 becomes time level n, and time level i becomes time level i-1.
!
!-----------------------------------------------------------------------
   recursive subroutine mpas_pool_shift_time_levels(inPool)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool


      integer :: i, j
      type (mpas_pool_member_type), pointer :: ptr
      type (mpas_pool_data_type), pointer :: dptr


      do i=1,inPool % size

         ptr => inPool % table(i) % head
         do while(associated(ptr))

            if (ptr % contentsType == MPAS_POOL_FIELD) then

               dptr => ptr % data

               if (associated(dptr % r0a)) then
                  call mpas_shift_time_levs(dptr % r0a)
               else if (associated(dptr % r1a)) then
                  call mpas_shift_time_levs(dptr % r1a)
               else if (associated(dptr % r2a)) then
                  call mpas_shift_time_levs(dptr % r2a)
               else if (associated(dptr % r3a)) then
                  call mpas_shift_time_levs(dptr % r3a)
               else if (associated(dptr % r4a)) then
                  call mpas_shift_time_levs(dptr % r4a)
               else if (associated(dptr % r5a)) then
                  call mpas_shift_time_levs(dptr % r5a)
               else if (associated(dptr % i0a)) then
                  call mpas_shift_time_levs(dptr % i0a)
               else if (associated(dptr % i1a)) then
                  call mpas_shift_time_levs(dptr % i1a)
               else if (associated(dptr % i2a)) then
                  call mpas_shift_time_levs(dptr % i2a)
               else if (associated(dptr % i3a)) then
                  call mpas_shift_time_levs(dptr % i3a)
               else if (associated(dptr % c0a)) then
                  call mpas_shift_time_levs(dptr % c0a)
               else if (associated(dptr % c1a)) then
                  call mpas_shift_time_levs(dptr % c1a)
               else if (associated(dptr % l0a)) then
                  call mpas_shift_time_levs(dptr % l0a)
               end if

            else if (ptr % contentsType == MPAS_POOL_SUBPOOL) then

               call mpas_pool_shift_time_levels(ptr % data % p)

            end if

            ptr => ptr % next
         end do

      end do

   end subroutine mpas_pool_shift_time_levels!}}}


!!!!!!!!!! Private subroutines !!!!!!!!!!

   logical function pool_add_member(inPool, key, newmem)

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      type (mpas_pool_member_type), pointer :: newmem

      integer :: hash
      type (mpas_pool_member_type), pointer :: ptr


      call pool_hash(hash, newmem % key, newmem % keylen)

      hash = mod(hash, inPool % size) + 1

      pool_add_member = .true.

      if (.not. associated(inPool % table(hash) % head)) then
         inPool % table(hash) % head => newmem
      else
         ptr => inPool % table(hash) % head
         do while (associated(ptr % next))

! TODO: should we never add a duplicate member, even if error level is set to silent?
            if (currentErrorLevel == MPAS_POOL_WARN .or. &
                currentErrorLevel == MPAS_POOL_FATAL) then
               if (ptr % contentsType == newmem % contentsType .and. &
                   ptr % keyLen == newmem % keyLen) then
                  if (ptr % key(1:ptr%keyLen) == newmem % key(1:newmem%keyLen)) then
                     pool_add_member = .false.
                     call pool_mesg('Error: Field '//trim(key)//' already exists in pool.')
                     return
                  end if
               end if
            end if

            ptr => ptr % next
         end do
         ptr % next => newmem
      end if

   end function pool_add_member


   function pool_get_member(inPool, key, memType)

      implicit none

      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in) :: memType

      type (mpas_pool_data_type), pointer :: pool_get_member

      integer :: hash, endl
      type (mpas_pool_member_type), pointer :: ptr


      nullify(pool_get_member)

      endl = len_trim(key)
      call pool_hash(hash, key, endl)

      hash = mod(hash, inPool % size) + 1

      ptr => inPool % table(hash) % head
      do while (associated(ptr))
         if (ptr % contentsType == memType) then
            if (endl == ptr % keyLen) then
               if (key(1:endl) == ptr % key(1:endl)) then
                  pool_get_member => ptr % data
                  exit
               end if
            end if
         end if
         ptr => ptr % next
      end do

   end function pool_get_member


   logical function pool_remove_member(inPool, key, memType)!{{{

      implicit none

      type (mpas_pool_type), intent(inout) :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in) :: memType

      integer :: hash, endl
      type (mpas_pool_member_type), pointer :: ptr, ptr_prev


      endl = len_trim(key)
      call pool_hash(hash, key, endl)

      hash = mod(hash, inPool % size) + 1

      if (associated(inPool % table(hash) % head)) then

         ! Is the sub-pool at the head of the list?
         ptr_prev => inPool % table(hash) % head
         if (ptr_prev % contentsType == memType) then
            if (endl == ptr_prev % keyLen) then
               if (key(1:endl) == ptr_prev % key(1:endl)) then
                  inPool % table(hash) % head => ptr_prev % next
!TODO: are there cases where we need to delete more data here?
                  deallocate(ptr_prev)
                  pool_remove_member = .true.
                  return
               end if
            end if
         end if

         ! Possibly later in the list?
         ptr => ptr_prev % next
         do while (associated(ptr))
            if (ptr % contentsType == memType) then
               if (endl == ptr % keyLen) then
                  if (key(1:endl) == ptr % key(1:endl)) then
                     ptr_prev % next => ptr % next
!TODO: are there cases where we need to delete more data here?
                     deallocate(ptr)
                     pool_remove_member = .true.
                     return
                  end if
               end if
            end if
            ptr => ptr % next
            ptr_prev => ptr_prev % next
         end do

      end if

      pool_remove_member = .false.

   end function pool_remove_member!}}}


   subroutine pool_mesg(mesg)

      implicit none

      character (len=*), intent(in) :: mesg

      if (currentErrorLevel == MPAS_POOL_WARN) then
         write(stderrUnit,*) trim(mesg)
      else if (currentErrorLevel == MPAS_POOL_FATAL) then
!TODO: using dmpar routines here will create circular module dependencies...
!         call mpas_dmpar_global_abort(trim(mesg))
      end if

   end subroutine pool_mesg


   subroutine pool_print_table_size(pool)

      implicit none

      type (mpas_pool_type), intent(in) :: pool

      integer :: i, head_size, total_size
      type (mpas_pool_member_type), pointer :: ptr


      total_size = 0
      do i=1,pool % size
         head_size = 0
         ptr => pool % table(i) % head
         do while (associated(ptr))
            head_size = head_size + 1
            ptr => ptr % next
         end do
         write(stderrUnit,*) 'List ', i, ' : ', head_size
         total_size = total_size + head_size
      end do
      write(stderrUnit,*) '----------------'
      write(stderrUnit,*) 'Total: ', total_size

   end subroutine pool_print_table_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           END POOL SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
