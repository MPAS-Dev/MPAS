! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

#define COMMA ,
#define FILTER_DEBUG_WRITE(M) ! write(stderrUnit, *) M
#define FILTER_WARN_WRITE(M) write(stderrUnit, *) 'WARNING: '//M
#define FILTER_ERROR_WRITE(M) write(stderrUnit, *) 'ERROR: '//M

!***********************************************************************
!
!  mpas_filter_list
!
!> \brief   MPAS Filter List Module
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This module defines the routines that can be used to create, use, destroy,
!> and query an mpas_filter_list type.
!
!-----------------------------------------------------------------------
module mpas_filter_list
   use mpas_kind_types
   use mpas_derived_types
   use mpas_dmpar
   use mpas_io_units
   use mpas_pool_routines

   implicit none
   private

   public :: mpas_filter_list_create_list, mpas_filter_list_destroy_list
   public :: mpas_filter_list_register_filter, mpas_filter_list_remove_filter
   public :: mpas_filter_list_get_filter_init, mpas_filter_list_get_filter_compute
   public :: mpas_filter_list_get_filter_restart, mpas_filter_list_get_filter_finalize
   public :: mpas_filter_list_begin_iteration, mpas_filter_list_get_next_filter, mpas_filter_list_get_filter_names

   contains

!***********************************************************************
!
!  routine mpas_filter_list_create_list
!
!> \brief   MPAS Filter List creation routine
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This routine will create a filter list.
!
!-----------------------------------------------------------------------
   subroutine mpas_filter_list_create_list(inList)!{{{
      type (mpas_filter_list_type), pointer :: inList

      inList % nameLength = -1
      inList % filterName = 'HEAD'
   end subroutine mpas_filter_list_create_list!}}}

!***********************************************************************
!
!  routine mpas_filter_list_destroy_list
!
!> \brief   MPAS Filter List destruction routine
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This routine will destroy a filter list.
!
!-----------------------------------------------------------------------
   subroutine mpas_filter_list_destroy_list(inList)!{{{
      type (mpas_filter_list_type), pointer :: inList

      type (mpas_filter_list_type), pointer :: listCursor

      do while (associated(inList % next))
         listCursor => inList % next

         if ( associated(listCursor % next) ) then
            inList % next => listCursor % next
         else
            nullify(inList % next)
         end if

         deallocate(listCursor)
      end do
   end subroutine mpas_filter_list_destroy_list!}}}

!***********************************************************************
!
!  routine mpas_filter_list_register_filter
!
!> \brief   MPAS Filter list filter registration routine
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This routine will add a new filter to a filter list.
!
!-----------------------------------------------------------------------
   subroutine mpas_filter_list_register_filter(inList, filterName, initFunc, computeFunc, restartFunc, finalizeFunc)!{{{
      type (mpas_filter_list_type), pointer :: inList
      character (len=*), intent(in) :: filterName
      procedure (mpas_filter_list_init_function), pointer :: initFunc
      procedure (mpas_filter_list_compute_function), pointer :: computeFunc
      procedure (mpas_filter_list_restart_function), pointer :: restartFunc
      procedure (mpas_filter_list_finalize_function), pointer :: finalizeFunc

      type (mpas_filter_list_type), pointer :: listCursor
      type (mpas_filter_list_type), pointer :: listCursorPrev
      integer :: nameLength

      nameLength = len_trim(filterName)
      listCursor => inList
      nullify(listCursorPrev)

      do while (associated(listCursor))
         if ( listCursor % nameLength == nameLength ) then
            if ( listCursor % filterName(1:nameLength) == filterName(1:nameLength) ) then
               call mpas_dmpar_global_abort('ERROR: Filter '''//filterName(1:nameLength)//''' already exists in filter list. Exiting...')
            end if
         end if
       
         listCursorPrev => listCursor
         listCursor => listCursor % next
      end do

      listCursor => listCursorPrev

      allocate(listCursor % next)
      listCursor => listCursor % next
      listCursor % nameLength = nameLength
      listCursor % filterName = filterName
      listCursor % initFunc => initFunc
      listCursor % computeFunc => computeFunc
      listCursor % restartFunc => restartFunc
      listCursor % finalizeFunc => finalizeFunc
      nullify(listCursor % next)

   end subroutine mpas_filter_list_register_filter!}}}

!***********************************************************************
!
!  routine mpas_filter_list_remove_filter
!
!> \brief   MPAS Filter list filter removal routine
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This routine will remove a filter from a filter list.
!
!-----------------------------------------------------------------------
   subroutine mpas_filter_list_remove_filter(inList, filterName)!{{{
      type (mpas_filter_list_type), pointer :: inList
      character (len=*), intent(in) :: filterName

      type (mpas_filter_list_type), pointer :: listCursor
      type (mpas_filter_list_type), pointer :: listCursorPrev
      integer :: nameLength

      nameLength = len_trim(filterName)

      listCursor => inList
      nullify(listCursorPrev)

      do while (associated(listCursor))
         if ( listCursor % nameLength == nameLength ) then
            if ( listCursor % filterName(1:nameLength) == filterName(1:nameLength) ) then
               if ( associated(listCursor % next) ) then
                  listCursorPrev % next => listCursor % next
                  deallocate(listCursor)
                  FILTER_DEBUG_WRITE(' -- Removing filter '''//filterName(1:nameLength)//''' from filter list.')
                  return
               end if
            end if
         end if
       
         listCursorPrev => listCursor
         listCursor => listCursor % next
      end do

      call mpas_dmpar_global_abort('ERROR: Filter '''//filterName(1:nameLength)//''' does not exist in filter list to remove. Exiting...')
   end subroutine mpas_filter_list_remove_filter!}}}

!***********************************************************************
!
!  routine mpas_filter_list_get_filter_init
!
!> \brief   MPAS Filter list filter init routine query routine
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This routine will find a filter in a filter list, and return a pointer to
!> the init function for that filter.
!
!-----------------------------------------------------------------------
   subroutine mpas_filter_list_get_filter_init(inList, filterName, outFunc, iErr)!{{{
      type (mpas_filter_list_type), pointer :: inList
      character (len=*), intent(in) :: filterName
      procedure (mpas_filter_list_init_function), pointer :: outFunc
      integer, intent(out) :: iErr

      type (mpas_filter_list_type), pointer :: listCursor
      integer :: nameLength

      iErr = MPAS_FILTER_NOERR

      nameLength = len_trim(filterName)

      listCursor => inList
      do while (associated(listCursor))
         if ( listCursor % nameLength == nameLength ) then
            if ( listCursor % filterName(1:nameLength) == filterName ) then
               FILTER_DEBUG_WRITE(' -- Found filter '''//filterName(1:nameLength)//''' returning init function pointer.')
               outFunc => listCursor % initFunc
               return
            end if
         end if
         listCursor => listCursor % next
      end do

      FILTER_ERROR_WRITE('ERROR: Filter '''//filterName(1:nameLength)//''' not found in filter list.')
      nullify(outFunc)
      iErr = MPAS_FILTER_NOTFOUND
   end subroutine mpas_filter_list_get_filter_init!}}}

!***********************************************************************
!
!  routine mpas_filter_list_get_filter_compute
!
!> \brief   MPAS Filter list filter compute routine query routine
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This routine will find a filter in a filter list, and return a pointer to
!> the compute function for that filter.
!
!-----------------------------------------------------------------------
   subroutine mpas_filter_list_get_filter_compute(inList, filterName, outFunc, iErr)!{{{
      type (mpas_filter_list_type), pointer :: inList
      character (len=*), intent(in) :: filterName
      procedure (mpas_filter_list_compute_function), pointer :: outFunc
      integer, intent(out) :: iErr

      type (mpas_filter_list_type), pointer :: listCursor
      integer :: nameLength

      iErr = MPAS_FILTER_NOERR

      nameLength = len_trim(filterName)

      listCursor => inList
      do while (associated(listCursor))
         if ( listCursor % nameLength == nameLength ) then
            if ( listCursor % filterName(1:nameLength) == filterName ) then
               FILTER_DEBUG_WRITE(' -- Found filter '''//filterName(1:nameLength)//''' returning compute function pointer.')
               outFunc => listCursor % computeFunc
               return
            end if
         end if
         listCursor => listCursor % next
      end do

      FILTER_ERROR_WRITE('ERROR: Filter '''//filterName(1:nameLength)//''' not found in filter list.')
      nullify(outFunc)
      iErr = MPAS_FILTER_NOTFOUND
   end subroutine mpas_filter_list_get_filter_compute!}}}

!***********************************************************************
!
!  routine mpas_filter_list_get_filter_restart
!
!> \brief   MPAS Filter list filter restart routine query routine
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This routine will find a filter in a filter list, and return a pointer to
!> the restart function for that filter.
!
!-----------------------------------------------------------------------
   subroutine mpas_filter_list_get_filter_restart(inList, filterName, outFunc, iErr)!{{{
      type (mpas_filter_list_type), pointer :: inList
      character (len=*), intent(in) :: filterName
      procedure (mpas_filter_list_restart_function), pointer :: outFunc
      integer, intent(out) :: iErr

      type (mpas_filter_list_type), pointer :: listCursor
      integer :: nameLength

      iErr = MPAS_FILTER_NOERR

      nameLength = len_trim(filterName)

      listCursor => inList
      do while (associated(listCursor))
         if ( listCursor % nameLength == nameLength ) then
            if ( listCursor % filterName(1:nameLength) == filterName ) then
               FILTER_DEBUG_WRITE(' -- Found filter '''//filterName(1:nameLength)//''' returning restart function pointer.')
               outFunc => listCursor % restartFunc
               return
            end if
         end if
         listCursor => listCursor % next
      end do

      FILTER_ERROR_WRITE('ERROR: Filter '''//filterName(1:nameLength)//''' not found in filter list.')
      nullify(outFunc)
      iErr = MPAS_FILTER_NOTFOUND
   end subroutine mpas_filter_list_get_filter_restart!}}}

!***********************************************************************
!
!  routine mpas_filter_list_get_filter_finalize
!
!> \brief   MPAS Filter list filter finalize routine query routine
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This routine will find a filter in a filter list, and return a pointer to
!> the finalize function for that filter.
!
!-----------------------------------------------------------------------
   subroutine mpas_filter_list_get_filter_finalize(inList, filterName, outFunc, iErr)!{{{
      type (mpas_filter_list_type), pointer :: inList
      character (len=*), intent(in) :: filterName
      procedure (mpas_filter_list_finalize_function), pointer :: outFunc
      integer, intent(out) :: iErr

      type (mpas_filter_list_type), pointer :: listCursor
      integer :: nameLength

      iErr = MPAS_FILTER_NOERR

      nameLength = len_trim(filterName)

      listCursor => inList
      do while (associated(listCursor))
         if ( listCursor % nameLength == nameLength ) then
            if ( listCursor % filterName(1:nameLength) == filterName ) then
               FILTER_DEBUG_WRITE(' -- Found filter '''//filterName(1:nameLength)//''' returning finalize function pointer.')
               outFunc => listCursor % finalizeFunc
               return
            end if
         end if
         listCursor => listCursor % next
      end do

      FILTER_ERROR_WRITE('ERROR: Filter '''//filterName(1:nameLength)//''' not found in filter list.')
      nullify(outFunc)
      iErr = MPAS_FILTER_NOTFOUND
   end subroutine mpas_filter_list_get_filter_finalize!}}}


!***********************************************************************
!
!  routine mpas_filter_list_begin_iteration
!
!> \brief   MPAS Filter list iteration begin routine
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This routine begins the iteration process for an mpas_filter_list_type.
!
!-----------------------------------------------------------------------
   subroutine mpas_filter_list_begin_iteration(inList)!{{{
      type (mpas_filter_list_type), pointer :: inList
      
      if ( associated(inList % next) ) then
         inList % currentFilter => inList % next
      else
         nullify(inList % currentFilter)
      end if
   end subroutine mpas_filter_list_begin_iteration!}}}


!***********************************************************************
!
!  function mpas_filter_list_get_next_filter
!
!> \brief   MPAS Filter list get next filter function
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This function takes as input a filter list, and returns a logical stating if
!> the current filter in the iteration process is valid or not. A return value
!> of false means there are not more filters to iterate through, while a return
!> value of true means the current filter is valid.
!> 
!> Upon returning, there are optional arguments for the different functions a
!> filter can provide, to easily query these while iterating through the list.
!
!-----------------------------------------------------------------------
   function mpas_filter_list_get_next_filter(inList, filterName, initFunc, computeFunc, restartFunc, finalizeFunc) result(validFilter)!{{{
      type (mpas_filter_list_type), pointer :: inList
      character (len=StrKIND), intent(out) :: filterName
      procedure (mpas_filter_list_init_function), pointer, optional :: initFunc
      procedure (mpas_filter_list_compute_function), pointer, optional :: computeFunc
      procedure (mpas_filter_list_restart_function), pointer, optional :: restartFunc
      procedure (mpas_filter_list_finalize_function), pointer, optional :: finalizeFunc
      logical :: validFilter

      validFilter = .false.
      filterName = 'garbage'
      if ( present(initFunc) ) then
         nullify(initFunc)
      end if

      if ( present(computeFunc) ) then
         nullify(computeFunc)
      end if

      if ( present(restartFunc) ) then
         nullify(restartFunc)
      end if

      if ( present(finalizeFunc) ) then
         nullify(finalizeFunc)
      end if

      if ( associated(inList % currentFilter) ) then
         validFilter = .true.
         filterName = trim(inList % currentFilter % filterName)

         if ( present(initFunc) ) then
            initFunc => inList % currentFilter % initFunc
         end if

         if ( present(computeFunc) ) then
            computeFunc => inList % currentFilter % computeFunc
         end if

         if ( present(restartFunc) ) then
            restartFunc => inList % currentFilter % restartFunc
         end if

         if ( present(finalizeFunc) ) then
            finalizeFunc => inList % currentFilter % finalizeFunc
         end if

         if ( associated(inList % currentFilter % next) ) then
            inList % currentFilter => inList % currentFilter % next
         else
            nullify(inList % currentFilter)
         end if
      end if
   end function mpas_filter_list_get_next_filter!}}}

!***********************************************************************
!
!  routine mpas_filter_list_get_filter_names
!
!> \brief   MPAS Filter list get filter names routine
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This routine will take a filter list as input, and add their names as configs to 
!> the pool passed in as outPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_filter_list_get_filter_names(inList, outPool)!{{{
      type (mpas_filter_list_type), pointer :: inList
      type (mpas_pool_type), pointer :: outPool

      character (len=StrKIND) :: filterName
      
      call mpas_filter_list_begin_iteration(inList)
      do while ( mpas_filter_list_get_next_filter(inList, filterName) )
         call mpas_pool_add_config(outPool, filterName, 1)
      end do

   end subroutine mpas_filter_list_get_filter_names!}}}

end module mpas_filter_list
