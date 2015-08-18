! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_filter_test_compute_interval
!
!> \brief   MPAS Test compute interval filter module
!> \author  Doug Jacobsen
!> \date    07/11/2015
!> \details
!> This module defines a filter that can be used to test the compute interval
!> portion of the filter infrastructure.
!
!-----------------------------------------------------------------------
module mpas_filter_test_compute_interval
   use mpas_io_units
   use mpas_kind_types
   use mpas_derived_types
   use mpas_stream_manager
   use mpas_pool_routines

   implicit none
   private

   public :: mpas_filter_test_compute_interval_init
   public :: mpas_filter_test_compute_interval_compute
   public :: mpas_filter_test_compute_interval_restart
   public :: mpas_filter_test_compute_interval_finalize

   contains

   function mpas_filter_test_compute_interval_init(domain, filterName) result(iErr)!{{{
      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: filterName
      integer :: iErr

      character (len=StrKIND) :: configName, fieldName
      character (len=StrKIND), pointer :: filterStreamName

      type (block_type), pointer :: block_ptr
      type (field0DInteger), pointer :: computeIndex, computeIndexPrev
      integer :: nameLength

      integer :: err_local

      iErr = 0

      nameLength = len_trim(filterName)

      configName = 'config_' // filterName(1:nameLength) // '_stream_name'
      call mpas_pool_get_config(domain % configs, configName, filterStreamName)

      if ( .not. associated(filterStreamName) ) then
         write(stderrUnit, *) 'ERROR: Namelist option ''' // trim(configName) // ''' not found for filter ''' // filterName(1:nameLength) // ''''
         iErr = 1
         return
      end if

      nullify(computeIndexPrev)

      ! Create field for counting the number of computations, and add it to the allFields pool
      fieldName = 'computeIndex_' // filterName(1:nameLength)
      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         allocate(computeIndex)

         computeIndex % block => block_ptr
         computeIndex % fieldName = 'computeIndex'
         computeIndex % isDecomposed = .false.
         computeIndex % hasTimeDimension = .true.
         computeIndex % isActive = .true.
         computeIndex % isVarArray = .false.
         computeIndex % scalar = 0
         nullify(computeIndex % sendList)
         nullify(computeIndex % recvList)
         nullify(computeIndex % copyList)

         if ( associated(computeIndexPrev) ) then
            computeIndex % prev => computeIndexPrev
            computeIndexPrev % next => computeIndex
         end if

         call mpas_pool_add_field(block_ptr % allFields, fieldName, computeIndex)

         computeIndexPrev => computeIndex
         block_ptr => block_ptr % next
      end do

      ! Add the computeIndex field to the stream
      call mpas_stream_mgr_add_field(domain % streamManager, streamID=filterStreamName, fieldName=fieldName, iErr=err_local)

   end function mpas_filter_test_compute_interval_init!}}}

   function mpas_filter_test_compute_interval_compute(domain, filterName, timeLevel) result(iErr)!{{{
      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: filterName
      integer, intent(in), optional :: timeLevel
      integer :: iErr

      character (len=StrKIND) :: fieldName
      type (block_type), pointer :: block_ptr
      integer, pointer :: computeIndex

      iErr = 0

      fieldName = 'computeIndex_' // trim(filterName)

      block_ptr => domain % blocklist
      do while ( associated(block_ptr) )
         call mpas_pool_get_array(block_ptr % allFields, fieldName, computeIndex)

         computeIndex = computeIndex + 1

         block_ptr => block_ptr % next
      end do
   end function mpas_filter_test_compute_interval_compute!}}}

   function mpas_filter_test_compute_interval_restart(domain, filterName) result(iErr)!{{{
      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: filterName
      integer :: iErr

      iErr = 0
   end function mpas_filter_test_compute_interval_restart!}}}

   function mpas_filter_test_compute_interval_finalize(domain, filterName) result(iErr)!{{{
      type (domain_type), intent(inout) :: domain
      character (len=*), intent(in) :: filterName
      integer :: iErr

      iErr = 0
   end function mpas_filter_test_compute_interval_finalize!}}}

end module mpas_filter_test_compute_interval
