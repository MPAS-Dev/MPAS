   integer, parameter :: MPAS_FILTER_NOERR = 0, &
                         MPAS_FILTER_NOTFOUND = 1

   abstract interface
      function mpas_filter_list_init_function(domain, filterName) result(iErr)
         import domain_type

         type (domain_type), intent(inout) :: domain
         character (len=*), intent(in) :: filterName
         integer :: iErr
      end function
   end interface

   abstract interface
      function mpas_filter_list_compute_function(domain, filterName, timeLevel) result(iErr)
         import domain_type

         type (domain_type), intent(inout) :: domain
         character (len=*), intent(in) :: filterName
         integer, intent(in) :: timeLevel
         integer :: iErr
      end function
   end interface

   abstract interface
      function mpas_filter_list_restart_function(domain, filterName) result(iErr)
         import domain_type

         type (domain_type), intent(inout) :: domain
         character (len=*), intent(in) :: filterName
         integer :: iErr
      end function
   end interface

   abstract interface
      function mpas_filter_list_finalize_function(domain, filterName) result(iErr)
         import domain_type

         type (domain_type), intent(inout) :: domain
         character (len=*), intent(in) :: filterName
         integer :: iErr
      end function
   end interface

   type mpas_filter_list_type
      integer :: nameLength
      character (len=StrKIND) :: filterName
      procedure (mpas_filter_list_init_function), pointer, nopass :: initFunc => null()
      procedure (mpas_filter_list_compute_function), pointer, nopass :: computeFunc => null()
      procedure (mpas_filter_list_restart_function), pointer, nopass :: restartFunc => null()
      procedure (mpas_filter_list_finalize_function), pointer, nopass :: finalizeFunc => null()

	  type (mpas_filter_list_type), pointer :: currentFilter
	  type (mpas_filter_list_type), pointer :: next
   end type mpas_filter_list_type
