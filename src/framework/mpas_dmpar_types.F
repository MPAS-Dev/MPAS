!-----------------------------------------------------------------------
!  mpas_dmpar_types
!
!> \brief MPAS Communication Type Definitions
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This module defines all communication related derived data types
!
!-----------------------------------------------------------------------
module mpas_dmpar_types

   use mpas_kind_types

   type dm_info
     integer :: nprocs, my_proc_id, comm, info
     logical :: using_external_comm
   end type dm_info


   type exchange_list
     integer :: procID
     integer :: blockID
     integer :: nlist
     integer, dimension(:), pointer :: list
     type (mpas_exchange_list), pointer :: next
     real (kind=RKIND), dimension(:), pointer :: rbuffer
     integer, dimension(:), pointer           :: ibuffer
     integer :: reqID

   end type exchange_list

   type mpas_exchange_list
     integer :: endPointID
     integer :: nlist
     integer, dimension(:), pointer :: srcList
     integer, dimension(:), pointer :: destList
     type (mpas_exchange_list), pointer :: next

   end type mpas_exchange_list

   type mpas_exchange_list_pointer
     type (mpas_exchange_list), pointer :: exchList
   end type mpas_exchange_list_pointer

   type mpas_multihalo_exchange_list
     type (mpas_exchange_list_pointer), dimension(:), pointer :: halos
   end type mpas_multihalo_exchange_list


   type mpas_communication_list
     integer :: procID
     integer :: nlist
     real (kind=RKIND), dimension(:), pointer :: rbuffer
     integer, dimension(:), pointer :: ibuffer
     integer :: reqID
     type (mpas_communication_list), pointer :: next

   end type mpas_communication_list


end module mpas_dmpar_types
