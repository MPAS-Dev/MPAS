! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_grid_types
!
!> \brief   MPAS Grid and field type defintion module
!> \author  Michael Duda, Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This module defines derived data types related to fields, and variable structures.
!> It also includes routines for allocating and deallocating these types.
!
!-----------------------------------------------------------------------
module mpas_grid_types

   use mpas_kind_types
   use mpas_dmpar_types
   use mpas_attlist
   use mpas_packages

   integer, parameter :: nTimeLevs = 2

   ! Derived type describing info for doing I/O specific to a field
   type io_info
      character (len=StrKIND) :: fieldName
      character (len=StrKIND) :: units
      character (len=StrKIND) :: description
      integer, dimension(4) :: start
      integer, dimension(4) :: count
      logical :: input
      logical :: sfc
      logical :: restart
      logical :: output
   end type io_info

   ! Derived type for storing fields
   type field5DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(5) :: dimNames
      integer, dimension(5) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field5DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field5DReal


   ! Derived type for storing fields
   type field4DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(4) :: dimNames
      integer, dimension(4) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field4DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field4DReal



   ! Derived type for storing fields
   type field3DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(3) :: dimNames
      integer, dimension(3) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field3DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field3DReal


   ! Derived type for storing fields
   type field2DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(2) :: dimNames
      integer, dimension(2) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field2DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field2DReal


   ! Derived type for storing fields
   type field1DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DReal


   ! Derived type for storing fields
   type field0DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND) :: scalar

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DReal


   ! Derived type for storing fields
   type field3DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:,:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(3) :: dimNames
      integer, dimension(3) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field3DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field3DInteger


   ! Derived type for storing fields
   type field2DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(2) :: dimNames
      integer, dimension(2) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field2DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field2DInteger


   ! Derived type for storing fields
   type field1DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DInteger


   ! Derived type for storing fields
   type field0DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer :: scalar

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DInteger


   ! Derived type for storing fields
   type field1DChar
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      character (len=StrKIND), dimension(:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DChar), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DChar


   ! Derived type for storing fields
   type field0DChar
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      character (len=StrKIND) :: scalar

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DChar), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DChar


   ! Derived type for storing fields
   type field0DLogical
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      logical :: scalar

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DLogical), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DLogical

   !!! TYPES FOR MPAS_POOL
   integer, parameter :: MPAS_POOL_SILENT = 1, &
                         MPAS_POOL_WARN = 2, &
                         MPAS_POOL_FATAL = 3

   integer, parameter :: MPAS_POOL_FIELD = 1, &
                         MPAS_POOL_CONFIG = 2, &
                         MPAS_POOL_DIMENSION = 3, &
                         MPAS_POOL_SUBPOOL = 4, &
                         MPAS_POOL_PACKAGE = 5

   integer, parameter :: MPAS_POOL_REAL = 1, &
                         MPAS_POOL_INTEGER = 2, &
                         MPAS_POOL_LOGICAL = 3, &
                         MPAS_POOL_CHARACTER = 4

   type mpas_pool_data_type
      integer :: contentsType
      integer :: contentsDims
      integer :: contentsTimeLevs
      type (field0DReal), dimension(:), pointer :: r0 => null()
      type (field1DReal), dimension(:), pointer :: r1 => null()
      type (field2DReal), dimension(:), pointer :: r2 => null()
      type (field3DReal), dimension(:), pointer :: r3 => null()
      type (field4DReal), dimension(:), pointer :: r4 => null()
      type (field5DReal), dimension(:), pointer :: r5 => null()
      type (field0DInteger), dimension(:), pointer :: i0 => null()
      type (field1DInteger), dimension(:), pointer :: i1 => null()
      type (field2DInteger), dimension(:), pointer :: i2 => null()
      type (field3DInteger), dimension(:), pointer :: i3 => null()
      type (field0DChar), dimension(:), pointer :: c0 => null()
      type (field1DChar), dimension(:), pointer :: c1 => null()
      type (field0DLogical), dimension(:), pointer :: l0 => null()
   end type mpas_pool_data_type

   type mpas_pool_member_type
      character (len=StrKIND) :: key
      integer :: keyLen
      integer contentsType
      type (mpas_pool_data_type), pointer :: data => null()
      type (mpas_pool_type), pointer :: pool => null()
      type (mpas_pool_member_type), pointer :: next => null()
   end type mpas_pool_member_type

   type mpas_pool_head_type
      type (mpas_pool_member_type), pointer :: head => null()
   end type mpas_pool_head_type

   type mpas_pool_type
      integer :: size
      integer :: iteratorIndex
      type (mpas_pool_head_type), pointer :: head => null()
      type (mpas_pool_member_type), pointer :: iterator => null()
   end type mpas_pool_type

   type mpas_pool_iterator_type
      character (len=StrKIND) :: memberName
      integer :: memberType
      integer :: dataType
      integer :: nDims
      integer :: nTimeLevels
   end type mpas_pool_iterator_type

   type mpas_pool_field_info_type
      integer :: fieldType
      integer :: nDims
      integer :: nTimeLevels
   end type mpas_pool_field_info_type
   !!! END TYPES FOR MPAS POOL



   ! Derived type for storing grid meta-data
   type mesh_type

      type (block_type), pointer :: block

#include "field_dimensions.inc"

      logical :: on_a_sphere
      real (kind=RKIND) :: sphere_radius

#include "time_invariant_fields.inc"

   end type mesh_type


#include "variable_groups.inc"


   ! Type for storing (possibly architecture specific) information concerning to parallelism
   type parallel_info
      type (mpas_multihalo_exchange_list), pointer :: cellsToSend            ! List of types describing which cells to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: cellsToRecv            ! List of types describing which cells to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: cellsToCopy            ! List of types describing which cells to copy from other blocks

      type (mpas_multihalo_exchange_list), pointer :: edgesToSend            ! List of types describing which edges to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: edgesToRecv            ! List of types describing which edges to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: edgesToCopy            ! List of types describing which edges to copy from other blocks

      type (mpas_multihalo_exchange_list), pointer :: verticesToSend         ! List of types describing which vertices to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: verticesToRecv         ! List of types describing which vertices to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: verticesToCopy         ! List of types describing which vertices to copy from other blocks
   end type parallel_info


   ! Derived type for storing part of a domain; used as a basic unit of work for a process
   type block_type

#include "block_group_members.inc"

      integer :: blockID   ! Unique global ID number for this block
      integer :: localBlockID  ! Unique local ID number for this block

      type (domain_type), pointer :: domain

      type (parallel_info), pointer :: parinfo

      type (block_type), pointer :: prev, next

      type (mpas_pool_type), pointer :: structs, dimensions, configs, packages
   end type block_type


   ! Derived type for storing list of blocks from a domain to be handled by a process
   type domain_type
      type (block_type), pointer :: blocklist
      type (mpas_pool_type), pointer :: configs, packages
   
      ! Also store parallelization info here
      type (dm_info), pointer :: dminfo
#include "model_variables.inc"
      character (len=StrKIND*2) :: history !< History attribute, read in from input file.
      character (len=StrKIND) :: mesh_id !< mesh_id attribute, randomly generated
      character (len=StrKIND) :: Conventions !< Conventions attribute, read in from input file.
      character (len=StrKIND) :: source !< source attribute, read in from input file.
      character (len=StrKIND) :: mesh_spec !< mesh_spec attribute, read in from input file.
      character (len=StrKIND) :: parent_id !< parent_id attribute, read in from input file.
   end type domain_type

   interface mpas_allocate_scratch_field
     module procedure mpas_allocate_scratch_field1d_integer
     module procedure mpas_allocate_scratch_field2d_integer
     module procedure mpas_allocate_scratch_field3d_integer
     module procedure mpas_allocate_scratch_field1d_real
     module procedure mpas_allocate_scratch_field2d_real
     module procedure mpas_allocate_scratch_field3d_real
     module procedure mpas_allocate_scratch_field4d_real
     module procedure mpas_allocate_scratch_field5d_real
     module procedure mpas_allocate_scratch_field1d_char
   end interface

   interface mpas_deallocate_scratch_field
     module procedure mpas_deallocate_scratch_field1d_integer
     module procedure mpas_deallocate_scratch_field2d_integer
     module procedure mpas_deallocate_scratch_field3d_integer
     module procedure mpas_deallocate_scratch_field1d_real
     module procedure mpas_deallocate_scratch_field2d_real
     module procedure mpas_deallocate_scratch_field3d_real
     module procedure mpas_deallocate_scratch_field4d_real
     module procedure mpas_deallocate_scratch_field5d_real
     module procedure mpas_deallocate_scratch_field1d_char
   end interface

   interface mpas_deallocate_field
     module procedure mpas_deallocate_field0d_integer
     module procedure mpas_deallocate_field1d_integer
     module procedure mpas_deallocate_field2d_integer
     module procedure mpas_deallocate_field3d_integer
     module procedure mpas_deallocate_field0d_real
     module procedure mpas_deallocate_field1d_real
     module procedure mpas_deallocate_field2d_real
     module procedure mpas_deallocate_field3d_real
     module procedure mpas_deallocate_field4d_real
     module procedure mpas_deallocate_field5d_real
     module procedure mpas_deallocate_field0d_char
     module procedure mpas_deallocate_field1d_char
   end interface

   interface mpas_pool_add_field
      module procedure mpas_pool_add_field_0d_real
      module procedure mpas_pool_add_field_1d_real
      module procedure mpas_pool_add_field_2d_real
      module procedure mpas_pool_add_field_3d_real
      module procedure mpas_pool_add_field_4d_real
      module procedure mpas_pool_add_field_5d_real
      module procedure mpas_pool_add_field_0d_int
      module procedure mpas_pool_add_field_1d_int
      module procedure mpas_pool_add_field_2d_int
      module procedure mpas_pool_add_field_3d_int
      module procedure mpas_pool_add_field_0d_char
      module procedure mpas_pool_add_field_1d_char
      module procedure mpas_pool_add_field_0d_reals
      module procedure mpas_pool_add_field_1d_reals
      module procedure mpas_pool_add_field_2d_reals
      module procedure mpas_pool_add_field_3d_reals
      module procedure mpas_pool_add_field_4d_reals
      module procedure mpas_pool_add_field_5d_reals
      module procedure mpas_pool_add_field_0d_ints
      module procedure mpas_pool_add_field_1d_ints
      module procedure mpas_pool_add_field_2d_ints
      module procedure mpas_pool_add_field_3d_ints
      module procedure mpas_pool_add_field_0d_chars
      module procedure mpas_pool_add_field_1d_chars
   end interface

   interface mpas_pool_get_field
      module procedure mpas_pool_get_field_0d_real
      module procedure mpas_pool_get_field_1d_real
      module procedure mpas_pool_get_field_2d_real
      module procedure mpas_pool_get_field_3d_real
      module procedure mpas_pool_get_field_4d_real
      module procedure mpas_pool_get_field_5d_real
      module procedure mpas_pool_get_field_0d_int
      module procedure mpas_pool_get_field_1d_int
      module procedure mpas_pool_get_field_2d_int
      module procedure mpas_pool_get_field_3d_int
      module procedure mpas_pool_get_field_0d_char
      module procedure mpas_pool_get_field_1d_char
   end interface

   interface mpas_pool_get_array
      module procedure mpas_pool_get_array_0d_real
      module procedure mpas_pool_get_array_1d_real
      module procedure mpas_pool_get_array_2d_real
      module procedure mpas_pool_get_array_3d_real
      module procedure mpas_pool_get_array_4d_real
      module procedure mpas_pool_get_array_5d_real
      module procedure mpas_pool_get_array_0d_int
      module procedure mpas_pool_get_array_1d_int
      module procedure mpas_pool_get_array_2d_int
      module procedure mpas_pool_get_array_3d_int
      module procedure mpas_pool_get_array_0d_char
      module procedure mpas_pool_get_array_1d_char
   end interface

   interface mpas_pool_add_config
      module procedure mpas_pool_add_config_real
      module procedure mpas_pool_add_config_int
      module procedure mpas_pool_add_config_char
      module procedure mpas_pool_add_config_logical
   end interface

   interface mpas_pool_get_config
      module procedure mpas_pool_get_config_real
      module procedure mpas_pool_get_config_int
      module procedure mpas_pool_get_config_char
      module procedure mpas_pool_get_config_logical
   end interface

   interface mpas_pool_add_dimension
      module procedure mpas_pool_add_dimension_0d
      module procedure mpas_pool_add_dimension_1d
   end interface

   interface mpas_pool_get_dimension
      module procedure mpas_pool_get_dimension_0d
      module procedure mpas_pool_get_dimension_1d
   end interface

   integer :: currentErrorLevel = MPAS_POOL_SILENT



   contains

!***********************************************************************
!
!  routine mpas_allocate_domain
!
!> \brief   MPAS Domain allocation routine
!> \author  Michael Duda
!> \date    04/02/13
!> \details 
!> This routine allocates a domain structure.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_domain(dom, dminfo)!{{{

      implicit none

      type (domain_type), pointer :: dom !< Input/Output: Domain structure
      type (dm_info), pointer :: dminfo !< Input: Domain Information

      nullify(dom % blocklist)
      dom % dminfo => dminfo

      allocate(dom % configs)
      allocate(dom % packages)

      call mpas_pool_create_pool(dom % configs)
      call mpas_pool_create_pool(dom % packages)

   end subroutine mpas_allocate_domain!}}}

!***********************************************************************
!
!  routine mpas_allocate_block
!
!> \brief   MPAS Block allocation routine
!> \author  Michael Duda
!> \date    04/02/13
!> \details 
!> This routine allocates a block structure. It calls routines to allocate the variable structures
!> that are members of the block type.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_block(nHaloLayers, b, dom, blockID, &!{{{
#include "dim_dummy_args.inc"
                            )

      implicit none

      integer, intent(in) :: nHaloLayers !< Input: Number of halo laters
      type (block_type), pointer :: b !< Input/Output: Block structure
      type (domain_type), pointer :: dom !< Input: Domain structure
      integer, intent(in) :: blockID !< Input: Global ID of block
#include "dim_dummy_decls.inc"

      integer :: i

      b % blockID = blockID

      allocate(b % parinfo)

      b % domain => dom

      allocate(b % structs)
      allocate(b % dimensions)
      call mpas_pool_create_pool(b % structs)
      call mpas_pool_create_pool(b % dimensions)

      b % configs => dom % configs
      b % packages => dom % packages


#include "block_allocs.inc"

   end subroutine mpas_allocate_block!}}}


#include "group_alloc_routines.inc"

#include "provis_alloc_routines.inc"

!***********************************************************************
!
!  routine mpas_deallocate_domain
!
!> \brief   MPAS Domain deallocation routine
!> \author  Michael Duda
!> \date    04/02/13
!> \details 
!> This routine deallocates a domain structure. 
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_domain(dom)!{{{

      implicit none

      type (domain_type), pointer :: dom !< Input/Output: Domain to deallocate

      type (block_type), pointer :: block_ptr

      block_ptr => dom % blocklist
      do while (associated(block_ptr))
         call mpas_deallocate_block(block_ptr)
         block_ptr => block_ptr % next
      end do

   end subroutine mpas_deallocate_domain!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field1d_integer
!
!> \brief   MPAS 1D Scratch integer allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 1D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field1d_integer(f, single_block_in)!{{{
       type (field1dInteger), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated or all blocks.
       logical :: single_block
       type (field1dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field1d_integer!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field2d_integer
!
!> \brief   MPAS 2D Scratch integer allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 2D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field2d_integer(f, single_block_in)!{{{
       type (field2dInteger), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field2dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field2d_integer!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field3d_integer
!
!> \brief   MPAS 3D Scratch integer allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 3D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field3d_integer(f, single_block_in)!{{{
       type (field3dInteger), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field3dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field3d_integer!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field1d_real
!
!> \brief   MPAS 1D Scratch real allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 1D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field1d_real(f, single_block_in)!{{{
       type (field1dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field1dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field1d_real!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field2d_real
!
!> \brief   MPAS 2D Scratch real allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 2D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field2d_real(f, single_block_in)!{{{
       type (field2dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field2dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field2d_real!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field3d_real
!
!> \brief   MPAS 3D Scratch real allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 3D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field3d_real(f, single_block_in)!{{{
       type (field3dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field3dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field3d_real!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field4D_real
!
!> \brief   MPAS 4D Scratch real allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 4D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field4d_real(f, single_block_in)!{{{
       type (field4dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field4dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3), f_cursor % dimSizes(4)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3), f % dimSizes(4)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field4d_real!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field5D_real
!
!> \brief   MPAS 5D Scratch real allocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 5D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field5d_real(f, single_block_in)!{{{
       type (field5dReal), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field5dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3), f_cursor % dimSizes(4), f_cursor % dimSizes(5)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3), f % dimSizes(4), f % dimSizes(5)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field5d_real!}}}

!***********************************************************************
!
!  routine mpas_allocate_scratch_field1D_char
!
!> \brief   MPAS 1D Scratch character deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine allocates a 1D scratch character field.
!
!-----------------------------------------------------------------------
   subroutine mpas_allocate_scratch_field1d_char(f, single_block_in)!{{{
       type (field1dChar), pointer :: f !< Input: Field to allocate
       logical, intent(in), optional :: single_block_in !< Input: Logical flag that determines if a single block should be allocated, or all blocks.
       logical :: single_block
       type (field1dChar), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field1d_char!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field1D_integer
!
!> \brief   MPAS 1D Scratch integer deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field1d_integer(f, single_block_in)!{{{
       type (field1dInteger), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field1dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field1d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field2D_integer
!
!> \brief   MPAS 2D Scratch integer deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 2D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field2d_integer(f, single_block_in)!{{{
       type (field2dInteger), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field2dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field2d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field3D_integer
!
!> \brief   MPAS 3D Scratch integer deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 3D scratch integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field3d_integer(f, single_block_in)!{{{
       type (field3dInteger), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field3dInteger), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field3d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field1D_real
!
!> \brief   MPAS 1D Scratch real deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field1d_real(f, single_block_in)!{{{
       type (field1dReal), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field1dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field1d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field2D_real
!
!> \brief   MPAS 2D Scratch real deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 2D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field2d_real(f, single_block_in)!{{{
       type (field2dReal), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field2dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field2d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field3D_real
!
!> \brief   MPAS 3D Scratch real deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 3D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field3d_real(f, single_block_in)!{{{
       type (field3dReal), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field3dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field3d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field4D_real
!
!> \brief   MPAS 4D Scratch real deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 4D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field4d_real(f, single_block_in)!{{{
       type (field4dReal), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field4dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field4d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field5D_real
!
!> \brief   MPAS 5D Scratch real deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 5D scratch real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field5d_real(f, single_block_in)!{{{
       type (field5dReal), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field5dReal), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field5d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_scratch_field1D_char
!
!> \brief   MPAS 1D Scratch character deallocation rotuine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D scratch character field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_scratch_field1d_char(f, single_block_in)!{{{
       type (field1dChar), pointer :: f !< Input: Field to deallocate
       logical, intent(in), optional :: single_block_in !< Input: Logical that determines if a single block should be deallocated, or all blocks.
       logical :: single_block
       type (field1dChar), pointer :: f_cursor

       if(f % isPersistent) then
          return
       end if

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field1d_char!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field0d_integer
!
!> \brief   MPAS 0D integer deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 0D integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field0d_integer(f)!{{{
       type (field0dInteger), pointer :: f !< Input: Field to deallocate
       type (field0dInteger), pointer :: f_cursor

       f_cursor => f

       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         deallocate(f_cursor)
         f_cursor => f
       end do

   end subroutine mpas_deallocate_field0d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field1D_integer
!
!> \brief   MPAS 1D integer deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field1d_integer(f)!{{{
       type (field1dInteger), pointer :: f !< Input: Field to deallocate
       type (field1dInteger), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field1d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field2D_integer
!
!> \brief   MPAS 2D integer deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 2D integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field2d_integer(f)!{{{
       type (field2dInteger), pointer :: f !< Input: Field to deallocate
       type (field2dInteger), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field2d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field3D_integer
!
!> \brief   MPAS 3D integer deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 3D integer field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field3d_integer(f)!{{{
       type (field3dInteger), pointer :: f !< Input: Field to deallocate
       type (field3dInteger), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field3d_integer!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field0d_real
!
!> \brief   MPAS 0D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 0D real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field0d_real(f)!{{{
       type (field0dReal), pointer :: f !< Input: Field to deallocate
       type (field0dReal), pointer :: f_cursor

       f_cursor => f

       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field0d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field1D_real
!
!> \brief   MPAS 1D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field1d_real(f)!{{{
       type (field1dReal), pointer :: f !< Input: Field to deallocate
       type (field1dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field1d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field2D_real
!
!> \brief   MPAS 2D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 2D real field.
!

   subroutine mpas_deallocate_field2d_real(f)!{{{
       type (field2dReal), pointer :: f !< Input: Field to deallocate
       type (field2dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field2d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field3D_real
!
!> \brief   MPAS 3D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 3D real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field3d_real(f)!{{{
       type (field3dReal), pointer :: f !< Input: Field to deallocate
       type (field3dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field3d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field4D_real
!
!> \brief   MPAS 4D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 4D real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field4d_real(f)!{{{
       type (field4dReal), pointer :: f !< Input: Field to deallocate
       type (field4dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field4d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field5D_real
!
!> \brief   MPAS 5D real deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 5D real field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field5d_real(f)!{{{
       type (field5dReal), pointer :: f !< Input: Field to deallocate
       type (field5dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field5d_real!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field0D_char
!
!> \brief   MPAS 0D character deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 0D character field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field0d_char(f)!{{{
       type (field0dChar), pointer :: f !< Input: Field to deallocate
       type (field0dChar), pointer :: f_cursor

       f_cursor => f

       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         deallocate(f_cursor)
         f_cursor => f
       end do

   end subroutine mpas_deallocate_field0d_char!}}}

!***********************************************************************
!
!  routine mpas_deallocate_field1D_char
!
!> \brief   MPAS 1D character deallocation routine.
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a 1D character field.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_field1d_char(f)!{{{
       type (field1dChar), pointer :: f !< Input: Field to deallocate
       type (field1dChar), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field1d_char!}}}

!***********************************************************************
!
!  routine mpas_deallocate_block
!
!> \brief   MPAS Block deallocation routine
!> \author  Doug Jacobsen
!> \date    04/02/13
!> \details 
!> This routine deallocates a block structure.
!
!-----------------------------------------------------------------------
   subroutine mpas_deallocate_block(b)!{{{
 
      implicit none

      type (block_type), intent(inout) :: b !< Input/Output: Block to be deallocated.

      integer :: i

      ! BUG: It seems like we should be deallocating the exchange lists before we 
      !      deallocate the array of head pointers and the parinfo type...
      !      It also seems like these deallocations should happen with mpas_dmpar_destroy_multihalo_exchange_list

      deallocate(b % parinfo % cellsToSend)
      deallocate(b % parinfo % cellsToRecv)
      deallocate(b % parinfo % cellsToCopy)

      deallocate(b % parinfo % edgesToSend)
      deallocate(b % parinfo % edgesToRecv)
      deallocate(b % parinfo % edgesToCopy)

      deallocate(b % parinfo % verticesToSend)
      deallocate(b % parinfo % verticesToRecv)
      deallocate(b % parinfo % verticesToCopy)

      deallocate(b % parinfo)

#include "block_deallocs.inc"

   end subroutine mpas_deallocate_block!}}}


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
   subroutine mpas_pool_create_pool(newPool)!{{{
      type (mpas_pool_type), pointer :: newPool
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
   subroutine mpas_pool_destroy_pool(inPool)!{{{
      type (mpas_pool_type), pointer :: inPool
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
      type (mpas_pool_type), pointer :: inPool
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
   subroutine mpas_pool_clone_pool(srcPool, destPool)!{{{
      type (mpas_pool_type), pointer :: srcPool
      type (mpas_pool_type), pointer :: destPool
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
   subroutine mpas_pool_copy_pool(srcPool, destPool)!{{{
      type (mpas_pool_type), pointer :: srcPool
      type (mpas_pool_type), pointer :: destPool
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field0DReal), pointer :: field
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field1DReal), pointer :: field
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field2DReal), pointer :: field
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field3DReal), pointer :: field
   end subroutine mpas_pool_add_field_3d_real!}}}

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
   subroutine mpas_pool_add_field_4d_real(inPool, key, field)!{{{
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field4DReal), pointer :: field
   end subroutine mpas_pool_add_field_4d_real!}}}

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
   subroutine mpas_pool_add_field_5d_real(inPool, key, field)!{{{
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field5DReal), pointer :: field
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field0DInteger), pointer :: field
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field1DInteger), pointer :: field
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field2DInteger), pointer :: field
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field3DInteger), pointer :: field
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field0DChar), pointer :: field
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field1DChar), pointer :: field
   end subroutine mpas_pool_add_field_1d_char!}}}

!-----------------------------------------------------------------------
!  routine mpas_pool_add_field_0d_reals
!
!> \brief MPAS Pool 0D Multi-level Real field add routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine inserts fields into inPool when fields is a multi-level  0D real field
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_add_field_0d_reals(inPool, key, fields)!{{{
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field0DReal), dimension(:), pointer :: fields
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field1DReal), dimension(:), pointer :: fields
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field2DReal), dimension(:), pointer :: fields
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field3DReal), dimension(:), pointer :: fields
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field4DReal), dimension(:), pointer :: fields
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field5DReal), dimension(:), pointer :: fields
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field0DInteger), dimension(:), pointer :: fields
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field1DInteger), dimension(:), pointer :: fields
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field2DInteger), dimension(:), pointer :: fields
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field3DInteger), dimension(:), pointer :: fields
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field0DChar), dimension(:), pointer :: fields
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
      type (mpas_pool_type), pointer, intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (field1DChar), dimension(:), pointer :: fields
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (mpas_pool_field_info_type), intent(out) :: info
   end subroutine mpas_pool_get_field_info!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_0d_real
!
!> \brief MPAS Pool 0D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_0d_real(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field0DReal), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_0d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_1d_real
!
!> \brief MPAS Pool 1D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_1d_real(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field1DReal), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_1d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_2d_real
!
!> \brief MPAS Pool 2D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_2d_real(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field2DReal), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_2d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_3d_real
!
!> \brief MPAS Pool 3D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_3d_real(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field3DReal), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_3d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_4d_real
!
!> \brief MPAS Pool 4D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_4d_real(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field4DReal), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_4d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_5d_real
!
!> \brief MPAS Pool 5D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_5d_real(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field5DReal), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_5d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_0d_int
!
!> \brief MPAS Pool 0D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_0d_int(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field0DInteger), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_0d_int!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_1d_int
!
!> \brief MPAS Pool 1D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_1d_int(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field1DInteger), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_1d_int!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_2d_int
!
!> \brief MPAS Pool 2D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_2d_int(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field2DInteger), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_2d_int!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_3d_int
!
!> \brief MPAS Pool 3D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_3d_int(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field3DInteger), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_3d_int!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_0d_char
!
!> \brief MPAS Pool 0D Character field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_0d_char(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field0DChar), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_0d_char!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_field_1d_char
!
!> \brief MPAS Pool 1D Character field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the field associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_field_1d_char(inPool, key, field, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      type (field1DChar), intent(out), pointer :: field
   end subroutine mpas_pool_get_field_1d_char!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_0d_real
!
!> \brief MPAS Pool 0D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_0d_real(inPool, key, scalar, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      real (kind=RKIND), intent(out), pointer :: scalar
   end subroutine mpas_pool_get_array_0d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_1d_real
!
!> \brief MPAS Pool 1D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_1d_real(inPool, key, array, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      real (kind=RKIND), intent(out), dimension(:), pointer :: array
   end subroutine mpas_pool_get_array_1d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_2d_real
!
!> \brief MPAS Pool 2D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_2d_real(inPool, key, array, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      real (kind=RKIND), intent(out), dimension(:,:), pointer :: array
   end subroutine mpas_pool_get_array_2d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_3d_real
!
!> \brief MPAS Pool 3D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_3d_real(inPool, key, array, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      real (kind=RKIND), intent(out), dimension(:,:,:), pointer :: array
   end subroutine mpas_pool_get_array_3d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_4d_real
!
!> \brief MPAS Pool 4D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_4d_real(inPool, key, array, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      real (kind=RKIND), intent(out), dimension(:,:,:,:), pointer :: array
   end subroutine mpas_pool_get_array_4d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_5d_real
!
!> \brief MPAS Pool 5D Real field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_5d_real(inPool, key, array, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      real (kind=RKIND), intent(out), dimension(:,:,:,:,:), pointer :: array
   end subroutine mpas_pool_get_array_5d_real!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_0d_int
!
!> \brief MPAS Pool 0D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_0d_int(inPool, key, scalar, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      integer, intent(out), pointer :: scalar
   end subroutine mpas_pool_get_array_0d_int!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_1d_int
!
!> \brief MPAS Pool 1D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_1d_int(inPool, key, array, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      integer, intent(out), dimension(:), pointer :: array
   end subroutine mpas_pool_get_array_1d_int!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_2d_int
!
!> \brief MPAS Pool 2D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_2d_int(inPool, key, array, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      integer, intent(out), dimension(:,:), pointer :: array
   end subroutine mpas_pool_get_array_2d_int!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_3d_int
!
!> \brief MPAS Pool 3D Integer field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_3d_int(inPool, key, array, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      integer, intent(out), dimension(:,:,:), pointer :: array
   end subroutine mpas_pool_get_array_3d_int!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_0d_char
!
!> \brief MPAS Pool 0D Character field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_0d_char(inPool, key, string, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      character (len=StrKIND), intent(out), pointer :: string
   end subroutine mpas_pool_get_array_0d_char!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_array_1d_char
!
!> \brief MPAS Pool 1D Character field get subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine returns a pointer to the array associted with key in inPool.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_array_1d_char(inPool, key, array, timeLevel)!{{{
      type (mpas_pool_type), intent(in), pointer :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in), optional :: timeLevel
      character (len=StrKIND), intent(out), dimension(:), pointer :: array
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      real (kind=RKIND), intent(in) :: value
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in) :: value
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      character (len=*), intent(in) :: value
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      logical, intent(in) :: value
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      real (kind=RKIND), pointer, intent(out) :: value
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, pointer, intent(out) :: value
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      character (len=StrKIND), pointer, intent(out) :: value
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      logical, pointer, intent(out) :: value
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, intent(in) :: dim
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, dimension(:), intent(in) :: dims
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, intent(out) :: dim
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      integer, intent(out), dimension(:) :: dims
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (mpas_pool_type), intent(in) :: subPool
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      type (mpas_pool_type), intent(out), pointer :: subPool
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      logical, intent(in) :: value
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
      logical, pointer, intent(out) :: package
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
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
      type (mpas_pool_type), intent(in) :: inPool
      character (len=*), intent(in) :: key
   end subroutine mpas_pool_remove_package!}}}

!-----------------------------------------------------------------------
!  routine mpas_pool_begin_iteration
!
!> \brief MPAS Pool Begin Iteration Routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This routine sets up the pools internal iteartor to iterate over fields.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_begin_iteration(inPool)!{{{
      type (mpas_pool_type), intent(in) :: inPool
   end subroutine mpas_pool_begin_iteration!}}}

!-----------------------------------------------------------------------
!  subroutine mpas_pool_get_next_member
!
!> \brief MPAS Pool Iterate To Next Member subroutine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/27/2014
!> \details
!> This subroutine advances the internal iterator to the next member in the pool,
!>  and returns an iterator type for the current member.
!
!-----------------------------------------------------------------------
   subroutine mpas_pool_get_next_member(inPool, iterator)!{{{
      type (mpas_pool_type), intent(in) :: inPool
      type (mpas_pool_iterator_type), intent(out), pointer :: iterator
   end subroutine mpas_pool_get_next_member!}}}

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
   subroutine mpas_pool_shift_time_levels(inPool)!{{{
      type (mpas_pool_type), intent(in) :: inPool
   end subroutine mpas_pool_shift_time_levels!}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           END POOL SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#include "group_dealloc_routines.inc"


#include "group_copy_routines.inc"


#include "group_shift_level_routines.inc"


#include "field_links.inc"


end module mpas_grid_types
