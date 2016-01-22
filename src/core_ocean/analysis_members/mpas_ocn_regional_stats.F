! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ocn_regional_stats
!
!> \brief MPAS ocean analysis core member: regional_stats
!> \author Jon Woodring
!> \date   December 17, 2015
!> \details
!>  Flexible regional averaging, mins, and maxes of fields.
!-----------------------------------------------------------------------
module ocn_regional_stats
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_dmpar
  use mpas_timekeeping
  use mpas_stream_manager

  use ocn_constants
  use ocn_diagnostics_routines

  implicit none
  private
  save

  ! Public parameters
  !--------------------------------------------------------------------

  ! Public member functions
  !--------------------------------------------------------------------
  public :: &
         ocn_packages_regional_stats, &
         ocn_init_regional_stats, &
         ocn_compute_regional_stats, &
         ocn_restart_regional_stats, &
         ocn_finalize_regional_stats

  ! Private module variables
  !--------------------------------------------------------------------

  type regional_variable_type
    ! state per instance variable, stored in framework
    character (len=StrKIND), pointer :: input_name

    ! generated on every instance call
    character (len=StrKIND), dimension(:), allocatable :: output_names
  end type regional_variable_type

  type regional_type
    ! dimensions
    integer, pointer :: number_of_regions
    integer, pointer :: max_regions_per
    integer, pointer :: number_of_groups

    ! state per instance, stored in framework
    integer, pointer :: number_of_variables

    integer, pointer :: operation
    integer, pointer :: region_element
    integer, pointer :: weighting_function 

    ! looked up on every instance call
    character (len=StrKIND), pointer :: weighting_field
    type (mpas_pool_field_info_type) :: weighting_info
    integer, dimension(:), pointer :: counts, num_regions
    integer, dimension(:, :), pointer :: groups
    integer, pointer :: group_index

    ! generated on every instance call
    character (len=StrKIND) :: mask_pool, mask_field

    ! allocated on every instance call
    type (regional_variable_type), dimension(:), allocatable :: variables
  end type regional_type

  ! enum of ops and types
  integer, parameter :: AVG_OP = 1
  integer, parameter :: MIN_OP = 2
  integer, parameter :: MAX_OP = 3

  integer, parameter :: ID_FUNC = 4
  integer, parameter :: MUL_FUNC = 5
  integer, parameter :: DIV_FUNC = 6

  integer, parameter :: CELL_REGION = 7
  integer, parameter :: VERTEX_REGION = 8
  integer, parameter :: EDGE_REGION = 9

  ! namelist operators and identifiers for unique names
  character (len=3), parameter :: AVG_TOKEN = 'avg'
  character (len=3), parameter :: MIN_TOKEN = 'min'
  character (len=3), parameter :: MAX_TOKEN = 'max'

  ! namelist operators
  character (len=2), parameter :: ID_TOKEN  = 'id'
  character (len=3), parameter :: MUL_TOKEN = 'mul'
  character (len=3), parameter :: DIV_TOKEN = 'div'

  ! namelist operators and identifiers for mask-struct data
  character (len=5), parameter :: CELL_TOKEN = 'cell'
  character (len=8), parameter :: VERTEX_TOKEN = 'vertex'
  character (len=5), parameter :: EDGE_TOKEN = 'edge'

  ! canonical dimension names
  character (len=6), parameter :: CELL_DIM = 'nCells'
  character (len=9), parameter :: VERTEX_DIM = 'nVertices'
  character (len=6), parameter :: EDGE_DIM = 'nEdges'

  ! canonical solve names
  character (len=11), parameter :: CELL_SOLVE = 'nCellsSolve'
  character (len=14), parameter :: VERTEX_SOLVE = 'nVerticesSolve'
  character (len=11), parameter :: EDGE_SOLVE = 'nEdgesSolve'

  ! canonical mins and maxes
  real (kind=RKIND), parameter :: DEFAULT_MPAS_MIN_VAL = -1.0e34
  real (kind=RKIND), parameter :: DEFAULT_MPAS_MAX_VAL =  1.0e34

  ! xtime
  character (len=5), parameter :: TIME_STREAM = 'xtime'

  ! none
  character (len=4), parameter :: NONE_TOKEN = 'none'

  ! memory names for duplication 
  character (len=StrKIND), parameter :: REGIONAL_STATS_POOL = &
    'regionalStatsAM'
  character (len=StrKIND), parameter :: ONE_INTEGER_MEMORY = &
    'regionalStatsOneInteger'
  character (len=StrKIND), parameter :: ONE_STRING_MEMORY = &
    'regionalStatsOneString'
 
  ! mask-struct array suffixes
  character (len=StrKIND), parameter :: MASK_POOL_SUFFIX = 'Regions'
  character (len=StrKIND), parameter :: MASKS_DATA_SUFFIX = 'RegionMasks'
  character (len=StrKIND), parameter :: COUNTS_DATA_SUFFIX = 'NElementsInRegion'
  character (len=StrKIND), parameter :: GROUPS_DATA_SUFFIX = 'RegionsInGroup'
  character (len=StrKIND), parameter :: REGCOUNTS_DATA_SUFFIX = &
    'NRegionsInGroup'
  character (len=6), parameter :: REGION_PREFIX = 'region'

  ! mask dimension names
  character (len=StrKIND), parameter :: NUM_REGIONS_SUFFIX = 'nRegions'
  character (len=StrKIND), parameter :: NUM_GROUPS_SUFFIX = 'nRegionGroups'
  character (len=StrKIND), parameter :: MAX_REGIONS_SUFFIX = &
    'maxRegionsInGroup'

  ! prefixes
  character (len=StrKIND), parameter :: CONFIG_PREFIX = &
    'config_AM_regionalStats'
  character (len=StrKIND), parameter :: FRAMEWORK_PREFIX = 'regionalStats'

  ! namelist-only suffixes
  character (len=StrKIND), parameter :: OUTPUT_STREAM_SUFFIX = '_output_stream'
  character (len=StrKIND), parameter :: WEIGHTING_FIELD_SUFFIX = &
    '_weighting_field'
  character (len=StrKIND), parameter :: REGION_GROUP_SUFFIX = '_region_group'
  character (len=StrKIND), parameter :: INPUT_STREAM_SUFFIX = '_input_stream'
  character (len=StrKIND), parameter :: RESTART_STREAM_SUFFIX = &
    '_restart_stream'

  ! namelist and instance suffixes
  character (len=StrKIND), parameter :: OPERATION_SUFFIX = '_operation'
  character (len=StrKIND), parameter :: REGION_TYPE_SUFFIX = '_region_type'
  character (len=StrKIND), parameter :: WEIGHTING_FUNCTION_SUFFIX = &
    '_weighting_function'

  ! instance-only suffixes
  character (len=StrKIND), parameter :: INPUT_NAME_SUFFIX = '_input_name'
  character (len=StrKIND), parameter :: NUMBER_OF_VARIABLES_SUFFIX = &
    '_number_of_variables'

  ! package suffixes
  character (len=StrKIND), parameter :: PACKAGE_ACTIVE_SUFFIX = &
    'RegionsPKGActive'

 
!***********************************************************************
contains

!***********************************************************************
! routine ocn_packages_regional_stats
!
!> \brief Enable additional packages for an MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts all additional package enables required for the
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_packages_regional_stats(configPool, packagePool)!{{{
  ! input variables

  ! input/output variables
  type (mpas_pool_type), intent(inout) :: configPool, packagePool

  ! output variables

  ! local variables
  character (len=StrKIND) :: instance ! TODO intent(in)
  logical, pointer :: active
  character (len=StrKIND) :: field_name, namelist_prefix
  character (len=StrKIND), pointer :: config_results

  ! start procedure

  ! TODO placeholder for some unique ID if this code is replicated
  instance = '' ! TODO to be passed in

  namelist_prefix = trim(CONFIG_PREFIX) // trim(instance)
  field_name = trim(namelist_prefix) // trim(REGION_TYPE_SUFFIX)
  call mpas_pool_get_config(configPool, field_name, config_results)

  ! error if unknown element type
  if ((config_results /= CELL_TOKEN) .and. &
      (config_results /= VERTEX_TOKEN) .and. &
      (config_results /= EDGE_TOKEN)) then
    call mpas_dmpar_global_abort('ERROR: unknown region element type in ' // &
      'regional stats analysis member configuration.')
  end if

  ! activate the package
  field_name = trim(config_results) // trim(PACKAGE_ACTIVE_SUFFIX)
  call mpas_pool_get_package(packagePool, field_name, active)
  if (.not. associated(active)) then
    call mpas_dmpar_global_abort('ERROR: unable to activate region ' // &
      'package ' // trim(field_name) // ' which is required for ' // &
      'the regional stats AM.')
  end if
  active = .true.

end subroutine ocn_packages_regional_stats!}}}



!***********************************************************************
! routine ocn_init_regional_stats
!
!> \brief Initialize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_init_regional_stats(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: v
  character (len=StrKIND) :: instance ! TODO intent(in)

  ! start procedure
  err = 0

  ! TODO placeholder for some unique ID if this code is replicated
  instance = '' ! TODO to be passed in

  ! create all of the state for this instance
  call start_state(domain, instance, err)

end subroutine ocn_init_regional_stats!}}}



!***********************************************************************
! routine ocn_compute_regional_stats
!
!> \brief Compute MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_compute_regional_stats(domain, timeLevel, err)!{{{
  ! input variables
  integer, intent(in) :: timeLevel

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  character (len=StrKIND) :: instance ! TODO intent(in)
  integer :: v
  type (regional_type) :: regions

  ! start procedure
  err = 0

  ! TODO placeholder for some unique ID if this code is replicated
  instance = '' ! TODO to be passed in

  ! get all of the state for this instance to be able to compute
  call get_state(domain, instance, regions)

  ! do all region reductions for each variable
  do v = 1, regions % number_of_variables
    call typed_operate(domain % dminfo, domain % blocklist, &
      regions, regions % variables(v))
  end do

  ! clean up the instance memory
  do v = 1, regions % number_of_variables
    deallocate(regions % variables(v) % output_names)
  end do
  deallocate(regions % variables)
end subroutine ocn_compute_regional_stats!}}}



!***********************************************************************
! routine ocn_restart_regional_stats
!
!> \brief Save restart for MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_restart_regional_stats(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

end subroutine ocn_restart_regional_stats!}}}



!***********************************************************************
! routine ocn_finalize_regional_stats
!
!> \brief Finalize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_finalize_regional_stats(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

end subroutine ocn_finalize_regional_stats!}}}

!
! local subroutines
!

!***********************************************************************
! routine add_new_string
!
!> \brief Allocate an string in the MPAS framework for this AM
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Allocate a new integer in the AM pool and return a pointer to it.
!-----------------------------------------------------------------------
subroutine add_new_string(pool, field_name, target_ptr)
  ! input variables
  character (len=StrKIND) :: field_name

  ! input/output variables
  type (mpas_pool_type), pointer, intent(inout) :: pool

  ! output variables
  character (len=StrKIND), pointer :: target_ptr

  ! local variables
  type (field0DChar), pointer :: srcString, dstString

  call mpas_pool_get_field(pool, ONE_STRING_MEMORY, srcString, 1)
  call mpas_duplicate_field(srcString, dstString)
  dstString % fieldName = field_name
  call mpas_pool_add_field(pool, dstString % fieldName, dstString)
  call mpas_pool_get_array(pool, dstString % fieldName, target_ptr, 1)
end subroutine add_new_string



!***********************************************************************
! routine add_new_integer
!
!> \brief Allocate an integer in the MPAS framework for this AM
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Allocate a new integer in the AM pool and return a pointer to it.
!-----------------------------------------------------------------------
subroutine add_new_integer(pool, field_name, target_ptr)
  ! input variables
  character (len=StrKIND) :: field_name

  ! input/output variables
  type (mpas_pool_type), pointer, intent(inout) :: pool

  ! output variables
  integer, pointer :: target_ptr

  ! local variables
  type (field0DInteger), pointer :: srcInteger, dstInteger

  call mpas_pool_get_field(pool, ONE_INTEGER_MEMORY, srcInteger, 1)
  call mpas_duplicate_field(srcInteger, dstInteger)
  dstInteger % fieldName = field_name
  call mpas_pool_add_field(pool, dstInteger % fieldName, dstInteger)
  call mpas_pool_get_array(pool, dstInteger % fieldName, target_ptr, 1)
end subroutine add_new_integer



!***********************************************************************
! function output_naming
!
!> \brief Given an input name, create a corresponding output name
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Code to create consistent output names from input names.
!-----------------------------------------------------------------------
character (len=StrKIND) function output_naming &
(op_name, input_name, region_identifier, instance)
  character (len=StrKIND), intent(in) :: op_name, &
    input_name, region_identifier, instance

  output_naming = trim(instance) // &
    trim(region_identifier) // '_' // trim(op_name) // '_' // &
    trim(input_name) 
end function output_naming



!***********************************************************************
! function operator_naming
!
!> \brief Given an operator enum, create a corresponding operator name
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Code to create consistent operator names from operator enums.
!-----------------------------------------------------------------------
character (len=StrKIND) function operator_naming(enum)
  integer, intent(in) :: enum

  if (enum == AVG_OP) then
    operator_naming = AVG_TOKEN 
  else if (enum == MIN_OP) then
    operator_naming = MIN_TOKEN 
  else if (enum == MAX_OP) then
    operator_naming = MAX_TOKEN 
  else
    call mpas_dmpar_global_abort('ERROR: The impossible happened. ' // &
      'Tried creating an operator of unknown kind in regional stats AM.')
  end if 
end function operator_naming



!***********************************************************************
! function element_naming
!
!> \brief Given an element enum, create a corresponding element name
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Code to create consistent operator names from operator enums.
!-----------------------------------------------------------------------
character (len=StrKIND) function element_naming(enum)
  integer, intent(in) :: enum

  if (enum == CELL_REGION) then
    element_naming = CELL_TOKEN 
  else if (enum == VERTEX_REGION) then
    element_naming = VERTEX_TOKEN
  else if (enum == EDGE_REGION) then
    element_naming = EDGE_TOKEN
  else
    call mpas_dmpar_global_abort('ERROR: The impossible happened. ' // &
      'Tried creating an element of unknown kind in regional stats AM.')
  end if 
end function element_naming



!***********************************************************************
! function dimension_naming
!
!> \brief Given a region type, return its canonical dimension string name
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Returns the canonical MPAS dimension name given a region element dimension.
!-----------------------------------------------------------------------
character (len=StrKIND) function dimension_naming(elem_type)
  integer, intent(in) :: elem_type

  if (elem_type == CELL_REGION) then
    dimension_naming = CELL_DIM
  else if (elem_type == VERTEX_REGION) then
    dimension_naming = VERTEX_DIM
  else if (elem_type == EDGE_REGION) then
    dimension_naming = EDGE_DIM
  else
    call mpas_dmpar_global_abort('ERROR: The impossible happened. ' // &
      'Tried creating a element dimension name of unknown kind in ' // &
      'regional stats AM.')
  end if
end function dimension_naming



!***********************************************************************
! function check_element_dim
!
!> \brief Return true if the last dimension matches the element name.
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Given an element dimension name and list of dim names, return true 
!> if the last dimension is the element, otherwise return false.
!-----------------------------------------------------------------------
integer function check_element_dim(dim_names, elem_name)
  character (len=StrKIND), dimension(:), intent(in) :: dim_names
  character (len=StrKIND), intent(in) :: elem_name

  integer :: last

  last = size(dim_names)
  check_element_dim = (trim(dim_names(last)) == trim(elem_name)) 

end function check_element_dim



!***********************************************************************
! routine debug_state
!
!> \brief Print all of the state for this instance.
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Given the current state in regions, print it all for debugging.
!-----------------------------------------------------------------------
subroutine debug_state(regions)
  type (regional_type), intent(in) :: regions
  integer :: v, b, last, m

  write(*,*) 'num vars: ', regions % number_of_variables
  write(*,*) 'num regions: ', regions % number_of_regions
  write(*,*) 'num groups: ', regions % number_of_groups
  write(*,*) 'max regions per: ', regions % max_regions_per
 
  write (*,*) 'op: ', regions % operation
  write (*,*) 'elem: ', regions % region_element
  write (*,*) 'func: ', regions % weighting_function
  write (*,*) 'weight: ', trim(regions % weighting_field)

  if (regions % weighting_field /= NONE_TOKEN) then
    write (*,*) 'type: ', regions % weighting_info % fieldType
  end if
  write (*,*) 'pool: ', trim(regions % mask_pool)
  write (*,*) 'mask: ', trim(regions % mask_field)

  write(*,*) 'regions in groups:'
  do v = 1, regions % number_of_groups
    write(*,*) 'group: ', v
    m = regions % num_regions(v)
    write(*,*) 'num regions: ', m
    do b = 1, regions % max_regions_per
      write(*,*) 'index: ', b
      write(*,*) 'region index: ', regions % groups(b, v)
    end do
  end do

  write(*,*) 'elements in region:'
  do v = 1, regions % number_of_regions
    write(*,*) 'region: ', v
    write(*,*) 'num elem: ', regions % counts(v)
  end do

  write(*,*) 'selected group: ', regions % group_index
  last = regions % num_regions(regions % group_index)
  write(*,*) 'num regions: ', last
  do v = 1, regions % number_of_variables
    write(*,*) 'var: ', trim(regions % variables(v) % input_name)
    do b = 1, last
      write(*,*) 'index: ', b
      m = regions % groups(b, regions % group_index)
      write(*,*) 'region index: ', m
      write(*,*) 'active elems: ', regions % counts(m)
      write(*,*) 'out var: ', trim(regions % variables(v) % output_names(b))
    end do
  end do
end subroutine debug_state 



!***********************************************************************
! routine get_state
!
!> \brief Get all of the state for this instance.
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> This will allocate and fetch all of the state necessary for this
!> instance that is being run (prevalidated data from start_state).
!-----------------------------------------------------------------------
subroutine get_state(domain, instance, regions)
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  type (regional_type), intent(out) :: regions

  ! local variables
  integer :: v, b, last
  character (len=StrKIND) :: namelist_prefix, storage_prefix, field_name, &
    element_prefix, var_identifier, op_name
  type (mpas_pool_type), pointer :: maskPool, amPool

  ! start procedure
  namelist_prefix = trim(CONFIG_PREFIX) // trim(instance)
  storage_prefix = trim(FRAMEWORK_PREFIX) // trim(instance)

  call mpas_pool_get_subpool(domain % blocklist % structs, &
    REGIONAL_STATS_POOL, amPool)

  ! 
  ! get dimensions
  !

  call mpas_pool_get_dimension(domain % blocklist % dimensions, &
    NUM_REGIONS_SUFFIX, regions % number_of_regions)
  call mpas_pool_get_dimension(domain % blocklist % dimensions, &
    NUM_GROUPS_SUFFIX, regions % number_of_groups)
  call mpas_pool_get_dimension(domain % blocklist % dimensions, &
    MAX_REGIONS_SUFFIX, regions % max_regions_per)

  !
  ! get ones we stored in framework for this instance
  !

  ! number_of_variables
  field_name = trim(storage_prefix) // trim(NUMBER_OF_VARIABLES_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, regions % number_of_variables, 1)

  ! operation
  field_name = trim(storage_prefix) // trim(OPERATION_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, regions % operation, 1)

  op_name = operator_naming(regions % operation)

  ! mask type
  field_name = trim(storage_prefix) // trim(REGION_TYPE_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, regions % region_element, 1)

  element_prefix = element_naming(regions % region_element)

  ! weighting function
  field_name = trim(storage_prefix) // trim(WEIGHTING_FUNCTION_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, regions % weighting_function, 1)

  ! get the input names for variables
  allocate(regions % variables(regions % number_of_variables))
  do v = 1, regions % number_of_variables
    ! identifier
    write(var_identifier, '(I0)') v
    field_name = trim(storage_prefix) // '_' // trim(var_identifier) // &
      trim(INPUT_NAME_SUFFIX)

    call mpas_pool_get_array(amPool, field_name, &
      regions % variables(v) % input_name, 1)
  end do

  ! 
  ! get ones that already exist from namelist/stream
  !

  ! weighting field & info
  field_name = trim(namelist_prefix) // trim(WEIGHTING_FIELD_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % weighting_field)

  if (regions % weighting_function /= ID_FUNC) then
    call mpas_pool_get_field_info(domain % blocklist % allFields, &
      regions % weighting_field, regions % weighting_info)
  end if

  ! mask field & pool
  regions % mask_pool = trim(element_prefix) // trim(MASK_POOL_SUFFIX)
  regions % mask_field = trim(element_prefix) // trim(MASKS_DATA_SUFFIX)
  call mpas_pool_get_subpool(domain % blocklist % structs, &
    regions % mask_pool, maskPool)

  ! counts
  field_name = trim(element_prefix) // trim(COUNTS_DATA_SUFFIX)
  call mpas_pool_get_array(maskPool, field_name, regions % counts, 1)

  ! groups
  field_name = trim(element_prefix) // trim(GROUPS_DATA_SUFFIX)
  call mpas_pool_get_array(maskPool, field_name, regions % groups, 1)

  ! num_regions
  field_name = trim(element_prefix) // trim(REGCOUNTS_DATA_SUFFIX)
  call mpas_pool_get_array(maskPool, field_name, regions % num_regions, 1)

  ! group index
  field_name = trim(namelist_prefix) // trim(REGION_GROUP_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % group_index)

  ! 
  ! generate output names
  !

  do v = 1, regions % number_of_variables
    ! allocate output names
    last = regions % num_regions(regions % group_index)
    allocate(regions % variables(v) % output_names(last))

    do b = 1, last
      ! region name
      write(field_name, '(I0)') regions % groups(b, regions % group_index)
      field_name = REGION_PREFIX // trim(field_name)

      regions % variables(v) % output_names(b) = output_naming(op_name, &
        regions % variables(v) % input_name, field_name, instance)
    end do
  end do

end subroutine get_state


!***********************************************************************
! routine start_state
!
!> \brief Begin the initialization of this analysis member
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> All the necessary details to initialize this analysis member
!> instance.
!-----------------------------------------------------------------------
subroutine start_state(domain, instance, err)
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  type (regional_type) :: regions
  character (len=StrKIND), pointer :: config_results, output_stream_name, &
    config_results_2
  character (len=StrKIND) :: namelist_prefix, storage_prefix, field_name, &
    element_prefix, elem_name, op_name, var_identifier
  integer :: v, b, last
  type (mpas_pool_type), pointer :: maskPool, amPool
  type (field1DReal), pointer :: real_field
  type (field1DInteger), pointer :: integer_field
  type (mpas_pool_field_info_type) :: info

  ! start procedure
  err = 0

  namelist_prefix = trim(CONFIG_PREFIX) // trim(instance)
  storage_prefix = trim(FRAMEWORK_PREFIX) // trim(instance)

  call mpas_pool_get_subpool(domain % blocklist % structs, &
    REGIONAL_STATS_POOL, amPool)

  ! 
  ! get dimensions
  !

  call mpas_pool_get_dimension(domain % blocklist % dimensions, &
    NUM_REGIONS_SUFFIX, regions % number_of_regions)
  call mpas_pool_get_dimension(domain % blocklist % dimensions, &
    NUM_GROUPS_SUFFIX, regions % number_of_groups)
  call mpas_pool_get_dimension(domain % blocklist % dimensions, &
    MAX_REGIONS_SUFFIX, regions % max_regions_per)

  !
  ! allocate some framework memory for instance state, and assign values
  !

  ! number of variables done in modify stream

  ! operation
  field_name = trim(storage_prefix) // trim(OPERATION_SUFFIX)
  call add_new_integer(amPool, field_name, regions % operation)

  ! get our operation
  field_name = trim(namelist_prefix) // trim(OPERATION_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, config_results)
  if (config_results == AVG_TOKEN) then
    regions % operation = AVG_OP
  else if (config_results == MIN_TOKEN) then
    regions % operation = MIN_OP
  else if (config_results == MAX_TOKEN) then
    regions % operation = MAX_OP
  else
    ! error if unknown operation
    call mpas_dmpar_global_abort('ERROR: unknown operation in ' // &
      'regional stats analysis member configuration.')
  end if

  op_name = operator_naming(regions % operation)

  ! element type
  field_name = trim(storage_prefix) // trim(REGION_TYPE_SUFFIX)
  call add_new_integer(amPool, field_name, regions % region_element)

  ! get the element type
  field_name = trim(namelist_prefix) // trim(REGION_TYPE_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, config_results)

  if (config_results == CELL_TOKEN) then
    regions % region_element = CELL_REGION
  else if (config_results == VERTEX_TOKEN) then
    regions % region_element = VERTEX_REGION
  else if (config_results == EDGE_TOKEN) then
    regions % region_element = EDGE_REGION
  else
    ! error if unknown element type
    call mpas_dmpar_global_abort('ERROR: unknown region element type in ' // &
      'regional stats analysis member configuration.')
  end if

  element_prefix = element_naming(regions % region_element)
  elem_name = dimension_naming(regions % region_element)

  ! weighting function
  field_name = trim(storage_prefix) // trim(WEIGHTING_FUNCTION_SUFFIX)
  call add_new_integer(amPool, field_name, regions % weighting_function)

  field_name = trim(namelist_prefix) // trim(WEIGHTING_FUNCTION_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, config_results)
  if (config_results == ID_TOKEN) then
    regions % weighting_function = ID_FUNC
  else if (config_results == MUL_TOKEN) then
    regions % weighting_function = MUL_FUNC
  else if (config_results == DIV_TOKEN) then
    regions % weighting_function = DIV_FUNC
  else
    ! error if unknown operation
    call mpas_dmpar_global_abort('ERROR: unknown weighting function in ' // &
      'regional stats analysis member configuration.')
  end if


  !
  ! verify the other data is OK before trying to use it later
  !

  ! validate input stream is not none
  field_name = trim(namelist_prefix) // trim(INPUT_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, config_results)

  if (config_results == NONE_TOKEN) then
    call mpas_dmpar_global_abort('ERROR: input stream for regional ' // &
      'stats AM cannot be "none". It needs to point to the region/mask stream.')
  end if

  ! validate the restart stream is not none
  field_name = trim(namelist_prefix) // trim(RESTART_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, config_results_2)

  if (config_results == NONE_TOKEN) then
    call mpas_dmpar_global_abort('ERROR: restart stream for regional ' // &
      'stats AM cannot be "none". It needs to point to the region/mask stream.')
  end if

  ! make sure input and restart are the same
  if (config_results /= config_results_2) then
    call mpas_dmpar_global_abort('ERROR: the input and restart stream ' // &
      'for regional stats AM need to be the same to ensure that it ' // &
      'has the same behavior on start and restart.')
  end if

  ! validate weighting field
  field_name = trim(namelist_prefix) // trim(WEIGHTING_FIELD_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % weighting_field)

  if (regions % weighting_function == ID_FUNC) then
    if (regions % weighting_field /= NONE_TOKEN) then
      call mpas_dmpar_global_abort('ERROR: weighting field "' // &
        trim(regions % weighting_field) // '" is not "none"' // &
        'when the weighting function is set to "id" in ' // &
        'regional stats AM.')
    end if
  else
    ! weighting field & info
    if (regions % weighting_field == NONE_TOKEN) then
      call mpas_dmpar_global_abort('ERROR: weighting field "' // &
        trim(regions % weighting_field) // '" is set to "none"' // &
        'when the weighting function is not set to "id" in ' // &
        'regional stats AM.')
    else
      call mpas_pool_get_field_info(domain % blocklist % allFields, &
        regions % weighting_field, regions % weighting_info)
  
      ! check if we can handle it
      if (regions % weighting_info % fieldType == MPAS_POOL_REAL) then
        call mpas_pool_get_field(domain % blocklist % allFields, &
          regions % weighting_field, real_field)

        if (.not. check_element_dim(real_field % dimNames, elem_name)) then 
          call mpas_dmpar_global_abort('ERROR: weighting field ' // &
            regions % weighting_field // &
            ' in regional stats AM needs to only have ' // &
            elem_name // ' dimension (rank-1).')
        end if
      else if (regions % weighting_info % fieldType == MPAS_POOL_INTEGER) then
        call mpas_pool_get_field(domain % blocklist % allFields, &
          regions % weighting_field, integer_field)

        if (.not. check_element_dim(integer_field % dimNames, elem_name)) then
          call mpas_dmpar_global_abort('ERROR: weighting field ' // &
            regions % weighting_field // &
            ' in regional stats AM needs to only have ' // &
            elem_name // ' dimension (rank-1).')
        end if
      else
        call mpas_dmpar_global_abort('ERROR: weighting field "' // &
          trim(regions % weighting_field) // '" listed in the ' // &
          'namelist, for regional stats analysis member ' // &
          'stream, is not real or integer.')
      end if
    end if
  end if

  ! mask field & pool
  regions % mask_pool = trim(element_prefix) // trim(MASK_POOL_SUFFIX)
  regions % mask_field = trim(element_prefix) // trim(MASKS_DATA_SUFFIX)
  call mpas_pool_get_subpool(domain % blocklist % structs, &
    trim(element_prefix) // trim(MASK_POOL_SUFFIX), maskPool)

  ! counts
  field_name = trim(element_prefix) // trim(COUNTS_DATA_SUFFIX)
  call mpas_pool_get_array(maskPool, field_name, regions % counts, 1)

  ! groups
  field_name = trim(element_prefix) // trim(GROUPS_DATA_SUFFIX)
  call mpas_pool_get_array(maskPool, field_name, regions % groups, 1)

  ! num_regions
  field_name = trim(element_prefix) // trim(REGCOUNTS_DATA_SUFFIX)
  call mpas_pool_get_array(maskPool, field_name, regions % num_regions, 1)

  ! group index
  field_name = trim(namelist_prefix) // trim(REGION_GROUP_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % group_index)

  ! verify selection is OK

  ! error if out of bounds
  v = regions % group_index
  if ((v < 1) .or. (v > regions % number_of_groups)) then
    call mpas_dmpar_global_abort('ERROR: the group id for ' // & 
      'regional stats analysis member is out of bounds for the size ' // &
      'of the region groups array found in masks input.')
  end if

  last = regions % num_regions(regions % group_index)
  if ((last < 1) .or. (last > regions % max_regions_per)) then
     call mpas_dmpar_global_abort('ERROR: the number of regions ' // &
       'for the selected group is out of bounds for the size of ' // &
       'the actual number of regions found in the masks input.')
  end if

  do b = 1, last
    v = regions % groups(b, regions % group_index)
    if ((v < 1) .or. (v > regions % number_of_regions)) then
      call mpas_dmpar_global_abort('ERROR: a region index found ' // &
        'in the list of regions for the selected group in the ' // &
        'regional stats analysis member is out of bounds for the size ' // &
        'of the regions found in the masks input.')
    end if
  end do

  !
  ! modify the stream to create destination output variables and remove inputs
  !

  ! get the output stream name
  field_name = trim(namelist_prefix) // trim(OUTPUT_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, output_stream_name)

  if (output_stream_name == NONE_TOKEN) then
    call mpas_dmpar_global_abort('ERROR: ouput stream cannot be "none" ' // &
      'for regional stats.')
  end if

  ! number_of_variables
  field_name = trim(storage_prefix) // trim(NUMBER_OF_VARIABLES_SUFFIX) 
  call add_new_integer(amPool, field_name, regions % number_of_variables)

  ! count the number of variables
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  regions % number_of_variables = 0
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    output_stream_name, field_name))
    regions % number_of_variables = regions % number_of_variables + 1
  end do

  ! create the memory
  allocate(regions % variables(regions % number_of_variables))

  ! create input variable name space
  do v = 1, regions % number_of_variables
    ! identifier
    write(var_identifier, '(I0)') v

    ! create input name space
    field_name = trim(storage_prefix) // '_' // trim(var_identifier) // &
      trim(INPUT_NAME_SUFFIX) 
    call add_new_string(amPool, field_name, regions % variables(v) % input_name)
  end do

  ! get the old field names and assign to input name
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  v = 1
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    output_stream_name, field_name))
    regions % variables(v) % input_name = field_name
    v = v + 1
  end do

  ! remove the old ones from the stream
  do v = 1, regions % number_of_variables
    call mpas_stream_mgr_remove_field(domain % streamManager, &
      output_stream_name, regions % variables(v) % input_name)
  end do

  !
  ! modify the stream
  !

  ! add xtime to the output stream
  call mpas_stream_mgr_add_field(domain % streamManager, &
    output_stream_name, TIME_STREAM, ierr=err)

  ! set up the variables and create output memory in stream
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)

  do v = 1, regions % number_of_variables
    ! get the info of the field
    call mpas_pool_get_field_info(domain % blocklist % allFields, &
      regions % variables(v) % input_name, info)

    ! check if we can handle it
    if(.not. ((info % fieldType == MPAS_POOL_REAL) &
       .or. (info % fieldType == MPAS_POOL_INTEGER))) then
      call mpas_dmpar_global_abort('ERROR: field "' // &
        trim(regions % variables(v) % input_name) // '" listed in the ' // &
        'output stream, for regional stats analysis member ' // &
        'stream, is not real or integer.')
    end if

    ! allocate output names
    last = regions % num_regions(regions % group_index)
    allocate(regions % variables(v) % output_names(last))

    ! allocate output fields if the mask is active for this group
    do b = 1, last
      ! region name
      write(field_name, '(I0)') regions % groups(b, regions % group_index)
      field_name = REGION_PREFIX // trim(field_name)

      regions % variables(v) % output_names(b) = output_naming(op_name, &
        regions % variables(v) % input_name, field_name, instance)

      ! create the field and add to pool
      call add_new_field(info, domain % blocklist, &
        regions % variables(v) % input_name, &
        regions % variables(v) % output_names(b), elem_name)

      ! add the field to the output stream
      call mpas_stream_mgr_add_field(domain % streamManager, &
        output_stream_name, regions % variables(v) % output_names(b), &
        ierr=err)
    end do
  end do ! number_of_variables

  ! clean up the instance memory
  do v = 1, regions % number_of_variables
    deallocate(regions % variables(v) % output_names)
  end do
  deallocate(regions % variables)
end subroutine start_state



!***********************************************************************
! routine add_new_field
!
!> \brief Function to create a new field from an existing field
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts all initializations required for
!>  duplicating a field and adding it to the allFields pool.
!-----------------------------------------------------------------------
subroutine add_new_field(info, block, inname, outname, elem_name)!{{{
  ! input variables
  type (mpas_pool_field_info_type), intent(in) :: info
  type (block_type), pointer, intent(inout) :: block
  character (len=StrKIND), intent(in) :: inname, outname, elem_name

  ! input/output variables

  ! output variables

  ! local variables

  ! duplicate field and add new field to pool
  if (info % fieldType == MPAS_POOL_REAL) then
    if (info % nDims == 0) then
      call mpas_dmpar_global_abort('ERROR: input field ' // inname // &
        'to regional stats AM is 0 dimensional. It needs to have at least' &
        // elem_name // ' dimensions.')
    else if (info % nDims == 1) then
      call copy_field_1r(block, inname, outname, elem_name)
    else if (info % nDims == 2) then
      call copy_field_2r(block, inname, outname, elem_name)
    else if (info % nDims == 3) then
      call copy_field_3r(block, inname, outname, elem_name)
    else if (info % nDims == 4) then
      call copy_field_4r(block, inname, outname, elem_name)
    else
      call copy_field_5r(block, inname, outname, elem_name)
    end if
  else
    if (info % nDims == 0) then
      call mpas_dmpar_global_abort('ERROR: input field ' // inname // &
        'to regional stats AM is 0 dimensional. It needs to have at least ' &
        // elem_name // ' dimensions.')
    else if (info % nDims == 1) then
      call copy_field_1i(block, inname, outname, elem_name)
    else if (info % nDims == 2) then
      call copy_field_2i(block, inname, outname, elem_name)
    else
      call copy_field_3i(block, inname, outname, elem_name)
    end if
  end if

end subroutine add_new_field!}}}



!***********************************************************************
! routine typed_operate
!
!> \brief Do the operation, but switch on run-time type
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  Since we don't know the type of the array, we need to do some
!>  run-time type switching based on the type of the array.
!-----------------------------------------------------------------------
subroutine typed_operate(dminfo, blocklist, regions, variable) !{{{

  ! input variables
  type (dm_info), pointer, intent(in) :: dminfo
  type (block_type), pointer, intent(in) :: blocklist
  type (regional_type), intent(in) :: regions
  type (regional_variable_type), intent(in) :: variable

  ! input/output variables

  ! output variables

  ! local variables
  type (mpas_pool_field_info_type) :: info

  ! get the info
  call mpas_pool_get_field_info(blocklist % allFields, &
    variable % input_name, info)

  ! switch based on the type, dimensionality, and operation
  if (info % fieldType == MPAS_POOL_REAL) then
    if (info % nDims == 1) then
      if (regions % operation == AVG_OP) then
        call operate1r_avg(dminfo, blocklist, regions, variable)
      else if (regions % operation == MIN_OP) then
        call operate1r_min(dminfo, blocklist, regions, variable)
      else
        call operate1r_max(dminfo, blocklist, regions, variable)
      end if
    else if (info % nDims == 2) then
      if (regions % operation == AVG_OP) then
        call operate2r_avg(dminfo, blocklist, regions, variable)
      else if (regions % operation == MIN_OP) then
        call operate2r_min(dminfo, blocklist, regions, variable)
      else
        call operate2r_max(dminfo, blocklist, regions, variable)
      end if
    else if (info % nDims == 3) then
      if (regions % operation == AVG_OP) then
        call operate3r_avg(dminfo, blocklist, regions, variable)
      else if (regions % operation == MIN_OP) then
        call operate3r_min(dminfo, blocklist, regions, variable)
      else
        call operate3r_max(dminfo, blocklist, regions, variable)
      end if
    else if (info % nDims == 4) then
      if (regions % operation == AVG_OP) then
        call operate4r_avg(dminfo, blocklist, regions, variable)
      else if (regions % operation == MIN_OP) then
        call operate4r_min(dminfo, blocklist, regions, variable)
      else
        call operate4r_max(dminfo, blocklist, regions, variable)
      end if
    else
      if (regions % operation == AVG_OP) then
        call operate5r_avg(dminfo, blocklist, regions, variable)
      else if (regions % operation == MIN_OP) then
        call operate5r_min(dminfo, blocklist, regions, variable)
      else
        call operate5r_max(dminfo, blocklist, regions, variable)
      end if
    end if
  else
    if (info % nDims == 1) then
      if (regions % operation == AVG_OP) then
        call operate1i_avg(dminfo, blocklist, regions, variable)
      else if (regions % operation == MIN_OP) then
        call operate1i_min(dminfo, blocklist, regions, variable)
      else
        call operate1i_max(dminfo, blocklist, regions, variable)
      end if
    else if (info % nDims == 2) then
      if (regions % operation == AVG_OP) then
        call operate2i_avg(dminfo, blocklist, regions, variable)
      else if (regions % operation == MIN_OP) then
        call operate2i_min(dminfo, blocklist, regions, variable)
      else
        call operate2i_max(dminfo, blocklist, regions, variable)
      end if
    else
      if (regions % operation == AVG_OP) then
        call operate3i_avg(dminfo, blocklist, regions, variable)
      else if (regions % operation == MIN_OP) then
        call operate3i_min(dminfo, blocklist, regions, variable)
      else
        call operate3i_max(dminfo, blocklist, regions, variable)
      end if
    end if
  end if
end subroutine typed_operate!}}}



!***********************************************************************
! routine copy_field_X
!
!> \brief Functions to create a new N-1 D field from an existing ND field
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts initializations required for
!>  duplicating a ND field and adding a N-1D to the AM's pool.
!>  This will only create an output in the first block of a blocklist.
!-----------------------------------------------------------------------

subroutine copy_field_1r(block, inname, outname, elem_name)!{{{
  type (block_type), pointer, intent(inout) :: block
  character (len=StrKIND), intent(in) :: inname, outname, elem_name

  type (field1DReal), pointer :: src
  type (field0DReal), pointer :: dst
  type (mpas_pool_type), pointer :: pool

#include "regional_stats_inc/regional_field_1d_to_0d.F"
end subroutine copy_field_1r!}}}


subroutine copy_field_2r(block, inname, outname, elem_name)!{{{
  type (block_type), pointer, intent(inout) :: block
  character (len=StrKIND), intent(in) :: inname, outname, elem_name

  type (field2DReal), pointer :: src
  type (field1DReal), pointer :: dst
  type (mpas_pool_type), pointer :: pool
  integer :: i
  integer, dimension(2) :: src_dims

#include "regional_stats_inc/regional_field_nd_to_n-1d-1.F"
  allocate(dst % array(src_dims(1)))
#include "regional_stats_inc/regional_field_nd_to_n-1d-2.F"
end subroutine copy_field_2r!}}}


subroutine copy_field_3r(block, inname, outname, elem_name)!{{{
  type (block_type), pointer, intent(inout) :: block
  character (len=StrKIND), intent(in) :: inname, outname, elem_name

  type (field3DReal), pointer :: src
  type (field2DReal), pointer :: dst
  type (mpas_pool_type), pointer :: pool
  integer :: i
  integer, dimension(3) :: src_dims

#include "regional_stats_inc/regional_field_nd_to_n-1d-1.F"
  allocate(dst % array(src_dims(1), src_dims(2)))
#include "regional_stats_inc/regional_field_nd_to_n-1d-2.F"
end subroutine copy_field_3r!}}}


subroutine copy_field_4r(block, inname, outname, elem_name)!{{{
  type (block_type), pointer, intent(inout) :: block
  character (len=StrKIND), intent(in) :: inname, outname, elem_name

  type (field4DReal), pointer :: src
  type (field3DReal), pointer :: dst
  type (mpas_pool_type), pointer :: pool
  integer :: i
  integer, dimension(4) :: src_dims

#include "regional_stats_inc/regional_field_nd_to_n-1d-1.F"
  allocate(dst % array(src_dims(1), src_dims(2), src_dims(3)))
#include "regional_stats_inc/regional_field_nd_to_n-1d-2.F"
end subroutine copy_field_4r!}}}


subroutine copy_field_5r(block, inname, outname, elem_name)!{{{
  type (block_type), pointer, intent(inout) :: block
  character (len=StrKIND), intent(in) :: inname, outname, elem_name 

  type (field5DReal), pointer :: src
  type (field4DReal), pointer :: dst
  type (mpas_pool_type), pointer :: pool
  integer :: i
  integer, dimension(5) :: src_dims

#include "regional_stats_inc/regional_field_nd_to_n-1d-1.F"
  allocate(dst % array(src_dims(1), src_dims(2), src_dims(3), src_dims(4)))
#include "regional_stats_inc/regional_field_nd_to_n-1d-2.F"
end subroutine copy_field_5r!}}}


subroutine copy_field_1i(block, inname, outname, elem_name)!{{{
  type (block_type), pointer, intent(inout) :: block
  character (len=StrKIND), intent(in) :: inname, outname, elem_name

  type (field1DInteger), pointer :: src
  type (field0DReal), pointer :: dst
  type (mpas_pool_type), pointer :: pool

#include "regional_stats_inc/regional_field_1d_to_0d.F"
end subroutine copy_field_1i!}}}


subroutine copy_field_2i(block, inname, outname, elem_name)!{{{
  type (block_type), pointer, intent(inout) :: block
  character (len=StrKIND), intent(in) :: inname, outname, elem_name

  type (field2DInteger), pointer :: src
  type (field1DReal), pointer :: dst
  type (mpas_pool_type), pointer :: pool
  integer :: i
  integer, dimension(2) :: src_dims

#include "regional_stats_inc/regional_field_nd_to_n-1d-1.F"
  allocate(dst % array(src_dims(1)))
#include "regional_stats_inc/regional_field_nd_to_n-1d-2.F"
end subroutine copy_field_2i!}}}


subroutine copy_field_3i(block, inname, outname, elem_name)!{{{
  type (block_type), pointer, intent(inout) :: block
  character (len=StrKIND), intent(in) :: inname, outname, elem_name

  type (field3DInteger), pointer :: src
  type (field2DReal), pointer :: dst
  type (mpas_pool_type), pointer :: pool
  integer :: i
  integer, dimension(3) :: src_dims

#include "regional_stats_inc/regional_field_nd_to_n-1d-1.F"
  allocate(dst % array(src_dims(1), src_dims(2)))
#include "regional_stats_inc/regional_field_nd_to_n-1d-2.F"
end subroutine copy_field_3i!}}}



!***********************************************************************
! routine operateX_Y
!
!> \brief Series of subroutines to support operations on run-time types
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  These subroutines encapsulate the different operations that can occur
!>  based on the run-time types. (This would likely be
!>  instantiated generics/templates in other languages.)
!>
!>  We will do the reductions using slicing operations rather
!>  than for loops, and then copy to and from a flattened array,
!>  during distributed communication.
!>
!>  This will be slower than for loops using flattened array output
!>  indexing on the first copy, because we're going to have to copy twice
!>  on the distributed communication. (Because the framework reduceAll
!>  only works for rank-1 arrays.) The reason to do this is that
!>  to preserve the dimensionality semantics for legibility and
!>  for fewer bugs in the reduction. Basically, we flatten and unflatten,
!>  and could remove the flatten, and reduce in a for loop with indexing
!>  on a flattened array, and then reshape. While it could be possible to
!>  do manual for loops without slicing, and increasing the code,
!>  I doubt the performance gain for doing one less copy of a small
!>  array will be that significant.
!>
!>  Also, there are inner if statements, per block, to switch based
!>  on the weighting function and weighting field type. This could
!>  be moved to an outer loop in typed_operate, but that would increase
!>  the number of specialized functions by a factor of 5.
!-----------------------------------------------------------------------

!!
!! first
!!
!! copy this one for making 3->2, 1->0 avg
!! copy this one for making 2->1 min
!!
subroutine operate2r_avg (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flattened, rather than out_array because it's 1D target already
  allocate(flattened(size(out_array)))
  ! zero the array
  flattened = 0

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flattened = flattened + in_array(:, i)
#include "regional_stats_inc/regional_op_mulreal.F"
  flattened = flattened + in_array(:, i) * real_weights(i)
#include "regional_stats_inc/regional_op_divreal.F"
  flattened = flattened + in_array(:, i) / real_weights(i)
#include "regional_stats_inc/regional_op_mulint.F"
  flattened = flattened + in_array(:, i) * integer_weights(i)
#include "regional_stats_inc/regional_op_divint.F"
  flattened = flattened + in_array(:, i) / integer_weights(i)
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! sum across processors and divide by total
  ! we can directly copy into out_array since it is 1D already
  call mpas_dmpar_sum_real_array(dminfo, size(flattened), &
    flattened, out_array)
  out_array = out_array / regions % counts(m)
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate2r_avg


!!
!! copied from 2->1 avg
!!
!! copy this one for making 4->3, 5->4 avg
!! copy this one for making 3->2 min
!!
subroutine operate3r_avg (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = 0

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = out_array + in_array(:, :, i)
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = out_array + in_array(:, :, i) * real_weights(i)
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = out_array + in_array(:, :, i) / real_weights(i)
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = out_array + in_array(:, :, i) * integer_weights(i)
#include "regional_stats_inc/regional_op_divint.F"
  out_array = out_array + in_array(:, :, i) / integer_weights(i)
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! sum across processors and divide by total
  allocate(flattened(size(out_array)))
  call mpas_dmpar_sum_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  flattened = flattened / regions % counts(m)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate3r_avg


!!
!! copied from 3->2 avg
!!
!! copy this one for making 4->3 min
!!
subroutine operate4r_avg (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = 0

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = out_array + in_array(:, :, :, i)
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = out_array + in_array(:, :, :, i) * real_weights(i)
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = out_array + in_array(:, :, :, i) / real_weights(i)
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = out_array + in_array(:, :, :, i) * integer_weights(i)
#include "regional_stats_inc/regional_op_divint.F"
  out_array = out_array + in_array(:, :, :, i) / integer_weights(i)
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! sum across processors and divide by total
  allocate(flattened(size(out_array)))
  call mpas_dmpar_sum_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  flattened = flattened / regions % counts(m)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate4r_avg


!!
!! copied from 3->2 avg
!!
!! copy this one for making 5->4 min
!!
subroutine operate5r_avg (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = 0

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = out_array + in_array(:, :, :, :, i)
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = out_array + in_array(:, :, :, :, i) * real_weights(i)
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = out_array + in_array(:, :, :, :, i) / real_weights(i)
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = out_array + in_array(:, :, :, :, i) * integer_weights(i)
#include "regional_stats_inc/regional_op_divint.F"
  out_array = out_array + in_array(:, :, :, :, i) / integer_weights(i)
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! sum across processors and divide by total
  allocate(flattened(size(out_array)))
  call mpas_dmpar_sum_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  flattened = flattened / regions % counts(m)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate5r_avg


!!
!! copied from 2->1 avg
!!
!! copy this one for making 1->0 min
!!
subroutine operate1r_avg (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:), pointer :: in_array
  real (kind=RKIND), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flat_scalar, rather than out_array because it's 0D target already
  ! zero the array
  flat_scalar = 0

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flat_scalar = flat_scalar + in_array(i)
#include "regional_stats_inc/regional_op_mulreal.F"
  flat_scalar = flat_scalar + in_array(i) * real_weights(i)
#include "regional_stats_inc/regional_op_divreal.F"
  flat_scalar = flat_scalar + in_array(i) / real_weights(i)
#include "regional_stats_inc/regional_op_mulint.F"
  flat_scalar = flat_scalar + in_array(i) * integer_weights(i)
#include "regional_stats_inc/regional_op_divint.F"
  flat_scalar = flat_scalar + in_array(i) / integer_weights(i)
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! sum across processors and divide by total
  ! we can directly copy into out_array since it is 0D already
  call mpas_dmpar_sum_real(dminfo, flat_scalar, out_array)
  out_array = out_array / regions % counts(m)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate1r_avg


!!
!! copied from 2->1 avg
!!
!! copy this one for making 2->1 max
!!
subroutine operate2r_min (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flattened, rather than out_array because it's 1D target already
  allocate(flattened(size(out_array)))
  ! zero the array
  flattened = DEFAULT_MPAS_MAX_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flattened = min(flattened, in_array(:, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  flattened = min(flattened, in_array(:, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  flattened = min(flattened, in_array(:, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  flattened = min(flattened, in_array(:, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  flattened = min(flattened, in_array(:, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! min across processors
  ! we can directly copy into out_array since it is 1D already
  call mpas_dmpar_min_real_array(dminfo, size(flattened), &
    flattened, out_array)
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate2r_min


!!
!! copied from 3->2 avg
!!
!! copy this one for making 3->2 max
!!
subroutine operate3r_min (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MAX_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = min(out_array, in_array(:, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = min(out_array, in_array(:, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = min(out_array, in_array(:, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = min(out_array, in_array(:, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = min(out_array, in_array(:, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! min across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_min_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate3r_min


!!
!! copied from 4->3 avg
!!
!! copy this one for making 4->3 max
!!
subroutine operate4r_min (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MAX_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = min(out_array, in_array(:, :, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = min(out_array, in_array(:, :, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = min(out_array, in_array(:, :, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = min(out_array, in_array(:, :, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = min(out_array, in_array(:, :, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! min across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_min_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate4r_min


!!
!! copied from 5->4 avg
!!
!! copy this one for making 5->4 max
!!
subroutine operate5r_min (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MAX_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = min(out_array, in_array(:, :, :, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = min(out_array, in_array(:, :, :, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = min(out_array, in_array(:, :, :, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = min(out_array, in_array(:, :, :, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = min(out_array, in_array(:, :, :, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! min across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_min_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate5r_min


!!
!! copied from 1->0 avg
!!
!! copy this one for making 1->0 max
!!
subroutine operate1r_min (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:), pointer :: in_array
  real (kind=RKIND), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flat_scalar, rather than out_array because it's 0D target already
  ! zero the array
  flat_scalar = DEFAULT_MPAS_MAX_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flat_scalar = min(flat_scalar, in_array(i))
#include "regional_stats_inc/regional_op_mulreal.F"
  flat_scalar = min(flat_scalar, in_array(i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  flat_scalar = min(flat_scalar, in_array(i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  flat_scalar = min(flat_scalar, in_array(i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  flat_scalar = min(flat_scalar, in_array(i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! min across processors
  ! we can directly copy into out_array since it is 0D already
  call mpas_dmpar_min_real(dminfo, flat_scalar, out_array)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate1r_min


!!
!! copied from 2->1 min
!!
subroutine operate2r_max (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flattened, rather than out_array because it's 1D target already
  allocate(flattened(size(out_array)))
  ! zero the array
  flattened = DEFAULT_MPAS_MIN_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flattened = max(flattened, in_array(:, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  flattened = max(flattened, in_array(:, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  flattened = max(flattened, in_array(:, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  flattened = max(flattened, in_array(:, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  flattened = max(flattened, in_array(:, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! max across processors
  ! we can directly copy into out_array since it is 1D already
  call mpas_dmpar_max_real_array(dminfo, size(flattened), &
    flattened, out_array)
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate2r_max


!!
!! copied from 3->2 min
!!
subroutine operate3r_max (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MIN_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = max(out_array, in_array(:, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = max(out_array, in_array(:, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = max(out_array, in_array(:, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = max(out_array, in_array(:, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = max(out_array, in_array(:, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! max across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_max_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate3r_max


!!
!! copied from 4->3 min
!!
subroutine operate4r_max (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MIN_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = max(out_array, in_array(:, :, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = max(out_array, in_array(:, :, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = max(out_array, in_array(:, :, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = max(out_array, in_array(:, :, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = max(out_array, in_array(:, :, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! max across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_max_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate4r_max


!!
!! copied from 5->4 min
!!
subroutine operate5r_max (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MIN_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = max(out_array, in_array(:, :, :, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = max(out_array, in_array(:, :, :, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = max(out_array, in_array(:, :, :, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = max(out_array, in_array(:, :, :, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = max(out_array, in_array(:, :, :, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! max across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_max_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate5r_max


!!
!! copied from 1->0 min
!!
subroutine operate1r_max (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:), pointer :: in_array
  real (kind=RKIND), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flat_scalar, rather than out_array because it's 0D target already
  ! zero the array
  flat_scalar = DEFAULT_MPAS_MIN_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flat_scalar = max(flat_scalar, in_array(i))
#include "regional_stats_inc/regional_op_mulreal.F"
  flat_scalar = max(flat_scalar, in_array(i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  flat_scalar = max(flat_scalar, in_array(i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  flat_scalar = max(flat_scalar, in_array(i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  flat_scalar = max(flat_scalar, in_array(i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! max across processors
  ! we can directly copy into out_array since it is 0D already
  call mpas_dmpar_max_real(dminfo, flat_scalar, out_array)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate1r_max

!!
!! all the integer versions are copied from real versions
!!

subroutine operate2i_avg (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flattened, rather than out_array because it's 1D target already
  allocate(flattened(size(out_array)))
  ! zero the array
  flattened = 0

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flattened = flattened + in_array(:, i)
#include "regional_stats_inc/regional_op_mulreal.F"
  flattened = flattened + in_array(:, i) * real_weights(i)
#include "regional_stats_inc/regional_op_divreal.F"
  flattened = flattened + in_array(:, i) / real_weights(i)
#include "regional_stats_inc/regional_op_mulint.F"
  flattened = flattened + in_array(:, i) * integer_weights(i)
#include "regional_stats_inc/regional_op_divint.F"
  flattened = flattened + in_array(:, i) / integer_weights(i)
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! sum across processors and divide by total
  ! we can directly copy into out_array since it is 1D already
  call mpas_dmpar_sum_real_array(dminfo, size(flattened), &
    flattened, out_array)
  out_array = out_array / regions % counts(m)
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate2i_avg


subroutine operate3i_avg (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = 0

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = out_array + in_array(:, :, i)
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = out_array + in_array(:, :, i) * real_weights(i)
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = out_array + in_array(:, :, i) / real_weights(i)
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = out_array + in_array(:, :, i) * integer_weights(i)
#include "regional_stats_inc/regional_op_divint.F"
  out_array = out_array + in_array(:, :, i) / integer_weights(i)
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! sum across processors and divide by total
  allocate(flattened(size(out_array)))
  call mpas_dmpar_sum_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  flattened = flattened / regions % counts(m)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate3i_avg


subroutine operate4i_avg (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = 0

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = out_array + in_array(:, :, :, i)
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = out_array + in_array(:, :, :, i) * real_weights(i)
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = out_array + in_array(:, :, :, i) / real_weights(i)
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = out_array + in_array(:, :, :, i) * integer_weights(i)
#include "regional_stats_inc/regional_op_divint.F"
  out_array = out_array + in_array(:, :, :, i) / integer_weights(i)
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! sum across processors and divide by total
  allocate(flattened(size(out_array)))
  call mpas_dmpar_sum_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  flattened = flattened / regions % counts(m)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate4i_avg


subroutine operate5i_avg (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = 0

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = out_array + in_array(:, :, :, :, i)
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = out_array + in_array(:, :, :, :, i) * real_weights(i)
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = out_array + in_array(:, :, :, :, i) / real_weights(i)
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = out_array + in_array(:, :, :, :, i) * integer_weights(i)
#include "regional_stats_inc/regional_op_divint.F"
  out_array = out_array + in_array(:, :, :, :, i) / integer_weights(i)
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! sum across processors and divide by total
  allocate(flattened(size(out_array)))
  call mpas_dmpar_sum_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  flattened = flattened / regions % counts(m)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate5i_avg


subroutine operate1i_avg (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:), pointer :: in_array
  real (kind=RKIND), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flat_scalar, rather than out_array because it's 0D target already
  ! zero the array
  flat_scalar = 0

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flat_scalar = flat_scalar + in_array(i)
#include "regional_stats_inc/regional_op_mulreal.F"
  flat_scalar = flat_scalar + in_array(i) * real_weights(i)
#include "regional_stats_inc/regional_op_divreal.F"
  flat_scalar = flat_scalar + in_array(i) / real_weights(i)
#include "regional_stats_inc/regional_op_mulint.F"
  flat_scalar = flat_scalar + in_array(i) * integer_weights(i)
#include "regional_stats_inc/regional_op_divint.F"
  flat_scalar = flat_scalar + in_array(i) / integer_weights(i)
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! sum across processors and divide by total
  ! we can directly copy into out_array since it is 0D already
  call mpas_dmpar_sum_real(dminfo, flat_scalar, out_array)
  out_array = out_array / regions % counts(m)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate1i_avg


subroutine operate2i_min (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flattened, rather than out_array because it's 1D target already
  allocate(flattened(size(out_array)))
  ! zero the array
  flattened = DEFAULT_MPAS_MAX_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flattened = min(flattened, in_array(:, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  flattened = min(flattened, in_array(:, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  flattened = min(flattened, in_array(:, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  flattened = min(flattened, in_array(:, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  flattened = min(flattened, in_array(:, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! min across processors
  ! we can directly copy into out_array since it is 1D already
  call mpas_dmpar_min_real_array(dminfo, size(flattened), &
    flattened, out_array)
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate2i_min


subroutine operate3i_min (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MAX_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = min(out_array, in_array(:, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = min(out_array, in_array(:, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = min(out_array, in_array(:, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = min(out_array, in_array(:, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = min(out_array, in_array(:, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! min across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_min_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate3i_min


subroutine operate4i_min (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MAX_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = min(out_array, in_array(:, :, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = min(out_array, in_array(:, :, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = min(out_array, in_array(:, :, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = min(out_array, in_array(:, :, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = min(out_array, in_array(:, :, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! min across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_min_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate4i_min


subroutine operate5i_min (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MAX_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = min(out_array, in_array(:, :, :, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = min(out_array, in_array(:, :, :, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = min(out_array, in_array(:, :, :, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = min(out_array, in_array(:, :, :, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = min(out_array, in_array(:, :, :, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! min across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_min_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate5i_min


subroutine operate1i_min (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:), pointer :: in_array
  real (kind=RKIND), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flat_scalar, rather than out_array because it's 0D target already
  ! zero the array
  flat_scalar = DEFAULT_MPAS_MAX_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flat_scalar = min(flat_scalar, in_array(i))
#include "regional_stats_inc/regional_op_mulreal.F"
  flat_scalar = min(flat_scalar, in_array(i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  flat_scalar = min(flat_scalar, in_array(i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  flat_scalar = min(flat_scalar, in_array(i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  flat_scalar = min(flat_scalar, in_array(i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! min across processors
  ! we can directly copy into out_array since it is 0D already
  call mpas_dmpar_min_real(dminfo, flat_scalar, out_array)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate1i_min


subroutine operate2i_max (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flattened, rather than out_array because it's 1D target already
  allocate(flattened(size(out_array)))
  ! zero the array
  flattened = DEFAULT_MPAS_MIN_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flattened = max(flattened, in_array(:, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  flattened = max(flattened, in_array(:, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  flattened = max(flattened, in_array(:, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  flattened = max(flattened, in_array(:, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  flattened = max(flattened, in_array(:, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! max across processors
  ! we can directly copy into out_array since it is 1D already
  call mpas_dmpar_max_real_array(dminfo, size(flattened), &
    flattened, out_array)
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate2i_max


subroutine operate3i_max (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MIN_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = max(out_array, in_array(:, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = max(out_array, in_array(:, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = max(out_array, in_array(:, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = max(out_array, in_array(:, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = max(out_array, in_array(:, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! max across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_max_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate3i_max


subroutine operate4i_max (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MIN_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = max(out_array, in_array(:, :, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = max(out_array, in_array(:, :, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = max(out_array, in_array(:, :, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = max(out_array, in_array(:, :, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = max(out_array, in_array(:, :, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! max across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_max_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate4i_max


subroutine operate5i_max (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! zero the array
  out_array = DEFAULT_MPAS_MIN_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  out_array = max(out_array, in_array(:, :, :, :, i))
#include "regional_stats_inc/regional_op_mulreal.F"
  out_array = max(out_array, in_array(:, :, :, :, i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  out_array = max(out_array, in_array(:, :, :, :, i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  out_array = max(out_array, in_array(:, :, :, :, i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  out_array = max(out_array, in_array(:, :, :, :, i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! max across processors
  allocate(flattened(size(out_array)))
  call mpas_dmpar_max_real_array(dminfo, size(flattened), &
    reshape(out_array, shape(flattened)), flattened)
  out_array = reshape(flattened, shape(out_array))
  deallocate(flattened)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate5i_max


subroutine operate1i_max (dminfo, start_block, regions, variable)

! declarations for each version
#include "regional_stats_inc/regional_op_header.F"
! dimension specific local declarations

  real (kind=RKIND), dimension(:), pointer :: in_array
  real (kind=RKIND), pointer :: out_array

! iterate over masks
#include "regional_stats_inc/regional_op_1.F"
! do your mask initializations

  ! target flat_scalar, rather than out_array because it's 0D target already
  ! zero the array
  flat_scalar = DEFAULT_MPAS_MIN_VAL

! iterate over blocks
#include "regional_stats_inc/regional_op_2.F"
! do your kernel per block

! operations per masked element
#include "regional_stats_inc/regional_op_id.F"
  flat_scalar = max(flat_scalar, in_array(i))
#include "regional_stats_inc/regional_op_mulreal.F"
  flat_scalar = max(flat_scalar, in_array(i) * real_weights(i))
#include "regional_stats_inc/regional_op_divreal.F"
  flat_scalar = max(flat_scalar, in_array(i) / real_weights(i))
#include "regional_stats_inc/regional_op_mulint.F"
  flat_scalar = max(flat_scalar, in_array(i) * integer_weights(i))
#include "regional_stats_inc/regional_op_divint.F"
  flat_scalar = max(flat_scalar, in_array(i) / integer_weights(i))
#include "regional_stats_inc/regional_op_end.F"

! close the block loop (started in regional_op_2.F)
#include "regional_stats_inc/regional_op_3.F"
! do your distributed communication per mask

  ! max across processors
  ! we can directly copy into out_array since it is 0D already
  call mpas_dmpar_max_real(dminfo, flat_scalar, out_array)

! close the mask loop (started in regional_op_1.F)
#include "regional_stats_inc/regional_op_4.F"
end subroutine operate1i_max


end module ocn_regional_stats
! vim: foldmethod=marker
