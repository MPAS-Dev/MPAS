  type (dm_info), pointer, intent(in) :: dminfo
  type (block_type), pointer, intent(in) :: start_block
  type (regional_type), intent(in) :: regions
  type (regional_variable_type), intent(in) :: variable

  integer :: b, m, i, last
  integer, pointer :: solve
  type (block_type), pointer :: block
  type (mpas_pool_type), pointer :: amPool, maskPool
  real (kind=RKIND), dimension(:), pointer :: real_weights
  integer, dimension(:), pointer :: integer_weights
  integer, dimension(:,:), pointer :: mask
  real (kind=RKIND), dimension(:), allocatable :: flattened
  real (kind=RKIND) :: flat_scalar
