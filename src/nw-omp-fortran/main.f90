program nw_omp_fortran
  use, intrinsic :: iso_fortran_env, only: int32, int64, real64
  use omp_lib, only: omp_get_wtime
  implicit none

  integer, parameter :: block_size = 16
  integer, parameter :: warmup = 100
  integer, parameter :: blosum62(0:23,0:23) = reshape([ &
     4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4, &
    -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4, &
    -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4, &
    -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4, &
     0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4, &
    -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4, &
    -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4, &
     0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4, &
    -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4, &
    -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4, &
    -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4, &
    -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4, &
    -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4, &
    -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4, &
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4, &
     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4, &
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4, &
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4, &
    -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4, &
     0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4, &
    -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4, &
    -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4, &
     0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4, &
    -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1], &
    [24, 24], order=[2, 1])

  integer :: argc, dim, penalty, repeat, ios
  integer :: max_rows, max_cols, n_items, i, j
  integer(int32), allocatable :: reference(:), host_items(:), device_items(:)
  character(len=128) :: arg, program_name
  real(real64) :: start_time, end_time, total_time
  logical :: ok

  write(*, '(A,I0,A)') 'WG size of kernel = ', block_size, ' '

  argc = command_argument_count()
  if (argc /= 3) call usage()

  call get_command_argument(1, arg)
  read(arg, *, iostat=ios) dim
  if (ios /= 0) call usage()
  call get_command_argument(2, arg)
  read(arg, *, iostat=ios) penalty
  if (ios /= 0) call usage()
  call get_command_argument(3, arg)
  read(arg, *, iostat=ios) repeat
  if (ios /= 0 .or. repeat <= 0) call usage()

  if (mod(dim, block_size) /= 0) then
    write(0, '(A)') 'The dimension values must be a multiple of 16'
    stop 1
  end if

  max_rows = dim + 1
  max_cols = dim + 1
  n_items = max_rows * max_cols
  allocate(reference(0:n_items - 1), host_items(0:n_items - 1), device_items(0:n_items - 1))
  reference = 0_int32
  host_items = 0_int32
  device_items = 0_int32

  call init_inputs(reference, host_items, device_items, max_rows, max_cols, penalty)

  !$omp target data map(tofrom: device_items(0:n_items - 1)) map(to: reference(0:n_items - 1))
    do i = 1, warmup + repeat
      if (i == warmup + 1) start_time = omp_get_wtime()
      call nw_device(device_items, reference, max_cols, penalty)
    end do
    end_time = omp_get_wtime()
  !$omp end target data

  total_time = (end_time - start_time) / real(repeat, real64)
  write(*, '(A,F8.6,A)') 'Total kernel execution time: ', total_time, ' (s)'

  call nw_host(host_items, reference, max_cols, penalty)
  ok = all(device_items == host_items)
  if (ok) then
    write(*, '(A)') 'PASS'
  else
    write(*, '(A)') 'FAIL'
    stop 1
  end if

  deallocate(reference, host_items, device_items)

contains

  subroutine usage()
    call get_command_argument(0, program_name)
    write(0, '(A,A,A)') 'Usage: ', trim(program_name), ' <max_rows/max_cols> <penalty> <repeat>'
    write(0, '(A)') char(9)//'<dimension>  - x and y dimensions'
    write(0, '(A)') char(9)//'<penalty> - penalty(positive integer)'
    write(0, '(A)') char(9)//'<repeat> - the number of kernel executions'
    write(0, '(A)') char(9)//'<file> - filename'
    stop 1
  end subroutine usage

  subroutine init_inputs(reference, host_items, device_items, max_rows, max_cols, penalty)
    integer(int32), intent(out) :: reference(0:), host_items(0:), device_items(0:)
    integer, intent(in) :: max_rows, max_cols, penalty
    integer :: i, j, seed

    seed = 7
    reference = 0_int32
    host_items = 0_int32
    device_items = 0_int32

    do i = 1, max_rows - 1
      host_items(i * max_cols) = next_rand10(seed)
      device_items(i * max_cols) = host_items(i * max_cols)
    end do

    do j = 1, max_cols - 1
      host_items(j) = next_rand10(seed)
      device_items(j) = host_items(j)
    end do

    do i = 1, max_cols - 1
      do j = 1, max_rows - 1
        reference(i * max_cols + j) = blosum62(device_items(i * max_cols), device_items(j))
      end do
    end do

    do i = 1, max_rows - 1
      host_items(i * max_cols) = -i * penalty
      device_items(i * max_cols) = host_items(i * max_cols)
    end do
    do j = 1, max_cols - 1
      host_items(j) = -j * penalty
      device_items(j) = host_items(j)
    end do
  end subroutine init_inputs

  integer(int32) function next_rand10(seed)
    integer, intent(inout) :: seed
    seed = int(mod(1103515245_int64 * int(seed, int64) + 12345_int64, 2147483648_int64))
    next_rand10 = int(mod(seed / 65536, 10) + 1, int32)
  end function next_rand10

  subroutine nw_device(items, reference, max_cols, penalty)
    integer(int32), intent(inout) :: items(0:)
    integer(int32), intent(in) :: reference(0:)
    integer, intent(in) :: max_cols, penalty
    integer :: blk, bx, block_width

    block_width = (max_cols - 1) / block_size

    do blk = 1, block_width
      !$omp target teams distribute parallel do thread_limit(128)
      do bx = 0, blk - 1
        call compute_block(items, reference, max_cols, penalty, bx, blk - 1 - bx)
      end do
      !$omp end target teams distribute parallel do
    end do

    do blk = block_width - 1, 1, -1
      !$omp target teams distribute parallel do thread_limit(128)
      do bx = 0, blk - 1
        call compute_block(items, reference, max_cols, penalty, bx + block_width - blk, block_width - bx - 1)
      end do
      !$omp end target teams distribute parallel do
    end do
  end subroutine nw_device

  subroutine nw_host(items, reference, max_cols, penalty)
    integer(int32), intent(inout) :: items(0:)
    integer(int32), intent(in) :: reference(0:)
    integer, intent(in) :: max_cols, penalty
    integer :: blk, bx, block_width

    block_width = (max_cols - 1) / block_size
    do blk = 1, block_width
      do bx = 0, blk - 1
        call compute_block(items, reference, max_cols, penalty, bx, blk - 1 - bx)
      end do
    end do
    do blk = 2, block_width
      do bx = blk - 1, block_width - 1
        call compute_block(items, reference, max_cols, penalty, bx, block_width + blk - 2 - bx)
      end do
    end do
  end subroutine nw_host

  subroutine compute_block(items, reference, max_cols, penalty, b_index_x, b_index_y)
    integer(int32), intent(inout) :: items(0:)
    integer(int32), intent(in) :: reference(0:)
    integer, intent(in) :: max_cols, penalty, b_index_x, b_index_y
    integer :: ii, jj, global_i, global_j, base
    integer(int32) :: local_items(0:block_size,0:block_size)
    integer(int32) :: local_ref(1:block_size,1:block_size)

    base = max_cols * block_size * b_index_y + block_size * b_index_x
    do ii = 0, block_size
      do jj = 0, block_size
        local_items(ii, jj) = items(base + ii * max_cols + jj)
      end do
    end do

    do ii = 1, block_size
      do jj = 1, block_size
        global_i = b_index_y * block_size + ii
        global_j = b_index_x * block_size + jj
        local_ref(ii, jj) = reference(global_i * max_cols + global_j)
      end do
    end do

    do ii = 1, block_size
      do jj = 1, block_size
        local_items(ii, jj) = max(local_items(ii - 1, jj - 1) + local_ref(ii, jj), &
                                  local_items(ii, jj - 1) - penalty, &
                                  local_items(ii - 1, jj) - penalty)
      end do
    end do

    do ii = 1, block_size
      do jj = 1, block_size
        items(base + ii * max_cols + jj) = local_items(ii, jj)
      end do
    end do
  end subroutine compute_block

end program nw_omp_fortran
