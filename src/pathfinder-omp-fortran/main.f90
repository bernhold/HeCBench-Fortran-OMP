program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use omp_lib
  implicit none

  integer, parameter :: m_seed = 9
  integer :: rows, cols, pyramid_height, total, row, col
  integer :: left_col, right_col, active_parity, mismatches
  integer(int32), allocatable :: wall(:), path_a(:), path_b(:), reference(:)
  real(real64) :: offload_start, offload_end, kernel_start, kernel_end

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    integer(c_int) function c_rand() bind(C, name="rand")
      import :: c_int
    end function c_rand
  end interface

  if (command_argument_count() == 3) then
    cols = read_arg(1)
    rows = read_arg(2)
    pyramid_height = read_arg(3)
  else
    print '(A)', 'Usage: ./main <column length> <row length> <pyramid_height>'
    stop 0
  end if

  if (cols <= 0 .or. rows <= 0 .or. pyramid_height <= 0) then
    error stop 'pathfinder dimensions and pyramid_height must be positive'
  end if

  total = rows * cols
  allocate(wall(total), path_a(cols), path_b(cols), reference(cols))
  call initialize_wall(wall, rows, cols)
  path_a = wall(1:cols)
  path_b = 0_int32
  reference = path_a

  offload_start = omp_get_wtime()
  kernel_start = omp_get_wtime()

  !$omp target data map(to: wall(1:total)) map(tofrom: path_a(1:cols), path_b(1:cols))
  do row = 2, rows
    if (mod(row, 2) == 0) then
      !$omp target teams distribute parallel do thread_limit(256) private(left_col, right_col)
      do col = 1, cols
        left_col = max(1, col - 1)
        right_col = min(cols, col + 1)
        path_b(col) = wall((row - 1) * cols + col) + min(path_a(left_col), min(path_a(col), path_a(right_col)))
      end do
      !$omp end target teams distribute parallel do
    else
      !$omp target teams distribute parallel do thread_limit(256) private(left_col, right_col)
      do col = 1, cols
        left_col = max(1, col - 1)
        right_col = min(cols, col + 1)
        path_a(col) = wall((row - 1) * cols + col) + min(path_b(left_col), min(path_b(col), path_b(right_col)))
      end do
      !$omp end target teams distribute parallel do
    end if
  end do
  !$omp end target data

  kernel_end = omp_get_wtime()
  offload_end = omp_get_wtime()

  call cpu_reference(wall, rows, cols, reference)
  active_parity = mod(rows, 2)
  if (active_parity == 1) then
    mismatches = count(path_a /= reference)
  else
    mismatches = count(path_b /= reference)
  end if
  if (mismatches /= 0) then
    write(*, '(A,I0)') 'FAIL: mismatches=', mismatches
    error stop 'pathfinder validation failed'
  end if

  write(*, '(A,F0.6,A)') 'Total kernel execution time: ', kernel_end - kernel_start, ' (s)'
  write(*, '(A,F0.6,A)') 'Device offloading time = ', offload_end - offload_start, ' (s)'

  deallocate(wall, path_a, path_b, reference)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine initialize_wall(wall, rows, cols)
    integer(int32), intent(out) :: wall(:)
    integer, intent(in) :: rows, cols
    integer :: i, j, idx
    call c_srand(int(m_seed, c_int))
    do j = 1, rows
      do i = 1, cols
        idx = (j - 1) * cols + i
        wall(idx) = int(mod(c_rand(), 10_c_int), int32)
      end do
    end do
  end subroutine initialize_wall

  subroutine cpu_reference(wall, rows, cols, result)
    integer(int32), intent(in) :: wall(:)
    integer, intent(in) :: rows, cols
    integer(int32), intent(out) :: result(:)
    integer(int32), allocatable :: scratch(:)
    integer :: row, col, left_col, right_col

    allocate(scratch(cols))
    result = wall(1:cols)
    scratch = 0_int32
    do row = 2, rows
      do col = 1, cols
        left_col = max(1, col - 1)
        right_col = min(cols, col + 1)
        scratch(col) = wall((row - 1) * cols + col) + min(result(left_col), min(result(col), result(right_col)))
      end do
      result = scratch
    end do
    deallocate(scratch)
  end subroutine cpu_reference

end program main
