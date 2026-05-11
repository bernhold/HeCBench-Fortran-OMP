module michalewicz_mod
  use iso_fortran_env, only: int32, int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  real(real32), parameter :: pi32 = acos(-1.0_real32)

contains

  pure integer(int32) function next_rand(state) result(value)
    integer(int32), intent(inout) :: state
    integer(int64) :: tmp

    tmp = mod(1664525_int64 * int(state, int64) + 1013904223_int64, 2147483648_int64)
    state = int(tmp, int32)
    value = iand(state, int(z'7fffffff', int32))
  end function next_rand

  subroutine fill_values(values)
    real(real32), intent(out) :: values(0:)
    integer(int32) :: state
    integer(int64) :: i

    state = 19937_int32
    do i = 0, int(size(values), int64) - 1
      values(i) = 4.0_real32 * real(next_rand(state), real32) / real(huge(1_int32), real32)
    end do
  end subroutine fill_values

  pure real(real32) function michalewicz_cpu(values, offset, dim) result(value)
    real(real32), intent(in) :: values(0:)
    integer(int64), intent(in) :: offset
    integer, intent(in) :: dim
    integer :: k
    real(real32) :: x, a, b, c

    value = 0.0_real32
    do k = 0, dim - 1
      x = values(offset + k)
      a = sin(x)
      b = sin((real(k + 1, real32) * x * x) / pi32)
      c = b ** 20
      value = value + a * c
    end do
    value = -value
  end function michalewicz_cpu

  real(real32) function cpu_minimum(values, n, dim) result(min_value)
    real(real32), intent(in) :: values(0:)
    integer(int64), intent(in) :: n
    integer, intent(in) :: dim
    integer(int64) :: j

    min_value = 0.0_real32
    do j = 0, n - 1
      min_value = min(min_value, michalewicz_cpu(values, j * dim, dim))
    end do
  end function cpu_minimum

  subroutine print_error(value, dim)
    real(real32), intent(in) :: value
    integer, intent(in) :: dim
    real(real32) :: true_min
    character(len=32) :: value_text, error_text

    value_text = fixed6(value)
    write(*,'(A,A)') 'Global minima = ', trim(value_text)
    true_min = 0.0_real32
    if (dim == 2) then
      true_min = -1.8013_real32
    else if (dim == 5) then
      true_min = -4.687658_real32
    else if (dim == 10) then
      true_min = -9.66015_real32
    end if
    error_text = fixed6(abs(true_min - value))
    write(*,'(A,A)') 'Error = ', trim(error_text)
  end subroutine print_error

  pure character(len=32) function fixed6(value) result(text)
    real(real32), intent(in) :: value
    character(len=32) :: tmp

    write(tmp, '(F0.6)') value
    tmp = adjustl(tmp)
    if (tmp(1:1) == '.') then
      text = '0' // trim(tmp)
    else if (tmp(1:2) == '-.') then
      text = '-0' // trim(tmp(2:))
    else
      text = trim(tmp)
    end if
  end function fixed6

  subroutine run_dimension(n, repeat, dim)
    integer(int64), intent(in) :: n
    integer, intent(in) :: repeat, dim
    integer(int64) :: total_size, j
    integer :: iter, k
    real(real32), allocatable :: values(:)
    real(real32) :: min_value, reference_value, x, a, b, c, candidate
    real(real64) :: start_time, elapsed_ns

    total_size = n * int(dim, int64)
    allocate(values(0:total_size - 1))
    call fill_values(values)

    min_value = 0.0_real32

    !$omp target data map(to: values(0:total_size - 1)) map(tofrom: min_value)
    start_time = omp_get_wtime()
    do iter = 1, repeat
      !$omp target teams distribute parallel do thread_limit(block_size) reduction(min: min_value) &
      !$omp& private(k, x, a, b, c, candidate)
      do j = 0, n - 1
        candidate = 0.0_real32
        do k = 0, dim - 1
          x = values(j * dim + k)
          a = sin(x)
          b = sin((real(k + 1, real32) * x * x) / pi32)
          c = b ** 20
          candidate = candidate + a * c
        end do
        candidate = -candidate
        min_value = min(min_value, candidate)
      end do
      !$omp end target teams distribute parallel do
    end do
    elapsed_ns = (omp_get_wtime() - start_time) * 1.0e9_real64
    write(*,'(A,I0,A,F0.6,A)') 'Average execution time of kernel (dim = ', dim, '): ', &
      elapsed_ns * 1.0e-3_real64 / real(repeat, real64), ' (us)'
    !$omp end target data

    reference_value = cpu_minimum(values, n, dim)
    if (abs(reference_value - min_value) > 1.0e-4_real32) then
      write(*,'(A,I0,A,F0.6,A,F0.6)') 'FAIL dim ', dim, ': CPU reference ', reference_value, &
        ' GPU result ', min_value
      stop 2
    end if

    call print_error(min_value, dim)
    deallocate(values)
  end subroutine run_dimension

end module michalewicz_mod

program main
  use iso_fortran_env, only: int64
  use michalewicz_mod
  implicit none

  integer :: argc, repeat, status, d
  integer(int64) :: n
  integer, parameter :: dims(3) = [2, 5, 10]
  character(len=256) :: arg, prog

  argc = command_argument_count()
  if (argc /= 2) then
    call get_command_argument(0, prog)
    write(*,'(A,A,A)') 'Usage: ', trim(prog), ' <number of vectors> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *, iostat=status) n
  if (status /= 0) stop 1
  call get_command_argument(2, arg)
  read(arg, *, iostat=status) repeat
  if (status /= 0) stop 1

  do d = 1, size(dims)
    call run_dimension(n, repeat, dims(d))
  end do
end program main
