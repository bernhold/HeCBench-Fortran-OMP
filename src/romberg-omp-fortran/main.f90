program main
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  integer, parameter :: row_size = 17
  real(real64), parameter :: lower_limit = 0.0_real64
  real(real64), parameter :: upper_limit = 15.0_real64
  real(real64), parameter :: eps = 1.0e-7_real64

  character(len=128) :: arg1, arg2, arg3
  integer :: nwg, wgs, repeat, iter, k
  real(real64), allocatable :: result(:)
  real(real64) :: start_time, end_time, elapsed_s, d_sum, ref_sum

  if (command_argument_count() /= 3) then
    print '(A)', 'Usage: ./main <number of work-groups> <work-group size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) nwg
  read(arg2, *) wgs
  read(arg3, *) repeat
  if (nwg <= 0 .or. wgs <= 0 .or. repeat <= 0) stop 1

  allocate(result(nwg))
  result = 0.0_real64
  d_sum = 0.0_real64

  !$omp target data map(from: result(1:nwg))
  start_time = omp_get_wtime()
  do iter = 1, repeat
    call compute_segments(result, nwg, wgs)
    !$omp target update from(result(1:nwg))
    d_sum = 0.0_real64
    do k = 1, nwg
      d_sum = d_sum + result(k)
    end do
  end do
  end_time = omp_get_wtime()
  elapsed_s = (end_time - start_time) / real(repeat, real64)
  print '(A,F0.6,A)', 'Average kernel execution time: ', elapsed_s, ' (s)'
  !$omp end target data

  ref_sum = romberg_reference(lower_limit, upper_limit, row_size, eps)
  if (abs(d_sum - ref_sum) > eps) then
    print '(A)', 'FAIL'
    stop 1
  else
    print '(A)', 'PASS'
  end if

  deallocate(result)

contains

  subroutine compute_segments(result, nwg, wgs)
    real(real64), intent(inout) :: result(:)
    integer, intent(in) :: nwg, wgs
    integer :: block
    real(real64) :: width, a, b

    width = (upper_limit - lower_limit) / real(nwg, real64)
    !$omp target teams distribute parallel do num_teams(nwg) thread_limit(wgs) private(block, a, b) firstprivate(nwg, width)
    do block = 0, nwg - 1
      a = lower_limit + real(block, real64) * width
      b = a + width
      result(block + 1) = romberg_reference(a, b, row_size, eps)
    end do
    !$omp end target teams distribute parallel do
  end subroutine compute_segments

  real(real64) function integrand(x)
    real(real64), intent(in) :: x

    integrand = exp(x) * sin(x)
  end function integrand

  real(real64) function romberg_reference(a, b, max_steps, acc)
    real(real64), intent(in) :: a, b, acc
    integer, intent(in) :: max_steps
    real(real64) :: r_prev(row_size), r_curr(row_size)
    real(real64) :: h, c, n_k
    integer :: i, j, ep

    r_prev = 0.0_real64
    r_curr = 0.0_real64
    h = b - a
    r_prev(1) = (integrand(a) + integrand(b)) * h * 0.5_real64

    do i = 2, max_steps
      h = h * 0.5_real64
      c = 0.0_real64
      ep = ishft(1, i - 2)
      do j = 1, ep
        c = c + integrand(a + real(2 * j - 1, real64) * h)
      end do
      r_curr(1) = h * c + 0.5_real64 * r_prev(1)

      do j = 2, i
        n_k = real(ishft(1, 2 * (j - 1)), real64)
        r_curr(j) = (n_k * r_curr(j - 1) - r_prev(j - 1)) / (n_k - 1.0_real64)
      end do

      if (i > 3 .and. abs(r_prev(i - 1) - r_curr(i)) < acc) then
        romberg_reference = r_curr(i - 1)
        return
      end if

      r_prev = r_curr
    end do

    romberg_reference = r_prev(max_steps)
  end function romberg_reference

end program main
