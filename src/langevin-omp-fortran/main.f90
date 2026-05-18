program langevin_main
  use, intrinsic :: iso_fortran_env, only: real32
  use, intrinsic :: ieee_arithmetic
  use omp_lib
  implicit none

  integer :: n, repeat, i
  real(real32), allocatable :: a(:), ref(:), o0(:), o1(:), o2(:)
  real(real32) :: x, x2, x4, x6
  real(real32) :: err0, err1, err2
  real(8) :: start_time, end_time

  call parse_args(n, repeat)

  allocate(a(n), ref(n), o0(n), o1(n), o2(n))

  do i = 1, n
    a(i) = -1.8_real32 + real(i - 1, real32) * (1.79999_real32 / real(n, real32))
  end do

  !$omp target data map(to: a) map(from: o0, o1, o2)
  start_time = omp_get_wtime()
  do i = 1, repeat
    call k0_device(a, o0, n)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') "Average execution time of k0: ", &
       real((end_time - start_time) / dble(repeat), real32), " (s)"

  start_time = omp_get_wtime()
  do i = 1, repeat
    call k1_device(a, o1, n)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') "Average execution time of k1: ", &
       real((end_time - start_time) / dble(repeat), real32), " (s)"

  start_time = omp_get_wtime()
  do i = 1, repeat
    call k2_device(a, o2, n)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') "Average execution time of k2: ", &
       real((end_time - start_time) / dble(repeat), real32), " (s)"
  !$omp end target data

  do i = 1, n
    x = a(i)
    x2 = x * x
    x4 = x2 * x2
    x6 = x4 * x2
    ref(i) = x * (1.0_real32 / 3.0_real32 - x2 / 45.0_real32 + &
         2.0_real32 * x4 / 945.0_real32 - x6 / 4725.0_real32)
  end do

  if (any(.not. ieee_is_finite(o0)) .or. any(.not. ieee_is_finite(o1)) .or. &
      any(.not. ieee_is_finite(o2))) then
    write(*,'(A)') "FAIL: non-finite device result"
    error stop 1
  end if

  err0 = sqrt(sum((ref - o0) * (ref - o0)))
  err1 = sqrt(sum((ref - o1) * (ref - o1)))
  err2 = sqrt(sum((ref - o2) * (ref - o2)))

  write(*,*)
  write(*,'(A)') "Error statistics for the kernels:"
  if (max(err0, err1, err2) < 10.0_real32) then
    write(*,'(F8.6,1X,F8.6,1X,F8.6,A)') err0, err1, err2, " "
  else
    write(*,'(F0.6,1X,F0.6,1X,F0.6,A)') err0, err1, err2, " "
  end if

contains

  subroutine parse_args(n, repeat)
    integer, intent(out) :: n, repeat
    character(len=64) :: arg
    integer :: status

    if (command_argument_count() /= 2) then
      call get_command_argument(0, arg)
      write(*,'(A,A,A)') "Usage ", trim(arg), " <n> <repeat>"
      stop 1
    end if

    call get_command_argument(1, arg)
    read(arg, *, iostat=status) n
    if (status /= 0 .or. n <= 0) error stop "invalid n"

    call get_command_argument(2, arg)
    read(arg, *, iostat=status) repeat
    if (status /= 0 .or. repeat <= 0) error stop "invalid repeat"
  end subroutine parse_args

  subroutine k0_device(a, o, n)
    integer, intent(in) :: n
    real(real32), intent(in) :: a(n)
    real(real32), intent(out) :: o(n)
    integer :: t
    real(real32) :: x

    !$omp target teams distribute parallel do thread_limit(256) private(x)
    do t = 1, n
      x = a(t)
      o(t) = cosh(x) / sinh(x) - 1.0_real32 / x
    end do
    !$omp end target teams distribute parallel do
  end subroutine k0_device

  subroutine k1_device(a, o, n)
    integer, intent(in) :: n
    real(real32), intent(in) :: a(n)
    real(real32), intent(out) :: o(n)
    integer :: t
    real(real32) :: x

    !$omp target teams distribute parallel do thread_limit(256) private(x)
    do t = 1, n
      x = a(t)
      o(t) = 1.0_real32 / tanh(x) - 1.0_real32 / x
    end do
    !$omp end target teams distribute parallel do
  end subroutine k1_device

  subroutine k2_device(a, o, n)
    integer, intent(in) :: n
    real(real32), intent(in) :: a(n)
    real(real32), intent(out) :: o(n)
    integer :: t
    real(real32) :: x, s, r

    !$omp target teams distribute parallel do thread_limit(256) private(x, s, r)
    do t = 1, n
      x = a(t)
      s = x * x
      r = 7.70960469e-8_real32
      r = r * s - 1.65101926e-6_real32
      r = r * s + 2.03457112e-5_real32
      r = r * s - 2.10521728e-4_real32
      r = r * s + 2.11580913e-3_real32
      r = r * s - 2.22220998e-2_real32
      r = r * s + 8.33333284e-2_real32
      r = r * x + 0.25_real32 * x
      o(t) = r
    end do
    !$omp end target teams distribute parallel do
  end subroutine k2_device

end program langevin_main
