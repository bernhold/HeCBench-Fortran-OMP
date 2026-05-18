program rodrigues_main
  use, intrinsic :: iso_fortran_env, only: real32
  use, intrinsic :: iso_c_binding, only: c_int
  use omp_lib
  implicit none

  integer :: n, repeat, i
  real(real32), allocatable :: x(:), y(:), z(:), x_ref(:), y_ref(:), z_ref(:)
  real(real32), allocatable :: x4(:), y4(:), z4(:), w4(:), x4_ref(:), y4_ref(:), z4_ref(:), w4_ref(:)
  real(real32) :: wx, wy, wz, norm, angle, max_error
  real(real32), parameter :: tolerance = 5.0e-4_real32
  real(real32) :: a, b, c, d
  real(8) :: start_time, end_time

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() result(value) bind(C, name="rand")
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  call parse_args(n, repeat)

  wx = -0.3_real32
  wy = -0.6_real32
  wz = 0.15_real32
  norm = 1.0_real32 / sqrt(wx * wx + wy * wy + wz * wz)
  wx = wx * norm
  wy = wy * norm
  wz = wz * norm
  angle = 0.5_real32

  allocate(x(n), y(n), z(n), x_ref(n), y_ref(n), z_ref(n))
  allocate(x4(n), y4(n), z4(n), w4(n), x4_ref(n), y4_ref(n), z4_ref(n), w4_ref(n))

  call c_srand(123_c_int)
  do i = 1, n
    a = real(c_rand(), real32)
    b = real(c_rand(), real32)
    c = real(c_rand(), real32)
    d = sqrt(a * a + b * b + c * c)
    x(i) = a / d
    y(i) = b / d
    z(i) = c / d
    x4(i) = x(i)
    y4(i) = y(i)
    z4(i) = z(i)
    w4(i) = 0.0_real32
  end do

  x_ref = x
  y_ref = y
  z_ref = z
  x4_ref = x4
  y4_ref = y4
  z4_ref = z4
  w4_ref = w4

  !$omp target data map(tofrom: x, y, z)
  start_time = omp_get_wtime()
  do i = 1, repeat
    call rotate3_device(n, angle, wx, wy, wz, x, y, z)
  end do
  end_time = omp_get_wtime()
  !$omp end target data
  write(*,'(A,F0.6,A)') "Average kernel execution time (float3): ", &
       real((end_time - start_time) * 1.0d6 / dble(repeat), real32), " (us)"

  !$omp target data map(tofrom: x4, y4, z4, w4)
  start_time = omp_get_wtime()
  do i = 1, repeat
    call rotate4_device(n, angle, wx, wy, wz, x4, y4, z4, w4)
  end do
  end_time = omp_get_wtime()
  !$omp end target data
  write(*,'(A,F0.6,A)') "Average kernel execution time (float4): ", &
       real((end_time - start_time) * 1.0d6 / dble(repeat), real32), " (us)"

  do i = 1, repeat
    call rotate3_cpu(n, angle, wx, wy, wz, x_ref, y_ref, z_ref)
    call rotate4_cpu(n, angle, wx, wy, wz, x4_ref, y4_ref, z4_ref, w4_ref)
  end do

  max_error = maxval(abs(x - x_ref))
  max_error = max(max_error, maxval(abs(y - y_ref)))
  max_error = max(max_error, maxval(abs(z - z_ref)))
  max_error = max(max_error, maxval(abs(x4 - x4_ref)))
  max_error = max(max_error, maxval(abs(y4 - y4_ref)))
  max_error = max(max_error, maxval(abs(z4 - z4_ref)))
  max_error = max(max_error, maxval(abs(w4 - w4_ref)))
  if (max_error > tolerance) then
    write(*,'(A,ES12.4)') "FAIL: max rotation error = ", max_error
    error stop 1
  end if

contains

  subroutine parse_args(n, repeat)
    integer, intent(out) :: n, repeat
    character(len=64) :: arg
    integer :: argc, status

    argc = command_argument_count()
    if (argc /= 2) then
      call get_command_argument(0, arg)
      write(*,'(A,A,A)') "Usage: ", trim(arg), " <number of points> <repeat>"
      stop 1
    end if

    call get_command_argument(1, arg)
    read(arg, *, iostat=status) n
    if (status /= 0 .or. n <= 0) error stop "invalid number of points"

    call get_command_argument(2, arg)
    read(arg, *, iostat=status) repeat
    if (status /= 0 .or. repeat <= 0) error stop "invalid repeat"
  end subroutine parse_args

  subroutine rotate3_device(n, angle, wx, wy, wz, x, y, z)
    integer, intent(in) :: n
    real(real32), intent(in) :: angle, wx, wy, wz
    real(real32), intent(inout) :: x(n), y(n), z(n)
    integer :: i
    real(real32) :: s, c, mc, px, py, pz
    real(real32) :: m1, m2, m3, m4, m5, m6, m7, m8, m9

    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(s, c, mc, px, py, pz, m1, m2, m3, m4, m5, m6, m7, m8, m9)
    do i = 1, n
      s = sin(angle)
      c = cos(angle)
      mc = 1.0_real32 - c
      px = x(i)
      py = y(i)
      pz = z(i)
      m1 = c + wx * wx * mc
      m2 = wz * s + wx * wy * mc
      m3 = -wy * s + wx * wz * mc
      m4 = -wz * s + wx * wy * mc
      m5 = c + wy * wy * mc
      m6 = wx * s + wy * wz * mc
      m7 = wy * s + wx * wz * mc
      m8 = -wx * s + wy * wz * mc
      m9 = c + wz * wz * mc
      x(i) = px * m1 + py * m2 + pz * m3
      y(i) = px * m4 + py * m5 + pz * m6
      z(i) = px * m7 + py * m8 + pz * m9
    end do
    !$omp end target teams distribute parallel do
  end subroutine rotate3_device

  subroutine rotate4_device(n, angle, wx, wy, wz, x, y, z, w)
    integer, intent(in) :: n
    real(real32), intent(in) :: angle, wx, wy, wz
    real(real32), intent(inout) :: x(n), y(n), z(n), w(n)
    integer :: i
    real(real32) :: s, c, mc, px, py, pz
    real(real32) :: m1, m2, m3, m4, m5, m6, m7, m8, m9

    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(s, c, mc, px, py, pz, m1, m2, m3, m4, m5, m6, m7, m8, m9)
    do i = 1, n
      s = sin(angle)
      c = cos(angle)
      mc = 1.0_real32 - c
      px = x(i)
      py = y(i)
      pz = z(i)
      m1 = c + wx * wx * mc
      m2 = wz * s + wx * wy * mc
      m3 = -wy * s + wx * wz * mc
      m4 = -wz * s + wx * wy * mc
      m5 = c + wy * wy * mc
      m6 = wx * s + wy * wz * mc
      m7 = wy * s + wx * wz * mc
      m8 = -wx * s + wy * wz * mc
      m9 = c + wz * wz * mc
      x(i) = px * m1 + py * m2 + pz * m3
      y(i) = px * m4 + py * m5 + pz * m6
      z(i) = px * m7 + py * m8 + pz * m9
      w(i) = 0.0_real32
    end do
    !$omp end target teams distribute parallel do
  end subroutine rotate4_device

  subroutine rotate3_cpu(n, angle, wx, wy, wz, x, y, z)
    integer, intent(in) :: n
    real(real32), intent(in) :: angle, wx, wy, wz
    real(real32), intent(inout) :: x(n), y(n), z(n)
    integer :: i

    do i = 1, n
      call rotate_one(angle, wx, wy, wz, x(i), y(i), z(i))
    end do
  end subroutine rotate3_cpu

  subroutine rotate4_cpu(n, angle, wx, wy, wz, x, y, z, w)
    integer, intent(in) :: n
    real(real32), intent(in) :: angle, wx, wy, wz
    real(real32), intent(inout) :: x(n), y(n), z(n), w(n)
    integer :: i

    do i = 1, n
      call rotate_one(angle, wx, wy, wz, x(i), y(i), z(i))
      w(i) = 0.0_real32
    end do
  end subroutine rotate4_cpu

  subroutine rotate_one(angle, wx, wy, wz, x, y, z)
    real(real32), intent(in) :: angle, wx, wy, wz
    real(real32), intent(inout) :: x, y, z
    real(real32) :: s, c, mc, px, py, pz
    real(real32) :: m1, m2, m3, m4, m5, m6, m7, m8, m9

    s = sin(angle)
    c = cos(angle)
    mc = 1.0_real32 - c
    px = x
    py = y
    pz = z
    m1 = c + wx * wx * mc
    m2 = wz * s + wx * wy * mc
    m3 = -wy * s + wx * wz * mc
    m4 = -wz * s + wx * wy * mc
    m5 = c + wy * wy * mc
    m6 = wx * s + wy * wz * mc
    m7 = wy * s + wx * wz * mc
    m8 = -wx * s + wy * wz * mc
    m9 = c + wz * wz * mc
    x = px * m1 + py * m2 + pz * m3
    y = px * m4 + py * m5 + pz * m6
    z = px * m7 + py * m8 + pz * m9
  end subroutine rotate_one

end program rodrigues_main
