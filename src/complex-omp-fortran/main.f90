program main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2
  integer :: n, repeat, iter
  integer, allocatable :: checksum(:)
  real(real64) :: start_time, end_time
  logical :: complex_float_check, complex_double_check
  integer(int64), parameter :: lcg_a = 2806196910506780709_int64
  integer(int64), parameter :: lcg_c = 1_int64
  integer(int64), parameter :: lcg_mask = huge(0_int64)
  real(real64), parameter :: lcg_scale = 1.0_real64 / 9223372036854775808.0_real64

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *) n
  read(arg2, *) repeat
  if (n <= 0 .or. repeat <= 0) stop 1

  allocate(checksum(n))

  !$omp target data map(alloc: checksum(1:n))
  call complex_float(checksum, n)
  call complex_double(checksum, n)

  start_time = omp_get_wtime()
  do iter = 1, repeat
    call complex_float(checksum, n)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average kernel execution time (float) ', &
    (end_time - start_time) / real(repeat, real64), ' (s)'

  !$omp target update from(checksum(1:n))
  complex_float_check = all(checksum == 5)

  start_time = omp_get_wtime()
  do iter = 1, repeat
    call complex_double(checksum, n)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average kernel execution time (double) ', &
    (end_time - start_time) / real(repeat, real64), ' (s)'

  !$omp target update from(checksum(1:n))
  complex_double_check = all(checksum == 5)
  !$omp end target data

  if (complex_float_check .and. complex_double_check) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(checksum)

contains

  subroutine complex_float(checksum, n)
    integer, intent(out) :: checksum(n)
    integer, intent(in) :: n
    integer :: i, s
    real(real32) :: r1, r2, r3, r4
    complex(real32) :: z1, z2

    !$omp target teams distribute parallel do thread_limit(256)
    do i = 1, n
      call random_float4(i - 1, r1, r2, r3, r4)
      z1 = cmplx(r1, r2, kind=real32)
      z2 = cmplx(r3, r4, kind=real32)

      s = 0
      if (abs(abs(z1 * z2) - abs(z1) * abs(z2)) < 1.0e-3_real32) s = s + 1
      if (abs(abs(z1 + z2) * abs(z1 + z2) - real((z1 + z2) * conjg(z1 + z2), real32)) &
          < 1.0e-3_real32) s = s + 1
      if (abs(abs(z1 - z2) * abs(z1 - z2) - real((z1 - z2) * conjg(z1 - z2), real32)) &
          < 1.0e-3_real32) s = s + 1
      if (abs(real(z1 * conjg(z2) + z2 * conjg(z1), real32) - &
          2.0_real32 * (real(z1, real32) * real(z2, real32) + aimag(z1) * aimag(z2))) &
          < 1.0e-3_real32) s = s + 1
      if (abs(abs(conjg(z1) / z2) - abs(conjg(z1) / conjg(z2))) < 1.0e-3_real32) s = s + 1

      checksum(i) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine complex_float

  subroutine complex_double(checksum, n)
    integer, intent(out) :: checksum(n)
    integer, intent(in) :: n
    integer :: i, s
    real(real64) :: r1, r2, r3, r4
    complex(real64) :: z1, z2

    !$omp target teams distribute parallel do thread_limit(256)
    do i = 1, n
      call random_double4(i - 1, r1, r2, r3, r4)
      z1 = cmplx(r1, r2, kind=real64)
      z2 = cmplx(r3, r4, kind=real64)

      s = 0
      if (abs(abs(z1 * z2) - abs(z1) * abs(z2)) < 1.0e-3_real64) s = s + 1
      if (abs(abs(z1 + z2) * abs(z1 + z2) - real((z1 + z2) * conjg(z1 + z2), real64)) &
          < 1.0e-3_real64) s = s + 1
      if (abs(abs(z1 - z2) * abs(z1 - z2) - real((z1 - z2) * conjg(z1 - z2), real64)) &
          < 1.0e-3_real64) s = s + 1
      if (abs(real(z1 * conjg(z2) + z2 * conjg(z1), real64) - &
          2.0_real64 * (real(z1, real64) * real(z2, real64) + aimag(z1) * aimag(z2))) &
          < 1.0e-3_real64) s = s + 1
      if (abs(abs(conjg(z1) / z2) - abs(conjg(z1) / conjg(z2))) < 1.0e-3_real64) s = s + 1

      checksum(i) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine complex_double

  subroutine random_float4(i0, r1, r2, r3, r4)
    integer, intent(in) :: i0
    real(real32), intent(out) :: r1, r2, r3, r4
    integer(int64) :: seed

    seed = fast_forward_lcg(1_int64, int(i0, int64))
    r1 = real(lcg_random_double(seed), real32)
    r2 = real(lcg_random_double(seed), real32)
    r3 = real(lcg_random_double(seed), real32)
    r4 = real(lcg_random_double(seed), real32)
  end subroutine random_float4

  subroutine random_double4(i0, r1, r2, r3, r4)
    integer, intent(in) :: i0
    real(real64), intent(out) :: r1, r2, r3, r4
    integer(int64) :: seed

    seed = fast_forward_lcg(1_int64, int(i0, int64))
    r1 = lcg_random_double(seed)
    r2 = lcg_random_double(seed)
    r3 = lcg_random_double(seed)
    r4 = lcg_random_double(seed)
  end subroutine random_double4

  real(real64) function lcg_random_double(seed)
    integer(int64), intent(inout) :: seed

    seed = iand(lcg_a * seed + lcg_c, lcg_mask)
    lcg_random_double = real(seed, real64) * lcg_scale
  end function lcg_random_double

  integer(int64) function fast_forward_lcg(seed, n)
    integer(int64), intent(in) :: seed, n
    integer(int64) :: a, c, a_new, c_new, n_work

    a = lcg_a
    c = lcg_c
    a_new = 1_int64
    c_new = 0_int64
    n_work = iand(n, lcg_mask)

    do while (n_work > 0_int64)
      if (iand(n_work, 1_int64) /= 0_int64) then
        a_new = iand(a_new * a, lcg_mask)
        c_new = iand(c_new * a + c, lcg_mask)
      end if

      c = iand(c * (a + 1_int64), lcg_mask)
      a = iand(a * a, lcg_mask)
      n_work = shiftr(n_work, 1)
    end do

    fast_forward_lcg = iand(a_new * seed + c_new, lcg_mask)
  end function fast_forward_lcg

end program main
