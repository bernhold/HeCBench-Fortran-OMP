program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2
  integer :: n, repeat, iter
  integer, allocatable :: checksum(:)
  real(real64) :: start_time, end_time
  logical :: complex_float_check, complex_double_check

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
      r1 = random_float(i, 1)
      r2 = random_float(i, 2)
      r3 = random_float(i, 3)
      r4 = random_float(i, 4)
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
      r1 = random_double(i, 1)
      r2 = random_double(i, 2)
      r3 = random_double(i, 3)
      r4 = random_double(i, 4)
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

  real(real32) function random_float(i, lane)
    integer, intent(in) :: i, lane
    integer :: value

    value = mod(i * (lane * 37 + 11) + lane * 97, 1000)
    random_float = 0.125_real32 + real(value, real32) / 2048.0_real32
  end function random_float

  real(real64) function random_double(i, lane)
    integer, intent(in) :: i, lane
    integer :: value

    value = mod(i * (lane * 37 + 11) + lane * 97, 1000)
    random_double = 0.125_real64 + real(value, real64) / 2048.0_real64
  end function random_double

end program main
