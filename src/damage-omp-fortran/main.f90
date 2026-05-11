program main
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  character(len=256) :: arg0, arg
  integer :: n, repeat, m, i, j, s
  integer(int64) :: seed
  integer, allocatable :: nlist(:), family(:), n_neigh(:)
  real(real64), allocatable :: damage(:)
  real(real64) :: start_time, elapsed

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <number of points> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) n
  call get_command_argument(2, arg); read(arg, *) repeat
  if (n <= 0 .or. repeat <= 0) stop 1

  m = (n + block_size - 1) / block_size
  allocate(nlist(n), family(m), n_neigh(m), damage(m))

  seed = 123_int64
  do i = 1, n
    if (lcg_random_double(seed) > 0.5_real64) then
      nlist(i) = 1
    else
      nlist(i) = -1
    end if
  end do

  do i = 1, m
    s = 0
    do j = (i - 1) * block_size + 1, min(i * block_size, n)
      if (nlist(j) /= -1) s = s + 1
    end do
    family(i) = int(real(s + 1, real64) + real(s, real64) * lcg_random_double(seed))
  end do

  !$omp target data map(to: nlist(1:n), family(1:m)) map(from: n_neigh(1:m), damage(1:m))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call damage_kernel(n, m, nlist, family, n_neigh, damage)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data
  write(*,'(A,F0.6,A)') 'Average kernel execution time ', elapsed / real(repeat, real64), ' (s)'

  call validate(n, m, nlist, family, n_neigh, damage)

  !$omp target data map(to: nlist(1:n), family(1:m)) map(from: n_neigh(1:m), damage(1:m))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call damage_kernel(n, m, nlist, family, n_neigh, damage)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data
  write(*,'(A,F0.6,A)') 'Average kernel execution time ', elapsed / real(repeat, real64), ' (s)'

  call validate(n, m, nlist, family, n_neigh, damage)

  deallocate(nlist, family, n_neigh, damage)

contains

  real(real64) function lcg_random_double(seed)
    integer(int64), intent(inout) :: seed
    integer(int64), parameter :: a = 1103515245_int64
    integer(int64), parameter :: c = 12345_int64
    integer(int64), parameter :: modulus = 2147483647_int64

    seed = modulo(a * seed + c, modulus)
    lcg_random_double = real(seed, real64) / real(modulus, real64)
  end function lcg_random_double

  subroutine damage_kernel(n, m, nlist, family, n_neigh, damage)
    integer, intent(in) :: n, m
    integer, intent(in) :: nlist(n), family(m)
    integer, intent(out) :: n_neigh(m)
    real(real64), intent(out) :: damage(m)
    integer :: gid, j, lower, upper, sum

    !$omp target teams distribute parallel do private(lower, upper, sum, j) thread_limit(block_size)
    do gid = 1, m
      lower = (gid - 1) * block_size + 1
      upper = min(gid * block_size, n)
      sum = 0
      do j = lower, upper
        if (nlist(j) /= -1) sum = sum + 1
      end do
      n_neigh(gid) = sum
      damage(gid) = 1.0_real64 - real(sum, real64) / real(family(gid), real64)
    end do
    !$omp end target teams distribute parallel do
  end subroutine damage_kernel

  subroutine validate(n, m, nlist, family, n_neigh, damage)
    integer, intent(in) :: n, m
    integer, intent(in) :: nlist(n), family(m), n_neigh(m)
    real(real64), intent(in) :: damage(m)
    integer :: gid, j, lower, upper, sum
    real(real64) :: damage_ref
    logical :: ok

    ok = .true.
    do gid = 1, m
      lower = (gid - 1) * block_size + 1
      upper = min(gid * block_size, n)
      sum = 0
      do j = lower, upper
        if (nlist(j) /= -1) sum = sum + 1
      end do
      damage_ref = 1.0_real64 - real(sum, real64) / real(family(gid), real64)
      if (n_neigh(gid) /= sum .or. abs(damage(gid) - damage_ref) > 1.0e-6_real64) then
        ok = .false.
        exit
      end if
    end do

    if (ok) then
      write(*,'(A)') 'PASS'
    else
      write(*,'(A)') 'FAIL'
    end if
  end subroutine validate

end program main
