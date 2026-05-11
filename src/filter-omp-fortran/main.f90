program main
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2, arg3
  integer :: num_elems, block_size, repeat, i
  integer :: nres(1)
  integer, allocatable :: input(:), output(:)
  integer :: expected_count, got_count
  integer(int64) :: expected_sum, got_sum
  real(real64) :: start_time, end_time, avg_ms
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    print '(3A)', 'Usage: ', trim(arg0), ' <number of elements> <block size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) num_elems
  read(arg2, *) block_size
  read(arg3, *) repeat

  if (num_elems <= 0 .or. block_size <= 0 .or. repeat <= 0) stop 1

  allocate(input(num_elems), output(num_elems))
  do i = 1, num_elems
    input(i) = (i - 1) - num_elems / 2
  end do
  output = 0
  nres = 0

  !$omp target data map(to: input(1:num_elems)) map(tofrom: nres(1:1)) map(from: output(1:num_elems))
  start_time = omp_get_wtime()
  do i = 1, repeat
    nres(1) = 0
    !$omp target update to(nres(1:1))
    call filter_positive(input, output, nres, num_elems, block_size)
  end do
  end_time = omp_get_wtime()
  avg_ms = (end_time - start_time) * 1.0e3_real64 / real(repeat, real64)
  print '(A,F0.6,A)', 'Average kernel execution time ', avg_ms, ' (ms)'
  !$omp end target data

  expected_count = 0
  expected_sum = 0_int64
  do i = 1, num_elems
    if (input(i) > 0) then
      expected_count = expected_count + 1
      expected_sum = expected_sum + int(input(i), int64)
    end if
  end do

  got_count = nres(1)
  got_sum = 0_int64
  do i = 1, got_count
    got_sum = got_sum + int(output(i), int64)
  end do

  ok = (got_count == expected_count) .and. (got_sum == expected_sum)

  print '(A)'
  if (ok) then
    print '(A)', 'Filter using shared memory PASS '
  else
    print '(A)', 'Filter using shared memory FAIL '
    stop 1
  end if

  deallocate(input, output)

contains

  subroutine filter_positive(input, output, nres, num_elems, block_size)
    integer, intent(in) :: input(:), num_elems, block_size
    integer, intent(out) :: output(:)
    integer, intent(inout) :: nres(:)
    integer :: idx, pos, value, teams

    teams = (num_elems + block_size - 1) / block_size
    !$omp target teams distribute parallel do num_teams(teams) thread_limit(block_size) private(value, pos)
    do idx = 1, num_elems
      value = input(idx)
      if (value > 0) then
        !$omp atomic capture
        pos = nres(1)
        nres(1) = nres(1) + 1
        !$omp end atomic
        output(pos + 1) = value
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine filter_positive

end program main
