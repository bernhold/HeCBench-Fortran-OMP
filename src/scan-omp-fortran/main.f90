program main
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0
  integer :: n, repeat

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage: ', trim(arg0), ' <number of elements> <repeat>'
    stop 1
  end if

  n = read_arg(1)
  repeat = read_arg(2)
  if (n <= 0 .or. repeat <= 0) stop 1

  call run_block_size(128, n, repeat)
  call run_block_size(256, n, repeat)
  call run_block_size(512, n, repeat)
  call run_block_size(1024, n, repeat)
  call run_block_size(2048, n, repeat)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer

    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine run_block_size(block_elems, n, repeat)
    integer, intent(in) :: block_elems, n, repeat
    integer :: phase
    logical :: timing

    do phase = 0, 1
      timing = phase > 0
      print '(A)'
      print '(A,I0)', 'The number of elements to scan in a thread block: ', block_elems
      call run_test(block_elems, 1, n, repeat, timing)
      call run_test(block_elems, 2, n, repeat, timing)
      call run_test(block_elems, 4, n, repeat, timing)
      call run_test(block_elems, 8, n, repeat, timing)
    end do
  end subroutine run_block_size

  subroutine run_test(block_elems, elem_size, n, repeat, timing)
    integer, intent(in) :: block_elems, elem_size, n, repeat
    logical, intent(in) :: timing
    integer :: num_blocks, nelems, i
    integer(int64), allocatable :: input(:), output(:), ref_output(:)
    real(real64) :: start_time, end_time, base_us, bcao_us, reduction

    num_blocks = (n + block_elems - 1) / block_elems
    nelems = num_blocks * block_elems
    allocate(input(nelems), output(nelems), ref_output(nelems))

    do i = 1, nelems
      input(i) = int(mod(17 * (i - 1) + 3, 5) + 1, int64)
    end do
    call scan_reference(input, ref_output, num_blocks, block_elems)
    output = 0_int64

    !$omp target data map(to: input(1:nelems)) map(alloc: output(1:nelems))
    start_time = omp_get_wtime()
    do i = 1, repeat
      call scan_device(input, output, num_blocks, block_elems)
    end do
    end_time = omp_get_wtime()
    base_us = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)

    !$omp target update from(output(1:nelems))
    if (timing) then
      print '(A,I0,A,F0.6,A)', 'Element size in bytes is ', elem_size, &
        '. Average execution time of scan (w/  bank conflicts): ', base_us, ' (us)'
    else
      call verify(ref_output, output, nelems)
    end if

    start_time = omp_get_wtime()
    do i = 1, repeat
      call scan_device(input, output, num_blocks, block_elems)
    end do
    end_time = omp_get_wtime()
    bcao_us = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)

    !$omp target update from(output(1:nelems))
    if (timing) then
      if (base_us /= 0.0_real64) then
        reduction = (base_us - bcao_us) / base_us * 100.0_real64
      else
        reduction = 0.0_real64
      end if
      write(*, '(A,I0,A,F0.6,A,F0.1,A)') 'Element size in bytes is ', elem_size, &
        '. Average execution time of scan (w/o bank conflicts): ', bcao_us, ' (us). Reduce the time by ', reduction, '%'
    else
      call verify(ref_output, output, nelems)
    end if
    !$omp end target data

    deallocate(input, output, ref_output)
  end subroutine run_test

  subroutine scan_reference(input, output, num_blocks, block_elems)
    integer(int64), intent(in) :: input(:)
    integer(int64), intent(out) :: output(:)
    integer, intent(in) :: num_blocks, block_elems
    integer :: bid, j, base

    do bid = 0, num_blocks - 1
      base = bid * block_elems
      output(base + 1) = 0_int64
      do j = 2, block_elems
        output(base + j) = output(base + j - 1) + input(base + j - 1)
      end do
    end do
  end subroutine scan_reference

  subroutine scan_device(input, output, num_blocks, block_elems)
    integer(int64), intent(in) :: input(:)
    integer(int64), intent(out) :: output(:)
    integer, intent(in) :: num_blocks, block_elems
    integer :: bid, j, base

    !$omp target teams distribute parallel do thread_limit(256) private(j, base)
    do bid = 0, num_blocks - 1
      base = bid * block_elems
      output(base + 1) = 0_int64
      do j = 2, block_elems
        output(base + j) = output(base + j - 1) + input(base + j - 1)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine scan_device

  subroutine verify(ref_output, output, nelems)
    integer(int64), intent(in) :: ref_output(:), output(:)
    integer, intent(in) :: nelems
    integer :: i
    logical :: ok

    ok = .true.
    do i = 1, nelems
      if (ref_output(i) /= output(i)) then
        ok = .false.
        exit
      end if
    end do

    if (ok) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
      stop 1
    end if
  end subroutine verify

end program main
