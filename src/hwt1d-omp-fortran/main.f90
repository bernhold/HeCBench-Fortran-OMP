program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  real(real32), parameter :: tolerance = 0.1_real32
  real(real32), parameter :: rsqrt_two = 0.7071_real32
  integer, parameter :: thread_count = 256

  integer(int64) :: signal_length_arg
  integer :: signal_length, iterations, levels, i
  real(real32), allocatable :: in_data(:), device_out(:), work(:), host_out(:)
  real(real64) :: start_time, end_time
  logical :: ok

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name='rand') result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  if (command_argument_count() /= 2) then
    write(*,'(A,A)') 'Usage: ', './main <signal length> <repeat>'
    stop 1
  end if

  signal_length_arg = read_int64_arg(1)
  iterations = read_int_arg(2)
  signal_length = round_to_power_of_2(signal_length_arg)

  if (get_levels(signal_length, levels) /= 0) then
    write(*,'(A)') 'signalLength > 2 ^ 23 not supported'
    stop 1
  end if

  allocate(in_data(signal_length), device_out(signal_length), work(signal_length), host_out(signal_length))

  call initialize_input(in_data)
  device_out = 0.0_real32
  work = 0.0_real32
  host_out = 0.0_real32

  write(*,'(A,I0,A)') 'Executing kernel for ', iterations, ' iterations'
  write(*,'(A)') '-------------------------------------------'

  !$omp target data map(to: in_data(1:signal_length)) map(tofrom: work(1:signal_length), device_out(1:signal_length))
  start_time = omp_get_wtime()
  do i = 1, iterations
    call hwt_device(in_data, work, device_out, signal_length)
  end do
  end_time = omp_get_wtime()
  !$omp end target data

  write(*,'(A,ES12.5,A)') 'Average device offload time ', real((end_time - start_time) / real(iterations, real64), real32), ' (s)'

  call hwt_host(in_data, host_out, signal_length)
  ok = compare_outputs(device_out, host_out)

  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(in_data, device_out, work, host_out)

contains

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  integer(int64) function read_int64_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int64_arg

  integer function round_to_power_of_2(value) result(rounded)
    integer(int64), intent(in) :: value
    integer(int64) :: current
    if (value <= 1_int64) then
      rounded = 1
      return
    end if
    current = 1_int64
    do while (current < value)
      current = current * 2_int64
    end do
    rounded = int(current)
  end function round_to_power_of_2

  integer function get_levels(length, levels) result(status)
    integer, intent(in) :: length
    integer, intent(out) :: levels
    integer :: idx
    status = 1
    levels = 0
    do idx = 0, 23
      if (length == 2 ** idx) then
        levels = idx
        status = 0
        exit
      end if
    end do
  end function get_levels

  subroutine initialize_input(data)
    real(real32), intent(out) :: data(:)
    integer :: idx
    call c_srand(2_c_int)
    do idx = 1, size(data)
      data(idx) = real(mod(c_rand(), 10_c_int), real32)
    end do
  end subroutine initialize_input

  subroutine hwt_device(in_data, work, out_data, signal_length)
    real(real32), intent(in) :: in_data(:)
    real(real32), intent(inout) :: work(:), out_data(:)
    integer, intent(in) :: signal_length
    integer :: idx, length
    real(real32) :: norm, data0, data1

    norm = 1.0_real32 / sqrt(real(signal_length, real32))
    !$omp target teams distribute parallel do thread_limit(thread_count)
    do idx = 1, signal_length
      work(idx) = in_data(idx) * norm
      out_data(idx) = 0.0_real32
    end do
    !$omp end target teams distribute parallel do

    length = signal_length
    do while (length > 1)
      !$omp target teams distribute parallel do thread_limit(thread_count) private(data0, data1)
      do idx = 1, length / 2
        data0 = work(2 * idx - 1)
        data1 = work(2 * idx)
        out_data(idx) = (data0 + data1) * rsqrt_two
        out_data(length / 2 + idx) = (data0 - data1) * rsqrt_two
      end do
      !$omp end target teams distribute parallel do

      !$omp target teams distribute parallel do thread_limit(thread_count)
      do idx = 1, length
        work(idx) = out_data(idx)
      end do
      !$omp end target teams distribute parallel do
      length = length / 2
    end do
  end subroutine hwt_device

  subroutine hwt_host(in_data, out_data, signal_length)
    real(real32), intent(in) :: in_data(:)
    real(real32), intent(out) :: out_data(:)
    integer, intent(in) :: signal_length
    real(real32), allocatable :: temp(:)
    integer :: idx, length
    real(real32) :: data0, data1

    allocate(temp(signal_length))
    temp = in_data / sqrt(real(signal_length, real32))
    out_data = 0.0_real32

    length = signal_length
    do while (length > 1)
      do idx = 1, length / 2
        data0 = temp(2 * idx - 1)
        data1 = temp(2 * idx)
        out_data(idx) = (data0 + data1) * rsqrt_two
        out_data(length / 2 + idx) = (data0 - data1) * rsqrt_two
      end do
      temp(1:length) = out_data(1:length)
      length = length / 2
    end do

    deallocate(temp)
  end subroutine hwt_host

  logical function compare_outputs(device_out, host_out) result(ok)
    real(real32), intent(in) :: device_out(:), host_out(:)
    integer :: idx
    ok = .true.
    do idx = 1, size(device_out)
      if (abs(device_out(idx) - host_out(idx)) > tolerance) then
        ok = .false.
        exit
      end if
    end do
  end function compare_outputs

end program main
