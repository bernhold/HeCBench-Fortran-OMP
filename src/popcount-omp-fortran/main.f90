program popcount_main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name="rand") result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  integer, parameter :: block_size = 256
  integer(int64), parameter :: m1 = int(z'5555555555555555', int64)
  integer(int64), parameter :: m2 = int(z'3333333333333333', int64)
  integer(int64), parameter :: m4 = int(z'0f0f0f0f0f0f0f0f', int64)
  integer :: length, repeat, mode
  integer(int64), allocatable :: data(:)
  integer(int32), allocatable :: result(:)
  integer :: i
  real(real64) :: start_time, elapsed_us
  character(len=64) :: arg

  if (command_argument_count() /= 2) then
    write(*,'("Usage: ./main <length> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) length
  call get_command_argument(2, arg)
  read(arg, *) repeat

  allocate(data(0:length-1), result(0:length-1))

  call c_srand(2_c_int)
  do i = 0, length - 1
    data(i) = ior(ishft(int(c_rand(), int64), 32), int(c_rand(), int64))
  end do

  !$omp target data map(to: data) map(alloc: result)
  do mode = 1, 6
    start_time = omp_get_wtime()
    do i = 1, repeat
      call popcount_kernel(data, result, length, mode)
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
    !$omp target update from(result)
    write(*,'("Average kernel execution time (pc",I0,"): ",F0.6," (us)")') mode, elapsed_us
    call check_results(data, result, length)
  end do
  !$omp end target data

contains

  subroutine popcount_kernel(data, result, length, mode)
    integer(int64), intent(in) :: data(0:)
    integer(int32), intent(out) :: result(0:)
    integer, intent(in) :: length, mode
    integer :: i, bit_index
    integer(int64) :: x
    integer(int32) :: count_value

    !$omp target teams distribute parallel do thread_limit(block_size) private(x, bit_index, count_value)
    do i = 0, length - 1
      x = data(i)
      select case (mode)
      case (1)
        x = x - iand(shiftr(x, 1), m1)
        x = iand(x, m2) + iand(shiftr(x, 2), m2)
        x = iand(x + shiftr(x, 4), m4)
        x = x + shiftr(x, 8)
        x = x + shiftr(x, 16)
        x = x + shiftr(x, 32)
        count_value = int(iand(x, 127_int64), int32)
      case (2)
        count_value = int(popcnt(x), int32)
      case (3)
        count_value = 0_int32
        do while (x /= 0_int64)
          count_value = count_value + 1_int32
          x = iand(x, x - 1_int64)
        end do
      case (4)
        count_value = 0_int32
        do bit_index = 0, 63
          count_value = count_value + int(iand(x, 1_int64), int32)
          x = shiftr(x, 1)
        end do
      case (5)
        count_value = 0_int32
        do bit_index = 0, 56, 8
          count_value = count_value + int(popcnt(iand(shiftr(x, bit_index), 255_int64)), int32)
        end do
      case default
        count_value = int(popcnt(x), int32)
      end select
      result(i) = count_value
    end do
    !$omp end target teams distribute parallel do
  end subroutine popcount_kernel

  subroutine check_results(data, result, length)
    integer(int64), intent(in) :: data(0:)
    integer(int32), intent(in) :: result(0:)
    integer, intent(in) :: length
    integer :: i
    logical :: error

    error = .false.
    do i = 0, length - 1
      if (popcnt(data(i)) /= result(i)) then
        error = .true.
        exit
      end if
    end do

    if (error) then
      write(*,'("Fail")')
    else
      write(*,'("Success")')
    end if
  end subroutine check_results

end program popcount_main
