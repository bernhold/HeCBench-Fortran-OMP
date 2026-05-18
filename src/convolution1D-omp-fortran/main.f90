program main
  use, intrinsic :: iso_fortran_env, only : int16, int32, int64, real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer, parameter :: max_mask_width = 10
  integer, parameter :: max_block_size = 1024

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

  integer :: input_width, repeat, mask_width

  if (command_argument_count() /= 2) then
    print '(A)', 'Usage: ./main <input_width> <repeat>'
    stop 1
  end if

  input_width = read_arg(1)
  input_width = ((input_width + max_block_size - 1) / max_block_size) * max_block_size
  repeat = read_arg(2)

  do mask_width = 3, max_mask_width - 1, 2
    print '(A)', ''
    print '(A)', '---------------------'
    print '(A,I0)', 'Mask width: ', mask_width

    print '(A)', '1D convolution (FP64)'
    call conv1d_real64(input_width, mask_width, repeat)

    print '(A)', '1D convolution (FP32)'
    call conv1d_real32(input_width, mask_width, repeat)

    print '(A)', '1D convolution (INT16)'
    call conv1d_int16(input_width, mask_width, repeat)
  end do

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine conv1d_real64(input_width, mask_width, repeat)
    integer, intent(in) :: input_width, mask_width, repeat
    real(real64), allocatable :: a(:), b(:), mask(:)
    integer :: i, bs

    allocate(a(input_width), b(input_width), mask(max_mask_width))
    mask = 1.0_real64
    call c_srand(123_c_int)
    do i = 1, input_width
      a(i) = real(mod(c_rand(), 256_c_int), real64)
    end do

    !$omp target data map(to: a(1:input_width), mask(1:mask_width)) map(alloc: b(1:input_width))
    do bs = 64, max_block_size, 64
      if (iand(bs, bs - 1) == 0) call run_real64_case('conv1d kernel', bs, repeat, a, b, mask, input_width, mask_width)
    end do
    do bs = 64, max_block_size, 64
      if (iand(bs, bs - 1) == 0) call run_real64_case('conv1d-tiled kernel', bs, repeat, a, b, mask, input_width, mask_width)
    end do
    do bs = 64, max_block_size, 64
      if (iand(bs, bs - 1) == 0) call run_real64_case('conv1d-tiled-caching kernel', bs, repeat, a, b, mask, input_width, mask_width)
    end do
    !$omp end target data
    deallocate(a, b, mask)
  end subroutine conv1d_real64

  subroutine conv1d_real32(input_width, mask_width, repeat)
    integer, intent(in) :: input_width, mask_width, repeat
    real(real32), allocatable :: a(:), b(:), mask(:)
    integer :: i, bs

    allocate(a(input_width), b(input_width), mask(max_mask_width))
    mask = 1.0_real32
    call c_srand(123_c_int)
    do i = 1, input_width
      a(i) = real(mod(c_rand(), 256_c_int), real32)
    end do

    !$omp target data map(to: a(1:input_width), mask(1:mask_width)) map(alloc: b(1:input_width))
    do bs = 64, max_block_size, 64
      if (iand(bs, bs - 1) == 0) call run_real32_case('conv1d kernel', bs, repeat, a, b, mask, input_width, mask_width)
    end do
    do bs = 64, max_block_size, 64
      if (iand(bs, bs - 1) == 0) call run_real32_case('conv1d-tiled kernel', bs, repeat, a, b, mask, input_width, mask_width)
    end do
    do bs = 64, max_block_size, 64
      if (iand(bs, bs - 1) == 0) call run_real32_case('conv1d-tiled-caching kernel', bs, repeat, a, b, mask, input_width, mask_width)
    end do
    !$omp end target data
    deallocate(a, b, mask)
  end subroutine conv1d_real32

  subroutine conv1d_int16(input_width, mask_width, repeat)
    integer, intent(in) :: input_width, mask_width, repeat
    integer(int16), allocatable :: a(:), b(:), mask(:)
    integer :: i, bs

    allocate(a(input_width), b(input_width), mask(max_mask_width))
    mask = 1_int16
    call c_srand(123_c_int)
    do i = 1, input_width
      a(i) = int(mod(c_rand(), 256_c_int), int16)
    end do

    !$omp target data map(to: a(1:input_width), mask(1:mask_width)) map(alloc: b(1:input_width))
    do bs = 64, max_block_size, 64
      if (iand(bs, bs - 1) == 0) call run_int16_case('conv1d kernel', bs, repeat, a, b, mask, input_width, mask_width)
    end do
    do bs = 64, max_block_size, 64
      if (iand(bs, bs - 1) == 0) call run_int16_case('conv1d-tiled kernel', bs, repeat, a, b, mask, input_width, mask_width)
    end do
    do bs = 64, max_block_size, 64
      if (iand(bs, bs - 1) == 0) call run_int16_case('conv1d-tiled-caching kernel', bs, repeat, a, b, mask, input_width, mask_width)
    end do
    !$omp end target data
    deallocate(a, b, mask)
  end subroutine conv1d_int16

  subroutine run_real64_case(label, bs, repeat, a, b, mask, input_width, mask_width)
    character(len=*), intent(in) :: label
    integer, intent(in) :: bs, repeat, input_width, mask_width
    real(real64), intent(in) :: a(:), mask(:)
    real(real64), intent(inout) :: b(:)
    integer :: rep
    real(real64) :: start_time, end_time

    start_time = omp_get_wtime()
    do rep = 1, repeat
      call conv_kernel_real64(mask, a, b, input_width, mask_width, bs)
    end do
    end_time = omp_get_wtime()
    print '(A,A,A,F0.6,A)', 'Average kernel execution time of ', label, ': ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(b(1:input_width))
    call reference_real64(a, b, mask, input_width, mask_width)
  end subroutine run_real64_case

  subroutine run_real32_case(label, bs, repeat, a, b, mask, input_width, mask_width)
    character(len=*), intent(in) :: label
    integer, intent(in) :: bs, repeat, input_width, mask_width
    real(real32), intent(in) :: a(:), mask(:)
    real(real32), intent(inout) :: b(:)
    integer :: rep
    real(real64) :: start_time, end_time

    start_time = omp_get_wtime()
    do rep = 1, repeat
      call conv_kernel_real32(mask, a, b, input_width, mask_width, bs)
    end do
    end_time = omp_get_wtime()
    print '(A,A,A,F0.6,A)', 'Average kernel execution time of ', label, ': ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(b(1:input_width))
    call reference_real32(a, b, mask, input_width, mask_width)
  end subroutine run_real32_case

  subroutine run_int16_case(label, bs, repeat, a, b, mask, input_width, mask_width)
    character(len=*), intent(in) :: label
    integer, intent(in) :: bs, repeat, input_width, mask_width
    integer(int16), intent(in) :: a(:), mask(:)
    integer(int16), intent(inout) :: b(:)
    integer :: rep
    real(real64) :: start_time, end_time

    start_time = omp_get_wtime()
    do rep = 1, repeat
      call conv_kernel_int16(mask, a, b, input_width, mask_width, bs)
    end do
    end_time = omp_get_wtime()
    print '(A,A,A,F0.6,A)', 'Average kernel execution time of ', label, ': ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(b(1:input_width))
    call reference_int16(a, b, mask, input_width, mask_width)
  end subroutine run_int16_case

  subroutine conv_kernel_real64(mask, a, b, input_width, mask_width, bs)
    real(real64), intent(in) :: mask(:), a(:)
    real(real64), intent(inout) :: b(:)
    integer, intent(in) :: input_width, mask_width, bs
    integer :: i, j, start_idx, in_idx
    real(real64) :: s
    !$omp target teams distribute parallel do num_threads(bs) private(j, start_idx, in_idx, s)
    do i = 1, input_width
      start_idx = i - 1 - mask_width / 2
      s = 0.0_real64
      do j = 1, mask_width
        in_idx = start_idx + j
        if (in_idx >= 1 .and. in_idx <= input_width) s = s + a(in_idx) * mask(j)
      end do
      b(i) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine conv_kernel_real64

  subroutine conv_kernel_real32(mask, a, b, input_width, mask_width, bs)
    real(real32), intent(in) :: mask(:), a(:)
    real(real32), intent(inout) :: b(:)
    integer, intent(in) :: input_width, mask_width, bs
    integer :: i, j, start_idx, in_idx
    real(real32) :: s
    !$omp target teams distribute parallel do num_threads(bs) private(j, start_idx, in_idx, s)
    do i = 1, input_width
      start_idx = i - 1 - mask_width / 2
      s = 0.0_real32
      do j = 1, mask_width
        in_idx = start_idx + j
        if (in_idx >= 1 .and. in_idx <= input_width) s = s + a(in_idx) * mask(j)
      end do
      b(i) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine conv_kernel_real32

  subroutine conv_kernel_int16(mask, a, b, input_width, mask_width, bs)
    integer(int16), intent(in) :: mask(:), a(:)
    integer(int16), intent(inout) :: b(:)
    integer, intent(in) :: input_width, mask_width, bs
    integer :: i, j, start_idx, in_idx
    integer(int32) :: s
    !$omp target teams distribute parallel do num_threads(bs) private(j, start_idx, in_idx, s)
    do i = 1, input_width
      start_idx = i - 1 - mask_width / 2
      s = 0_int32
      do j = 1, mask_width
        in_idx = start_idx + j
        if (in_idx >= 1 .and. in_idx <= input_width) s = s + int(a(in_idx), int32) * int(mask(j), int32)
      end do
      b(i) = int(s, int16)
    end do
    !$omp end target teams distribute parallel do
  end subroutine conv_kernel_int16

  subroutine reference_real64(a, b, mask, input_width, mask_width)
    real(real64), intent(in) :: a(:), b(:), mask(:)
    integer, intent(in) :: input_width, mask_width
    integer :: i, j, start_idx, in_idx
    real(real64) :: s
    logical :: ok
    ok = .true.
    do i = 1, input_width
      start_idx = i - 1 - mask_width / 2
      s = 0.0_real64
      do j = 1, mask_width
        in_idx = start_idx + j
        if (in_idx >= 1 .and. in_idx <= input_width) s = s + a(in_idx) * mask(j)
      end do
      if (abs(s - b(i)) > 1.0e-3_real64) ok = .false.
    end do
    call print_status(ok)
  end subroutine reference_real64

  subroutine reference_real32(a, b, mask, input_width, mask_width)
    real(real32), intent(in) :: a(:), b(:), mask(:)
    integer, intent(in) :: input_width, mask_width
    integer :: i, j, start_idx, in_idx
    real(real32) :: s
    logical :: ok
    ok = .true.
    do i = 1, input_width
      start_idx = i - 1 - mask_width / 2
      s = 0.0_real32
      do j = 1, mask_width
        in_idx = start_idx + j
        if (in_idx >= 1 .and. in_idx <= input_width) s = s + a(in_idx) * mask(j)
      end do
      if (abs(s - b(i)) > 1.0e-3_real32) ok = .false.
    end do
    call print_status(ok)
  end subroutine reference_real32

  subroutine reference_int16(a, b, mask, input_width, mask_width)
    integer(int16), intent(in) :: a(:), b(:), mask(:)
    integer, intent(in) :: input_width, mask_width
    integer :: i, j, start_idx, in_idx
    integer(int32) :: s
    logical :: ok
    ok = .true.
    do i = 1, input_width
      start_idx = i - 1 - mask_width / 2
      s = 0_int32
      do j = 1, mask_width
        in_idx = start_idx + j
        if (in_idx >= 1 .and. in_idx <= input_width) s = s + int(a(in_idx), int32) * int(mask(j), int32)
      end do
      if (int(b(i), int32) /= s) ok = .false.
    end do
    call print_status(ok)
  end subroutine reference_int16

  subroutine print_status(ok)
    logical, intent(in) :: ok
    if (ok) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
    end if
  end subroutine print_status

end program main
