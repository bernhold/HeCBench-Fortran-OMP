program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0
  integer, parameter :: n_vec = 4
  integer, parameter :: vec_len(n_vec) = [1, 2, 4, 8]
  integer :: count, repeat, vl, iv, i
  real(real32), allocatable :: h_gradient(:), h_feature(:), h_backprop(:), r_backprop(:)
  integer(int32), allocatable :: h_in(:), h_out(:), r_out(:)
  real(real64) :: start_time, end_time, avg_us
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage: ', trim(arg0), ' <count> <repeat>'
    stop 1
  end if

  count = read_arg(1)
  repeat = read_arg(2)
  if (count <= 0 .or. repeat <= 0) stop 1

  allocate(h_gradient(count), h_feature(count), h_backprop(count), r_backprop(count))

  do i = 1, count
    h_feature(i) = deterministic_signed_unit_float(i)
    h_gradient(i) = 1.0_real32
  end do

  call relu_grad_reference(count, h_gradient, h_feature, r_backprop)

  !$omp target data map(to: h_gradient(1:count), h_feature(1:count)) map(alloc: h_backprop(1:count))
  do iv = 1, n_vec
    vl = vec_len(iv)
    start_time = omp_get_wtime()
    do i = 1, repeat
      call relu_grad_device(h_gradient, h_feature, h_backprop, count, vl)
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)
    print '(A,I0,A,F0.6,A)', 'Average execution time of ReluGrad_impl (VL=', vl, '): ', avg_us, ' (us)'

    !$omp target update from(h_backprop(1:count))
    ok = .true.
    do i = 1, count
      if (abs(h_backprop(i) - r_backprop(i)) > 1.0e-3_real32) then
        ok = .false.
        exit
      end if
    end do
    call print_status(ok)
  end do
  !$omp end target data

  deallocate(h_gradient, h_feature, h_backprop, r_backprop)

  allocate(h_in(count), h_out(count), r_out(count))
  do i = 1, count
    h_in(i) = make_packed_input(i)
  end do

  call relu_reference(count, h_in, r_out)

  !$omp target data map(to: h_in(1:count)) map(alloc: h_out(1:count))
  do iv = 1, n_vec
    vl = vec_len(iv)
    start_time = omp_get_wtime()
    do i = 1, repeat
      call relu_device(count, h_in, h_out, vl)
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)
    print '(A,I0,A,F0.6,A)', 'Average execution time of Relu_impl (VL=', vl, '): ', avg_us, ' (us)'

    !$omp target update from(h_out(1:count))
    ok = .true.
    do i = 1, count
      if (h_out(i) /= r_out(i)) then
        ok = .false.
        exit
      end if
    end do
    call print_status(ok)
  end do
  !$omp end target data

  deallocate(h_in, h_out, r_out)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer

    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine print_status(ok)
    logical, intent(in) :: ok

    if (ok) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
      stop 1
    end if
  end subroutine print_status

  real(real32) function deterministic_signed_unit_float(i)
    integer, intent(in) :: i
    integer(int64) :: value

    value = modulo(1103515245_int64 * int(i, int64) + 12345_int64, 2147483647_int64)
    deterministic_signed_unit_float = 2.0_real32 * (real(value, real32) / 2147483647.0_real32) - 1.0_real32
  end function deterministic_signed_unit_float

  integer(int32) function make_packed_input(i)
    integer, intent(in) :: i
    integer(int32) :: b0, b1, b2, b3

    b0 = int(mod(17 * i + 3, 256), int32)
    b1 = int(mod(31 * i + 5, 256), int32)
    b2 = int(mod(47 * i + 7, 256), int32)
    b3 = int(mod(61 * i + 11, 256), int32)
    make_packed_input = ior(ior(b0, ishft(b1, 8)), ior(ishft(b2, 16), ishft(b3, 24)))
  end function make_packed_input

  integer(int32) function relu_byte(byte_value)
    integer(int32), intent(in) :: byte_value
    integer(int32) :: signed_value

    signed_value = iand(byte_value, int(z'000000FF', int32))
    if (signed_value >= 128_int32) signed_value = signed_value - 256_int32
    if (signed_value > 0_int32) then
      relu_byte = signed_value
    else
      relu_byte = 0_int32
    end if
  end function relu_byte

  subroutine relu_grad_reference(count, gradient, feature, backprop)
    integer, intent(in) :: count
    real(real32), intent(in) :: gradient(:), feature(:)
    real(real32), intent(out) :: backprop(:)
    integer :: i

    do i = 1, count
      if (feature(i) > 0.0_real32) then
        backprop(i) = gradient(i)
      else
        backprop(i) = 0.0_real32
      end if
    end do
  end subroutine relu_grad_reference

  subroutine relu_grad_device(gradient, feature, backprop, count, vl)
    real(real32), intent(in) :: gradient(:), feature(:)
    real(real32), intent(out) :: backprop(:)
    integer, intent(in) :: count, vl
    integer :: v_count, index, j, base

    v_count = count / vl
    !$omp target teams distribute parallel do thread_limit(256) private(j, base)
    do index = 0, v_count - 1
      base = index * vl
      do j = 1, vl
        if (feature(base + j) > 0.0_real32) then
          backprop(base + j) = gradient(base + j)
        else
          backprop(base + j) = 0.0_real32
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine relu_grad_device

  subroutine relu_reference(count, input, output)
    integer, intent(in) :: count
    integer(int32), intent(in) :: input(:)
    integer(int32), intent(out) :: output(:)
    integer :: i

    do i = 1, count
      output(i) = relu_packed(input(i))
    end do
  end subroutine relu_reference

  integer(int32) function relu_packed(input)
    integer(int32), intent(in) :: input
    integer(int32) :: r0, r1, r2, r3

    r0 = relu_byte(input)
    r1 = relu_byte(ishft(input, -8))
    r2 = relu_byte(ishft(input, -16))
    r3 = relu_byte(ishft(input, -24))
    relu_packed = ior(ior(r0, ishft(r1, 8)), ior(ishft(r2, 16), ishft(r3, 24)))
  end function relu_packed

  subroutine relu_device(count, input, output, vl)
    integer, intent(in) :: count, vl
    integer(int32), intent(in) :: input(:)
    integer(int32), intent(out) :: output(:)
    integer :: v_count, index, j, base

    v_count = count / vl
    !$omp target teams distribute parallel do thread_limit(256) private(j, base)
    do index = 0, v_count - 1
      base = index * vl
      do j = 1, vl
        output(base + j) = relu_packed(input(base + j))
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine relu_device

end program main
