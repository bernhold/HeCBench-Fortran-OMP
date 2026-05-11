program main
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use omp_lib, only: omp_get_wtime
  implicit none

  integer, parameter :: kernel_radius = 8
  integer, parameter :: kernel_length = 2 * kernel_radius + 1

  integer :: argc, image_w, image_h, repeat
  integer(int64) :: n
  character(len=256) :: arg0
  real(real32), allocatable :: kernel(:), input(:), buffer(:), output_cpu(:), output_gpu(:)
  real(real64) :: sum_ref, delta, l2norm
  integer :: i

  call get_command_argument(0, arg0)
  argc = command_argument_count()
  if (argc /= 3) then
    write(*,'("Usage: ",A," <image width> <image height> <repeat>")') trim(arg0)
    stop 1
  end if

  image_w = read_arg(1)
  image_h = read_arg(2)
  repeat = read_arg(3)
  if (image_w <= 0 .or. image_h <= 0 .or. repeat <= 0) then
    write(*,'("Usage: ",A," <image width> <image height> <repeat>")') trim(arg0)
    stop 1
  end if

  n = int(image_w, int64) * int(image_h, int64)
  if (n > huge(i)) then
    write(*,'("Image is too large")')
    stop 1
  end if

  allocate(kernel(kernel_length), input(n), buffer(n), output_cpu(n), output_gpu(n))
  call initialize_inputs(kernel, input, int(n))

  !$omp target data map(to: kernel(1:kernel_length), input(1:n)) &
  !$omp& map(alloc: buffer(1:n)) map(from: output_gpu(1:n))
    call convolution_rows_device(buffer, input, kernel, image_w, image_h)
    call convolution_columns_device(output_gpu, buffer, kernel, image_w, image_h)

    call run_timed_convolution(output_gpu, buffer, input, kernel, image_w, image_h, repeat)
  !$omp end target data

  write(*,'("Comparing against Host/C++ computation...")')
  call convolution_row_host(buffer, input, kernel, image_w, image_h)
  call convolution_column_host(output_cpu, buffer, kernel, image_w, image_h)

  delta = 0.0_real64
  sum_ref = 0.0_real64
  do i = 1, int(n)
    delta = delta + real(output_cpu(i) - output_gpu(i), real64) * real(output_cpu(i) - output_gpu(i), real64)
    sum_ref = sum_ref + real(output_cpu(i), real64) * real(output_cpu(i), real64)
  end do
  if (sum_ref > 0.0_real64) then
    l2norm = sqrt(delta / sum_ref)
  else
    l2norm = sqrt(delta)
  end if
  write(*,'("Relative L2 norm: ",ES10.3)') l2norm
  write(*,'()')

  if (l2norm < 1.0e-6_real64) then
    write(*,'("PASS")')
  else
    write(*,'("FAIL")')
    stop 1
  end if

  deallocate(kernel, input, buffer, output_cpu, output_gpu)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer_arg

    call get_command_argument(position, buffer_arg)
    read(buffer_arg, *) read_arg
  end function read_arg

  subroutine initialize_inputs(kernel, input, n)
    real(real32), intent(out) :: kernel(:), input(:)
    integer, intent(in) :: n
    integer(int64) :: state
    integer :: idx

    state = 2009_int64
    do idx = 1, kernel_length
      kernel(idx) = real(c_rand_mod(state, 16), real32)
    end do

    do idx = 1, n
      input(idx) = real(c_rand_mod(state, 16), real32)
    end do
  end subroutine initialize_inputs

  integer function c_rand_mod(state, modulus)
    integer(int64), intent(inout) :: state
    integer, intent(in) :: modulus

    state = modulo(1103515245_int64 * state + 12345_int64, 2147483648_int64)
    c_rand_mod = int(modulo(state / 65536_int64, int(modulus, int64)))
  end function c_rand_mod

  subroutine run_timed_convolution(output, buffer, input, kernel, image_w, image_h, repeat)
    real(real32), intent(inout) :: output(:), buffer(:)
    real(real32), intent(in) :: input(:), kernel(:)
    integer, intent(in) :: image_w, image_h, repeat
    integer :: iter
    real(real64) :: start_time, end_time

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call convolution_rows_device(buffer, input, kernel, image_w, image_h)
      call convolution_columns_device(output, buffer, kernel, image_w, image_h)
    end do
    end_time = omp_get_wtime()

    write(*,'("Average kernel execution time ",F0.6," (s)")') (end_time - start_time) / real(repeat, real64)
  end subroutine run_timed_convolution

  subroutine convolution_rows_device(dst, src, kernel, image_w, image_h)
    real(real32), intent(inout) :: dst(:)
    real(real32), intent(in) :: src(:), kernel(:)
    integer, intent(in) :: image_w, image_h
    integer :: x, y, k, d, idx
    real(real32) :: accum

    !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(k, d, idx, accum)
    do y = 0, image_h - 1
      do x = 0, image_w - 1
        accum = 0.0_real32
        do k = -kernel_radius, kernel_radius
          d = x + k
          if (d >= 0 .and. d < image_w) then
            accum = accum + src(y * image_w + d + 1) * kernel(kernel_radius - k + 1)
          end if
        end do
        idx = y * image_w + x + 1
        dst(idx) = accum
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine convolution_rows_device

  subroutine convolution_columns_device(dst, src, kernel, image_w, image_h)
    real(real32), intent(inout) :: dst(:)
    real(real32), intent(in) :: src(:), kernel(:)
    integer, intent(in) :: image_w, image_h
    integer :: x, y, k, d, idx
    real(real32) :: accum

    !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(k, d, idx, accum)
    do y = 0, image_h - 1
      do x = 0, image_w - 1
        accum = 0.0_real32
        do k = -kernel_radius, kernel_radius
          d = y + k
          if (d >= 0 .and. d < image_h) then
            accum = accum + src(d * image_w + x + 1) * kernel(kernel_radius - k + 1)
          end if
        end do
        idx = y * image_w + x + 1
        dst(idx) = accum
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine convolution_columns_device

  subroutine convolution_row_host(dst, src, kernel, image_w, image_h)
    real(real32), intent(out) :: dst(:)
    real(real32), intent(in) :: src(:), kernel(:)
    integer, intent(in) :: image_w, image_h
    integer :: x, y, k, d
    real(real64) :: accum

    do y = 0, image_h - 1
      do x = 0, image_w - 1
        accum = 0.0_real64
        do k = -kernel_radius, kernel_radius
          d = x + k
          if (d >= 0 .and. d < image_w) then
            accum = accum + real(src(y * image_w + d + 1), real64) * real(kernel(kernel_radius - k + 1), real64)
          end if
        end do
        dst(y * image_w + x + 1) = real(accum, real32)
      end do
    end do
  end subroutine convolution_row_host

  subroutine convolution_column_host(dst, src, kernel, image_w, image_h)
    real(real32), intent(out) :: dst(:)
    real(real32), intent(in) :: src(:), kernel(:)
    integer, intent(in) :: image_w, image_h
    integer :: x, y, k, d
    real(real64) :: accum

    do y = 0, image_h - 1
      do x = 0, image_w - 1
        accum = 0.0_real64
        do k = -kernel_radius, kernel_radius
          d = y + k
          if (d >= 0 .and. d < image_h) then
            accum = accum + real(src(d * image_w + x + 1), real64) * real(kernel(kernel_radius - k + 1), real64)
          end if
        end do
        dst(y * image_w + x + 1) = real(accum, real32)
      end do
    end do
  end subroutine convolution_column_host

end program main
