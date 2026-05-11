program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg
  integer :: width, height, block_size, repeat
  integer :: image_size, i, iter, errors, max_error, err
  integer, allocatable :: image(:), reference(:), pixel(:)
  real(real64) :: start_time, total_time

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <image width> <image height> <block size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) width
  call get_command_argument(2, arg); read(arg, *) height
  call get_command_argument(3, arg); read(arg, *) block_size
  call get_command_argument(4, arg); read(arg, *) repeat

  if (width <= 0 .or. height <= 0 .or. block_size <= 0 .or. repeat <= 0) stop 1

  image_size = width * height
  allocate(image(image_size), reference(image_size), pixel(image_size))

  call fill_fractal(width, height, image)
  do i = 1, image_size
    reference(i) = gamma_value(image(i))
  end do

  pixel = image
  total_time = 0.0_real64

  !$omp target data map(tofrom: pixel(1:image_size))
  do iter = 1, repeat
    pixel = image
    !$omp target update to(pixel(1:image_size))

    start_time = omp_get_wtime()
    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 1, image_size
      pixel(i) = int(255.0_real32 * &
        (real(pixel(i), real32) / 255.0_real32) * &
        (real(pixel(i), real32) / 255.0_real32))
      if (pixel(i) > 255) pixel(i) = 255
    end do
    !$omp end target teams distribute parallel do
    total_time = total_time + (omp_get_wtime() - start_time)
  end do
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average kernel execution time ', total_time / real(repeat, real64), ' (s)'

  errors = 0
  max_error = 0
  do i = 1, image_size
    err = abs(pixel(i) - reference(i))
    if (err /= 0) then
      errors = errors + 1
      if (err > max_error) max_error = err
    end if
  end do

  if (errors == 0) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A,I0,A,I0)') 'FAIL errors=', errors, ' max_error=', max_error
  end if

  deallocate(image, reference, pixel)

contains

  subroutine fill_fractal(width, height, image)
    integer, intent(in) :: width, height
    integer, intent(out) :: image(:)
    integer :: x, y, idx
    real(real64) :: fractal_pixel

    do y = 0, height - 1
      do x = 0, width - 1
        idx = y * width + x + 1
        fractal_pixel = fractal_value(x, y, width, height)
        if (fractal_pixel < 0.0_real64) fractal_pixel = 0.0_real64
        if (fractal_pixel > 255.0_real64) fractal_pixel = 255.0_real64
        image(idx) = int(fractal_pixel)
      end do
    end do
  end subroutine fill_fractal

  real(real64) function fractal_value(x, y, width, height)
    integer, intent(in) :: x, y, width, height
    integer :: iter
    real(real64) :: fx, fy, res, nx, ny, val
    real(real64), parameter :: cx = -0.7436_real64
    real(real64), parameter :: cy = 0.1319_real64
    real(real64), parameter :: magn = 2000000.0_real64

    fx = (real(x, real64) - real(width, real64) / 2.0_real64) * (1.0_real64 / magn) + cx
    fy = (real(y, real64) - real(height, real64) / 2.0_real64) * (1.0_real64 / magn) + cy

    res = 0.0_real64
    nx = 0.0_real64
    ny = 0.0_real64
    val = 0.0_real64
    iter = 0
    do while (nx * nx + ny * ny <= 4.0_real64 .and. iter < 1000)
      val = nx * nx - ny * ny + fx
      ny = 2.0_real64 * nx * ny + fy
      nx = val
      res = res + exp(-sqrt(nx * nx + ny * ny))
      iter = iter + 1
    end do

    fractal_value = res
  end function fractal_value

  integer function gamma_value(value)
    integer, intent(in) :: value
    real(real32) :: v

    v = real(value, real32) / 255.0_real32
    gamma_value = int(255.0_real32 * v * v)
    if (gamma_value > 255) gamma_value = 255
  end function gamma_value

end program main
