program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2, arg3, arg4
  integer :: i_img_width, i_img_height, i_img_count, repeat
  integer :: hstride, vstride, o_img_width, o_img_height
  integer :: size_image, size_output, total_input, total_output
  integer :: i
  integer(c_int) :: c_rand_value
  real(real32), allocatable :: h_image(:), h_output(:), d_output(:)
  real(real64) :: start_time, end_time, avg_time
  logical :: ok

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name="rand") result(r)
      import :: c_int
      integer(c_int) :: r
    end function c_rand
  end interface

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    print '(3A)', 'Usage: ', trim(arg0), ' <image width> <image height> <image count> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  call get_command_argument(4, arg4)
  read(arg1, *) i_img_width
  read(arg2, *) i_img_height
  read(arg3, *) i_img_count
  read(arg4, *) repeat

  if (i_img_width <= 0 .or. i_img_height <= 0 .or. i_img_count <= 0 .or. repeat <= 0) stop 1

  hstride = 2
  vstride = 2
  o_img_width = i_img_width / hstride
  o_img_height = i_img_height / vstride

  print '(A,I0,A,I0)', 'input image width ', i_img_width, ' Hstride ', hstride
  print '(A,I0,A,I0)', 'input image height ', i_img_height, ' Vstride ', vstride
  print '(A,I0)', 'output image width ', o_img_width
  print '(A,I0)', 'output image height ', o_img_height

  size_image = i_img_width * i_img_height
  size_output = o_img_width * o_img_height
  total_input = size_image * i_img_count
  total_output = size_output * i_img_count

  allocate(h_image(total_input), h_output(total_output), d_output(total_output))

  call c_srand(2_c_int)
  do i = 1, total_input
    c_rand_value = c_rand()
    h_image(i) = real(mod(c_rand_value, 256_c_int), real32) / 255.0_real32
  end do
  h_output = 0.0_real32
  d_output = 0.0_real32

  !$omp target data map(to: h_image(1:total_input)) map(from: d_output(1:total_output))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call maxpool_device(h_image, d_output, i_img_width, i_img_height, i_img_count, o_img_width, o_img_height)
  end do
  end_time = omp_get_wtime()
  avg_time = (end_time - start_time) / real(repeat, real64)
  print '(A,F8.6,A)', 'Average kernel execution time: ', avg_time, ' (s)'
  !$omp end target data

  call maxpool_host(h_image, h_output, i_img_width, i_img_height, i_img_count, o_img_width, o_img_height)

  ok = .true.
  do i = 1, total_output
    if (h_output(i) /= d_output(i)) then
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

  deallocate(h_image, h_output, d_output)

contains

  subroutine maxpool_device(h_image, d_output, i_img_width, i_img_height, i_img_count, o_img_width, o_img_height)
    real(real32), intent(in) :: h_image(:)
    real(real32), intent(out) :: d_output(:)
    integer, intent(in) :: i_img_width, i_img_height, i_img_count, o_img_width, o_img_height
    integer :: z, y, x, r, c, xidx, yidx, idx_in, idx_out, size_image, size_output
    real(real32) :: maxval

    size_image = i_img_width * i_img_height
    size_output = o_img_width * o_img_height

    !$omp target teams distribute parallel do collapse(3) thread_limit(256) private(xidx, yidx, idx_in, idx_out, r, c, maxval)
    do z = 0, i_img_count - 1
      do y = 0, o_img_height - 1
        do x = 0, o_img_width - 1
          xidx = 2 * x
          yidx = 2 * y
          maxval = 0.0_real32
          do r = 0, 1
            do c = 0, 1
              idx_in = z * size_image + (yidx + r) * i_img_width + xidx + c + 1
              maxval = max(maxval, h_image(idx_in))
            end do
          end do
          idx_out = z * size_output + y * o_img_width + x + 1
          d_output(idx_out) = maxval
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine maxpool_device

  subroutine maxpool_host(h_image, h_output, i_img_width, i_img_height, i_img_count, o_img_width, o_img_height)
    real(real32), intent(in) :: h_image(:)
    real(real32), intent(out) :: h_output(:)
    integer, intent(in) :: i_img_width, i_img_height, i_img_count, o_img_width, o_img_height
    integer :: z, y, x, r, c, xidx, yidx, idx_in, idx_out, size_image, size_output
    real(real32) :: maxval

    size_image = i_img_width * i_img_height
    size_output = o_img_width * o_img_height

    do z = 0, i_img_count - 1
      do y = 0, o_img_height - 1
        do x = 0, o_img_width - 1
          xidx = 2 * x
          yidx = 2 * y
          maxval = 0.0_real32
          do r = 0, 1
            do c = 0, 1
              idx_in = z * size_image + (yidx + r) * i_img_width + xidx + c + 1
              maxval = max(maxval, h_image(idx_in))
            end do
          end do
          idx_out = z * size_output + y * o_img_width + x + 1
          h_output(idx_out) = maxval
        end do
      end do
    end do
  end subroutine maxpool_host

end program main
