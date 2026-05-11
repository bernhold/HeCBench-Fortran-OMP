program main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg
  integer :: width, height, img_size, num_detections, i
  real(real32), allocatable :: input_x(:), input_y(:), input_z(:)
  real(real32), allocatable :: output_x(:), output_y(:), output_z(:)
  real(real32), allocatable :: ref_x(:), ref_y(:), ref_z(:)
  integer, allocatable :: box_width(:), box_height(:), box_left(:), box_top(:)
  real(real64) :: start_time, end_time
  logical :: ok

  if (command_argument_count() /= 2) then
    call get_command_argument(0, arg)
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg)
    write(*,'(A)') ' <width> <height>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) width
  call get_command_argument(2, arg)
  read(arg, *) height
  if (width <= 64 .or. height <= 64) stop 1

  img_size = width * height
  num_detections = int(real(img_size, real32) * 0.8_real32)

  allocate(input_x(0:img_size - 1), input_y(0:img_size - 1), input_z(0:img_size - 1), &
           output_x(0:img_size - 1), output_y(0:img_size - 1), output_z(0:img_size - 1), &
           ref_x(0:img_size - 1), ref_y(0:img_size - 1), ref_z(0:img_size - 1), &
           box_width(0:num_detections - 1), box_height(0:num_detections - 1), &
           box_left(0:num_detections - 1), box_top(0:num_detections - 1))

  call initialize_inputs(input_x, input_y, input_z, output_x, output_y, output_z, &
                         ref_x, ref_y, ref_z, box_width, box_height, box_left, box_top, &
                         img_size, width, height, num_detections)

  !$omp target data map(to: input_x(0:img_size - 1), input_y(0:img_size - 1), input_z(0:img_size - 1)) &
  !$omp& map(tofrom: output_x(0:img_size - 1), output_y(0:img_size - 1), output_z(0:img_size - 1))
  start_time = omp_get_wtime()
  do i = 0, num_detections - 1
    call detection_overlay_box(input_x, input_y, input_z, output_x, output_y, output_z, width, height, &
                               box_left(i), box_top(i), box_width(i), box_height(i))
  end do
  end_time = omp_get_wtime()
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Total kernel execution time: ', end_time - start_time, ' (s)'

  call reference_overlay(input_x, input_y, input_z, ref_x, ref_y, ref_z, width, height, &
                         box_width, box_height, box_left, box_top, num_detections)

  ok = .true.
  do i = 0, img_size - 1
    if (abs(ref_x(i) - output_x(i)) > 1.0e-3_real32 .or. &
        abs(ref_y(i) - output_y(i)) > 1.0e-3_real32 .or. &
        abs(ref_z(i) - output_z(i)) > 1.0e-3_real32) then
      write(*,'(A,I0)') 'Error at index ', i
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(input_x, input_y, input_z, output_x, output_y, output_z, ref_x, ref_y, ref_z, &
             box_width, box_height, box_left, box_top)

contains

  subroutine initialize_inputs(input_x, input_y, input_z, output_x, output_y, output_z, &
                               ref_x, ref_y, ref_z, box_width, box_height, box_left, box_top, &
                               img_size, width, height, num_detections)
    real(real32), intent(out) :: input_x(0:), input_y(0:), input_z(0:)
    real(real32), intent(out) :: output_x(0:), output_y(0:), output_z(0:)
    real(real32), intent(out) :: ref_x(0:), ref_y(0:), ref_z(0:)
    integer, intent(out) :: box_width(0:), box_height(0:), box_left(0:), box_top(0:)
    integer, intent(in) :: img_size, width, height, num_detections
    integer :: i
    integer(int64) :: state

    state = 123_int64
    do i = 0, img_size - 1
      input_x(i) = real(next_rand(state, 256), real32)
      input_y(i) = real(next_rand(state, 256), real32)
      input_z(i) = real(next_rand(state, 256), real32)
      output_x(i) = input_x(i)
      output_y(i) = input_y(i)
      output_z(i) = input_z(i)
      ref_x(i) = input_x(i)
      ref_y(i) = input_y(i)
      ref_z(i) = input_z(i)
    end do

    do i = 0, num_detections - 1
      box_width(i) = 64 + next_rand(state, 128)
      box_height(i) = 64 + next_rand(state, 128)
      box_left(i) = next_rand(state, width - 64)
      box_top(i) = next_rand(state, height - 64)
    end do
  end subroutine initialize_inputs

  integer function next_rand(state, modulus)
    integer(int64), intent(inout) :: state
    integer, intent(in) :: modulus

    state = mod(1103515245_int64 * state + 12345_int64, 2147483648_int64)
    next_rand = int(mod(state, int(modulus, int64)))
  end function next_rand

  subroutine detection_overlay_box(input_x, input_y, input_z, output_x, output_y, output_z, &
                                   img_width, img_height, x0, y0, box_width, box_height)
    real(real32), intent(in) :: input_x(0:), input_y(0:), input_z(0:)
    real(real32), intent(inout) :: output_x(0:), output_y(0:), output_z(0:)
    integer, intent(in) :: img_width, img_height, x0, y0, box_width, box_height
    integer :: box_x, box_y, x, y, idx
    real(real32), parameter :: color_x = 255.0_real32, color_y = 204.0_real32, color_z = 203.0_real32
    real(real32), parameter :: alpha = 1.0_real32 / 255.0_real32, ialph = 1.0_real32 - alpha

    !$omp target teams distribute parallel do collapse(2) thread_limit(64) &
    !$omp& private(box_x, box_y, x, y, idx)
    do box_y = 0, box_height - 1
      do box_x = 0, box_width - 1
        x = box_x + x0
        y = box_y + y0
        if (x < img_width .and. y < img_height) then
          idx = y * img_width + x
          output_x(idx) = alpha * color_x + ialph * input_x(idx)
          output_y(idx) = alpha * color_y + ialph * input_y(idx)
          output_z(idx) = alpha * color_z + ialph * input_z(idx)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine detection_overlay_box

  subroutine reference_overlay(input_x, input_y, input_z, output_x, output_y, output_z, &
                               img_width, img_height, box_width, box_height, box_left, box_top, num_detections)
    real(real32), intent(in) :: input_x(0:), input_y(0:), input_z(0:)
    real(real32), intent(inout) :: output_x(0:), output_y(0:), output_z(0:)
    integer, intent(in) :: img_width, img_height, num_detections
    integer, intent(in) :: box_width(0:), box_height(0:), box_left(0:), box_top(0:)
    integer :: n, box_x, box_y, x, y, idx
    real(real32), parameter :: color_x = 255.0_real32, color_y = 204.0_real32, color_z = 203.0_real32
    real(real32), parameter :: alpha = 1.0_real32 / 255.0_real32, ialph = 1.0_real32 - alpha

    do n = 0, num_detections - 1
      do box_y = 0, box_height(n) - 1
        do box_x = 0, box_width(n) - 1
          x = box_x + box_left(n)
          y = box_y + box_top(n)
          if (x < img_width .and. y < img_height) then
            idx = y * img_width + x
            output_x(idx) = alpha * color_x + ialph * input_x(idx)
            output_y(idx) = alpha * color_y + ialph * input_y(idx)
            output_z(idx) = alpha * color_z + ialph * input_z(idx)
          end if
        end do
      end do
    end do
  end subroutine reference_overlay

end program main
