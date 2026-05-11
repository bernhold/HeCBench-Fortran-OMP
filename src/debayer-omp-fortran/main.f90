program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use omp_lib
  implicit none

  integer, parameter :: tile_rows = 5
  integer, parameter :: tile_cols = 32
  integer, parameter :: kernel_size = 5
  integer, parameter :: rggb = 0

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

  character(len=256) :: arg0, arg
  integer :: width, height, repeat, num_pix, input_image_pitch, output_image_pitch
  integer, allocatable :: input(:), output(:), reference(:)
  integer :: i
  integer(int64) :: checksum
  real(real64) :: start_time, elapsed

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <width> <height> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) width
  call get_command_argument(2, arg); read(arg, *) height
  call get_command_argument(3, arg); read(arg, *) repeat
  if (width <= kernel_size .or. height <= kernel_size .or. repeat <= 0) stop 1

  input_image_pitch = width
  output_image_pitch = width * 4
  num_pix = width * height
  allocate(input(num_pix), output(4 * num_pix), reference(4 * num_pix))

  call c_srand(123_c_int)
  do i = 1, num_pix
    input(i) = modulo(c_rand(), 256_c_int)
  end do
  output = 0
  reference = 0

  !$omp target data map(to: input(1:num_pix)) map(from: output(1:4*num_pix))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call malvar_he_cutler_demosaic(height, width, input, input_image_pitch, output, &
      output_image_pitch, rggb)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  call reference_demosaic(height, width, input, input_image_pitch, reference, output_image_pitch, rggb)
  if (any(output /= reference)) error stop 'Fortran reference validation failed'

  write(*,'(A,F0.6,A)') 'Average kernel execution time ', elapsed / real(repeat, real64), ' (s)'

  checksum = 0_int64
  do i = 1, num_pix
    checksum = checksum + int(output(i), int64)
  end do
  write(*,'(A,I0)') 'Checksum: ', checksum

  deallocate(input, output, reference)

contains

  subroutine malvar_he_cutler_demosaic(height, width, input_image, input_pitch, output_image, &
      output_pitch, bayer_pattern)
    integer, intent(in) :: height, width, input_pitch, output_pitch, bayer_pattern
    integer, intent(in) :: input_image(:)
    integer, intent(inout) :: output_image(:)
    integer :: row, col

    !$omp target teams distribute parallel do collapse(2) thread_limit(tile_rows * tile_cols) &
    !$omp& map(to: input_image(1:height*width)) map(tofrom: output_image(1:4*height*width))
    do row = 0, height - 1
      do col = 0, width - 1
        call demosaic_pixel(row, col, height, width, input_image, input_pitch, output_image, &
          output_pitch, bayer_pattern)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine malvar_he_cutler_demosaic

  subroutine reference_demosaic(height, width, input_image, input_pitch, output_image, &
      output_pitch, bayer_pattern)
    integer, intent(in) :: height, width, input_pitch, output_pitch, bayer_pattern
    integer, intent(in) :: input_image(:)
    integer, intent(inout) :: output_image(:)
    integer :: row, col

    do row = 0, height - 1
      do col = 0, width - 1
        call demosaic_pixel(row, col, height, width, input_image, input_pitch, output_image, &
          output_pitch, bayer_pattern)
      end do
    end do
  end subroutine reference_demosaic

  subroutine demosaic_pixel(g_r, g_c, height, width, input_image, input_pitch, output_image, &
      output_pitch, bayer_pattern)
    integer, intent(in) :: g_r, g_c, height, width, input_pitch, output_pitch, bayer_pattern
    integer, intent(in) :: input_image(:)
    integer, intent(inout) :: output_image(:)
    integer :: f_ij, r1, r2, r3, r4
    integer :: green_at_red_or_blue, red_at_green_in_red, red_at_green_in_blue
    integer :: blue_at_green_in_red, blue_at_green_in_blue, red_at_blue, blue_at_red
    integer :: r_mod_2, c_mod_2, red_col, red_row, blue_col, blue_row
    logical :: in_red_row, in_blue_row, is_red_pixel, is_blue_pixel, is_green_pixel
    integer :: red_value, green_value, blue_value, out_idx

    f_ij = sample_pixel(input_image, height, width, input_pitch, g_r, g_c)

    r1 = (4 * f(g_r, g_c, input_image, height, width, input_pitch) + &
      2 * (f(g_r, g_c - 1, input_image, height, width, input_pitch) + &
           f(g_r - 1, g_c, input_image, height, width, input_pitch) + &
           f(g_r, g_c + 1, input_image, height, width, input_pitch) + &
           f(g_r + 1, g_c, input_image, height, width, input_pitch)) - &
      f(g_r, g_c - 2, input_image, height, width, input_pitch) - &
      f(g_r, g_c + 2, input_image, height, width, input_pitch) - &
      f(g_r - 2, g_c, input_image, height, width, input_pitch) - &
      f(g_r + 2, g_c, input_image, height, width, input_pitch)) / 8

    r2 = (8 * (f(g_r, g_c - 1, input_image, height, width, input_pitch) + &
               f(g_r, g_c + 1, input_image, height, width, input_pitch)) + &
      10 * f(g_r, g_c, input_image, height, width, input_pitch) + &
      f(g_r - 2, g_c, input_image, height, width, input_pitch) + &
      f(g_r + 2, g_c, input_image, height, width, input_pitch) - &
      2 * (f(g_r - 1, g_c - 1, input_image, height, width, input_pitch) + &
           f(g_r - 1, g_c + 1, input_image, height, width, input_pitch) + &
           f(g_r + 1, g_c - 1, input_image, height, width, input_pitch) + &
           f(g_r + 1, g_c + 1, input_image, height, width, input_pitch) + &
           f(g_r, g_c - 2, input_image, height, width, input_pitch) + &
           f(g_r, g_c + 2, input_image, height, width, input_pitch))) / 16

    r3 = (8 * (f(g_r - 1, g_c, input_image, height, width, input_pitch) + &
               f(g_r + 1, g_c, input_image, height, width, input_pitch)) + &
      10 * f(g_r, g_c, input_image, height, width, input_pitch) + &
      f(g_r, g_c - 2, input_image, height, width, input_pitch) + &
      f(g_r, g_c + 2, input_image, height, width, input_pitch) - &
      2 * (f(g_r - 1, g_c - 1, input_image, height, width, input_pitch) + &
           f(g_r - 1, g_c + 1, input_image, height, width, input_pitch) + &
           f(g_r + 1, g_c - 1, input_image, height, width, input_pitch) + &
           f(g_r + 1, g_c + 1, input_image, height, width, input_pitch) + &
           f(g_r - 2, g_c, input_image, height, width, input_pitch) + &
           f(g_r + 2, g_c, input_image, height, width, input_pitch))) / 16

    r4 = (12 * f(g_r, g_c, input_image, height, width, input_pitch) - &
      3 * (f(g_r, g_c - 2, input_image, height, width, input_pitch) + &
           f(g_r, g_c + 2, input_image, height, width, input_pitch) + &
           f(g_r - 2, g_c, input_image, height, width, input_pitch) + &
           f(g_r + 2, g_c, input_image, height, width, input_pitch)) + &
      4 * (f(g_r - 1, g_c - 1, input_image, height, width, input_pitch) + &
           f(g_r - 1, g_c + 1, input_image, height, width, input_pitch) + &
           f(g_r + 1, g_c - 1, input_image, height, width, input_pitch) + &
           f(g_r + 1, g_c + 1, input_image, height, width, input_pitch))) / 16

    green_at_red_or_blue = r1
    red_at_green_in_red = r2
    blue_at_green_in_blue = r2
    red_at_green_in_blue = r3
    blue_at_green_in_red = r3
    red_at_blue = r4
    blue_at_red = r4

    r_mod_2 = iand(g_r, 1)
    c_mod_2 = iand(g_c, 1)
    red_col = merge(1, 0, bayer_pattern == 1 .or. bayer_pattern == 3)
    red_row = merge(1, 0, bayer_pattern == 2 .or. bayer_pattern == 3)
    blue_col = 1 - red_col
    blue_row = 1 - red_row

    in_red_row = r_mod_2 == red_row
    in_blue_row = r_mod_2 == blue_row
    is_red_pixel = (r_mod_2 == red_row) .and. (c_mod_2 == red_col)
    is_blue_pixel = (r_mod_2 == blue_row) .and. (c_mod_2 == blue_col)
    is_green_pixel = .not. (is_red_pixel .or. is_blue_pixel)

    red_value = bool_int(is_red_pixel) * f_ij + bool_int(is_blue_pixel) * red_at_blue + &
      bool_int(is_green_pixel .and. in_red_row) * red_at_green_in_red + &
      bool_int(is_green_pixel .and. in_blue_row) * red_at_green_in_blue
    blue_value = bool_int(is_blue_pixel) * f_ij + bool_int(is_red_pixel) * blue_at_red + &
      bool_int(is_green_pixel .and. in_red_row) * blue_at_green_in_red + &
      bool_int(is_green_pixel .and. in_blue_row) * blue_at_green_in_blue
    green_value = bool_int(is_green_pixel) * f_ij + bool_int(.not. is_green_pixel) * green_at_red_or_blue

    out_idx = g_r * output_pitch + g_c * 4 + 1
    output_image(out_idx) = saturate_uchar(red_value)
    output_image(out_idx + 1) = saturate_uchar(green_value)
    output_image(out_idx + 2) = saturate_uchar(blue_value)
    output_image(out_idx + 3) = 0
  end subroutine demosaic_pixel

  integer function f(row, col, input_image, height, width, pitch) result(value)
    integer, intent(in) :: row, col, height, width, pitch
    integer, intent(in) :: input_image(:)
    value = sample_pixel(input_image, height, width, pitch, row, col)
  end function f

  integer function sample_pixel(input_image, height, width, pitch, row, col) result(value)
    integer, intent(in) :: input_image(:), height, width, pitch, row, col
    integer :: rr, cc

    rr = reflect_exclusive(row, height)
    cc = reflect_exclusive(col, width)
    value = input_image(rr * pitch + cc + 1)
  end function sample_pixel

  integer function reflect_exclusive(coord, limit) result(value)
    integer, intent(in) :: coord, limit

    value = coord
    if (value < 0) then
      value = -value
    end if
    if (value >= limit) then
      value = limit - (value - limit) - 2
    end if
  end function reflect_exclusive

  integer function saturate_uchar(value) result(out)
    integer, intent(in) :: value

    if (value > 255) then
      out = 255
    else if (value < 0) then
      out = 0
    else
      out = value
    end if
  end function saturate_uchar

  integer function bool_int(value) result(out)
    logical, intent(in) :: value
    out = merge(1, 0, value)
  end function bool_int

end program main
