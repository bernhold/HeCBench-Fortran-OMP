program main
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
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

  type :: properties_t
    real(real32) :: k
    real(real32) :: center_x
    real(real32) :: center_y
    integer :: width
    integer :: height
    real(real32) :: thresh
    real(real32) :: xscale
    real(real32) :: yscale
    real(real32) :: xshift
    real(real32) :: yshift
  end type properties_t

  character(len=256) :: arg0, arg
  integer :: width, height, repeat_count, image_size
  real(real32) :: coeff, new_center_x, new_center_y, xshift_2, yshift_2
  type(properties_t) :: prop
  integer(int32), allocatable :: src_r(:), src_g(:), src_b(:)
  integer(int32), allocatable :: dst_r(:), dst_g(:), dst_b(:)
  integer(int32), allocatable :: ref_r(:), ref_g(:), ref_b(:)
  integer :: ex, ey, ez

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), &
      '<input image width> <input image height> <coefficient of distortion> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) width
  call get_command_argument(2, arg); read(arg, *) height
  call get_command_argument(3, arg); read(arg, *) coeff
  call get_command_argument(4, arg); read(arg, *) repeat_count
  if (width <= 0 .or. height <= 0 .or. repeat_count <= 0) stop 1

  prop%k = coeff
  prop%center_x = real(width / 2, real32)
  prop%center_y = real(height / 2, real32)
  prop%width = width
  prop%height = height
  prop%thresh = 1.0_real32
  prop%xshift = calc_shift(0.0_real32, prop%center_x - 1.0_real32, prop%center_x, prop%k, prop%thresh)
  new_center_x = real(prop%width, real32) - prop%center_x
  xshift_2 = calc_shift(0.0_real32, new_center_x - 1.0_real32, new_center_x, prop%k, prop%thresh)
  prop%yshift = calc_shift(0.0_real32, prop%center_y - 1.0_real32, prop%center_y, prop%k, prop%thresh)
  new_center_y = real(prop%height, real32) - prop%center_y
  yshift_2 = calc_shift(0.0_real32, new_center_y - 1.0_real32, new_center_y, prop%k, prop%thresh)
  prop%xscale = (real(prop%width, real32) - prop%xshift - xshift_2) / real(prop%width, real32)
  prop%yscale = (real(prop%height, real32) - prop%yshift - yshift_2) / real(prop%height, real32)

  image_size = width * height
  allocate(src_r(image_size), src_g(image_size), src_b(image_size))
  allocate(dst_r(image_size), dst_g(image_size), dst_b(image_size))
  allocate(ref_r(image_size), ref_g(image_size), ref_b(image_size))

  call c_srand(123_c_int)
  call fill_image(src_r, src_g, src_b)

  call run_distort(src_r, src_g, src_b, dst_r, dst_g, dst_b, prop, repeat_count)
  call reference(src_r, src_g, src_b, ref_r, ref_g, ref_b, prop)
  call max_error(dst_r, dst_g, dst_b, ref_r, ref_g, ref_b, ex, ey, ez)
  write(*,'(A,I0,1X,I0,1X,I0)') 'Max error of each channel: ', ex, ey, ez

  deallocate(src_r, src_g, src_b, dst_r, dst_g, dst_b, ref_r, ref_g, ref_b)

contains

  recursive real(real32) function calc_shift(x1, x2, cx, k, thresh) result(shift)
    real(real32), intent(in) :: x1, x2, cx, k, thresh
    real(real32) :: x3, result1, result3

    x3 = x1 + (x2 - x1) * 0.5_real32
    result1 = x1 + ((x1 - cx) * k * ((x1 - cx) * (x1 - cx)))
    result3 = x3 + ((x3 - cx) * k * ((x3 - cx) * (x3 - cx)))
    if (result1 > -thresh .and. result1 < thresh) then
      shift = x1
    else if (result3 < 0.0_real32) then
      shift = calc_shift(x3, x2, cx, k, thresh)
    else
      shift = calc_shift(x1, x3, cx, k, thresh)
    end if
  end function calc_shift

  subroutine fill_image(src_r, src_g, src_b)
    integer(int32), intent(out) :: src_r(:), src_g(:), src_b(:)
    integer :: i

    do i = 1, size(src_r)
      src_r(i) = next_rand_mod(256)
      src_g(i) = next_rand_mod(256)
      src_b(i) = next_rand_mod(256)
    end do
  end subroutine fill_image

  integer function next_rand_mod(divisor)
    integer, intent(in) :: divisor
    integer(c_int) :: value

    value = c_rand()
    next_rand_mod = modulo(value, divisor)
  end function next_rand_mod

  subroutine run_distort(src_r, src_g, src_b, dst_r, dst_g, dst_b, prop, repeat_count)
    integer(int32), intent(in) :: src_r(:), src_g(:), src_b(:)
    integer(int32), intent(out) :: dst_r(:), dst_g(:), dst_b(:)
    type(properties_t), intent(in) :: prop
    integer, intent(in) :: repeat_count
    integer :: iter, image_size
    real(real64) :: start_time, elapsed_ms

    image_size = prop%width * prop%height
    !$omp target data map(to: src_r(1:image_size), src_g(1:image_size), src_b(1:image_size), prop) &
    !$omp& map(from: dst_r(1:image_size), dst_g(1:image_size), dst_b(1:image_size))
    start_time = omp_get_wtime()
    do iter = 1, repeat_count
      call barrel_distort(src_r, src_g, src_b, dst_r, dst_g, dst_b, prop)
    end do
    elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64 / real(repeat_count, real64)
    write(*,'(A,F0.6,A)') 'Average kernel execution time: ', elapsed_ms, ' (ms)'
    !$omp end target data
  end subroutine run_distort

  subroutine barrel_distort(src_r, src_g, src_b, dst_r, dst_g, dst_b, prop)
    integer(int32), intent(in) :: src_r(:), src_g(:), src_b(:)
    integer(int32), intent(out) :: dst_r(:), dst_g(:), dst_b(:)
    type(properties_t), intent(in) :: prop
    integer :: row, col, idx
    real(real32) :: radial_x, radial_y
    integer(int32) :: rr, gg, bb

    !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(idx, radial_x, radial_y, rr, gg, bb)
    do row = 0, prop%height - 1
      do col = 0, prop%width - 1
        radial_x = get_radial_x(real(col, real32), real(row, real32), prop)
        radial_y = get_radial_y(real(col, real32), real(row, real32), prop)
        call sample_image(src_r, src_g, src_b, radial_y, radial_x, rr, gg, bb, prop)
        idx = row * prop%width + col + 1
        dst_r(idx) = rr
        dst_g(idx) = gg
        dst_b(idx) = bb
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine barrel_distort

  subroutine reference(src_r, src_g, src_b, dst_r, dst_g, dst_b, prop)
    integer(int32), intent(in) :: src_r(:), src_g(:), src_b(:)
    integer(int32), intent(out) :: dst_r(:), dst_g(:), dst_b(:)
    type(properties_t), intent(in) :: prop
    integer :: row, col, idx
    real(real32) :: radial_x, radial_y
    integer(int32) :: rr, gg, bb

    do row = 0, prop%height - 1
      do col = 0, prop%width - 1
        radial_x = get_radial_x(real(col, real32), real(row, real32), prop)
        radial_y = get_radial_y(real(col, real32), real(row, real32), prop)
        call sample_image(src_r, src_g, src_b, radial_y, radial_x, rr, gg, bb, prop)
        idx = row * prop%width + col + 1
        dst_r(idx) = rr
        dst_g(idx) = gg
        dst_b(idx) = bb
      end do
    end do
  end subroutine reference

  real(real32) function get_radial_x(x, y, prop)
    real(real32), intent(in) :: x, y
    type(properties_t), intent(in) :: prop
    real(real32) :: scaled_x, scaled_y

    scaled_x = x * prop%xscale + prop%xshift
    scaled_y = y * prop%yscale + prop%yshift
    get_radial_x = scaled_x + ((scaled_x - prop%center_x) * prop%k * &
      ((scaled_x - prop%center_x) * (scaled_x - prop%center_x) + &
       (scaled_y - prop%center_y) * (scaled_y - prop%center_y)))
  end function get_radial_x

  real(real32) function get_radial_y(x, y, prop)
    real(real32), intent(in) :: x, y
    type(properties_t), intent(in) :: prop
    real(real32) :: scaled_x, scaled_y

    scaled_x = x * prop%xscale + prop%xshift
    scaled_y = y * prop%yscale + prop%yshift
    get_radial_y = scaled_y + ((scaled_y - prop%center_y) * prop%k * &
      ((scaled_x - prop%center_x) * (scaled_x - prop%center_x) + &
       (scaled_y - prop%center_y) * (scaled_y - prop%center_y)))
  end function get_radial_y

  subroutine sample_image(src_r, src_g, src_b, idx0, idx1, rr, gg, bb, prop)
    integer(int32), intent(in) :: src_r(:), src_g(:), src_b(:)
    real(real32), intent(in) :: idx0, idx1
    integer(int32), intent(out) :: rr, gg, bb
    type(properties_t), intent(in) :: prop
    integer :: idx0_floor, idx0_ceil, idx1_floor, idx1_ceil
    integer :: i1, i2, i3, i4
    real(real32) :: x, y, r_value, g_value, b_value

    if (idx0 < 0.0_real32 .or. idx1 < 0.0_real32 .or. &
        idx0 > real(prop%height - 1, real32) .or. idx1 > real(prop%width - 1, real32)) then
      rr = 0
      gg = 0
      bb = 0
      return
    end if

    idx0_floor = floor_int(idx0)
    idx0_ceil = ceiling_int(idx0)
    idx1_floor = floor_int(idx1)
    idx1_ceil = ceiling_int(idx1)
    i1 = idx0_floor * prop%width + idx1_floor + 1
    i2 = idx0_floor * prop%width + idx1_ceil + 1
    i3 = idx0_ceil * prop%width + idx1_ceil + 1
    i4 = idx0_ceil * prop%width + idx1_floor + 1
    x = idx0 - real(idx0_floor, real32)
    y = idx1 - real(idx1_floor, real32)

    r_value = real(src_r(i1), real32) * (1.0_real32 - x) * (1.0_real32 - y) + &
      real(src_r(i2), real32) * (1.0_real32 - x) * y + real(src_r(i3), real32) * x * y + &
      real(src_r(i4), real32) * x * (1.0_real32 - y)
    g_value = real(src_g(i1), real32) * (1.0_real32 - x) * (1.0_real32 - y) + &
      real(src_g(i2), real32) * (1.0_real32 - x) * y + real(src_g(i3), real32) * x * y + &
      real(src_g(i4), real32) * x * (1.0_real32 - y)
    b_value = real(src_b(i1), real32) * (1.0_real32 - x) * (1.0_real32 - y) + &
      real(src_b(i2), real32) * (1.0_real32 - x) * y + real(src_b(i3), real32) * x * y + &
      real(src_b(i4), real32) * x * (1.0_real32 - y)
    rr = int(r_value, int32)
    gg = int(g_value, int32)
    bb = int(b_value, int32)
  end subroutine sample_image

  integer function floor_int(value)
    real(real32), intent(in) :: value

    floor_int = int(floor(value))
  end function floor_int

  integer function ceiling_int(value)
    real(real32), intent(in) :: value

    ceiling_int = int(ceiling(value))
  end function ceiling_int

  subroutine max_error(dst_r, dst_g, dst_b, ref_r, ref_g, ref_b, ex, ey, ez)
    integer(int32), intent(in) :: dst_r(:), dst_g(:), dst_b(:), ref_r(:), ref_g(:), ref_b(:)
    integer, intent(out) :: ex, ey, ez
    integer :: i

    ex = 0
    ey = 0
    ez = 0
    do i = 1, size(dst_r)
      ex = max(abs(dst_r(i) - ref_r(i)), ex)
      ey = max(abs(dst_g(i) - ref_g(i)), ey)
      ez = max(abs(dst_b(i) - ref_b(i)), ez)
    end do
  end subroutine max_error

end program main
