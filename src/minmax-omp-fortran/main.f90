program minmax
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  integer :: repeat, n
  real(real32) :: box_size
  character(len=64) :: arg
  real(real32), allocatable :: points_x(:), points_y(:)
  real(real32) :: min_value(2), max_value(2), ref_min, ref_max
  real(real64) :: start_time, elapsed_us
  logical :: ok

  if (command_argument_count() /= 2) then
    write(*,'("Usage: ./main <bounding-box size> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) box_size
  call get_command_argument(2, arg)
  read(arg, *) repeat

  call generate_points(box_size, points_x, points_y, n)
  write(*,'("Total number of points: ",I0)') n

  !$omp target data map(to: points_x, points_y)
  start_time = omp_get_wtime()
  call run_separate(points_x, points_y, n, repeat, min_value(1), max_value(1))
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'("Average execution time of omp:min() + omp:max(): ",F0.6," (us)")') elapsed_us

  start_time = omp_get_wtime()
  call run_combined(points_x, points_y, n, repeat, min_value(2), max_value(2))
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'("Average execution time of omp:minmax(): ",F0.6," (us)")') elapsed_us
  !$omp end target data

  call cpu_reference(points_x, points_y, n, ref_min, ref_max)
  ok = abs(min_value(1) - ref_min) <= 1.0e-5_real32 .and. &
       abs(max_value(1) - ref_max) <= 1.0e-5_real32 .and. &
       abs(min_value(2) - ref_min) <= 1.0e-5_real32 .and. &
       abs(max_value(2) - ref_max) <= 1.0e-5_real32

  if (ok) then
    write(*,'("PASS")')
  else
    write(*,'("FAIL")')
    stop 1
  end if

contains

  subroutine run_separate(px, py, count, repeat, min_out, max_out)
    real(real32), intent(in) :: px(:), py(:)
    integer, intent(in) :: count, repeat
    real(real32), intent(out) :: min_out, max_out
    integer :: r, i
    real(real32) :: min_val, max_val, val

    do r = 1, repeat
      min_val = huge(1.0_real32)
      !$omp target teams distribute parallel do reduction(min:min_val)
      do i = 1, count
        val = px(i) * px(i) + py(i) * py(i)
        min_val = min(min_val, val)
      end do
      !$omp end target teams distribute parallel do

      max_val = -huge(1.0_real32)
      !$omp target teams distribute parallel do reduction(max:max_val)
      do i = 1, count
        val = px(i) * px(i) + py(i) * py(i)
        max_val = max(max_val, val)
      end do
      !$omp end target teams distribute parallel do
    end do

    min_out = min_val
    max_out = max_val
  end subroutine run_separate

  subroutine run_combined(px, py, count, repeat, min_out, max_out)
    real(real32), intent(in) :: px(:), py(:)
    integer, intent(in) :: count, repeat
    real(real32), intent(out) :: min_out, max_out
    integer :: r, i
    real(real32) :: min_val, max_val, val

    do r = 1, repeat
      min_val = huge(1.0_real32)
      max_val = -huge(1.0_real32)
      !$omp target teams distribute parallel do reduction(min:min_val) reduction(max:max_val)
      do i = 1, count
        val = px(i) * px(i) + py(i) * py(i)
        min_val = min(min_val, val)
        max_val = max(max_val, val)
      end do
      !$omp end target teams distribute parallel do
    end do

    min_out = min_val
    max_out = max_val
  end subroutine run_combined

  subroutine cpu_reference(px, py, count, min_out, max_out)
    real(real32), intent(in) :: px(:), py(:)
    integer, intent(in) :: count
    real(real32), intent(out) :: min_out, max_out
    integer :: i
    real(real32) :: val

    min_out = px(1) * px(1) + py(1) * py(1)
    max_out = min_out
    do i = 2, count
      val = px(i) * px(i) + py(i) * py(i)
      min_out = min(min_out, val)
      max_out = max(max_out, val)
    end do
  end subroutine cpu_reference

  subroutine generate_points(size, px, py, total_points)
    real(real32), intent(in) :: size
    real(real32), allocatable, intent(out) :: px(:), py(:)
    integer, intent(out) :: total_points
    integer, allocatable :: rect_counts(:)
    real(real32), allocatable :: tlx(:), tly(:), brx(:), bry(:)
    real(real32) :: cur_tlx, cur_tly, cur_brx, cur_bry, area_x, area_y, phi
    integer :: rect_capacity, rects, nrect_points, offset, i, j
    integer(int64) :: seed

    phi = (1.0_real32 + sqrt(5.0_real32)) * 0.5_real32
    cur_tlx = 0.0_real32
    cur_tly = 0.0_real32
    cur_brx = size
    cur_bry = size
    area_x = cur_brx - cur_tlx
    area_y = cur_bry - cur_tly

    rect_capacity = 64
    allocate(rect_counts(rect_capacity), tlx(rect_capacity), tly(rect_capacity), brx(rect_capacity), bry(rect_capacity))
    rects = 0
    total_points = 0

    do while (area_x > 1.0_real32 .and. area_y > 1.0_real32)
      select case (mod(rects, 4))
      case (0)
        cur_brx = cur_tlx - (cur_tlx - cur_brx) / phi
      case (1)
        cur_bry = cur_tly - (cur_tly - cur_bry) / phi
      case (2)
        cur_tlx = cur_tlx + (cur_brx - cur_tlx) / phi
      case default
        cur_tly = cur_tly + (cur_bry - cur_tly) / phi
      end select

      area_x = cur_brx - cur_tlx
      area_y = cur_bry - cur_tly
      nrect_points = int(sqrt(area_x * area_y * 1000000.0_real32))

      rects = rects + 1
      if (rects > rect_capacity) stop 1
      rect_counts(rects) = nrect_points
      tlx(rects) = cur_tlx
      tly(rects) = cur_tly
      brx(rects) = cur_brx
      bry(rects) = cur_bry
      total_points = total_points + nrect_points
    end do

    allocate(px(total_points), py(total_points))
    seed = 123_int64
    offset = 0
    do i = 1, rects
      do j = 1, rect_counts(i)
        offset = offset + 1
        px(offset) = tlx(i) + (brx(i) - tlx(i)) * next_unit(seed)
        py(offset) = tly(i) + (bry(i) - tly(i)) * next_unit(seed)
      end do
    end do
  end subroutine generate_points

  real(real32) function next_unit(seed)
    integer(int64), intent(inout) :: seed

    seed = mod(seed * 1103515245_int64 + 12345_int64, 2147483648_int64)
    next_unit = real(seed, real32) / 2147483647.0_real32
  end function next_unit

end program minmax
