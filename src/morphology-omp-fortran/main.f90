program main
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: black = 0_int32
  integer(int32), parameter :: white = 255_int32
  character(len=256) :: arg0, arg1, arg2, arg3, arg4, arg5
  integer(int32) :: hsize, vsize, width, height, repeat
  integer(int32), allocatable :: src_img(:), tmp_img(:)
  integer(int32) :: i, j, n, total
  real(real64) :: dilate_time, erode_time

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 5) then
    print '(3A)', 'Usage: ', trim(arg0), ' <kernel width> <kernel height> ', &
      '<image width> <image height> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  call get_command_argument(4, arg4)
  call get_command_argument(5, arg5)
  read(arg1, *) hsize
  read(arg2, *) vsize
  read(arg3, *) width
  read(arg4, *) height
  read(arg5, *) repeat

  if (hsize <= 0_int32 .or. vsize <= 0_int32 .or. width <= 0_int32 .or. height <= 0_int32) stop 1
  if (repeat <= 0_int32) stop 1

  allocate(src_img(width * height), tmp_img(width * height))
  src_img = black
  tmp_img = black

  do i = 1, height
    do j = 1, width
      if (i == height / 2 .and. j == width / 2) then
        src_img(index_1d(i, j, width)) = white
      end if
    end do
  end do

  dilate_time = 0.0_real64
  erode_time = 0.0_real64

  !$omp target data map(tofrom: src_img(1:width * height)) map(alloc: tmp_img(1:width * height))
  do n = 1, repeat
    dilate_time = dilate_time + morphology_pass(src_img, tmp_img, width, height, hsize, vsize, .true.)
    erode_time = erode_time + morphology_pass(src_img, tmp_img, width, height, hsize, vsize, .false.)
  end do
  !$omp end target data

  print '(A,F0.6,A)', 'Average kernel execution time (dilate): ', dilate_time / real(repeat, real64), ' (s)'
  print '(A,F0.6,A)', 'Average kernel execution time (erode): ', erode_time / real(repeat, real64), ' (s)'

  total = sum(src_img)
  if (total == white) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(src_img, tmp_img)

contains

  integer(int32) function index_1d(row, col, width_value) result(idx)
    integer(int32), intent(in) :: row, col, width_value

    idx = (row - 1_int32) * width_value + col
  end function index_1d

  real(real64) function morphology_pass(img, tmp, width, height, hsize, vsize, is_dilate) result(elapsed)
    integer(int32), intent(inout) :: img(:), tmp(:)
    integer(int32), intent(in) :: width, height, hsize, vsize
    logical, intent(in) :: is_dilate
    real(real64) :: start_time, end_time

    call clear_image(tmp, width, height)
    start_time = omp_get_wtime()
    call horizontal_pass(img, tmp, width, height, hsize, is_dilate)
    call vertical_pass(tmp, img, width, height, vsize, is_dilate)
    end_time = omp_get_wtime()
    elapsed = end_time - start_time
  end function morphology_pass

  subroutine clear_image(img, width, height)
    integer(int32), intent(inout) :: img(:)
    integer(int32), intent(in) :: width, height
    integer(int32) :: p

    !$omp target teams distribute parallel do thread_limit(256)
    do p = 1, width * height
      img(p) = black
    end do
    !$omp end target teams distribute parallel do
  end subroutine clear_image

  subroutine horizontal_pass(src, dst, width, height, hsize, is_dilate)
    integer(int32), intent(in) :: src(:)
    integer(int32), intent(inout) :: dst(:)
    integer(int32), intent(in) :: width, height, hsize
    logical, intent(in) :: is_dilate
    integer(int32) :: row, col, k, half, value, sample

    half = hsize / 2_int32
    !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(k, value, sample)
    do row = 1, height
      do col = 1, width
        if (col <= half .or. col > width - half) cycle
        if (is_dilate) then
          value = black
          do k = -half, half
            sample = src(index_1d(row, col + k, width))
            if (sample > value) value = sample
          end do
        else
          value = white
          do k = -half, half
            sample = src(index_1d(row, col + k, width))
            if (sample < value) value = sample
          end do
        end if
        dst(index_1d(row, col, width)) = value
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine horizontal_pass

  subroutine vertical_pass(src, dst, width, height, vsize, is_dilate)
    integer(int32), intent(in) :: src(:)
    integer(int32), intent(inout) :: dst(:)
    integer(int32), intent(in) :: width, height, vsize
    logical, intent(in) :: is_dilate
    integer(int32) :: row, col, k, half, value, sample, border_value

    half = vsize / 2_int32
    border_value = merge(white, black, is_dilate)
    !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(k, value, sample)
    do row = 1, height
      do col = 1, width
        if (row <= half .or. row > height - half) then
          dst(index_1d(row, col, width)) = border_value
        else
          if (is_dilate) then
            value = black
            do k = -half, half
              sample = src(index_1d(row + k, col, width))
              if (sample > value) value = sample
            end do
          else
            value = white
            do k = -half, half
              sample = src(index_1d(row + k, col, width))
              if (sample < value) value = sample
            end do
          end if
          dst(index_1d(row, col, width)) = value
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine vertical_pass

end program main
