program main
  use, intrinsic :: iso_fortran_env, only : int8, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: max_detections = 4096
  integer, parameter :: n_partitions = 32

  character(len=512) :: input_file, output_file
  integer :: repeat, ndetections, limit, threads, totaldets
  real(real32), allocatable :: px(:), py(:), pz(:), ps(:)
  integer(int8), allocatable :: pointsbitmap(:), nmsbitmap(:)
  real(real64) :: start_time, end_time

  if (command_argument_count() /= 3) then
    call print_help()
    stop
  end if

  call get_command_argument(1, input_file)
  call get_command_argument(2, output_file)
  repeat = read_arg(3)

  allocate(px(max_detections), py(max_detections), pz(max_detections), ps(max_detections))
  allocate(pointsbitmap(max_detections), nmsbitmap(max_detections * max_detections))
  px = 0.0_real32
  py = 0.0_real32
  pz = 0.0_real32
  ps = 0.0_real32
  pointsbitmap = 0_int8
  nmsbitmap = 1_int8

  call read_points(trim(input_file), px, py, pz, ps, ndetections)
  print '(A,A,A,I0)', 'Number of detections read from input file (', trim(input_file), '): ', ndetections

  limit = get_upper_limit(ndetections, 16)
  threads = get_optimal_dim(limit) * get_optimal_dim(limit)

  !$omp target data map(to: px(1:max_detections), py(1:max_detections), pz(1:max_detections), ps(1:max_detections), &
  !$omp& nmsbitmap(1:max_detections * max_detections)) map(tofrom: pointsbitmap(1:max_detections))
  start_time = omp_get_wtime()
  call generate_nms_bitmap(px, py, pz, ps, nmsbitmap, limit, repeat, threads)
  end_time = omp_get_wtime()
  print '(A,F0.6,A)', 'Average kernel execution time (generate_nms_bitmap): ', &
    (end_time - start_time) / real(repeat, real64), ' (s)'

  start_time = omp_get_wtime()
  call reduce_nms_bitmap(nmsbitmap, pointsbitmap, ndetections, repeat)
  end_time = omp_get_wtime()
  print '(A,F0.6,A)', 'Average kernel execution time (reduce_nms_bitmap): ', &
    (end_time - start_time) / real(repeat, real64), ' (s)'
  !$omp end target data

  call write_points(trim(output_file), px, py, pz, ps, pointsbitmap, ndetections, totaldets)
  print '(A,I0)', 'Detections after NMS: ', totaldets

  deallocate(px, py, pz, ps, pointsbitmap, nmsbitmap)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine print_help()
    print '(A)', ''
    print '(A)', 'Usage: nmstest  <detections.txt>  <output.txt>'
    print '(A)', ''
    print '(A)', '               detections.txt -> Input file containing the coordinates, width, and scores of detected objects'
    print '(A)', '               output.txt     -> Output file after performing NMS'
    print '(A)', '               repeat         -> Kernel execution count'
    print '(A)', ''
  end subroutine print_help

  integer function get_optimal_dim(val)
    integer, intent(in) :: val
    integer :: div, neg, cntneg, cntpos, i
    neg = 1
    div = 16
    cntneg = div
    cntpos = div
    do i = 1, 5
      if (mod(val, div) == 0) then
        get_optimal_dim = div
        return
      end if
      if (neg /= 0) then
        cntneg = cntneg - 1
        div = cntneg
        neg = 0
      else
        cntpos = cntpos + 1
        div = cntpos
        neg = 1
      end if
    end do
    get_optimal_dim = 16
  end function get_optimal_dim

  integer function get_upper_limit(val, mul)
    integer, intent(in) :: val, mul
    integer :: cnt
    cnt = mul
    do while (cnt < val)
      cnt = cnt + mul
    end do
    if (cnt > max_detections) cnt = max_detections
    get_upper_limit = cnt
  end function get_upper_limit

  subroutine read_points(path, px, py, pz, ps, ndetections)
    character(len=*), intent(in) :: path
    real(real32), intent(inout) :: px(:), py(:), pz(:), ps(:)
    integer, intent(out) :: ndetections
    character(len=256) :: line
    integer :: unit, ios, k, x, y, w
    real(real32) :: score

    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print '(A,A,A)', 'Error: Unable to open file ', path, ' for input detection coordinates.'
      stop 1
    end if

    ndetections = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      do k = 1, len_trim(line)
        if (line(k:k) == ',') line(k:k) = ' '
      end do
      read(line, *, iostat=ios) x, y, w, score
      if (ios /= 0) then
        print '(A,I0,A,A)', 'Error: Invalid file format in line ', ndetections, ' when reading ', path
        stop 1
      end if
      ndetections = ndetections + 1
      if (ndetections > max_detections) exit
      px(ndetections) = real(x, real32)
      py(ndetections) = real(y, real32)
      pz(ndetections) = real(w, real32)
      ps(ndetections) = score
    end do
    close(unit)
  end subroutine read_points

  subroutine generate_nms_bitmap(px, py, pz, ps, nmsbitmap, limit, repeat, threads)
    real(real32), intent(in) :: px(:), py(:), pz(:), ps(:)
    integer(int8), intent(inout) :: nmsbitmap(:)
    integer, intent(in) :: limit, repeat, threads
    integer :: rep, i, j, idx
    real(real32) :: area, overlap_w, overlap_h

    do rep = 1, repeat
      !$omp target teams distribute parallel do collapse(2) thread_limit(threads) private(idx, area, overlap_w, overlap_h)
      do i = 1, limit
        do j = 1, limit
          if (ps(i) < ps(j)) then
            area = (pz(j) + 1.0_real32) * (pz(j) + 1.0_real32)
            overlap_w = max(0.0_real32, min(px(i) + pz(i), px(j) + pz(j)) - max(px(i), px(j)) + 1.0_real32)
            overlap_h = max(0.0_real32, min(py(i) + pz(i), py(j) + pz(j)) - max(py(i), py(j)) + 1.0_real32)
            idx = (i - 1) * max_detections + j
            if (((overlap_w * overlap_h) / area) < 0.3_real32 .and. pz(j) /= 0.0_real32) then
              nmsbitmap(idx) = 1_int8
            else
              nmsbitmap(idx) = 0_int8
            end if
          end if
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
  end subroutine generate_nms_bitmap

  subroutine reduce_nms_bitmap(nmsbitmap, pointsbitmap, ndetections, repeat)
    integer(int8), intent(in) :: nmsbitmap(:)
    integer(int8), intent(inout) :: pointsbitmap(:)
    integer, intent(in) :: ndetections, repeat
    integer :: rep, bid, j, accum

    do rep = 1, repeat
      !$omp target teams distribute parallel do thread_limit(max_detections / n_partitions) private(j, accum)
      do bid = 1, ndetections
        accum = 1
        do j = 1, max_detections
          accum = iand(accum, int(nmsbitmap((bid - 1) * max_detections + j)))
        end do
        pointsbitmap(bid) = int(accum, int8)
      end do
      !$omp end target teams distribute parallel do
    end do
  end subroutine reduce_nms_bitmap

  subroutine write_points(path, px, py, pz, ps, pointsbitmap, ndetections, totaldets)
    character(len=*), intent(in) :: path
    real(real32), intent(in) :: px(:), py(:), pz(:), ps(:)
    integer(int8), intent(in) :: pointsbitmap(:)
    integer, intent(in) :: ndetections
    integer, intent(out) :: totaldets
    integer :: unit, ios, i

    open(newunit=unit, file=path, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      print '(A,A,A)', 'Error: Unable to open file ', path, ' for detection outcome.'
      stop 1
    end if
    totaldets = 0
    do i = 1, ndetections
      if (pointsbitmap(i) /= 0_int8) then
        write(unit,'(I0,A,I0,A,I0,A,F0.6)') int(px(i)), ',', int(py(i)), ',', int(pz(i)), ',', ps(i)
        totaldets = totaldets + 1
      end if
    end do
    close(unit)
  end subroutine write_points

end program main
