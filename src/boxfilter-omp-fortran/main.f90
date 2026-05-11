program main
  use, intrinsic :: iso_fortran_env, only : int8, int32, int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: radius = 10
  real(real32), parameter :: scale = 1.0_real32 / real(2 * radius + 1, real32)

  character(len=512) :: image_path
  integer :: repeat, width, height, n_pixels, i, error
  integer(int32), allocatable :: input(:), tmp(:), dev_output(:), host_output(:)

  if (command_argument_count() /= 2) then
    print '(A)', 'Usage ./main <PPM image> <repeat>'
    stop 1
  end if

  call get_command_argument(1, image_path)
  repeat = read_arg(2)
  if (repeat <= 0) error stop 'repeat must be positive'

  call load_ppm(trim(image_path), input, width, height)
  n_pixels = width * height
  allocate(tmp(n_pixels), dev_output(n_pixels), host_output(n_pixels))
  tmp = 0_int32
  dev_output = 0_int32
  host_output = 0_int32

  print '(A,I0,A,I0,A,I0,A,I0)', 'Image Width = ', width, ', Height = ', height, &
      ', bpp = ', 32, ', Mask Radius = ', radius
  print '(A)', 'Using Local Memory for Row Processing'
  print '(A)'

  !$omp target data map(to: input(1:n_pixels)) map(tofrom: tmp(1:n_pixels), dev_output(1:n_pixels))
  print '(A)', 'Warmup..'
  call box_filter_device(input, tmp, dev_output, width, height, repeat)
  print '(A)'
  print '(A,I0,A)', 'Running BoxFilterGPU for ', repeat, ' cycles...'
  print '(A)'
  call box_filter_device(input, tmp, dev_output, width, height, repeat)
  !$omp end target data

  call box_filter_host(input, tmp, host_output, width, height)

  error = 0
  do i = radius * width + 1, n_pixels - radius * width
    if (dev_output(i) /= host_output(i)) then
      write(*, '(I0,1X,Z8.8,1X,Z8.8)') i - 1, dev_output(i), host_output(i)
      error = 1
      exit
    end if
  end do
  if (error == 0) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(input, tmp, dev_output, host_output)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine load_ppm(path, image, width, height)
    character(len=*), intent(in) :: path
    integer(int32), allocatable, intent(out) :: image(:)
    integer, intent(out) :: width, height
    integer :: unit, file_size, pos, maxval, pixel, base
    character(len=1), allocatable :: bytes(:)
    character(len=64) :: token

    inquire(file=path, size=file_size)
    if (file_size <= 0) error stop 'invalid PPM input'
    allocate(bytes(file_size))
    open(newunit=unit, file=path, access='stream', form='unformatted', status='old', action='read')
    read(unit) bytes
    close(unit)

    pos = 1
    call next_token(bytes, file_size, pos, token)
    if (trim(token) /= 'P6') error stop 'expected binary P6 PPM'
    call next_token(bytes, file_size, pos, token)
    read(token, *) width
    call next_token(bytes, file_size, pos, token)
    read(token, *) height
    call next_token(bytes, file_size, pos, token)
    read(token, *) maxval
    if (maxval /= 255) error stop 'expected 8-bit PPM'

    allocate(image(width * height))
    do pixel = 1, width * height
      base = pos + (pixel - 1) * 3
      image(pixel) = pack_rgba(byte_value(bytes(base)), byte_value(bytes(base + 1)), byte_value(bytes(base + 2)), 0)
    end do
    deallocate(bytes)
  end subroutine load_ppm

  subroutine next_token(bytes, file_size, pos, token)
    character(len=1), intent(in) :: bytes(:)
    integer, intent(in) :: file_size
    integer, intent(inout) :: pos
    character(len=*), intent(out) :: token
    integer :: out_pos, ch
    token = ''
    do while (pos <= file_size)
      ch = iachar(bytes(pos))
      if (ch == iachar('#')) then
        do while (pos <= file_size .and. iachar(bytes(pos)) /= 10)
          pos = pos + 1
        end do
      else if (.not. is_space(ch)) then
        exit
      end if
      pos = pos + 1
    end do
    out_pos = 1
    do while (pos <= file_size)
      ch = iachar(bytes(pos))
      if (is_space(ch) .or. ch == iachar('#')) exit
      if (out_pos <= len(token)) token(out_pos:out_pos) = bytes(pos)
      out_pos = out_pos + 1
      pos = pos + 1
    end do
    do while (pos <= file_size .and. is_space(iachar(bytes(pos))))
      pos = pos + 1
    end do
  end subroutine next_token

  logical function is_space(ch)
    integer, intent(in) :: ch
    is_space = ch == 9 .or. ch == 10 .or. ch == 13 .or. ch == 32
  end function is_space

  integer function byte_value(ch)
    character(len=1), intent(in) :: ch
    byte_value = iachar(ch)
  end function byte_value

  integer(int32) function pack_rgba(r, g, b, a)
    integer, intent(in) :: r, g, b, a
    pack_rgba = int(ior(ior(ior(iand(r, 255), ishft(iand(g, 255), 8)), &
        ishft(iand(b, 255), 16)), ishft(iand(a, 255), 24)), int32)
  end function pack_rgba

  subroutine unpack_rgba(pixel, r, g, b, a)
    integer(int32), intent(in) :: pixel
    real(real32), intent(out) :: r, g, b, a
    r = real(iand(pixel, int(z'000000ff', int32)), real32)
    g = real(iand(ishft(pixel, -8), int(z'000000ff', int32)), real32)
    b = real(iand(ishft(pixel, -16), int(z'000000ff', int32)), real32)
    a = real(iand(ishft(pixel, -24), int(z'000000ff', int32)), real32)
  end subroutine unpack_rgba

  integer(int32) function pack_scaled(r, g, b, a)
    real(real32), intent(in) :: r, g, b, a
    pack_scaled = pack_rgba(int(r * scale), int(g * scale), int(b * scale), int(a * scale))
  end function pack_scaled

  subroutine box_filter_device(input, tmp, output, width, height, repeat)
    integer(int32), intent(in) :: input(:)
    integer(int32), intent(inout) :: tmp(:), output(:)
    integer, intent(in) :: width, height, repeat
    integer :: cycle, idx, x, y, offset, yy
    real(real32) :: rs, gs, bs, as, r, g, b, a
    real(real64) :: start_time, end_time, avg_us

    start_time = omp_get_wtime()
    do cycle = 1, repeat
      !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(idx, offset, rs, gs, bs, as, r, g, b, a)
      do y = 1, height
        do x = 1, width
          rs = 0.0_real32
          gs = 0.0_real32
          bs = 0.0_real32
          as = 0.0_real32
          do offset = -radius, radius
            if (x + offset >= 1 .and. x + offset <= width) then
              call unpack_rgba(input((y - 1) * width + x + offset), r, g, b, a)
              rs = rs + r
              gs = gs + g
              bs = bs + b
              as = as + a
            end if
          end do
          idx = (y - 1) * width + x
          tmp(idx) = pack_scaled(rs, gs, bs, as)
        end do
      end do
      !$omp end target teams distribute parallel do

      !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(idx, yy, rs, gs, bs, as, r, g, b, a)
      do y = 1, height
        do x = 1, width
          rs = 0.0_real32
          gs = 0.0_real32
          bs = 0.0_real32
          as = 0.0_real32
          do offset = -radius, radius
            yy = min(height, max(1, y + offset))
            call unpack_rgba(tmp((yy - 1) * width + x), r, g, b, a)
            rs = rs + r
            gs = gs + g
            bs = bs + b
            as = as + a
          end do
          idx = (y - 1) * width + x
          output(idx) = pack_scaled(rs, gs, bs, as)
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)
    write(*, '(A,F0.6,A)') 'Average kernel execution time ', avg_us, ' (us)'
  end subroutine box_filter_device

  subroutine box_filter_host(input, tmp, output, width, height)
    integer(int32), intent(in) :: input(:)
    integer(int32), intent(inout) :: tmp(:), output(:)
    integer, intent(in) :: width, height
    integer :: idx, x, y, offset, yy
    real(real32) :: rs, gs, bs, as, r, g, b, a
    do y = 1, height
      do x = 1, width
        rs = 0.0_real32
        gs = 0.0_real32
        bs = 0.0_real32
        as = 0.0_real32
        do offset = -radius, radius
          if (x + offset >= 1 .and. x + offset <= width) then
            call unpack_rgba(input((y - 1) * width + x + offset), r, g, b, a)
            rs = rs + r
            gs = gs + g
            bs = bs + b
            as = as + a
          end if
        end do
        idx = (y - 1) * width + x
        tmp(idx) = pack_scaled(rs, gs, bs, as)
      end do
    end do
    do y = 1, height
      do x = 1, width
        rs = 0.0_real32
        gs = 0.0_real32
        bs = 0.0_real32
        as = 0.0_real32
        do offset = -radius, radius
          yy = min(height, max(1, y + offset))
          call unpack_rgba(tmp((yy - 1) * width + x), r, g, b, a)
          rs = rs + r
          gs = gs + g
          bs = bs + b
          as = as + a
        end do
        idx = (y - 1) * width + x
        output(idx) = pack_scaled(rs, gs, bs, as)
      end do
    end do
  end subroutine box_filter_host

end program main
