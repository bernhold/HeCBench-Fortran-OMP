program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
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

  character(len=256) :: arg0, arg
  integer :: width, height, repeat, num_pts, i
  real(real32), allocatable :: px(:), py(:), pz(:)
  real(real32), allocatable :: nx(:), ny(:), nz(:), nw(:)
  real(real32), allocatable :: rx(:), ry(:), rz(:), rw(:)
  real(real32) :: sx, sy, sz, sw
  real(real64) :: start_time, elapsed

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <width> <height> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) width
  call get_command_argument(2, arg); read(arg, *) height
  call get_command_argument(3, arg); read(arg, *) repeat
  if (width <= 0 .or. height <= 0 .or. repeat <= 0) stop 1

  num_pts = width * height
  allocate(px(num_pts), py(num_pts), pz(num_pts))
  allocate(nx(num_pts), ny(num_pts), nz(num_pts), nw(num_pts))
  allocate(rx(num_pts), ry(num_pts), rz(num_pts), rw(num_pts))

  call c_srand(123_c_int)
  do i = 1, num_pts
    px(i) = real(modulo(c_rand(), width), real32)
    py(i) = real(modulo(c_rand(), height), real32)
    pz(i) = real(modulo(c_rand(), 256_c_int), real32)
  end do
  nx = 0.0_real32; ny = 0.0_real32; nz = 0.0_real32; nw = 0.0_real32
  rx = 0.0_real32; ry = 0.0_real32; rz = 0.0_real32; rw = 0.0_real32

  !$omp target data map(to: px(1:num_pts), py(1:num_pts), pz(1:num_pts)) &
  !$omp& map(from: nx(1:num_pts), ny(1:num_pts), nz(1:num_pts), nw(1:num_pts))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call estimate_normals(px, py, pz, nx, ny, nz, nw, width, height, num_pts)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', &
    elapsed / real(repeat, real64), ' (s)'

  call estimate_normals_cpu(px, py, pz, rx, ry, rz, rw, width, height, num_pts)
  if (maxval(abs(nx - rx)) > 1.0e-5_real32 .or. maxval(abs(ny - ry)) > 1.0e-5_real32 .or. &
      maxval(abs(nz - rz)) > 1.0e-5_real32 .or. maxval(abs(nw - rw)) > 1.0e-5_real32) then
    error stop 'Fortran reference validation failed'
  end if

  sx = sum(nx)
  sy = sum(ny)
  sz = sum(nz)
  sw = sum(nw)
  write(*,'(A,F0.6,A,F0.6,A,F0.6,A,F0.6)') 'Checksum: x=', sx, ' y=', sy, &
    ' z=', sz, ' w=', sw

  deallocate(px, py, pz, nx, ny, nz, nw, rx, ry, rz, rw)

contains

  subroutine estimate_normals(px, py, pz, nx, ny, nz, nw, width, height, num_pts)
    integer, intent(in) :: width, height, num_pts
    real(real32), intent(in) :: px(:), py(:), pz(:)
    real(real32), intent(out) :: nx(:), ny(:), nz(:), nw(:)
    integer :: idx

    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& map(to: px(1:num_pts), py(1:num_pts), pz(1:num_pts)) &
    !$omp& map(from: nx(1:num_pts), ny(1:num_pts), nz(1:num_pts), nw(1:num_pts))
    do idx = 0, num_pts - 1
      call normal_estimate(px, py, pz, idx, width, height, nx(idx + 1), ny(idx + 1), &
        nz(idx + 1), nw(idx + 1))
    end do
    !$omp end target teams distribute parallel do
  end subroutine estimate_normals

  subroutine estimate_normals_cpu(px, py, pz, nx, ny, nz, nw, width, height, num_pts)
    integer, intent(in) :: width, height, num_pts
    real(real32), intent(in) :: px(:), py(:), pz(:)
    real(real32), intent(out) :: nx(:), ny(:), nz(:), nw(:)
    integer :: idx

    do idx = 0, num_pts - 1
      call normal_estimate(px, py, pz, idx, width, height, nx(idx + 1), ny(idx + 1), &
        nz(idx + 1), nw(idx + 1))
    end do
  end subroutine estimate_normals_cpu

  subroutine normal_estimate(px, py, pz, idx, width, height, ox, oy, oz, ow)
    real(real32), intent(in) :: px(:), py(:), pz(:)
    integer, intent(in) :: idx, width, height
    real(real32), intent(out) :: ox, oy, oz, ow
    integer :: p, x_idx, y_idx
    logical :: west_valid, east_valid, north_valid, south_valid
    real(real32) :: qx, qy, qz, hx, hy, hz, vx, vy, vz
    real(real32) :: nx0, ny0, nz0, curvature, len, inv_len, dotv

    p = idx + 1
    qx = px(p); qy = py(p); qz = pz(p)
    if (ieee_is_nan(qz)) then
      ox = 0.0_real32; oy = 0.0_real32; oz = 0.0_real32; ow = 0.0_real32
      return
    end if

    x_idx = modulo(idx, width)
    y_idx = idx / width

    west_valid = .false.
    east_valid = .false.
    north_valid = .false.
    south_valid = .false.
    if (x_idx > 1) west_valid = (.not. ieee_is_nan(pz(p - 1))) .and. abs(pz(p - 1) - qz) < 200.0_real32
    if (x_idx < width - 1) east_valid = (.not. ieee_is_nan(pz(p + 1))) .and. abs(pz(p + 1) - qz) < 200.0_real32
    if (y_idx > 1) north_valid = (.not. ieee_is_nan(pz(p - width))) .and. &
      abs(pz(p - width) - qz) < 200.0_real32
    if (y_idx < height - 1) south_valid = (.not. ieee_is_nan(pz(p + width))) .and. &
      abs(pz(p + width) - qz) < 200.0_real32

    if (west_valid .and. east_valid) then
      hx = px(p + 1) - px(p - 1)
      hy = py(p + 1) - py(p - 1)
      hz = pz(p + 1) - pz(p - 1)
    else if (west_valid .and. .not. east_valid) then
      hx = qx - px(p - 1)
      hy = qy - py(p - 1)
      hz = qz - pz(p - 1)
    else if (.not. west_valid .and. east_valid) then
      hx = px(p + 1) - qx
      hy = py(p + 1) - qy
      hz = pz(p + 1) - qz
    else
      ox = 0.0_real32; oy = 0.0_real32; oz = 0.0_real32; ow = 1.0_real32
      return
    end if

    if (south_valid .and. north_valid) then
      vx = px(p - width) - px(p + width)
      vy = py(p - width) - py(p + width)
      vz = pz(p - width) - pz(p + width)
    else if (south_valid .and. .not. north_valid) then
      vx = qx - px(p + width)
      vy = qy - py(p + width)
      vz = qz - pz(p + width)
    else if (.not. south_valid .and. north_valid) then
      vx = px(p - width) - qx
      vy = py(p - width) - qy
      vz = pz(p - width) - qz
    else
      ox = 0.0_real32; oy = 0.0_real32; oz = 0.0_real32; ow = 1.0_real32
      return
    end if

    nx0 = hy * vz - hz * vy
    ny0 = hz * vx - hx * vz
    nz0 = hx * vy - hy * vx
    len = sqrt(nx0 * nx0 + ny0 * ny0 + nz0 * nz0)
    curvature = merge(1.0_real32, 0.0_real32, abs(hz) > 0.04_real32 .or. &
      abs(vz) > 0.04_real32 .or. .not. west_valid .or. .not. east_valid .or. &
      .not. north_valid .or. .not. south_valid)

    inv_len = 1.0_real32 / len
    ox = nx0 * inv_len
    oy = ny0 * inv_len
    oz = nz0 * inv_len
    dotv = qx * ox + qy * oy + qz * oz
    if (dotv > 0.0_real32) then
      ox = -ox
      oy = -oy
      oz = -oz
    end if
    ow = curvature
  end subroutine normal_estimate

end program main
