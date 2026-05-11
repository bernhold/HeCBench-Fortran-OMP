program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: number_par_per_box = 100
  integer, parameter :: number_threads = 128
  real(real32), parameter :: alpha = 0.5_real32
  real(real32), parameter :: tolerance = 2.0e-3_real32

  integer :: boxes1d, number_boxes, space_elem, i
  integer, allocatable :: box_offset(:), box_nn(:), box_nei(:,:)
  real(real32), allocatable :: rv_v(:), rv_x(:), rv_y(:), rv_z(:), qv(:)
  real(real32), allocatable :: fv_v(:), fv_x(:), fv_y(:), fv_z(:)
  real(real32), allocatable :: ref_v(:), ref_x(:), ref_y(:), ref_z(:)
  integer(int64) :: seed
  real(real64) :: start_total, end_total, start_kernel, end_kernel
  logical :: ok
  character(len=256) :: flag

  write(*,'(A,I0,1X)') 'WG size of kernel = ', number_threads

  if (command_argument_count() /= 2) then
    write(*,'(A)', advance='no') 'Provide boxes1d argument, example: -boxes1d 16'
    stop 0
  end if
  call get_command_argument(1, flag)
  if (trim(flag) /= '-boxes1d') then
    write(*,'(A)') 'ERROR: Unknown argument'
    stop 0
  end if
  boxes1d = read_int_arg(2)
  if (boxes1d <= 0) then
    write(*,'(A)') 'ERROR: Wrong value to -boxes1d argument, cannot be <=0'
    stop 0
  end if

  write(*,'(A,I0,A,I0,A,I0)') 'Configuration used: arch = ', 0, ', cores = ', 1, ', boxes1d = ', boxes1d

  number_boxes = boxes1d * boxes1d * boxes1d
  space_elem = number_boxes * number_par_per_box
  allocate(box_offset(number_boxes), box_nn(number_boxes), box_nei(26, number_boxes))
  allocate(rv_v(space_elem), rv_x(space_elem), rv_y(space_elem), rv_z(space_elem), qv(space_elem))
  allocate(fv_v(space_elem), fv_x(space_elem), fv_y(space_elem), fv_z(space_elem))
  allocate(ref_v(space_elem), ref_x(space_elem), ref_y(space_elem), ref_z(space_elem))

  call initialize_boxes(boxes1d, number_boxes, box_offset, box_nn, box_nei)
  seed = 2_int64
  do i = 1, space_elem
    rv_v(i) = lava_random(seed)
    rv_x(i) = lava_random(seed)
    rv_y(i) = lava_random(seed)
    rv_z(i) = lava_random(seed)
  end do
  do i = 1, space_elem
    qv(i) = lava_random(seed)
  end do
  fv_v = 0.0_real32
  fv_x = 0.0_real32
  fv_y = 0.0_real32
  fv_z = 0.0_real32
  ref_v = 0.0_real32
  ref_x = 0.0_real32
  ref_y = 0.0_real32
  ref_z = 0.0_real32

  start_total = omp_get_wtime()
  !$omp target data map(to: box_offset(1:number_boxes), box_nn(1:number_boxes), box_nei(1:26,1:number_boxes), &
  !$omp& rv_v(1:space_elem), rv_x(1:space_elem), rv_y(1:space_elem), rv_z(1:space_elem), qv(1:space_elem)) &
  !$omp& map(tofrom: fv_v(1:space_elem), fv_x(1:space_elem), fv_y(1:space_elem), fv_z(1:space_elem))
  start_kernel = omp_get_wtime()
  call compute_forces_device(number_boxes, box_offset, box_nn, box_nei, rv_v, rv_x, rv_y, rv_z, qv, &
      fv_v, fv_x, fv_y, fv_z)
  end_kernel = omp_get_wtime()
  !$omp end target data
  end_total = omp_get_wtime()

  call compute_forces_host(number_boxes, box_offset, box_nn, box_nei, rv_v, rv_x, rv_y, rv_z, qv, &
      ref_v, ref_x, ref_y, ref_z)
  ok = compare_forces(space_elem, fv_v, fv_x, fv_y, fv_z, ref_v, ref_x, ref_y, ref_z)

  write(*,'(A)') 'Device offloading time:'
  write(*,'(F0.12,A)') real(end_total - start_total, real32), ' s'
  write(*,'(A)') 'Kernel execution time:'
  write(*,'(F0.12,A)') real(end_kernel - start_kernel, real32), ' s'
  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(box_offset, box_nn, box_nei, rv_v, rv_x, rv_y, rv_z, qv, fv_v, fv_x, fv_y, fv_z)
  deallocate(ref_v, ref_x, ref_y, ref_z)

contains

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  real(real32) function lava_random(seed) result(value)
    integer(int64), intent(inout) :: seed
    seed = mod(seed * 1103515245_int64 + 12345_int64, 2147483648_int64)
    value = real(mod(seed, 10_int64) + 1_int64, real32) / 10.0_real32
  end function lava_random

  subroutine initialize_boxes(boxes1d, number_boxes, box_offset, box_nn, box_nei)
    integer, intent(in) :: boxes1d, number_boxes
    integer, intent(out) :: box_offset(:), box_nn(:), box_nei(:,:)
    integer :: i, j, k, l, m, n, nh, neighbor
    box_nei = 0
    nh = 0
    do i = 0, boxes1d - 1
      do j = 0, boxes1d - 1
        do k = 0, boxes1d - 1
          nh = nh + 1
          box_offset(nh) = (nh - 1) * number_par_per_box + 1
          box_nn(nh) = 0
          do l = -1, 1
            do m = -1, 1
              do n = -1, 1
                if ((i + l) >= 0 .and. (j + m) >= 0 .and. (k + n) >= 0 .and. &
                    (i + l) < boxes1d .and. (j + m) < boxes1d .and. (k + n) < boxes1d .and. &
                    .not. (l == 0 .and. m == 0 .and. n == 0)) then
                  box_nn(nh) = box_nn(nh) + 1
                  neighbor = ((i + l) * boxes1d * boxes1d) + ((j + m) * boxes1d) + (k + n) + 1
                  box_nei(box_nn(nh), nh) = neighbor
                end if
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine initialize_boxes

  subroutine compute_forces_device(number_boxes, box_offset, box_nn, box_nei, rv_v, rv_x, rv_y, rv_z, qv, &
      fv_v, fv_x, fv_y, fv_z)
    integer, intent(in) :: number_boxes, box_offset(:), box_nn(:), box_nei(:,:)
    real(real32), intent(in) :: rv_v(:), rv_x(:), rv_y(:), rv_z(:), qv(:)
    real(real32), intent(inout) :: fv_v(:), fv_x(:), fv_y(:), fv_z(:)
    integer :: bx, pi, k, j, first_i, first_j, pointer, home_idx, nei_idx
    real(real32) :: a2, r2, u2, vij, fs, dx, dy, dz

    a2 = 2.0_real32 * alpha * alpha
    !$omp target teams distribute parallel do collapse(2) thread_limit(number_threads) &
    !$omp& private(first_i, first_j, pointer, home_idx, nei_idx, k, j, r2, u2, vij, fs, dx, dy, dz)
    do bx = 1, number_boxes
      do pi = 1, number_par_per_box
        first_i = box_offset(bx)
        home_idx = first_i + pi - 1
        do k = 0, box_nn(bx)
          if (k == 0) then
            pointer = bx
          else
            pointer = box_nei(k, bx)
          end if
          first_j = box_offset(pointer)
          do j = 1, number_par_per_box
            nei_idx = first_j + j - 1
            r2 = rv_v(home_idx) + rv_v(nei_idx) - &
                (rv_x(home_idx) * rv_x(nei_idx) + rv_y(home_idx) * rv_y(nei_idx) + rv_z(home_idx) * rv_z(nei_idx))
            u2 = a2 * r2
            vij = exp(-u2)
            fs = 2.0_real32 * vij
            dx = rv_x(home_idx) - rv_x(nei_idx)
            dy = rv_y(home_idx) - rv_y(nei_idx)
            dz = rv_z(home_idx) - rv_z(nei_idx)
            fv_v(home_idx) = fv_v(home_idx) + qv(nei_idx) * vij
            fv_x(home_idx) = fv_x(home_idx) + qv(nei_idx) * fs * dx
            fv_y(home_idx) = fv_y(home_idx) + qv(nei_idx) * fs * dy
            fv_z(home_idx) = fv_z(home_idx) + qv(nei_idx) * fs * dz
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine compute_forces_device

  subroutine compute_forces_host(number_boxes, box_offset, box_nn, box_nei, rv_v, rv_x, rv_y, rv_z, qv, &
      fv_v, fv_x, fv_y, fv_z)
    integer, intent(in) :: number_boxes, box_offset(:), box_nn(:), box_nei(:,:)
    real(real32), intent(in) :: rv_v(:), rv_x(:), rv_y(:), rv_z(:), qv(:)
    real(real32), intent(inout) :: fv_v(:), fv_x(:), fv_y(:), fv_z(:)
    integer :: bx, pi, k, j, first_i, first_j, pointer, home_idx, nei_idx
    real(real32) :: a2, r2, u2, vij, fs, dx, dy, dz

    a2 = 2.0_real32 * alpha * alpha
    do bx = 1, number_boxes
      first_i = box_offset(bx)
      do pi = 1, number_par_per_box
        home_idx = first_i + pi - 1
        do k = 0, box_nn(bx)
          if (k == 0) then
            pointer = bx
          else
            pointer = box_nei(k, bx)
          end if
          first_j = box_offset(pointer)
          do j = 1, number_par_per_box
            nei_idx = first_j + j - 1
            r2 = rv_v(home_idx) + rv_v(nei_idx) - &
                (rv_x(home_idx) * rv_x(nei_idx) + rv_y(home_idx) * rv_y(nei_idx) + rv_z(home_idx) * rv_z(nei_idx))
            u2 = a2 * r2
            vij = exp(-u2)
            fs = 2.0_real32 * vij
            dx = rv_x(home_idx) - rv_x(nei_idx)
            dy = rv_y(home_idx) - rv_y(nei_idx)
            dz = rv_z(home_idx) - rv_z(nei_idx)
            fv_v(home_idx) = fv_v(home_idx) + qv(nei_idx) * vij
            fv_x(home_idx) = fv_x(home_idx) + qv(nei_idx) * fs * dx
            fv_y(home_idx) = fv_y(home_idx) + qv(nei_idx) * fs * dy
            fv_z(home_idx) = fv_z(home_idx) + qv(nei_idx) * fs * dz
          end do
        end do
      end do
    end do
  end subroutine compute_forces_host

  logical function compare_forces(n, fv_v, fv_x, fv_y, fv_z, ref_v, ref_x, ref_y, ref_z) result(ok)
    integer, intent(in) :: n
    real(real32), intent(in) :: fv_v(:), fv_x(:), fv_y(:), fv_z(:), ref_v(:), ref_x(:), ref_y(:), ref_z(:)
    integer :: i
    ok = .true.
    do i = 1, n
      if (abs(fv_v(i) - ref_v(i)) > tolerance .or. abs(fv_x(i) - ref_x(i)) > tolerance .or. &
          abs(fv_y(i) - ref_y(i)) > tolerance .or. abs(fv_z(i) - ref_z(i)) > tolerance) then
        ok = .false.
        return
      end if
    end do
  end function compare_forces

end program main
