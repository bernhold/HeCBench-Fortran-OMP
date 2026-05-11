program dslash_main
  use iso_fortran_env, only: int64, real64, output_unit
  use omp_lib, only: omp_get_wtime
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: ldim = 32
  integer, parameter :: iterations = 100
  integer, parameter :: warmups = 1
  integer, parameter :: total_sites = ldim * ldim * ldim * ldim
  integer, parameter :: even_sites = total_sites / 2
  real(dp), parameter :: eps = 2.0e-6_dp

  integer :: workgroup_size, arg_status
  character(len=64) :: arg
  real(dp), allocatable :: src_r(:,:), src_i(:,:), dst_r(:,:), dst_i(:,:)
  real(dp), allocatable :: chk_r(:,:), chk_i(:,:)
  real(dp), allocatable :: fat_r(:,:,:,:), fat_i(:,:,:,:)
  real(dp), allocatable :: lng_r(:,:,:,:), lng_i(:,:,:,:)
  real(dp), allocatable :: fatbck_r(:,:,:,:), fatbck_i(:,:,:,:)
  real(dp), allocatable :: lngbck_r(:,:,:,:), lngbck_i(:,:,:,:)
  integer(int64), allocatable :: fwd(:,:), bck(:,:), fwd3(:,:), bck3(:,:)
  real(dp) :: total_time, tflop, memory_usage, memory_allocated
  real(dp) :: max_error

  if (command_argument_count() < 1) then
    write(*,'(A)') "Usage <workgroup size>"
    stop 1
  end if

  call get_command_argument(1, arg, status=arg_status)
  if (arg_status /= 0) stop 1
  read(arg, *) workgroup_size

  allocate(src_r(3,total_sites), src_i(3,total_sites))
  allocate(dst_r(3,total_sites), dst_i(3,total_sites))
  allocate(chk_r(3,total_sites), chk_i(3,total_sites))
  allocate(fat_r(3,3,4,total_sites), fat_i(3,3,4,total_sites))
  allocate(lng_r(3,3,4,total_sites), lng_i(3,3,4,total_sites))
  allocate(fatbck_r(3,3,4,total_sites), fatbck_i(3,3,4,total_sites))
  allocate(lngbck_r(3,3,4,total_sites), lngbck_i(3,3,4,total_sites))
  allocate(fwd(4,total_sites), bck(4,total_sites), fwd3(4,total_sites), bck3(4,total_sites))

  dst_r = 0.0_dp
  dst_i = 0.0_dp
  chk_r = 0.0_dp
  chk_i = 0.0_dp
  fatbck_r = 0.0_dp
  fatbck_i = 0.0_dp
  lngbck_r = 0.0_dp
  lngbck_i = 0.0_dp

  call set_neighbors(fwd, bck, fwd3, bck3)
  call make_data(src_r, src_i, fat_r, fat_i, lng_r, lng_i)

  write(*,'(A,I0,A)') "Number of sites = ", ldim, "^4"
  write(*,'(A,I0,A,I0,A)') "Executing ", iterations, " iterations with ", warmups, " warmups"
  if (workgroup_size /= 0) write(*,'(A,I0)') "Threads per group = ", workgroup_size
  flush(output_unit)

  total_time = dslash_fn(src_r, src_i, dst_r, dst_i, fat_r, fat_i, lng_r, lng_i, &
                         fatbck_r, fatbck_i, lngbck_r, lngbck_i, &
                         fwd, bck, fwd3, bck3, workgroup_size)
  write(*,'(A,ES15.8,A)') "Total execution time = ", total_time, " secs"

  write(*,'(A)') "Validating the result"
  flush(output_unit)
  call dslash_host(src_r, src_i, chk_r, chk_i, fat_r, fat_i, lng_r, lng_i, &
                   fatbck_r, fatbck_i, lngbck_r, lngbck_i, fwd, bck, fwd3, bck3)
  max_error = max(maxval(abs(dst_r(:,1:even_sites) - chk_r(:,1:even_sites))), &
                  maxval(abs(dst_i(:,1:even_sites) - chk_i(:,1:even_sites))))
  if (max_error >= eps) then
    write(*,'(A,ES15.8)') "Validation failed, max error = ", max_error
    stop 2
  end if

  tflop = real(iterations, dp) * real(even_sites, dp) * 1182.0_dp
  write(*,'(A,ES15.8)') "Total GFLOP/s = ", tflop / total_time / 1.0e9_dp

  memory_usage = real(even_sites, dp) * &
      (real(144 * 4 * 4, dp) + real(48 * 16, dp) + real(8 * 16, dp) + real(48, dp))
  write(*,'(A,ES15.8)') "Total GByte/s (GPU memory) = ", &
      real(iterations, dp) * memory_usage / total_time / 1.0e9_dp

  memory_allocated = real(total_sites, dp) * &
      (real(144 * 4 * 4, dp) + real(48 * 2, dp) + real(8 * 4 * 4, dp))
  write(*,'(A,ES15.8)') "Total allocation for matrices = ", memory_allocated / 1048576.0_dp
  write(*,'(A,ES15.8)') "Approximate memory usage = ", 0.0_dp

contains

  integer function node_index(x, y, z, t) result(idx)
    integer, intent(in) :: x, y, z, t
    integer :: xr, yr, zr, tr, linear

    xr = modulo(x + ldim, ldim)
    yr = modulo(y + ldim, ldim)
    zr = modulo(z + ldim, ldim)
    tr = modulo(t + ldim, ldim)
    linear = xr + ldim * (yr + ldim * (zr + ldim * tr))
    if (modulo(x + y + z + t, 2) == 0) then
      idx = linear / 2 + 1
    else
      idx = (linear + total_sites) / 2 + 1
    end if
  end function node_index

  subroutine set_neighbors(fwd, bck, fwd3, bck3)
    integer(int64), intent(out) :: fwd(4,total_sites), bck(4,total_sites)
    integer(int64), intent(out) :: fwd3(4,total_sites), bck3(4,total_sites)
    integer :: x, y, z, t, idx

    do t = 0, ldim - 1
      do z = 0, ldim - 1
        do y = 0, ldim - 1
          do x = 0, ldim - 1
            idx = node_index(x, y, z, t)
            fwd(1,idx) = node_index(x + 1, y, z, t)
            bck(1,idx) = node_index(x - 1, y, z, t)
            fwd(2,idx) = node_index(x, y + 1, z, t)
            bck(2,idx) = node_index(x, y - 1, z, t)
            fwd(3,idx) = node_index(x, y, z + 1, t)
            bck(3,idx) = node_index(x, y, z - 1, t)
            fwd(4,idx) = node_index(x, y, z, t + 1)
            bck(4,idx) = node_index(x, y, z, t - 1)
            fwd3(1,idx) = node_index(x + 3, y, z, t)
            bck3(1,idx) = node_index(x - 3, y, z, t)
            fwd3(2,idx) = node_index(x, y + 3, z, t)
            bck3(2,idx) = node_index(x, y - 3, z, t)
            fwd3(3,idx) = node_index(x, y, z + 3, t)
            bck3(3,idx) = node_index(x, y, z - 3, t)
            fwd3(4,idx) = node_index(x, y, z, t + 3)
            bck3(4,idx) = node_index(x, y, z, t - 3)
          end do
        end do
      end do
    end do
  end subroutine set_neighbors

  subroutine make_data(src_r, src_i, fat_r, fat_i, lng_r, lng_i)
    real(dp), intent(out) :: src_r(3,total_sites), src_i(3,total_sites)
    real(dp), intent(out) :: fat_r(3,3,4,total_sites), fat_i(3,3,4,total_sites)
    real(dp), intent(out) :: lng_r(3,3,4,total_sites), lng_i(3,3,4,total_sites)
    real(dp) :: values(2)
    integer :: site, dir

    call random_seed()
    do site = 1, total_sites
      call random_number(values)
      src_r(:,site) = 2.0_dp * values(1) - 1.0_dp
      src_i(:,site) = 2.0_dp * values(2) - 1.0_dp
      do dir = 1, 4
        call random_number(values)
        fat_r(:,:,dir,site) = 2.0_dp * values(1) - 1.0_dp
        fat_i(:,:,dir,site) = 2.0_dp * values(2) - 1.0_dp
        call random_number(values)
        lng_r(:,:,dir,site) = 2.0_dp * values(1) - 1.0_dp
        lng_i(:,:,dir,site) = 2.0_dp * values(2) - 1.0_dp
      end do
    end do
  end subroutine make_data

  real(dp) function dslash_fn(src_r, src_i, dst_r, dst_i, fat_r, fat_i, lng_r, lng_i, &
                              fatbck_r, fatbck_i, lngbck_r, lngbck_i, &
                              fwd, bck, fwd3, bck3, wgsize) result(ttotal)
    real(dp), intent(in) :: src_r(3,total_sites), src_i(3,total_sites)
    real(dp), intent(inout) :: dst_r(3,total_sites), dst_i(3,total_sites)
    real(dp), intent(in) :: fat_r(3,3,4,total_sites), fat_i(3,3,4,total_sites)
    real(dp), intent(in) :: lng_r(3,3,4,total_sites), lng_i(3,3,4,total_sites)
    real(dp), intent(inout) :: fatbck_r(3,3,4,total_sites), fatbck_i(3,3,4,total_sites)
    real(dp), intent(inout) :: lngbck_r(3,3,4,total_sites), lngbck_i(3,3,4,total_sites)
    integer(int64), intent(in) :: fwd(4,total_sites), bck(4,total_sites)
    integer(int64), intent(in) :: fwd3(4,total_sites), bck3(4,total_sites)
    integer, intent(in) :: wgsize

    integer :: my_site, dir, row, col, k, iter
    integer(int64) :: src_site
    real(dp) :: start_time, stop_time
    real(dp) :: acc_r, acc_i, v_r(3), v_i(3)

    !$omp target data map(to: src_r, src_i, fat_r, fat_i, lng_r, lng_i, fwd, bck, fwd3, bck3) &
    !$omp& map(tofrom: dst_r, dst_i, fatbck_r, fatbck_i, lngbck_r, lngbck_i)

    !$omp target teams distribute parallel do thread_limit(1) &
    !$omp& private(my_site, dir, row, col, src_site)
    do my_site = 1, even_sites
      do dir = 1, 4
        src_site = bck(dir,my_site)
        do row = 1, 3
          do col = 1, 3
            fatbck_r(row,col,dir,my_site) = fat_r(col,row,dir,src_site)
            fatbck_i(row,col,dir,my_site) = -fat_i(col,row,dir,src_site)
          end do
        end do
        src_site = bck3(dir,my_site)
        do row = 1, 3
          do col = 1, 3
            lngbck_r(row,col,dir,my_site) = lng_r(col,row,dir,src_site)
            lngbck_i(row,col,dir,my_site) = -lng_i(col,row,dir,src_site)
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

    write(*,'(A)') "Running dslash loop"
    write(*,'(A,I0)') "Setting number of work items to ", even_sites
    write(*,'(A,I0)') "Setting workgroup size to ", wgsize
    flush(output_unit)

    start_time = omp_get_wtime()
    do iter = 0, iterations + warmups - 1
      if (iter == warmups) start_time = omp_get_wtime()
      !$omp target teams distribute parallel do thread_limit(wgsize) &
      !$omp& private(my_site, dir, row, col, k, src_site, acc_r, acc_i, v_r, v_i)
      do my_site = 1, even_sites
        do row = 1, 3
          acc_r = 0.0_dp
          acc_i = 0.0_dp
          src_site = fwd(1,my_site)
          do col = 1, 3
            acc_r = acc_r + fat_r(row,col,1,my_site) * src_r(col,src_site) - &
                            fat_i(row,col,1,my_site) * src_i(col,src_site)
            acc_i = acc_i + fat_r(row,col,1,my_site) * src_i(col,src_site) + &
                            fat_i(row,col,1,my_site) * src_r(col,src_site)
          end do
          dst_r(row,my_site) = acc_r
          dst_i(row,my_site) = acc_i
        end do
        do dir = 2, 4
          src_site = fwd(dir,my_site)
          do row = 1, 3
            acc_r = 0.0_dp
            acc_i = 0.0_dp
            do col = 1, 3
              acc_r = acc_r + fat_r(row,col,dir,my_site) * src_r(col,src_site) - &
                              fat_i(row,col,dir,my_site) * src_i(col,src_site)
              acc_i = acc_i + fat_r(row,col,dir,my_site) * src_i(col,src_site) + &
                              fat_i(row,col,dir,my_site) * src_r(col,src_site)
            end do
            dst_r(row,my_site) = dst_r(row,my_site) + acc_r
            dst_i(row,my_site) = dst_i(row,my_site) + acc_i
          end do
        end do

        do row = 1, 3
          v_r(row) = 0.0_dp
          v_i(row) = 0.0_dp
        end do
        do dir = 1, 4
          src_site = fwd3(dir,my_site)
          do row = 1, 3
            acc_r = 0.0_dp
            acc_i = 0.0_dp
            do col = 1, 3
              acc_r = acc_r + lng_r(row,col,dir,my_site) * src_r(col,src_site) - &
                              lng_i(row,col,dir,my_site) * src_i(col,src_site)
              acc_i = acc_i + lng_r(row,col,dir,my_site) * src_i(col,src_site) + &
                              lng_i(row,col,dir,my_site) * src_r(col,src_site)
            end do
            v_r(row) = v_r(row) + acc_r
            v_i(row) = v_i(row) + acc_i
          end do
        end do
        do row = 1, 3
          dst_r(row,my_site) = dst_r(row,my_site) + v_r(row)
          dst_i(row,my_site) = dst_i(row,my_site) + v_i(row)
        end do

        do row = 1, 3
          v_r(row) = 0.0_dp
          v_i(row) = 0.0_dp
        end do
        do dir = 1, 4
          src_site = bck(dir,my_site)
          do row = 1, 3
            acc_r = 0.0_dp
            acc_i = 0.0_dp
            do col = 1, 3
              acc_r = acc_r + fatbck_r(row,col,dir,my_site) * src_r(col,src_site) - &
                              fatbck_i(row,col,dir,my_site) * src_i(col,src_site)
              acc_i = acc_i + fatbck_r(row,col,dir,my_site) * src_i(col,src_site) + &
                              fatbck_i(row,col,dir,my_site) * src_r(col,src_site)
            end do
            v_r(row) = v_r(row) + acc_r
            v_i(row) = v_i(row) + acc_i
          end do
        end do
        do row = 1, 3
          dst_r(row,my_site) = dst_r(row,my_site) - v_r(row)
          dst_i(row,my_site) = dst_i(row,my_site) - v_i(row)
        end do

        do row = 1, 3
          v_r(row) = 0.0_dp
          v_i(row) = 0.0_dp
        end do
        do dir = 1, 4
          src_site = bck3(dir,my_site)
          do row = 1, 3
            acc_r = 0.0_dp
            acc_i = 0.0_dp
            do col = 1, 3
              acc_r = acc_r + lngbck_r(row,col,dir,my_site) * src_r(col,src_site) - &
                              lngbck_i(row,col,dir,my_site) * src_i(col,src_site)
              acc_i = acc_i + lngbck_r(row,col,dir,my_site) * src_i(col,src_site) + &
                              lngbck_i(row,col,dir,my_site) * src_r(col,src_site)
            end do
            v_r(row) = v_r(row) + acc_r
            v_i(row) = v_i(row) + acc_i
          end do
        end do
        do row = 1, 3
          dst_r(row,my_site) = dst_r(row,my_site) - v_r(row)
          dst_i(row,my_site) = dst_i(row,my_site) - v_i(row)
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    stop_time = omp_get_wtime()

    !$omp end target data
    ttotal = stop_time - start_time
  end function dslash_fn

  subroutine dslash_host(src_r, src_i, dst_r, dst_i, fat_r, fat_i, lng_r, lng_i, &
                         fatbck_r, fatbck_i, lngbck_r, lngbck_i, fwd, bck, fwd3, bck3)
    real(dp), intent(in) :: src_r(3,total_sites), src_i(3,total_sites)
    real(dp), intent(out) :: dst_r(3,total_sites), dst_i(3,total_sites)
    real(dp), intent(in) :: fat_r(3,3,4,total_sites), fat_i(3,3,4,total_sites)
    real(dp), intent(in) :: lng_r(3,3,4,total_sites), lng_i(3,3,4,total_sites)
    real(dp), intent(in) :: fatbck_r(3,3,4,total_sites), fatbck_i(3,3,4,total_sites)
    real(dp), intent(in) :: lngbck_r(3,3,4,total_sites), lngbck_i(3,3,4,total_sites)
    integer(int64), intent(in) :: fwd(4,total_sites), bck(4,total_sites)
    integer(int64), intent(in) :: fwd3(4,total_sites), bck3(4,total_sites)
    integer :: my_site, dir, row, col
    integer(int64) :: src_site
    real(dp) :: acc_r, acc_i, v_r(3), v_i(3)

    dst_r = 0.0_dp
    dst_i = 0.0_dp
    !$omp parallel do private(my_site, dir, row, col, src_site, acc_r, acc_i, v_r, v_i)
    do my_site = 1, even_sites
      do row = 1, 3
        acc_r = 0.0_dp
        acc_i = 0.0_dp
        src_site = fwd(1,my_site)
        do col = 1, 3
          acc_r = acc_r + fat_r(row,col,1,my_site) * src_r(col,src_site) - &
                          fat_i(row,col,1,my_site) * src_i(col,src_site)
          acc_i = acc_i + fat_r(row,col,1,my_site) * src_i(col,src_site) + &
                          fat_i(row,col,1,my_site) * src_r(col,src_site)
        end do
        dst_r(row,my_site) = acc_r
        dst_i(row,my_site) = acc_i
      end do
      do dir = 2, 4
        src_site = fwd(dir,my_site)
        do row = 1, 3
          acc_r = 0.0_dp
          acc_i = 0.0_dp
          do col = 1, 3
            acc_r = acc_r + fat_r(row,col,dir,my_site) * src_r(col,src_site) - &
                            fat_i(row,col,dir,my_site) * src_i(col,src_site)
            acc_i = acc_i + fat_r(row,col,dir,my_site) * src_i(col,src_site) + &
                            fat_i(row,col,dir,my_site) * src_r(col,src_site)
          end do
          dst_r(row,my_site) = dst_r(row,my_site) + acc_r
          dst_i(row,my_site) = dst_i(row,my_site) + acc_i
        end do
      end do

      v_r = 0.0_dp
      v_i = 0.0_dp
      do dir = 1, 4
        src_site = fwd3(dir,my_site)
        do row = 1, 3
          acc_r = 0.0_dp
          acc_i = 0.0_dp
          do col = 1, 3
            acc_r = acc_r + lng_r(row,col,dir,my_site) * src_r(col,src_site) - &
                            lng_i(row,col,dir,my_site) * src_i(col,src_site)
            acc_i = acc_i + lng_r(row,col,dir,my_site) * src_i(col,src_site) + &
                            lng_i(row,col,dir,my_site) * src_r(col,src_site)
          end do
          v_r(row) = v_r(row) + acc_r
          v_i(row) = v_i(row) + acc_i
        end do
      end do
      dst_r(:,my_site) = dst_r(:,my_site) + v_r
      dst_i(:,my_site) = dst_i(:,my_site) + v_i

      v_r = 0.0_dp
      v_i = 0.0_dp
      do dir = 1, 4
        src_site = bck(dir,my_site)
        do row = 1, 3
          acc_r = 0.0_dp
          acc_i = 0.0_dp
          do col = 1, 3
            acc_r = acc_r + fatbck_r(row,col,dir,my_site) * src_r(col,src_site) - &
                            fatbck_i(row,col,dir,my_site) * src_i(col,src_site)
            acc_i = acc_i + fatbck_r(row,col,dir,my_site) * src_i(col,src_site) + &
                            fatbck_i(row,col,dir,my_site) * src_r(col,src_site)
          end do
          v_r(row) = v_r(row) + acc_r
          v_i(row) = v_i(row) + acc_i
        end do
      end do
      dst_r(:,my_site) = dst_r(:,my_site) - v_r
      dst_i(:,my_site) = dst_i(:,my_site) - v_i

      v_r = 0.0_dp
      v_i = 0.0_dp
      do dir = 1, 4
        src_site = bck3(dir,my_site)
        do row = 1, 3
          acc_r = 0.0_dp
          acc_i = 0.0_dp
          do col = 1, 3
            acc_r = acc_r + lngbck_r(row,col,dir,my_site) * src_r(col,src_site) - &
                            lngbck_i(row,col,dir,my_site) * src_i(col,src_site)
            acc_i = acc_i + lngbck_r(row,col,dir,my_site) * src_i(col,src_site) + &
                            lngbck_i(row,col,dir,my_site) * src_r(col,src_site)
          end do
          v_r(row) = v_r(row) + acc_r
          v_i(row) = v_i(row) + acc_i
        end do
      end do
      dst_r(:,my_site) = dst_r(:,my_site) - v_r
      dst_i(:,my_site) = dst_i(:,my_site) - v_i
    end do
    !$omp end parallel do
  end subroutine dslash_host

end program dslash_main
