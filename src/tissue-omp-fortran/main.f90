program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg
  integer, parameter :: nnt_dev = 32 * 32 * 32
  integer, parameter :: nsp = 2
  integer, parameter :: step = 4
  integer, parameter :: isp = 2
  integer :: dim, repeat, nnt, i, iter
  integer(int32) :: seed
  integer, allocatable :: tisspoints(:)
  real(real32), allocatable :: gtt(:), gbartt(:), ct(:), ctprev(:), qt(:), ct_gold(:)
  real(real64) :: start_time, elapsed
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <dimension of a 3D grid> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) dim
  if (dim > 32) then
    write(*,'(A)') 'Maximum dimension is 32'
    stop 1
  end if
  call get_command_argument(2, arg); read(arg, *) repeat
  if (dim <= 0 .or. repeat <= 0) stop 1

  nnt = dim * dim * dim
  allocate(tisspoints(3 * nnt_dev))
  allocate(gtt(nsp * nnt_dev), gbartt(nsp * nnt_dev))
  allocate(ct(nnt_dev), ctprev(nnt_dev), qt(nnt_dev), ct_gold(nnt_dev))

  seed = 1_int32
  do i = 1, 3 * nnt_dev
    tisspoints(i) = modulo(c_rand(seed), nnt_dev / 3)
  end do
  do i = 1, nsp * nnt_dev
    gtt(i) = real(c_rand(seed), real32) / 32767.0_real32
    gbartt(i) = real(c_rand(seed), real32) / 32767.0_real32
  end do
  do i = 1, nnt_dev
    ct(i) = 0.0_real32
    ct_gold(i) = 0.0_real32
    ctprev(i) = real(c_rand(seed), real32) / 32767.0_real32
    qt(i) = real(c_rand(seed), real32) / 32767.0_real32
  end do

  !$omp target data map(to: tisspoints(1:3*nnt_dev), gtt(1:nsp*nnt_dev), gbartt(1:nsp*nnt_dev), &
  !$omp& ctprev(1:nnt_dev), qt(1:nnt_dev)) map(tofrom: ct(1:nnt_dev))
  do i = 1, 2
    call tissue(tisspoints, gtt, gbartt, ct, ctprev, qt, nnt, nnt_dev, isp)
  end do

  do i = 1, 2
    call reference(tisspoints, gtt, gbartt, ct_gold, ctprev, qt, nnt, nnt_dev, isp)
  end do

  !$omp target update from(ct(1:nnt_dev))
  ok = .true.
  do i = 1, nnt_dev
    if (abs(ct(i) - ct_gold(i)) > 1.0e-1_real32) then
      write(*,'(A,I0,A,F0.6,1X,F0.6)') '@', i - 1, ': ', ct(i), ct_gold(i)
      ok = .false.
      exit
    end if
  end do
  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  start_time = omp_get_wtime()
  do iter = 1, repeat
    call tissue(tisspoints, gtt, gbartt, ct, ctprev, qt, nnt, nnt_dev, isp)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', elapsed / real(repeat, real64), ' (s)'

  deallocate(tisspoints, gtt, gbartt, ct, ctprev, qt, ct_gold)

contains

  subroutine tissue(tisspoints, gtt, gbartt, ct, ctprev, qt, nnt, nnt_dev, isp)
    integer, intent(in) :: nnt, nnt_dev, isp
    integer, intent(in) :: tisspoints(3 * nnt_dev)
    real(real32), intent(in) :: gtt(nsp * nnt_dev), gbartt(nsp * nnt_dev), ctprev(nnt_dev), qt(nnt_dev)
    real(real32), intent(out) :: ct(nnt_dev)
    integer :: itp, jtp, ix, iy, iz, jx, jy, jz, ixyz
    real(real32) :: p

    !$omp target teams distribute parallel do private(jtp, ix, iy, iz, jx, jy, jz, ixyz, p) thread_limit(256)
    do itp = 0, nnt - 1
      ix = tisspoints(itp + 1)
      iy = tisspoints(itp + nnt + 1)
      iz = tisspoints(itp + 2 * nnt + 1)
      p = 0.0_real32
      do jtp = 0, nnt - 1
        jx = tisspoints(jtp + 1)
        jy = tisspoints(jtp + nnt + 1)
        jz = tisspoints(jtp + 2 * nnt + 1)
        ixyz = abs(jx - ix) + abs(jy - iy) + abs(jz - iz) + (isp - 1) * nnt_dev
        p = p + gtt(ixyz + 1) * ctprev(jtp + 1) + gbartt(ixyz + 1) * qt(jtp + 1)
      end do
      ct(itp + 1) = p
    end do
    !$omp end target teams distribute parallel do
  end subroutine tissue

  subroutine reference(tisspoints, gtt, gbartt, ct, ctprev, qt, nnt, nnt_dev, isp)
    integer, intent(in) :: nnt, nnt_dev, isp
    integer, intent(in) :: tisspoints(3 * nnt_dev)
    real(real32), intent(in) :: gtt(nsp * nnt_dev), gbartt(nsp * nnt_dev), ctprev(nnt_dev), qt(nnt_dev)
    real(real32), intent(out) :: ct(nnt_dev)
    integer :: itp, jtp, ix, iy, iz, jx, jy, jz, ixyz
    real(real32) :: p

    do itp = 0, nnt - 1
      ix = tisspoints(itp + 1)
      iy = tisspoints(itp + nnt + 1)
      iz = tisspoints(itp + 2 * nnt + 1)
      p = 0.0_real32
      do jtp = 0, nnt - 1
        jx = tisspoints(jtp + 1)
        jy = tisspoints(jtp + nnt + 1)
        jz = tisspoints(jtp + 2 * nnt + 1)
        ixyz = abs(jx - ix) + abs(jy - iy) + abs(jz - iz) + (isp - 1) * nnt_dev
        p = p + gtt(ixyz + 1) * ctprev(jtp + 1) + gbartt(ixyz + 1) * qt(jtp + 1)
      end do
      ct(itp + 1) = p
    end do
  end subroutine reference

  integer(int32) function c_rand(seed)
    integer(int32), intent(inout) :: seed
    integer(int64) :: next_value

    next_value = mod(1103515245_int64 * int(seed, int64) + 12345_int64, 2147483648_int64)
    seed = int(next_value, int32)
    c_rand = iand(seed / 65536_int32, 32767_int32)
  end function c_rand

end program main
