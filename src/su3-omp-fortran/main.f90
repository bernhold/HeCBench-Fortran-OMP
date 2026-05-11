program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  integer(int64) :: iterations, ldim, threads_per_group, warmups
  integer(int32) :: verbose, device
  integer :: argc

  iterations = 100_int64
  ldim = 32_int64
  threads_per_group = 128_int64
  warmups = 1_int64
  verbose = 1_int32
  device = -1_int32

  argc = command_argument_count()
  call parse_args(argc, iterations, ldim, threads_per_group, verbose, device, warmups)
  call run_benchmark(iterations, ldim, threads_per_group, verbose, warmups)

contains

  subroutine parse_args(argc, iterations, ldim, threads_per_group, verbose, device, warmups)
    integer, intent(in) :: argc
    integer(int64), intent(inout) :: iterations, ldim, threads_per_group, warmups
    integer(int32), intent(inout) :: verbose, device
    character(len=256) :: arg, value
    integer :: i

    i = 1
    do while (i <= argc)
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('-h')
        write(*,'(A)') 'Usage: ./main [-i iterations] [-l lattice dimension] [-t threads per workgroup] [-d device] [-v verbosity level [0,1,2,3]] [-w warmups]'
        stop 1
      case ('-i')
        call get_command_argument(i + 1, value)
        read(value, *) iterations
        i = i + 1
      case ('-l')
        call get_command_argument(i + 1, value)
        read(value, *) ldim
        i = i + 1
      case ('-t')
        call get_command_argument(i + 1, value)
        read(value, *) threads_per_group
        i = i + 1
      case ('-v')
        call get_command_argument(i + 1, value)
        read(value, *) verbose
        i = i + 1
      case ('-d')
        call get_command_argument(i + 1, value)
        read(value, *) device
        i = i + 1
      case ('-w')
        call get_command_argument(i + 1, value)
        read(value, *) warmups
        i = i + 1
      case ('-n')
        i = i + 1
      end select
      i = i + 1
    end do
  end subroutine parse_args

  subroutine run_benchmark(iterations, ldim, threads_per_group, verbose, warmups)
    integer(int64), intent(in) :: iterations, ldim, warmups
    integer(int64), intent(inout) :: threads_per_group
    integer(int32), intent(in) :: verbose
    integer(int64) :: total_sites
    real(real64), allocatable :: a_real(:,:,:,:), a_imag(:,:,:,:)
    real(real64), allocatable :: b_real(:,:,:), b_imag(:,:,:)
    real(real64), allocatable :: c_real(:,:,:,:), c_imag(:,:,:,:)
    real(real64) :: ttotal, tflop, memory_usage, checksum

    if (threads_per_group == 0_int64) threads_per_group = 36_int64
    total_sites = ldim * ldim * ldim * ldim

    allocate(a_real(total_sites,4,3,3), a_imag(total_sites,4,3,3))
    allocate(c_real(total_sites,4,3,3), c_imag(total_sites,4,3,3))
    allocate(b_real(4,3,3), b_imag(4,3,3))

    a_real = 1.0_real64
    a_imag = 0.0_real64
    b_real = 1.0_real64 / 3.0_real64
    b_imag = 0.0_real64
    c_real = 0.0_real64
    c_imag = 0.0_real64

    if (verbose >= 1_int32) then
      write(*,'(A,I0,A)') 'Number of sites = ', ldim, '^4'
      write(*,'(A,I0,A,I0,A)') 'Executing ', iterations, ' iterations with ', warmups, ' warmups'
      if (threads_per_group /= 0_int64) write(*,'(A,I0)') 'Threads per group = ', threads_per_group
    end if

    ttotal = su3_mat_nn(a_real, a_imag, b_real, b_imag, c_real, c_imag, total_sites, iterations, &
      threads_per_group, warmups, verbose, checksum)

    if (verbose >= 1_int32) write(*,'(A,F0.6,A)') 'Total kernel execution time = ', ttotal, ' (s)'
    tflop = real(iterations, real64) * real(total_sites, real64) * 864.0_real64
    write(*,'(A,F0.3)') 'Total GFLOP/s = ', tflop / ttotal / 1.0e9_real64
    memory_usage = real(2_int64 * total_sites, real64) * 4.0_real64 * 3.0_real64 * 3.0_real64 * 2.0_real64 * &
      real(storage_size(1.0_real64) / 8, real64) + 4.0_real64 * 3.0_real64 * 3.0_real64 * 2.0_real64 * &
      real(storage_size(1.0_real64) / 8, real64)
    write(*,'(A,F0.3)') 'Total GByte/s (GPU memory)  = ', &
      real(iterations, real64) * memory_usage / ttotal / 1.0e9_real64

    if (abs(checksum - 36.0_real64) > 1.0e-6_real64) error stop 'verification failed'

    deallocate(b_imag, b_real, c_imag, c_real, a_imag, a_real)
  end subroutine run_benchmark

  real(real64) function su3_mat_nn(a_real, a_imag, b_real, b_imag, c_real, c_imag, total_sites, iterations, &
      threads_per_group, warmups, verbose, checksum) result(ttotal)
    real(real64), intent(in) :: a_real(:,:,:,:), a_imag(:,:,:,:), b_real(:,:,:), b_imag(:,:,:)
    real(real64), intent(inout) :: c_real(:,:,:,:), c_imag(:,:,:,:)
    integer(int64), intent(in) :: total_sites, iterations, threads_per_group, warmups
    integer(int32), intent(in) :: verbose
    real(real64), intent(out) :: checksum
    integer(int64) :: num_work_items, id, i, j, k, l, m, iters
    real(real64) :: cc_real, cc_imag, tstart

    num_work_items = total_sites * threads_per_group
    if (verbose >= 1_int32) then
      write(*,'(A,I0)') 'Number of teams = ', total_sites
      write(*,'(A,I0)') 'Threads per team = ', threads_per_group
      write(*,'(A,I0)') 'Number of work items = ', num_work_items
    end if

    !$omp target data map(to: a_real(1:total_sites,1:4,1:3,1:3), a_imag(1:total_sites,1:4,1:3,1:3), &
    !$omp& b_real(1:4,1:3,1:3), b_imag(1:4,1:3,1:3)) &
    !$omp& map(from: c_real(1:total_sites,1:4,1:3,1:3), c_imag(1:total_sites,1:4,1:3,1:3))
    do iters = 0_int64, iterations + warmups - 1_int64
      if (iters == warmups) tstart = omp_get_wtime()
      !$omp target teams distribute parallel do num_teams(total_sites) thread_limit(threads_per_group) &
      !$omp& private(i,j,k,l,m,cc_real,cc_imag)
      do id = 0_int64, num_work_items - 1_int64
        i = id / 36_int64 + 1_int64
        if (i <= total_sites) then
          j = mod(id, 36_int64) / 9_int64 + 1_int64
          k = mod(id, 9_int64) / 3_int64 + 1_int64
          l = mod(id, 3_int64) + 1_int64
          cc_real = 0.0_real64
          cc_imag = 0.0_real64
          do m = 1_int64, 3_int64
            cc_real = cc_real + a_real(i,j,k,m) * b_real(j,m,l) - a_imag(i,j,k,m) * b_imag(j,m,l)
            cc_imag = cc_imag + a_real(i,j,k,m) * b_imag(j,m,l) + a_imag(i,j,k,m) * b_real(j,m,l)
          end do
          c_real(i,j,k,l) = cc_real
          c_imag(i,j,k,l) = cc_imag
        end if
      end do
      !$omp end target teams distribute parallel do
    end do
    !$omp end target data

    ttotal = omp_get_wtime() - tstart

    checksum = 0.0_real64
    do i = 1_int64, total_sites
      do j = 1_int64, 4_int64
        do k = 1_int64, 3_int64
          do l = 1_int64, 3_int64
            cc_real = 0.0_real64
            cc_imag = 0.0_real64
            do m = 1_int64, 3_int64
              cc_real = cc_real + a_real(i,j,k,m) * b_real(j,m,l) - a_imag(i,j,k,m) * b_imag(j,m,l)
              cc_imag = cc_imag + a_real(i,j,k,m) * b_imag(j,m,l) + a_imag(i,j,k,m) * b_real(j,m,l)
            end do
            if (abs(c_real(i,j,k,l) - cc_real) > 1.0e-6_real64 .or. &
                abs(c_imag(i,j,k,l) - cc_imag) > 1.0e-6_real64) error stop 'verification failed'
            checksum = checksum + c_real(i,j,k,l)
          end do
        end do
      end do
    end do
    checksum = checksum / real(total_sites, real64)

    if (abs(checksum - 36.0_real64) < 1.0e-6_real64) then
      write(*,'(A,I0)') 'Checksum SUCCESS... though please be diligent and check the following value is not NaN: checksum=', nint(checksum)
    else
      write(*,'(A)') 'Checksum FAILURE'
    end if
  end function su3_mat_nn

end program main
