program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_sizes(6) = [32, 64, 128, 256, 512, 1024]
  real(real32), parameter :: eps = 1.0e-5_real32

  character(len=256) :: arg0, arg
  integer :: rows, cols, repeat, i, block_size
  integer(int64) :: total_size, memory_ops
  integer(int32) :: seed
  real(real32), allocatable :: inp(:), gamma(:), out(:), d_out(:)
  real(real64) :: start_time, elapsed_ms, bandwidth

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <number of rows> <number of columns> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) rows
  call get_command_argument(2, arg); read(arg, *) cols
  call get_command_argument(3, arg); read(arg, *) repeat
  if (rows <= 0 .or. cols <= 0 .or. repeat <= 0) stop 1

  total_size = int(rows, int64) * int(cols, int64)
  allocate(inp(total_size), gamma(cols), out(total_size), d_out(total_size))

  seed = 0_int32
  call fill_random(inp, seed)
  call fill_random(gamma, seed)

  !$omp target data map(to: inp(1:total_size), gamma(1:cols)) map(tofrom: d_out(1:total_size))
  call rmsnorm_forward_cpu(out, inp, gamma, rows, cols)

  do i = 1, size(block_sizes)
    block_size = block_sizes(i)
    write(*,'(A,I0,A)') 'Checking block size ', block_size, '.'

    call rmsnorm_forward_kernel(inp, gamma, d_out, rows, cols, block_size)
    !$omp target update from(d_out(1:total_size))
    call validate_result(d_out, out, 'out', total_size, 1.0e-5_real32)

    call rmsnorm_forward_kernel2(inp, gamma, d_out, rows, cols, block_size)
    !$omp target update from(d_out(1:total_size))
    call validate_result(d_out, out, 'out', total_size, 1.0e-5_real32)
  end do

  write(*,'(A)') 'All results match. Starting benchmarks.'
  write(*,*)

  do i = 1, size(block_sizes)
    block_size = block_sizes(i)
    start_time = omp_get_wtime()
    call run_repeated(1, inp, gamma, d_out, rows, cols, block_size, repeat)
    elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64 / real(repeat, real64)
    memory_ops = (2_int64 * total_size + int(cols, int64)) * 4_int64
    bandwidth = real(memory_ops, real64) / elapsed_ms / 1.0e6_real64
    write(*,'(A,I4,A,F0.4,A,F0.2,A)') 'block_size ', block_size, ' | time ', elapsed_ms, &
      ' ms | bandwidth ', bandwidth, ' GB/s'
  end do

  write(*,*)
  do i = 1, size(block_sizes)
    block_size = block_sizes(i)
    start_time = omp_get_wtime()
    call run_repeated(2, inp, gamma, d_out, rows, cols, block_size, repeat)
    elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64 / real(repeat, real64)
    memory_ops = (2_int64 * total_size + int(cols, int64)) * 4_int64
    bandwidth = real(memory_ops, real64) / elapsed_ms / 1.0e6_real64
    write(*,'(A,I4,A,F0.4,A,F0.2,A)') 'block_size ', block_size, ' | time ', elapsed_ms, &
      ' ms | bandwidth ', bandwidth, ' GB/s'
  end do
  !$omp end target data

  deallocate(inp, gamma, out, d_out)

contains

  subroutine fill_random(values, seed)
    real(real32), intent(out) :: values(:)
    integer(int32), intent(inout) :: seed
    integer :: i

    do i = 1, size(values)
      values(i) = real(c_rand(seed), real32) / 32767.0_real32 * 2.0_real32 - 1.0_real32
    end do
  end subroutine fill_random

  subroutine rmsnorm_forward_kernel(inp, gamma, out, rows, cols, block_size)
    real(real32), intent(in) :: inp(:), gamma(:)
    real(real32), intent(out) :: out(:)
    integer, intent(in) :: rows, cols, block_size
    integer :: t, j, base
    real(real32) :: m, s

    !$omp target teams distribute parallel do private(j, base, m, s) thread_limit(block_size)
    do t = 0, rows - 1
      base = t * cols
      m = 0.0_real32
      do j = 1, cols
        m = m + inp(base + j) * inp(base + j)
      end do
      m = m / real(cols, real32)
      s = 1.0_real32 / sqrt(m + eps)
      do j = 1, cols
        out(base + j) = inp(base + j) * s * gamma(j)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine rmsnorm_forward_kernel

  subroutine rmsnorm_forward_kernel2(inp, gamma, out, rows, cols, block_size)
    real(real32), intent(in) :: inp(:), gamma(:)
    real(real32), intent(out) :: out(:)
    integer, intent(in) :: rows, cols, block_size
    integer :: t, j, base
    real(real32) :: m, s

    !$omp target teams distribute parallel do private(j, base, m, s) thread_limit(block_size)
    do t = 0, rows - 1
      base = t * cols
      m = 0.0_real32
      do j = 1, cols
        m = m + inp(base + j) * inp(base + j)
      end do
      s = 1.0_real32 / sqrt(m / real(cols, real32) + eps)
      do j = 1, cols
        out(base + j) = inp(base + j) * s * gamma(j)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine rmsnorm_forward_kernel2

  subroutine run_repeated(which_kernel, inp, gamma, out, rows, cols, block_size, repeat)
    integer, intent(in) :: which_kernel, rows, cols, block_size, repeat
    real(real32), intent(in) :: inp(:), gamma(:)
    real(real32), intent(out) :: out(:)
    integer :: iter

    do iter = 1, repeat
      if (which_kernel == 1) then
        call rmsnorm_forward_kernel(inp, gamma, out, rows, cols, block_size)
      else
        call rmsnorm_forward_kernel2(inp, gamma, out, rows, cols, block_size)
      end if
    end do
  end subroutine run_repeated

  subroutine rmsnorm_forward_cpu(out, inp, gamma, rows, cols)
    real(real32), intent(out) :: out(:)
    real(real32), intent(in) :: inp(:), gamma(:)
    integer, intent(in) :: rows, cols
    integer :: t, j, base
    real(real32) :: m, s

    do t = 0, rows - 1
      base = t * cols
      m = 0.0_real32
      do j = 1, cols
        m = m + inp(base + j) * inp(base + j)
      end do
      m = m / real(cols, real32)
      s = 1.0_real32 / sqrt(m + eps)
      do j = 1, cols
        out(base + j) = inp(base + j) * s * gamma(j)
      end do
    end do
  end subroutine rmsnorm_forward_cpu

  subroutine validate_result(device_result, cpu_reference, name, num_elements, tolerance)
    real(real32), intent(in) :: device_result(:), cpu_reference(:), tolerance
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: num_elements
    integer(int64) :: i
    integer :: nfaults
    real(real32) :: fp_epsilon, t_eff

    nfaults = 0
    fp_epsilon = epsilon(1.0_real32)
    do i = 1, num_elements
      if (cpu_reference(i) /= cpu_reference(i) .or. abs(cpu_reference(i)) > huge(cpu_reference(i))) cycle
      t_eff = tolerance + abs(cpu_reference(i)) * fp_epsilon
      if (abs(cpu_reference(i) - device_result(i)) > t_eff) then
        write(*,'(A,A,A,I0,A,F0.6,A,F0.6)') 'Mismatch of ', trim(name), ' at ', i - 1, &
          ': CPU_ref: ', cpu_reference(i), ' vs GPU: ', device_result(i)
        nfaults = nfaults + 1
        if (nfaults >= 10) stop 1
      end if
    end do
    if (nfaults > 0) stop 1
  end subroutine validate_result

  integer(int32) function c_rand(seed)
    integer(int32), intent(inout) :: seed
    integer(int64) :: next_value

    next_value = mod(1103515245_int64 * int(seed, int64) + 12345_int64, 2147483648_int64)
    seed = int(next_value, int32)
    c_rand = iand(seed / 65536_int32, 32767_int32)
  end function c_rand

end program main
