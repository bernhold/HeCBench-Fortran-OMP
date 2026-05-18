program main
  use, intrinsic :: iso_c_binding, only : c_double, c_long
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  integer, parameter :: num_floats = 2 * 1024 * 1024
  integer :: repeat

  interface
    subroutine c_srand48(seed) bind(C, name="srand48")
      import :: c_long
      integer(c_long), value :: seed
    end subroutine c_srand48

    real(c_double) function c_drand48() bind(C, name="drand48")
      import :: c_double
    end function c_drand48
  end interface

  if (command_argument_count() /= 1) then
    call print_usage()
    stop 1
  end if

  repeat = read_int_arg(1)
  if (repeat < 0) stop 1

  write(*,'(A)') '=== Single-precision floating-point kernels ==='
  call test_sp(repeat, num_floats)
  write(*,'(A)') '=== Double-precision floating-point kernels ==='
  call test_dp(repeat, num_floats)
  write(*,'(A)') 'PASS'

contains

  subroutine print_usage()
    character(len=256) :: arg0
    call get_command_argument(0, arg0)
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <repeat>'
  end subroutine print_usage

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  subroutine test_sp(repeat, n)
    integer, intent(in) :: repeat, n
    real(real32), allocatable :: data(:)
    integer :: j
    real(real64) :: t0, elapsed
    allocate(data(n))
    call c_srand48(123_c_long)
    do j = 1, n / 2
      data(j) = real(c_drand48() * 10.0_c_double, real32)
      data(n - j + 1) = data(j)
    end do
    !$omp target data map(alloc: data(1:n))
    do j = 1, 4
      call kernel_sp(data, n, repeat, 1, 1)
      call kernel_sp(data, n, repeat, 1, 2)
      call kernel_sp(data, n, repeat, 1, 4)
      call kernel_sp(data, n, repeat, 1, 8)
    end do
    call timed_sp(data, n, repeat, 1, 1, 'Add1')
    call timed_sp(data, n, repeat, 1, 2, 'Add2')
    call timed_sp(data, n, repeat, 1, 4, 'Add4')
    call timed_sp(data, n, repeat, 1, 8, 'Add8')
    do j = 1, 4
      call kernel_sp(data, n, repeat, 2, 1)
      call kernel_sp(data, n, repeat, 2, 2)
      call kernel_sp(data, n, repeat, 2, 4)
      call kernel_sp(data, n, repeat, 2, 8)
    end do
    call timed_sp(data, n, repeat, 2, 1, 'Mul1')
    call timed_sp(data, n, repeat, 2, 2, 'Mul2')
    call timed_sp(data, n, repeat, 2, 4, 'Mul4')
    call timed_sp(data, n, repeat, 2, 8, 'Mul8')
    do j = 1, 4
      call kernel_sp(data, n, repeat, 3, 1)
      call kernel_sp(data, n, repeat, 3, 2)
      call kernel_sp(data, n, repeat, 3, 4)
      call kernel_sp(data, n, repeat, 3, 8)
    end do
    call timed_sp(data, n, repeat, 3, 1, 'MAdd1')
    call timed_sp(data, n, repeat, 3, 2, 'MAdd2')
    call timed_sp(data, n, repeat, 3, 4, 'MAdd4')
    call timed_sp(data, n, repeat, 3, 8, 'MAdd8')
    do j = 1, 4
      call kernel_sp(data, n, repeat, 4, 1)
      call kernel_sp(data, n, repeat, 4, 2)
      call kernel_sp(data, n, repeat, 4, 4)
      call kernel_sp(data, n, repeat, 4, 8)
    end do
    call timed_sp(data, n, repeat, 4, 1, 'MulMAdd1')
    call timed_sp(data, n, repeat, 4, 2, 'MulMAdd2')
    call timed_sp(data, n, repeat, 4, 4, 'MulMAdd4')
    call timed_sp(data, n, repeat, 4, 8, 'MulMAdd8')
    !$omp end target data
    deallocate(data)
  end subroutine test_sp

  subroutine timed_sp(data, n, repeat, op, lanes, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat, op, lanes
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call kernel_sp(data, n, repeat, op, lanes)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_sp

  subroutine kernel_sp(data, n, repeat, op, lanes)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat, op, lanes
    integer :: gid, j, k, reps
    real(real32) :: s, s2, s3, s4, s5, s6, s7, s8
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,reps,s,s2,s3,s4,s5,s6,s7,s8)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real32 - s
      s3 = 9.0_real32 - s
      s4 = 9.0_real32 - s2
      s5 = 8.0_real32 - s
      s6 = 8.0_real32 - s2
      s7 = 7.0_real32 - s
      s8 = 7.0_real32 - s2
      do j = 1, repeat
        select case (op)
        case (1)
          reps = merge(240, merge(120, merge(60, 30, lanes == 4), lanes == 2), lanes == 1)
          do k = 1, reps
            s = 10.0_real32 - s
            if (lanes >= 2) s2 = 10.0_real32 - s2
            if (lanes >= 4) then
              s3 = 10.0_real32 - s3
              s4 = 10.0_real32 - s4
            end if
            if (lanes >= 8) then
              s5 = 10.0_real32 - s5
              s6 = 10.0_real32 - s6
              s7 = 10.0_real32 - s7
              s8 = 10.0_real32 - s8
            end if
          end do
        case (2)
          reps = merge(200, merge(100, merge(50, 25, lanes == 4), lanes == 2), lanes == 1)
          s = data(gid) - data(gid) + 0.999_real32
          s2 = s - 0.0001_real32
          s3 = s - 0.0002_real32
          s4 = s - 0.0003_real32
          s5 = s - 0.0004_real32
          s6 = s - 0.0005_real32
          s7 = s - 0.0006_real32
          s8 = s - 0.0007_real32
          do k = 1, reps
            s = s * s * 1.01_real32
            if (lanes >= 2) s2 = s2 * s2 * 1.01_real32
            if (lanes >= 4) then
              s3 = s3 * s3 * 1.01_real32
              s4 = s4 * s4 * 1.01_real32
            end if
            if (lanes >= 8) then
              s5 = s5 * s5 * 1.01_real32
              s6 = s6 * s6 * 1.01_real32
              s7 = s7 * s7 * 1.01_real32
              s8 = s8 * s8 * 1.01_real32
            end if
          end do
        case (3)
          reps = merge(240, merge(120, merge(60, 30, lanes == 4), lanes == 2), lanes == 1)
          do k = 1, reps
            s = 10.0_real32 - s * 0.9899_real32
            if (lanes >= 2) s2 = 10.0_real32 - s2 * 0.9899_real32
            if (lanes >= 4) then
              s3 = 10.0_real32 - s3 * 0.9899_real32
              s4 = 10.0_real32 - s4 * 0.9899_real32
            end if
            if (lanes >= 8) then
              s5 = 10.0_real32 - s5 * 0.9899_real32
              s6 = 10.0_real32 - s6 * 0.9899_real32
              s7 = 10.0_real32 - s7 * 0.9899_real32
              s8 = 10.0_real32 - s8 * 0.9899_real32
            end if
          end do
        case (4)
          reps = merge(240, merge(120, merge(60, 30, lanes == 4), lanes == 2), lanes == 1)
          do k = 1, reps
            s = (3.75_real32 - 0.355_real32 * s) * s
            if (lanes >= 2) s2 = (3.75_real32 - 0.355_real32 * s2) * s2
            if (lanes >= 4) then
              s3 = (3.75_real32 - 0.355_real32 * s3) * s3
              s4 = (3.75_real32 - 0.355_real32 * s4) * s4
            end if
            if (lanes >= 8) then
              s5 = (3.75_real32 - 0.355_real32 * s5) * s5
              s6 = (3.75_real32 - 0.355_real32 * s6) * s6
              s7 = (3.75_real32 - 0.355_real32 * s7) * s7
              s8 = (3.75_real32 - 0.355_real32 * s8) * s8
            end if
          end do
        end select
      end do
      data(gid) = merge(s, merge(s + s2, merge((s+s2)+(s3+s4), ((s+s2)+(s3+s4))+((s5+s6)+(s7+s8)), lanes == 4), lanes == 2), lanes == 1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel_sp

  subroutine test_dp(repeat, n)
    integer, intent(in) :: repeat, n
    real(real64), allocatable :: data(:)
    integer :: j
    allocate(data(n))
    call c_srand48(123_c_long)
    do j = 1, n / 2
      data(j) = real(c_drand48() * 10.0_c_double, real64)
      data(n - j + 1) = data(j)
    end do
    !$omp target data map(alloc: data(1:n))
    call timed_dp_set(data, n, repeat)
    !$omp end target data
    deallocate(data)
  end subroutine test_dp

  subroutine timed_dp_set(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: op, lanes_idx
    character(len=8), parameter :: op_names(4) = [character(len=8) :: 'Add', 'Mul', 'MAdd', 'MulMAdd']
    integer, parameter :: lane_values(4) = [1, 2, 4, 8]
    do op = 1, 4
      do lanes_idx = 1, 4
        call timed_dp(data, n, repeat, op, lane_values(lanes_idx), trim(op_names(op)) // trim(int_to_string(lane_values(lanes_idx))))
      end do
    end do
  end subroutine timed_dp_set

  subroutine timed_dp(data, n, repeat, op, lanes, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat, op, lanes
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call kernel_dp(data, n, repeat, op, lanes)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_dp

  subroutine kernel_dp(data, n, repeat, op, lanes)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat, op, lanes
    integer :: gid, j, k, reps
    real(real64) :: s, s2, s3, s4, s5, s6, s7, s8
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,reps,s,s2,s3,s4,s5,s6,s7,s8)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real64 - s
      s3 = 9.0_real64 - s
      s4 = 9.0_real64 - s2
      s5 = 8.0_real64 - s
      s6 = 8.0_real64 - s2
      s7 = 7.0_real64 - s
      s8 = 7.0_real64 - s2
      do j = 1, repeat
        reps = merge(240, merge(120, merge(60, 30, lanes == 4), lanes == 2), lanes == 1)
        if (op == 2) reps = merge(200, merge(100, merge(50, 25, lanes == 4), lanes == 2), lanes == 1)
        do k = 1, reps
          select case (op)
          case (1)
            s = 10.0_real64 - s
            if (lanes >= 2) s2 = 10.0_real64 - s2
            if (lanes >= 4) then
              s3 = 10.0_real64 - s3
              s4 = 10.0_real64 - s4
            end if
            if (lanes >= 8) then
              s5 = 10.0_real64 - s5; s6 = 10.0_real64 - s6
              s7 = 10.0_real64 - s7; s8 = 10.0_real64 - s8
            end if
          case (2)
            s = s * s * 1.01_real64
            if (lanes >= 2) s2 = s2 * s2 * 1.01_real64
            if (lanes >= 4) then
              s3 = s3 * s3 * 1.01_real64
              s4 = s4 * s4 * 1.01_real64
            end if
            if (lanes >= 8) then
              s5 = s5 * s5 * 1.01_real64; s6 = s6 * s6 * 1.01_real64
              s7 = s7 * s7 * 1.01_real64; s8 = s8 * s8 * 1.01_real64
            end if
          case (3)
            s = 10.0_real64 - s * 0.9899_real64
            if (lanes >= 2) s2 = 10.0_real64 - s2 * 0.9899_real64
            if (lanes >= 4) then
              s3 = 10.0_real64 - s3 * 0.9899_real64
              s4 = 10.0_real64 - s4 * 0.9899_real64
            end if
            if (lanes >= 8) then
              s5 = 10.0_real64 - s5 * 0.9899_real64; s6 = 10.0_real64 - s6 * 0.9899_real64
              s7 = 10.0_real64 - s7 * 0.9899_real64; s8 = 10.0_real64 - s8 * 0.9899_real64
            end if
          case (4)
            s = (3.75_real64 - 0.355_real64 * s) * s
            if (lanes >= 2) s2 = (3.75_real64 - 0.355_real64 * s2) * s2
            if (lanes >= 4) then
              s3 = (3.75_real64 - 0.355_real64 * s3) * s3
              s4 = (3.75_real64 - 0.355_real64 * s4) * s4
            end if
            if (lanes >= 8) then
              s5 = (3.75_real64 - 0.355_real64 * s5) * s5; s6 = (3.75_real64 - 0.355_real64 * s6) * s6
              s7 = (3.75_real64 - 0.355_real64 * s7) * s7; s8 = (3.75_real64 - 0.355_real64 * s8) * s8
            end if
          end select
        end do
      end do
      data(gid) = merge(s, merge(s + s2, merge((s+s2)+(s3+s4), ((s+s2)+(s3+s4))+((s5+s6)+(s7+s8)), lanes == 4), lanes == 2), lanes == 1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel_dp

  character(len=2) function int_to_string(value) result(text)
    integer, intent(in) :: value
    write(text,'(I0)') value
  end function int_to_string

end program main
