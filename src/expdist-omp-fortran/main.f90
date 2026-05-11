program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2
  integer :: size, repeat

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage ./', trim(arg0), ' <size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *) size
  read(arg2, *) repeat
  if (size <= 0 .or. repeat <= 0) stop 1

  print '(A)', 'Test single precision'
  call test_real32(size, repeat)

  print '(A)', 'Test double precision'
  call test_real64(size, repeat)

contains

  subroutine test_real32(size, repeat)
    integer, intent(in) :: size, repeat
    integer :: i
    real(real32), allocatable :: a(:), b(:), scale_a(:), scale_b(:)
    real(real32) :: output, host_output
    real(real64) :: start_time, end_time, avg_time, analytical

    allocate(a(size * 2), b(size * 2), scale_a(size), scale_b(size))
    a = 1.0_real32
    b = 0.0_real32
    scale_a = 1.0_real32
    scale_b = 1.0_real32

    !$omp target data map(to: a(1:size*2), b(1:size*2), scale_a(1:size), scale_b(1:size))
    start_time = omp_get_wtime()
    do i = 1, repeat
      output = distance_real32(a, b, scale_a, scale_b, size)
    end do
    end_time = omp_get_wtime()
    avg_time = (end_time - start_time) / real(repeat, real64)
    print '(A,F8.6,A)', 'Average kernel execution time ', avg_time, ' (s)'
    print '(A,F0.6)', '    device result: ', real(output, real64)

    host_output = host_cost_real32(a, b, scale_a, scale_b, size)
    print '(A,F0.6)', '      host result: ', real(host_output, real64)
    analytical = real(size, real64) * real(size, real64) * exp(-1.0_real64)
    print '(A,F0.6)', 'analytical result: ', analytical
    print '(A)'
    !$omp end target data

    deallocate(a, b, scale_a, scale_b)
  end subroutine test_real32

  subroutine test_real64(size, repeat)
    integer, intent(in) :: size, repeat
    integer :: i
    real(real64), allocatable :: a(:), b(:), scale_a(:), scale_b(:)
    real(real64) :: output, host_output, start_time, end_time, avg_time, analytical

    allocate(a(size * 2), b(size * 2), scale_a(size), scale_b(size))
    a = 1.0_real64
    b = 0.0_real64
    scale_a = 1.0_real64
    scale_b = 1.0_real64

    !$omp target data map(to: a(1:size*2), b(1:size*2), scale_a(1:size), scale_b(1:size))
    start_time = omp_get_wtime()
    do i = 1, repeat
      output = distance_real64(a, b, scale_a, scale_b, size)
    end do
    end_time = omp_get_wtime()
    avg_time = (end_time - start_time) / real(repeat, real64)
    print '(A,F8.6,A)', 'Average kernel execution time ', avg_time, ' (s)'
    print '(A,F0.6)', '    device result: ', output

    host_output = host_cost_real64(a, b, scale_a, scale_b, size)
    print '(A,F0.6)', '      host result: ', host_output
    analytical = real(size, real64) * real(size, real64) * exp(-1.0_real64)
    print '(A,F0.6)', 'analytical result: ', analytical
    print '(A)'
    !$omp end target data

    deallocate(a, b, scale_a, scale_b)
  end subroutine test_real64

  function distance_real32(a, b, scale_a, scale_b, size) result(total)
    real(real32), intent(in) :: a(:), b(:), scale_a(:), scale_b(:)
    integer, intent(in) :: size
    real(real32) :: total
    integer :: i, j
    real(real32) :: dist

    total = 0.0_real32
    !$omp target teams distribute parallel do collapse(2) reduction(+:total) private(dist) thread_limit(256)
    do j = 1, size
      do i = 1, size
        dist = (a(i) - b(j)) * (a(i) - b(j)) + &
               (a(i + size) - b(j + size)) * (a(i + size) - b(j + size))
        total = total + exp(-dist / (scale_a(i) + scale_b(j)))
      end do
    end do
    !$omp end target teams distribute parallel do
  end function distance_real32

  function distance_real64(a, b, scale_a, scale_b, size) result(total)
    real(real64), intent(in) :: a(:), b(:), scale_a(:), scale_b(:)
    integer, intent(in) :: size
    real(real64) :: total
    integer :: i, j
    real(real64) :: dist

    total = 0.0_real64
    !$omp target teams distribute parallel do collapse(2) reduction(+:total) private(dist) thread_limit(256)
    do j = 1, size
      do i = 1, size
        dist = (a(i) - b(j)) * (a(i) - b(j)) + &
               (a(i + size) - b(j + size)) * (a(i + size) - b(j + size))
        total = total + exp(-dist / (scale_a(i) + scale_b(j)))
      end do
    end do
    !$omp end target teams distribute parallel do
  end function distance_real64

  function host_cost_real32(a, b, scale_a, scale_b, size) result(total32)
    real(real32), intent(in) :: a(:), b(:), scale_a(:), scale_b(:)
    integer, intent(in) :: size
    real(real32) :: total32
    real(real64) :: total64
    integer :: i, j
    real(real32) :: dist

    total64 = 0.0_real64
    do i = 1, size
      do j = 1, size
        dist = (a(i) - b(j)) * (a(i) - b(j)) + &
               (a(i + size) - b(j + size)) * (a(i + size) - b(j + size))
        total64 = total64 + real(exp(-dist / (scale_a(i) + scale_b(j))), real64)
      end do
    end do
    total32 = real(total64, real32)
  end function host_cost_real32

  function host_cost_real64(a, b, scale_a, scale_b, size) result(total)
    real(real64), intent(in) :: a(:), b(:), scale_a(:), scale_b(:)
    integer, intent(in) :: size
    real(real64) :: total
    integer :: i, j
    real(real64) :: dist

    total = 0.0_real64
    do i = 1, size
      do j = 1, size
        dist = (a(i) - b(j)) * (a(i) - b(j)) + &
               (a(i + size) - b(j + size)) * (a(i + size) - b(j + size))
        total = total + exp(-dist / (scale_a(i) + scale_b(j)))
      end do
    end do
  end function host_cost_real64

end program main
