program threadfence_main
  use, intrinsic :: iso_fortran_env, only: real32
  use omp_lib
  implicit none

  integer :: repeat, n, blocks, grids, iter
  real(real32), allocatable :: h_array(:)
  real(real32) :: result
  real(8) :: elapsed, start_time, end_time
  logical :: ok

  call parse_args(repeat, n)

  blocks = 256
  grids = (n + blocks - 1) / blocks

  allocate(h_array(n))
  h_array = -1.0_real32
  elapsed = 0.0_8
  ok = .true.

  !$omp target data map(to: h_array)
  do iter = 1, repeat
    start_time = omp_get_wtime()
    call sum_device(h_array, n, result)
    end_time = omp_get_wtime()
    elapsed = elapsed + (end_time - start_time)

    if (result /= -1.0_real32 * real(n, real32)) then
      ok = .false.
      exit
    end if
  end do
  !$omp end target data

  if (ok) then
    write(*,'(A,F0.6,A)') "Average kernel execution time: ", &
         real(elapsed * 1.0d3 / dble(repeat), real32), " (ms)"
  end if

  write(*,'(A)') merge("PASS", "FAIL", ok)

contains

  subroutine parse_args(repeat, n)
    integer, intent(out) :: repeat, n
    character(len=64) :: arg
    integer :: status

    if (command_argument_count() /= 2) then
      call get_command_argument(0, arg)
      write(*,'(A,A,A)') "Usage: ", trim(arg), " <repeat> <array length>"
      stop 1
    end if

    call get_command_argument(1, arg)
    read(arg, *, iostat=status) repeat
    if (status /= 0 .or. repeat <= 0) error stop "invalid repeat"

    call get_command_argument(2, arg)
    read(arg, *, iostat=status) n
    if (status /= 0 .or. n <= 0) error stop "invalid array length"
  end subroutine parse_args

  subroutine sum_device(array, n, result)
    integer, intent(in) :: n
    real(real32), intent(in) :: array(n)
    real(real32), intent(out) :: result
    integer :: i
    real(real32) :: total

    total = 0.0_real32
    !$omp target teams distribute parallel do thread_limit(256) reduction(+:total)
    do i = 1, n
      total = total + array(i)
    end do
    !$omp end target teams distribute parallel do
    result = total
  end subroutine sum_device

end program threadfence_main
