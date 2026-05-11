program main
  use, intrinsic :: iso_c_binding, only : c_signed_char
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1
  integer :: repeat
  integer(c_signed_char) :: small_args(16), medium_args(256), large_args(4096)
  integer(c_signed_char) :: sink

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    print '(3A)', 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) repeat
  if (repeat <= 0) stop 1

  small_args = 1_c_signed_char
  medium_args = 1_c_signed_char
  large_args = 1_c_signed_char
  sink = 0_c_signed_char

  call run_case('kernelWithSmallArgs', small_args, repeat, sink)
  call run_case('kernelWithMediumArgs', medium_args, repeat, sink)
  call run_case('kernelWithLargeArgs', large_args, repeat, sink)

contains

  subroutine run_case(label, args, repeat, sink)
    character(len=*), intent(in) :: label
    integer(c_signed_char), intent(in) :: args(:)
    integer, intent(in) :: repeat
    integer(c_signed_char), intent(inout) :: sink
    integer :: i
    real(real64) :: start_time, end_time, elapsed_us

    do i = 1, repeat
      call launch_payload(args, sink)
    end do

    start_time = omp_get_wtime()
    do i = 1, repeat
      call launch_payload(args, sink)
    end do
    end_time = omp_get_wtime()

    elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
    print '(A,A,A,F0.6,A)', 'Average execution time of ', trim(label), ': ', elapsed_us, ' (us)'
  end subroutine run_case

  subroutine launch_payload(args, sink)
    integer(c_signed_char), intent(in) :: args(:)
    integer(c_signed_char), intent(inout) :: sink
    integer :: n

    n = size(args)
    !$omp target map(to: args(1:n)) map(tofrom: sink)
    sink = args(1)
    !$omp end target
  end subroutine launch_payload

end program main
