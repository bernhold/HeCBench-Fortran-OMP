program main
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use omp_lib
  implicit none

  integer, parameter :: d1 = 41, d2 = 13, d3 = 11, d4 = 9, d5 = 76, d6 = 50
  integer(int64), parameter :: data_size = int(d1, int64) * int(d2, int64) * &
    int(d3, int64) * int(d4, int64) * int(d5, int64) * int(d6, int64)
  character(len=256) :: arg0, arg1
  integer :: repeat, iter
  integer(int64) :: idx
  real(real64), allocatable :: input(:), output(:)
  real(real64) :: start_time, end_time

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) repeat
  if (repeat <= 0) stop 1

  allocate(input(0:data_size - 1), output(0:data_size - 1))
  !$omp parallel do
  do idx = 0, data_size - 1
    input(idx) = real(idx, real64)
    output(idx) = 0.0_real64
  end do
  !$omp end parallel do

  !$omp target data map(to: input(0:data_size - 1)) map(from: output(0:data_size - 1))
  start_time = omp_get_wtime()
  do iter = 1, repeat
    call tensor_transpose(input, output)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', &
    (end_time - start_time) * 1.0e3_real64 / real(repeat, real64), ' (ms)'
  !$omp end target data

  call verify(input, output)
  deallocate(input, output)

contains

  subroutine tensor_transpose(input, output)
    real(real64), intent(in) :: input(0:)
    real(real64), intent(out) :: output(0:)
    integer(int64) :: idx, it
    integer :: i1, i2, i3, i4, i5, i6
    integer(int64) :: out_idx

    !$omp target teams distribute parallel do private(idx, it, i1, i2, i3, i4, i5, i6, out_idx)
    do idx = 0, data_size - 1
      it = idx
      i1 = int(mod(it, int(d1, int64)))
      it = it / int(d1, int64)
      i2 = int(mod(it, int(d2, int64)))
      it = it / int(d2, int64)
      i3 = int(mod(it, int(d3, int64)))
      it = it / int(d3, int64)
      i4 = int(mod(it, int(d4, int64)))
      it = it / int(d4, int64)
      i5 = int(mod(it, int(d5, int64)))
      i6 = int(it / int(d5, int64))

      out_idx = int(i2, int64) + int(d2, int64) * (int(i3, int64) + &
        int(d3, int64) * (int(i4, int64) + int(d4, int64) * (int(i6, int64) + &
        int(d6, int64) * (int(i1, int64) + int(d1, int64) * int(i5, int64)))))
      output(out_idx) = input(idx)
    end do
    !$omp end target teams distribute parallel do
  end subroutine tensor_transpose

  subroutine verify(input, output)
    real(real64), intent(in) :: input(0:), output(0:)
    integer(int64) :: input_offset, output_offset, step_input, step_output
    integer :: i
    logical :: error

    input_offset = 2_int64 + int(d1, int64) * (2_int64 + int(d2, int64) * &
      (2_int64 + int(d3, int64) * (2_int64 + int(d4, int64) * &
      (0_int64 + 2_int64 * int(d5, int64)))))
    output_offset = 2_int64 + int(d2, int64) * (2_int64 + int(d3, int64) * &
      (2_int64 + int(d4, int64) * (2_int64 + int(d6, int64) * &
      (2_int64 + 0_int64 * int(d1, int64)))))
    step_input = int(d1, int64) * int(d2, int64) * int(d3, int64) * int(d4, int64)
    step_output = int(d2, int64) * int(d3, int64) * int(d4, int64) * int(d6, int64) * int(d1, int64)

    error = .false.
    do i = 0, d5 - 1
      if (input(input_offset + int(i, int64) * step_input) /= &
          output(output_offset + int(i, int64) * step_output)) then
        write(*,'(A)') 'FAIL'
        error = .true.
        exit
      end if
    end do
    if (.not. error) write(*,'(A)') 'PASS'
  end subroutine verify

end program main
