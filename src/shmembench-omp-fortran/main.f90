program main
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  integer, parameter :: vector_size = 1024 * 1024
  integer, parameter :: total_iterations = 1024
  integer, parameter :: block_size = 256
  character(len=256) :: arg0, arg
  integer :: repeat, i, iter
  real(real64), allocatable :: c(:)
  real(real64) :: start_time, time_ms, time_ns
  integer(8) :: datasize, operations_bytes, operations_128bit

  call get_command_argument(0, arg0)
  write(*,'(A)') 'Shared memory bandwidth microbenchmark'

  if (command_argument_count() /= 1) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) repeat
  if (repeat <= 0) stop 1

  datasize = int(vector_size, 8) * 8_8
  write(*,'(A,I0,A)') 'Buffer sizes: ', datasize / (1024_8 * 1024_8), 'MB'

  allocate(c(vector_size))
  c = 0.0_real64

  !$omp target data map(from: c(1:vector_size))
  start_time = omp_get_wtime()
  do iter = 1, repeat
    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 1, vector_size
      c(i) = real(i + total_iterations, real64)
    end do
    !$omp end target teams distribute parallel do
  end do
  time_ns = ((omp_get_wtime() - start_time) * 1.0e9_real64) / real(repeat, real64)
  !$omp end target data

  time_ms = time_ns * 1.0e-6_real64
  write(*,'(A,F0.6,A)') 'Average kernel execution time : ', time_ms, ' (ms)'
  write(*,'(A)') 'Memory throughput'

  operations_bytes = (6_8 + 4_8 * 5_8 * int(total_iterations, 8) + 6_8) * &
    int(vector_size, 8) * 4_8
  operations_128bit = (6_8 + 4_8 * 5_8 * int(total_iterations, 8) + 6_8) * &
    int(vector_size, 8) / 4_8

  write(*,'(A,F12.2,A,F10.2,A)') achar(9)//'using 128bit operations : ', &
    real(operations_bytes, real64) / time_ns, ' GB/sec (', &
    real(operations_128bit, real64) / time_ns, ' billion accesses/sec)'

  deallocate(c)
end program main
