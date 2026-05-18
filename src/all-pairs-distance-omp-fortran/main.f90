program main
  use, intrinsic :: iso_c_binding, only : c_int, c_long
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use omp_lib
  implicit none

  integer, parameter :: instances = 224, attributes = 4096
  integer :: iterations, i, attr, instance_id
  integer(int32), allocatable :: data(:), cpu_distance(:), gpu_distance(:)
  real(real64) :: start_time, end_time, elapsed_us

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_random() bind(C, name='random') result(value)
      import :: c_long
      integer(c_long) :: value
    end function c_random
  end interface

  if (command_argument_count() /= 1) then
    print '(A)', 'Usage: ./main <iterations>'
    stop 1
  end if
  iterations = read_arg(1)
  if (iterations <= 0) stop 1

  allocate(data(instances * attributes), cpu_distance(instances * instances), gpu_distance(instances * instances))
  call c_srand(2_c_int)
  do attr = 1, attributes
    do instance_id = 1, instances
      data(attr + attributes * (instance_id - 1)) = int(mod(c_random(), 3_c_long), int32)
    end do
  end do

  cpu_distance = 0_int32
  start_time = omp_get_wtime()
  call cpu_reference(data, cpu_distance)
  end_time = omp_get_wtime()
  print '(A,F0.6,A)', 'CPU time: ', (end_time - start_time) * 1.0e6_real64, ' (us)'

  gpu_distance = 0_int32
  !$omp target data map(to: data(1:instances*attributes)) map(alloc: gpu_distance(1:instances*instances))
  elapsed_us = 0.0_real64
  do i = 1, iterations
    gpu_distance = 0_int32
    !$omp target update to(gpu_distance(1:instances*instances))
    start_time = omp_get_wtime()
    call distance_device(data, gpu_distance)
    end_time = omp_get_wtime()
    elapsed_us = elapsed_us + (end_time - start_time) * 1.0e6_real64
  end do
  !$omp target update from(gpu_distance(1:instances*instances))
  print '(A,F0.6,A)', 'Average kernel execution time (w/o shared memory): ', elapsed_us / real(iterations, real64), ' (us)'
  call print_status(cpu_distance, gpu_distance)

  elapsed_us = 0.0_real64
  do i = 1, iterations
    gpu_distance = 0_int32
    !$omp target update to(gpu_distance(1:instances*instances))
    start_time = omp_get_wtime()
    call distance_device(data, gpu_distance)
    end_time = omp_get_wtime()
    elapsed_us = elapsed_us + (end_time - start_time) * 1.0e6_real64
  end do
  !$omp target update from(gpu_distance(1:instances*instances))
  print '(A,F0.6,A)', 'Average kernel execution time (w/ shared memory): ', elapsed_us / real(iterations, real64), ' (us)'
  call print_status(cpu_distance, gpu_distance)
  !$omp end target data

  deallocate(data, cpu_distance, gpu_distance)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine cpu_reference(data, distance)
    integer(int32), intent(in) :: data(:)
    integer(int32), intent(inout) :: distance(:)
    integer :: i, j, k, count
    !$omp parallel do collapse(2) private(k, count)
    do j = 1, instances
      do i = 1, instances
        count = 0
        do k = 1, attributes
          if (data((i - 1) * attributes + k) /= data((j - 1) * attributes + k)) count = count + 1
        end do
        distance(i + instances * (j - 1)) = count
      end do
    end do
    !$omp end parallel do
  end subroutine cpu_reference

  subroutine distance_device(data, distance)
    integer(int32), intent(in) :: data(:)
    integer(int32), intent(out) :: distance(:)
    integer :: i, j, k, count
    !$omp target teams distribute parallel do collapse(2) thread_limit(128) private(k, count)
    do j = 1, instances
      do i = 1, instances
        count = 0
        do k = 1, attributes
          if (data((i - 1) * attributes + k) /= data((j - 1) * attributes + k)) count = count + 1
        end do
        distance(i + instances * (j - 1)) = count
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine distance_device

  subroutine print_status(expected, actual)
    integer(int32), intent(in) :: expected(:), actual(:)
    integer :: i
    logical :: ok
    ok = .true.
    do i = 1, instances * instances
      if (expected(i) /= actual(i)) then
        ok = .false.
        exit
      end if
    end do
    if (ok) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
      stop 1
    end if
  end subroutine print_status

end program main
