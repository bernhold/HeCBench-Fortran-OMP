program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name="rand") result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  character(len=256) :: arg0, arg
  integer :: nelem, repeat, i
  integer, allocatable :: grad_in(:), grad_out(:)
  integer(int64) :: grad_in_dim_size
  real(real64) :: start_time, elapsed
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <number of elements> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) nelem
  call get_command_argument(2, arg); read(arg, *) repeat
  if (nelem <= 0 .or. repeat <= 0) stop 1

  grad_in_dim_size = int(nelem, int64)
  allocate(grad_in(nelem), grad_out(nelem))

  call c_srand(123_c_int)
  do i = 1, nelem
    grad_in(i) = modulo(c_rand(), 256_c_int)
  end do
  grad_out = 0

  !$omp target data map(to: grad_in(1:nelem)) map(tofrom: grad_out(1:nelem))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call unfold_backward_internal_kernel(grad_out, grad_in, grad_in_dim_size)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average execution time of unfold backward kernel: ', &
    (elapsed * 1000000.0_real64) / real(repeat, real64), ' (us)'

  ok = .true.
  do i = 1, nelem
    if (repeat * grad_in(i) /= grad_out(i)) then
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(grad_in, grad_out)

contains

  subroutine unfold_backward_internal_kernel(grad_out, grad_in, grad_in_dim_size)
    integer, intent(inout) :: grad_out(:)
    integer, intent(in) :: grad_in(:)
    integer(int64), intent(in) :: grad_in_dim_size
    integer, parameter :: n_threads = 64
    integer, parameter :: n_elems_per_thread = 4
    integer, parameter :: total_work_block = n_threads * n_elems_per_thread
    integer :: grid

    grid = int((grad_in_dim_size + total_work_block - 1_int64) / total_work_block)
    call unfold_backward_elementwise_kernel(grid, int(grad_in_dim_size), grad_out, grad_in)
  end subroutine unfold_backward_internal_kernel

  subroutine unfold_backward_elementwise_kernel(grid, total_n_elems, grad_out, grad_in)
    integer, intent(in) :: grid, total_n_elems
    integer, intent(inout) :: grad_out(:)
    integer, intent(in) :: grad_in(:)
    integer, parameter :: n_threads = 64
    integer, parameter :: n_elems_per_thread = 4
    integer, parameter :: total_work_block = n_threads * n_elems_per_thread
    integer :: bid, tid, elem, j

    !$omp target teams distribute parallel do collapse(2) num_teams(grid) thread_limit(n_threads) &
    !$omp& private(elem, j) map(tofrom: grad_out(1:total_n_elems)) map(to: grad_in(1:total_n_elems))
    do bid = 0, grid - 1
      do tid = 0, n_threads - 1
        elem = bid * total_work_block + tid
        do j = 1, n_elems_per_thread
          if (elem < total_n_elems) then
            grad_out(elem + 1) = grad_out(elem + 1) + grad_in(elem + 1)
            elem = elem + n_threads
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine unfold_backward_elementwise_kernel

end program main
