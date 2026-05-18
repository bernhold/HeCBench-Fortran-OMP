program nlll
  use iso_c_binding, only: c_int
  use iso_fortran_env, only: int32, int64, real32, real64
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

  integer :: argc
  character(len=64) :: arg
  integer(int64) :: nframe, n_classes
  integer :: repeat

  argc = command_argument_count()
  if (argc /= 3) then
    call get_command_argument(0, arg)
    print '(A,A,A)', 'Usage: ', trim(arg), ' <minibatch size> <number of classes> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) nframe
  call get_command_argument(2, arg)
  read(arg, *) n_classes
  call get_command_argument(3, arg)
  read(arg, *) repeat

  print '(A)', '=========== Data type is FP32 =========='
  call driver(nframe, n_classes, repeat)

contains

  subroutine driver(nframe, n_classes, repeat)
    integer(int64), intent(in) :: nframe, n_classes
    integer, intent(in) :: repeat

    integer(int64) :: input_size, weights_size, target_size
    integer(int64) :: i
    real(real32), allocatable :: input(:), weights(:)
    integer(int32), allocatable :: target(:)
    real(real32) :: r_output, r_total_weight
    logical :: size_average
    integer(int64) :: ignore_index

    input_size = nframe * n_classes
    weights_size = nframe
    target_size = nframe

    allocate(input(input_size), weights(weights_size), target(target_size))

    call c_srand(123_c_int)

    print '(A)', 'Initialization of input data may take a while..'
    do i = 1_int64, input_size
      input(i) = rand_real()
    end do
    do i = 1_int64, weights_size
      weights(i) = rand_real()
    end do
    do i = 1_int64, target_size
      target(i) = int(modulo(int(c_rand(), int64), n_classes), int32) + 1_int32
    end do

    size_average = .true.
    ignore_index = n_classes / 2_int64 + 1_int64

    call reference_nll(r_output, r_total_weight, input, target, weights, size_average, &
                       nframe, n_classes, ignore_index)

    call eval_nll(64, nframe, n_classes, size_average, ignore_index, r_output, &
                  r_total_weight, input, weights, target, repeat)
    call eval_nll(128, nframe, n_classes, size_average, ignore_index, r_output, &
                  r_total_weight, input, weights, target, repeat)
    call eval_nll(256, nframe, n_classes, size_average, ignore_index, r_output, &
                  r_total_weight, input, weights, target, repeat)
    call eval_nll(512, nframe, n_classes, size_average, ignore_index, r_output, &
                  r_total_weight, input, weights, target, repeat)
    call eval_nll(1024, nframe, n_classes, size_average, ignore_index, r_output, &
                  r_total_weight, input, weights, target, repeat)

    deallocate(input, weights, target)
  end subroutine driver

  real(real32) function rand_real() result(value)
    integer(int64), parameter :: rand_max = 2147483647_int64

    value = 2.0_real32 * (real(c_rand(), real32) / real(rand_max, real32)) - 1.0_real32
  end function rand_real

  subroutine reference_nll(output, total_weight, input, target, weights, size_average, &
                           nframe, kdim, ignore_index)
    real(real32), intent(out) :: output, total_weight
    real(real32), intent(in) :: input(:), weights(:)
    integer(int32), intent(in) :: target(:)
    logical, intent(in) :: size_average
    integer(int64), intent(in) :: nframe, kdim, ignore_index

    integer(int64) :: i, t, input_index
    real(real32) :: output_acc, total_weight_acc, cur_weight

    output_acc = 0.0_real32
    total_weight_acc = 0.0_real32
    do i = 1_int64, nframe
      t = int(target(i), int64)
      if (t /= ignore_index) then
        cur_weight = weights(t)
        input_index = (i - 1_int64) * kdim + t
        output_acc = output_acc - input(input_index) * cur_weight
        total_weight_acc = total_weight_acc + cur_weight
      end if
    end do

    total_weight = total_weight_acc
    if (size_average) then
      output = output_acc / total_weight_acc
    else
      output = output_acc
    end if
  end subroutine reference_nll

  subroutine eval_nll(gpu_threads, nframe, n_classes, size_average, ignore_index, &
                      r_output, r_total_weight, input, weights, target, repeat)
    integer, intent(in) :: gpu_threads, repeat
    integer(int64), intent(in) :: nframe, n_classes, ignore_index
    logical, intent(in) :: size_average
    real(real32), intent(in) :: r_output, r_total_weight
    real(real32), intent(in) :: input(:), weights(:)
    integer(int32), intent(in) :: target(:)

    real(real32) :: output(1), total_weight(1)
    real(real32), allocatable :: partial_output(:), partial_weight(:)
    integer(int64) :: input_size, weights_size, target_size
    integer :: iter
    real(real64) :: start_time, end_time, average_us
    logical :: ok

    input_size = nframe * n_classes
    weights_size = nframe
    target_size = nframe
    allocate(partial_output(gpu_threads), partial_weight(gpu_threads))

    start_time = omp_get_wtime()
    !$omp target data map(to: input(1:input_size), weights(1:weights_size), target(1:target_size)) &
    !$omp& map(alloc: partial_output(1:gpu_threads), partial_weight(1:gpu_threads)) &
    !$omp& map(from: output(1:1), total_weight(1:1))
    do iter = 1, repeat
      call nll_loss_kernel(output, total_weight, input, target, weights, partial_output, &
                           partial_weight, size_average, nframe, n_classes, ignore_index, &
                           gpu_threads)
    end do
    !$omp end target data
    end_time = omp_get_wtime()

    average_us = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)
    print *
    print '(A,I0)', 'Thread block size: ', gpu_threads
    print '(A,F0.6,A)', 'Average execution time of nll loss forward kernel: ', average_us, ' (us)'

    ok = abs(output(1) - r_output) <= 1.0e-1_real32 .and. &
         abs(total_weight(1) - r_total_weight) <= 1.0e-1_real32
    if (.not. ok) then
      print '(4(F0.6,1X))', output(1), r_output, total_weight(1), r_total_weight
    end if
    if (ok) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
    end if

    deallocate(partial_output, partial_weight)
  end subroutine eval_nll

  subroutine nll_loss_kernel(output, total_weight, input, target, weights, partial_output, &
                             partial_weight, size_average, nframe, kdim, ignore_index, gpu_threads)
    real(real32), intent(out) :: output(:), total_weight(:)
    real(real32), intent(in) :: input(:), weights(:)
    integer(int32), intent(in) :: target(:)
    real(real32), intent(inout) :: partial_output(:), partial_weight(:)
    logical, intent(in) :: size_average
    integer(int64), intent(in) :: nframe, kdim, ignore_index
    integer, intent(in) :: gpu_threads

    integer :: tid, nthreads, slot
    integer(int64) :: i, t, input_index
    real(real32) :: cur_weight, output_acc, weight_acc

    !$omp target teams num_teams(1) thread_limit(gpu_threads)
    !$omp parallel private(tid,nthreads,slot,i,t,input_index,cur_weight,output_acc,weight_acc)
      tid = omp_get_thread_num()
      nthreads = omp_get_num_threads()
      slot = tid + 1
      output_acc = 0.0_real32
      weight_acc = 0.0_real32

      do i = int(slot, int64), nframe, int(nthreads, int64)
        t = int(target(i), int64)
        if (t /= ignore_index) then
          cur_weight = weights(t)
          input_index = (i - 1_int64) * kdim + t
          output_acc = output_acc - input(input_index) * cur_weight
          weight_acc = weight_acc + cur_weight
        end if
      end do
      partial_output(slot) = output_acc
      partial_weight(slot) = weight_acc

      !$omp barrier

      if (tid == 0) then
        output_acc = 0.0_real32
        weight_acc = 0.0_real32
        do slot = 1, nthreads
          output_acc = output_acc + partial_output(slot)
          weight_acc = weight_acc + partial_weight(slot)
        end do
        total_weight(1) = weight_acc
        if (size_average) then
          output(1) = output_acc / weight_acc
        else
          output(1) = output_acc
        end if
      end if
    !$omp end parallel
    !$omp end target teams
  end subroutine nll_loss_kernel

end program nlll
