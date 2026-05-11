program main
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: num_threads = 256_int32

  character(len=256) :: arg0, arg1, arg2, arg3, arg4
  integer(int32) :: nrows, ndims, top_k, repeat, data_size
  integer(int32), allocatable :: labels(:)
  real(real32), allocatable :: data(:)
  integer(int32) :: count_ref, ngrid, iter, count_device
  real(real64) :: start_time, end_time, avg_us

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    print '(3A)', 'Usage: ', trim(arg0), ' <number of rows> <number of columns> <top K> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  call get_command_argument(4, arg4)
  read(arg1, *) nrows
  read(arg2, *) ndims
  read(arg3, *) top_k
  read(arg4, *) repeat
  if (nrows <= 0_int32 .or. ndims <= 0_int32 .or. top_k <= 0_int32 .or. repeat <= 0_int32) stop 1

  data_size = nrows * ndims
  allocate(labels(0:nrows - 1), data(0:data_size - 1))
  call initialize_inputs(nrows, ndims, labels, data)

  count_ref = reference_count(nrows, ndims, top_k, data, labels)

  !$omp target data map(to: labels(0:nrows - 1), data(0:data_size - 1))
  do ngrid = nrows / 4, nrows, nrows / 4
    print '(A,I0)', 'Grid size is ', ngrid
    start_time = omp_get_wtime()
    do iter = 1, repeat
      count_device = 0_int32
      call accuracy_kernel(nrows, ndims, top_k, ngrid, data, labels, count_device)
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0d6) / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of accuracy kernel: ', avg_us, ' (us)'
    if (count_device == count_ref) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
      stop 1
    end if
  end do
  !$omp end target data

  deallocate(labels, data)

contains

  subroutine initialize_inputs(rows, dims, labels, data)
    integer(int32), intent(in) :: rows, dims
    integer(int32), intent(out) :: labels(0:)
    real(real32), intent(out) :: data(0:)
    integer(int32) :: row, col, idx

    do row = 0, rows - 1
      labels(row) = modulo(row * 17_int32 + 3_int32, dims)
    end do

    do row = 0, rows - 1
      do col = 0, dims - 1
        idx = row * dims + col
        data(idx) = real(modulo(idx * 37_int32 + row * 13_int32 + col * 7_int32 + 123_int32, 4096_int32), real32) / 4096.0_real32
      end do
    end do
  end subroutine initialize_inputs

  integer(int32) function reference_count(rows, dims, topk, data, labels) result(count)
    integer(int32), intent(in) :: rows, dims, topk
    real(real32), intent(in) :: data(0:)
    integer(int32), intent(in) :: labels(0:)
    integer(int32) :: row, col, label_value, ngt
    real(real32) :: label_pred, pred

    count = 0_int32
    do row = 0, rows - 1
      label_value = labels(row)
      label_pred = data(row * dims + label_value)
      ngt = 0_int32
      do col = 0, dims - 1
        pred = data(row * dims + col)
        if (pred > label_pred .or. (pred == label_pred .and. col <= label_value)) ngt = ngt + 1_int32
      end do
      if (ngt <= topk) count = count + 1_int32
    end do
  end function reference_count

  subroutine accuracy_kernel(rows, dims, topk, teams, data, labels, count)
    integer(int32), intent(in) :: rows, dims, topk, teams
    real(real32), intent(in) :: data(0:)
    integer(int32), intent(in) :: labels(0:)
    integer(int32), intent(out) :: count
    integer(int32) :: row, col, label_value, ngt
    real(real32) :: label_pred, pred

    count = 0_int32
    !$omp target teams distribute parallel do num_teams(teams) num_threads(num_threads) reduction(+:count) &
    !$omp& private(col, label_value, label_pred, pred, ngt)
    do row = 0, rows - 1
      label_value = labels(row)
      label_pred = data(row * dims + label_value)
      ngt = 0_int32
      do col = 0, dims - 1
        pred = data(row * dims + col)
        if (pred > label_pred .or. (pred == label_pred .and. col <= label_value)) ngt = ngt + 1_int32
      end do
      if (ngt <= topk) count = count + 1_int32
    end do
    !$omp end target teams distribute parallel do
  end subroutine accuracy_kernel

end program main
