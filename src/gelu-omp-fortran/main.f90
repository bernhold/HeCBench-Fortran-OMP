program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name='rand') result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  character(len=256) :: arg0, arg
  integer :: batch_size, seq_len, hidden_dim, repeat, block_size, i
  integer(int64) :: src_size
  real(real32), allocatable :: output(:), output_ref(:), bias(:)
  real(real64) :: start_time, end_time
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <batch> <sequence length> <hidden dimension> <repeat>'
    write(*,'(A)') 'The hidden dimension is a multiple of two'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) batch_size
  call get_command_argument(2, arg)
  read(arg, *) seq_len
  call get_command_argument(3, arg)
  read(arg, *) hidden_dim
  call get_command_argument(4, arg)
  read(arg, *) repeat
  if (batch_size <= 0 .or. seq_len <= 0 .or. hidden_dim <= 0 .or. repeat <= 0) stop 1
  if (mod(hidden_dim, 2) /= 0) stop 1

  src_size = int(batch_size, int64) * int(seq_len, int64) * int(hidden_dim, int64)
  allocate(output(0:src_size - 1), output_ref(0:src_size - 1), bias(0:hidden_dim - 1))
  call initialize_inputs(output, output_ref, bias, src_size, hidden_dim)

  if (hidden_dim >= 4096) then
    block_size = 512
  else if (hidden_dim >= 2048) then
    block_size = 256
  else
    block_size = 128
  end if

  call gelu_bias_loop_cpu(output_ref, bias, batch_size, hidden_dim, seq_len)

  !$omp target data map(to: bias(0:hidden_dim - 1)) map(tofrom: output(0:src_size - 1))
  call gelu_bias_loop(output, bias, batch_size, hidden_dim, seq_len, block_size)
  !$omp target update from(output(0:src_size - 1))

  ok = .true.
  do i = 0, int(src_size - 1_int64)
    if (abs(output_ref(i) - output(i)) > 1.0e-3_real32) then
      ok = .false.
      exit
    end if
  end do
  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  start_time = omp_get_wtime()
  do i = 1, repeat
    call gelu_bias_loop(output, bias, batch_size, hidden_dim, seq_len, block_size)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average execution time of vectorized kernel ', &
    (end_time - start_time) * 1.0e3_real64 / real(repeat, real64), ' (ms)'

  start_time = omp_get_wtime()
  do i = 1, repeat
    call gelu_bias_loop_base(output, bias, batch_size, hidden_dim, seq_len, block_size)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average execution time of baseline kernel ', &
    (end_time - start_time) * 1.0e3_real64 / real(repeat, real64), ' (ms)'
  !$omp end target data

  deallocate(output, output_ref, bias)

contains

  subroutine initialize_inputs(output, output_ref, bias, src_size, hidden_dim)
    real(real32), intent(out) :: output(0:), output_ref(0:), bias(0:)
    integer(int64), intent(in) :: src_size
    integer, intent(in) :: hidden_dim
    integer(int64) :: i
    integer :: j

    call c_srand(123_c_int)
    do i = 0_int64, src_size - 1_int64
      output(i) = real(c_rand(), real32) / real(huge(0_c_int), real32)
      output_ref(i) = output(i)
    end do
    do j = 0, hidden_dim - 1
      bias(j) = real(-6_c_int + mod(c_rand(), 12_c_int), real32)
    end do
  end subroutine initialize_inputs

  subroutine gelu_bias_loop(src, bias, batch_size, width, height, block_size)
    real(real32), intent(inout) :: src(0:)
    real(real32), intent(in) :: bias(0:)
    integer, intent(in) :: batch_size, width, height, block_size
    integer :: batch, x, y
    integer(int64) :: base, left, right
    real(real32) :: tx, ty

    !$omp target teams distribute collapse(2) num_teams(batch_size * height) private(batch, x, y, base, left, right, tx, ty)
    do batch = 0, batch_size - 1
      do x = 0, height - 1
        base = int(batch, int64) * width * height + int(x, int64) * width
        !$omp parallel do num_threads(block_size) private(y, left, right, tx, ty)
        do y = 0, width / 2 - 1
          left = base + int(2 * y, int64)
          right = left + 1_int64
          tx = src(left) + bias(2 * y)
          ty = src(right) + bias(2 * y + 1)
          src(left) = gelu_value(tx)
          src(right) = gelu_value(ty)
        end do
        !$omp end parallel do
      end do
    end do
    !$omp end target teams distribute
  end subroutine gelu_bias_loop

  subroutine gelu_bias_loop_base(src, bias, batch_size, width, height, block_size)
    real(real32), intent(inout) :: src(0:)
    real(real32), intent(in) :: bias(0:)
    integer, intent(in) :: batch_size, width, height, block_size
    integer :: batch, x, y
    integer(int64) :: base, idx
    real(real32) :: t

    !$omp target teams distribute collapse(2) num_teams(batch_size * height) private(batch, x, y, base, idx, t)
    do batch = 0, batch_size - 1
      do x = 0, height - 1
        base = int(batch, int64) * width * height + int(x, int64) * width
        !$omp parallel do num_threads(block_size) private(y, idx, t)
        do y = 0, width - 1
          idx = base + int(y, int64)
          t = src(idx) + bias(y)
          src(idx) = gelu_value(t)
        end do
        !$omp end parallel do
      end do
    end do
    !$omp end target teams distribute
  end subroutine gelu_bias_loop_base

  subroutine gelu_bias_loop_cpu(src, bias, batch_size, width, height)
    real(real32), intent(inout) :: src(0:)
    real(real32), intent(in) :: bias(0:)
    integer, intent(in) :: batch_size, width, height
    integer :: batch, x, y
    integer(int64) :: idx
    real(real32) :: t

    do batch = 0, batch_size - 1
      do x = 0, height - 1
        do y = 0, width - 1
          idx = int(batch, int64) * width * height + int(x, int64) * width + int(y, int64)
          t = src(idx) + bias(y)
          src(idx) = gelu_value(t)
        end do
      end do
    end do
  end subroutine gelu_bias_loop_cpu

  pure real(real32) function gelu_value(t)
    real(real32), intent(in) :: t

    gelu_value = 0.5_real32 * t * (1.0_real32 + tanh(0.79788456_real32 * (t + 0.044715_real32 * t * t * t)))
  end function gelu_value

end program main
