program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
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

  type params_t
    integer :: n_gpu_threads = 64
    integer :: n_gpu_blocks = 16
    integer :: n_warmup = 10
    integer :: n_reps = 100
    integer :: m = 197
    integer :: n = 35588
    integer :: s = 32
  end type params_t

  type(params_t) :: p
  integer :: tiled_n, padded_n, in_size, status
  real(real32), allocatable :: input_backup(:), device_output(:), host_output(:)
  real(real64) :: elapsed

  call parse_args(p)
  if (p%n_gpu_threads > 256) error stop 'The thread block size is greater than the maximum thread block size that can be used on this device'
  if (p%n_gpu_threads <= 0 .or. p%n_gpu_blocks <= 0 .or. p%n_warmup < 0 .or. p%n_reps <= 0) error stop 'invalid launch or repeat count'
  if (p%m <= 0 .or. p%n <= 0 .or. p%s <= 0) error stop 'invalid matrix dimensions'

  tiled_n = divceil(p%n, p%s)
  padded_n = tiled_n * p%s
  in_size = p%m * padded_n

  allocate(input_backup(in_size), device_output(in_size), host_output(in_size))
  call read_input(input_backup)
  device_output = 0.0_real32
  host_output = 0.0_real32

  call asta_device(input_backup, device_output, p, tiled_n, padded_n, in_size, elapsed)
  write(*, '(A,F0.6,A)') 'Average kernel execution time ', elapsed / real(p%n_reps, real64), ' (s)'

  call cpu_soa_asta(input_backup, host_output, padded_n, p%m, p%s)
  status = compare_output(device_output, host_output)
  if (status == 0) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(input_backup, device_output, host_output)

contains

  integer function divceil(n, m)
    integer, intent(in) :: n, m
    divceil = (n - 1) / m + 1
  end function divceil

  subroutine parse_args(p)
    type(params_t), intent(inout) :: p
    integer :: argc, pos
    character(len=128) :: arg

    argc = command_argument_count()
    pos = 1
    do while (pos <= argc)
      call get_command_argument(pos, arg)
      select case (trim(arg))
      case ('-h')
        call usage()
        stop
      case ('-i')
        p%n_gpu_threads = read_int_arg(pos, argc)
      case ('-g')
        p%n_gpu_blocks = read_int_arg(pos, argc)
      case ('-w')
        p%n_warmup = read_int_arg(pos, argc)
      case ('-r')
        p%n_reps = read_int_arg(pos, argc)
      case ('-m')
        p%m = read_int_arg(pos, argc)
      case ('-n')
        p%n = read_int_arg(pos, argc)
      case ('-s')
        p%s = read_int_arg(pos, argc)
      case default
        write(*, '(A)') ''
        write(*, '(A)') 'Unrecognized option!'
        call usage()
        stop 1
      end select
      pos = pos + 2
    end do
  end subroutine parse_args

  integer function read_int_arg(pos, argc)
    integer, intent(in) :: pos, argc
    character(len=128) :: value
    if (pos + 1 > argc) then
      write(*, '(A)') ''
      write(*, '(A)') 'Unrecognized option!'
      call usage()
      stop 1
    end if
    call get_command_argument(pos + 1, value)
    read(value, *) read_int_arg
  end function read_int_arg

  subroutine usage()
    write(*, '(A)') ''
    write(*, '(A)') 'Usage:  ./main [options]'
    write(*, '(A)') ''
    write(*, '(A)') 'General options:'
    write(*, '(A)') '    -h        help'
    write(*, '(A)') '    -i <I>    # of device threads per block (default=64)'
    write(*, '(A)') '    -g <G>    # of device blocks (default=16)'
    write(*, '(A)') '    -w <W>    # of warmup iterations (default=10)'
    write(*, '(A)') '    -r <R>    # of repetition iterations (default=100)'
    write(*, '(A)') ''
    write(*, '(A)') 'Benchmark-specific options:'
    write(*, '(A)') '    -m <M>    matrix height (default=197)'
    write(*, '(A)') '    -n <N>    matrix width (default=35588)'
    write(*, '(A)') '    -s <M>    super-element size (default=32)'
    write(*, '(A)') ''
  end subroutine usage

  subroutine read_input(x)
    real(real32), intent(out) :: x(:)
    integer :: i
    integer(c_int) :: value

    call c_srand(5432_c_int)
    do i = 1, size(x)
      value = c_rand()
      x(i) = real(mod(value, 100_c_int), real32) / 100.0_real32
    end do
  end subroutine read_input

  subroutine asta_device(input, output, p, tiled_n, padded_n, in_size, elapsed)
    real(real32), intent(in) :: input(:)
    real(real32), intent(inout) :: output(:)
    type(params_t), intent(in) :: p
    integer, intent(in) :: tiled_n, padded_n, in_size
    real(real64), intent(out) :: elapsed
    integer :: rep, tile, row, elem, src_idx, dst_idx
    real(real64) :: start_time, end_time

    start_time = 0.0_real64
    !$omp target data map(to: input(1:in_size)) map(from: output(1:in_size))
    do rep = 1, p%n_warmup + p%n_reps
      if (rep == p%n_warmup + 1) start_time = omp_get_wtime()
      !$omp target teams distribute parallel do collapse(3) num_teams(p%n_gpu_blocks) thread_limit(p%n_gpu_threads) &
      !$omp& private(src_idx, dst_idx)
      do tile = 1, tiled_n
        do row = 1, p%m
          do elem = 1, p%s
            src_idx = (row - 1) * padded_n + (tile - 1) * p%s + elem
            dst_idx = (tile - 1) * p%m * p%s + (row - 1) * p%s + elem
            output(dst_idx) = input(src_idx)
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    elapsed = end_time - start_time
  end subroutine asta_device

  subroutine cpu_soa_asta(src, dst, height, width, tile_size)
    real(real32), intent(in) :: src(:)
    real(real32), intent(out) :: dst(:)
    integer, intent(in) :: height, width, tile_size
    integer :: k, i, j, src_idx, dst_idx

    dst = 0.0_real32
    if ((height / tile_size) * tile_size == height) then
      do k = 1, width
        do i = 1, height / tile_size
          do j = 1, tile_size
            src_idx = (k - 1) * height + (i - 1) * tile_size + j
            dst_idx = (i - 1) * width * tile_size + (k - 1) * tile_size + j
            dst(dst_idx) = src(src_idx)
          end do
        end do
      end do
    end if
  end subroutine cpu_soa_asta

  integer function compare_output(output, ref)
    real(real32), intent(in) :: output(:), ref(:)
    integer :: i
    real(real32) :: diff

    compare_output = 0
    do i = 1, size(output)
      diff = abs(ref(i) - output(i))
      if ((diff - 0.0_real32) > 0.00001_real32 .and. diff > 0.01_real32 * abs(ref(i))) then
        write(*, '(A,I0,A,F0.6,A,F0.6,A,F0.6)') 'Failed at line: ', i - 1, ' ref: ', ref(i), ' actual: ', output(i), ' diff: ', diff
        compare_output = 1
        exit
      end if
    end do
  end function compare_output

end program main
