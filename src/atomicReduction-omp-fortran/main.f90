program atomic_reduction_main
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: num_block_sizes = 4
  integer(int32), parameter :: block_sizes(num_block_sizes) = [128_int32, 256_int32, 512_int32, 1024_int32]
  integer(int32) :: array_length, repeats, checksum
  integer(int32), allocatable :: array(:)
  integer :: i, k, threads, blocks
  integer(int32) :: warmup_sum
  real(real32) :: gb
  character(len=64) :: arg

  array_length = 52428800_int32
  repeats = 100_int32
  if (command_argument_count() == 2) then
    call get_command_argument(1, arg)
    read(arg, *) array_length
    call get_command_argument(2, arg)
    read(arg, *) repeats
  end if

  write(*,'("Array size: ",G0," MB")') real(array_length * 4_int32, real64) / 1024.0_real64 / 1024.0_real64
  write(*,'("Repeat the kernel execution: ",I0," times")') repeats

  allocate(array(0:array_length-1))
  checksum = 0_int32
  do i = 0, array_length - 1
    array(i) = int(mod(i, 2), int32)
    checksum = checksum + array(i)
  end do

  gb = real(array_length, real32) * 4.0_real32 * real(repeats, real32)

  !$omp target data map(to: array)
  do i = 1, repeats
    call reduction_variant(array, array_length, 1_int32, 2048, 256, warmup_sum)
  end do

  do k = 1, num_block_sizes
    threads = block_sizes(k)
    blocks = min((array_length + threads - 1) / threads, 2048)
    call benchmark_variant(array, array_length, repeats, threads, blocks, 1_int32, checksum, gb, .true.)
    call benchmark_variant(array, array_length, repeats, threads, blocks / 2, 2_int32, checksum, gb, .false.)
    call benchmark_variant(array, array_length, repeats, threads, blocks / 4, 4_int32, checksum, gb, .false.)
    call benchmark_variant(array, array_length, repeats, threads, blocks / 8, 8_int32, checksum, gb, .false.)
    call benchmark_variant(array, array_length, repeats, threads, blocks / 16, 16_int32, checksum, gb, .false.)
  end do
  !$omp end target data

contains

  subroutine benchmark_variant(array, array_length, repeats, threads, blocks, width, checksum, gb, print_sum)
    integer(int32), intent(in) :: array(0:), array_length, repeats, threads, blocks, width, checksum
    real(real32), intent(in) :: gb
    logical, intent(in) :: print_sum
    integer(int32) :: sum
    integer :: n
    real(real64) :: start_time, elapsed

    start_time = omp_get_wtime()
    do n = 1, repeats
      call reduction_variant(array, array_length, width, blocks, threads, sum)
    end do
    elapsed = omp_get_wtime() - start_time
    write(*,'("Thread block size: ",I0,", The average performance of reduction is ",G0," GBytes/sec")') &
        threads, 1.0e-09_real64 * real(gb, real64) / elapsed
    if (print_sum) then
      write(*,'(I0,1X,I0)') sum, checksum
    end if
    if (sum == checksum) then
      write(*,'("VERIFICATION: PASS")')
      write(*,*)
    else
      write(*,'("VERIFICATION: FAIL!!")')
      write(*,*)
    end if
  end subroutine benchmark_variant

  subroutine reduction_variant(array, array_length, width, blocks, threads, sum)
    integer(int32), intent(in) :: array(0:), array_length, width
    integer, intent(in) :: blocks, threads
    integer(int32), intent(out) :: sum
    integer(int32) :: local_sum
    integer :: i

    local_sum = 0_int32
    select case (width)
    case (1_int32)
      !$omp target teams distribute parallel do num_teams(blocks) thread_limit(threads) reduction(+:local_sum)
      do i = 0, array_length - 1
        local_sum = local_sum + array(i)
      end do
      !$omp end target teams distribute parallel do
    case (2_int32)
      !$omp target teams distribute parallel do num_teams(blocks) thread_limit(threads) reduction(+:local_sum)
      do i = 0, array_length / 2 - 1
        local_sum = local_sum + array(i * 2) + array(i * 2 + 1)
      end do
      !$omp end target teams distribute parallel do
    case (4_int32)
      !$omp target teams distribute parallel do num_teams(blocks) thread_limit(threads) reduction(+:local_sum)
      do i = 0, array_length / 4 - 1
        local_sum = local_sum + array(i * 4) + array(i * 4 + 1) + array(i * 4 + 2) + array(i * 4 + 3)
      end do
      !$omp end target teams distribute parallel do
    case (8_int32)
      !$omp target teams distribute parallel do num_teams(blocks) thread_limit(threads) reduction(+:local_sum)
      do i = 0, array_length / 8 - 1
        local_sum = local_sum + array(i * 8) + array(i * 8 + 1) + array(i * 8 + 2) + array(i * 8 + 3) + &
            array(i * 8 + 4) + array(i * 8 + 5) + array(i * 8 + 6) + array(i * 8 + 7)
      end do
      !$omp end target teams distribute parallel do
    case (16_int32)
      !$omp target teams distribute parallel do num_teams(blocks) thread_limit(threads) reduction(+:local_sum)
      do i = 0, array_length / 16 - 1
        local_sum = local_sum + array(i * 16) + array(i * 16 + 1) + array(i * 16 + 2) + array(i * 16 + 3) + &
            array(i * 16 + 4) + array(i * 16 + 5) + array(i * 16 + 6) + array(i * 16 + 7) + &
            array(i * 16 + 8) + array(i * 16 + 9) + array(i * 16 + 10) + array(i * 16 + 11) + &
            array(i * 16 + 12) + array(i * 16 + 13) + array(i * 16 + 14) + array(i * 16 + 15)
      end do
      !$omp end target teams distribute parallel do
    end select

    sum = local_sum
  end subroutine reduction_variant

end program atomic_reduction_main
