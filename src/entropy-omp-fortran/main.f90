program entropy_main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  integer :: width, height, repeat, i
  integer :: input_size
  integer, allocatable :: input(:)
  real(real32), allocatable :: output(:), output_ref(:)
  real(real32) :: log_table(0:25)
  real(real64) :: start_time, elapsed
  integer(int64) :: seed
  character(len=64) :: arg
  logical :: ok

  if (command_argument_count() /= 3) then
    write(*,'("Usage: ./main <width> <height> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) width
  call get_command_argument(2, arg)
  read(arg, *) height
  call get_command_argument(3, arg)
  read(arg, *) repeat

  input_size = width * height
  allocate(input(0:input_size-1), output(0:input_size-1), output_ref(0:input_size-1))

  do i = 0, 25
    if (i <= 1) then
      log_table(i) = 0.0_real32
    else
      log_table(i) = real(i, real32) * log2_real(real(i, real32))
    end if
  end do

  seed = 123_int64
  do i = 0, input_size - 1
    input(i) = next_mod_16(seed)
  end do

  !$omp target data map(to: input, log_table) map(from: output)
  start_time = omp_get_wtime()
  do i = 1, repeat
    call entropy_baseline(output, input, height, width)
  end do
  elapsed = (omp_get_wtime() - start_time) / real(repeat, real64)
  write(*,'("Average kernel (baseline) execution time ",F0.6," (s)")') elapsed

  start_time = omp_get_wtime()
  do i = 1, repeat
    call entropy_optimized(output, input, log_table, height, width)
  end do
  elapsed = (omp_get_wtime() - start_time) / real(repeat, real64)
  !$omp end target data
  write(*,'("Average kernel (optimized) execution time ",F0.6," (s)")') elapsed

  call reference_entropy(output_ref, input, height, width)
  ok = check_entropy(output, output_ref, input_size)
  if (ok) then
    write(*,'("PASS")')
  else
    write(*,'("FAIL")')
    stop 1
  end if

contains

  subroutine entropy_baseline(entropy_out, values, rows, cols)
    real(real32), intent(out) :: entropy_out(0:)
    integer, intent(in) :: values(0:)
    integer, intent(in) :: rows, cols
    integer :: y, x, dy, dx, xx, yy, k, total
    integer :: counts(0:15)
    real(real32) :: p, s

    !$omp target teams distribute parallel do collapse(2) private(dy, dx, xx, yy, k, total, counts, p, s) thread_limit(256)
    do y = 0, rows - 1
      do x = 0, cols - 1
        counts = 0
        total = 0
        do dy = -2, 2
          do dx = -2, 2
            xx = x + dx
            yy = y + dy
            if (xx >= 0 .and. yy >= 0 .and. yy < rows .and. xx < cols) then
              counts(values(yy * cols + xx)) = counts(values(yy * cols + xx)) + 1
              total = total + 1
            end if
          end do
        end do
        if (total < 1) total = 1
        s = 0.0_real32
        do k = 0, 15
          if (counts(k) > 0) then
            p = real(counts(k), real32) / real(total, real32)
            s = s - p * (log(p) / log(2.0_real32))
          end if
        end do
        entropy_out(y * cols + x) = s
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine entropy_baseline

  subroutine entropy_optimized(entropy_out, values, log_table, rows, cols)
    real(real32), intent(out) :: entropy_out(0:)
    integer, intent(in) :: values(0:)
    real(real32), intent(in) :: log_table(0:)
    integer, intent(in) :: rows, cols
    integer :: y, x, dy, dx, xx, yy, k, total
    integer :: counts(0:15)
    real(real32) :: s

    !$omp target teams distribute parallel do collapse(2) private(dy, dx, xx, yy, k, total, counts, s) thread_limit(256)
    do y = 0, rows - 1
      do x = 0, cols - 1
        counts = 0
        total = 0
        do dy = -2, 2
          do dx = -2, 2
            xx = x + dx
            yy = y + dy
            if (xx >= 0 .and. yy >= 0 .and. yy < rows .and. xx < cols) then
              counts(values(yy * cols + xx)) = counts(values(yy * cols + xx)) + 1
              total = total + 1
            end if
          end do
        end do
        if (total < 1) total = 1
        s = 0.0_real32
        do k = 0, 15
          s = s - log_table(counts(k))
        end do
        entropy_out(y * cols + x) = s / real(total, real32) + &
            log(real(total, real32)) / log(2.0_real32)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine entropy_optimized

  subroutine reference_entropy(entropy_out, values, rows, cols)
    real(real32), intent(out) :: entropy_out(0:)
    integer, intent(in) :: values(0:)
    integer, intent(in) :: rows, cols
    integer :: y, x, dy, dx, xx, yy, k, total
    integer :: counts(0:15)
    real(real32) :: p, s

    do y = 0, rows - 1
      do x = 0, cols - 1
        counts = 0
        total = 0
        do dx = -2, 2
          do dy = -2, 2
            xx = x + dx
            yy = y + dy
            if (xx >= 0 .and. yy >= 0 .and. yy < rows .and. xx < cols) then
              counts(values(yy * cols + xx)) = counts(values(yy * cols + xx)) + 1
              total = total + 1
            end if
          end do
        end do
        if (total < 1) total = 1
        s = 0.0_real32
        do k = 0, 15
          if (counts(k) > 0) then
            p = real(counts(k), real32) / real(total, real32)
            s = s - p * log2_real(p)
          end if
        end do
        entropy_out(y * cols + x) = s
      end do
    end do
  end subroutine reference_entropy

  logical function check_entropy(device_values, reference_values, size)
    real(real32), intent(in) :: device_values(0:), reference_values(0:)
    integer, intent(in) :: size
    integer :: idx

    check_entropy = .true.
    do idx = 0, size - 1
      if (abs(device_values(idx) - reference_values(idx)) > 1.0e-3_real32) then
        check_entropy = .false.
        exit
      end if
    end do
  end function check_entropy

  real(real32) function log2_real(value)
    real(real32), intent(in) :: value
    log2_real = log(value) / log(2.0_real32)
  end function log2_real

  integer function next_mod_16(seed)
    integer(int64), intent(inout) :: seed
    seed = mod(seed * 1103515245_int64 + 12345_int64, 2147483648_int64)
    next_mod_16 = int(mod(seed, 16_int64))
  end function next_mod_16

end program entropy_main
