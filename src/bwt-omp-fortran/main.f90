module bwt_mod
  use iso_fortran_env, only: int32, int64, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  integer, parameter :: etx = 0

contains

  pure integer(int32) function next_rand(state) result(value)
    integer(int32), intent(inout) :: state
    integer(int64) :: tmp

    tmp = mod(1103515245_int64 * int(state, int64) + 12345_int64, 2147483648_int64)
    state = int(tmp, int32)
    value = iand(ishft(state, -16), int(z'7fff', int32))
  end function next_rand

  subroutine generate_sequence(sequence, n)
    integer(int32), intent(out) :: sequence(0:)
    integer, intent(in) :: n
    integer(int32) :: rng
    integer :: i, pick
    integer(int32), parameter :: alphabet(0:3) = [iachar('A'), iachar('T'), iachar('C'), iachar('G')]

    rng = 123_int32
    do i = 0, n - 1
      pick = mod(next_rand(rng), 4)
      sequence(i) = alphabet(pick)
    end do
    sequence(n) = etx
  end subroutine generate_sequence

  subroutine generate_table(table, table_size, n)
    integer(int32), intent(out) :: table(0:)
    integer, intent(in) :: table_size, n
    integer :: i

    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 0, table_size - 1
      if (i < n) then
        table(i) = i
      else
        table(i) = -1
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine generate_table

  subroutine bitonic_sort_step(table, table_size, j, k, genome, n)
    integer(int32), intent(inout) :: table(0:)
    integer(int32), intent(in) :: genome(0:)
    integer, intent(in) :: table_size, j, k, n
    integer :: i, ixj, t1, t2, cmp_a, cmp_b, offset
    logical :: forward, swap_entries

    !$omp target teams distribute parallel do thread_limit(block_size) private(ixj, forward, t1, t2, cmp_a, cmp_b, offset, swap_entries)
    do i = 0, table_size - 1
      ixj = ieor(i, j)
      if (i < ixj) then
        forward = iand(i, k) == 0
        t1 = table(i)
        t2 = table(ixj)
        cmp_a = merge(t2, t1, forward)
        cmp_b = merge(t1, t2, forward)
        if (cmp_a < 0) then
          swap_entries = .false.
        else if (cmp_b < 0) then
          swap_entries = .true.
        else
          swap_entries = .false.
          do offset = 0, n - 1
            if (genome(mod(cmp_a + offset, n)) /= genome(mod(cmp_b + offset, n))) then
              swap_entries = genome(mod(cmp_a + offset, n)) < genome(mod(cmp_b + offset, n))
              exit
            end if
          end do
        end if
        if (swap_entries) then
          table(i) = t2
          table(ixj) = t1
        end if
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine bitonic_sort_step

  logical function compare_rotations_cpu(a, b, genome, n) result(is_less)
    integer, intent(in) :: a, b, n
    integer(int32), intent(in) :: genome(0:)
    integer :: i

    if (a < 0) then
      is_less = .false.
      return
    end if
    if (b < 0) then
      is_less = .true.
      return
    end if

    do i = 0, n - 1
      if (genome(mod(a + i, n)) /= genome(mod(b + i, n))) then
        is_less = genome(mod(a + i, n)) < genome(mod(b + i, n))
        return
      end if
    end do
    is_less = .false.
  end function compare_rotations_cpu

  subroutine reconstruct_sequence(table, sequence, transformed, n)
    integer(int32), intent(in) :: table(0:)
    integer(int32), intent(in) :: sequence(0:)
    integer(int32), intent(out) :: transformed(0:)
    integer, intent(in) :: n
    integer :: i

    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 0, n - 1
      transformed(i) = sequence(mod(n + table(i) - 1, n))
    end do
    !$omp end target teams distribute parallel do
  end subroutine reconstruct_sequence

  subroutine bwt_gpu(sequence, transformed, suffix_table, n)
    integer(int32), intent(in) :: sequence(0:)
    integer(int32), intent(out) :: transformed(0:)
    integer(int32), allocatable, intent(out) :: suffix_table(:)
    integer, intent(in) :: n
    integer :: table_size, j, k

    table_size = 1
    do while (table_size < n)
      table_size = table_size * 2
    end do

    allocate(suffix_table(0:table_size - 1))

    !$omp target data map(from: suffix_table(0:table_size - 1), transformed(0:n - 1)) map(to: sequence(0:n - 1))
    call generate_table(suffix_table, table_size, n)

    k = 2
    do while (k <= table_size)
      j = k / 2
      do while (j > 0)
        call bitonic_sort_step(suffix_table, table_size, j, k, sequence, n)
        j = j / 2
      end do
      k = k * 2
    end do

    call reconstruct_sequence(suffix_table, sequence, transformed, n)
    !$omp end target data
  end subroutine bwt_gpu

  subroutine bwt_cpu(sequence, transformed, n)
    integer(int32), intent(in) :: sequence(0:)
    integer(int32), intent(out) :: transformed(0:)
    integer, intent(in) :: n
    integer(int32), allocatable :: table(:)
    integer :: i, j, tmp

    allocate(table(0:n - 1))
    do i = 0, n - 1
      table(i) = i
    end do

    do i = 1, n - 1
      tmp = table(i)
      j = i - 1
      do while (j >= 0 .and. compare_rotations_cpu(tmp, table(j), sequence, n))
        table(j + 1) = table(j)
        j = j - 1
      end do
      table(j + 1) = tmp
    end do

    do i = 0, n - 1
      transformed(i) = sequence(mod(n + table(i) - 1, n))
    end do

    deallocate(table)
  end subroutine bwt_cpu

  logical function arrays_equal(a, b, n) result(equal)
    integer(int32), intent(in) :: a(0:), b(0:)
    integer, intent(in) :: n
    integer :: i

    equal = .true.
    do i = 0, n - 1
      if (a(i) /= b(i)) then
        equal = .false.
        return
      end if
    end do
  end function arrays_equal

end module bwt_mod

program main
  use iso_fortran_env, only: int32, real64
  use omp_lib
  use bwt_mod
  implicit none

  integer :: n, arg_count
  character(len=64) :: arg
  integer(int32), allocatable :: sequence(:), cpu_seq(:), gpu_seq(:), suffix_table(:)
  real(real64) :: start_time, cpu_ms, gpu_ms

  arg_count = command_argument_count()
  if (arg_count > 0) then
    call get_command_argument(1, arg)
    read(arg, *) n
  else
    n = 1000000
  end if

  write(*,'(A,I0)') 'running a sample sequence of length ', n

  allocate(sequence(0:n), cpu_seq(0:n), gpu_seq(0:n))
  call generate_sequence(sequence, n)

  start_time = omp_get_wtime()
  call bwt_cpu(sequence, cpu_seq, n)
  cpu_ms = (omp_get_wtime() - start_time) * 1000.0_real64

  start_time = omp_get_wtime()
  call bwt_gpu(sequence, gpu_seq, suffix_table, n)
  gpu_ms = (omp_get_wtime() - start_time) * 1000.0_real64

  write(*,'(A,I0,A)') 'Host time: ', nint(cpu_ms), ' ms'
  write(*,'(A,I0,A)') 'Device time: ', nint(gpu_ms), ' ms'

  if (arrays_equal(cpu_seq, gpu_seq, n)) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(sequence, cpu_seq, gpu_seq, suffix_table)
end program main
