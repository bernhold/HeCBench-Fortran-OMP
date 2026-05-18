program main
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer, parameter :: chunk_size = 512
  integer, parameter :: value_count = 16

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
  integer :: nkeys, repeat_count
  integer(int32), allocatable :: keys(:), out(:)
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <number of keys> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) nkeys
  call get_command_argument(2, arg)
  read(arg, *) repeat_count
  if (nkeys <= 0 .or. repeat_count <= 0 .or. modulo(nkeys, chunk_size) /= 0) stop 1

  allocate(keys(nkeys), out(nkeys))
  call c_srand(512_c_int)
  call initialize(keys)
  out = keys

  call split_sort(out, repeat_count)
  ok = verify(out, keys)
  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
    stop 1
  end if

  deallocate(keys, out)

contains

  subroutine initialize(keys)
    integer(int32), intent(out) :: keys(:)
    integer :: i

    do i = 1, size(keys)
      ! Preserve the C++ original's srand(512)/rand() input stream.
      keys(i) = int(modulo(c_rand(), value_count), int32)
    end do
  end subroutine initialize

  subroutine split_sort(out, repeat_count)
    integer(int32), intent(inout) :: out(:)
    integer, intent(in) :: repeat_count
    integer :: iter, nkeys
    real(real64) :: start_time, elapsed_us

    nkeys = size(out)
    !$omp target data map(tofrom: out(1:nkeys))
    start_time = omp_get_wtime()
    do iter = 1, repeat_count
      call sort_chunks(out, nkeys)
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat_count, real64)
    write(*,'(A,F0.6,A)') 'Average kernel execution time: ', elapsed_us, ' (us)'
    !$omp end target data
  end subroutine split_sort

  subroutine sort_chunks(out, nkeys)
    integer(int32), intent(inout) :: out(:)
    integer, intent(in) :: nkeys
    integer :: chunk, base, i, value, pos
    integer :: hist0, hist1, hist2, hist3, hist4, hist5, hist6, hist7
    integer :: hist8, hist9, hist10, hist11, hist12, hist13, hist14, hist15

    !$omp target teams distribute parallel do thread_limit(1) private(base, i, value, pos, hist0, hist1, hist2, hist3, &
    !$omp& hist4, hist5, hist6, hist7, hist8, hist9, hist10, hist11, hist12, hist13, hist14, hist15)
    do chunk = 0, nkeys / chunk_size - 1
      base = chunk * chunk_size
      hist0 = 0; hist1 = 0; hist2 = 0; hist3 = 0
      hist4 = 0; hist5 = 0; hist6 = 0; hist7 = 0
      hist8 = 0; hist9 = 0; hist10 = 0; hist11 = 0
      hist12 = 0; hist13 = 0; hist14 = 0; hist15 = 0
      do i = 1, chunk_size
        value = out(base + i)
        select case (value)
        case (0); hist0 = hist0 + 1
        case (1); hist1 = hist1 + 1
        case (2); hist2 = hist2 + 1
        case (3); hist3 = hist3 + 1
        case (4); hist4 = hist4 + 1
        case (5); hist5 = hist5 + 1
        case (6); hist6 = hist6 + 1
        case (7); hist7 = hist7 + 1
        case (8); hist8 = hist8 + 1
        case (9); hist9 = hist9 + 1
        case (10); hist10 = hist10 + 1
        case (11); hist11 = hist11 + 1
        case (12); hist12 = hist12 + 1
        case (13); hist13 = hist13 + 1
        case (14); hist14 = hist14 + 1
        case (15); hist15 = hist15 + 1
        end select
      end do
      pos = base + 1
      call fill_value(out, pos, 0, hist0)
      call fill_value(out, pos, 1, hist1)
      call fill_value(out, pos, 2, hist2)
      call fill_value(out, pos, 3, hist3)
      call fill_value(out, pos, 4, hist4)
      call fill_value(out, pos, 5, hist5)
      call fill_value(out, pos, 6, hist6)
      call fill_value(out, pos, 7, hist7)
      call fill_value(out, pos, 8, hist8)
      call fill_value(out, pos, 9, hist9)
      call fill_value(out, pos, 10, hist10)
      call fill_value(out, pos, 11, hist11)
      call fill_value(out, pos, 12, hist12)
      call fill_value(out, pos, 13, hist13)
      call fill_value(out, pos, 14, hist14)
      call fill_value(out, pos, 15, hist15)
    end do
    !$omp end target teams distribute parallel do
  end subroutine sort_chunks

  subroutine fill_value(out, pos, value, count_value)
    integer(int32), intent(inout) :: out(:)
    integer, intent(inout) :: pos
    integer, intent(in) :: value, count_value
    integer :: i

    do i = 1, count_value
      out(pos) = int(value, int32)
      pos = pos + 1
    end do
  end subroutine fill_value

  logical function verify(sorted_keys, keys)
    integer(int32), intent(in) :: sorted_keys(:), keys(:)
    integer :: nkeys, base, i, value
    integer :: hist_original(0:value_count - 1), hist_sorted(0:value_count - 1)

    nkeys = size(sorted_keys)
    verify = .true.
    do base = 1, nkeys, chunk_size
      do i = 0, chunk_size - 2
        if (sorted_keys(base + i) > sorted_keys(base + i + 1)) then
          verify = .false.
          return
        end if
      end do
    end do
    do i = 1, nkeys
      if (sorted_keys(i) < 0 .or. sorted_keys(i) >= value_count) then
        verify = .false.
        return
      end if
    end do
    do base = 1, nkeys, chunk_size
      hist_original = 0
      hist_sorted = 0
      do i = 0, chunk_size - 1
        value = keys(base + i)
        hist_original(value) = hist_original(value) + 1
        value = sorted_keys(base + i)
        hist_sorted(value) = hist_sorted(value) + 1
      end do
      do value = 0, value_count - 1
        if (hist_original(value) /= hist_sorted(value)) then
          verify = .false.
          return
        end if
      end do
    end do
  end function verify

end program main
