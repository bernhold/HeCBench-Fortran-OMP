module crc64_mod
  use iso_fortran_env, only: int32, int64, real64
  implicit none

  integer(int64), parameter :: crc64_poly = int(z'c96c5795d7870f42', int64)
  integer(int64), parameter :: all_ones = int(z'ffffffffffffffff', int64)
  integer(int64), parameter :: high_bit = int(z'8000000000000000', int64)
  integer(int64), parameter :: x_term = int(z'4000000000000000', int64)
  integer(int64), parameter :: mask48 = int(z'0000ffffffffffff', int64)
  integer(int64), parameter :: lcg_a = 25214903917_int64
  integer(int64), parameter :: lcg_c = 11_int64
  integer(int64), parameter :: base24 = 16777216_int64
  real(real64), parameter :: two48 = 281474976710656.0_real64

contains

  subroutine srand48_state(seed, state)
    integer(int32), intent(in) :: seed
    integer(int64), intent(out) :: state

    state = iand(ior(shiftl(int(seed, int64), 16), int(z'330e', int64)), mask48)
  end subroutine srand48_state

  subroutine next_drand48(state, value)
    integer(int64), intent(inout) :: state
    real(real64), intent(out) :: value
    integer(int64) :: a0, a1, x0, x1, mid

    a0 = iand(lcg_a, base24 - 1_int64)
    a1 = shiftr(lcg_a, 24)
    x0 = iand(state, base24 - 1_int64)
    x1 = shiftr(state, 24)
    mid = iand(a0 * x1 + a1 * x0, base24 - 1_int64)
    state = iand(a0 * x0 + shiftl(mid, 24) + lcg_c, mask48)
    value = real(state, real64) / two48
  end subroutine next_drand48

  pure integer(int64) function crc64_serial(buffer, first, last) result(cs_out)
    integer(int32), intent(in) :: buffer(:)
    integer(int64), intent(in) :: first, last
    integer(int64) :: cs, idx, pos
    integer :: bit

    cs = all_ones
    do pos = first, last
      cs = ieor(cs, int(iand(buffer(pos), 255_int32), int64))
      do bit = 1, 8
        if (iand(cs, 1_int64) /= 0_int64) then
          cs = ieor(shiftr(cs, 1), crc64_poly)
        else
          cs = shiftr(cs, 1)
        end if
      end do
    end do
    cs_out = ieor(cs, all_ones)
  end function crc64_serial

  pure integer(int64) function crc64_multiply(a_in, b_in) result(r)
    integer(int64), intent(in) :: a_in, b_in
    integer(int64) :: a, b
    integer :: bit

    a = a_in
    b = b_in
    r = 0_int64
    do bit = 1, 64
      if (iand(a, high_bit) /= 0_int64) r = ieor(r, b)
      a = shiftl(a, 1)
      if (iand(b, 1_int64) /= 0_int64) then
        b = ieor(shiftr(b, 1), crc64_poly)
      else
        b = shiftr(b, 1)
      end if
    end do
  end function crc64_multiply

  subroutine crc64_pow_table(pow2)
    integer(int64), intent(out) :: pow2(64)
    integer :: i

    pow2(1) = x_term
    do i = 2, 64
      pow2(i) = crc64_multiply(pow2(i - 1), pow2(i - 1))
    end do
  end subroutine crc64_pow_table

  integer(int64) function crc64_x_pow_n(n) result(r)
    integer(int64), intent(in) :: n
    integer(int64) :: pow2(64), work
    integer :: i

    call crc64_pow_table(pow2)
    r = high_bit
    work = n
    i = 1
    do while (work /= 0_int64)
      if (iand(work, 1_int64) /= 0_int64) r = crc64_multiply(r, pow2(i))
      work = shiftr(work, 1)
      i = i + 1
    end do
  end function crc64_x_pow_n

  integer(int64) function crc64_combine(cs1, cs2, nbytes2) result(cs)
    integer(int64), intent(in) :: cs1, cs2, nbytes2

    cs = ieor(cs2, crc64_multiply(cs1, crc64_x_pow_n(8_int64 * nbytes2)))
  end function crc64_combine

  subroutine crc64_invert(cs, check_bytes)
    integer(int64), intent(in) :: cs
    integer(int32), intent(out) :: check_bytes(8)
    integer(int64) :: work
    integer :: i

    work = ieor(cs, all_ones)
    do i = 1, 8
      check_bytes(i) = int(iand(shiftr(work, 8 * (i - 1)), 255_int64), int32)
    end do
  end subroutine crc64_invert

  integer(int64) function crc64_omp(buffer, nbytes) result(cs)
    integer(int32), intent(in) :: buffer(:)
    integer(int64), intent(in) :: nbytes
    integer(int64), allocatable :: chunk_cs(:), chunk_sz(:)
    integer(int64) :: bpt, first, last, local_cs
    integer :: tid, nthreads, bit
    integer(int64) :: pos

    if (nbytes <= 2048_int64) then
      cs = crc64_serial(buffer, 1_int64, nbytes)
      return
    end if

    nthreads = 96 * 8 * 32
    if (nbytes < int(nthreads, int64) * 1024_int64) then
      nthreads = max(1, int(nbytes / 1024_int64))
    end if

    allocate(chunk_cs(nthreads), chunk_sz(nthreads))
    bpt = nbytes / int(nthreads, int64)

    !$omp target teams distribute parallel do num_teams(max(1, nthreads / 64)) thread_limit(64) &
    !$omp& map(to: buffer(1:nbytes)) map(from: chunk_cs(1:nthreads), chunk_sz(1:nthreads)) &
    !$omp& private(first, last, local_cs, pos, bit)
    do tid = 1, nthreads
      first = int(tid - 1, int64) * bpt + 1_int64
      if (tid /= nthreads) then
        last = first + bpt - 1_int64
      else
        last = nbytes
      end if

      local_cs = all_ones
      do pos = first, last
        local_cs = ieor(local_cs, int(iand(buffer(pos), 255_int32), int64))
        do bit = 1, 8
          if (iand(local_cs, 1_int64) /= 0_int64) then
            local_cs = ieor(shiftr(local_cs, 1), crc64_poly)
          else
            local_cs = shiftr(local_cs, 1)
          end if
        end do
      end do
      chunk_cs(tid) = ieor(local_cs, all_ones)
      chunk_sz(tid) = last - first + 1_int64
    end do
    !$omp end target teams distribute parallel do

    cs = chunk_cs(1)
    do tid = 2, nthreads
      cs = crc64_combine(cs, chunk_cs(tid), chunk_sz(tid))
    end do

    deallocate(chunk_cs, chunk_sz)
  end function crc64_omp

end module crc64_mod

program main
  use iso_fortran_env, only: int32, int64, real64
  use omp_lib, only: omp_get_wtime
  use crc64_mod
  implicit none

  integer :: argc, ntests, seed, max_test_length, ntest
  integer(int64) :: rng_state, test_length, div_pt, tlend
  integer(int32), allocatable :: buffer(:)
  integer(int32) :: check_bytes(8)
  integer(int64) :: cs, csc, cs1, cs2
  real(real64) :: rnd, start_time, end_time, b_time, tot_time, tot_bytes
  character(len=64) :: arg
  character(len=4) :: check1, check2
  integer(int64) :: i

  ntests = 10
  seed = 5
  max_test_length = 2097152

  argc = command_argument_count()
  if (argc > 0) then
    call get_command_argument(1, arg)
    read(arg, *) ntests
  end if
  if (argc > 1) then
    call get_command_argument(2, arg)
    read(arg, *) seed
  end if
  if (argc > 2) then
    call get_command_argument(3, arg)
    read(arg, *) max_test_length
  end if

  write(*,'("Running ",I0," tests with seed ",I0)') ntests, seed

  call srand48_state(int(seed, int32), rng_state)
  tot_time = 0.0_real64
  tot_bytes = 0.0_real64
  tlend = 8_int64

  do ntest = 1, ntests
    call next_drand48(rng_state, rnd)
    test_length = int(real(max_test_length, real64) * (rnd + 1.0_real64), int64)
    allocate(buffer(test_length + tlend))

    do i = 1, test_length
      call next_drand48(rng_state, rnd)
      buffer(i) = int(255.0_real64 * rnd, int32)
    end do

    start_time = omp_get_wtime()
    cs = crc64_omp(buffer, test_length)
    end_time = omp_get_wtime()
    b_time = end_time - start_time

    if (ntest > 1) then
      tot_time = tot_time + b_time
      tot_bytes = tot_bytes + real(test_length, real64)
    end if

    buffer(test_length + 1:test_length + tlend) = 0_int32
    call crc64_invert(cs, check_bytes)
    buffer(test_length + 1:test_length + tlend) = check_bytes

    csc = crc64_serial(buffer, 1_int64, test_length + tlend)
    if (csc == all_ones) then
      check1 = "pass"
    else
      check1 = "fail"
    end if

    call next_drand48(rng_state, rnd)
    div_pt = int(real(test_length, real64) * rnd, int64)
    if (div_pt > 0_int64) then
      cs1 = crc64_serial(buffer, 1_int64, div_pt)
    else
      cs1 = crc64_serial(buffer, 1_int64, 0_int64)
    end if
    cs2 = crc64_serial(buffer, div_pt + 1_int64, test_length)
    csc = crc64_combine(cs1, cs2, test_length - div_pt)
    if (csc == cs) then
      check2 = "pass"
    else
      check2 = "fail"
    end if

    write(*,'(I0,1X,I0,1X,A,1X,A)') ntest, test_length, trim(check1), trim(check2)
    deallocate(buffer)
  end do

  if (tot_time > 0.0_real64) then
    write(*,'(G0,1X,"MB/s")') (tot_bytes / (1024.0_real64 * 1024.0_real64)) / tot_time
  else
    write(*,'("inf MB/s")')
  end if

end program main
