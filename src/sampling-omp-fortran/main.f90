program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: ncases = 3
  integer, parameter :: exact_cases(ncases) = [1000, 0, 1000]
  integer, parameter :: sampled_cases(ncases) = [0, 1000, 1000]
  integer, parameter :: ncols_cases(ncases) = [2000, 2000, 2000]
  integer, parameter :: background_cases(ncases) = [10, 10, 10]
  integer, parameter :: max_samples_cases(ncases) = [11, 11, 11]
  integer :: repeat, case_id, r
  real(real64) :: total_time

  if (command_argument_count() /= 1) then
    print '(A)', 'Usage: ./main <repeat>'
    stop 1
  end if
  repeat = read_arg(1)
  if (repeat <= 0) then
    print '(A)', 'invalid repeat'
    stop 1
  end if

  do case_id = 1, ncases
    total_time = 0.0_real64
    do r = 1, repeat
      call run_case(exact_cases(case_id), sampled_cases(case_id), ncols_cases(case_id), &
                    background_cases(case_id), max_samples_cases(case_id), total_time)
    end do
    write(*, '(A,F0.6,A)') 'Average execution time of kernels: ', (total_time * 1.0e6_real64) / real(repeat, real64), ' (us)'
  end do

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=128) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine run_case(nrows_exact, nrows_sampled, ncols, nrows_background, max_samples, total_time)
    integer, intent(in) :: nrows_exact, nrows_sampled, ncols, nrows_background, max_samples
    real(real64), intent(inout) :: total_time
    integer :: nrows_x, i, j
    real(real32) :: sent_value
    real(real32), allocatable :: background(:), observation(:), x(:), dataset(:)
    integer, allocatable :: nsamples(:)
    real(real64) :: start_time, end_time

    nrows_x = nrows_exact + nrows_sampled
    allocate(background(nrows_background * ncols), observation(ncols), x(nrows_x * ncols))
    allocate(dataset(nrows_x * nrows_background * ncols))
    allocate(nsamples(max(1, nrows_sampled / 2)))

    sent_value = real(nrows_x * nrows_background * ncols * 100, real32)
    observation = sent_value
    do i = 0, nrows_background - 1
      do j = 0, ncols - 1
        background(i * ncols + j + 1) = real(i * 2 + 1, real32)
      end do
    end do

    x = 0.0_real32
    do i = 0, nrows_exact - 1
      x(i * ncols + i + 1) = 1.0_real32
      x(i * ncols + i + 2) = 1.0_real32
    end do

    do i = 0, nrows_sampled / 2 - 1
      nsamples(i + 1) = max_samples - mod(i, 2)
    end do
    dataset = 0.0_real32

    start_time = omp_get_wtime()
    !$omp target data map(to: background(1:nrows_background*ncols), observation(1:ncols), nsamples(1:max(1,nrows_sampled/2))) &
    !$omp& map(tofrom: x(1:nrows_x*ncols)) map(from: dataset(1:nrows_x*nrows_background*ncols))
      if (nrows_exact > 0) then
        call fill_exact(background, observation, x, dataset, nrows_exact, ncols, nrows_background)
      end if
      if (nrows_sampled > 0) then
        call fill_sampled(background, observation, nsamples, x, dataset, nrows_exact, nrows_sampled, ncols, nrows_background)
      end if
    !$omp end target data
    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)

    call validate_case(x, dataset, nsamples, sent_value, nrows_exact, nrows_sampled, ncols, nrows_background)
    deallocate(background, observation, x, dataset, nsamples)
  end subroutine run_case

  subroutine fill_exact(background, observation, x, dataset, nrows_exact, ncols, nrows_background)
    real(real32), intent(in) :: background(:), observation(:), x(:)
    real(real32), intent(out) :: dataset(:)
    integer, intent(in) :: nrows_exact, ncols, nrows_background
    integer :: row, bg, col, curr_x, out_idx

    !$omp target teams distribute parallel do collapse(3) thread_limit(256) private(row, bg, col, curr_x, out_idx)
    do row = 0, nrows_exact - 1
      do bg = 0, nrows_background - 1
        do col = 0, ncols - 1
          curr_x = int(x(row * ncols + col + 1))
          out_idx = (row * nrows_background + bg) * ncols + col + 1
          if (curr_x == 0) then
            dataset(out_idx) = background(bg * ncols + col + 1)
          else
            dataset(out_idx) = observation(col + 1)
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine fill_exact

  subroutine fill_sampled(background, observation, nsamples, x, dataset, nrows_exact, nrows_sampled, ncols, nrows_background)
    real(real32), intent(in) :: background(:), observation(:)
    integer, intent(in) :: nsamples(:)
    real(real32), intent(inout) :: x(:)
    real(real32), intent(inout) :: dataset(:)
    integer, intent(in) :: nrows_exact, nrows_sampled, ncols, nrows_background
    integer :: pair_id, bg, col, row0, row1, sample_count, curr_x, out_idx

    !$omp target teams distribute parallel do collapse(3) thread_limit(256) private(pair_id, bg, col, row0, row1, sample_count, curr_x, out_idx)
    do pair_id = 0, nrows_sampled / 2 - 1
      do bg = 0, nrows_background - 1
        do col = 0, ncols - 1
          row0 = nrows_exact + 2 * pair_id
          row1 = row0 + 1
          sample_count = nsamples(pair_id + 1)
          if (col < sample_count) then
            x(row0 * ncols + col + 1) = 1.0_real32
            x(row1 * ncols + col + 1) = 0.0_real32
          else
            x(row0 * ncols + col + 1) = 0.0_real32
            x(row1 * ncols + col + 1) = 1.0_real32
          end if

          curr_x = int(x(row0 * ncols + col + 1))
          out_idx = (row0 * nrows_background + bg) * ncols + col + 1
          if (curr_x == 0) then
            dataset(out_idx) = background(bg * ncols + col + 1)
          else
            dataset(out_idx) = observation(col + 1)
          end if

          curr_x = int(x(row1 * ncols + col + 1))
          out_idx = (row1 * nrows_background + bg) * ncols + col + 1
          if (curr_x == 1) then
            dataset(out_idx) = observation(col + 1)
          else
            dataset(out_idx) = background(bg * ncols + col + 1)
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine fill_sampled

  subroutine validate_case(x, dataset, nsamples, sent_value, nrows_exact, nrows_sampled, ncols, nrows_background)
    real(real32), intent(in) :: x(:), dataset(:), sent_value
    integer, intent(in) :: nsamples(:), nrows_exact, nrows_sampled, ncols, nrows_background
    logical :: test_sampled_x, test_scatter_exact, test_scatter_sampled
    integer :: i, j, k, counter, sample_idx, compliment_ctr

    test_sampled_x = .true.
    sample_idx = 1
    do i = nrows_exact * ncols + 1, (nrows_exact + nrows_sampled) * ncols, 2 * ncols
      if (nrows_sampled <= 0) exit
      counter = 0
      do k = i, i + ncols - 1
        if (x(k) == 1.0_real32) counter = counter + 1
      end do
      test_sampled_x = test_sampled_x .and. (counter == nsamples(sample_idx))
      counter = 0
      do k = i + ncols, i + 2 * ncols - 1
        if (x(k) == 1.0_real32) counter = counter + 1
      end do
      test_sampled_x = test_sampled_x .and. (counter == ncols - nsamples(sample_idx))
      sample_idx = sample_idx + 1
    end do

    test_scatter_exact = .true.
    do i = 0, nrows_exact - 1
      do j = i * nrows_background * ncols + 1, (i + 1) * nrows_background * ncols, ncols
        counter = 0
        do k = j, j + ncols - 1
          if (dataset(k) == sent_value) counter = counter + 1
        end do
        test_scatter_exact = test_scatter_exact .and. (counter == 2)
        if (.not. test_scatter_exact) then
          write(*, '(A,I0,A)') 'test_scatter_exact counter failed with: ', counter, ', expected value was 2.'
          exit
        end if
      end do
      if (.not. test_scatter_exact) exit
    end do

    test_scatter_sampled = .true.
    compliment_ctr = 0
    do i = nrows_exact, nrows_exact + nrows_sampled / 2 - 1
      do j = (i + compliment_ctr) * nrows_background * ncols + 1, (i + compliment_ctr + 1) * nrows_background * ncols, ncols
        counter = 0
        do k = j, j + ncols - 1
          if (dataset(k) == sent_value) counter = counter + 1
        end do
        test_scatter_sampled = test_scatter_sampled .and. (counter == nsamples(i - nrows_exact + 1))
        if (.not. test_scatter_sampled) then
          write(*, '(A,I0,A,I0,A)') 'test_scatter_sampled counter failed with: ', counter, &
            ', expected value was ', nsamples(i - nrows_exact + 1), '.'
          exit
        end if
      end do

      compliment_ctr = compliment_ctr + 1
      do j = (i + compliment_ctr) * nrows_background * ncols + 1, (i + compliment_ctr + 1) * nrows_background * ncols, ncols
        counter = 0
        do k = j, j + ncols - 1
          if (dataset(k) == sent_value) counter = counter + 1
        end do
        test_scatter_sampled = test_scatter_sampled .and. (counter == ncols - nsamples(i - nrows_exact + 1))
        if (.not. test_scatter_sampled) then
          write(*, '(A,I0,A,I0,A)') 'test_scatter_sampled counter failed with: ', counter, &
            ', expected value was ', ncols - nsamples(i - nrows_exact + 1), '.'
          exit
        end if
      end do
      if (.not. test_scatter_sampled) exit
    end do

    if (.not. test_sampled_x) print '(A)', 'test_sampled_X failed'
  end subroutine validate_case

end program main
