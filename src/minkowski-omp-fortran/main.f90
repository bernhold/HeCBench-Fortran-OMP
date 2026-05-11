program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: m_size = 512_int32 * 8_int32
  integer(int32), parameter :: m_dim = m_size / 8_int32
  integer(int32), parameter :: n_dim = m_size / 4_int32
  integer(int32), parameter :: k_dim = m_size / 2_int32
  real(real32), parameter :: value_tolerance = 1.0e-5_real32

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

  character(len=256) :: arg0, arg1
  integer(int32) :: repeat, i, j, m, iter
  real(real32), allocatable :: a_host(:, :), b_host(:, :), c_host(:, :), c_back(:, :)
  real(real32) :: p_value, one_over_p
  real(real64) :: start_time, end_time, avg_seconds

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    print '(3A)', 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) repeat
  if (repeat <= 0_int32) then
    stop 1
  end if

  allocate(a_host(1:m_dim, 1:n_dim))
  allocate(b_host(1:n_dim, 1:k_dim))
  allocate(c_host(1:m_dim, 1:k_dim))
  allocate(c_back(1:m_dim, 1:k_dim))

  a_host = 1.0_real32 / real(n_dim, real32)

  call c_srand(123_c_int)
  do i = 1, n_dim
    do j = 1, k_dim
      b_host(i, j) = real(modulo(c_rand(), 256_c_int), real32)
    end do
  end do

  call normalize_columns(b_host)

  write(*, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'Problem size: c(', m_dim, ',', k_dim, ') = a(', m_dim, ',', &
                                                 n_dim, ') * b(', n_dim, ',', k_dim, ')'

  !$omp target data map(to: a_host(1:m_dim, 1:n_dim), b_host(1:n_dim, 1:k_dim)) &
  !$omp& map(alloc: c_back(1:m_dim, 1:k_dim))
  do m = 1, 4
    print '(A,I0)', 'Minkowski distance with p = ', m
    p_value = real(m, real32)
    one_over_p = 1.0_real32 / p_value

    start_time = omp_get_wtime()
    do iter = 1, repeat
      !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(j)
      do i = 1, m_dim
        do j = 1, k_dim
          c_back(i, j) = minkowski_element(a_host, b_host, i, j, p_value, one_over_p)
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()

    avg_seconds = (end_time - start_time) / real(repeat, real64)
    print '(A,F0.6,A)', 'Average kernel execution time: ', avg_seconds, ' (s)'

    !$omp target update from(c_back(1:m_dim, 1:k_dim))
    call verify_result(a_host, b_host, c_host, c_back, p_value, one_over_p)
  end do
  !$omp end target data

contains

  subroutine normalize_columns(matrix)
    real(real32), intent(inout) :: matrix(1:, 1:)
    integer(int32) :: row, col
    real(real32) :: sum_value

    do col = 1, size(matrix, 2)
      sum_value = 0.0_real32
      do row = 1, size(matrix, 1)
        sum_value = sum_value + matrix(row, col)
      end do
      do row = 1, size(matrix, 1)
        matrix(row, col) = matrix(row, col) / sum_value
      end do
    end do
  end subroutine normalize_columns

  real(real32) function minkowski_element(a_matrix, b_matrix, row_idx, col_idx, p_value, one_over_p) result(value)
    real(real32), intent(in) :: a_matrix(1:, 1:), b_matrix(1:, 1:)
    integer(int32), intent(in) :: row_idx, col_idx
    real(real32), intent(in) :: p_value, one_over_p
    integer(int32) :: k
    real(real32) :: sum_value

    sum_value = 0.0_real32
    do k = 1, size(a_matrix, 2)
      sum_value = sum_value + abs(a_matrix(row_idx, k) - b_matrix(k, col_idx)) ** p_value
    end do
    value = sum_value ** one_over_p
  end function minkowski_element

  logical function value_same(lhs, rhs) result(ok)
    real(real32), intent(in) :: lhs, rhs
    ok = abs(lhs - rhs) <= value_tolerance
  end function value_same

  subroutine verify_result(a_matrix, b_matrix, c_matrix, c_device, p_value, one_over_p)
    real(real32), intent(in) :: a_matrix(1:, 1:), b_matrix(1:, 1:)
    real(real32), intent(inout) :: c_matrix(1:, 1:)
    real(real32), intent(in) :: c_device(1:, 1:)
    real(real32), intent(in) :: p_value, one_over_p
    integer(int32) :: i, j, k, print_count
    logical :: mismatch_found

    c_matrix = 0.0_real32
    do i = 1, size(c_matrix, 1)
      do k = 1, size(a_matrix, 2)
        do j = 1, size(c_matrix, 2)
          c_matrix(i, j) = c_matrix(i, j) + abs(a_matrix(i, k) - b_matrix(k, j)) ** p_value
        end do
      end do
    end do

    do i = 1, size(c_matrix, 1)
      do j = 1, size(c_matrix, 2)
        c_matrix(i, j) = c_matrix(i, j) ** one_over_p
      end do
    end do

    mismatch_found = .false.
    print_count = 0

    do i = 1, size(c_matrix, 1)
      do j = 1, size(c_matrix, 2)
        if (.not. value_same(c_device(i, j), c_matrix(i, j))) then
          write(*, '(A,I0,A,I0,A,F0.6,A,F0.6)') 'Fail - The result is incorrect for element: [', i - 1, ', ', &
                                                 j - 1, '], expected: ', c_matrix(i, j), ', but found: ', c_device(i, j)
          mismatch_found = .true.
          print_count = print_count + 1
          if (print_count == 5) exit
        end if
      end do
      if (print_count == 5) exit
    end do

    if (.not. mismatch_found) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
    end if
  end subroutine verify_result

end program main
