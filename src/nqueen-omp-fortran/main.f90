program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  integer, parameter :: queens_block_size = 128
  integer, parameter :: empty = -1
  integer, parameter :: max_depth = 20
  integer, parameter :: max_board = 12

  integer :: size, initial_depth, repeat
  integer(int64) :: tree_size, n_explorers, qtd_sols_global
  integer(int64), allocatable :: controls(:), vector_of_tree_size(:), solutions(:)
  integer(int32), allocatable :: boards(:, :)
  integer(int64) :: i

  if (command_argument_count() /= 3) then
    print '(A)', 'Usage: ./main <size> <initial depth> <repeat>'
    stop 1
  end if

  size = read_arg(1)
  initial_depth = read_arg(2)
  repeat = read_arg(3)

  write(*,'(/,A,I0,A,I0,A)', advance='no') '### Initial depth: ', initial_depth, ' - Size: ', size, ':'

  if (size > max_depth .or. initial_depth > max_board) then
    print '(A)', ''
    print '(A)', 'Error: unsupported board size or initial depth'
    stop 1
  end if

  call count_prefixes(size, initial_depth, tree_size, n_explorers)
  allocate(controls(n_explorers), boards(0:max_board - 1, n_explorers))
  allocate(vector_of_tree_size(n_explorers), solutions(n_explorers))
  boards = empty
  controls = 0_int64
  vector_of_tree_size = 0_int64
  solutions = 0_int64

  call build_prefixes(size, initial_depth, tree_size, controls, boards)
  call nqueens(size, initial_depth, int(n_explorers), controls, boards, vector_of_tree_size, solutions, repeat)

  print '(A,I0)', ''
  print '(A,I0)', 'Tree size: ', tree_size

  qtd_sols_global = 0_int64
  do i = 1, n_explorers
    if (solutions(i) > 0_int64) qtd_sols_global = qtd_sols_global + solutions(i)
    if (vector_of_tree_size(i) > 0_int64) tree_size = tree_size + vector_of_tree_size(i)
  end do

  print '(A,I0,/,A,I0)', 'Number of solutions found: ', qtd_sols_global, 'Tree size: ', tree_size
  if (size == 15 .and. initial_depth == 7) then
    if (qtd_sols_global == 2279184_int64 .and. tree_size == 171129071_int64) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
    end if
  end if

  deallocate(controls, boards, vector_of_tree_size, solutions)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  logical function host_still_legal(board, r)
    integer(int32), intent(in) :: board(0:)
    integer, intent(in) :: r
    integer :: i, ld, rd

    host_still_legal = .true.
    do i = 0, r - 1
      if (board(i) == board(r)) then
        host_still_legal = .false.
        return
      end if
    end do

    ld = board(r)
    rd = board(r)
    do i = r - 1, 0, -1
      ld = ld - 1
      rd = rd + 1
      if (board(i) == ld .or. board(i) == rd) then
        host_still_legal = .false.
        return
      end if
    end do
  end function host_still_legal

  subroutine count_prefixes(size, initial_depth, tree_size, num_sol)
    integer, intent(in) :: size, initial_depth
    integer(int64), intent(out) :: tree_size, num_sol
    integer(int64) :: flag, bit_test
    integer(int32) :: vertice(0:max_depth - 1)
    integer :: nivel

    flag = 0_int64
    tree_size = 0_int64
    num_sol = 0_int64
    vertice = empty
    nivel = 0

    do while (nivel >= 0)
      vertice(nivel) = vertice(nivel) + 1
      bit_test = shiftl(1_int64, vertice(nivel))
      if (vertice(nivel) == size) then
        vertice(nivel) = empty
      else if (host_still_legal(vertice, nivel) .and. (iand(flag, bit_test) == 0_int64)) then
        flag = ior(flag, bit_test)
        nivel = nivel + 1
        tree_size = tree_size + 1_int64
        if (nivel == initial_depth) then
          num_sol = num_sol + 1_int64
        else
          cycle
        end if
      else
        cycle
      end if

      nivel = nivel - 1
      if (nivel >= 0) flag = iand(flag, not(shiftl(1_int64, vertice(nivel))))
    end do
  end subroutine count_prefixes

  subroutine build_prefixes(size, initial_depth, tree_size, controls, boards)
    integer, intent(in) :: size, initial_depth
    integer(int64), intent(out) :: tree_size
    integer(int64), intent(out) :: controls(:)
    integer(int32), intent(out) :: boards(0:, :)
    integer(int64) :: flag, bit_test, num_sol
    integer(int32) :: vertice(0:max_depth - 1)
    integer :: nivel, j

    flag = 0_int64
    tree_size = 0_int64
    num_sol = 0_int64
    vertice = empty
    nivel = 0

    do while (nivel >= 0)
      vertice(nivel) = vertice(nivel) + 1
      bit_test = shiftl(1_int64, vertice(nivel))
      if (vertice(nivel) == size) then
        vertice(nivel) = empty
      else if (host_still_legal(vertice, nivel) .and. (iand(flag, bit_test) == 0_int64)) then
        flag = ior(flag, bit_test)
        nivel = nivel + 1
        tree_size = tree_size + 1_int64
        if (nivel == initial_depth) then
          num_sol = num_sol + 1_int64
          controls(num_sol) = flag
          do j = 0, initial_depth - 1
            boards(j, num_sol) = vertice(j)
          end do
        else
          cycle
        end if
      else
        cycle
      end if

      nivel = nivel - 1
      if (nivel >= 0) flag = iand(flag, not(shiftl(1_int64, vertice(nivel))))
    end do
  end subroutine build_prefixes

  subroutine nqueens(size, initial_depth, n_explorers, controls, boards, vector_of_tree_size, solutions, repeat)
    integer, intent(in) :: size, initial_depth, n_explorers, repeat
    integer(int64), intent(in) :: controls(:)
    integer(int32), intent(in) :: boards(0:, :)
    integer(int64), intent(out) :: vector_of_tree_size(:), solutions(:)
    integer :: rep
    real(real64) :: start_time, end_time

    print '(A)', ''
    print '(A)', '### Regular BP-DFS search. ###'

    !$omp target data map(to: controls(1:n_explorers), boards(0:max_board - 1, 1:n_explorers)) &
    !$omp& map(from: vector_of_tree_size(1:n_explorers), solutions(1:n_explorers))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      call bp_queens_root_dfs(size, n_explorers, initial_depth, controls, boards, vector_of_tree_size, solutions)
    end do
    end_time = omp_get_wtime()
    print '(A,F0.6,A)', 'Average kernel execution time: ', (end_time - start_time) / real(repeat, real64), ' (s)'
    !$omp end target data
  end subroutine nqueens

  subroutine bp_queens_root_dfs(size, n_explorers, initial_depth, controls, boards, vector_of_tree_size, solutions)
    integer, intent(in) :: size, n_explorers, initial_depth
    integer(int64), intent(in) :: controls(:)
    integer(int32), intent(in) :: boards(0:, :)
    integer(int64), intent(out) :: vector_of_tree_size(:), solutions(:)
    integer :: idx, i, depth, depth_global, ld, rd
    integer(int64) :: flag, bit_test, qtd_solutions_thread, tree_size
    integer(int32) :: vertice(0:max_depth - 1)
    logical :: safe

    !$omp target teams distribute parallel do thread_limit(queens_block_size) &
    !$omp& private(idx, i, depth, depth_global, ld, rd, flag, bit_test, qtd_solutions_thread, tree_size, vertice, safe)
    do idx = 1, n_explorers
      flag = controls(idx)
      vertice = empty
      do i = 0, initial_depth - 1
        vertice(i) = boards(i, idx)
      end do

      depth_global = initial_depth
      depth = depth_global
      qtd_solutions_thread = 0_int64
      tree_size = 0_int64

      do while (depth >= depth_global)
        vertice(depth) = vertice(depth) + 1
        bit_test = shiftl(1_int64, vertice(depth))
        if (vertice(depth) == size) then
          vertice(depth) = empty
        else if (iand(flag, bit_test) == 0_int64) then
          safe = .true.
          do i = 0, depth - 1
            if (vertice(i) == vertice(depth)) safe = .false.
          end do
          ld = vertice(depth)
          rd = vertice(depth)
          do i = depth - 1, 0, -1
            ld = ld - 1
            rd = rd + 1
            if (vertice(i) == ld .or. vertice(i) == rd) safe = .false.
          end do

          if (safe) then
            tree_size = tree_size + 1_int64
            flag = ior(flag, bit_test)
            depth = depth + 1
            if (depth == size) then
              qtd_solutions_thread = qtd_solutions_thread + 1_int64
            else
              cycle
            end if
          else
            cycle
          end if
        else
          cycle
        end if

        depth = depth - 1
        if (depth >= 0) flag = iand(flag, not(shiftl(1_int64, vertice(depth))))
      end do

      solutions(idx) = qtd_solutions_thread
      vector_of_tree_size(idx) = tree_size
    end do
    !$omp end target teams distribute parallel do
  end subroutine bp_queens_root_dfs

end program main
