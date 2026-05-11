program main
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  integer :: no_of_nodes, edge_list_size, source
  integer, allocatable :: node_start(:), node_edges(:), graph_edges(:)
  integer, allocatable :: graph_mask(:), updating_graph_mask(:), graph_visited(:)
  integer, allocatable :: cost(:), cost_ref(:)
  character(len=4096) :: input_file

  if (command_argument_count() /= 1) then
    call usage()
    stop 0
  end if

  call get_command_argument(1, input_file)
  print '(A)', 'Reading File'
  call read_graph(trim(input_file), no_of_nodes, edge_list_size, source, node_start, node_edges, graph_edges)

  allocate(graph_mask(no_of_nodes), updating_graph_mask(no_of_nodes), graph_visited(no_of_nodes))
  allocate(cost(no_of_nodes), cost_ref(no_of_nodes))

  graph_mask = 0
  updating_graph_mask = 0
  graph_visited = 0
  cost = -1
  cost_ref = -1
  source = 0
  graph_mask(source + 1) = 1
  graph_visited(source + 1) = 1
  cost(source + 1) = 0
  cost_ref(source + 1) = 0

  print '(A,I0,A)', 'run bfs (#nodes = ', no_of_nodes, ') on device'
  call run_bfs_gpu(no_of_nodes, node_start, node_edges, edge_list_size, graph_edges, &
      graph_mask, updating_graph_mask, graph_visited, cost)

  print '(A,I0,A)', 'run bfs (#nodes = ', no_of_nodes, ') on host (cpu) '
  graph_mask = 0
  updating_graph_mask = 0
  graph_visited = 0
  graph_mask(source + 1) = 1
  graph_visited(source + 1) = 1
  call run_bfs_cpu(no_of_nodes, node_start, node_edges, graph_edges, &
      graph_mask, updating_graph_mask, graph_visited, cost_ref)

  call compare_results(cost_ref, cost)

  deallocate(node_start, node_edges, graph_edges)
  deallocate(graph_mask, updating_graph_mask, graph_visited, cost, cost_ref)

contains

  subroutine usage()
    character(len=4096) :: program_name
    call get_command_argument(0, program_name)
    write(0, '(A,A,A)') 'Usage: ', trim(program_name), ' <input_file>'
  end subroutine usage

  subroutine read_graph(path, no_of_nodes, edge_list_size, source, node_start, node_edges, graph_edges)
    character(len=*), intent(in) :: path
    integer, intent(out) :: no_of_nodes, edge_list_size, source
    integer, allocatable, intent(out) :: node_start(:), node_edges(:), graph_edges(:)
    integer :: unit, ios, idx, edge_id, edge_cost

    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print '(A,A)', 'Error Reading graph file ', trim(path)
      stop 1
    end if

    read(unit, *) no_of_nodes
    allocate(node_start(no_of_nodes), node_edges(no_of_nodes))
    do idx = 1, no_of_nodes
      read(unit, *) node_start(idx), node_edges(idx)
    end do
    read(unit, *) source
    source = 0
    read(unit, *) edge_list_size
    allocate(graph_edges(edge_list_size))
    do idx = 1, edge_list_size
      read(unit, *) edge_id, edge_cost
      graph_edges(idx) = edge_id
    end do
    close(unit)
  end subroutine read_graph

  subroutine run_bfs_cpu(no_of_nodes, node_start, node_edges, graph_edges, graph_mask, updating_graph_mask, graph_visited, cost_ref)
    integer, intent(in) :: no_of_nodes
    integer, intent(in) :: node_start(:), node_edges(:), graph_edges(:)
    integer, intent(inout) :: graph_mask(:), updating_graph_mask(:), graph_visited(:), cost_ref(:)
    integer :: tid, edge_idx, graph_id, stop_flag
    do
      stop_flag = 0
      do tid = 0, no_of_nodes - 1
        if (graph_mask(tid + 1) == 1) then
          graph_mask(tid + 1) = 0
          do edge_idx = node_start(tid + 1), node_start(tid + 1) + node_edges(tid + 1) - 1
            graph_id = graph_edges(edge_idx + 1)
            if (graph_visited(graph_id + 1) == 0) then
              cost_ref(graph_id + 1) = cost_ref(tid + 1) + 1
              updating_graph_mask(graph_id + 1) = 1
            end if
          end do
        end if
      end do
      do tid = 0, no_of_nodes - 1
        if (updating_graph_mask(tid + 1) == 1) then
          graph_mask(tid + 1) = 1
          graph_visited(tid + 1) = 1
          stop_flag = 1
          updating_graph_mask(tid + 1) = 0
        end if
      end do
      if (stop_flag == 0) exit
    end do
  end subroutine run_bfs_cpu

  subroutine run_bfs_gpu(no_of_nodes, node_start, node_edges, edge_list_size, graph_edges, graph_mask, updating_graph_mask, graph_visited, cost)
    integer, intent(in) :: no_of_nodes, edge_list_size
    integer, intent(in) :: node_start(:), node_edges(:), graph_edges(:)
    integer, intent(inout) :: graph_mask(:), updating_graph_mask(:), graph_visited(:), cost(:)
    integer :: tid, edge_idx, graph_id
    integer :: over(1)
    real(real64) :: start_time, end_time, elapsed

    elapsed = 0.0_real64
    !$omp target data map(to: node_start(1:no_of_nodes), node_edges(1:no_of_nodes), graph_edges(1:edge_list_size)) &
    !$omp& map(to: graph_mask(1:no_of_nodes), updating_graph_mask(1:no_of_nodes), graph_visited(1:no_of_nodes)) &
    !$omp& map(alloc: over(1:1)) map(tofrom: cost(1:no_of_nodes))
    do
      over(1) = 0
      !$omp target update to(over(1:1))
      start_time = omp_get_wtime()

      !$omp target teams distribute parallel do thread_limit(256)
      do tid = 0, no_of_nodes - 1
        if (graph_mask(tid + 1) /= 0) then
          graph_mask(tid + 1) = 0
          do edge_idx = node_start(tid + 1), node_start(tid + 1) + node_edges(tid + 1) - 1
            graph_id = graph_edges(edge_idx + 1)
            if (graph_visited(graph_id + 1) == 0) then
              cost(graph_id + 1) = cost(tid + 1) + 1
              updating_graph_mask(graph_id + 1) = 1
            end if
          end do
        end if
      end do
      !$omp end target teams distribute parallel do

      !$omp target teams distribute parallel do thread_limit(256)
      do tid = 0, no_of_nodes - 1
        if (updating_graph_mask(tid + 1) /= 0) then
          graph_mask(tid + 1) = 1
          graph_visited(tid + 1) = 1
          over(1) = 1
          updating_graph_mask(tid + 1) = 0
        end if
      end do
      !$omp end target teams distribute parallel do

      end_time = omp_get_wtime()
      elapsed = elapsed + end_time - start_time
      !$omp target update from(over(1:1))
      if (over(1) == 0) exit
    end do
    !$omp end target data

    print '(A,F0.6,A)', 'Total kernel execution time : ', elapsed * 1.0e6_real64, ' (us)'
  end subroutine run_bfs_gpu

  subroutine compare_results(cpu_results, gpu_results)
    integer, intent(in) :: cpu_results(:), gpu_results(:)
    if (all(cpu_results == gpu_results)) then
      print '(A)', 'Passed'
    else
      print '(A)', 'Failed'
    end if
  end subroutine compare_results

end program main
