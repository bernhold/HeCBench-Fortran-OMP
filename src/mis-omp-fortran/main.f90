program mis
  use iso_fortran_env, only: int8, int32, int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: threads_per_block = 256
  integer, parameter :: blocks = 24
  integer(int32), parameter :: in_status = int(z'000000FE', int32)
  integer(int32), parameter :: out_status = 0_int32

  character(len=512) :: graph_file
  character(len=64) :: arg
  integer :: argc, repeat
  integer(int32) :: nodes, edges
  integer(int32), allocatable :: nidx(:), nlist(:), eweight(:)
  integer(int32), allocatable :: nstatus(:)

  argc = command_argument_count()
  print '(A)', 'ECL-MIS v1.3 (main.cpp)'
  print '(A)', 'Copyright 2017-2020 Texas State University'

  if (argc /= 2) then
    call get_command_argument(0, arg)
    write(*, '(A,A,A)') 'USAGE: ', trim(arg), ' <input_file_name> <repeat>'
    print *
    stop 255
  end if

  call get_command_argument(1, graph_file)
  call get_command_argument(2, arg)
  read(arg, *) repeat

  call read_ecl_graph(trim(graph_file), nodes, edges, nidx, nlist, eweight)
  print '(A,I0,A,I0,A,A,A)', 'configuration: ', nodes, ' nodes and ', edges, &
      ' edges (', trim(graph_file), ')'
  print '(A,F0.2,A)', 'average degree: ', real(edges, real64) / real(nodes, real64), &
      ' edges per node'

  allocate(nstatus(nodes))
  call compute_mis(repeat, nodes, edges, nidx, nlist, nstatus)
  call verify_mis(nodes, nidx, nlist, nstatus)

  deallocate(nidx, nlist, nstatus)
  if (allocated(eweight)) deallocate(eweight)

contains

  subroutine read_ecl_graph(path, nodes, edges, nidx, nlist, eweight)
    character(len=*), intent(in) :: path
    integer(int32), intent(out) :: nodes, edges
    integer(int32), allocatable, intent(out) :: nidx(:), nlist(:), eweight(:)

    integer :: unit, ios
    integer(int64) :: file_size, pos_after_graph, remaining_bytes

    open(newunit=unit, file=path, access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(0, '(A,A)') 'ERROR: could not open file ', trim(path)
      print *
      stop 255
    end if

    inquire(unit=unit, size=file_size)
    read(unit) nodes
    read(unit) edges
    if (nodes < 1 .or. edges < 0) then
      write(0, '(A)') 'ERROR: node or edge count too low'
      print *
      close(unit)
      stop 255
    end if

    allocate(nidx(nodes + 1), nlist(edges))
    read(unit, iostat=ios) nidx
    if (ios /= 0) write(0, '(A)') 'ERROR: failed to read neighbor index list'
    read(unit, iostat=ios) nlist
    if (ios /= 0) write(0, '(A)') 'ERROR: failed to read neighbor list'

    inquire(unit=unit, pos=pos_after_graph)
    remaining_bytes = file_size - pos_after_graph + 1_int64
    if (remaining_bytes >= int(edges, int64) * 4_int64) then
      allocate(eweight(edges))
      read(unit, iostat=ios) eweight
      if (ios /= 0) write(0, '(A)') 'ERROR: failed to read edge weights'
    else
      allocate(eweight(0))
    end if
    close(unit)
  end subroutine read_ecl_graph

  subroutine compute_mis(repeat, nodes, edges, nidx, nlist, nstatus)
    integer, intent(in) :: repeat
    integer(int32), intent(in) :: nodes, edges
    integer(int32), intent(in) :: nidx(:), nlist(:)
    integer(int32), intent(out) :: nstatus(:)

    real(real32) :: avg, scaledavg, runtime
    real(real64) :: start_time, end_time
    integer :: pass

    avg = real(edges, real32) / real(nodes, real32)
    scaledavg = real((in_status / 2_int32) - 1_int32, real32) * avg

    start_time = omp_get_wtime()
    !$omp target data map(to: nidx(1:nodes+1), nlist(1:edges)) &
    !$omp& map(from: nstatus(1:nodes))
    do pass = 1, 100
      call initialize_status(nodes, nidx, nstatus, avg, scaledavg)
      call converge_status(nodes, nidx, nlist, nstatus)
    end do
    !$omp end target data
    end_time = omp_get_wtime()

    runtime = real(end_time - start_time, real32) / real(repeat, real32)
    print '(A,F0.6,A)', 'compute time: ', runtime, ' s'
    print '(A,F0.6,A)', 'throughput: ', real(nodes, real32) * 0.000001_real32 / runtime, &
        ' Mnodes/s'
    print '(A,F0.6,A)', 'throughput: ', real(edges, real32) * 0.000001_real32 / runtime, &
        ' Medges/s'
  end subroutine compute_mis

  subroutine initialize_status(nodes, nidx, nstatus, avg, scaledavg)
    integer(int32), intent(in) :: nodes
    integer(int32), intent(in) :: nidx(:)
    integer(int32), intent(inout) :: nstatus(:)
    real(real32), intent(in) :: avg, scaledavg

    integer(int32) :: i, degree, res
    integer(int64) :: h
    real(real32) :: x

    !$omp target teams distribute parallel do num_teams(blocks) thread_limit(threads_per_block) &
    !$omp& private(i,degree,res,h,x)
    do i = 1, nodes
      nstatus(i) = in_status
      degree = nidx(i + 1) - nidx(i)
      if (degree > 0) then
        h = modulo(int(i - 1, int64), 4294967296_int64)
        h = modulo(ieor(ishft(h, -16), h) * int(z'045D9F3B', int64), 4294967296_int64)
        h = modulo(ieor(ishft(h, -16), h) * int(z'045D9F3B', int64), 4294967296_int64)
        h = modulo(ieor(ishft(h, -16), h), 4294967296_int64)
        x = real(degree, real32) - real(h, real32) * 0.00000000023283064365386962890625_real32
        res = int(scaledavg / (avg + x), int32)
        nstatus(i) = ior(res + res, 1_int32)
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine initialize_status

  subroutine converge_status(nodes, nidx, nlist, nstatus)
    integer(int32), intent(in) :: nodes
    integer(int32), intent(in) :: nidx(:), nlist(:)
    integer(int32), intent(inout) :: nstatus(:)

    integer(int32) :: from, incr, missing, v, i, neighbor, nv
    logical :: blocked

    !$omp target teams num_teams(blocks) thread_limit(threads_per_block)
    !$omp parallel private(from,incr,missing,v,i,neighbor,nv,blocked)
      from = omp_get_thread_num() + omp_get_team_num() * threads_per_block + 1
      incr = omp_get_num_teams() * threads_per_block

      do
        missing = 0_int32
        do v = from, nodes, incr
          nv = nstatus(v)
          if (iand(nv, 1_int32) /= 0_int32) then
            blocked = .false.
            i = nidx(v) + 1
            do while (i <= nidx(v + 1))
              neighbor = nlist(i) + 1
              if (.not. ((nv > nstatus(neighbor)) .or. &
                  ((nv == nstatus(neighbor)) .and. ((v - 1) > nlist(i))))) then
                blocked = .true.
                exit
              end if
              i = i + 1
            end do
            if (blocked) then
              missing = 1_int32
            else
              do i = nidx(v) + 1, nidx(v + 1)
                nstatus(nlist(i) + 1) = out_status
              end do
              nstatus(v) = in_status
            end if
          end if
        end do
        if (missing == 0_int32) exit
      end do
    !$omp end parallel
    !$omp end target teams
  end subroutine converge_status

  subroutine verify_mis(nodes, nidx, nlist, nstatus)
    integer(int32), intent(in) :: nodes
    integer(int32), intent(in) :: nidx(:), nlist(:), nstatus(:)

    integer(int32) :: v, i, flag

    do v = 1, nodes
      if ((nstatus(v) /= in_status) .and. (nstatus(v) /= out_status)) then
        write(0, '(A)') 'ERROR: found unprocessed node in graph'
        print *
        exit
      end if
      if (nstatus(v) == in_status) then
        do i = nidx(v) + 1, nidx(v + 1)
          if (nstatus(nlist(i) + 1) == in_status) then
            write(0, '(A)') 'ERROR: found adjacent nodes in MIS'
            print *
            exit
          end if
        end do
      else
        flag = 0_int32
        do i = nidx(v) + 1, nidx(v + 1)
          if (nstatus(nlist(i) + 1) == in_status) flag = 1_int32
        end do
        if (flag == 0_int32) then
          write(0, '(A)') 'ERROR: set is not maximal'
          print *
          exit
        end if
      end if
    end do
  end subroutine verify_mis

end program mis
