program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

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

  integer, parameter :: x_dim = 11
  integer, parameter :: y_dim = 40
  integer, parameter :: batch = 4
  integer, parameter :: states = 3
  integer, parameter :: nstate = states - 1

  character(len=256) :: arg0, arg
  integer :: repeat, count, i, j, ii, jj, b, s
  real(real64), allocatable :: cur(:,:,:,:), next(:,:,:,:), cpu_cur(:,:,:,:), cpu_next(:,:,:,:)
  real(real64), allocatable :: emis(:,:,:,:), trans(:,:,:,:), like(:,:,:,:), start(:,:)
  real(real64) :: start_time, elapsed_ms, checksum
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if
  call get_command_argument(1, arg)
  read(arg, *) repeat

  allocate(cur(0:x_dim,0:y_dim,0:batch-1,0:nstate-1))
  allocate(next(0:x_dim,0:y_dim,0:batch-1,0:nstate-1))
  allocate(cpu_cur(0:x_dim,0:y_dim,0:batch-1,0:nstate-1))
  allocate(cpu_next(0:x_dim,0:y_dim,0:batch-1,0:nstate-1))
  allocate(emis(0:x_dim,0:y_dim,0:batch-1,0:nstate-1))
  allocate(trans(0:x_dim,0:batch-1,0:nstate-1,0:states-1))
  allocate(like(0:1,0:1,0:batch-1,0:nstate-1))
  allocate(start(0:batch-1,0:nstate-1))

  call initialize(cur, emis, trans, like, start)
  next = 0.0_real64
  cpu_cur = cur
  cpu_next = 0.0_real64

  do count = 1, repeat
    do i = 1, x_dim
      do j = 1, y_dim
        call pair_hmm_forward_host(i, j, cpu_cur, trans, emis, like, start, cpu_next)
        cpu_cur = cpu_next
      end do
    end do
  end do

  !$omp target data map(to: emis, trans, like, start) map(tofrom: cur, next)
  start_time = omp_get_wtime()
  do count = 1, repeat
    do i = 1, x_dim
      do j = 1, y_dim
        call pair_hmm_forward_device(i, j, cur, trans, emis, like, start, next)
        !$omp target teams distribute parallel do collapse(4) private(ii, jj, b, s)
        do ii = 0, x_dim
          do jj = 0, y_dim
            do b = 0, batch - 1
              do s = 0, nstate - 1
                cur(ii,jj,b,s) = next(ii,jj,b,s)
              end do
            end do
          end do
        end do
        !$omp end target teams distribute parallel do
      end do
    end do
  end do
  elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64
  !$omp end target data

  ok = all(abs(cpu_cur - cur) <= 1.0e-9_real64)
  if (.not. ok) stop 1

  checksum = sum(cur)
  write(*,'(A,F0.6,A)') 'Total execution time ', elapsed_ms, ' milliseconds'
  write(*,'(A,F0.6)') 'Checksum ', checksum

  deallocate(cur, next, cpu_cur, cpu_next, emis, trans, like, start)

contains

  subroutine initialize(cur, emis, trans, like, start)
    real(real64), intent(out) :: cur(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(out) :: emis(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(out) :: trans(0:x_dim,0:batch-1,0:nstate-1,0:states-1)
    real(real64), intent(out) :: like(0:1,0:1,0:batch-1,0:nstate-1)
    real(real64), intent(out) :: start(0:batch-1,0:nstate-1)
    integer :: i, j, b, s, t, a, c

    call c_srand(123_c_int)
    do i = 0, x_dim
      do j = 0, y_dim
        do b = 0, batch - 1
          do s = 0, nstate - 1
            cur(i,j,b,s) = next_random()
            emis(i,j,b,s) = next_random()
          end do
        end do
      end do
    end do
    do i = 0, x_dim
      do b = 0, batch - 1
        do s = 0, nstate - 1
          do t = 0, states - 1
            trans(i,b,s,t) = next_random()
          end do
        end do
      end do
    end do
    do b = 0, batch - 1
      do s = 0, nstate - 1
        start(b,s) = next_random()
      end do
    end do
    do a = 0, 1
      do c = 0, 1
        do b = 0, batch - 1
          do s = 0, nstate - 1
            like(a,c,b,s) = next_random()
          end do
        end do
      end do
    end do
  end subroutine initialize

  real(real64) function next_random()
    next_random = real(c_rand(), real64) / 2147483647.0_real64
  end function next_random

  subroutine pair_hmm_forward_host(cur_i, cur_j, forward_in, transitions, emissions, likelihood, start_transitions, forward_out)
    integer, intent(in) :: cur_i, cur_j
    real(real64), intent(in) :: forward_in(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: transitions(0:x_dim,0:batch-1,0:nstate-1,0:states-1)
    real(real64), intent(in) :: emissions(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: likelihood(0:1,0:1,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: start_transitions(0:batch-1,0:nstate-1)
    real(real64), intent(inout) :: forward_out(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    integer :: b, s

    do b = 0, batch - 1
      do s = 0, nstate - 1
        forward_out(cur_i,cur_j,b,s) = pair_hmm_value(cur_i, cur_j, b, s, forward_in, transitions, emissions, likelihood, start_transitions)
      end do
    end do
  end subroutine pair_hmm_forward_host

  subroutine pair_hmm_forward_device(cur_i, cur_j, forward_in, transitions, emissions, likelihood, start_transitions, forward_out)
    integer, intent(in) :: cur_i, cur_j
    real(real64), intent(in) :: forward_in(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: transitions(0:x_dim,0:batch-1,0:nstate-1,0:states-1)
    real(real64), intent(in) :: emissions(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: likelihood(0:1,0:1,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: start_transitions(0:batch-1,0:nstate-1)
    real(real64), intent(inout) :: forward_out(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    integer :: b, s

    !$omp target teams distribute parallel do collapse(2) num_teams(batch) thread_limit(nstate)
    do b = 0, batch - 1
      do s = 0, nstate - 1
        forward_out(cur_i,cur_j,b,s) = pair_hmm_value(cur_i, cur_j, b, s, forward_in, transitions, emissions, likelihood, start_transitions)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine pair_hmm_forward_device

  real(real64) function pair_hmm_value(cur_i, cur_j, b, s, forward_in, transitions, emissions, likelihood, start_transitions)
    integer, intent(in) :: cur_i, cur_j, b, s
    real(real64), intent(in) :: forward_in(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: transitions(0:x_dim,0:batch-1,0:nstate-1,0:states-1)
    real(real64), intent(in) :: emissions(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: likelihood(0:1,0:1,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: start_transitions(0:batch-1,0:nstate-1)
    integer :: k
    real(real64) :: total

    if (cur_i == 1 .and. cur_j == 0) then
      pair_hmm_value = start_transitions(b,s) * emissions(cur_i,cur_j,b,s)
    else if (cur_i > 0 .and. cur_j == 0) then
      total = 0.0_real64
      do k = 0, nstate - 1
        total = total + forward_in(cur_i - 1,cur_j,b,k) * transitions(cur_i - 1,b,k,s)
      end do
      pair_hmm_value = total * emissions(cur_i,cur_j,b,s) * likelihood(0,1,b,s)
    else
      total = 0.0_real64
      do k = 0, nstate - 1
        total = total + forward_in(cur_i - 1,cur_j - 1,b,k) * transitions(cur_i - 1,b,k,s) * likelihood(0,0,b,s)
        total = total + forward_in(cur_i - 1,cur_j,b,k) * transitions(cur_i - 1,b,k,s) * likelihood(0,1,b,s)
        total = total + forward_in(cur_i,cur_j - 1,b,k) * transitions(cur_i,b,k,s) * likelihood(1,0,b,s)
        total = total + forward_in(cur_i,cur_j,b,k) * transitions(cur_i,b,k,s) * likelihood(1,1,b,s)
      end do
      pair_hmm_value = total * emissions(cur_i,cur_j,b,s)
    end if
  end function pair_hmm_value

end program main
