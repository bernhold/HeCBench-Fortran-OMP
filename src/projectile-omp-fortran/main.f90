program main
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: num_elements = 10000000
  integer, parameter :: block_size = 256
  real(real32), parameter :: k_pi_value = 3.1415_real32
  real(real32), parameter :: k_g_value = 9.81_real32

  character(len=256) :: arg0, arg
  integer :: repeat, i, iter, errors
  integer(int32) :: seed
  real(real32), allocatable :: angle(:), velocity(:)
  real(real32), allocatable :: range_out(:), time_out(:), height_out(:)
  real(real32), allocatable :: range_ref(:), time_ref(:), height_ref(:)
  real(real64) :: start_time, elapsed

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) repeat
  if (repeat <= 0) stop 1

  allocate(angle(num_elements), velocity(num_elements))
  allocate(range_out(num_elements), time_out(num_elements), height_out(num_elements))
  allocate(range_ref(num_elements), time_ref(num_elements), height_ref(num_elements))

  seed = 2_int32
  do i = 1, num_elements
    angle(i) = real(mod(c_rand(seed), 90_int32) + 10_int32, real32)
    velocity(i) = real(mod(c_rand(seed), 400_int32) + 10_int32, real32)
  end do

  range_out = 0.0_real32
  time_out = 0.0_real32
  height_out = 0.0_real32

  !$omp target data map(to: angle(1:num_elements), velocity(1:num_elements)) &
  !$omp& map(from: range_out(1:num_elements), time_out(1:num_elements), height_out(1:num_elements))
  start_time = omp_get_wtime()
  do iter = 1, repeat
    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 1, num_elements
      call compute_projectile(angle(i), velocity(i), range_out(i), time_out(i), height_out(i))
    end do
    !$omp end target teams distribute parallel do
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', elapsed / real(repeat, real64), ' (s)'

  do i = 1, num_elements
    call compute_projectile(angle(i), velocity(i), range_ref(i), time_ref(i), height_ref(i))
  end do

  errors = 0
  do i = 1, num_elements
    if (abs(angle(i) - angle(i)) > 1.0_real32 .or. &
        abs(velocity(i) - velocity(i)) > 1.0_real32 .or. &
        abs(range_out(i) - range_ref(i)) > 1.0_real32 .or. &
        abs(time_out(i) - time_ref(i)) > 1.0_real32 .or. &
        abs(height_out(i) - height_ref(i)) > 1.0_real32) then
      errors = errors + 1
      exit
    end if
  end do

  if (errors == 0) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(angle, velocity, range_out, time_out, height_out, range_ref, time_ref, height_ref)

contains

  integer(int32) function c_rand(seed)
    integer(int32), intent(inout) :: seed
    integer, parameter :: int64 = selected_int_kind(18)
    integer(int64) :: next_value

    next_value = mod(1103515245_int64 * int(seed, int64) + 12345_int64, 2147483648_int64)
    seed = int(next_value, int32)
    c_rand = iand(seed / 65536_int32, 32767_int32)
  end function c_rand

  subroutine compute_projectile(proj_angle, proj_vel, max_range, total_time, max_height)
    real(real32), intent(in) :: proj_angle, proj_vel
    real(real32), intent(out) :: max_range, total_time, max_height
    real(real32) :: sin_value, cos_value

    sin_value = sin(proj_angle * k_pi_value / 180.0_real32)
    cos_value = cos(proj_angle * k_pi_value / 180.0_real32)
    total_time = abs(2.0_real32 * proj_vel * sin_value) / k_g_value
    max_range = abs(proj_vel * total_time * cos_value)
    max_height = (proj_vel * proj_vel * sin_value * sin_value) / 2.0_real32 * k_g_value
  end subroutine compute_projectile

end program main
