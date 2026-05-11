program main
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2, arg3, arg4
  integer :: dim_x, dim_y, dim_z, repeat, size, i
  real(real64), allocatable :: ksat(:), psi(:), c(:), theta(:), k(:)
  real(real64), allocatable :: c_ref(:), theta_ref(:), k_ref(:)
  real(real64) :: start_time, end_time, avg_time
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    print '(3A)', 'Usage: ./', trim(arg0), ' <dimX> <dimY> <dimZ> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  call get_command_argument(4, arg4)
  read(arg1, *) dim_x
  read(arg2, *) dim_y
  read(arg3, *) dim_z
  read(arg4, *) repeat

  if (dim_x <= 0 .or. dim_y <= 0 .or. dim_z <= 0 .or. repeat <= 0) stop 1
  size = dim_x * dim_y * dim_z

  allocate(ksat(size), psi(size), c(size), theta(size), k(size))
  allocate(c_ref(size), theta_ref(size), k_ref(size))

  do i = 1, size
    ksat(i) = 1.0e-6_real64 + (1.0_real64 - 1.0e-6_real64) * real(i - 1, real64) / real(size, real64)
    psi(i) = -100.0_real64 + 101.0_real64 * real(i - 1, real64) / real(size, real64)
  end do

  call reference(ksat, psi, c_ref, theta_ref, k_ref, size)

  !$omp target data map(to: ksat(1:size), psi(1:size)) map(from: c(1:size), theta(1:size), k(1:size))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call van_genuchten(ksat, psi, c, theta, k, size)
  end do
  end_time = omp_get_wtime()
  avg_time = (end_time - start_time) / real(repeat, real64)
  print '(A,F8.6,A)', 'Average kernel execution time: ', avg_time, ' (s)'
  !$omp end target data

  ok = .true.
  do i = 1, size
    if (abs(c(i) - c_ref(i)) > 1.0e-3_real64 .or. &
        abs(theta(i) - theta_ref(i)) > 1.0e-3_real64 .or. &
        abs(k(i) - k_ref(i)) > 1.0e-3_real64) then
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(ksat, psi, c, theta, k, c_ref, theta_ref, k_ref)

contains

  subroutine van_genuchten(ksat, psi, c, theta, k, size)
    real(real64), intent(in) :: ksat(:), psi(:)
    real(real64), intent(out) :: c(:), theta(:), k(:)
    integer, intent(in) :: size
    integer :: idx
    real(real64) :: se, theta_value, psi_value, lambda, m, t
    real(real64), parameter :: alpha = 0.02_real64
    real(real64), parameter :: theta_s = 0.45_real64
    real(real64), parameter :: theta_r = 0.1_real64
    real(real64), parameter :: n_value = 1.8_real64

    lambda = n_value - 1.0_real64
    m = lambda / n_value

    !$omp target teams distribute parallel do thread_limit(256) private(se, theta_value, psi_value, t)
    do idx = 1, size
      psi_value = psi(idx) * 100.0_real64
      if (psi_value < 0.0_real64) then
        theta_value = (theta_s - theta_r) / &
            (1.0_real64 + (alpha * (-psi_value)) ** n_value) ** m + theta_r
      else
        theta_value = theta_s
      end if

      theta(idx) = theta_value
      se = (theta_value - theta_r) / (theta_s - theta_r)
      t = 1.0_real64 - (1.0_real64 - se ** (1.0_real64 / m)) ** m
      k(idx) = ksat(idx) * sqrt(se) * t * t

      if (psi_value < 0.0_real64) then
        c(idx) = 100.0_real64 * alpha * n_value * (1.0_real64 / n_value - 1.0_real64) * &
            (alpha * abs(psi_value)) ** (n_value - 1.0_real64) * (theta_r - theta_s) * &
            ((alpha * abs(psi_value)) ** n_value + 1.0_real64) ** (1.0_real64 / n_value - 2.0_real64)
      else
        c(idx) = 0.0_real64
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine van_genuchten

  subroutine reference(ksat, psi, c, theta, k, size)
    real(real64), intent(in) :: ksat(:), psi(:)
    real(real64), intent(out) :: c(:), theta(:), k(:)
    integer, intent(in) :: size
    integer :: idx
    real(real64) :: se, theta_value, psi_value, lambda, m, t
    real(real64), parameter :: alpha = 0.02_real64
    real(real64), parameter :: theta_s = 0.45_real64
    real(real64), parameter :: theta_r = 0.1_real64
    real(real64), parameter :: n_value = 1.8_real64

    lambda = n_value - 1.0_real64
    m = lambda / n_value

    do idx = 1, size
      psi_value = psi(idx) * 100.0_real64
      if (psi_value < 0.0_real64) then
        theta_value = (theta_s - theta_r) / &
            (1.0_real64 + (alpha * (-psi_value)) ** n_value) ** m + theta_r
      else
        theta_value = theta_s
      end if

      theta(idx) = theta_value
      se = (theta_value - theta_r) / (theta_s - theta_r)
      t = 1.0_real64 - (1.0_real64 - se ** (1.0_real64 / m)) ** m
      k(idx) = ksat(idx) * sqrt(se) * t * t

      if (psi_value < 0.0_real64) then
        c(idx) = 100.0_real64 * alpha * n_value * (1.0_real64 / n_value - 1.0_real64) * &
            (alpha * abs(psi_value)) ** (n_value - 1.0_real64) * (theta_r - theta_s) * &
            ((alpha * abs(psi_value)) ** n_value + 1.0_real64) ** (1.0_real64 / n_value - 2.0_real64)
      else
        c(idx) = 0.0_real64
      end if
    end do
  end subroutine reference

end program main
