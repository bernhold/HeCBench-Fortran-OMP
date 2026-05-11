program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: repeat = 1000
  integer, parameter :: n = 370
  real(real32), parameter :: alpha_re = 3.14_real32
  real(real32), parameter :: alpha_im = 1.59_real32
  real(real32), parameter :: beta_re = 2.71_real32
  real(real32), parameter :: beta_im = 8.28_real32

  real(real32) :: at_re(0:n-1, 0:n-1), at_im(0:n-1, 0:n-1)
  real(real32) :: x_re(0:n-1), x_im(0:n-1)
  real(real32) :: y0_re(0:n-1), y0_im(0:n-1)
  real(real32) :: y_cpu_re(0:n-1), y_cpu_im(0:n-1)
  real(real32) :: y_gpu_re(0:n-1), y_gpu_im(0:n-1)
  real(real32) :: tol_re, tol_im
  integer :: i, j
  real(real64) :: start_time, elapsed_us
  logical :: ok

  do i = 0, n - 1
    x_re(i) = real(i + 5, real32)
    x_im(i) = real(i * 2, real32)
    y0_re(i) = real(i * 3, real32)
    y0_im(i) = real(i + 7, real32)
    do j = 0, n - 1
      at_re(i, j) = real(i + j, real32)
      at_im(i, j) = real(i + 3, real32)
    end do
  end do

  call chemv_cpu(at_re, at_im, x_re, x_im, y0_re, y0_im, y_cpu_re, y_cpu_im)
  y_gpu_re = y0_re
  y_gpu_im = y0_im

  !$omp target data map(to: at_re, at_im, x_re, x_im, y0_re, y0_im) map(from: y_gpu_re, y_gpu_im)
  start_time = omp_get_wtime()
  do i = 1, repeat
    call chemv_kernel(at_re, at_im, x_re, x_im, y0_re, y0_im, y_gpu_re, y_gpu_im)
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'(A,F0.6,A)') 'Average execution time of chemv kernels: ', elapsed_us, ' (us)'
  !$omp end target data

  ok = .true.
  do i = 0, n - 1
    tol_re = max(1.0e-3_real32, 1.0e-4_real32 * abs(y_cpu_re(i)))
    tol_im = max(1.0e-3_real32, 1.0e-4_real32 * abs(y_cpu_im(i)))
    if ((abs(y_cpu_re(i) - y_gpu_re(i)) > tol_re) .or. &
        (abs(y_cpu_im(i) - y_gpu_im(i)) > tol_im)) then
      ok = .false.
      write(*,'(I0,1X,F0.6,1X,F0.6)') i, y_cpu_re(i), y_gpu_re(i)
      exit
    end if
  end do

  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
    stop 1
  end if

contains

  subroutine chemv_cpu(at_re, at_im, x_re, x_im, y_in_re, y_in_im, y_re, y_im)
    real(real32), intent(in) :: at_re(0:n-1, 0:n-1), at_im(0:n-1, 0:n-1)
    real(real32), intent(in) :: x_re(0:n-1), x_im(0:n-1)
    real(real32), intent(in) :: y_in_re(0:n-1), y_in_im(0:n-1)
    real(real32), intent(out) :: y_re(0:n-1), y_im(0:n-1)
    integer :: row

    do row = 0, n - 1
      call chemv_row(row, at_re, at_im, x_re, x_im, y_in_re, y_in_im, y_re(row), y_im(row))
    end do
  end subroutine chemv_cpu

  subroutine chemv_kernel(at_re, at_im, x_re, x_im, y_in_re, y_in_im, y_re, y_im)
    real(real32), intent(in) :: at_re(0:n-1, 0:n-1), at_im(0:n-1, 0:n-1)
    real(real32), intent(in) :: x_re(0:n-1), x_im(0:n-1)
    real(real32), intent(in) :: y_in_re(0:n-1), y_in_im(0:n-1)
    real(real32), intent(inout) :: y_re(0:n-1), y_im(0:n-1)
    integer :: row

    !$omp target teams distribute parallel do thread_limit(256)
    do row = 0, n - 1
      call chemv_row(row, at_re, at_im, x_re, x_im, y_in_re, y_in_im, y_re(row), y_im(row))
    end do
    !$omp end target teams distribute parallel do
  end subroutine chemv_kernel

  subroutine chemv_row(row, at_re, at_im, x_re, x_im, y_in_re, y_in_im, out_re, out_im)
    integer, intent(in) :: row
    real(real32), intent(in) :: at_re(0:n-1, 0:n-1), at_im(0:n-1, 0:n-1)
    real(real32), intent(in) :: x_re(0:n-1), x_im(0:n-1)
    real(real32), intent(in) :: y_in_re(0:n-1), y_in_im(0:n-1)
    real(real32), intent(out) :: out_re, out_im
    integer :: col
    real(real32) :: ar, ai, mr, mi

    out_re = y_in_re(row) * beta_re - y_in_im(row) * beta_im
    out_im = y_in_im(row) * beta_re + y_in_re(row) * beta_im

    ar = alpha_re * at_re(row, row)
    ai = alpha_im * at_re(row, row)
    out_re = out_re + ar * x_re(row) - ai * x_im(row)
    out_im = out_im + ai * x_re(row) + ar * x_im(row)

    do col = row + 1, n - 1
      ar = alpha_re * at_re(row, col) + alpha_im * at_im(row, col)
      ai = alpha_im * at_re(row, col) - alpha_re * at_im(row, col)
      mr = ar * x_re(col) - ai * x_im(col)
      mi = ai * x_re(col) + ar * x_im(col)
      out_re = out_re + mr
      out_im = out_im + mi
    end do

    do col = 0, row - 1
      ar = alpha_re * at_re(col, row) - alpha_im * at_im(col, row)
      ai = alpha_im * at_re(col, row) + alpha_re * at_im(col, row)
      mr = ar * x_re(col) - ai * x_im(col)
      mi = ai * x_re(col) + ar * x_im(col)
      out_re = out_re + mr
      out_im = out_im + mi
    end do
  end subroutine chemv_row

end program main
