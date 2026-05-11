  subroutine pair32(pr, pi, a, b, idx, peer)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: idx, peer
    real(real32) :: tmp_real, tmp_imag
    tmp_real = pr(idx)
    tmp_imag = pi(idx)
    pr(idx) = a * tmp_real - b * pi(peer)
    pi(idx) = a * tmp_imag + b * pr(peer)
    pr(peer) = a * pr(peer) - b * tmp_imag
    pi(peer) = a * pi(peer) + b * tmp_real
  end subroutine pair32

  subroutine pair64(pr, pi, a, b, idx, peer)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: idx, peer
    real(real64) :: tmp_real, tmp_imag
    tmp_real = pr(idx)
    tmp_imag = pi(idx)
    pr(idx) = a * tmp_real - b * pi(peer)
    pi(idx) = a * tmp_imag + b * pr(peer)
    pr(peer) = a * pr(peer) - b * tmp_imag
    pi(peer) = a * pi(peer) + b * tmp_real
  end subroutine pair64

  subroutine kernel1_32(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    do x = 0, width - 1, 2
      call pair32(pr, pi, a, b, x, width + x)
    end do
    do y = 1, height - 2
      do idx = y * width + mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel1_32

  subroutine kernel2_32(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    do y = 0, height - 1
      do idx = y * width + mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel2_32

  subroutine kernel3_32(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    do x = 1, width - 1, 2
      call pair32(pr, pi, a, b, x, width + x)
    end do
    do y = 1, height - 2
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel3_32

  subroutine kernel4_32(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    do y = 0, height - 1
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel4_32

  subroutine kernel1_64(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    do x = 0, width - 1, 2
      call pair64(pr, pi, a, b, x, width + x)
    end do
    do y = 1, height - 2
      do idx = y * width + mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel1_64

  subroutine kernel2_64(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    do y = 0, height - 1
      do idx = y * width + mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel2_64

  subroutine kernel3_64(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    do x = 1, width - 1, 2
      call pair64(pr, pi, a, b, x, width + x)
    end do
    do y = 1, height - 2
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel3_64

  subroutine kernel4_64(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    do y = 0, height - 1
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel4_64

  subroutine kernel1_32_device(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    !$omp target teams distribute parallel do thread_limit(256)
    do x = 0, width - 1, 2
      call pair32(pr, pi, a, b, x, width + x)
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 1, height - 2
      do idx = y * width + mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel1_32_device

  subroutine kernel2_32_device(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 0, height - 1
      do idx = y * width + mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel2_32_device

  subroutine kernel3_32_device(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    !$omp target teams distribute parallel do thread_limit(256)
    do x = 1, width - 1, 2
      call pair32(pr, pi, a, b, x, width + x)
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 1, height - 2
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel3_32_device

  subroutine kernel4_32_device(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 0, height - 1
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel4_32_device

  subroutine kernel1_64_device(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    !$omp target teams distribute parallel do thread_limit(256)
    do x = 0, width - 1, 2
      call pair64(pr, pi, a, b, x, width + x)
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 1, height - 2
      do idx = y * width + mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel1_64_device

  subroutine kernel2_64_device(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 0, height - 1
      do idx = y * width + mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel2_64_device

  subroutine kernel3_64_device(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    !$omp target teams distribute parallel do thread_limit(256)
    do x = 1, width - 1, 2
      call pair64(pr, pi, a, b, x, width + x)
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 1, height - 2
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel3_64_device

  subroutine kernel4_64_device(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 0, height - 1
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel4_64_device
