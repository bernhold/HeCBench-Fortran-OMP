program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  type :: particle
    real(real32) :: pos(3)
    real(real32) :: vel(3)
    real(real32) :: acc(3)
    real(real32) :: mass
  end type particle

  integer(int32) :: npart, nsteps
  character(len=256) :: arg

  npart = 16000_int32
  nsteps = 10_int32

  if (command_argument_count() > 0) then
    call get_command_argument(1, arg)
    read(arg, *) npart
    if (command_argument_count() == 2) then
      call get_command_argument(2, arg)
      read(arg, *) nsteps
      if (nsteps < 3_int32) then
        write(*,'(A)') 'The number of integration steps should be at least 3'
        stop 1
      end if
    end if
  end if

  call run_simulation(npart, nsteps)

contains

  subroutine run_simulation(npart, nsteps)
    integer(int32), intent(in) :: npart, nsteps
    type(particle), allocatable :: particles(:), ref_particles(:)
    real(real32), allocatable :: energy(:), ref_energy(:)
    real(real32) :: kenergy, ref_kenergy
    real(real64) :: total_time, total_flops
    logical :: ok

    write(*,'(A)') '==============================='
    write(*,'(A)') ' Initialize Gravity Simulation'

    allocate(particles(npart), ref_particles(npart), energy(npart), ref_energy(npart))
    call initialize_particles(particles)
    ref_particles = particles
    energy = 0.0_real32
    ref_energy = 0.0_real32

    call start_device(particles, energy, npart, nsteps, kenergy, total_time, total_flops)
    call start_reference(ref_particles, ref_energy, npart, nsteps, ref_kenergy)

    ok = abs(kenergy - ref_kenergy) < 1.0e-3_real32
    write(*,*)
    write(*,'(A)') merge('PASS', 'FAIL', ok)

    deallocate(ref_energy, energy, ref_particles, particles)
  end subroutine run_simulation

  subroutine initialize_particles(particles)
    type(particle), intent(out) :: particles(:)
    integer(int32) :: i
    integer(int64) :: state_pos, state_vel, state_mass

    state_pos = 42_int64
    state_vel = 42_int64
    state_mass = 42_int64

    do i = 1, size(particles)
      particles(i)%pos(1) = rand_uniform(state_pos, 0.0_real32, 1.0_real32)
      particles(i)%pos(2) = rand_uniform(state_pos, 0.0_real32, 1.0_real32)
      particles(i)%pos(3) = rand_uniform(state_pos, 0.0_real32, 1.0_real32)
    end do

    do i = 1, size(particles)
      particles(i)%vel(1) = rand_uniform(state_vel, -1.0_real32, 1.0_real32) * 1.0e-3_real32
      particles(i)%vel(2) = rand_uniform(state_vel, -1.0_real32, 1.0_real32) * 1.0e-3_real32
      particles(i)%vel(3) = rand_uniform(state_vel, -1.0_real32, 1.0_real32) * 1.0e-3_real32
    end do

    do i = 1, size(particles)
      particles(i)%acc = 0.0_real32
    end do

    do i = 1, size(particles)
      particles(i)%mass = real(size(particles), real32) * rand_uniform(state_mass, 0.0_real32, 1.0_real32)
    end do
  end subroutine initialize_particles

  real(real32) function rand_uniform(state, low, high) result(value)
    integer(int64), intent(inout) :: state
    real(real32), intent(in) :: low, high
    integer(int64), parameter :: a = 6364136223846793005_int64
    integer(int64), parameter :: c = 1442695040888963407_int64
    integer(int64) :: bits
    real(real32) :: unit

    state = a * state + c
    bits = iand(ishft(state, -40), int(z'FFFFFF', int64))
    unit = real(bits, real32) / real(int(z'1000000', int64), real32)
    value = low + (high - low) * unit
  end function rand_uniform

  subroutine start_device(particles, energy, n, nsteps, kenergy, total_time, total_flops)
    type(particle), intent(inout) :: particles(:)
    real(real32), intent(inout) :: energy(:)
    integer(int32), intent(in) :: n, nsteps
    real(real32), intent(out) :: kenergy
    real(real64), intent(out) :: total_time, total_flops
    real(real32), parameter :: dt = 0.1_real32
    real(real32), parameter :: softening_squared = 1.0e-3_real32
    real(real32), parameter :: grav_const = 6.67259e-11_real32
    real(real64) :: gflops, av, dev, elapsed_seconds, ts0, t0
    integer(int32) :: s, nf

    gflops = 1.0e-9_real64 * ((11.0_real64 + 18.0_real64) * real(n, real64) * real(n, real64) + &
      real(n, real64) * 19.0_real64)
    nf = 0_int32
    av = 0.0_real64
    dev = 0.0_real64
    kenergy = 0.0_real32
    total_time = 0.0_real64
    total_flops = 0.0_real64

    t0 = omp_get_wtime()
    !$omp target data map(to: particles(1:n)) map(alloc: energy(1:n))
    do s = 1, nsteps
      ts0 = omp_get_wtime()
      call accelerate_particles_device(particles, n, softening_squared, grav_const)
      call update_particles_device(particles, energy, n, dt)
      call accumulate_energy_device(energy, n)
      elapsed_seconds = omp_get_wtime() - ts0

      !$omp target update from(energy(1:1))
      kenergy = 0.5_real32 * energy(1)
      energy(1) = 0.0_real32

      nf = nf + 1_int32
      if (nf > 2_int32) then
        av = av + gflops / elapsed_seconds
        dev = dev + gflops * gflops / (elapsed_seconds * elapsed_seconds)
      end if
    end do
    !$omp end target data

    total_time = omp_get_wtime() - t0
    total_flops = gflops * real(nsteps, real64)
    av = av / real(nf - 2_int32, real64)
    if (nf == 3_int32) then
      dev = 0.0_real64
    else
      dev = sqrt(dev / real(nf - 2_int32, real64) - av * av)
    end if

    write(*,*)
    write(*,'(A,1X,ES16.8)') '# Total Energy        :', kenergy
    write(*,'(A,1X,ES16.8)') '# Total Time (s)      :', total_time
    write(*,'(A,1X,ES16.8,A,1X,ES16.8)') '# Average Performance :', av, ' +-', dev
    write(*,'(A)') '==============================='
  end subroutine start_device

  subroutine accelerate_particles_device(particles, n, softening_squared, grav_const)
    type(particle), intent(inout) :: particles(:)
    integer(int32), intent(in) :: n
    real(real32), intent(in) :: softening_squared, grav_const
    integer(int32) :: i, j
    type(particle) :: pi, pj
    real(real32) :: acc0, acc1, acc2, dx, dy, dz, distance_sqr, distance_inv, scale

    !$omp target teams distribute parallel do thread_limit(256) private(pi,pj,acc0,acc1,acc2,dx,dy,dz,distance_sqr,distance_inv,scale)
    do i = 1, n
      pi = particles(i)
      acc0 = pi%acc(1)
      acc1 = pi%acc(2)
      acc2 = pi%acc(3)
      do j = 1, n
        pj = particles(j)
        dx = pj%pos(1) - pi%pos(1)
        dy = pj%pos(2) - pi%pos(2)
        dz = pj%pos(3) - pi%pos(3)
        distance_sqr = dx * dx + dy * dy + dz * dz + softening_squared
        distance_inv = 1.0_real32 / sqrt(distance_sqr)
        scale = grav_const * pj%mass * distance_inv * distance_inv * distance_inv
        acc0 = acc0 + dx * scale
        acc1 = acc1 + dy * scale
        acc2 = acc2 + dz * scale
      end do
      pi%acc(1) = acc0
      pi%acc(2) = acc1
      pi%acc(3) = acc2
      particles(i) = pi
    end do
    !$omp end target teams distribute parallel do
  end subroutine accelerate_particles_device

  subroutine update_particles_device(particles, energy, n, dt)
    type(particle), intent(inout) :: particles(:)
    real(real32), intent(inout) :: energy(:)
    integer(int32), intent(in) :: n
    real(real32), intent(in) :: dt
    integer(int32) :: i
    type(particle) :: pi

    !$omp target teams distribute parallel do thread_limit(256) private(pi)
    do i = 1, n
      pi = particles(i)
      pi%vel(1) = pi%vel(1) + pi%acc(1) * dt
      pi%vel(2) = pi%vel(2) + pi%acc(2) * dt
      pi%vel(3) = pi%vel(3) + pi%acc(3) * dt
      pi%pos(1) = pi%pos(1) + pi%vel(1) * dt
      pi%pos(2) = pi%pos(2) + pi%vel(2) * dt
      pi%pos(3) = pi%pos(3) + pi%vel(3) * dt
      pi%acc = 0.0_real32
      energy(i) = pi%mass * (pi%vel(1) * pi%vel(1) + pi%vel(2) * pi%vel(2) + pi%vel(3) * pi%vel(3))
      particles(i) = pi
    end do
    !$omp end target teams distribute parallel do
  end subroutine update_particles_device

  subroutine accumulate_energy_device(energy, n)
    real(real32), intent(inout) :: energy(:)
    integer(int32), intent(in) :: n
    integer(int32) :: i

    !$omp target
    do i = 2, n
      energy(1) = energy(1) + energy(i)
    end do
    !$omp end target
  end subroutine accumulate_energy_device

  subroutine start_reference(particles, energy, n, nsteps, ref_kenergy)
    type(particle), intent(inout) :: particles(:)
    real(real32), intent(inout) :: energy(:)
    integer(int32), intent(in) :: n, nsteps
    real(real32), intent(out) :: ref_kenergy
    real(real32), parameter :: dt = 0.1_real32
    real(real32), parameter :: softening_squared = 1.0e-3_real32
    real(real32), parameter :: grav_const = 6.67259e-11_real32
    integer(int32) :: s

    ref_kenergy = 0.0_real32
    do s = 1, nsteps
      call accelerate_particles_host(particles, n, softening_squared, grav_const)
      call update_particles_host(particles, energy, n, dt)
      call accumulate_energy_host(energy, n)
      ref_kenergy = 0.5_real32 * energy(1)
      energy(1) = 0.0_real32
    end do
  end subroutine start_reference

  subroutine accelerate_particles_host(particles, n, softening_squared, grav_const)
    type(particle), intent(inout) :: particles(:)
    integer(int32), intent(in) :: n
    real(real32), intent(in) :: softening_squared, grav_const
    integer(int32) :: i, j
    type(particle) :: pi, pj
    real(real32) :: acc0, acc1, acc2, dx, dy, dz, distance_sqr, distance_inv, scale

    do i = 1, n
      pi = particles(i)
      acc0 = pi%acc(1)
      acc1 = pi%acc(2)
      acc2 = pi%acc(3)
      do j = 1, n
        pj = particles(j)
        dx = pj%pos(1) - pi%pos(1)
        dy = pj%pos(2) - pi%pos(2)
        dz = pj%pos(3) - pi%pos(3)
        distance_sqr = dx * dx + dy * dy + dz * dz + softening_squared
        distance_inv = 1.0_real32 / sqrt(distance_sqr)
        scale = grav_const * pj%mass * distance_inv * distance_inv * distance_inv
        acc0 = acc0 + dx * scale
        acc1 = acc1 + dy * scale
        acc2 = acc2 + dz * scale
      end do
      pi%acc(1) = acc0
      pi%acc(2) = acc1
      pi%acc(3) = acc2
      particles(i) = pi
    end do
  end subroutine accelerate_particles_host

  subroutine update_particles_host(particles, energy, n, dt)
    type(particle), intent(inout) :: particles(:)
    real(real32), intent(inout) :: energy(:)
    integer(int32), intent(in) :: n
    real(real32), intent(in) :: dt
    integer(int32) :: i
    type(particle) :: pi

    do i = 1, n
      pi = particles(i)
      pi%vel(1) = pi%vel(1) + pi%acc(1) * dt
      pi%vel(2) = pi%vel(2) + pi%acc(2) * dt
      pi%vel(3) = pi%vel(3) + pi%acc(3) * dt
      pi%pos(1) = pi%pos(1) + pi%vel(1) * dt
      pi%pos(2) = pi%pos(2) + pi%vel(2) * dt
      pi%pos(3) = pi%pos(3) + pi%vel(3) * dt
      pi%acc = 0.0_real32
      energy(i) = pi%mass * (pi%vel(1) * pi%vel(1) + pi%vel(2) * pi%vel(2) + pi%vel(3) * pi%vel(3))
      particles(i) = pi
    end do
  end subroutine update_particles_host

  subroutine accumulate_energy_host(energy, n)
    real(real32), intent(inout) :: energy(:)
    integer(int32), intent(in) :: n
    integer(int32) :: i

    do i = 2, n
      energy(1) = energy(1) + energy(i)
    end do
  end subroutine accumulate_energy_host

end program main
