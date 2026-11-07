PROGRAM mixc

  USE prms
  USE data
  USE derivatives
  USE rungekutta
  USE initialcond
  USE out
  USE mesh
  USE bc
  USE filter
  USE bspline
  USE col_linalg
  USE test_collocation

  IMPLICIT none

  call init
!  call solve
CONTAINS

  SUBROUTINE init

    Pi = 4.*ATAN(1.)


    call initprms; print *, "INITALIZED PRMS"
    call initdata; print *, "INITALIZED DATA"
    call initbspline; print *, "INITALIZED BSPLINES"
    call initmesh; print *, "INITALIZED MESH"
!    call initders; print *, "INITALIZED DERS"

!    call test_GCV_vector
    call restart_vector_GCV_channel_solve(103620)
!    call mms_solve
!    call test_ns_rhs
!    call initRK
!    call initFilters; print *, "INITIALIZED FILTERS"
!    call ic; print *, "INITIAL CONDITIONS"
!    call initbc(Q); print *, "INITIALIZED BCS"
!    call output(Q,0); print *, "INITIAL OUTPUT"

  END SUBROUTINE init

  SUBROUTINE restart_vector_GCV_channel_solve(itime_start)
    real, allocatable, dimension(:, :) :: tempx, tempy, tempz
    real, allocatable, dimension(:, :, :) :: rho, T
    real, allocatable, dimension(:, :, :, :) :: u
    real, allocatable, dimension(:, :, :) :: alpha_rho, alpha_T, alpha_rho_T
    real, allocatable, dimension(:, :, :) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, allocatable, dimension(:, :, :) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, allocatable, dimension(:, :, :, :) :: alpha_u, alpha_rho_u, alpha_rho_uT
    real, allocatable, dimension(:, :, :) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, allocatable, dimension(:, :, :) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, allocatable, dimension(:, :, :, :) :: output
    real, allocatable, dimension(:, :, :) :: temp1, temp2
    real :: sum_rho, sum_rho_u
    real :: gamma, Mach, Re, Pr, mu
    real    :: time
    integer :: itime
    integer :: itime_start

    real, dimension(xNc, Ny, zNc) :: bodyforce

    call channel_init(output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    sum_rho = SUM(rho)
    sum_rho_u = SUM(rho * u(:,:,:,1))

    time_start = 0.001 * itime_start
    call channel_read(rho, u, T, itime_start)

    do itime = itime_start+1,itime_start+Ntime

       time = time_start + REAL(itime-itime_start-1)*dt

       if (mod(itime, 20) .eq. 0) then
          global_rho = rho; global_u = u; global_T = T
          call vector_GCV(rho, u, T, alpha_rho, alpha_u, alpha_T, temp1, tempx, tempy, tempz)
       end if

       call vector_GCV_channel_rkw3(time, itime, bodyforce, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sum_rho, sum_rho_u)

    end do

  END SUBROUTINE restart_vector_GCV_channel_solve

  SUBROUTINE restart_GCV_channel_solve(itime_start)
    real, allocatable, dimension(:, :) :: tempx, tempy, tempz
    real, allocatable, dimension(:, :, :) :: rho, T
    real, allocatable, dimension(:, :, :, :) :: u
    real, allocatable, dimension(:, :, :) :: alpha_rho, alpha_T, alpha_rho_T
    real, allocatable, dimension(:, :, :) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, allocatable, dimension(:, :, :) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, allocatable, dimension(:, :, :, :) :: alpha_u, alpha_rho_u, alpha_rho_uT
    real, allocatable, dimension(:, :, :) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, allocatable, dimension(:, :, :) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, allocatable, dimension(:, :, :, :) :: output
    real, allocatable, dimension(:, :, :) :: temp1, temp2
    real :: sum_rho, sum_rho_u
    real :: gamma, Mach, Re, Pr, mu
    real    :: time
    integer :: itime
    integer :: itime_start

    real, dimension(xNc, Ny, zNc) :: bodyforce

    call channel_init(output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    sum_rho = SUM(rho)
    sum_rho_u = SUM(rho * u(:,:,:,1))

    time_start = 0.001 * itime_start
    call channel_read(rho, u, T, itime_start)

    do itime = itime_start+1,itime_start+Ntime

       time = time_start + REAL(itime-itime_start-1)*dt

       global_rho = rho; global_u = u; global_T = T
       call GCV(rho, u, T, alpha_rho, alpha_u, alpha_T, temp1, tempx, tempy, tempz)

       call GCV_channel_rkw3(time, itime, bodyforce, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sum_rho, sum_rho_u)

    end do

  END SUBROUTINE restart_GCV_channel_solve

  SUBROUTINE restart_pchannel_solve(itime_start)
    real, allocatable, dimension(:, :) :: tempx, tempy, tempz
    real, allocatable, dimension(:, :, :) :: rho, T
    real, allocatable, dimension(:, :, :, :) :: u
    real, allocatable, dimension(:, :, :) :: alpha_rho, alpha_T, alpha_rho_T
    real, allocatable, dimension(:, :, :) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, allocatable, dimension(:, :, :) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, allocatable, dimension(:, :, :, :) :: alpha_u, alpha_rho_u, alpha_rho_uT
    real, allocatable, dimension(:, :, :) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, allocatable, dimension(:, :, :) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, allocatable, dimension(:, :, :, :) :: output
    real, allocatable, dimension(:, :, :) :: temp1, temp2
    real :: sum_rho, sum_rho_u
    real :: gamma, Mach, Re, Pr, mu
    real    :: time
    integer :: itime
    integer :: itime_start

    real, dimension(xNc, Ny, zNc) :: bodyforce

    call channel_init(output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    sum_rho = SUM(rho)
    sum_rho_u = SUM(rho * u(:,:,:,1))

    time_start = 0.001 * itime_start
    call channel_read(rho, u, T, itime_start)

    ! re-construct penalty matrices
    lamda = 1.0e-5
    call re_operator

    do itime = itime_start+1,itime_start+Ntime

       time = time_start + REAL(itime-itime_start-1)*dt

       call pchannel_rkw3(time, itime, bodyforce, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sum_rho, sum_rho_u)

    end do

  END SUBROUTINE restart_pchannel_solve

  SUBROUTINE restart_channel_solve(itime_start)
    real, allocatable, dimension(:, :) :: tempx, tempy, tempz
    real, allocatable, dimension(:, :, :) :: rho, T
    real, allocatable, dimension(:, :, :, :) :: u
    real, allocatable, dimension(:, :, :) :: alpha_rho, alpha_T, alpha_rho_T
    real, allocatable, dimension(:, :, :) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, allocatable, dimension(:, :, :) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, allocatable, dimension(:, :, :, :) :: alpha_u, alpha_rho_u, alpha_rho_uT
    real, allocatable, dimension(:, :, :) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, allocatable, dimension(:, :, :) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, allocatable, dimension(:, :, :, :) :: output
    real, allocatable, dimension(:, :, :) :: temp1, temp2
    real :: sum_rho, sum_rho_u
    real :: gamma, Mach, Re, Pr, mu
    real    :: time
    integer :: itime
    integer :: itime_start

    real, dimension(xNc, Ny, zNc) :: bodyforce

    call channel_init(output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    sum_rho = SUM(rho)
    sum_rho_u = SUM(rho * u(:,:,:,1))

    time_start = 0.001 * itime_start
    call channel_read(rho, u, T, itime_start)

    do itime = itime_start+1,itime_start+Ntime

       time = time_start + REAL(itime-itime_start-1)*dt

       call channel_rkw3(time, itime, bodyforce, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sum_rho, sum_rho_u)

    end do

  END SUBROUTINE restart_channel_solve

  SUBROUTINE GCV_channel_solve
    real, allocatable, dimension(:, :) :: tempx, tempy, tempz
    real, allocatable, dimension(:, :, :) :: rho, T
    real, allocatable, dimension(:, :, :, :) :: u
    real, allocatable, dimension(:, :, :) :: alpha_rho, alpha_T, alpha_rho_T
    real, allocatable, dimension(:, :, :) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, allocatable, dimension(:, :, :) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, allocatable, dimension(:, :, :, :) :: alpha_u, alpha_rho_u, alpha_rho_uT
    real, allocatable, dimension(:, :, :) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, allocatable, dimension(:, :, :) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, allocatable, dimension(:, :, :, :) :: output
    real, allocatable, dimension(:, :, :) :: temp1, temp2
    real :: sum_rho, sum_rho_u
    real :: gamma, Mach, Re, Pr, mu
    real    :: time
    integer :: itime

    real, dimension(xNc, Ny, zNc) :: bodyforce

    call channel_init(output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    sum_rho = SUM(rho)
    sum_rho_u = SUM(rho * u(:,:,:,1))

    alpha_u = u
    call solve3d(alpha_u(:,:,:,1), tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
    temp1 = mu * temp1
    call channel_output(rho, u, T, temp1, 0)

    do itime = itime_start+1,itime_start+Ntime

       time = time_start + REAL(itime-1)*dt

       global_rho = rho; global_u = u; global_T = T
       call GCV(rho, u, T, alpha_rho, alpha_u, alpha_T, temp1, tempx, tempy, tempz)

       call GCV_channel_rkw3(time, itime, bodyforce, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sum_rho, sum_rho_u)

    end do

  END SUBROUTINE GCV_channel_solve

  SUBROUTINE pchannel_solve
    real, allocatable, dimension(:, :) :: tempx, tempy, tempz
    real, allocatable, dimension(:, :, :) :: rho, T
    real, allocatable, dimension(:, :, :, :) :: u
    real, allocatable, dimension(:, :, :) :: alpha_rho, alpha_T, alpha_rho_T
    real, allocatable, dimension(:, :, :) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, allocatable, dimension(:, :, :) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, allocatable, dimension(:, :, :, :) :: alpha_u, alpha_rho_u, alpha_rho_uT
    real, allocatable, dimension(:, :, :) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, allocatable, dimension(:, :, :) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, allocatable, dimension(:, :, :, :) :: output
    real, allocatable, dimension(:, :, :) :: temp1, temp2
    real :: sum_rho, sum_rho_u
    real :: gamma, Mach, Re, Pr, mu
    real    :: time
    integer :: itime

    real, dimension(xNc, Ny, zNc) :: bodyforce

    call channel_init(output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    sum_rho = SUM(rho)
    sum_rho_u = SUM(rho * u(:,:,:,1))

    alpha_u = u
    call solve3d(alpha_u(:,:,:,1), tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
    temp1 = mu * temp1
    call channel_output(rho, u, T, temp1, 0)

    ! re-construct penalty matrices
    lamda = 1.0e-5
    call re_operator

    do itime = itime_start+1,itime_start+Ntime

       time = time_start + REAL(itime-1)*dt

       call pchannel_rkw3(time, itime, bodyforce, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sum_rho, sum_rho_u)

    end do

  END SUBROUTINE pchannel_solve

  SUBROUTINE channel_solve
    real, allocatable, dimension(:, :) :: tempx, tempy, tempz
    real, allocatable, dimension(:, :, :) :: rho, T
    real, allocatable, dimension(:, :, :, :) :: u
    real, allocatable, dimension(:, :, :) :: alpha_rho, alpha_T, alpha_rho_T
    real, allocatable, dimension(:, :, :) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, allocatable, dimension(:, :, :) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, allocatable, dimension(:, :, :, :) :: alpha_u, alpha_rho_u, alpha_rho_uT
    real, allocatable, dimension(:, :, :) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, allocatable, dimension(:, :, :) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, allocatable, dimension(:, :, :, :) :: output
    real, allocatable, dimension(:, :, :) :: temp1, temp2
    real :: sum_rho, sum_rho_u
    real :: gamma, Mach, Re, Pr, mu
    real    :: time
    integer :: itime

    real, dimension(xNc, Ny, zNc) :: bodyforce

    call channel_init(output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    sum_rho = SUM(rho)
    sum_rho_u = SUM(rho * u(:,:,:,1))

    alpha_u = u
    call solve3d(alpha_u(:,:,:,1), tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
    temp1 = mu * temp1
    call channel_output(rho, u, T, temp1, 0)

    do itime = itime_start+1,itime_start+Ntime

       time = time_start + REAL(itime-1)*dt

       call channel_rkw3(time, itime, bodyforce, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sum_rho, sum_rho_u)

    end do

  END SUBROUTINE channel_solve

  SUBROUTINE mms_solve
    real, allocatable, dimension(:, :) :: tempx, tempy, tempz
    real, allocatable, dimension(:, :, :) :: rho, T
    real, allocatable, dimension(:, :, :, :) :: u
    real, allocatable, dimension(:, :, :) :: alpha_rho, alpha_T, alpha_rho_T
    real, allocatable, dimension(:, :, :) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, allocatable, dimension(:, :, :) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, allocatable, dimension(:, :, :, :) :: alpha_u, alpha_rho_u, alpha_rho_uT
    real, allocatable, dimension(:, :, :) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, allocatable, dimension(:, :, :) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, allocatable, dimension(:, :, :, :) :: sol, output
    real, allocatable, dimension(:, :, :) :: temp1, temp2
    real :: gamma, Mach, Re, Pr, mu
    real    :: time
    integer :: itime

    real, dimension(xNc, Ny, zNc) :: sol_rho, sol_T
    real, dimension(xNc, Ny, zNc, 3) :: sol_u

    call mms_init(output, sol, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    sol_rho = rho
    sol_u = u
    sol_T = T

    do itime = itime_start+1,itime_start+Ntime

       time = time_start + REAL(itime-1)*dt

       call mms_rkw3(time, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sol_rho, sol_u, sol_T)

    end do

  END SUBROUTINE mms_solve

  SUBROUTINE final

    ! do nothing

  END SUBROUTINE final

  SUBROUTINE valid
    real, dimension(Nx,Ny,Nz,5) :: dQ
    integer                   :: n


    call rhs(Q,dQ)
    open(1,file='l.in')
    write(1,"(F20.10)")X(Nx/3)
    write(1,"(F20.10)")Y(Ny/3)
    write(1,"(F20.10)")Z(Nz/3)
    close(1)

    write(6,"(9F20.10)")(dQ(Nx/3,Ny/3,Nz/3,n),n=1,5)

    stop 'VALIDATION DONE'

  END SUBROUTINE valid


  SUBROUTINE testfilter

    real, dimension(Nx,Ny,Nz,5) :: Q
    integer :: i,n

    do i = 1,Nx*Ny*Nz*5

       call RANDOM_NUMBER(Q(i,1,1,1))

    end do

    do n = 1,1
       call filterQ(Q)
    end do

    call output(Q,0)

  END SUBROUTINE testfilter

  SUBROUTINE test

    real,dimension(Nx,Ny,Nz) :: F,dF
    real,dimension(Ny)    :: Fbs, dFbs
    real                  :: exact, Linf
    integer               :: i,j,k

    do i = 1,Nx
       F(i,:,:) = EXP(SIN(X(i)) + COS(X(i)))
    end do

    call derivativeX(F,dF,1)

    do i = 1,Nx
       exact = (COS(X(i)) - SIN(X(i)))*EXP(SIN(X(i)) + COS(X(i)))
       write(88,"(6E20.10)") X(i),exact,dF(i,Ny/2,Nz/2),F(i,Ny/2,Nz/2)
    end do

    do k = 1,Nz
       do j = 1,Ny
          F(:,j,k) = Z(k)*EXP(SIN(Y(j)) + COS(Y(j)))
       end do
    end do

    print *,"x"
    call derivativeY(F,dF)
    print *,"y"

    do j = 1,Ny
       exact = Z(Nz/2)*(COS(Y(j)) - SIN(Y(j)))*EXP(SIN(Y(j)) + COS(Y(j)))
       write(89,"(6E20.10)") Y(j),exact,dF(Nx/2,j,Nz/2),F(Nx/2,j,Nz/2)
    end do



    do k = 1,Nz
       do j = 1,Ny
          F(:,j,k) = Z(k)*EXP(SIN(Yc(j)) + COS(Yc(j)))
       end do
    end do

!    call derivativeYBS(F,dF)
!
!    do j = 1,Ny
!       exact = Z(Nz/2)*(COS(Yc(j)) - SIN(Yc(j)))*EXP(SIN(Yc(j)) + COS(Yc(j)))
!       write(30,"(6E20.10)") Yc(j),exact,dF(Nx/2,j,Nz/2),F(Nx/2,j,Nz/2)
!    end do


    do j = 1,Ny
       Fbs(j) = EXP(SIN(Yc(j)) + COS(Yc(j)))
    end do

    print *,"x"
    call derivativeYBS1(1,Fbs,dFbs)
    print *,"y"

    Linf = 0.0
    do j = 1,Ny
       exact = (COS(Yc(j)) - SIN(Yc(j)))*EXP(SIN(Yc(j)) + COS(Yc(j)))
       write(20,"(6E20.10)") Yc(j),exact,dFbs(j),Fbs(j)

       Linf = MAX(ABS(exact - dFbs(j)), Linf)
    end do

    print *, Linf

    call derivativeYBS1(2,Fbs,dFbs)
    print *,"y"

    Linf = 0.0
    do j = 1,Ny
       exact = (1.0 - COS(Yc(j)) - SIN(Yc(j)) - SIN(2.0*Yc(j)))*EXP(SIN(Yc(j)) + COS(Yc(j)))
       write(20,"(6E20.10)") Yc(j),exact,dFbs(j),Fbs(j)

       Linf = MAX(ABS(exact - dFbs(j)), Linf)
    end do

    print *, Linf

    do k = 1,Nz
       do i = 1,Nx
          F(i,:,k) = X(i)*EXP(SIN(Z(k)*2.*Pi/Lz) + COS(Z(k)*2.*Pi/Lz))
       end do
    end do
    call derivativeZ(F,dF,1)

    do k = 1,Nz
       exact = X(Nx)*(COS(Z(k)*2.*Pi/Lz) - SIN(Z(k)*2.*Pi/Lz))*EXP(SIN(Z(k)*2.*Pi/Lz) + COS(Z(k)*2.*Pi/Lz))
       write(90,"(6E20.10)") Z(k),exact*2.*Pi/Lz,dF(Nx,Ny/2,k),F(Nx,Ny/2,k)
    end do

    do i = 1,Nx
       F(i,:,:) = EXP(SIN(X(i)) + COS(X(i)))
    end do
    call derivative2X(F,dF)
    do i = 1,Nx
       exact = (-SIN(X(i)) - COS(X(i)))*EXP(SIN(X(i)) + COS(X(i))) &
            + ((COS(X(i)) - SIN(X(i)))**2)*EXP(SIN(X(i)) + COS(X(i)))
       write(78,"(6E20.10)") X(i),exact,dF(i,Ny/2,Nz/2),F(i,Ny/2,Nz/2)
    end do


    do k = 1,Nz
       do j = 1,Ny
          F(:,j,k) = Z(k)*EXP(SIN(Y(j)) + COS(Y(j)))
       end do
    end do
    call derivative2Y(F,dF)
    do j = 1,Ny
       exact = Z(Nz/2)*((-SIN(Y(j)) - COS(Y(j)))*EXP(SIN(Y(j)) + COS(Y(j))) &
            + ((COS(Y(j)) - SIN(Y(j)))**2)*EXP(SIN(Y(j)) + COS(Y(j))) )
       write(79,"(6E20.10)") Y(j),exact,dF(Nx/2,j,Nz/2),F(Nx/2,j,Nz/2)
    end do


    do k = 1,Nz
       do j = 1,Ny
          F(:,j,k) = Z(k)*EXP(SIN(Yc(j)) + COS(Yc(j)))
       end do
    end do
!    call derivative2YBS(F,dF)
!    do j = 1,Ny
!       exact = Z(Nz/2)*((-SIN(Yc(j)) - COS(Yc(j)))*EXP(SIN(Yc(j)) + COS(Yc(j))) &
!            + ((COS(Yc(j)) - SIN(Yc(j)))**2)*EXP(SIN(Yc(j)) + COS(Yc(j))) )
!       write(46,"(6E20.10)") Yc(j),exact,dF(Nx/2,j,Nz/2),F(Nx/2,j,Nz/2)
!    end do

    do k = 1,Nz
       do i = 1,Nx
          F(i,:,k) = X(i)*EXP(SIN(Z(k)*2.*Pi/Lz) + COS(Z(k)*2.*Pi/Lz))
       end do
    end do
    call derivative2Z(F,dF)
    do k = 1,Nz
       exact = X(Nx)*((-SIN(Z(k)*2.*Pi/Lz) - COS(Z(k)*2.*Pi/Lz))*EXP(SIN(Z(k)*2.*Pi/Lz) + COS(Z(k)*2.*Pi/Lz)) &
            + ((COS(Z(k)*2.*Pi/Lz) - SIN(Z(k)*2.*Pi/Lz))**2)*EXP(SIN(Z(k)*2.*Pi/Lz) + COS(Z(k)*2.*Pi/Lz)))
       write(80,"(6E20.10)") Z(k),exact*(2.*Pi/Lz)**2,dF(Nx,Ny/2,k),F(Nx,Ny/2,k)
    end do

  END SUBROUTINE test
END PROGRAM mixc
