MODULE rungekutta

  USE prms
  USE data
  USE rhside
  USE bc
  USE col_linalg

  IMPLICIT none

  real                       :: fac_3,fac_6,fac_2

CONTAINS

  SUBROUTINE initRK

    fac_3 = 1./3.*dt
    fac_2 = .5*dt
    fac_6 = 1./6.*dt

  END SUBROUTINE initRK

  SUBROUTINE channel_conserve_mass(rho, u, sum_rho0, sum_rho_u0)
    real, dimension(xNc, Ny, zNc) :: rho
    real, dimension(xNc, Ny, zNc, 3) :: u
    real :: sum_rho0, sum_rho_u0

    real :: sum_rho1, sum_rho_u2, correction, sum_j
    integer :: i, j, k

    sum_rho1 = SUM(rho)
    sum_rho_u2 = SUM(rho * u(:,:,:,1))

    do k=1, zNc
       do j=1, Ny
          do i=1, xNc
             sum_j = SUM(rho(i, :, k))
             correction = sum_j - sum_rho0/sum_rho1 * sum_j
             rho(i,j,k) = rho(i,j,k) - rho(i,j,k)/sum_j * correction

             sum_j = SUM(rho(i, :, k) * u(i, :, k, 1))
             correction = sum_j - (sum_rho1/sum_rho0) * (sum_rho_u0/sum_rho_u2) * sum_j
             u(i,j,k,1) = u(i,j,k,1) - u(i,j,k,1)/sum_j * correction
          end do
       end do
    end do
  END SUBROUTINE channel_conserve_mass

  SUBROUTINE channel_read(rho, u, T, itime)
    real, dimension(xNc, Ny, zNc) :: rho, T
    real, dimension(xNc, Ny, zNc, 3) :: u
    integer :: itime

    integer :: i, j, k
    character(len=10) :: otime
    real, allocatable, dimension(:) :: xx, yy, zz
    character(len=:), allocatable :: f0, f1, f2, f3, f4

    allocate(xx(xNc), yy(Ny), zz(zNc))

    f0 = "channel/rho/"
    f1 = "channel/u/"
    f2 = "channel/v/"
    f3 = "channel/w/"
    f4 = "channel/T/"
    write(otime,"(I10.10)")itime

    open(0,file= f0//"rho"//"_"//otime//".out",form='FORMATTED')
    open(1,file= f1//"u"//"_"//otime//".out",form='FORMATTED')
    open(2,file= f2//"v"//"_"//otime//".out",form='FORMATTED')
    open(3,file= f3//"w"//"_"//otime//".out",form='FORMATTED')
    open(4,file= f4//"T"//"_"//otime//".out",form='FORMATTED')

    do k=1, zNc
       do j=1, Ny
          do i=1, xNc
             read(0,"(4F20.10)") xx(i), yy(j), zz(k), rho(i, j, k)
             read(1,"(4F20.10)") xx(i), yy(j), zz(k), u(i, j, k, 1)
             read(2,"(4F20.10)") xx(i), yy(j), zz(k), u(i, j, k, 2)
             read(3,"(4F20.10)") xx(i), yy(j), zz(k), u(i, j, k, 3)
             read(4,"(4F20.10)") xx(i), yy(j), zz(k), T(i, j, k)
          end do
       end do
    end do

    close(0)
    close(1)
    close(2)
    close(3)
    close(4)

    deallocate(xx, yy, zz)
  END SUBROUTINE channel_read

  SUBROUTINE channel_vector_lamda_output(itime)
    integer :: itime

    character(len=:), allocatable :: f0

    f0 = "channel/lamda/"

    open(0,file= f0//"lamda"//".out",form='FORMATTED',access='APPEND')

    write(0,"(4F20.10)") itime*0.001, vector_lamda(1), SUM(vector_lamda(3:Ny+2))/Ny, vector_lamda(2)

    close(0)
  END SUBROUTINE channel_vector_lamda_output

  SUBROUTINE channel_lamda_output(itime)
    integer :: itime

    character(len=:), allocatable :: f0

    f0 = "channel/lamda/"

    open(0,file= f0//"lamda"//".out",form='FORMATTED',access='APPEND')

    write(0,"(4F20.10)") itime*0.001, lamda(1), lamda(2), lamda(3)

    close(0)
  END SUBROUTINE channel_lamda_output

  SUBROUTINE channel_output(rho, u, T, tw, itime)
    real, dimension(xNc, Ny, zNc) :: rho, T, tw
    real, dimension(xNc, Ny, zNc, 3) :: u
    integer :: itime

    integer :: i, j, k
    character(len=10) :: otime
    character(len=:), allocatable :: f0, f1, f2, f3, f4, f5

    f0 = "channel/rho/"
    f1 = "channel/u/"
    f2 = "channel/v/"
    f3 = "channel/w/"
    f4 = "channel/T/"
    f5 = "channel/tw/"
    write(otime,"(I10.10)")itime

    open(0,file= f0//"rho"//"_"//otime//".out",form='FORMATTED')
    open(1,file= f1//"u"//"_"//otime//".out",form='FORMATTED')
    open(2,file= f2//"v"//"_"//otime//".out",form='FORMATTED')
    open(3,file= f3//"w"//"_"//otime//".out",form='FORMATTED')
    open(4,file= f4//"T"//"_"//otime//".out",form='FORMATTED')
    open(5,file= f5//"tw"//"_"//otime//".out",form='FORMATTED')

    do k=1, zNc
       do j=1, Ny
          do i=1, xNc
             write(0,"(4F20.10)") Xc(i), Yc(j), Zc(k), rho(i, j, k)
             write(1,"(4F20.10)") Xc(i), Yc(j), Zc(k), u(i, j, k, 1)
             write(2,"(4F20.10)") Xc(i), Yc(j), Zc(k), u(i, j, k, 2)
             write(3,"(4F20.10)") Xc(i), Yc(j), Zc(k), u(i, j, k, 3)
             write(4,"(4F20.10)") Xc(i), Yc(j), Zc(k), T(i, j, k)
             if (j .eq. 1) then
                write(5,"(4F20.10)") Xc(i), Yc(j), Zc(k), tw(i, j, k)
             end if
          end do
       end do
    end do

    close(0)
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
  END SUBROUTINE channel_output

  SUBROUTINE channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    real, dimension(xNc, Ny, zNc) :: bodyforce
    real, dimension(xNc, Ny, zNc) :: alpha_tau_12, temp1
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz
    real :: Re

    real :: up, down
    real :: N
    integer :: i, j, k

    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_12, temp1, tempx, tempy, tempz)

    N = 0.0
    up = 0.0
    down = 0.0
    do k=1, zNc
       do i=1, xNc
          j = 1
          down = down + temp1(i, j, k)
          j = Ny
          up = up + temp1(i, j, k)
          N = N + 1.0
       end do
    end do
    bodyforce = - (up/N - down/N) / (2.0 * Re)

  END SUBROUTINE channel_bodyforce

  SUBROUTINE channel_bc(u, T)
    real, dimension(xNc, Ny, zNc, 3) :: u
    real, dimension(xNc, Ny, zNc) :: T
    integer :: i, j, k

    do k=1, zNc
       do i=1, xNc
          j = 1
          u(i, j, k, 1) = 0.0
          u(i, j, k, 2) = 0.0
          u(i, j, k, 3) = 0.0
          T(i, j, k) = 1.0
          j = Ny
          u(i, j, k, 1) = 0.0
          u(i, j, k, 2) = 0.0
          u(i, j, k, 3) = 0.0
          T(i, j, k) = 1.0
       end do
    end do
  END SUBROUTINE channel_bc

  SUBROUTINE vector_GCV_channel_rkw3(time, itime, bodyforce, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sum_rho, sum_rho_u)
    real, parameter, dimension(3) :: a = (/ 8.0/15.0, 1.0/4.0, 5.0/12.0 /)
    real, parameter, dimension(3) :: b = (/ 1.0/4.0, 0.0, 3.0/4.0 /)
    real, parameter, dimension(3) :: c = (/ 0.0, 8.0/15.0, 2.0/3.0 /)

    real :: gamma, Mach, Re, Pr, mu
    real, dimension(xNc, Ny, zNc) :: rho, T
    real, dimension(xNc, Ny, zNc) :: alpha_rho, alpha_T
    real, dimension(xNc, Ny, zNc, 0:4) :: output
    real, dimension(xNc, Ny, zNc, 3) :: u, alpha_u, alpha_rho_u, alpha_rho_uT
    real, dimension(xNc, Ny, zNc) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, dimension(xNc, Ny, zNc) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, dimension(xNc, Ny, zNc) :: alpha_rho_T
    real, dimension(xNc, Ny, zNc) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, dimension(xNc, Ny, zNc) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, dimension(xNc, Ny, zNc) :: bodyforce, temp1, temp2
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz
    real :: sum_rho, sum_rho_u
    integer :: itime

    real :: time_local
    real :: time
    real, dimension(xNc, Ny, zNc, 0:4) :: rktempx, rktempy, f

    ! Stage 1
    time_local = time + c(1) * dt
    rktempx(:,:,:,0) = rho
    rktempy(:,:,:,0) = rho
    f(:,:,:,0) = rho

    rktempx(:,:,:,1:3) = u(:,:,:,1:3)
    rktempy(:,:,:,1:3) = u(:,:,:,1:3)
    f(:,:,:,1:3) = u(:,:,:,1:3)

    rktempx(:,:,:,4) = T
    rktempy(:,:,:,4) = T
    f(:,:,:,4) = T

    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rktempx = rktempx + b(1) * dt * f

    ! Stage 2
    time_local = time + c(2) * dt
    rktempy = rktempx + (a(1) - b(1)) * dt * f
    f = rktempy
    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rktempx = rktempx + b(2) * dt * f

    ! Stage 3
    time_local = time + c(3) * dt
    rktempy = rktempx  + (a(3) - b(2)) * dt * f
    f = rktempy
    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rho = rktempx(:,:,:,0) + b(3) * dt * f(:,:,:,0)
    u(:,:,:,1:3) = rktempx(:,:,:,1:3) + b(3) * dt * f(:,:,:,1:3)
    T = rktempx(:,:,:,4)  + b(3) * dt * f(:,:,:,4)

    ! set bc
    call channel_bc(u, T)

    ! conserve mass
    call channel_conserve_mass(rho, u, sum_rho, sum_rho_u)

    WRITE(*, '(A11, F8.4, A11, F7.5, A9, F7.5, A9, F7.5, A9, F7.5, A9, F7.5, A9, E9.2, A9, E9.2, A9, E9.2)') &
        'Time = ', time+dt, &
        '<Rho>= ', SUM(rho)/(xNc*Ny*zNc), &
        '<U>= ', SUM(u(:,:,:,1))/(xNc*Ny*zNc), &
        '<V>= ', SUM(u(:,:,:,2))/(xNc*Ny*zNc), &
        '<W>= ', SUM(u(:,:,:,3))/(xNc*Ny*zNc), &
        '<T>= ', SUM(T)/(xNc*Ny*zNc), &
        'lam1= ', vector_lamda(1), &
        'lam2= ', SUM(vector_lamda(3:Ny+2))/Ny, &
        'lam3= ', vector_lamda(2)

    call channel_vector_lamda_output(itime)
    if (mod(itime, 60) .eq. 0) then
       call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
       temp1 = mu * temp1
       call channel_output(rho, u, T, temp1, itime)
    end if
  END SUBROUTINE vector_GCV_channel_rkw3

  SUBROUTINE GCV_channel_rkw3(time, itime, bodyforce, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sum_rho, sum_rho_u)
    real, parameter, dimension(3) :: a = (/ 8.0/15.0, 1.0/4.0, 5.0/12.0 /)
    real, parameter, dimension(3) :: b = (/ 1.0/4.0, 0.0, 3.0/4.0 /)
    real, parameter, dimension(3) :: c = (/ 0.0, 8.0/15.0, 2.0/3.0 /)

    real :: gamma, Mach, Re, Pr, mu
    real, dimension(xNc, Ny, zNc) :: rho, T
    real, dimension(xNc, Ny, zNc) :: alpha_rho, alpha_T
    real, dimension(xNc, Ny, zNc, 0:4) :: output
    real, dimension(xNc, Ny, zNc, 3) :: u, alpha_u, alpha_rho_u, alpha_rho_uT
    real, dimension(xNc, Ny, zNc) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, dimension(xNc, Ny, zNc) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, dimension(xNc, Ny, zNc) :: alpha_rho_T
    real, dimension(xNc, Ny, zNc) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, dimension(xNc, Ny, zNc) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, dimension(xNc, Ny, zNc) :: bodyforce, temp1, temp2
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz
    real :: sum_rho, sum_rho_u
    integer :: itime

    real :: time_local
    real :: time
    real, dimension(xNc, Ny, zNc, 0:4) :: rktempx, rktempy, f

    ! Stage 1
    time_local = time + c(1) * dt
    rktempx(:,:,:,0) = rho
    rktempy(:,:,:,0) = rho
    f(:,:,:,0) = rho

    rktempx(:,:,:,1:3) = u(:,:,:,1:3)
    rktempy(:,:,:,1:3) = u(:,:,:,1:3)
    f(:,:,:,1:3) = u(:,:,:,1:3)

    rktempx(:,:,:,4) = T
    rktempy(:,:,:,4) = T
    f(:,:,:,4) = T

    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rktempx = rktempx + b(1) * dt * f

    ! Stage 2
    time_local = time + c(2) * dt
    rktempy = rktempx + (a(1) - b(1)) * dt * f
    f = rktempy
    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rktempx = rktempx + b(2) * dt * f

    ! Stage 3
    time_local = time + c(3) * dt
    rktempy = rktempx  + (a(3) - b(2)) * dt * f
    f = rktempy
    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rho = rktempx(:,:,:,0) + b(3) * dt * f(:,:,:,0)
    u(:,:,:,1:3) = rktempx(:,:,:,1:3) + b(3) * dt * f(:,:,:,1:3)
    T = rktempx(:,:,:,4)  + b(3) * dt * f(:,:,:,4)

    ! set bc
    call channel_bc(u, T)

    ! conserve mass
    call channel_conserve_mass(rho, u, sum_rho, sum_rho_u)

    WRITE(*, '(A11, F8.4, A11, F7.5, A9, F7.5, A9, F7.5, A9, F7.5, A9, F7.5, A9, E9.2, A9, E9.2, A9, E9.2)') &
        'Time = ', time+dt, &
        '<Rho>= ', SUM(rho)/(xNc*Ny*zNc), &
        '<U>= ', SUM(u(:,:,:,1))/(xNc*Ny*zNc), &
        '<V>= ', SUM(u(:,:,:,2))/(xNc*Ny*zNc), &
        '<W>= ', SUM(u(:,:,:,3))/(xNc*Ny*zNc), &
        '<T>= ', SUM(T)/(xNc*Ny*zNc), &
        'lam1= ', lamda(1), &
        'lam2= ', lamda(2), &
        'lam3= ', lamda(3)

    call channel_lamda_output(itime)
    if (mod(itime, 60) .eq. 0) then
       call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
       temp1 = mu * temp1
       call channel_output(rho, u, T, temp1, itime)
    end if
  END SUBROUTINE GCV_channel_rkw3

  SUBROUTINE pchannel_rkw3(time, itime, bodyforce, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sum_rho, sum_rho_u)
    real, parameter, dimension(3) :: a = (/ 8.0/15.0, 1.0/4.0, 5.0/12.0 /)
    real, parameter, dimension(3) :: b = (/ 1.0/4.0, 0.0, 3.0/4.0 /)
    real, parameter, dimension(3) :: c = (/ 0.0, 8.0/15.0, 2.0/3.0 /)

    real :: gamma, Mach, Re, Pr, mu
    real, dimension(xNc, Ny, zNc) :: rho, T
    real, dimension(xNc, Ny, zNc) :: alpha_rho, alpha_T
    real, dimension(xNc, Ny, zNc, 0:4) :: output
    real, dimension(xNc, Ny, zNc, 3) :: u, alpha_u, alpha_rho_u, alpha_rho_uT
    real, dimension(xNc, Ny, zNc) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, dimension(xNc, Ny, zNc) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, dimension(xNc, Ny, zNc) :: alpha_rho_T
    real, dimension(xNc, Ny, zNc) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, dimension(xNc, Ny, zNc) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, dimension(xNc, Ny, zNc) :: bodyforce, temp1, temp2
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz
    real :: sum_rho, sum_rho_u
    integer :: itime

    real :: time_local
    real :: time
    real, dimension(xNc, Ny, zNc, 0:4) :: rktempx, rktempy, f

    ! Stage 1
    time_local = time + c(1) * dt
    rktempx(:,:,:,0) = rho
    rktempy(:,:,:,0) = rho
    f(:,:,:,0) = rho

    rktempx(:,:,:,1:3) = u(:,:,:,1:3)
    rktempy(:,:,:,1:3) = u(:,:,:,1:3)
    f(:,:,:,1:3) = u(:,:,:,1:3)

    rktempx(:,:,:,4) = T
    rktempy(:,:,:,4) = T
    f(:,:,:,4) = T

    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call pns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rktempx = rktempx + b(1) * dt * f

    ! Stage 2
    time_local = time + c(2) * dt
    rktempy = rktempx + (a(1) - b(1)) * dt * f
    f = rktempy
    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call pns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rktempx = rktempx + b(2) * dt * f

    ! Stage 3
    time_local = time + c(3) * dt
    rktempy = rktempx  + (a(3) - b(2)) * dt * f
    f = rktempy
    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call pns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rho = rktempx(:,:,:,0) + b(3) * dt * f(:,:,:,0)
    u(:,:,:,1:3) = rktempx(:,:,:,1:3) + b(3) * dt * f(:,:,:,1:3)
    T = rktempx(:,:,:,4)  + b(3) * dt * f(:,:,:,4)

    ! set bc
    call channel_bc(u, T)

    ! conserve mass
    call channel_conserve_mass(rho, u, sum_rho, sum_rho_u)

    WRITE(*, '(A11, F8.4, A11, F7.5, A9, F7.5, A9, F7.5, A9, F7.5, A9, F7.5)') &
        'Time = ', time+dt, &
        '<Rho>= ', SUM(rho)/(xNc*Ny*zNc), &
        '<U>= ', SUM(u(:,:,:,1))/(xNc*Ny*zNc), &
        '<V>= ', SUM(u(:,:,:,2))/(xNc*Ny*zNc), &
        '<W>= ', SUM(u(:,:,:,3))/(xNc*Ny*zNc), &
        '<T>= ', SUM(T)/(xNc*Ny*zNc)

    if (mod(itime, 60) .eq. 0) then
       call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
       temp1 = mu * temp1
       call channel_output(rho, u, T, temp1, itime)
    end if
  END SUBROUTINE pchannel_rkw3

  SUBROUTINE channel_rkw3(time, itime, bodyforce, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sum_rho, sum_rho_u)
    real, parameter, dimension(3) :: a = (/ 8.0/15.0, 1.0/4.0, 5.0/12.0 /)
    real, parameter, dimension(3) :: b = (/ 1.0/4.0, 0.0, 3.0/4.0 /)
    real, parameter, dimension(3) :: c = (/ 0.0, 8.0/15.0, 2.0/3.0 /)

    real :: gamma, Mach, Re, Pr, mu
    real, dimension(xNc, Ny, zNc) :: rho, T
    real, dimension(xNc, Ny, zNc) :: alpha_rho, alpha_T
    real, dimension(xNc, Ny, zNc, 0:4) :: output
    real, dimension(xNc, Ny, zNc, 3) :: u, alpha_u, alpha_rho_u, alpha_rho_uT
    real, dimension(xNc, Ny, zNc) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, dimension(xNc, Ny, zNc) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, dimension(xNc, Ny, zNc) :: alpha_rho_T
    real, dimension(xNc, Ny, zNc) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, dimension(xNc, Ny, zNc) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, dimension(xNc, Ny, zNc) :: bodyforce, temp1, temp2
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz
    real :: sum_rho, sum_rho_u
    integer :: itime

    real :: time_local
    real :: time
    real, dimension(xNc, Ny, zNc, 0:4) :: rktempx, rktempy, f

    ! Stage 1
    time_local = time + c(1) * dt
    rktempx(:,:,:,0) = rho
    rktempy(:,:,:,0) = rho
    f(:,:,:,0) = rho

    rktempx(:,:,:,1:3) = u(:,:,:,1:3)
    rktempy(:,:,:,1:3) = u(:,:,:,1:3)
    f(:,:,:,1:3) = u(:,:,:,1:3)

    rktempx(:,:,:,4) = T
    rktempy(:,:,:,4) = T
    f(:,:,:,4) = T

    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rktempx = rktempx + b(1) * dt * f

    ! Stage 2
    time_local = time + c(2) * dt
    rktempy = rktempx + (a(1) - b(1)) * dt * f
    f = rktempy
    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rktempx = rktempx + b(2) * dt * f

    ! Stage 3
    time_local = time + c(3) * dt
    rktempy = rktempx  + (a(3) - b(2)) * dt * f
    f = rktempy
    ! set bc
    call channel_bc(f(:,:,:,1:3), f(:,:,:,4))

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    ! add source term
    call channel_bodyforce(bodyforce, alpha_tau_12, temp1, tempx, tempy, tempz, Re)
    f(:,:,:,1) = output(:,:,:,1) + bodyforce

    f(:,:,:,0) = output(:,:,:,0)
    f(:,:,:,2) = output(:,:,:,2)
    f(:,:,:,3) = output(:,:,:,3)
    f(:,:,:,4) = output(:,:,:,4)

    rho = rktempx(:,:,:,0) + b(3) * dt * f(:,:,:,0)
    u(:,:,:,1:3) = rktempx(:,:,:,1:3) + b(3) * dt * f(:,:,:,1:3)
    T = rktempx(:,:,:,4)  + b(3) * dt * f(:,:,:,4)

    ! set bc
    call channel_bc(u, T)

    ! conserve mass
    call channel_conserve_mass(rho, u, sum_rho, sum_rho_u)

    WRITE(*, '(A11, F8.4, A11, F7.5, A9, F7.5, A9, F7.5, A9, F7.5, A9, F7.5)') &
        'Time = ', time+dt, &
        '<Rho>= ', SUM(rho)/(xNc*Ny*zNc), &
        '<U>= ', SUM(u(:,:,:,1))/(xNc*Ny*zNc), &
        '<V>= ', SUM(u(:,:,:,2))/(xNc*Ny*zNc), &
        '<W>= ', SUM(u(:,:,:,3))/(xNc*Ny*zNc), &
        '<T>= ', SUM(T)/(xNc*Ny*zNc)

    if (mod(itime, 60) .eq. 0) then
       call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
       temp1 = mu * temp1
       call channel_output(rho, u, T, temp1, itime)
    end if
  END SUBROUTINE channel_rkw3

  SUBROUTINE mms_source(source_rho, source_u, source_T, t, gamma, Mach, Re, Pr, mu)
    real, dimension(xNc, Ny, zNc) :: source_rho, source_T
    real, dimension(xNc, Ny, zNc, 3) :: source_u
    real :: t, gamma, Mach, Re, Pr, mu
    integer :: i, j, k

    do k=1, zNc
       do j=1, Ny
          do i=1, xNc
             source_rho(i, j, k) =  &
             (R0*(-2*exp(t) + 2*W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(W3*Zc(k)) - &
             exp(t)*Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)) + &
             W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(W3*Zc(k))*Sin(R1*Xc(i))*Sin(R2*Yc(j))* &
             Sin(R3*Zc(k)) + R1*U0*Cos(R1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Sin(U1*Xc(i))* &
             Sin(R2*Yc(j))*Sin(R3*Zc(k)) + &
             U0*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))) + &
             V0*Cos(V1*Xc(i))*Cos(V3*Zc(k))* &
             (R2*Cos(R2*Yc(j))*Sin(R1*Xc(i))*Sin(V2*Yc(j))*Sin(R3*Zc(k)) + &
             V2*Cos(V2*Yc(j))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))) + &
             R3*W0*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(R3*Zc(k))*Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(W3*Zc(k)) &
             ))/exp(2*t)

             source_u(i,j,k,1) = &
             -((U0*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Sin(U1*Xc(i)))/exp(t)) - &
             (mu/Re*exp(t)*((-4*U0*U1**2*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Sin(U1*Xc(i)))/(3.*exp(t)) - &
             (U0*U2**2*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Sin(U1*Xc(i)))/exp(t) - &
             (U0*U3**2*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Sin(U1*Xc(i)))/exp(t) - &
             (V0*V1*V2*Cos(V2*Yc(j))*Cos(V3*Zc(k))*Sin(V1*Xc(i)))/exp(t) - &
             (W0*W1*W3*Cos(W2*Yc(j))*Cos(W3*Zc(k))*Sin(W1*Xc(i)))/exp(t) - &
             (2*(-((V0*V1*V2*Cos(V2*Yc(j))*Cos(V3*Zc(k))*Sin(V1*Xc(i)))/exp(t)) - &
             (W0*W1*W3*Cos(W2*Yc(j))*Cos(W3*Zc(k))*Sin(W1*Xc(i)))/exp(t)))/3.))/ &
             (R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))) + &
             (1.0/(gamma * Mach**2)*exp(t)*((R0*R1*T0*Cos(R1*Xc(i))* &
             (2 + Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))* &
             Sin(R2*Yc(j))*Sin(R3*Zc(k)))/exp(2*t) - &
             (R0*T0*T1*Cos(T2*Yc(j))*Cos(T3*Zc(k))*Sin(T1*Xc(i))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t)))/ &
             (R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))) + &
             (exp(t)*((R0*R1*U0**2*Cos(R1*Xc(i))*Cos(U2*Yc(j))**2*Cos(U3*Zc(k))**2*Sin(U1*Xc(i))**2* &
             Sin(R2*Yc(j))*Sin(R3*Zc(k)))/exp(3*t) + &
             (R0*R2*U0*V0*Cos(V1*Xc(i))*Cos(R2*Yc(j))*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Cos(V3*Zc(k))* &
             Sin(R1*Xc(i))*Sin(U1*Xc(i))*Sin(V2*Yc(j))*Sin(R3*Zc(k)))/exp(3*t) + &
             (3*R0*U0**2*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))**2*Cos(U3*Zc(k))**2*Sin(U1*Xc(i))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(3*t) + &
             (R0*U0*V0*V2*Cos(V1*Xc(i))*Cos(U2*Yc(j))*Cos(V2*Yc(j))*Cos(U3*Zc(k))*Cos(V3*Zc(k))* &
             Sin(U1*Xc(i))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(3*t) + &
             (R0*U0*W0*W3*Cos(W1*Xc(i))*Cos(U2*Yc(j))*Cos(W2*Yc(j))*Cos(U3*Zc(k))*Cos(W3*Zc(k))* &
             Sin(U1*Xc(i))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(3*t) - &
             (2*R0*U0*U2*V0*Cos(V1*Xc(i))*Cos(U3*Zc(k))*Cos(V3*Zc(k))*Sin(U1*Xc(i))* &
             Sin(U2*Yc(j))*Sin(V2*Yc(j))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/ &
             exp(3*t) + (R0*R3*U0*W0*Cos(W1*Xc(i))*Cos(U2*Yc(j))*Cos(W2*Yc(j))* &
             Cos(R3*Zc(k))*Cos(U3*Zc(k))*Sin(R1*Xc(i))*Sin(U1*Xc(i))*Sin(R2*Yc(j))*Sin(W3*Zc(k)))/ &
             exp(3*t) - (2*R0*U0*U3*W0*Cos(W1*Xc(i))*Cos(U2*Yc(j))*Cos(W2*Yc(j))* &
             Sin(U1*Xc(i))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))*Sin(U3*Zc(k))* &
             Sin(W3*Zc(k)))/exp(3*t) - &
             (U0*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Sin(U1*Xc(i))* &
             ((R0*R1*U0*Cos(R1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Sin(U1*Xc(i))*Sin(R2*Yc(j))* &
             Sin(R3*Zc(k)))/exp(2*t) + &
             (R0*R2*V0*Cos(V1*Xc(i))*Cos(R2*Yc(j))*Cos(V3*Zc(k))*Sin(R1*Xc(i))* &
             Sin(V2*Yc(j))*Sin(R3*Zc(k)))/exp(2*t) + &
             (R0*U0*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*V0*V2*Cos(V1*Xc(i))*Cos(V2*Yc(j))*Cos(V3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(W3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*R3*W0*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(R3*Zc(k))*Sin(R1*Xc(i))* &
             Sin(R2*Yc(j))*Sin(W3*Zc(k)))/exp(2*t)))/exp(t)))/ &
             (2.*R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))

             source_u(i,j,k,2) = &
             -((V0*Cos(V1*Xc(i))*Cos(V3*Zc(k))*Sin(V2*Yc(j)))/exp(t)) - &
             (mu/Re*exp(t)*(-((U0*U1*U2*Cos(U1*Xc(i))*Cos(U3*Zc(k))*Sin(U2*Yc(j)))/exp(t)) - &
             (V0*V1**2*Cos(V1*Xc(i))*Cos(V3*Zc(k))*Sin(V2*Yc(j)))/exp(t) - &
             (4*V0*V2**2*Cos(V1*Xc(i))*Cos(V3*Zc(k))*Sin(V2*Yc(j)))/(3.*exp(t)) - &
             (V0*V3**2*Cos(V1*Xc(i))*Cos(V3*Zc(k))*Sin(V2*Yc(j)))/exp(t) - &
             (W0*W2*W3*Cos(W1*Xc(i))*Cos(W3*Zc(k))*Sin(W2*Yc(j)))/exp(t) - &
             (2*(-((U0*U1*U2*Cos(U1*Xc(i))*Cos(U3*Zc(k))*Sin(U2*Yc(j)))/exp(t)) - &
             (W0*W2*W3*Cos(W1*Xc(i))*Cos(W3*Zc(k))*Sin(W2*Yc(j)))/exp(t)))/3.))/ &
             (R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))) + &
             (1.0/(gamma * Mach**2)*exp(t)*((R0*R2*T0*Cos(R2*Yc(j))*(2 + &
             Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))* &
             Sin(R1*Xc(i))*Sin(R3*Zc(k)))/exp(2*t) - &
             (R0*T0*T2*Cos(T1*Xc(i))*Cos(T3*Zc(k))*Sin(T2*Yc(j))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t)))/ &
             (R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))) + &
             (exp(t)*((R0*R1*U0*V0*Cos(R1*Xc(i))*Cos(V1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))* &
             Cos(V3*Zc(k))*Sin(U1*Xc(i))*Sin(R2*Yc(j))*Sin(V2*Yc(j))*Sin(R3*Zc(k)))/exp(3*t) &
             + (R0*R2*V0**2*Cos(V1*Xc(i))**2*Cos(R2*Yc(j))*Cos(V3*Zc(k))**2*Sin(R1*Xc(i))* &
             Sin(V2*Yc(j))**2*Sin(R3*Zc(k)))/exp(3*t) + &
             (R0*U0*U1*V0*Cos(U1*Xc(i))*Cos(V1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Cos(V3*Zc(k))* &
             Sin(V2*Yc(j))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(3*t) + &
             (3*R0*V0**2*V2*Cos(V1*Xc(i))**2*Cos(V2*Yc(j))*Cos(V3*Zc(k))**2*Sin(V2*Yc(j))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(3*t) + &
             (R0*V0*W0*W3*Cos(V1*Xc(i))*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(V3*Zc(k))*Cos(W3*Zc(k))* &
             Sin(V2*Yc(j))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(3*t) - &
             (2*R0*U0*V0*V1*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Cos(V3*Zc(k))*Sin(U1*Xc(i))* &
             Sin(V1*Xc(i))*Sin(V2*Yc(j))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/ &
             exp(3*t) + (R0*R3*V0*W0*Cos(V1*Xc(i))*Cos(W1*Xc(i))*Cos(W2*Yc(j))* &
             Cos(R3*Zc(k))*Cos(V3*Zc(k))*Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(V2*Yc(j))*Sin(W3*Zc(k)))/ &
             exp(3*t) - (2*R0*V0*V3*W0*Cos(V1*Xc(i))*Cos(W1*Xc(i))*Cos(W2*Yc(j))* &
             Sin(V2*Yc(j))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))*Sin(V3*Zc(k))* &
             Sin(W3*Zc(k)))/exp(3*t) - &
             (V0*Cos(V1*Xc(i))*Cos(V3*Zc(k))*Sin(V2*Yc(j))* &
             ((R0*R1*U0*Cos(R1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Sin(U1*Xc(i))*Sin(R2*Yc(j))* &
             Sin(R3*Zc(k)))/exp(2*t) + &
             (R0*R2*V0*Cos(V1*Xc(i))*Cos(R2*Yc(j))*Cos(V3*Zc(k))*Sin(R1*Xc(i))* &
             Sin(V2*Yc(j))*Sin(R3*Zc(k)))/exp(2*t) + &
             (R0*U0*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*V0*V2*Cos(V1*Xc(i))*Cos(V2*Yc(j))*Cos(V3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(W3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*R3*W0*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(R3*Zc(k))*Sin(R1*Xc(i))* &
             Sin(R2*Yc(j))*Sin(W3*Zc(k)))/exp(2*t)))/exp(t)))/ &
             (2.*R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))

             source_u(i,j,k,3) = &
             (1.0/(gamma * Mach**2)*exp(t)*((R0*R3*T0*Cos(R3*Zc(k))* &
             (2 + Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))* &
             Sin(R1*Xc(i))*Sin(R2*Yc(j)))/exp(2*t) - &
             (R0*T0*T3*Cos(T1*Xc(i))*Cos(T2*Yc(j))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))*Sin(T3*Zc(k)))/exp(2*t)))/ &
             (R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))) - &
             (W0*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Sin(W3*Zc(k)))/exp(t) - &
             (mu/Re*exp(t)*(-((U0*U1*U3*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Sin(U3*Zc(k)))/exp(t)) - &
             (V0*V2*V3*Cos(V1*Xc(i))*Cos(V2*Yc(j))*Sin(V3*Zc(k)))/exp(t) - &
             (2*(-((U0*U1*U3*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Sin(U3*Zc(k)))/exp(t)) - &
             (V0*V2*V3*Cos(V1*Xc(i))*Cos(V2*Yc(j))*Sin(V3*Zc(k)))/exp(t)))/3. - &
             (W0*W1**2*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Sin(W3*Zc(k)))/exp(t) - &
             (W0*W2**2*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Sin(W3*Zc(k)))/exp(t) - &
             (4*W0*W3**2*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Sin(W3*Zc(k)))/(3.*exp(t))))/ &
             (R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))) + &
             (exp(t)*((R0*R1*U0*W0*Cos(R1*Xc(i))*Cos(W1*Xc(i))*Cos(U2*Yc(j))*Cos(W2*Yc(j))* &
             Cos(U3*Zc(k))*Sin(U1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))*Sin(W3*Zc(k)))/exp(3*t) &
             + (R0*R2*V0*W0*Cos(V1*Xc(i))*Cos(W1*Xc(i))*Cos(R2*Yc(j))*Cos(W2*Yc(j))* &
             Cos(V3*Zc(k))*Sin(R1*Xc(i))*Sin(V2*Yc(j))*Sin(R3*Zc(k))*Sin(W3*Zc(k)))/exp(3*t) &
             + (R0*U0*U1*W0*Cos(U1*Xc(i))*Cos(W1*Xc(i))*Cos(U2*Yc(j))*Cos(W2*Yc(j))* &
             Cos(U3*Zc(k))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))*Sin(W3*Zc(k)))/ &
             exp(3*t) + (R0*V0*V2*W0*Cos(V1*Xc(i))*Cos(W1*Xc(i))*Cos(V2*Yc(j))* &
             Cos(W2*Yc(j))*Cos(V3*Zc(k))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))* &
             Sin(W3*Zc(k)))/exp(3*t) + &
             (3*R0*W0**2*W3*Cos(W1*Xc(i))**2*Cos(W2*Yc(j))**2*Cos(W3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))*Sin(W3*Zc(k)))/exp(3*t) - &
             (2*R0*U0*W0*W1*Cos(U2*Yc(j))*Cos(W2*Yc(j))*Cos(U3*Zc(k))*Sin(U1*Xc(i))* &
             Sin(W1*Xc(i))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))*Sin(W3*Zc(k)))/ &
             exp(3*t) - (2*R0*V0*W0*W2*Cos(V1*Xc(i))*Cos(W1*Xc(i))*Cos(V3*Zc(k))* &
             Sin(V2*Yc(j))*Sin(W2*Yc(j))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))* &
             Sin(W3*Zc(k)))/exp(3*t) + &
             (R0*R3*W0**2*Cos(W1*Xc(i))**2*Cos(W2*Yc(j))**2*Cos(R3*Zc(k))*Sin(R1*Xc(i))* &
             Sin(R2*Yc(j))*Sin(W3*Zc(k))**2)/exp(3*t) - &
             (W0*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Sin(W3*Zc(k))* &
             ((R0*R1*U0*Cos(R1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Sin(U1*Xc(i))*Sin(R2*Yc(j))* &
             Sin(R3*Zc(k)))/exp(2*t) + &
             (R0*R2*V0*Cos(V1*Xc(i))*Cos(R2*Yc(j))*Cos(V3*Zc(k))*Sin(R1*Xc(i))* &
             Sin(V2*Yc(j))*Sin(R3*Zc(k)))/exp(2*t) + &
             (R0*U0*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*V0*V2*Cos(V1*Xc(i))*Cos(V2*Yc(j))*Cos(V3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(W3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*R3*W0*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(R3*Zc(k))*Sin(R1*Xc(i))* &
             Sin(R2*Yc(j))*Sin(W3*Zc(k)))/exp(2*t)))/exp(t)))/ &
             (2.*R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))

             source_T(i,j,k) = &
             -((T0*(2 + Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k))))/exp(t)) + &
             ((gamma - 1.0)*T0*(2 + Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))* &
             ((U0*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k)))/exp(t) + &
             (V0*V2*Cos(V1*Xc(i))*Cos(V2*Yc(j))*Cos(V3*Zc(k)))/exp(t) + &
             (W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(W3*Zc(k)))/exp(t)))/exp(t) - &
             (exp(t)*((gamma * mu)/(Re * Pr))* &
             (-((T0*T1**2*Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))/exp(t)) - &
             (T0*T2**2*Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))/exp(t) - &
             (T0*T3**2*Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))/exp(t)))/ &
             (R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))) + &
             (exp(t)*((R0*R1*T0*U0*Cos(R1*Xc(i))*Cos(U2*Yc(j))* &
             (2 + Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))*Cos(U3*Zc(k))*Sin(U1*Xc(i))* &
             Sin(R2*Yc(j))*Sin(R3*Zc(k)))/exp(3*t) + &
             (R0*R2*T0*V0*Cos(V1*Xc(i))*Cos(R2*Yc(j))* &
             (2 + Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))*Cos(V3*Zc(k))*Sin(R1*Xc(i))* &
             Sin(V2*Yc(j))*Sin(R3*Zc(k)))/exp(3*t) + &
             (R0*T0*U0*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))* &
             (2 + Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))*Cos(U3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(3*t) + &
             (R0*T0*V0*V2*Cos(V1*Xc(i))*Cos(V2*Yc(j))* &
             (2 + Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))*Cos(V3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(3*t) + &
             (R0*T0*W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))* &
             (2 + Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))*Cos(W3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(3*t) - &
             (2*R0*T0*T1*U0*Cos(T2*Yc(j))*Cos(U2*Yc(j))*Cos(T3*Zc(k))*Cos(U3*Zc(k))* &
             Sin(T1*Xc(i))*Sin(U1*Xc(i))*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/ &
             exp(3*t) - (2*R0*T0*T2*V0*Cos(T1*Xc(i))*Cos(V1*Xc(i))*Cos(T3*Zc(k))* &
             Cos(V3*Zc(k))*Sin(T2*Yc(j))*Sin(V2*Yc(j))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(3*t) + &
             (R0*R3*T0*W0*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(R3*Zc(k))* &
             (2 + Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))*Sin(R1*Xc(i))*Sin(R2*Yc(j))* &
             Sin(W3*Zc(k)))/exp(3*t) - &
             (2*R0*T0*T3*W0*Cos(T1*Xc(i))*Cos(W1*Xc(i))*Cos(T2*Yc(j))*Cos(W2*Yc(j))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))*Sin(T3*Zc(k))*Sin(W3*Zc(k)))/ &
             exp(3*t) - (T0*(2 + Cos(T1*Xc(i))*Cos(T2*Yc(j))*Cos(T3*Zc(k)))* &
             ((R0*R1*U0*Cos(R1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))*Sin(U1*Xc(i))*Sin(R2*Yc(j))* &
             Sin(R3*Zc(k)))/exp(2*t) + &
             (R0*R2*V0*Cos(V1*Xc(i))*Cos(R2*Yc(j))*Cos(V3*Zc(k))*Sin(R1*Xc(i))* &
             Sin(V2*Yc(j))*Sin(R3*Zc(k)))/exp(2*t) + &
             (R0*U0*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*V0*V2*Cos(V1*Xc(i))*Cos(V2*Yc(j))*Cos(V3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(W3*Zc(k))* &
             (2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))/exp(2*t) + &
             (R0*R3*W0*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(R3*Zc(k))*Sin(R1*Xc(i))* &
             Sin(R2*Yc(j))*Sin(W3*Zc(k)))/exp(2*t)))/exp(t)))/ &
             (2.*R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k)))) - &
             ((mu * gamma * (gamma - 1.0) * Mach**2 / Re)*exp(t)* &
             ((W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(W3*Zc(k))* &
             ((-2*((U0*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k)))/exp(t) + &
             (V0*V2*Cos(V1*Xc(i))*Cos(V2*Yc(j))*Cos(V3*Zc(k)))/exp(t)))/3. + &
             (4*W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(W3*Zc(k)))/(3.*exp(t))))/exp(t) + &
             (V0*V2*Cos(V1*Xc(i))*Cos(V2*Yc(j))*Cos(V3*Zc(k))* &
             ((4*V0*V2*Cos(V1*Xc(i))*Cos(V2*Yc(j))*Cos(V3*Zc(k)))/(3.*exp(t)) - &
             (2*((U0*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k)))/exp(t) + &
             (W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(W3*Zc(k)))/exp(t)))/3.))/exp(t) &
             + (U0*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k))* &
             ((4*U0*U1*Cos(U1*Xc(i))*Cos(U2*Yc(j))*Cos(U3*Zc(k)))/(3.*exp(t)) - &
             (2*((V0*V2*Cos(V1*Xc(i))*Cos(V2*Yc(j))*Cos(V3*Zc(k)))/exp(t) + &
             (W0*W3*Cos(W1*Xc(i))*Cos(W2*Yc(j))*Cos(W3*Zc(k)))/exp(t)))/3.))/exp(t) &
             + (-((U0*U2*Cos(U3*Zc(k))*Sin(U1*Xc(i))*Sin(U2*Yc(j)))/exp(t)) - &
             (V0*V1*Cos(V3*Zc(k))*Sin(V1*Xc(i))*Sin(V2*Yc(j)))/exp(t))**2 + &
             (-((U0*U3*Cos(U2*Yc(j))*Sin(U1*Xc(i))*Sin(U3*Zc(k)))/exp(t)) - &
             (W0*W1*Cos(W2*Yc(j))*Sin(W1*Xc(i))*Sin(W3*Zc(k)))/exp(t))**2 + &
             (-((V0*V3*Cos(V1*Xc(i))*Sin(V2*Yc(j))*Sin(V3*Zc(k)))/exp(t)) - &
             (W0*W2*Cos(W1*Xc(i))*Sin(W2*Yc(j))*Sin(W3*Zc(k)))/exp(t))**2))/ &
             (R0*(2 + Sin(R1*Xc(i))*Sin(R2*Yc(j))*Sin(R3*Zc(k))))
          end do
       end do
    end do
  END SUBROUTINE mms_source

  SUBROUTINE mms_bc(rho, u, T, mms_t)
    real :: mms_t
    real, dimension(xNc, Ny, zNc) :: rho, T
    real, dimension(xNc, Ny, zNc, 3) :: u

    integer :: i, j, k

    do k=1, zNc
       do i=1, xNc
          j = 1
          rho(i, j, k)  = R0 * (SIN(R1*Xc(i)) * SIN(R2*Yc(j)) * SIN(R3*Zc(k)) + 2.0) * exp(-mms_t)
          u(i, j, k, 1) = U0 * (SIN(U1*Xc(i)) * COS(U2*Yc(j)) * COS(U3*Zc(k))) * exp(-mms_t)
          u(i, j, k, 2) = V0 * (COS(V1*Xc(i)) * SIN(V2*Yc(j)) * COS(V3*Zc(k))) * exp(-mms_t)
          u(i, j, k, 3) = W0 * (COS(W1*Xc(i)) * COS(W2*Yc(j)) * SIN(W3*Zc(k))) * exp(-mms_t)
          T(i, j, k)    = T0 * (COS(T1*Xc(i)) * COS(T2*Yc(j)) * COS(T3*Zc(k)) + 2.0) * exp(-mms_t)
          j = Ny
          rho(i, j, k)  = R0 * (SIN(R1*Xc(i)) * SIN(R2*Yc(j)) * SIN(R3*Zc(k)) + 2.0) * exp(-mms_t)
          u(i, j, k, 1) = U0 * (SIN(U1*Xc(i)) * COS(U2*Yc(j)) * COS(U3*Zc(k))) * exp(-mms_t)
          u(i, j, k, 2) = V0 * (COS(V1*Xc(i)) * SIN(V2*Yc(j)) * COS(V3*Zc(k))) * exp(-mms_t)
          u(i, j, k, 3) = W0 * (COS(W1*Xc(i)) * COS(W2*Yc(j)) * SIN(W3*Zc(k))) * exp(-mms_t)
          T(i, j, k)    = T0 * (COS(T1*Xc(i)) * COS(T2*Yc(j)) * COS(T3*Zc(k)) + 2.0) * exp(-mms_t)
       end do
    end do
  END SUBROUTINE mms_bc

  SUBROUTINE mms_rkw3(time, output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz, sol_rho, sol_u, sol_T)
    real, parameter, dimension(3) :: a = (/ 8.0/15.0, 1.0/4.0, 5.0/12.0 /)
    real, parameter, dimension(3) :: b = (/ 1.0/4.0, 0.0, 3.0/4.0 /)
    real, parameter, dimension(3) :: c = (/ 0.0, 8.0/15.0, 2.0/3.0 /)

    real :: gamma, Mach, Re, Pr, mu
    real, dimension(xNc, Ny, zNc) :: rho, T
    real, dimension(xNc, Ny, zNc) :: alpha_rho, alpha_T
    real, dimension(xNc, Ny, zNc, 0:4) :: output
    real, dimension(xNc, Ny, zNc, 3) :: u, alpha_u, alpha_rho_u, alpha_rho_uT
    real, dimension(xNc, Ny, zNc) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, dimension(xNc, Ny, zNc) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, dimension(xNc, Ny, zNc) :: alpha_rho_T
    real, dimension(xNc, Ny, zNc) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, dimension(xNc, Ny, zNc) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, dimension(xNc, Ny, zNc) :: temp1, temp2
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz
    real, dimension(xNc, Ny, zNc) :: sol_rho, sol_T
    real, dimension(xNc, Ny, zNc, 3) :: sol_u

    real :: time_local
    real :: time
    real, dimension(xNc, Ny, zNc, 0:4) :: rktempx, rktempy, f
    real, dimension(xNc, Ny, zNc) :: source_rho, source_T
    real, dimension(xNc, Ny, zNc, 3) :: source_u

    ! Stage 1
    time_local = time + c(1) * dt
    rktempx(:,:,:,0) = rho
    rktempy(:,:,:,0) = rho
    f(:,:,:,0) = rho

    rktempx(:,:,:,1:3) = u(:,:,:,1:3)
    rktempy(:,:,:,1:3) = u(:,:,:,1:3)
    f(:,:,:,1:3) = u(:,:,:,1:3)

    rktempx(:,:,:,4) = T
    rktempy(:,:,:,4) = T
    f(:,:,:,4) = T

    ! set bc
    call mms_bc(f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), time_local)

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    ! add source term
    call mms_source(source_rho, source_u, source_T, time_local, gamma, Mach, Re, Pr, mu)

    f(:,:,:,0) = output(:,:,:,0) + source_rho
    f(:,:,:,1:3) = output(:,:,:,1:3) + source_u(:,:,:,1:3)
    f(:,:,:,4) = output(:,:,:,4) + source_T
    rktempx = rktempx + b(1) * dt * f

    ! Stage 2
    time_local = time + c(2) * dt
    rktempy = rktempx + (a(1) - b(1)) * dt * f
    f = rktempy
    ! set bc
    call mms_bc(f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), time_local)

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    ! add source term
    call mms_source(source_rho, source_u, source_T, time_local, gamma, Mach, Re, Pr, mu)
    f(:,:,:,0) = output(:,:,:,0) + source_rho
    f(:,:,:,1:3) = output(:,:,:,1:3) + source_u(:,:,:,1:3)
    f(:,:,:,4) = output(:,:,:,4) + source_T
    rktempx = rktempx + b(2) * dt * f

    ! Stage 3
    time_local = time + c(3) * dt
    rktempy = rktempx  + (a(3) - b(2)) * dt * f
    f = rktempy
    ! set bc
    call mms_bc(f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), time_local)

    ! form rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, f(:,:,:,0), f(:,:,:,1:3), f(:,:,:,4), &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    ! add source term
    call mms_source(source_rho, source_u, source_T, time_local, gamma, Mach, Re, Pr, mu)
    f(:,:,:,0) = output(:,:,:,0) + source_rho
    rho = rktempx(:,:,:,0)  + b(3) * dt * f(:,:,:,0)
    f(:,:,:,1:3) = output(:,:,:,1:3) + source_u(:,:,:,1:3)
    u(:,:,:,1:3) = rktempx(:,:,:,1:3) + b(3) * dt * f(:,:,:,1:3)
    f(:,:,:,4) = output(:,:,:,4) + source_T
    T = rktempx(:,:,:,4)  + b(3) * dt * f(:,:,:,4)

    ! set bc
    call mms_bc(rho, u, T, time+dt)

    WRITE(*, '(A10, F8.5, A10, E10.3, A10, E10.3, A10, E10.3, A10, E10.3, A10, E10.3)') &
    'Time = ', time+dt, &
    'Rho = ', MAXVAL(ABS(sol_rho*exp(-(time+dt)) - rho)), &
    'U = ', MAXVAL(ABS(sol_u(:,:,:,1)*exp(-(time+dt)) - u(:,:,:,1))), &
    'V = ', MAXVAL(ABS(sol_u(:,:,:,2)*exp(-(time+dt)) - u(:,:,:,2))), &
    'W = ', MAXVAL(ABS(sol_u(:,:,:,3)*exp(-(time+dt)) - u(:,:,:,3))), &
    'T = ', MAXVAL(ABS(sol_T*exp(-(time+dt)) - T))
    print*, MAXLOC(ABS(sol_rho*exp(-(time+dt)) - rho))
    print*, MAXLOC(ABS(sol_u(:,:,:,1)*exp(-(time+dt)) - u(:,:,:,1)))
    print*, MAXLOC(ABS(sol_u(:,:,:,2)*exp(-(time+dt)) - u(:,:,:,2)))
    print*, MAXLOC(ABS(sol_u(:,:,:,3)*exp(-(time+dt)) - u(:,:,:,3)))
    print*, MAXLOC(ABS(sol_T*exp(-(time+dt)) - T))

  END SUBROUTINE mms_rkw3

!  SUBROUTINE rkw3Periodic(time)
!    real, parameter, dimension(3) :: a = (/ 8.0/15.0, 1.0/4.0, 5.0/12.0 /)
!    real, parameter, dimension(3) :: b = (/ 1.0/4.0, 0.0, 3.0/4.0 /)
!    real, parameter, dimension(3) :: c = (/ 0.0, 8.0/15.0, 2.0/3.0 /)
!    real :: time_local
!    real :: time
!    real, dimension(Nz) :: tempx, tempy
!    real, dimension(Nz) :: f
!
!    ! Stage 1
!    time_local = time + c(1) * dt
!    tempx = Gzrhs
!    tempy = Gzrhs
!    f = Gzrhs
!    call rhsPeriodicGalerkin(f, time_local)
!    tempx = tempx + b(1) * dt * f
!
!    ! Stage 2
!    time_local = time + c(2) * dt
!    tempy = tempx + (a(1) - b(1)) * dt * f
!    f = tempy
!    call rhsPeriodicGalerkin(f, time_local)
!    tempx = tempx + b(2) * dt * f
!
!    ! Stage 3
!    time_local = time + c(3) * dt
!    tempy = tempx  + (a(3) - b(2)) * dt * f
!    f = tempy
!    call rhsPeriodicGalerkin(f, time_local)
!    Gzrhs  = tempx  + b(3) * dt * f
!
!  END SUBROUTINE rkw3Periodic

!  SUBROUTINE rkw3Dirichlet(time)
!    real, parameter, dimension(3) :: a = (/ 8.0/15.0, 1.0/4.0, 5.0/12.0 /)
!    real, parameter, dimension(3) :: b = (/ 1.0/4.0, 0.0, 3.0/4.0 /)
!    real, parameter, dimension(3) :: c = (/ 0.0, 8.0/15.0, 2.0/3.0 /)
!    real :: time_local
!    real :: time
!    real, dimension(Nz) :: tempx, tempy
!    real, dimension(Nz) :: f
!
!    ! Stage 1
!    time_local = time + c(1) * dt
!    tempx = Gzrhs
!    tempy = Gzrhs
!    f = Gzrhs
!    call doGalerkinBC(f, f_mms_ptr(MINVAL(Zk), time_local), &
!                         f_mms_ptr(MAXVAL(Zk), time_local), Nz)
!    call rhsDirichletGalerkin(f, time_local)
!    tempx = tempx + b(1) * dt * f
!
!    ! Stage 2
!    time_local = time + c(2) * dt
!    tempy = tempx + (a(1) - b(1)) * dt * f
!    f = tempy
!    call doGalerkinBC(f, f_mms_ptr(MINVAL(Zk), time_local), &
!                         f_mms_ptr(MAXVAL(Zk), time_local), Nz)
!    call rhsDirichletGalerkin(f, time_local)
!    tempx = tempx + b(2) * dt * f
!
!    ! Stage 3
!    time_local = time + c(3) * dt
!    tempy = tempx  + (a(3) - b(2)) * dt * f
!    f = tempy
!    call doGalerkinBC(f, f_mms_ptr(MINVAL(Zk), time_local), &
!                         f_mms_ptr(MAXVAL(Zk), time_local), Nz)
!    call rhsDirichletGalerkin(f, time_local)
!    Gzrhs  = tempx  + b(3) * dt * f
!    call doGalerkinBC(Gzrhs, f_mms_ptr(MINVAL(Zk), time + dt), &
!                             f_mms_ptr(MAXVAL(Zk), time + dt), Nz)
!
!  END SUBROUTINE rkw3Dirichlet

!  SUBROUTINE rkw3Neumann(time)
!    real, parameter, dimension(3) :: a = (/ 8.0/15.0, 1.0/4.0, 5.0/12.0 /)
!    real, parameter, dimension(3) :: b = (/ 1.0/4.0, 0.0, 3.0/4.0 /)
!    real, parameter, dimension(3) :: c = (/ 0.0, 8.0/15.0, 2.0/3.0 /)
!    real :: time_local
!    real :: time
!    real, dimension(Ny) :: tempx, tempy
!    real, dimension(Ny) :: f
!
!    ! Stage 1
!    time_local = time + c(1) * dt
!    tempx = Gyrhs
!    tempy = Gyrhs
!    f = Gyrhs
!    call rhsNeumannGalerkin(f, time_local)
!    tempx = tempx + b(1) * dt * f
!
!    ! Stage 2
!    time_local = time + c(2) * dt
!    tempy = tempx + (a(1) - b(1)) * dt * f
!    f = tempy
!    call rhsNeumannGalerkin(f, time_local)
!    tempx = tempx + b(2) * dt * f
!
!    ! Stage 3
!    time_local = time + c(3) * dt
!    tempy = tempx  + (a(3) - b(2)) * dt * f
!    f = tempy
!    call rhsNeumannGalerkin(f, time_local)
!    Gyrhs  = tempx  + b(3) * dt * f
!  END SUBROUTINE rkw3Neumann

  SUBROUTINE rk(time)

    real,dimension(Nx,Ny,Nz,5) :: dQ,Q3,Q4
    real                       :: time_local
    real                       :: time

! Step 1 --> Euler Predictor
    time_local = time
    call rhs(Q,dQ)
    Q3 = Q + fac_2*dQ
    Q4 = Q + fac_6*dQ
    call dobc(Q3)

! Step 2 --> Euler Corrector
    time_local = time + 0.5*dt
    call rhs(Q3,dQ)
    Q3 = Q + fac_2*dQ
    Q4 = Q4 + fac_3*dQ
    call dobc(Q3)

! Step 3 --> Leapfrog Predictor
    time_local = time + 0.5*dt
    call rhs(Q3,dQ)
    Q3 = Q + dt*dQ
    Q4 = Q4 + fac_3*dQ
    call dobc(Q3)

! Step 4 --> Milne Corrector
    time_local = time + dt
    call rhs(Q3,dQ)
    Q = Q4 + fac_6*dQ
    call dobc(Q)

    time = time + dt

  END SUBROUTINE rk

END MODULE rungekutta
