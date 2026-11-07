MODULE col_linalg

  USE prms
  USE col_fft
  USE col_operator

  IMPLICIT none

  real :: R0, R1, R2, R3
  real :: U0, U1, U2, U3
  real :: V0, V1, V2, V3
  real :: W0, W1, W2, W3
  real :: T0, T1, T2, T3

CONTAINS

  SUBROUTINE GCV_lamda(x)
    real ( kind = 8 ) x(3)

    real :: tol
    integer :: i

    tol = 5.0e-5
    do i=1, 3
       if (x(i) .lt. 0 .or. x(i) .gt. tol) then
          lamda(i) = 0.0
       else
          lamda(i) = x(i)
       end if
    end do
  END SUBROUTINE GCV_lamda

  SUBROUTINE vector_GCV_lamda(x)
    real ( kind = 8 ) x(82)

    real :: tol
    integer :: i

    tol = 5.0e-5
    do i=1, 82
       if (x(i) .lt. 0 .or. x(i) .gt. tol) then
          vector_lamda(i) = 0.0
       else
          vector_lamda(i) = x(i)
       end if
    end do
  END SUBROUTINE vector_GCV_lamda

  SUBROUTINE vector_GCV(rho, u, T, alpha_rho, alpha_u, alpha_T, temp1, tempx, tempy, tempz)
    implicit none

    real, dimension(xNc, Ny, zNc) :: rho, T
    real, dimension(xNc, Ny, zNc) :: alpha_rho, alpha_T
    real, dimension(xNc, Ny, zNc, 3) :: u, alpha_u
    real, dimension(xNc, Ny, zNc) :: temp1
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz


    integer ( kind = 4 ), parameter :: n = 80 + 2

    real ( kind = 8 ) pr
    real ( kind = 8 ) praxis
    real ( kind = 8 ), external :: GCV_vector_u_praxis
    real ( kind = 8 ), external :: GCV_vector_v_praxis
    real ( kind = 8 ), external :: GCV_vector_w_praxis
    real ( kind = 8 ), external :: GCV_vector_rho_praxis
    real ( kind = 8 ), external :: GCV_vector_T_praxis
    real ( kind = 8 ) x(n)

!    ! rho
!    x(:) = 1e-5
!    pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_vector_rho_praxis )
!    call vector_GCV_lamda(x)
!    call vector_re_operator
!    alpha_rho = rho
!    call psolve3d(alpha_rho, temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
!    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_rho, rho, tempx, tempy, tempz)

    ! u
    x(:) = 1e-5
    pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_vector_u_praxis )
    call vector_GCV_lamda(x)
    call vector_re_operator
    alpha_u(:,:,:,1) = u(:,:,:,1)
    call psolve3d(alpha_u(:,:,:,1), temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,1), u(:,:,:,1), tempx, tempy, tempz)

!    ! v
!    if (SUM(global_u(:,:,:,2)) .gt. 1e-5) then
!       x(:) = 1e-5
!       pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_vector_v_praxis )
!       call vector_GCV_lamda(x)
!       call vector_re_operator
!       alpha_u(:,:,:,2) = u(:,:,:,2)
!       call psolve3d(alpha_u(:,:,:,2), temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
!       call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,2), u(:,:,:,2), tempx, tempy, tempz)
!    end if
!
!    ! w
!    if (SUM(global_u(:,:,:,3)) .gt. 1e-5) then
!       x(:) = 1e-5
!       pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_vector_w_praxis )
!       call vector_GCV_lamda(x)
!       call vector_re_operator
!       alpha_u(:,:,:,3) = u(:,:,:,3)
!       call psolve3d(alpha_u(:,:,:,3), temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
!       call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,3), u(:,:,:,3), tempx, tempy, tempz)
!    end if
!
!    ! T
!    x(:) = 1e-5
!    pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_vector_T_praxis )
!    call vector_GCV_lamda(x)
!    call vector_re_operator
!    alpha_T = T
!    call psolve3d(alpha_T, temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
!    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_T, T, tempx, tempy, tempz)

    return
  END SUBROUTINE vector_GCV

  SUBROUTINE GCV(rho, u, T, alpha_rho, alpha_u, alpha_T, temp1, tempx, tempy, tempz)
    implicit none

    real, dimension(xNc, Ny, zNc) :: rho, T
    real, dimension(xNc, Ny, zNc) :: alpha_rho, alpha_T
    real, dimension(xNc, Ny, zNc, 3) :: u, alpha_u
    real, dimension(xNc, Ny, zNc) :: temp1
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz


    integer ( kind = 4 ), parameter :: n = 3

    real ( kind = 8 ) pr
    real ( kind = 8 ) praxis
    real ( kind = 8 ), external :: GCV_u_praxis
    real ( kind = 8 ), external :: GCV_v_praxis
    real ( kind = 8 ), external :: GCV_w_praxis
    real ( kind = 8 ), external :: GCV_rho_praxis
    real ( kind = 8 ), external :: GCV_T_praxis
    real ( kind = 8 ) x(n)

!    ! rho
!    x(1) = 1e-5; x(2) = 1e-5; x(3) = 1e-5
!    pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_rho_praxis )
!    call GCV_lamda(x)
!    call re_operator
!    alpha_rho = rho
!    call psolve3d(alpha_rho, temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
!    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_rho, rho, tempx, tempy, tempz)

    ! u
    x(1) = 1e-5; x(2) = 1e-5; x(3) = 1e-5
    pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_u_praxis )
    call GCV_lamda(x)
    call re_operator
    alpha_u(:,:,:,1) = u(:,:,:,1)
    call psolve3d(alpha_u(:,:,:,1), temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,1), u(:,:,:,1), tempx, tempy, tempz)

    ! v
    if (SUM(global_u(:,:,:,2)) .gt. 1e-5) then
       x(1) = 1e-5; x(2) = 1e-5; x(3) = 1e-5
       pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_v_praxis )
       call GCV_lamda(x)
       call re_operator
       alpha_u(:,:,:,2) = u(:,:,:,2)
       call psolve3d(alpha_u(:,:,:,2), temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
       call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,2), u(:,:,:,2), tempx, tempy, tempz)
    end if

    ! w
    if (SUM(global_u(:,:,:,3)) .gt. 1e-5) then
       x(1) = 1e-5; x(2) = 1e-5; x(3) = 1e-5
       pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_w_praxis )
       call GCV_lamda(x)
       call re_operator
       alpha_u(:,:,:,3) = u(:,:,:,3)
       call psolve3d(alpha_u(:,:,:,3), temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
       call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,3), u(:,:,:,3), tempx, tempy, tempz)
    end if

!    ! T
!    x(1) = 1e-5; x(2) = 1e-5; x(3) = 1e-5
!    pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_T_praxis )
!    call GCV_lamda(x)
!    alpha_T = T
!    call psolve3d(alpha_T, temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
!    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_T, T, tempx, tempy, tempz)

    return
  END SUBROUTINE GCV

  SUBROUTINE nonlinear_alpha(mu, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, tempx, tempy, tempz)
    real :: mu
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz
    real, dimension(xNc, Ny, zNc) :: alpha_rho, alpha_T, alpha_rho_T
    real, dimension(xNc, Ny, zNc) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, dimension(xNc, Ny, zNc) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, dimension(xNc, Ny, zNc, 3) :: alpha_u, alpha_rho_u, alpha_rho_uT
    real, dimension(xNc, Ny, zNc) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, dimension(xNc, Ny, zNc) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, dimension(xNc, Ny, zNc) :: temp1
    integer :: i

    ! find alpha_rho_u/v/w
    do i=1, 3
       call alpha_2_to_1(alpha_rho_u(:,:,:,i), alpha_rho, alpha_u(:,:,:,i), &
             temp1, tempx, tempy, tempz)
       call alpha_3_to_1(alpha_rho_uT(:,:,:,i), alpha_rho, alpha_u(:,:,:,i), &
             alpha_T, temp1, tempx, tempy, tempz)
    end do

    call alpha_2_to_1(alpha_rho_T, alpha_rho, alpha_T, &
          temp1, tempx, tempy, tempz)

    ! find alpha_rho_uu/uv/uw/vv/vw/ww/
    call alpha_3_to_1(alpha_rho_uu, alpha_rho, alpha_u(:,:,:,1), alpha_u(:,:,:,1),&
          temp1, tempx, tempy, tempz)
    call alpha_3_to_1(alpha_rho_uv, alpha_rho, alpha_u(:,:,:,1), alpha_u(:,:,:,2),&
          temp1, tempx, tempy, tempz)
    call alpha_3_to_1(alpha_rho_uw, alpha_rho, alpha_u(:,:,:,1), alpha_u(:,:,:,3),&
          temp1, tempx, tempy, tempz)
    call alpha_3_to_1(alpha_rho_vv, alpha_rho, alpha_u(:,:,:,2), alpha_u(:,:,:,2),&
          temp1, tempx, tempy, tempz)
    call alpha_3_to_1(alpha_rho_vw, alpha_rho, alpha_u(:,:,:,2), alpha_u(:,:,:,3),&
          temp1, tempx, tempy, tempz)
    call alpha_3_to_1(alpha_rho_ww, alpha_rho, alpha_u(:,:,:,3), alpha_u(:,:,:,3),&
          temp1, tempx, tempy, tempz)

    ! find tau
    call initialize_tau(alpha_tau_11, alpha_tau_12, alpha_tau_13, alpha_tau_22, &
        alpha_tau_23, alpha_tau_33, alpha_u, mu, temp1, tempx, tempy, tempz)
  END SUBROUTINE nonlinear_alpha

  SUBROUTINE channel_init(output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    real :: gamma, Mach, Re, Pr, mu
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
    real :: random
    integer :: i, j, k

    allocate(alpha_rho(xNc, Ny, zNc), alpha_T(xNc, Ny, zNc), alpha_rho_T(xNc, Ny, zNc))
    allocate(alpha_tau_11(xNc, Ny, zNc), alpha_tau_12(xNc, Ny, zNc), alpha_tau_13(xNc, Ny, zNc))
    allocate(alpha_tau_22(xNc, Ny, zNc), alpha_tau_23(xNc, Ny, zNc), alpha_tau_33(xNc, Ny, zNc))
    allocate(alpha_u(xNc, Ny, zNc, 3), alpha_rho_uT(xNc, Ny, zNc, 3))
    allocate(alpha_rho_u(xNc, Ny, zNc, 3))
    allocate(rho(xNc, Ny, zNc), u(xNc, Ny, zNc, 3), T(xNc, Ny, zNc))
    allocate(alpha_rho_uu(xNc, Ny, zNc), alpha_rho_uv(xNc, Ny, zNc), alpha_rho_uw(xNc, Ny, zNc))
    allocate(alpha_rho_vv(xNc, Ny, zNc), alpha_rho_vw(xNc, Ny, zNc), alpha_rho_ww(xNc, Ny, zNc))
    allocate(tempx(xNc, 2), tempy(Ny, 2), tempz(zNc, 2))
    allocate(temp1(xNc, Ny, zNc), temp2(xNc, Ny, zNc))
    allocate(output(xNc, Ny, zNc, 0:4))

    gamma = 1.4; Mach = 1.5; Re=3000.0; Pr=0.72
    mu = (1.0 + 110.4/293.15)/(1.0 + 110.4/293.15)

    call generate_operator

    rho = 1.0
    u(:, :, :, 2:3) = 0.0

    call RANDOM_SEED

    ! create rhs
    do k=1, zNc
       do j=1, Ny
          do i=1, xNc
             call RANDOM_NUMBER(random)
             random = 2.0 * (random - 0.5) * 0.1
             u(i, j, k, 1) = 1.5 * (1.0 - (Yc(j) - 1.0)**2) * (1.0 + random)
             T(i, j, k) = 1.0 + ((gamma - 1.0) * Pr * Mach**2)/3.0 * &
                          1.5 * (1.0 - (Yc(j) - 1.0)**4)
          end do
       end do
    end do
  END SUBROUTINE channel_init

  SUBROUTINE mms_init(output, sol, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    real :: gamma, Mach, Re, Pr, mu
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

    integer :: i, j, k

    R0 = 1.0; R1 = 6.0; R2 = 6.0; R3 = 6.0
    U0 = 1.0; U1 = 8.0; U2 = 8.0; U3 = 8.0
    V0 = 1.0; V1 = 6.0; V2 = 6.0; V3 = 6.0
    W0 = 1.0; W1 = 8.0; W2 = 8.0; W3 = 8.0
    T0 = 1.0; T1 = 6.0; T2 = 6.0; T3 = 6.0

    allocate(alpha_rho(xNc, Ny, zNc), alpha_T(xNc, Ny, zNc), alpha_rho_T(xNc, Ny, zNc))
    allocate(alpha_tau_11(xNc, Ny, zNc), alpha_tau_12(xNc, Ny, zNc), alpha_tau_13(xNc, Ny, zNc))
    allocate(alpha_tau_22(xNc, Ny, zNc), alpha_tau_23(xNc, Ny, zNc), alpha_tau_33(xNc, Ny, zNc))
    allocate(alpha_u(xNc, Ny, zNc, 3), alpha_rho_uT(xNc, Ny, zNc, 3))
    allocate(alpha_rho_u(xNc, Ny, zNc, 3))
    allocate(rho(xNc, Ny, zNc), u(xNc, Ny, zNc, 3), T(xNc, Ny, zNc))
    allocate(alpha_rho_uu(xNc, Ny, zNc), alpha_rho_uv(xNc, Ny, zNc), alpha_rho_uw(xNc, Ny, zNc))
    allocate(alpha_rho_vv(xNc, Ny, zNc), alpha_rho_vw(xNc, Ny, zNc), alpha_rho_ww(xNc, Ny, zNc))
    allocate(tempx(xNc, 2), tempy(Ny, 2), tempz(zNc, 2))
    allocate(temp1(xNc, Ny, zNc), temp2(xNc, Ny, zNc))
    allocate(output(xNc, Ny, zNc, 0:4), sol(xNc, Ny, zNc, 0:4))

    gamma = 1.4; Mach = 1.5; Re=1.0; Pr=0.72
    mu = (1.0 + 110.4/293.15)/(1.0 + 110.4/293.15)

    call generate_operator

    ! create rhs
    do k=1, zNc
       do j=1, Ny
          do i=1, xNc
             rho(i, j, k)  = R0 * (SIN(R1*Xc(i)) * SIN(R2*Yc(j)) * SIN(R3*Zc(k)) + 2.0)
             u(i, j, k, 1) = U0 * (SIN(U1*Xc(i)) * COS(U2*Yc(j)) * COS(U3*Zc(k)))
             u(i, j, k, 2) = V0 * (COS(V1*Xc(i)) * SIN(V2*Yc(j)) * COS(V3*Zc(k)))
             u(i, j, k, 3) = W0 * (COS(W1*Xc(i)) * COS(W2*Yc(j)) * SIN(W3*Zc(k)))
             T(i, j, k)    = T0 * (COS(T1*Xc(i)) * COS(T2*Yc(j)) * COS(T3*Zc(k)) + 2.0)
          end do
       end do
    end do
  END SUBROUTINE mms_init

  SUBROUTINE pns_rhs(output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    real :: gamma, Mach, Re, Pr, mu
    real, dimension(xNc, Ny, zNc) :: rho, T, alpha_T, alpha_rho
    real, dimension(xNc, Ny, zNc, 0:4) :: output
    real, dimension(xNc, Ny, zNc, 3) :: alpha_u, alpha_rho_u, u, alpha_rho_uT
    real, dimension(xNc, Ny, zNc) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, dimension(xNc, Ny, zNc) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, dimension(xNc, Ny, zNc) :: alpha_rho_T
    real, dimension(xNc, Ny, zNc) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, dimension(xNc, Ny, zNc) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, dimension(xNc, Ny, zNc) :: temp1, temp2
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz

    integer :: i

    ! filter
    alpha_rho = rho
    alpha_u = u
    alpha_T = T

    call psolve3d(alpha_rho, temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_rho, rho, tempx, tempy, tempz)

    call psolve3d(alpha_T, temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_T, T, tempx, tempy, tempz)

    do i=1, 3
       call psolve3d(alpha_u(:,:,:,i), temp1, tempx(:, 1), tempy(:, 1), tempz(:, 1))
       call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
           alpha_u(:,:,:,i), u(:,:,:,i), tempx, tempy, tempz)
    end do

    ! initialize
    alpha_rho = rho
    alpha_u = u
    alpha_T = T
    alpha_rho_T = rho * T
    alpha_rho_u(:,:,:,1) = rho * u(:,:,:,1)
    alpha_rho_u(:,:,:,2) = rho * u(:,:,:,2)
    alpha_rho_u(:,:,:,3) = rho * u(:,:,:,3)
    alpha_rho_uT(:,:,:,1) = alpha_rho_u(:,:,:,1) * T
    alpha_rho_uT(:,:,:,2) = alpha_rho_u(:,:,:,2) * T
    alpha_rho_uT(:,:,:,3) = alpha_rho_u(:,:,:,3) * T
    alpha_rho_uu = alpha_rho_u(:,:,:,1) * u(:,:,:,1)
    alpha_rho_uv = alpha_rho_u(:,:,:,1) * u(:,:,:,2)
    alpha_rho_uw = alpha_rho_u(:,:,:,1) * u(:,:,:,3)
    alpha_rho_vv = alpha_rho_u(:,:,:,2) * u(:,:,:,2)
    alpha_rho_vw = alpha_rho_u(:,:,:,2) * u(:,:,:,3)
    alpha_rho_ww = alpha_rho_u(:,:,:,3) * u(:,:,:,3)

    ! solve
    call solve3d(alpha_rho, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_T, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_T, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_uu, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_uv, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_uw, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_vv, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_vw, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_ww, tempx(:, 1), tempy(:, 1), tempz(:, 1))

    do i=1, 3
       call solve3d(alpha_u(:,:,:,i), tempx(:, 1), tempy(:, 1), tempz(:, 1))
       call solve3d(alpha_rho_u(:,:,:,i), tempx(:, 1), tempy(:, 1), tempz(:, 1))
       call solve3d(alpha_rho_uT(:,:,:,i), tempx(:, 1), tempy(:, 1), tempz(:, 1))
    end do

    ! find tau
    call initialize_tau(alpha_tau_11, alpha_tau_12, alpha_tau_13, alpha_tau_22, &
        alpha_tau_23, alpha_tau_33, alpha_u, mu, temp1, tempx, tempy, tempz)

    output = 0.0
    !!! mass
    ! d/dx (alpha_rho_u) + d/dy (alpha_rho_v) +  d/dz (alpha_rho_w)
    call divergence(output(:,:,:,0), alpha_rho_u(:,:,:,1), alpha_rho_u(:,:,:,2), &
        alpha_rho_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,0) = - output(:,:,:,0)

    !!! momentum i=1
    ! convective term
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + rho * u(:,:,:,1) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + rho * u(:,:,:,2) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + rho * u(:,:,:,3) * temp1

    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_rho_uu, temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_rho_uv, temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_rho_uw, temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + temp1

    call divergence(temp2, alpha_rho_u(:,:,:,1), alpha_rho_u(:,:,:,2), &
        alpha_rho_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,1) = -(0.5 / rho) * (output(:,:,:,1) - u(:,:,:,1) * temp2)

    ! pressure
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_rho_T, temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) - (1.0/(gamma* Mach**2 * rho)) * temp1

    ! viscous
    call divergence(temp2, alpha_tau_11, alpha_tau_12, alpha_tau_13, temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + (1.0/(Re* rho)) * temp2

    !!! i=2
    ! convective term
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,2), temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + rho * u(:,:,:,1) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,2), temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + rho * u(:,:,:,2) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,2), temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + rho * u(:,:,:,3) * temp1

    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_rho_uv, temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_rho_vv, temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_rho_vw, temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + temp1

    call divergence(temp2, alpha_rho_u(:,:,:,1), alpha_rho_u(:,:,:,2), &
        alpha_rho_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,2) = -(0.5 / rho) * (output(:,:,:,2) - u(:,:,:,2) * temp2)

    ! pressure
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_rho_T, temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) - (1.0/(gamma* Mach**2 * rho)) * temp1

    ! viscous
    call divergence(temp2, alpha_tau_12, alpha_tau_22, alpha_tau_23, temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + (1.0/(Re* rho)) * temp2

    !!! i=3
    ! convective term
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + rho * u(:,:,:,1) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + rho * u(:,:,:,2) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + rho * u(:,:,:,3) * temp1

    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_rho_uw, temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_rho_vw, temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_rho_ww, temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + temp1

    call divergence(temp2, alpha_rho_u(:,:,:,1), alpha_rho_u(:,:,:,2), &
        alpha_rho_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,3) = -(0.5 / rho) * (output(:,:,:,3) - u(:,:,:,3) * temp2)

    ! pressure
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_rho_T, temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) - (1.0/(gamma* Mach**2 * rho)) * temp1

    ! viscous
    call divergence(temp2, alpha_tau_13, alpha_tau_23, alpha_tau_33, temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + (1.0/(Re* rho)) * temp2


    !!! Energy
    ! convective
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_T, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + rho * u(:,:,:,1) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_T, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + rho * u(:,:,:,2) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_T, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + rho * u(:,:,:,3) * temp1

    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_rho_uT(:,:,:,1), temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_rho_uT(:,:,:,2), temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_rho_uT(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + temp1

    call divergence(temp2, alpha_rho_u(:,:,:,1), alpha_rho_u(:,:,:,2), &
        alpha_rho_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,4) = -(0.5 / rho) * (output(:,:,:,4) - T * temp2)

    ! divergence
    call divergence(temp2, alpha_u(:,:,:,1), alpha_u(:,:,:,2), &
        alpha_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) - (gamma - 1.0) * (T * temp2)

    ! viscous
    ! tau_11 * dudx
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_11, temp2, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp2 * temp1)

    ! tau_22 * dvdy
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,2), temp1, tempx, tempy, tempz)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_22, temp2, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp2 * temp1)

    ! tau_33 * dwdz
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,3), temp1, tempx, tempy, tempz)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_33, temp2, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp2 * temp1)

    ! tau_12 * (dudy + dvdx)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_12, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp1**2/mu)

    ! tau_13 * (dudz + dwdx)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_13, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp1**2/mu)

    ! tau_23 * (dvdz + dwdy)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_23, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp1**2/mu)

    ! Laplacian
    call laplacian(temp2, alpha_T, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*mu)/(Re*Pr*rho) * temp2

  END SUBROUTINE pns_rhs

  SUBROUTINE ns_rhs(output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)
    real :: gamma, Mach, Re, Pr, mu
    real, dimension(xNc, Ny, zNc) :: rho, T, alpha_T, alpha_rho
    real, dimension(xNc, Ny, zNc, 0:4) :: output
    real, dimension(xNc, Ny, zNc, 3) :: alpha_u, alpha_rho_u, u, alpha_rho_uT
    real, dimension(xNc, Ny, zNc) :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
    real, dimension(xNc, Ny, zNc) :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
    real, dimension(xNc, Ny, zNc) :: alpha_rho_T
    real, dimension(xNc, Ny, zNc) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, dimension(xNc, Ny, zNc) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, dimension(xNc, Ny, zNc) :: temp1, temp2
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz

    integer :: i

    ! initialize
    alpha_rho = rho
    alpha_u = u
    alpha_T = T
    alpha_rho_T = rho * T
    alpha_rho_u(:,:,:,1) = rho * u(:,:,:,1)
    alpha_rho_u(:,:,:,2) = rho * u(:,:,:,2)
    alpha_rho_u(:,:,:,3) = rho * u(:,:,:,3)
    alpha_rho_uT(:,:,:,1) = alpha_rho_u(:,:,:,1) * T
    alpha_rho_uT(:,:,:,2) = alpha_rho_u(:,:,:,2) * T
    alpha_rho_uT(:,:,:,3) = alpha_rho_u(:,:,:,3) * T
    alpha_rho_uu = alpha_rho_u(:,:,:,1) * u(:,:,:,1)
    alpha_rho_uv = alpha_rho_u(:,:,:,1) * u(:,:,:,2)
    alpha_rho_uw = alpha_rho_u(:,:,:,1) * u(:,:,:,3)
    alpha_rho_vv = alpha_rho_u(:,:,:,2) * u(:,:,:,2)
    alpha_rho_vw = alpha_rho_u(:,:,:,2) * u(:,:,:,3)
    alpha_rho_ww = alpha_rho_u(:,:,:,3) * u(:,:,:,3)

    ! solve
    call solve3d(alpha_rho, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_T, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_T, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_uu, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_uv, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_uw, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_vv, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_vw, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call solve3d(alpha_rho_ww, tempx(:, 1), tempy(:, 1), tempz(:, 1))

    do i=1, 3
       call solve3d(alpha_u(:,:,:,i), tempx(:, 1), tempy(:, 1), tempz(:, 1))
       call solve3d(alpha_rho_u(:,:,:,i), tempx(:, 1), tempy(:, 1), tempz(:, 1))
       call solve3d(alpha_rho_uT(:,:,:,i), tempx(:, 1), tempy(:, 1), tempz(:, 1))
    end do

    ! find tau
    call initialize_tau(alpha_tau_11, alpha_tau_12, alpha_tau_13, alpha_tau_22, &
        alpha_tau_23, alpha_tau_33, alpha_u, mu, temp1, tempx, tempy, tempz)

    output = 0.0
    !!! mass
    ! d/dx (alpha_rho_u) + d/dy (alpha_rho_v) +  d/dz (alpha_rho_w)
    call divergence(output(:,:,:,0), alpha_rho_u(:,:,:,1), alpha_rho_u(:,:,:,2), &
        alpha_rho_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,0) = - output(:,:,:,0)

    !!! momentum i=1
    ! convective term
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + rho * u(:,:,:,1) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + rho * u(:,:,:,2) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + rho * u(:,:,:,3) * temp1

    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_rho_uu, temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_rho_uv, temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_rho_uw, temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + temp1

    call divergence(temp2, alpha_rho_u(:,:,:,1), alpha_rho_u(:,:,:,2), &
        alpha_rho_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,1) = -(0.5 / rho) * (output(:,:,:,1) - u(:,:,:,1) * temp2)

    ! pressure
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_rho_T, temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) - (1.0/(gamma* Mach**2 * rho)) * temp1

    ! viscous
    call divergence(temp2, alpha_tau_11, alpha_tau_12, alpha_tau_13, temp1, tempx, tempy, tempz)
    output(:,:,:,1) = output(:,:,:,1) + (1.0/(Re* rho)) * temp2

    !!! i=2
    ! convective term
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,2), temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + rho * u(:,:,:,1) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,2), temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + rho * u(:,:,:,2) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,2), temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + rho * u(:,:,:,3) * temp1

    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_rho_uv, temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_rho_vv, temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_rho_vw, temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + temp1

    call divergence(temp2, alpha_rho_u(:,:,:,1), alpha_rho_u(:,:,:,2), &
        alpha_rho_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,2) = -(0.5 / rho) * (output(:,:,:,2) - u(:,:,:,2) * temp2)

    ! pressure
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_rho_T, temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) - (1.0/(gamma* Mach**2 * rho)) * temp1

    ! viscous
    call divergence(temp2, alpha_tau_12, alpha_tau_22, alpha_tau_23, temp1, tempx, tempy, tempz)
    output(:,:,:,2) = output(:,:,:,2) + (1.0/(Re* rho)) * temp2

    !!! i=3
    ! convective term
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + rho * u(:,:,:,1) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + rho * u(:,:,:,2) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + rho * u(:,:,:,3) * temp1

    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_rho_uw, temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_rho_vw, temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_rho_ww, temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + temp1

    call divergence(temp2, alpha_rho_u(:,:,:,1), alpha_rho_u(:,:,:,2), &
        alpha_rho_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,3) = -(0.5 / rho) * (output(:,:,:,3) - u(:,:,:,3) * temp2)

    ! pressure
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_rho_T, temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) - (1.0/(gamma* Mach**2 * rho)) * temp1

    ! viscous
    call divergence(temp2, alpha_tau_13, alpha_tau_23, alpha_tau_33, temp1, tempx, tempy, tempz)
    output(:,:,:,3) = output(:,:,:,3) + (1.0/(Re* rho)) * temp2


    !!! Energy
    ! convective
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_T, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + rho * u(:,:,:,1) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_T, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + rho * u(:,:,:,2) * temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_T, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + rho * u(:,:,:,3) * temp1

    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_rho_uT(:,:,:,1), temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_rho_uT(:,:,:,2), temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + temp1
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_rho_uT(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + temp1

    call divergence(temp2, alpha_rho_u(:,:,:,1), alpha_rho_u(:,:,:,2), &
        alpha_rho_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,4) = -(0.5 / rho) * (output(:,:,:,4) - T * temp2)

    ! divergence
    call divergence(temp2, alpha_u(:,:,:,1), alpha_u(:,:,:,2), &
        alpha_u(:,:,:,3), temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) - (gamma - 1.0) * (T * temp2)

    ! viscous
    ! tau_11 * dudx
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,1), temp1, tempx, tempy, tempz)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_11, temp2, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp2 * temp1)

    ! tau_22 * dvdy
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,2), temp1, tempx, tempy, tempz)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_22, temp2, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp2 * temp1)

    ! tau_33 * dwdz
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,3), temp1, tempx, tempy, tempz)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_33, temp2, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp2 * temp1)

    ! tau_12 * (dudy + dvdx)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_12, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp1**2/mu)

    ! tau_13 * (dudz + dwdx)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_13, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp1**2/mu)

    ! tau_23 * (dvdz + dwdy)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha_tau_23, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*(gamma-1.0)*Mach**2/(Re*rho))*(temp1**2/mu)

    ! Laplacian
    call laplacian(temp2, alpha_T, temp1, tempx, tempy, tempz)
    output(:,:,:,4) = output(:,:,:,4) + (gamma*mu)/(Re*Pr*rho) * temp2

  END SUBROUTINE ns_rhs

  SUBROUTINE initialize_tau(alpha_tau_11, alpha_tau_12, alpha_tau_13, alpha_tau_22, &
        alpha_tau_23, alpha_tau_33, alpha_u, mu, temp, tempx, tempy, tempz)
    real, dimension(xNc, Ny, zNc) :: alpha_tau_11, alpha_tau_12, alpha_tau_13
    real, dimension(xNc, Ny, zNc) :: alpha_tau_22, alpha_tau_23, alpha_tau_33
    real, dimension(xNc, Ny, zNc, 3) :: alpha_u
    real, dimension(xNc, Ny, zNc) :: temp
    real, dimension(xNc, 2) :: tempx
    real, dimension(xNc, 2) :: tempy
    real, dimension(xNc, 2) :: tempz
    real :: mu

    ! tau_11 = mu * (4/3 dudx - 2/3 * (dvdy + dwdz))
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,1), temp, tempx, tempy, tempz)
    alpha_tau_11 = (4.0/3.0) * temp
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,2), temp, tempx, tempy, tempz)
    alpha_tau_11 = alpha_tau_11 - (2.0/3.0) * temp
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,3), temp, tempx, tempy, tempz)
    alpha_tau_11 = mu * (alpha_tau_11 - (2.0/3.0) * temp)
    call solve3d(alpha_tau_11, tempx(:,1), tempy(:,1), tempz(:,1))

    ! tau_12 = mu * (dudy + dvdx)
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,1), temp, tempx, tempy, tempz)
    alpha_tau_12 = temp
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,2), temp, tempx, tempy, tempz)
    alpha_tau_12 = mu * (alpha_tau_12 + temp)
    call solve3d(alpha_tau_12, tempx(:,1), tempy(:,1), tempz(:,1))

    ! tau_13 = mu * (dudz + dwdx)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,1), temp, tempx, tempy, tempz)
    alpha_tau_13 = temp
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,3), temp, tempx, tempy, tempz)
    alpha_tau_13 = mu * (alpha_tau_13 + temp)
    call solve3d(alpha_tau_13, tempx(:,1), tempy(:,1), tempz(:,1))

    ! tau_22 = mu * (4/3 dvdy - 2/3 * (dudx + dwdz))
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,2), temp, tempx, tempy, tempz)
    alpha_tau_22 = (4.0/3.0) * temp
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,1), temp, tempx, tempy, tempz)
    alpha_tau_22 = alpha_tau_22 - (2.0/3.0) * temp
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,3), temp, tempx, tempy, tempz)
    alpha_tau_22 = mu * (alpha_tau_22 - (2.0/3.0) * temp)
    call solve3d(alpha_tau_22, tempx(:,1), tempy(:,1), tempz(:,1))

    ! tau_23 = mu * (dvdz + dwdy)
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,2), temp, tempx, tempy, tempz)
    alpha_tau_23 = temp
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,3), temp, tempx, tempy, tempz)
    alpha_tau_23 = mu * (alpha_tau_23 + temp)
    call solve3d(alpha_tau_23, tempx(:,1), tempy(:,1), tempz(:,1))

    ! tau_33 = mu * (4/3 dwdz - 2/3 * (dudx + dvdy))
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha_u(:,:,:,3), temp, tempx, tempy, tempz)
    alpha_tau_33 = (4.0/3.0) * temp
    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha_u(:,:,:,1), temp, tempx, tempy, tempz)
    alpha_tau_33 = alpha_tau_33 - (2.0/3.0) * temp
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha_u(:,:,:,2), temp, tempx, tempy, tempz)
    alpha_tau_33 = mu * (alpha_tau_33 - (2.0/3.0) * temp)
    call solve3d(alpha_tau_33, tempx(:,1), tempy(:,1), tempz(:,1))

  END SUBROUTINE initialize_tau

  SUBROUTINE laplacian(output, alpha, temp, tempx, tempy, tempz)
    real, dimension(xNc, Ny, zNc) :: output, alpha, temp
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz

    call matvec(Bx(:, :, 2), dBy(:, :, 0), Bz(:, :, 0), alpha, temp, tempx, tempy, tempz)
    output = temp
    call matvec(Bx(:, :, 0), dBy(:, :, 2), Bz(:, :, 0), alpha, temp, tempx, tempy, tempz)
    output = output + temp
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 2), alpha, temp, tempx, tempy, tempz)
    output = output + temp
  END SUBROUTINE laplacian

  SUBROUTINE divergence(output, alpha1, alpha2, alpha3, temp, tempx, tempy, tempz)
    real, dimension(xNc, Ny, zNc) :: output, alpha1, alpha2, alpha3, temp
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz

    call matvec(Bx(:, :, 1), dBy(:, :, 0), Bz(:, :, 0), alpha1, temp, tempx, tempy, tempz)
    output = temp
    call matvec(Bx(:, :, 0), dBy(:, :, 1), Bz(:, :, 0), alpha2, temp, tempx, tempy, tempz)
    output = output + temp
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 1), alpha3, temp, tempx, tempy, tempz)
    output = output + temp
  END SUBROUTINE divergence

  SUBROUTINE alpha_3_to_1(alpha1, alpha2, alpha3, alpha4, temp, tempx, tempy, tempz)
    real, dimension(xNc, Ny, zNc) :: alpha1, alpha2, alpha3, alpha4, temp
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz

    ! B alpha2 + B alpha3 + B alpha4
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha2, temp, tempx, tempy, tempz)
    alpha1 = temp
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha3, temp, tempx, tempy, tempz)
    alpha1 = alpha1 * temp
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha4, temp, tempx, tempy, tempz)
    alpha1 = alpha1 * temp

    ! get combined alpha
    call solve3d(alpha1, tempx(:,1), tempy(:,1), tempz(:,1))
  END SUBROUTINE alpha_3_to_1

  SUBROUTINE alpha_2_to_1(alpha1, alpha2, alpha3, temp, tempx, tempy, tempz)
    real, dimension(xNc, Ny, zNc) :: alpha1, alpha2, alpha3, temp
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz

    ! B alpha1 * B alpha2
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha2, temp, tempx, tempy, tempz)
    alpha1 = temp
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha3, temp, tempx, tempy, tempz)
    alpha1 = alpha1 * temp

    ! get combined alpha
    call solve3d(alpha1, tempx(:,1), tempy(:,1), tempz(:,1))
  END SUBROUTINE alpha_2_to_1

  ! computes x = A^-1 b
  SUBROUTINE psolve3d(rhs, temp, tempx, tempy, tempz)
    real, dimension(xNc) :: tempx
    real, dimension(Ny) :: tempy
    real, dimension(zNc) :: tempz
    real, dimension(xNc, Ny, zNc) :: rhs, temp

    integer :: i, j, k

    temp = rhs

    call matvec(BTx, dBTy, BTz, temp, rhs, tempx, tempy, tempz)

    ! Solve in x dimension
    do k=1, zNc
       do j=1, Ny
          tempx = rhs(:, j, k)
          call pcirculant_xsolve(tempx)
          rhs(:, j, k) = tempx
       end do
    end do

    ! Solve in y dimension
    do k=1, zNc
       do i=1, xNc
          tempy = rhs(i, :, k)
          call DGBTRS('N', Ny, PKl, PKu, 1, PLUy, 2*PKl+PKu+1, ppivoty, tempy, Ny, pinfoy)
          rhs(i, :, k) = tempy
       end do
    end do

    ! Solve in z dimension
    do j=1, Ny
       do i=1, xNc
          tempz = rhs(i, j, :)
          call pcirculant_zsolve(tempz)
          rhs(i, j, :) = tempz
       end do
    end do
  END SUBROUTINE psolve3d

  ! computes x = A^-1 b
  SUBROUTINE solve3d(rhs, tempx, tempy, tempz)
    real, dimension(xNc) :: tempx
    real, dimension(Ny) :: tempy
    real, dimension(zNc) :: tempz
    real, dimension(xNc, Ny, zNc) :: rhs

    integer :: i, j, k

    ! Solve in x dimension
    do k=1, zNc
       do j=1, Ny
          tempx = rhs(:, j, k)
          call circulant_xsolve(tempx)
          rhs(:, j, k) = tempx
       end do
    end do

    ! Solve in y dimension
    do k=1, zNc
       do i=1, xNc
          tempy = rhs(i, :, k)
          call DGBTRS('N', Ny, BKl, BKu, 1, LUy, 2*BKl+BKu+1, pivoty, tempy, Ny, infoy)
          rhs(i, :, k) = tempy
       end do
    end do

    ! Solve in z dimension
    do j=1, Ny
       do i=1, xNc
          tempz = rhs(i, j, :)
          call circulant_zsolve(tempz)
          rhs(i, j, :) = tempz
       end do
    end do
  END SUBROUTINE solve3d

  ! computes Ax = b
  SUBROUTINE matvec(matX, matY, matZ, x, b, tempx, tempy, tempz)
    real, dimension(xNc, xNc) :: matX
    real, dimension(-BKl:BKu, Ny) :: matY
    real, dimension(zNc, zNc) :: matZ
    real, dimension(xNc, 2) :: tempx
    real, dimension(Ny, 2) :: tempy
    real, dimension(zNc, 2) :: tempz
    real, dimension(xNc, Ny, zNc) :: x, b

    integer :: i, j, k

    b = x

    ! x-dimension (periodic)
    do k=1, zNc
       do j=1, Ny
          tempx(:, 1) = b(:, j, k)
          call DGEMV('N', xNc, xNc, 1., matX, xNc, tempx(:, 1), 1, 0., tempx(:, 2), 1)
          b(:, j, k) = tempx(:, 2)
       end do
    end do

    ! y-dimension (dirichlet)
    do k=1, zNc
       do i=1, xNc
          tempy(:, 1) = b(i, :, k)
          call DGBMV('N', Ny, Ny, BKl, BKu, 1., matY, BKu+BKl+1, tempy(:, 1), 1, 0., tempy(:, 2), 1)
          b(i, :, k) = tempy(:, 2)
       end do
    end do

    ! z-dimension (periodic)
    do j=1, Ny
       do i=1, xNc
          tempz(:, 1) = b(i, j, :)
          call DGEMV('N', zNc, zNc, 1., matZ, zNc, tempz(:, 1), 1, 0., tempz(:, 2), 1)
          b(i, j, :) = tempz(:, 2)
       end do
    end do
  END SUBROUTINE matvec

END MODULE col_linalg
