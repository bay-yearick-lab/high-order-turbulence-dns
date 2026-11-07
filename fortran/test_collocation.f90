MODULE test_collocation
  USE prms
  USE bspline
  USE col_fft
  USE col_linalg
  USE col_operator
  USE rungekutta
  USE optimize

  IMPLICIT none

CONTAINS

  SUBROUTINE test_GCV_vector
    implicit none

    integer ( kind = 4 ), parameter :: n = 80 + 2

    real ( kind = 8 ) pr
    real ( kind = 8 ) praxis
    real ( kind = 8 ), external :: GCV_vector_u_praxis
    real ( kind = 8 ), external :: GCV_vector_v_praxis
    real ( kind = 8 ), external :: GCV_vector_w_praxis
    real ( kind = 8 ), external :: GCV_vector_rho_praxis
    real ( kind = 8 ), external :: GCV_vector_T_praxis
    real ( kind = 8 ) x(n)
    real ( kind = 8 ) tf, ts

    integer :: i, j, k
    real :: random

    call generate_operator

    call RANDOM_SEED

    ! create rhs
    do k=1, zNc
       do j=1, Ny
          do i=1, xNc
             call RANDOM_NUMBER(random)
             random = 2.0 * (random - 0.5) * 0.1
             global_u(i, j, k, 1) = 1.5 * (1.0 - (Yc(j) - 1.0)**2) * (1.0 + random)
             call RANDOM_NUMBER(random)
             random = 2.0 * (random - 0.5) * 0.1
             global_u(i, j, k, 2) = 1.5 * (1.0 - (Yc(j) - 1.0)**2) * (1.0 + random)
             call RANDOM_NUMBER(random)
             random = 2.0 * (random - 0.5) * 0.1
             global_u(i, j, k, 3) = 1.5 * (1.0 - (Yc(j) - 1.0)**2) * (1.0 + random)
             call RANDOM_NUMBER(random)
             random = 2.0 * (random - 0.5) * 0.1
             global_rho(i, j, k) = 1.5 * (1.0 - (Yc(j) - 1.0)**2) * (1.0 + random)
             call RANDOM_NUMBER(random)
             random = 2.0 * (random - 0.5) * 0.1
             global_T(i, j, k) = 1.5 * (1.0 - (Yc(j) - 1.0)**2) * (1.0 + random)
          end do
       end do
    end do

    ! Praxis
    x(:) = 1e-5
    write( *, '(A, I3)') '  Number of pts:', n
    write(*, '(A, E10.3)') '  Initial point (mean): ', sum(x)/n
    write ( *, '(a,g14.6)' ) '  Function value = ', GCV_vector_u_praxis ( x, n )
    call CPU_TIME(ts)
    pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_vector_u_praxis )
    call CPU_TIME(tf)
    print*, " Time taken: ", tf - ts
    write(*, '(A, E10.3)') '  Computed minimizer: ', sum(x)/n
    write ( *, '(a,g14.6)' ) '  Function value = ', pr
    print*, " "

    x(:) = 1e-5
    write( *, '(A, I3)') '  Number of pts:', n
    write(*, '(A, E10.3)') '  Initial point (mean): ', sum(x)/n
    write ( *, '(a,g14.6)' ) '  Function value = ', GCV_vector_v_praxis ( x, n )
    call CPU_TIME(ts)
    pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_vector_v_praxis )
    call CPU_TIME(tf)
    print*, " Time taken: ", tf - ts
    write(*, '(A, E10.3)') '  Computed minimizer: ', sum(x)/n
    write ( *, '(a,g14.6)' ) '  Function value = ', pr
    print*, " "

    x(:) = 1e-5
    write( *, '(A, I3)') '  Number of pts:', n
    write(*, '(A, E10.3)') '  Initial point (mean): ', sum(x)/n
    write ( *, '(a,g14.6)' ) '  Function value = ', GCV_vector_w_praxis ( x, n )
    call CPU_TIME(ts)
    pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_vector_w_praxis )
    call CPU_TIME(tf)
    print*, " Time taken: ", tf - ts
    write(*, '(A, E10.3)') '  Computed minimizer: ', sum(x)/n
    write ( *, '(a,g14.6)' ) '  Function value = ', pr
    print*, " "

    x(:) = 1e-5
    write( *, '(A, I3)') '  Number of pts:', n
    write(*, '(A, E10.3)') '  Initial point (mean): ', sum(x)/n
    write ( *, '(a,g14.6)' ) '  Function value = ', GCV_vector_rho_praxis ( x, n )
    call CPU_TIME(ts)
    pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_vector_rho_praxis )
    call CPU_TIME(tf)
    print*, " Time taken: ", tf - ts
    write(*, '(A, E10.3)') '  Computed minimizer: ', sum(x)/n
    write ( *, '(a,g14.6)' ) '  Function value = ', pr
    print*, " "

    x(:) = 1e-5
    write( *, '(A, I3)') '  Number of pts:', n
    write(*, '(A, E10.3)') '  Initial point (mean): ', sum(x)/n
    write ( *, '(a,g14.6)' ) '  Function value = ', GCV_vector_T_praxis ( x, n )
    call CPU_TIME(ts)
    pr = praxis ( 1.0e-5, 1.0e-7, n, 0, x, GCV_vector_T_praxis )
    call CPU_TIME(tf)
    print*, " Time taken: ", tf - ts
    write(*, '(A, E10.3)') '  Computed minimizer: ', sum(x)/n
    write ( *, '(a,g14.6)' ) '  Function value = ', pr
    print*, " "

    return
  END SUBROUTINE test_GCV_vector

  SUBROUTINE test_GCV
    implicit none

    integer ( kind = 4 ), parameter :: n = 3

    real ( kind = 8 ) pr
    real ( kind = 8 ) praxis
    real ( kind = 8 ), external :: GCV_u_praxis
    real ( kind = 8 ), external :: GCV_v_praxis
    real ( kind = 8 ), external :: GCV_w_praxis
    real ( kind = 8 ), external :: GCV_rho_praxis
    real ( kind = 8 ), external :: GCV_T_praxis
    real ( kind = 8 ) x(n)
    real ( kind = 8 ) tf, ts

    integer ( kind = 4 ) icount
    integer ( kind = 4 ) ifault
    integer ( kind = 4 ) kcount
    integer ( kind = 4 ) konvge
    integer ( kind = 4 ) numres
    real ( kind = 8 ) reqmin
    real ( kind = 8 ) start(n)
    real ( kind = 8 ) step(n)
    real ( kind = 8 ), external :: GCV_u_nelmin
    real ( kind = 8 ) xmin(n)
    real ( kind = 8 ) ynewlo

    integer :: i, j, k
    real :: random

    call generate_operator

    call RANDOM_SEED

    ! create rhs
    do k=1, zNc
       do j=1, Ny
          do i=1, xNc
             call RANDOM_NUMBER(random)
             random = 2.0 * (random - 0.5) * 0.1
             global_u(i, j, k, 1) = 1.5 * (1.0 - (Yc(j) - 1.0)**2) * (1.0 + random)
             call RANDOM_NUMBER(random)
             random = 2.0 * (random - 0.5) * 0.1
             global_u(i, j, k, 2) = 1.5 * (1.0 - (Yc(j) - 1.0)**2) * (1.0 + random)
             call RANDOM_NUMBER(random)
             random = 2.0 * (random - 0.5) * 0.1
             global_u(i, j, k, 3) = 1.5 * (1.0 - (Yc(j) - 1.0)**2) * (1.0 + random)
             call RANDOM_NUMBER(random)
             random = 2.0 * (random - 0.5) * 0.1
             global_rho(i, j, k) = 1.5 * (1.0 - (Yc(j) - 1.0)**2) * (1.0 + random)
             call RANDOM_NUMBER(random)
             random = 2.0 * (random - 0.5) * 0.1
             global_T(i, j, k) = 1.5 * (1.0 - (Yc(j) - 1.0)**2) * (1.0 + random)
          end do
       end do
    end do

    ! Praxis
    x(1) = 1e-5; x(2) = 1e-5; x(3) = 1e-5
    call r8vec_print ( n, x, '  Initial point:' )
    write ( *, '(a,g14.6)' ) '  Function value = ', GCV_u_praxis ( x, n )
    call CPU_TIME(ts)
    pr = praxis ( 1.0e-5, 1.0e-7, 3, 0, x, GCV_u_praxis )
    call CPU_TIME(tf)
    print*, "Time taken: ", tf - ts
    call r8vec_print ( n, x, '  Computed minimizer:' )
    write ( *, '(a,g14.6)' ) '  Function value = ', pr

    x(1) = 1e-5; x(2) = 1e-5; x(3) = 1e-5
    call r8vec_print ( n, x, '  Initial point:' )
    write ( *, '(a,g14.6)' ) '  Function value = ', GCV_v_praxis ( x, n )
    call CPU_TIME(ts)
    pr = praxis ( 1.0e-5, 1.0e-7, 3, 0, x, GCV_v_praxis )
    call CPU_TIME(tf)
    print*, "Time taken: ", tf - ts
    call r8vec_print ( n, x, '  Computed minimizer:' )
    write ( *, '(a,g14.6)' ) '  Function value = ', pr

    x(1) = 1e-5; x(2) = 1e-5; x(3) = 1e-5
    call r8vec_print ( n, x, '  Initial point:' )
    write ( *, '(a,g14.6)' ) '  Function value = ', GCV_w_praxis ( x, n )
    call CPU_TIME(ts)
    pr = praxis ( 1.0e-5, 1.0e-7, 3, 0, x, GCV_w_praxis )
    call CPU_TIME(tf)
    print*, "Time taken: ", tf - ts
    call r8vec_print ( n, x, '  Computed minimizer:' )
    write ( *, '(a,g14.6)' ) '  Function value = ', pr

    x(1) = 1e-5; x(2) = 1e-5; x(3) = 1e-5
    call r8vec_print ( n, x, '  Initial point:' )
    write ( *, '(a,g14.6)' ) '  Function value = ', GCV_rho_praxis ( x, n )
    call CPU_TIME(ts)
    pr = praxis ( 1.0e-5, 1.0e-7, 3, 0, x, GCV_rho_praxis )
    call CPU_TIME(tf)
    print*, "Time taken: ", tf - ts
    call r8vec_print ( n, x, '  Computed minimizer:' )
    write ( *, '(a,g14.6)' ) '  Function value = ', pr

    x(1) = 1e-5; x(2) = 1e-5; x(3) = 1e-5
    call r8vec_print ( n, x, '  Initial point:' )
    write ( *, '(a,g14.6)' ) '  Function value = ', GCV_T_praxis ( x, n )
    call CPU_TIME(ts)
    pr = praxis ( 1.0e-5, 1.0e-7, 3, 0, x, GCV_T_praxis )
    call CPU_TIME(tf)
    print*, "Time taken: ", tf - ts
    call r8vec_print ( n, x, '  Computed minimizer:' )
    write ( *, '(a,g14.6)' ) '  Function value = ', pr

    ! Nelder-Mead
    write ( *, '(a)' ) ' '
    start(1:n) = (/ 1.0e-5, 1.0e-5, 1.0e-5 /)
    reqmin = 1.0e-14
    step(1:n) = (/ 1.0e-5, 1.0e-5, 1.0e-5 /)
    konvge = 10
    kcount = 500
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Starting point X:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,g14.6)' ) start(i)
    end do
    ynewlo = GCV_u_nelmin ( start )
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo
    call CPU_TIME(ts)
    call nelmin ( GCV_u_nelmin, n, start, xmin, ynewlo, reqmin, step, &
      konvge, kcount, icount, numres, ifault )
    call CPU_TIME(tf)
    print*, "Time taken: ", tf - ts

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Return codeIFAULT = ', ifault
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Estimate of minimizing value X*:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,g24.10)' ) xmin(i)
    end do
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of iterations = ', icount
    write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

    return
  END SUBROUTINE test_GCV

  SUBROUTINE test_praxis
    implicit none

    integer ( kind = 4 ), parameter :: n = 2

    real ( kind = 8 ) h0
    real ( kind = 8 ) pr
    real ( kind = 8 ) praxis
    integer ( kind = 4 ) prin
    real ( kind = 8 ), external :: rosenbrock_praxis
    real ( kind = 8 ) t0
    real ( kind = 8 ) x(n)
    real ( kind = 8 ) finish, start

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROSENBROCK_TEST'
    write ( *, '(a)' ) '  The Rosenbrock function.'

    t0 = 1.0e-5
    h0 = 1.0
    prin = 0

    x(1) = -12.0; x(2) = 10.0

    call r8vec_print ( n, x, '  Initial point:' )

    write ( *, '(a,g14.6)' ) '  Function value = ', rosenbrock_praxis ( x, n )

    call CPU_TIME(start)
    pr = praxis ( t0, h0, n, prin, x, rosenbrock_praxis )
    call CPU_TIME(finish)
    print*, "Time taken: ", finish - start

    call r8vec_print ( n, x, '  Computed minimizer:' )

    write ( *, '(a,g14.6)' ) '  Function value = ', rosenbrock_praxis ( x, n )

    return
  END SUBROUTINE test_praxis

  SUBROUTINE test_3d_pspline
    real, allocatable, dimension(:, :) :: tempx, tempy, tempz
    real, allocatable, dimension(:, :, :) :: alpha, num, palpha, rhs
    integer :: i, j, k

    allocate(tempx(xNc, 2), tempy(Ny, 2), tempz(zNc, 2))
    allocate(palpha(xNc, Ny, zNc))
    allocate(rhs(xNc, Ny, zNc))
    allocate(num(xNc, Ny, zNc))
    allocate(alpha(xNc, Ny, zNc))

    call generate_operator

    ! create rhs
    do k=1, zNc
       do j=1, Ny
          do i=1, xNc
             rhs(i, j, k) = COS(2*Xc(i)) * COS(2*Yc(j)) * COS(2*Zc(k))
          end do
       end do
    end do

    alpha = rhs

    ! B-spline solve
    call solve3d(alpha, tempx(:, 1), tempy(:, 1), tempz(:, 1))

    ! P-spline solve
    lamda = 0.0
    call re_operator
    palpha = rhs
    call psolve3d(palpha, num, tempx(:, 1), tempy(:, 1), tempz(:, 1))
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), palpha, num, tempx, tempy, tempz)

    WRITE(*, '(A15, E10.3)') 'Solve error = ', MAXVAL(ABS(rhs - num))
    WRITE(*, '(A36, E10.3)') 'Difference between P and B solve = ', MAXVAL(ABS(alpha - palpha))

    call deallocate_operator
    deallocate(alpha, num, palpha, rhs)
    deallocate(tempx, tempy, tempz)
  END SUBROUTINE test_3d_pspline

  SUBROUTINE test_nelmin
    implicit none

    call test_nelmin1
    call test_nelmin2
    call test_nelmin3
  END SUBROUTINE test_nelmin

  SUBROUTINE test_nelmin3
    implicit none

    integer ( kind = 4 ), parameter :: n = 3

    integer ( kind = 4 ) i
    integer ( kind = 4 ) icount
    integer ( kind = 4 ) ifault
    integer ( kind = 4 ) kcount
    integer ( kind = 4 ) konvge
    integer ( kind = 4 ) numres
    real ( kind = 8 ) reqmin
    real ( kind = 8 ) start(n)
    real ( kind = 8 ) step(n)
    real ( kind = 8 ), external :: helical
    real ( kind = 8 ) xmin(n)
    real ( kind = 8 ) ynewlo

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HELICAL_TEST'
    write ( *, '(a)' ) '  Apply NELMIN to the HELICAL function.'

    start(1:n) = (/ 1.0, 1.0, 1.0 /)

    reqmin = 1.0e-14

    step(1:n) = (/ 1.0, 1.0, 1.0 /)

    konvge = 100
    kcount = 1000

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Starting point X:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,g14.6)' ) start(i)
    end do

    ynewlo = helical ( start )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

    call nelmin ( helical, n, start, xmin, ynewlo, reqmin, step, &
      konvge, kcount, icount, numres, ifault )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Return codeIFAULT = ', ifault
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Estimate of minimizing value X*:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,g24.10)' ) xmin(i)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of iterations = ', icount
    write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

    return
  END SUBROUTINE test_nelmin3

  SUBROUTINE test_nelmin2
    implicit none

    integer ( kind = 4 ), parameter :: n = 4

    integer ( kind = 4 ) i
    integer ( kind = 4 ) icount
    integer ( kind = 4 ) ifault
    integer ( kind = 4 ) kcount
    integer ( kind = 4 ) konvge
    integer ( kind = 4 ) numres
    real ( kind = 8 ) reqmin
    real ( kind = 8 ) start(n)
    real ( kind = 8 ) step(n)
    real ( kind = 8 ), external :: powell
    real ( kind = 8 ) xmin(n)
    real ( kind = 8 ) ynewlo

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POWELL_TEST'
    write ( *, '(a)' ) '  Apply NELMIN to the POWELL quartic function.'

    start(1:n) = (/ 3.0, -1.0, 0.0, 1.0 /)

    reqmin = 1.0e-14

    step(1:n) = (/ 1.0, 1.0, 1.0, 1.0 /)

    konvge = 400
    kcount = 1000

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Starting point X:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,g14.6)' ) start(i)
    end do

    ynewlo = powell ( start )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

    call nelmin ( powell, n, start, xmin, ynewlo, reqmin, step, &
      konvge, kcount, icount, numres, ifault )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Estimate of minimizing value X*:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,g24.10)' ) xmin(i)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of iterations = ', icount
    write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

    return
  END SUBROUTINE test_nelmin2

  SUBROUTINE test_nelmin1
    implicit none

    integer ( kind = 4 ), parameter :: n = 2

    integer ( kind = 4 ) i
    integer ( kind = 4 ) icount
    integer ( kind = 4 ) ifault
    integer ( kind = 4 ) kcount
    integer ( kind = 4 ) konvge
    integer ( kind = 4 ) numres
    real ( kind = 8 ) reqmin
    real ( kind = 8 ) start(n)
    real ( kind = 8 ) step(n)
    real ( kind = 8 ), external :: rosenbrock
    real ( kind = 8 ) xmin(n)
    real ( kind = 8 ) ynewlo
    real :: tstart, tfinish

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROSENBROCK_TEST'
    write ( *, '(a)' ) '  Apply NELMIN to the ROSENBROCK function.'

    start(1:n) = (/ -12.0, 10.0 /)

    reqmin = 1.0e-14

    step(1:n) = (/ 1.0, 1.0 /)

    konvge = 100
    kcount = 500

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Starting point X:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,g14.6)' ) start(i)
    end do

    ynewlo = rosenbrock ( start )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

    call CPU_TIME(tstart)
    call nelmin ( rosenbrock, n, start, xmin, ynewlo, reqmin, step, &
      konvge, kcount, icount, numres, ifault )
    call CPU_TIME(tfinish)
    print*, "Time taken: ", tfinish - tstart

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Estimate of minimizing value X*:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,g24.10)' ) xmin(i)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of iterations = ', icount
    write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

    return
  END SUBROUTINE test_nelmin1

  SUBROUTINE test_random
    real, dimension(Ny) :: random
    integer :: j

    call RANDOM_SEED

    do j=1, Ny
       call RANDOM_NUMBER(random(j))
       random(j) = 2.0 * (random(j) - 0.5)
       if (random(j) .gt. 1.0 .or. random(j) .lt. -1.0) then
          print*, "RANDOM NUMBER OUT OF [-1, 1]"
       end if
    end do
    print*, random
  END SUBROUTINE test_random

  SUBROUTINE test_ns_rhs
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

    call mms_init(output, sol, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    ! rhs
    call ns_rhs(output, gamma, Mach, Re, Pr, mu, rho, u, T, &
          alpha_rho, alpha_u, alpha_T, alpha_rho_u, alpha_rho_uT, &
          alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, &
          alpha_rho_vv, alpha_rho_vw, alpha_rho_ww, alpha_rho_T, &
          alpha_tau_11, alpha_tau_12, alpha_tau_13, &
          alpha_tau_22, alpha_tau_23, alpha_tau_33, &
          temp1, temp2, tempx, tempy, tempz)

    ! sol
    call mms_source(sol(:,:,:,0), sol(:,:,:,1:3), sol(:,:,:,4), 0.0, gamma, Mach, Re, Pr, mu)

    sol(:,:,:,0) = -(sol(:,:,:,0) + rho)
    sol(:,:,:,1:3) = -(sol(:,:,:,1:3) + u(:,:,:,1:3))
    sol(:,:,:,4) = -(sol(:,:,:,4) + T)

    WRITE(*, '(A26, E10.3)') 'Mass RHS error = ', &
        MAXVAL(ABS(sol(:,:,:,0) - output(:,:,:,0)))
    WRITE(*, '(A26, E10.3)') 'Momentum (i=1) RHS error = ', &
        MAXVAL(ABS(sol(:,:,:,1) - output(:,:,:,1)))
    WRITE(*, '(A26, E10.3)') 'Momentum (i=2) RHS error = ', &
        MAXVAL(ABS(sol(:,:,:,2) - output(:,:,:,2)))
    WRITE(*, '(A26, E10.3)') 'Momentum (i=3) RHS error = ', &
        MAXVAL(ABS(sol(:,:,:,3) - output(:,:,:,3)))
    WRITE(*, '(A26, E10.3)') 'Energy RHS error = ', &
        MAXVAL(ABS(sol(:,:,:,4) - output(:,:,:,4)))

    call deallocate_operator
    deallocate(alpha_rho_u, alpha_rho, alpha_u, sol)
    deallocate(temp1, temp2, tempx, tempy, tempz)

  END SUBROUTINE test_ns_rhs

  SUBROUTINE test_3d_bspline
    real, allocatable, dimension(:, :) :: tempx, tempy, tempz
    real, allocatable, dimension(:, :, :) :: alpha, num, rhs
    real, allocatable, dimension(:, :, :, :) :: sol
    integer :: i, j, k

    allocate(tempx(xNc, 2), tempy(Ny, 2), tempz(zNc, 2))
    allocate(rhs(xNc, Ny, zNc))
    allocate(sol(xNc, Ny, zNc, 0:2))
    allocate(num(xNc, Ny, zNc))
    allocate(alpha(xNc, Ny, zNc))

    call generate_operator

    ! create rhs
    do k=1, zNc
       do j=1, Ny
          do i=1, xNc
             sol(i, j, k, 0) = COS(2*Xc(i)) * COS(2*Yc(j)) * COS(2*Zc(k))
             sol(i, j, k, 1) = -2.0 * (COS(2*Xc(i)) * SIN(2*Yc(j)) * COS(2*Zc(k)) + &
                               COS(2*Yc(j)) * SIN(2.0 * (Xc(i) + Zc(k))))
             sol(i, j, k, 2) = -12.0 * COS(2*Xc(i)) * COS(2*Yc(j)) * COS(2*Zc(k))
          end do
       end do
    end do

    alpha = sol(:, :, :, 0)

    ! solve
    call solve3d(alpha, tempx(:, 1), tempy(:, 1), tempz(:, 1))

    ! matvec
    call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), alpha, num, tempx, tempy, tempz)
    WRITE(*, '(A18, E10.3)') 'Solve error = ', MAXVAL(ABS(sol(:,:,:,0) - num))

    call divergence(rhs, alpha, alpha, alpha, num, tempx, tempy, tempz)
    WRITE(*, '(A18, E10.3)') 'Divergence error = ', MAXVAL(ABS(sol(:,:,:,1) - rhs))

    call laplacian(rhs, alpha, num, tempx, tempy, tempz)
    WRITE(*, '(A18, E10.3)') 'Laplacian error = ', MAXVAL(ABS(sol(:,:,:,2) - rhs))

    call deallocate_operator
    deallocate(alpha, num, rhs, sol)
    deallocate(tempx, tempy, tempz)

  END SUBROUTINE test_3d_bspline

  SUBROUTINE test_1d_bspline
    real, allocatable, dimension(:) :: alphaz, numz
    real, allocatable, dimension(:) :: alphay, numy
    real, allocatable, dimension(:) :: alphax, numx

    allocate(alphax(xNc), alphay(Ny), alphaz(zNc))
    allocate(numx(xNc), numy(Ny), numz(zNc))

    call generate_operator

    ! set rhs
    alphax = SIN(2*Xc)

    call circulant_xsolve(alphax)

    ! error for zero derivative
    call DGEMV('N', xNc, xNc, 1., Bx(:, :, 0), xNc, alphax, 1, 0., numx, 1)
    WRITE(*, '(A36, E10.3)') 'X-direction 0th derivative error = ', MAXVAL(ABS(SIN(2*Xc) - numx))

    ! error for first derivative
    call DGEMV('N', xNc, xNc, 1., Bx(:, :, 1), xNc, alphax, 1, 0., numx, 1)
    WRITE(*, '(A36, E10.3)') 'X-direction 1st derivative error = ', MAXVAL(ABS(2*COS(2*Xc) - numx))

    ! error for second derivative
    call DGEMV('N', xNc, xNc, 1., Bx(:, :, 2), xNc, alphax, 1, 0., numx, 1)
    WRITE(*, '(A36, E10.3)') 'X-direction 2nd derivative error = ', MAXVAL(ABS(-4*SIN(2*Xc) - numx))

    ! set rhs
    alphay = SIN(2*Yc)

    ! solve
    call DGBTRS('N', Ny, BKl, BKu, 1, LUy, 2*BKl+BKu+1, pivoty, alphay, Ny, infoy)

    ! error for zero derivative
    call DGBMV('N', Ny, Ny, BKl, BKu, 1., dBy(:, :, 0), BKu+BKl+1, alphay, 1, 0., numy, 1)
    WRITE(*, '(A36, E10.3)') 'Y-direction 0th derivative error = ', MAXVAL(ABS(SIN(2*Yc) - numy))

    ! error for first derivative
    call DGBMV('N', Ny, Ny, BKl, BKu, 1., dBy(:, :, 1), BKu+BKl+1, alphay, 1, 0., numy, 1)
    WRITE(*, '(A36, E10.3)') 'Y-direction 1st derivative error = ', MAXVAL(ABS(2*COS(2*Yc) - numy))

    ! error for second derivative
    call DGBMV('N', Ny, Ny, BKl, BKu, 1., dBy(:, :, 2), BKu+BKl+1, alphay, 1, 0., numy, 1)
    WRITE(*, '(A36, E10.3)') 'Y-direction 2nd derivative error = ', MAXVAL(ABS(-4*SIN(2*Yc) - numy))

    ! set rhs
    alphaz = SIN(2*Zc)

    call circulant_zsolve(alphaz)

    ! error for zero derivative
    call DGEMV('N', zNc, zNc, 1., Bz(:, :, 0), zNc, alphaz, 1, 0., numz, 1)
    WRITE(*, '(A36, E10.3)') 'Z-direction 0th derivative error = ', MAXVAL(ABS(SIN(2*Zc) - numz))

    ! error for first derivative
    call DGEMV('N', zNc, zNc, 1., Bz(:, :, 1), zNc, alphaz, 1, 0., numz, 1)
    WRITE(*, '(A36, E10.3)') 'Z-direction 1st derivative error = ', MAXVAL(ABS(2*COS(2*Zc) - numz))

    ! error for second derivative
    call DGEMV('N', zNc, zNc, 1., Bz(:, :, 2), zNc, alphaz, 1, 0., numz, 1)
    WRITE(*, '(A36, E10.3)') 'Z-direction 2nd derivative error = ', MAXVAL(ABS(-4*SIN(2*Zc) - numz))

    call deallocate_operator
    deallocate(numx, numy, numz)
    deallocate(alphax, alphay, alphaz)
  END SUBROUTINE test_1d_bspline

END MODULE test_collocation
