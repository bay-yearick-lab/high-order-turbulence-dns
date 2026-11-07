FUNCTION rosenbrock_praxis(x, n)

  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  real ( kind = 8 ) rosenbrock_praxis
  real ( kind = 8 ) x(n)

  fx1 = x(2) - x(1)**2
  fx2 = 1.0 - x(1)

  fx = 100.0 * fx1 * fx1 &
     +         fx2 * fx2

  rosenbrock_praxis = fx
  return
END FUNCTION rosenbrock_praxis

FUNCTION GCV_vector_u_praxis(x, n)
  use col_linalg
  use col_operator
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) GCV_vector_u_praxis

  vector_lamda = x
  call vector_re_operator

  global_temp1 = global_u(:,:,:,1)

  call psolve3d(global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
      global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  GCV_vector_u_praxis = SUM(global_u(:,:,:,1) - global_temp2)**2
  return
END FUNCTION GCV_vector_u_praxis

FUNCTION GCV_vector_v_praxis(x, n)
  use col_linalg
  use col_operator
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) GCV_vector_v_praxis

  vector_lamda = x
  call vector_re_operator

  global_temp1 = global_u(:,:,:,2)

  call psolve3d(global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
      global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  GCV_vector_v_praxis = SUM(global_u(:,:,:,2) - global_temp2)**2
  return
END FUNCTION GCV_vector_v_praxis

FUNCTION GCV_vector_w_praxis(x, n)
  use col_linalg
  use col_operator
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) GCV_vector_w_praxis

  vector_lamda = x
  call vector_re_operator

  global_temp1 = global_u(:,:,:,3)

  call psolve3d(global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
      global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  GCV_vector_w_praxis = SUM(global_u(:,:,:,3) - global_temp2)**2
  return
END FUNCTION GCV_vector_w_praxis

FUNCTION GCV_vector_rho_praxis(x, n)
  use col_linalg
  use col_operator
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) GCV_vector_rho_praxis

  vector_lamda = x
  call vector_re_operator

  global_temp1 = global_rho

  call psolve3d(global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
      global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  GCV_vector_rho_praxis = SUM(global_rho - global_temp2)**2
  return
END FUNCTION GCV_vector_rho_praxis

FUNCTION GCV_vector_T_praxis(x, n)
  use col_linalg
  use col_operator
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) GCV_vector_T_praxis

  vector_lamda = x
  call vector_re_operator

  global_temp1 = global_T

  call psolve3d(global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
      global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  GCV_vector_T_praxis = SUM(global_T - global_temp2)**2
  return
END FUNCTION GCV_vector_T_praxis

FUNCTION GCV_u_praxis(x, n)
  use col_linalg
  use col_operator
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) GCV_u_praxis

  lamda = x
  call re_operator

  global_temp1 = global_u(:,:,:,1)

  call psolve3d(global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
      global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  GCV_u_praxis = SUM(global_u(:,:,:,1) - global_temp2)**2
  return
END FUNCTION GCV_u_praxis

FUNCTION GCV_v_praxis(x, n)
  use col_linalg
  use col_operator
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) GCV_v_praxis

  lamda = x
  call re_operator

  global_temp1 = global_u(:,:,:,2)

  call psolve3d(global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
      global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  GCV_v_praxis = SUM(global_u(:,:,:,2) - global_temp2)**2
  return
END FUNCTION GCV_v_praxis

FUNCTION GCV_w_praxis(x, n)
  use col_linalg
  use col_operator
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) GCV_w_praxis

  lamda = x
  call re_operator

  global_temp1 = global_u(:,:,:,3)

  call psolve3d(global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
      global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  GCV_w_praxis = SUM(global_u(:,:,:,3) - global_temp2)**2
  return
END FUNCTION GCV_w_praxis

FUNCTION GCV_rho_praxis(x, n)
  use col_linalg
  use col_operator
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) GCV_rho_praxis

  lamda = x
  call re_operator

  global_temp1 = global_rho

  call psolve3d(global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
      global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  GCV_rho_praxis = SUM(global_rho - global_temp2)**2
  return
END FUNCTION GCV_rho_praxis

FUNCTION GCV_T_praxis(x, n)
  use col_linalg
  use col_operator
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) GCV_T_praxis

  lamda = x
  call re_operator

  global_temp1 = global_T

  call psolve3d(global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
      global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  GCV_T_praxis = SUM(global_T - global_temp2)**2
  return
END FUNCTION GCV_T_praxis

FUNCTION GCV_u_nelmin(x)
  use col_linalg
  use col_operator
  implicit none

  real ( kind = 8 ) GCV_u_nelmin
  real ( kind = 8 ) x(3)

  lamda = x
  call re_operator

  global_temp1 = global_u(:,:,:,1)

  call psolve3d(global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  call matvec(Bx(:, :, 0), dBy(:, :, 0), Bz(:, :, 0), &
      global_temp1, global_temp2, global_tempx, global_tempy, global_tempz)

  GCV_u_nelmin = SUM(global_u(:,:,:,1) - global_temp2)**2

  return
END FUNCTION GCV_u_nelmin

FUNCTION rosenbrock(x)

  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  real ( kind = 8 ) rosenbrock
  real ( kind = 8 ) x(2)

  fx1 = x(2) - x(1)**2
  fx2 = 1.0 - x(1)

  fx = 100.0 * fx1 * fx1 &
     +         fx2 * fx2

  rosenbrock = fx

  return
END FUNCTION rosenbrock

FUNCTION powell(x)

  implicit none

  real ( kind = 8 ) powell
  real ( kind = 8 ) x(4)

  powell = (x(1) + 10.0 * x(2))**2 + &
           5.0 * (x(3) - x(4))**2 + &
           (x(2) -2.0 * x(3))**4 + &
           10.0 * (x(1) - x(4))**4

  return
END FUNCTION powell

FUNCTION helical(x)
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  real ( kind = 8 ) helical
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(3)

  if ( 0.0D+00 < x(1) ) then
    theta = atan2 ( x(2), x(1) ) / 2.0D+00 / r8_pi
  else if ( x(1) < 0.0D+00 ) then
    theta = 0.5D+00 + atan2 ( x(2), x(1) ) / 2.0D+00 / r8_pi
  else if ( x(1) == 0.0D+00 ) then
    theta = 0.25D+00
  end if

  fx1 = x(3) - 10.0 * theta
  fx2 = SQRT(x(1)**2 + x(2)**2) - 1.0

  helical = 100.0 * (fx1**2 + fx2**2) + x(3)**2

  return
END FUNCTION helical


