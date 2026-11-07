MODULE prms

  IMPLICIT none

  integer   ::  Nx              ! mesh size in x
  integer   ::  Ny              ! mesh size in y (B-spline knot points)
  integer   ::  Nz              ! mesh size in z

  integer   ::  Ntime           ! # timesteps
  real      ::  Lx              ! box in x
  real      ::  Ly              ! box in y
  real      ::  Lz              ! box in z

  integer   ::  xNc             ! B-spline collation points x
  integer   ::  yNc             ! B-spline collation points y
  integer   ::  zNc             ! B-spline collation points z
  integer   ::  d               ! B-spline order

  real      ::  dt              ! timestep
  integer   ::  out_q           ! output Q this often
  integer   ::  out_prof        ! output profile
  integer   ::  out_vort        ! output profile
  integer   ::  out_cfl         ! output cfl
  real      ::  epsilon1         ! wave ampl

  real      ::  Uo              ! reference velocity
  integer   ::  which_ic        ! select initial condition

  real      ::  Wxr             ! boundary zone widths
  real      ::  Wxl             ! boundary zone widths
  real      ::  Wy              !
  real      ::  Wz              !
  real      ::  sigmax          ! boundary zone max sig

  integer   ::  which_x_mesh
  real      ::  xmesh_sig
  real      ::  xmesh_sigl
  real      ::  xmesh_A
  real      ::  xmesh_Al

  integer   ::  which_y_mesh
  real      ::  ymesh_sig
  real      ::  ymesh_A

  integer   ::  which_z_mesh
  real      ::  zmesh_sig
  real      ::  zmesh_A


  real      ::  Xph             ! size of x physical domain
  real      ::  Yph             ! size of y physical domain
  real      ::  Zph             ! size of z physical domain

  logical   ::  periodic        ! periodic in x


  real      ::  Re              ! Reynolds number
  real      ::  Pr              ! Prandtl number

  integer          ::  itime_start  ! first timestep number
  real             ::  time_start   ! code time

  real             ::  mu       ! viscosity (* computed *)
  real             ::  kappa    ! thermal conductivity (* computed *)
  real             ::  Pi       ! cake
  real             ::  f4o3     ! four thirds
  real             ::  f1o3     ! one third
  real             ::  f2o3     ! two thirds
  real, parameter  ::  gamma=1.4  ! ratio of specific heats

CONTAINS

  SUBROUTINE initprms

    f4o3 = 4./3.
    f1o3 = 1./3.
    f2o3 = 2./3.
    call readprms

  END SUBROUTINE initprms

  SUBROUTINE readprms

    open(1,file='mix.in')

    read(1,*) Nx
    read(1,*) Ny
    read(1,*) Nz
    read(1,*) Ntime

    read(1,*) Lx
    read(1,*) Ly
    read(1,*) Lz
    read(1,*) dt

    read(1,*) d; xNc = Nx-d; yNc = Ny-d; zNc = Nz-d

    read(1,*) out_q
    read(1,*) out_prof
    read(1,*) out_vort
    read(1,*) out_cfl
    read(1,*) epsilon1
    read(1,*) Re ; mu = 1./Re
    read(1,*) Pr ; kappa = 1./Re/Pr

    read(1,*) Uo
    read(1,*) which_ic

    read(1,*) Wxl
    read(1,*) Wxr
    read(1,*) Wy
    read(1,*) Wz
    read(1,*) sigmax

    read(1,*)  which_x_mesh
    read(1,*)  xmesh_sig
    read(1,*)  xmesh_sigl
    read(1,*)  xmesh_A
    read(1,*)  xmesh_Al

    read(1,*)  which_y_mesh
    read(1,*)  ymesh_sig
    read(1,*)  ymesh_A

    read(1,*)  which_z_mesh
    read(1,*)  zmesh_sig
    read(1,*)  zmesh_A

    read(1,*)  Xph
    read(1,*)  Yph
    read(1,*)  Zph

    read(1,*)  periodic

    close(1)

  END SUBROUTINE readprms

END MODULE prms
