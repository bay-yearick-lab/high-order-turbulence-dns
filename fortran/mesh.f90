MODULE mesh

  USE prms
  USE bspline

  IMPLICIT none

  real,allocatable,dimension(:)       :: X,Y,Z
  real,allocatable,dimension(:,:)     :: sx,sy,sz

CONTAINS

  SUBROUTINE allocatemesh

    allocate(X(Nx),sx(Nx,3))
    allocate(Y(Ny),sy(Ny,3))
    allocate(Z(Nz),sz(Nz,3))

  END SUBROUTINE allocatemesh

  SUBROUTINE initmesh

    real, dimension(Nx,3) :: dx
    real, dimension(Ny,3) :: dy
    real, dimension(Nz,3) :: dz

    call allocatemesh

    select case(which_x_mesh)
    case(0)
       call uniformP(Nx,X,dx,Lx)
    case(1)
       call uniform(Nx,X,dx,Lx)
    case(2)
       call tanhmesh(Nx,X,dx,xmesh_sig,xmesh_sigl,xmesh_A,xmesh_Al,Lx)
    end select
    sx(:,1) = 1./dx(:,1)
    sx(:,2) = 1./dx(:,1)**2
    sx(:,3) = -dx(:,2)/dx(:,1)**3

!    call write1mesh(Nx,X,'x')

    select case(which_y_mesh)
    case(1)
       call uniform(Ny,Y,dy,Ly)
    case(2)
       call tanhmesh(Ny,Y,dy,ymesh_sig,ymesh_sig,ymesh_A,ymesh_A,Ly)
    end select
    Y = Y - Ly/2.
    sy(:,1) = 1./dy(:,1)
    sy(:,2) = 1./dy(:,1)**2
    sy(:,3) = -dy(:,2)/dy(:,1)**3

    Y = Yc
!    call write1mesh(Ny,Y,'y')

    select case(which_z_mesh)
    case(0)
       call uniformP(Nz,Z,dz,Lz)
    case(1)
       call uniform(Nz,Z,dz,Lz)
    case(2)
       call tanhmesh(Nz,Z,dz,zmesh_sig,zmesh_sig,zmesh_A,zmesh_A,Lz)
    end select
    Z = Z - Lz/2.
    sz(:,1) = 1./dz(:,1)
    sz(:,2) = 1./dz(:,1)**2
    sz(:,3) = -dz(:,2)/dz(:,1)**3

!    call write1mesh(Nz,Z,'z')


  END SUBROUTINE initmesh


  SUBROUTINE write1mesh(Na,A,scriptor)
    integer               :: Na
    real, dimension(Nx)   :: A
    character(1)          :: scriptor

    integer :: i

    open(1,file='mix.grid.'//scriptor)
    do i = 1,Na
       write(1,"(I5,4F20.10)")i,A(i),sx(i,:)
    end do
    close(1)

  END SUBROUTINE write1mesh

  SUBROUTINE uniformP(Na,A,da,Amax)
    integer               :: Na
    real, dimension(Na)   :: A
    real, dimension(Na,2) :: da
    real                  :: Amax

    integer :: i
    real    :: s

    ! ... s: 0 - > 1
    ! ... dx1 = dx/ds
    ! ... dx2 = d^2x/ds^2
    do i = 1,Na
       s       = real(i-1)/real(Na)
       A(i)    = Amax * s
       da(i,1) = Amax
       da(i,2) = 0.
    end do

  END SUBROUTINE uniformP

  SUBROUTINE uniform(Na,A,da,Amax)
    integer               :: Na
    real, dimension(Na)   :: A
    real, dimension(Na,2) :: da
    real                  :: Amax

    integer :: i
    real    :: s

    ! ... s: 0 - > 1
    ! ... dx1 = dx/ds
    ! ... dx2 = d^2x/ds^2
    do i = 1,Na
       s       = real(i-1)/real(Na-1)
       A(i)    = Amax * s
       da(i,1) = Amax
       da(i,2) = 0.
    end do

  END SUBROUTINE uniform

  SUBROUTINE tanhmesh(Na,A,da,sig1,sig2,Am1,Am2,Amax)
    integer             :: Na
    real, dimension(Na) :: A
    real, dimension(Na,2) ::da
    real :: sig1,Am1,sig2,Am2,Amax

    real :: s
    real :: go
    integer :: i

! ... s: 0 - > 1
! ... da1 = da/ds
! ... da2 = d^2a/ds^2
    go = (1.+Am1/2.+Am2/2.) - Am1/2./sig1*LOG(COSH(-sig1)) &
                         - Am2/2./sig2*LOG(COSH(sig2))
    do i = 1,Na
       s = REAL(i-1)/REAL(Na-1)

       A(i)    = s*(Amax/go*(1.+Am1/2.+Am2/2.)) &
            + Am1*Amax/2./go/sig1*(LOG(COSH(sig1*(s-1.)))-LOG(COSH(-sig1))) &
            - Am2/2.*Amax/go/sig2*LOG(COSH(sig2*s))

       da(i,1) = Amax/go*(1. + Am1/2.*(1.+TANH(sig1*(s-1.))) &
                              +Am2/2.*(1.-TANH(sig2*s))  )
       da(i,2) = Amax*Am1/2./go*sig1/COSH(sig1*(s-1.))**2 &
                -Amax*Am2/2./go*sig2/COSH(sig2*s)**2
    end do

  END SUBROUTINE tanhmesh

END MODULE mesh
