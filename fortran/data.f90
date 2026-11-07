MODULE data

  USE prms 

  IMPLICIT none 

  real,allocatable,dimension(:,:,:,:) :: Q

CONTAINS

  SUBROUTINE initdata

    allocate(Q(Nx,Ny,Nz,5))


  END SUBROUTINE initdata

! convert from primative to conservative variables
  SUBROUTINE CtoU(C,U)
    real, dimension(Nx,Ny,Nz,5) :: C
    real, dimension(Nx,Ny,Nz,3) :: U

    U(:,:,:,1) = C(:,:,:,1)/C(:,:,:,4)    ! rho u/rho -> u
    U(:,:,:,2) = C(:,:,:,2)/C(:,:,:,4)    ! rho v/rho -> v
    U(:,:,:,3) = C(:,:,:,3)/C(:,:,:,4)    ! rho w/rho -> w

  END SUBROUTINE CtoU

! convert from primative to conservative variables
! NOTE:  P and C can be same in memory
  SUBROUTINE Prim2Cons(P,C)
    real,dimension(Nx,Ny,Nz,5) :: C  ! rho u ; rho v ; rho w ;  rho  ; e
    real,dimension(Nx,Ny,Nz,5) :: P  ! u   ; v   ; w   ; rho ; p

    C(:,:,:,1) = P(:,:,:,1)*P(:,:,:,4)    ! u rho
    C(:,:,:,2) = P(:,:,:,2)*P(:,:,:,4)    ! v rho
    C(:,:,:,3) = P(:,:,:,3)*P(:,:,:,4)    ! w rho
    C(:,:,:,4) = P(:,:,:,4)            ! rho
    C(:,:,:,5) = P(:,:,:,5)/(gamma-1.) &  ! p/(gamma-1) + 0.5*rho*(u^2+v^2+w^2)
         + 0.5*(C(:,:,:,1)**2+C(:,:,:,2)**2+C(:,:,:,3)**2)/C(:,:,:,4)

  END SUBROUTINE Prim2Cons

! convert from conservative to primative variables
! NOTE:  P and C can be same in memory
  SUBROUTINE Cons2Prim(C,P)
    real,dimension(Nx,Ny,Nz,5) :: C  ! u/V ; v/V ; w/V ; V   ; e
    real,dimension(Nx,Ny,Nz,5) :: P  ! u   ; v   ; w   ; rho ; e

    P(:,:,:,4) = 1./C(:,:,:,4)              ! rho
    P(:,:,:,1) = C(:,:,:,1)/P(:,:,:,4)      ! u
    P(:,:,:,2) = C(:,:,:,2)/P(:,:,:,4)      ! v
    P(:,:,:,3) = C(:,:,:,3)/P(:,:,:,4)      ! w
    P(:,:,:,5) = (gamma-1.)*(C(:,:,:,5) - & ! p
         0.5*(P(:,:,:,1)**2+P(:,:,:,2)**2+P(:,:,:,3)**2)*P(:,:,:,4))

  END SUBROUTINE Cons2Prim

! compute pressure
  SUBROUTINE pressure(C,P)
    real,dimension(Nx,Ny,Nz,5) :: C  ! u/V ; v/V ; w/V ; V   ; e
    real,dimension(Nx,Ny,Nz)   :: P  ! pressure
    
! p = (gamma-1)*(e - 0.5 rho (u^2 + v^2 + w^2)

    P = (gamma-1.)*(C(:,:,:,5) - & ! p
         0.5*(C(:,:,:,1)**2+C(:,:,:,2)**2+C(:,:,:,3)**2)/C(:,:,:,4))

  END SUBROUTINE Pressure

END MODULE data
