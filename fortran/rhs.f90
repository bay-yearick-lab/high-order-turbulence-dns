MODULE rhside

  use prms
  use data
  use derivatives
  use bc

  IMPLICIT none

CONTAINS


  SUBROUTINE rhs(Q,dQ)

    real,dimension(Nx,Ny,Nz,5)   :: Q,dQ
    real,dimension(Nx,Ny,Nz)     :: F,P,dF
    real,dimension(Nx,Ny,Nz,3)   :: U
    real,dimension(Nx,Ny,Nz,3,3) :: A
    real,dimension(Nx,Ny,Nz,3)   :: V

    call pressure(Q,P)
    call CtoU(Q,U)

!  d rho u      d rho u u     d p
!  -------  = - ---------  -  ---
!    d t           d x        d x
    F = - Q(:,:,:,1)*U(:,:,:,1) - P
    call derivativeX(F,dQ(1,1,1,1),1)

!  d rho v      d rho u v
!  ------- = -  ---------
!    d t           d x
    F = - Q(:,:,:,1)*U(:,:,:,2)
    call derivativeX(F,dQ(1,1,1,2),1)

!  d rho w      d rho u w
!  ------- = -  ---------
!    d t           d x
    F = - Q(:,:,:,1)*U(:,:,:,3)
    call derivativeX(F,dQ(1,1,1,3),1)


!   d rho      d rho u
! - -----   =  -------
!    d t         d x
    call derivativeX(Q(1,1,1,1),dQ(1,1,1,4),1)

!  d e      d u(e + p)
!  ---  = - ---------
!  d t         d x
    F = -U(:,:,:,1)*(Q(:,:,:,5) + P)
    call derivativeX(F,dQ(1,1,1,5),1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  d rho u       d rho u v
!  ------- += -  ---------
!    d t            d y
    F = - Q(:,:,:,2)*U(:,:,:,1)
    call derivativeY(F,dF)
    dQ(:,:,:,1) = dQ(:,:,:,1) + dF

!  d rho v       d rho v v   d P
!  ------- += -  --------- - ---
!    d t            d y      d y

    F = - Q(:,:,:,2)*U(:,:,:,2) - P
    call derivativeY(F,dF)
    dQ(:,:,:,2) = dQ(:,:,:,2) + dF

!  d rho w      d rho v w
!  ------- += - ---------
!    d t           d y
    F = - Q(:,:,:,2)*U(:,:,:,3)
    call derivativeY(F,dF)
    dQ(:,:,:,3) = dQ(:,:,:,3) + dF

!   d rho        d rho v
!   ------  += - -------
!    d t           d y
    call derivativeY(Q(1,1,1,2),dF)
    dQ(:,:,:,4) = - dQ(:,:,:,4) - dF

!  d e      d v(e + p)
!  --- += - ----------
!  d t         d y
    F = -U(:,:,:,2)*(Q(:,:,:,5) + P)
    call derivativeY(F,dF)
    dQ(:,:,:,5) = dQ(:,:,:,5) + dF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  d rho u      d rho u w
!  ------- += - ---------
!    d t           d z
    F = - Q(:,:,:,3)*U(:,:,:,1)
    call derivativeZ(F,dF,1)
    dQ(:,:,:,1) = dQ(:,:,:,1) + dF

!  d rho v      d rho v w
!  ------- += - ---------
!    d t           d z
    F = - Q(:,:,:,3)*U(:,:,:,2)
    call derivativeZ(F,dF,1)
    dQ(:,:,:,2) = dQ(:,:,:,2) + dF


!  d rho w       d rho w w   d P
!  ------- += -  --------- - ---
!    d t            d z      d z
    F = - Q(:,:,:,3)*U(:,:,:,3) - P
    call derivativeZ(F,dF,1)
    dQ(:,:,:,3) = dQ(:,:,:,3) + dF

!   d rho        d rho w
!   ------  += - -------
!    d t           d z
    call derivativeZ(Q(1,1,1,3),dF,1)
    dQ(:,:,:,4) = dQ(:,:,:,4) - dF

!  d e      d w(e + p)
!  --- += - ----------
!  d t         d z
    F = -U(:,:,:,3)*(Q(:,:,:,5) + P)
    call derivativeZ(F,dF,1)
    dQ(:,:,:,5) = dQ(:,:,:,5) + dF


!   d e       d^2 T
!   ---  +=  -------
!   d t      dxi dxi                    !       p        p  gamma
    P = P/Q(:,:,:,4)*gamma/(gamma-1.)   ! T = -----  = ------------
                                        !     R rho    rho (gamma-1)
    call Lap(P,dF)
    dQ(:,:,:,5) = dQ(:,:,:,5) + kappa*dF

!
!  d e       [  (d u)^2     (d v)^2      (d w)^2
!  --- += mu [2 (---)   + 2 (---)    + 2 (---)
!  d t       [  (d x)       (d y)        (d z)
!
!               (d v   d u)^2   (d w   d v)^2   (d u   d w)^2
!             + (--- + ---)   + (--- + ---)   + (--- + ---)
!               (d x   d y)     (d y   d z)     (d z   d x)
!
!             2 (d u   d v   d w)^2 ]
!           - - (--- + --- + ---)   ]
!             3 (d x   d y   d z)   ]
!
!

    call getA(U,A,1)
    dQ(:,:,:,5) = dQ(:,:,:,5) + mu*( &
                    2.*(A(:,:,:,1,1)**2 +A(:,:,:,2,2)**2+A(:,:,:,3,3)**2) &
                    +  (A(:,:,:,2,1) + A(:,:,:,1,2))**2 &
                    +  (A(:,:,:,1,3) + A(:,:,:,3,1))**2 &
                    +  (A(:,:,:,2,3) + A(:,:,:,3,2))**2 &
                    -  f2o3*(A(:,:,:,1,1)+A(:,:,:,2,2)+A(:,:,:,3,3))**2 )

!       mu    ( d^2 v     d^2 w )
!  Vx = --- * ( -----  +  ----- )
!        3    ( dx dy     dx dz )
!
    F = mu*f1o3*(A(:,:,:,2,2) + A(:,:,:,3,3))
    call derivativeX(F,V(1,1,1,1),1)

!       mu    ( d^2 u     d^2 w )
!  Vy = --- * ( -----  +  ----- )
!        3    ( dx dy     dy dz )
!
    F = mu*f1o3*(A(:,:,:,1,1) + A(:,:,:,3,3))
    call derivativeY(F,V(1,1,1,2))

!       mu    ( d^2 u     d^2 v )
!  Vz = --- * ( -----  +  ----- )
!        3    ( dx dz     dy dz )
!
    F = mu*f1o3*(A(:,:,:,1,1) + A(:,:,:,2,2))
    call derivativeZ(F,V(1,1,1,3),1)


!                 ( 4 d^2 u   d^2 u   d^2 u )
!  Vx = Vx + mu * ( - ----- + ----- + ----- )
!                 ( 3 d x^2   d y^2   d z^2 )
!
!                 ( d^2 v   4 d^2 v   d^2 v )
!  Vy = Vy + mu * ( ----- + - ----- + ----- )
!                 ( d x^2   3 d y^2   d z^2 )
!
!                 ( d^2 w   d^2 w   4 d^2 w )
!  Vz = Vz + mu * ( ----- + ----- + - ----- )
!                 ( d x^2   d y^2   3 d z^2 )
!
!  d rho u            d rho v           d rho w
!  ------- +=  Vx  ;  ------- += Vy  ;  ------- += Vz
!    d t                d t               d t
!
!  d e
!  ---  += u Vx + v Vy + w Vz
!  d t
!
    call getA(U,A,2)
    V(:,:,:,1) =  V(:,:,:,1)  + mu*(f4o3*A(:,:,:,1,1)+A(:,:,:,1,2)+A(:,:,:,1,3))
    V(:,:,:,2) =  V(:,:,:,2)  + mu*(A(:,:,:,2,1)+f4o3*A(:,:,:,2,2)+A(:,:,:,2,3))
    V(:,:,:,3) =  V(:,:,:,3)  + mu*(A(:,:,:,3,1)+A(:,:,:,3,2)+f4o3*A(:,:,:,3,3))
    dQ(:,:,:,1:3) = dQ(:,:,:,1:3) + V
    dQ(:,:,:,5) = dQ(:,:,:,5) + U(:,:,:,1)*V(:,:,:,1) &
                       +  U(:,:,:,2)*V(:,:,:,2) + U(:,:,:,3)*V(:,:,:,3)

    call damp(Q,dQ)

  END SUBROUTINE rhs

  SUBROUTINE damp(Q,dQ)
    real, dimension(Nx,Ny,Nz,5)  :: Q,dQ

    integer        :: n,k,i,j

    do n = 1,5
       do k = 1,Nz
          do j = 1,Ny
             do i = 1,Nx
                dQ(i,j,k,n) = dQ(i,j,k,n) - (sigX(i)+sigY(j))*(Q(i,j,k,n)-Qtar(i,j,n))
             end do
          end do
       end do
    end do

  END SUBROUTINE damp

  SUBROUTINE getA(U,A,o)
    real, dimension(Nx,Ny,Nz,3)    :: U
    real, dimension(Nx,Ny,Nz,3,3)  :: A
    integer                        :: o
    integer                        :: i,j

!        d^o ui
!  Aij = ------
!        d xj^o
!
    do i = 1,3
       do j = 1,3
          call derivative(U(1,1,1,i),A(1,1,1,i,j),j,o)
       end do
    end do

  END SUBROUTINE getA

END MODULE rhside



