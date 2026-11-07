MODULE initialcond

  USE prms
  USE data
  USE mesh

  IMPLICIT none

CONTAINS

  SUBROUTINE ic

    real, dimension(Nx,Ny,Nz,5) :: P


    P = 1.  ! needed for N(1)+1 and N(1)+2 so don't get div by zero

    if (which_ic.eq.0) then

       call readRest(Q,itime_start,time_start)

    else

       itime_start = 0
       time_start = 0.

       select case (which_ic)
       case (-1)
          call analytic(P)
       case (1)
          call xsoundwave(P)
       case (2)
          call rpulse(P)
       case (3)
          call xper_mixing_layer(P)
       case (4)
          call nothing(P)
       case default
          stop 'illegal initial condition'
       end select

       call Prim2Cons(P,Q)

    end if

!    if (out_cfl.gt.0) call writeCFL(Q)

  END SUBROUTINE ic


  SUBROUTINE readRest(Q,itime,time)
    real, dimension(Nx,Ny,Nz,5) :: Q
    integer                     :: itime
    real                        :: time

    open(1,file='D/restart.in',form='UNFORMATTED')
    read(1)Q,time,itime
    close(1)

  END SUBROUTINE readRest

  SUBROUTINE nothing(P)
    real, dimension(Nx,Ny,Nz,5) :: P

    P(:,:,:,1:3) = 0.
    P(:,:,:,4) = 1.
    P(:,:,:,5) = 1./gamma

  END SUBROUTINE nothing

  SUBROUTINE analytic(P)
    real, dimension(Nx,Ny,Nz,5) :: P
    integer  :: i,j,k

    do k = 1,Nz
       do j = 1,Ny
          do i = 1,Nx
             P(i,j,k,1) = SIN(Y(j))*SIN(X(i))
             P(i,j,k,2) = EXP(COS(Z(k))*COS(X(i)))
             P(i,j,k,3) = EXP(SIN(X(i))*COS(Y(j))*SIN(Z(k)))
             P(i,j,k,4) = 2. + COS(X(i))*SIN(Z(k))
             P(i,j,k,5) = 2. + EXP(COS(Y(j))+SIN(Z(k)))
          end do
       end do
    end do

  END SUBROUTINE analytic

  SUBROUTINE xsoundwave(P)
    real, dimension(Nx,Ny,Nz,5) :: P
    integer  :: i

    do i = 1,Nx
       P(i,:,:,1) = epsilon1*SIN(X(i))                        !  u
       P(i,:,:,2) = 0.                                       !  v
       P(i,:,:,3) = 0.                                       !  w
       P(i,:,:,4) = 1. + epsilon1*SIN(X(i))                   !  rho
       P(i,:,:,5) = 1./gamma + epsilon1*SIN(X(i))             !  p
    end do

  END SUBROUTINE xsoundwave


  SUBROUTINE rpulse(P)
    real, dimension(Nx,Ny,Nz,5) :: P
    real     :: r,Xm,Ym,Zm
    integer  :: i,j,k

    Xm = X(Nx/2)
    Ym = Y(Ny/2)
    Zm = Z(Nz/2)
    do k = 1,Nz
       do j = 1,Ny
          do i = 1,Nx
             r = SQRT( (X(i)-Xm)*(X(i)-Xm) + (Y(j)-Ym)*(Y(j)-Ym) )
             P(i,j,k,1) = 0.
             P(i,j,k,2) = 0.
             P(i,j,k,3) = 0.
             P(i,j,k,4) = 1.
             P(i,j,k,5) = 1./gamma + epsilon1*EXP(-3.*r*r)
          end do
       end do
    end do

  END SUBROUTINE rpulse

  SUBROUTINE xper_mixing_layer(P)
    real, dimension(Nx,Ny,Nz,5) :: P
    real     :: sig,Ym
    integer  :: i,j,k

    sig = 6.
    Ym = Ly/2
    do k = 1,Nz
       do j = 1,Ny
          do i = 1,Nx
             P(i,j,k,1) =  0.5*Uo*(TANH(sig*Y(j))+1.)
             P(i,j,k,2) = (0.5*epsilon1*COS(2.*Pi*X(i)/X(Nx)) &
                  + epsilon1*SIN(4.*Pi*X(i)/Lx))*EXP(-sig*Y(j)**2/4.) &
                  + epsilon1*SIN(4.*Pi*Z(k)/Lz)
             P(i,j,k,3) = 0.
             P(i,j,k,4) = 1.
             P(i,j,k,5) = 1./gamma
          end do
       end do
    end do

  END SUBROUTINE xper_mixing_layer

END MODULE initialcond
