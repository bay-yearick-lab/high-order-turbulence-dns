MODULE bc

  USE prms
  USE mesh

  IMPLICIT none

  real, allocatable, dimension(:)     :: sigX,sigY,sigZ
  real, allocatable,dimension(:,:,:)  :: Qtar,Utar

CONTAINS

  SUBROUTINE initbc(Q)
    real, dimension(Nx,Ny,Nz,5) :: Q

    integer  :: i,j,k

    allocate (sigX(Nx),sigY(Ny),sigZ(Nz))

    sigX = 0.
    open(1,file='mix.sig.x')
    if (.not.periodic) then
       do i = 1,Nx
          if (X(i).gt.(Lx-Wxr)) then
             sigX(i) = sigmax*(X(i)-(Lx-Wxr))**2/Wxr**2
          else if (X(i).lt.Wxl) then
             sigX(i) = sigmax*(Wxl-X(i))**2/Wxl**2
          end if
          write(1,"(2F20.10)")X(i),sigX(i)
       end do
    end if
    close(1)

    sigY = 0.
    open(1,file='mix.sig.y')
    do j = 1,Ny
       if (ABS(Y(j)).gt.(Ly/2.-Wy)) then
          sigY(j) = sigmax*(ABS(Y(j))-(Ly/2.-Wy))**2/Wy**2
!!$       else if (Y(j).lt.Wy) then
!!$          sigY(j) = sigmax*(Wy-Y(j))**2/Wy**2
       end if
       write(1,"(2F20.10)")Y(j),sigY(j)
    end do
    close(1)

    sigZ = 0.
    open(1,file='mix.sig.z')
    do k = 1,Nz
!!$       if (ABS(Z(k)).gt.(Lz/2.-Wz)) then
!!$          sigZ(k) = sigmax*(ABS(Z(k))-(Lz/2.-Wz))**2/Wz**2
!!$       else if (Z(k).lt.Wz) then
!!$          sigZ(k) = sigmax*(Wz-Z(k))**2/Wz**2
!!$       end if
       write(1,"(2F20.10)")Z(k),sigZ(k)
    end do
    close(1)


    call initQtar(Q)

  END SUBROUTINE initbc

  SUBROUTINE initQtar(Q)
    real, dimension(Nx,Ny,Nz,5) :: Q

    allocate(Qtar(Nx,Ny,5),Utar(Nx,Ny,6))
    Qtar = Q(:,:,1,:)

    Utar(:,:,1) = Qtar(:,:,1)/Qtar(:,:,4)
    Utar(:,:,2) = Qtar(:,:,2)/Qtar(:,:,4)
    Utar(:,:,3) = Qtar(:,:,3)/Qtar(:,:,4)
    Utar(:,:,4) = Qtar(:,:,4)
    Utar(:,:,5) = (gamma-1.)*(Qtar(:,:,5) - 0.5*(Qtar(:,:,1)**2  &
         + Qtar(:,:,2)**2 + Qtar(:,:,3)**2)/Qtar(:,:,4))
    Utar(:,:,6) = SQRT(gamma*Utar(:,:,5)/Utar(:,:,4))

  END SUBROUTINE initQtar

  SUBROUTINE toChar(Q)
    real, dimension(Nx,Ny,Nz,5)  :: Q

    call XtoChar1(1,Q)
    call XtoChar1(Nx,Q)
    call YtoChar1(1,Q)
    call YtoChar1(Ny,Q)

  END SUBROUTINE toChar

  SUBROUTINE YtoChar1(j,Q)
    real, dimension(Nx,Ny,Nz,5)  :: Q
    integer  :: j

    integer  :: k
    real, dimension(Nx,Nz,5)  :: ux

    do k = 1,Nz
       ux(2:Nx-1,k,1) = Q(2:Nx-1,j,k,1) / Q(2:Nx-1,j,k,4)  - Utar(2:Nx-1,j,1)
       ux(2:Nx-1,k,2) = Q(2:Nx-1,j,k,2) / Q(2:Nx-1,j,k,4)  - Utar(2:Nx-1,j,2)
       ux(2:Nx-1,k,3) = Q(2:Nx-1,j,k,3) / Q(2:Nx-1,j,k,4)  - Utar(2:Nx-1,j,3)
       ux(2:Nx-1,k,4) = Q(2:Nx-1,j,k,4)                    - Utar(2:Nx-1,j,4)
       ux(2:Nx-1,k,5) = (gamma-1.)*(Q(2:Nx-1,j,k,5) - 0.5*(Q(2:Nx-1,j,k,1)**2  &
            + Q(2:Nx-1,j,k,2)**2 + Q(2:Nx-1,j,k,3)**2)/Q(2:Nx-1,j,k,4)) - Utar(2:Nx-1,j,5)
    end do

    Q(2:Nx-1,j,:,1) = ux(2:Nx-1,:,5) - Utar(2:Nx-1,:,6)**2*ux(2:Nx-1,:,4)
    Q(2:Nx-1,j,:,2) = Utar(2:Nx-1,:,4)*Utar(2:Nx-1,:,6)*ux(2:Nx-1,:,1)
    Q(2:Nx-1,j,:,3) = Utar(2:Nx-1,:,4)*Utar(2:Nx-1,:,6)*ux(2:Nx-1,:,3)
    Q(2:Nx-1,j,:,4) = ux(2:Nx-1,:,5) + Utar(2:Nx-1,:,4)*Utar(2:Nx-1,:,6)*ux(2:Nx-1,:,2)
    Q(2:Nx-1,j,:,5) = ux(2:Nx-1,:,5) - Utar(2:Nx-1,:,4)*Utar(2:Nx-1,:,6)*ux(2:Nx-1,:,2)

  END SUBROUTINE YtoChar1

  SUBROUTINE XtoChar1(i,Q)
    real, dimension(Nx,Ny,Nz,5)  :: Q
    integer  :: i

    integer  :: k
    real, dimension(Ny,Nz,5)  :: ux

    do k = 1,Nz
       ux(:,k,1) = Q(i,:,k,1) / Q(i,:,k,4)  - Utar(i,:,1)
       ux(:,k,2) = Q(i,:,k,2) / Q(i,:,k,4)  - Utar(i,:,2)
       ux(:,k,3) = Q(i,:,k,3) / Q(i,:,k,4)  - Utar(i,:,3)
       ux(:,k,4) = Q(i,:,k,4)               - Utar(i,:,4)
       ux(:,k,5) = (gamma-1.)*(Q(i,:,k,5) - 0.5*(Q(i,:,k,1)**2  &
            + Q(i,:,k,2)**2 + Q(i,:,k,3)**2)/Q(i,:,k,4)) - Utar(i,:,5)
    end do

    Q(i,:,:,1) = ux(:,:,5) - Utar(:,:,6)**2*ux(:,:,4)
    Q(i,:,:,2) = Utar(:,:,4)*Utar(:,:,6)*ux(:,:,2)
    Q(i,:,:,3) = Utar(:,:,4)*Utar(:,:,6)*ux(:,:,3)
    Q(i,:,:,4) = ux(:,:,5) + Utar(:,:,4)*Utar(:,:,6)*ux(:,:,1)
    Q(i,:,:,5) = ux(:,:,5) - Utar(:,:,4)*Utar(:,:,6)*ux(:,:,1)

  END SUBROUTINE XtoChar1

  SUBROUTINE applybc(Q)
    real, dimension(Nx,Ny,Nz,5)  :: Q

    if (.not.periodic) then
       Q(1,:,:,1) = 0.
       Q(1,:,:,2) = 0.
       Q(1,:,:,3) = 0.
       Q(1,:,:,4) = 0.

       Q(Nx,:,:,5) = 0.
    end if

    Q(:,1,:,1) = 0.
    Q(:,1,:,2) = 0.
    Q(:,1,:,3) = 0.
    Q(:,1,:,4) = 0.

    Q(:,Ny,:,1) = 0.
    Q(:,Ny,:,2) = 0.
    Q(:,Ny,:,3) = 0.
    Q(:,Ny,:,5) = 0.

  END SUBROUTINE applybc

  SUBROUTINE fromChar(Q)
    real, dimension(Nx,Ny,Nz,5)  :: Q

    call XfromChar1(1,Q)
    call XfromChar1(Nx,Q)
    call YfromChar1(1,Q)
    call YfromChar1(Ny,Q)

  END SUBROUTINE fromChar


  SUBROUTINE XfromChar1(i,Q)
    integer                      :: i
    real, dimension(Nx,Ny,Nz,5)  :: Q

    real, dimension(Ny,Nz,5)  :: ux
    integer  :: k

    do k = 1,Nz
       ux(:,k,1) = 0.5*(Q(i,:,k,4)-Q(i,:,k,5))/(Utar(i,:,4)*Utar(i,:,6)) + Utar(i,:,1)
       ux(:,k,2) = Q(i,:,k,2)/(Utar(i,:,4)*Utar(i,:,6)) + Utar(i,:,2)
       ux(:,k,3) = Q(i,:,k,3)/(Utar(i,:,4)*Utar(i,:,6)) + Utar(i,:,3)
       ux(:,k,4) = (0.5*(Q(i,:,k,4)+Q(i,:,k,5))-Q(i,:,k,1))/Utar(i,:,6)**2 + Utar(i,:,4)
       ux(:,k,5) = 0.5*(Q(i,:,k,4)+Q(i,:,k,5)) + Utar(i,:,5)
    end do

    Q(i,:,:,1) = ux(:,:,1)*ux(:,:,4)
    Q(i,:,:,2) = ux(:,:,2)*ux(:,:,4)
    Q(i,:,:,3) = ux(:,:,3)*ux(:,:,4)
    Q(i,:,:,4) = ux(:,:,4)
    Q(i,:,:,5) = ux(:,:,5)/(gamma-1.) &
         + 0.5*ux(:,:,4)*(ux(:,:,1)**2+ux(:,:,2)**2+ux(:,:,3)**2)

  END SUBROUTINE XfromChar1

  SUBROUTINE YfromChar1(j,Q)
    integer                      :: j
    real, dimension(Nx,Ny,Nz,5)  :: Q

    real, dimension(Nx,Nz,5)  :: ux
    integer  :: k

    do k = 1,Nz
       ux(2:Nx-1,k,2) = 0.5*(Q(2:Nx-1,j,k,4)-Q(2:Nx-1,j,k,5))/(Utar(2:Nx-1,j,4)*Utar(2:Nx-1,j,6)) + Utar(2:Nx-1,j,2)
       ux(2:Nx-1,k,1) = Q(2:Nx-1,j,k,2)/(Utar(2:Nx-1,j,4)*Utar(2:Nx-1,j,6)) + Utar(2:Nx-1,j,1)
       ux(2:Nx-1,k,3) = Q(2:Nx-1,j,k,3)/(Utar(2:Nx-1,j,4)*Utar(2:Nx-1,j,6)) + Utar(2:Nx-1,j,3)
       ux(2:Nx-1,k,4) = (0.5*(Q(2:Nx-1,j,k,4)+Q(2:Nx-1,j,k,5))-Q(2:Nx-1,j,k,1))/Utar(2:Nx-1,j,6)**2 + Utar(2:Nx-1,j,4)
       ux(2:Nx-1,k,5) = 0.5*(Q(2:Nx-1,j,k,4)+Q(2:Nx-1,j,k,5)) + Utar(2:Nx-1,j,5)
    end do

    Q(2:Nx-1,j,:,1) = ux(2:Nx-1,:,1)*ux(2:Nx-1,:,4)
    Q(2:Nx-1,j,:,2) = ux(2:Nx-1,:,2)*ux(2:Nx-1,:,4)
    Q(2:Nx-1,j,:,3) = ux(2:Nx-1,:,3)*ux(2:Nx-1,:,4)
    Q(2:Nx-1,j,:,4) = ux(2:Nx-1,:,4)
    Q(2:Nx-1,j,:,5) = ux(2:Nx-1,:,5)/(gamma-1.) &
         + 0.5*ux(2:Nx-1,:,4)*(ux(2:Nx-1,:,1)**2+ux(2:Nx-1,:,2)**2+ux(2:Nx-1,:,3)**2)

  END SUBROUTINE YfromChar1

  SUBROUTINE dobc(Q)
    real, dimension(Nx,Ny,Nz,5)  :: Q

    call toChar(Q)
    call applybc(Q)
    call fromChar(Q)

  END SUBROUTINE dobc

END MODULE bc
