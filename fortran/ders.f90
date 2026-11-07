MODULE derivatives

  USE prms
  USE mesh
  USE bspline

  IMPLICIT none


CONTAINS

  SUBROUTINE initders

    call initbsplineder

  END SUBROUTINE initders

  SUBROUTINE initbsplineder

    allocate(DM(-2*Kl:Ku,Ny,0:Mo))
    allocate(ipiv(Ny,0:Mo))

    DM(-Kl:Ku,:,:) = dB

    call doLUs

  END SUBROUTINE initbsplineder

  SUBROUTINE doLUs

    integer    :: m
    integer    :: info

    do m = 0,Mo
       call DGBTRF(Ny,Ny,Kl,Ku,DM(:,:,m),2*Kl+Ku+1,ipiv(:,m),info)
    end do

  END SUBROUTINE doLUs

  SUBROUTINE derivativeY(F,dF)
    real, dimension(Nx,Ny,Nz) :: F,dF

    real, dimension(Ny) :: a,ff,dff

    integer             :: i,k

    do k = 1,Nz
       do i = 1,Nx
          ff = F(i,:,k)
          call f2a(ff,a)
          call DGBMV('N',Ny,Ny,Kl,Ku,1.,dB(:,:,1),Ku+Kl+1,a,1,0.,dff,1)
          dF(i,:,k) = dff
       end do
    end do

  END SUBROUTINE derivativeY

  SUBROUTINE derivative2Y(F,dF)
    real, dimension(Nx,Ny,Nz) :: F,dF

    call derivativeY(F,dF)
    call derivativeY(df,dF)

  END SUBROUTINE derivative2Y

  SUBROUTINE derivativeYBS1(io,f,df)
    integer             :: io  ! derivative order
    real, dimension(Ny) :: f,df

    real, dimension(Ny) :: a

    if (io.gt.Mo) STOP "der do high order"
    call f2a(f,a)
    call DGBMV('N',Ny,Ny,Kl,Ku,1.,dB(:,:,io),Ku+Kl+1,a,1,0.,df,1)

  END SUBROUTINE derivativeYBS1

  SUBROUTINE derivative(F,dF,d,o)
    real,dimension(Nx,Ny,Nz) :: F,dF
    integer                  :: d,o

    if (o.eq.1) then
       if (d.eq.1) then
          call derivativeX(F,dF,1)
       else if (d.eq.2) then
          call derivativeY(F,dF)
       else if (d.eq.3) then
          call derivativeZ(F,dF,1)
       end if
    else if (o.eq.2) then
       if (d.eq.1) then
          call derivative2X(F,dF)
       else if (d.eq.2) then
          call derivative2Y(F,dF)
       else if (d.eq.3) then
          call derivative2Z(F,dF)
       end if
    end if

  END SUBROUTINE derivative

  SUBROUTINE Lap(F,LF)
    real,dimension(Nx,Ny,Nz) :: F,LF

    real, dimension(Nx,Ny,Nz) :: dF

    call derivative2X(F,LF)
    call derivative2Y(F,dF); LF = LF + dF
    call derivative2Z(F,dF); LF = LF + dF

  END SUBROUTINE Lap

  SUBROUTINE derivativeX(F,dF,in)
    real,dimension(Nx,Ny,Nz) :: F,dF
    integer                  :: in

    integer :: i,j,k

    real :: ds,fac1,fac2
    real :: a1,a2,a3

    ds = REAL(Nx-1)
    fac1 = 8./12.*ds
    fac2 = -1./12.*ds
!!$    a1 = ds*0.79926643E0 ; a2 =-ds*0.18941314E0; a3 = ds*0.0208434E0
!!$    a1 = ds*45.E0/60.E0 ;  a2 = -ds*9.E0/60.E0; a3 = ds*1.E0/60.E0
    a1 = 0.5*ds*1.6492796967641983044E0
    a2 = -0.5*ds*0.41942375741135864349E0
    a3 = 0.5*ds*0.063189272686172994206E0

    do i = 4,Nx-3
       dF(i,:,:) = a3*(F(i+3,:,:) - F(i-3,:,:)) &
            + a2*(F(i+2,:,:) - F(i-2,:,:)) &
            + a1*(F(i+1,:,:) - F(i-1,:,:))
    end do

    if (.not.periodic) then
       i = 3
       dF(i,:,:) = fac2*(F(i+2,:,:) - F(i-2,:,:)) &
            + fac1*(F(i+1,:,:) - F(i-1,:,:))

       i = 2
       dF(i,:,:) = 0.5*(F(i+1,:,:) - F(i-1,:,:))*ds

       i = 1
       dF(i,:,:) = (F(i+1,:,:) - F(i,:,:))*ds

       i = Nx-2
       dF(i,:,:) = fac2*(F(i+2,:,:) - F(i-2,:,:)) &
            &                      + fac1*(F(i+1,:,:) - F(i-1,:,:))

       i = Nx-1
       dF(i,:,:) = 0.5*(F(i+1,:,:) - F(i-1,:,:))*ds

       i = Nx
       dF(i,:,:) = (F(i,:,:) - F(i-1,:,:))*ds

    else if (periodic) then

       dF(3,:,:) = a3*(F(6,:,:) - F(Nx,:,:)) &
            + a2*(F(5,:,:) - F(1,:,:)) &
            + a1*(F(4,:,:) - F(2,:,:))

       dF(2,:,:) = a3*(F(5,:,:) - F(Nx-1,:,:)) &
            + a2*(F(4,:,:) - F(Nx,:,:)) &
            + a1*(F(3,:,:) - F(1,:,:))

       dF(1,:,:) = a3*(F(4,:,:) - F(Nx-2,:,:)) &
            + a2*(F(3,:,:) - F(Nx-1,:,:)) &
            + a1*(F(2,:,:) - F(Nx,:,:))

       dF(Nx,:,:) = a3*(F(3,:,:) - F(Nx-3,:,:)) &
            + a2*(F(2,:,:) - F(Nx-2,:,:)) &
            + a1*(F(1,:,:) - F(Nx-1,:,:))

       dF(Nx-1,:,:) = a3*(F(2,:,:) - F(Nx-4,:,:)) &
            + a2*(F(1,:,:) - F(Nx-3,:,:)) &
            + a1*(F(Nx,:,:) - F(Nx-2,:,:))

       dF(Nx-2,:,:) = a3*(F(1,:,:) - F(Nx-5,:,:)) &
            + a2*(F(Nx,:,:) - F(Nx-4,:,:)) &
            + a1*(F(Nx-1,:,:) - F(Nx-3,:,:))
    end if

    do k = 1,Nz
       do j = 1,Ny
          dF(:,j,k) = sx(:,in)*dF(:,j,k)
       end do
    end do


  END SUBROUTINE derivativeX

  SUBROUTINE derivativeYFD(F,dF,in)
    real, dimension(Nx,Ny,Nz) :: F,dF
    integer                   :: in

    integer :: i,j,k
    real :: ds,fac1,fac2
    real :: a1,a2,a3

    ds = REAL(Ny-1)
    fac1 = 8./12.*ds
    fac2 = -1./12.*ds
!!$    a1 = ds*0.79926643E0 ; a2 =-ds*0.18941314E0; a3 = ds*0.0208434E0
!!$    a1 = ds*45.E0/60.E0 ;  a2 = -ds*9.E0/60.E0; a3 = ds*1.E0/60.E0
    a1 = 0.5*ds*1.6492796967641983044E0
    a2 = -0.5*ds*0.41942375741135864349E0
    a3 = 0.5*ds*0.063189272686172994206E0

    do j = 4,Ny-3
       dF(:,j,:) = a3*(F(:,j+3,:) - F(:,j-3,:)) &
            + a2*(F(:,j+2,:) - F(:,j-2,:)) &
            + a1*(F(:,j+1,:) - F(:,j-1,:))
    end do

    j = 3
    dF(:,j,:) = fac2*(F(:,j+2,:) - F(:,j-2,:)) &
         + fac1*(F(:,j+1,:) - F(:,j-1,:))

    j = 2
    dF(:,j,:) = 0.5*(F(:,j+1,:) - F(:,j-1,:))*ds

    j = 1
    dF(:,j,:) = (F(:,j+1,:) - F(:,j,:))*ds

    j = Ny-2
    dF(:,j,:) = fac2*(F(:,j+2,:) - F(:,j-2,:)) &
     &                      + fac1*(F(:,j+1,:) - F(:,j-1,:))

    j = Ny-1
    dF(:,j,:) = 0.5*(F(:,j+1,:) - F(:,j-1,:))*ds

    j = Ny
    dF(:,j,:) = (F(:,j,:) - F(:,j-1,:))*ds

    do k = 1,Nz
       do j = 1,Ny
          do i = 1,Nx
             dF(i,j,k) = sy(j,in)*dF(i,j,k)
          end do
       end do
    end do

  END SUBROUTINE derivativeYFD


  SUBROUTINE derivativeZ(F,dF,in)
    real, dimension(Nx,Ny,Nz) :: F,dF
    integer                   :: in

    integer :: i,j,k
    real :: ds,fac1,fac2
    real :: a1,a2,a3

    ds = REAL(Nz-1)
    fac1 = 8./12.*ds
    fac2 = -1./12.*ds
!!$    a1 = ds*0.79926643E0 ; a2 =-ds*0.18941314E0; a3 = ds*0.0208434E0
!!$    a1 = ds*45.E0/60.E0 ;  a2 = -ds*9.E0/60.E0; a3 = ds*1.E0/60.E0
    a1 = 0.5*ds*1.6492796967641983044E0
    a2 = -0.5*ds*0.41942375741135864349E0
    a3 = 0.5*ds*0.063189272686172994206E0

    do k = 4,Nz-3
       dF(:,:,k) = a3*(F(:,:,k+3) - F(:,:,k-3)) &
            + a2*(F(:,:,k+2) - F(:,:,k-2)) &
            + a1*(F(:,:,k+1) - F(:,:,k-1))
    end do

    dF(:,:,3) = a3*(F(:,:,6) - F(:,:,Nz)) &
         + a2*(F(:,:,5) - F(:,:,1)) &
         + a1*(F(:,:,4) - F(:,:,2))

    dF(:,:,2) = a3*(F(:,:,5) - F(:,:,Nz-1)) &
         + a2*(F(:,:,4) - F(:,:,Nz)) &
         + a1*(F(:,:,3) - F(:,:,1))

    dF(:,:,1) = a3*(F(:,:,4) - F(:,:,Nz-2)) &
         + a2*(F(:,:,3) - F(:,:,Nz-1)) &
         + a1*(F(:,:,2) - F(:,:,Nz))

    dF(:,:,Nz) = a3*(F(:,:,3) - F(:,:,Nz-3)) &
         + a2*(F(:,:,2) - F(:,:,Nz-2)) &
         + a1*(F(:,:,1) - F(:,:,Nz-1))

    dF(:,:,Nz-1) = a3*(F(:,:,2) - F(:,:,Nz-4)) &
         + a2*(F(:,:,1) - F(:,:,Nz-3)) &
         + a1*(F(:,:,Nz) - F(:,:,Nz-2))

    dF(:,:,Nz-2) = a3*(F(:,:,1) - F(:,:,Nz-5)) &
         + a2*(F(:,:,Nz) - F(:,:,Nz-4)) &
         + a1*(F(:,:,Nz-1) - F(:,:,Nz-3))

    do k = 1,Nz
       do j = 1,Ny
          do i = 1,Nx
             dF(i,j,k) = sz(k,in)*dF(i,j,k)
          end do
       end do
    end do

  END SUBROUTINE derivativeZ


  SUBROUTINE derivative2X(F,dF)
    real, dimension(Nx,Ny,Nz) :: F,dF

    real, dimension(Nx,Ny,Nz) :: F1
    real :: a,b,c,ds
    integer :: i,j,k

    ds = DBLE(Nx-1)**2
    a = ds*6./4.
    b = -ds*3./5./4.
    c = ds*2./20./9.

    do i = 4,Nx-3
       dF(i,:,:) = a*(F(i+1,:,:) - 2.*F(i,:,:) + F(i-1,:,:)) &
            + b*(F(i+2,:,:) - 2.*F(i,:,:) + F(i-2,:,:)) &
            + c*(F(i+3,:,:) - 2.*F(i,:,:) + F(i-3,:,:))
    end do

    do i = 2,3
       dF(i,:,:) = ds*(F(i+1,:,:) - 2.*F(i,:,:) + F(i-1,:,:))
    end do
    do i = Nx-2,Nx-1
       dF(i,:,:) = ds*(F(i+1,:,:) - 2.*F(i,:,:) + F(i-1,:,:))
    end do

    dF(1,:,:) = dF(2,:,:)
    dF(Nx,:,:) = dF(Nx-1,:,:)

    call derivativeX(F,F1,3)
    do k = 1,Nz
       do j = 1,Ny
          dF(:,j,k) = sx(:,2)*dF(:,j,k) + F1(:,j,k)
       end do
    end do

  END SUBROUTINE derivative2X

  SUBROUTINE derivative2YFD(F,dF)
    real, dimension(Nx,Ny,Nz) :: F,dF

    real, dimension(Nx,Ny,Nz) :: F1

    real :: a,b,c,ds
    integer :: i,j,k

    ds = DBLE(Ny-1)**2
    a = ds*6./4.
    b = -ds*3./5./4.
    c = ds*2./20./9.

    do j = 4,Ny-3
       dF(:,j,:) = a*(F(:,j+1,:) - 2.*F(:,j,:) + F(:,j-1,:)) &
            + b*(F(:,j+2,:) - 2.*F(:,j,:) + F(:,j-2,:)) &
            + c*(F(:,j+3,:) - 2.*F(:,j,:) + F(:,j-3,:))
    end do

    do j = 2,3
       dF(:,j,:) = ds*(F(:,j+1,:) - 2.*F(:,j,:) + F(:,j-1,:))
    end do
    do j = Ny-2,Ny-1
       dF(:,j,:) = ds*(F(:,j+1,:) - 2.*F(:,j,:) + F(:,j-1,:))
    end do

    dF(:,1,:) = dF(:,2,:)
    dF(:,Ny,:) = dF(:,Ny-1,:)

    call derivativeYFD(F,F1,3)
    do k = 1,Nz
       do j = 1,Ny
          do i = 1,Nx
             dF(i,j,k) = sy(j,2)*dF(i,j,k) + F1(i,j,k)
          end do
       end do
    end do

  END SUBROUTINE derivative2YFD

  SUBROUTINE derivative2Z(F,dF)
    real, dimension(Nx,Ny,Nz) :: F,dF

    real, dimension(Nx,Ny,Nz) :: F1

    real :: a,b,c,ds
    integer :: i,j,k

    ds = DBLE(Nz-1)**2
    a = ds*6./4.
    b = -ds*3./5./4.
    c = ds*2./20./9.

    do k = 4,Nz-3
       dF(:,:,k) = a*(F(:,:,k+1) - 2.*F(:,:,k) + F(:,:,k-1)) &
            + b*(F(:,:,k+2) - 2.*F(:,:,k) + F(:,:,k-2)) &
            + c*(F(:,:,k+3) - 2.*F(:,:,k) + F(:,:,k-3))
    end do

    dF(:,:,3) = c*(F(:,:,6) - 2.*F(:,:,3) + F(:,:,Nz)) &
         + b*(F(:,:,5) - 2.*F(:,:,3) + F(:,:,1)) &
         + a*(F(:,:,4) - 2.*F(:,:,3) + F(:,:,2))

    dF(:,:,2) = c*(F(:,:,5) - 2.*F(:,:,2) +F(:,:,Nz-1)) &
         + b*(F(:,:,4) - 2.*F(:,:,2) + F(:,:,Nz)) &
         + a*(F(:,:,3) - 2.*F(:,:,2) +F(:,:,1))

    dF(:,:,1) = c*(F(:,:,4) - 2.*F(:,:,1) + F(:,:,Nz-2)) &
         + b*(F(:,:,3) - 2.*F(:,:,1) + F(:,:,Nz-1)) &
         + a*(F(:,:,2) - 2.*F(:,:,1) + F(:,:,Nz))

    dF(:,:,Nz) = c*(F(:,:,3) - 2.*F(:,:,Nz) + F(:,:,Nz-3)) &
         + b*(F(:,:,2) - 2.*F(:,:,Nz) + F(:,:,Nz-2)) &
         + a*(F(:,:,1) - 2.*F(:,:,Nz) + F(:,:,Nz-1))

    dF(:,:,Nz-1) = c*(F(:,:,2) - 2.*F(:,:,Nz-1) + F(:,:,Nz-4)) &
         + b*(F(:,:,1) - 2.*F(:,:,Nz-1) + F(:,:,Nz-3)) &
         + a*(F(:,:,Nz) - 2.*F(:,:,Nz-1) + F(:,:,Nz-2))

    dF(:,:,Nz-2) = c*(F(:,:,1) - 2.*F(:,:,Nz-2) + F(:,:,Nz-5)) &
         + b*(F(:,:,Nz) - 2.*F(:,:,Nz-2) + F(:,:,Nz-4)) &
         + a*(F(:,:,Nz-1) - 2.*F(:,:,Nz-2) + F(:,:,Nz-3))


    call derivativeZ(F,F1,3)
    do k = 1,Nz
       do j = 1,Ny
          do i = 1,Nx
             dF(i,j,k) = sz(k,2)*dF(i,j,k) + F1(i,j,k)
          end do
       end do
    end do

  END SUBROUTINE derivative2Z

END MODULE derivatives
