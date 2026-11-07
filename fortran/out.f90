MODULE out

  USE prms
  USE mesh
  USE data
  USE derivatives

  IMPLICIT none

CONTAINS

  SUBROUTINE writemesh

    integer :: i,j,k

    open(1,file='D/mesh.out',form='UNFORMATTED')
    write(1) Nx,Ny
    write(1) ((X(i),i=1,Nx),j=1,Ny),((Y(j),i=1,Nx),j=1,Ny)
    close(1)

    open(1,file='D/mesh3D.out',form='UNFORMATTED')
    write(1) Nx,Ny,Nz
    write(1) (((X(i),i=1,Nx),j=1,Ny),k=1,Nz), &
                       (((Y(j),i=1,Nx),j=1,Ny),k=1,Nz), &
                       (((Z(k),i=1,Nx),j=1,Ny),k=1,Nz)
    close(1)

  END SUBROUTINE writemesh

  SUBROUTINE output(Q,itime)
    real, dimension(Nx,Ny,Nz,5) :: Q
    integer :: itime

    if (out_vort.gt.0.and.MOD(itime,out_vort).eq.0) then
       call writeVort(Q,itime)
    end if
    if (out_q.gt.0.and.MOD(itime,out_q).eq.0) then
       call writeQ(Q,itime)
    end if
    if (out_prof.gt.0.and.MOD(itime,out_prof).eq.0) then
       call writeProf(Q,itime)
    end if
    if (out_cfl.gt.0.and.MOD(itime,out_cfl).eq.0) then
       call writeCFL(Q)
    end if

  END SUBROUTINE output

  SUBROUTINE writeRest(Qloc,itime,time)
    real,dimension(Nx,Ny,Nz,5) :: Qloc
    integer                    :: itime
    real                       :: time
    character(20)              :: otime

    write(otime,"(I10.10,'.t',F8.2)")itime,time
    open(1,file='D/restart.'//otime,form='UNFORMATTED')
    write(1)Qloc,time,itime
    close(1)

  END SUBROUTINE writeRest

  SUBROUTINE writeProf(Qloc,itime)
    real,dimension(Nx,Ny,Nz,5) :: Qloc
    integer                    :: itime
    integer                    :: i
    character(10)              :: otime

    write(otime,"(I10.10)")itime

    open(1,file='D/prof.out'//otime,form='FORMATTED')
    do i = 1,Nx
       write(1,"(2F20.10)") X(i),Qloc(i,1,1,4)
    end do
    close(1)

  END SUBROUTINE writeProf


  SUBROUTINE computeVort(Q,V)
    real, dimension(Nx,Ny,Nz,5) :: Q
    real, dimension(Nx,Ny,Nz,3) :: V

    real, dimension(Nx,Ny,Nz,3) :: U
    real, dimension(Nx,Ny,Nz)   :: dF

    call CtoU(Q,U)

    call derivativeZ(U(1,1,1,2),V(1,1,1,1),1)
    call derivativeY(U(1,1,1,3),dF)
    V(:,:,:,1) = V(:,:,:,1) - dF

    call derivativeX(U(1,1,1,3),V(1,1,1,2),1)
    call derivativeZ(U(1,1,1,1),dF,1)
    V(:,:,:,2) = V(:,:,:,2) - dF

    call derivativeY(U(1,1,1,1),V(1,1,1,3))
    call derivativeX(U(1,1,1,2),dF,1)
    V(:,:,:,3) = V(:,:,:,3) - dF

  END SUBROUTINE computeVort

  SUBROUTINE writeVort(Qloc,itime)

    real,dimension(Nx,Ny,Nz,5) :: Qloc
    integer                    :: itime
    integer                    :: i,j,k
    character(10)              :: otime

    real, dimension(Nx,Ny,Nz,3) :: V
    real, dimension(Nx,Ny,Nz)   :: Vm

    write(otime,"(I10.10)")itime

    call computeVort(Qloc,V)
    Vm = SQRT(V(:,:,:,1)**2 + V(:,:,:,2)**2 + V(:,:,:,3)**2)

    open(1,file='D/vort.out'//otime,form='UNFORMATTED')
    write(1) Nx,Ny,1
    write(1) ((Vm(i,j,1),i=1,Nx),j=1,Ny)
    close(1)

    open(1,file='D/vort3d.out'//otime,form='UNFORMATTED')
    write(1) Nx,Ny,Nz,1
    write(1)(((Vm(i,j,k),i=1,Nx),j=1,Ny),k=1,Nz)
    close(1)

  END SUBROUTINE writeVort

  SUBROUTINE writeQ(Qloc,itime)
    real,dimension(Nx,Ny,Nz,5) :: Qloc
    integer                    :: itime

    integer                    :: i,j,l
    character(10)              :: otime

    write(otime,"(I10.10)")itime

    open(1,file='D/q.out'//otime,form='UNFORMATTED')
    write(1) Nx,Ny,5
    write(1) (((Qloc(i,j,1,l),i=1,Nx),j=1,Ny),l=1,5)
    close(1)
    open(1,file='D/rho.out'//otime,form='UNFORMATTED')
    write(1) Nx,Ny,1
    write(1) ((Qloc(i,j,1,4),i=1,Nx),j=1,Ny)
    close(1)
    open(1,file='D/u.out'//otime,form='UNFORMATTED')
    write(1) Nx,Ny,1
    write(1) ((Qloc(i,j,1,1),i=1,Nx),j=1,Ny)
    close(1)

  END SUBROUTINE writeQ

  SUBROUTINE writeCFL(Q)
    real,dimension(Nx,Ny,Nz,5) :: Q

    integer                    :: i,j,k
    real                       :: c,cfl,cfl_l

    cfl = -1
    do k = 1,Nz
       do j = 1,Ny-1
          do i = 1,Nx-1

             c = SQRT( (gamma-1.)*gamma/Q(i,j,k,4)*(Q(i,j,k,5) &
                  - 0.5*(Q(i,j,k,1)**2+Q(i,j,k,2)**2+           &
                                           Q(i,j,k,3)**2)))/Q(i,j,k,4)
             cfl_l = (ABS(Q(i,j,k,1)/Q(i,j,k,4)) + c)*dt*(REAL(Nx)*sx(i,1)) &
                   + (ABS(Q(i,j,k,2)/Q(i,j,k,4)) + c)*dt*(REAL(Ny)*sy(j,1))! &
   !                + (ABS(Q(i,j,k,3)/Q(i,j,k,4)) + c)*dt/(Z(2)-Z(1))
             cfl = MAX(cfl_l,cfl)
          end do
       end do
    end do

    print *,"CFL = ",cfl

  END SUBROUTINE writeCFL

END MODULE out
