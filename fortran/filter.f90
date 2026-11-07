MODULE filter

  USE prms
  USE mesh
  
  IMPLICIT none

  real, allocatable, dimension(:) :: ffx,ffy
  real, allocatable, dimension(:,:) :: ff

CONTAINS
  
  SUBROUTINE filterQ(Q)
    real, dimension(Nx,Ny,Nz,5) :: Q
    
    real, dimension(Nx,Ny,Nz)   :: W

    call filterX(Q,W)
    call filterY(Q,W)
    if (Nz.gt.1) call filterZ(Q,W)

  END SUBROUTINE filterQ
           
  SUBROUTINE filterX(Q,W)
    real, dimension(Nx,Ny,Nz,5) :: Q
    real, dimension(Nx,Ny,Nz)   :: W

    integer           :: n

    do n = 1,5
       call filterX1(Q(:,:,:,n),W)
    end do

  END SUBROUTINE filterX

  SUBROUTINE filterY(Q,W)
    real, dimension(Nx,Ny,Nz,5) :: Q
    real, dimension(Nx,Ny,Nz)   :: W

    integer           :: n

    do n = 1,5
       call filterY1(Q(:,:,:,n),W)
    end do

  END SUBROUTINE filterY

  SUBROUTINE filterZ(Q,W)
    real, dimension(Nx,Ny,Nz,5) :: Q
    real, dimension(Nx,Ny,Nz)   :: W

    integer           :: n

    do n = 1,5
       call filterZ1(Q(:,:,:,n),W)
    end do

  END SUBROUTINE filterZ

  SUBROUTINE filterX1(Q1,W)
    real, dimension(Nx,Ny,Nz)   :: Q1,W

    integer :: i,j,k

    do k = 1,Nz
       do j = 1,Ny
          do i = 1,Nx
             W(i,j,k) = (1.-ff(i,j))*Q1(i,j,k) & 
                  + ff(i,j)*0.25*(Q1(i+1,j,k)+2.*Q1(i,j,k)+Q1(i-1,j,k))
          end do
          W(1,j,k) = Q1(1,j,k)
          W(Nx,j,k) = Q1(Nx,j,k)
       end do
    end do
    Q1 = W

  END SUBROUTINE filterX1
            
  SUBROUTINE filterY1(Q1,W)
    real, dimension(Nx,Ny,Nz)   :: Q1,W

    integer :: i,j,k

    do k = 1,Nz
       do j = 1,Ny
          do i = 1,Nx
             W(i,j,k) = (1.-ff(i,j))*Q1(i,j,k) & 
                  + ff(i,j)*0.25*(Q1(i,j+1,k)+2.*Q1(i,j,k)+Q1(i,j-1,k))
          end do

       end do
       W(:,1,k) = Q1(:,1,k)
       W(:,Ny,k) = Q1(:,Ny,k)
    end do
    Q1 = W

  END SUBROUTINE filterY1

  SUBROUTINE filterZ1(Q1,W)
    real, dimension(Nx,Ny,Nz)   :: Q1,W

    integer :: i,j,k

    do k = 2,Nz-1
       do j = 1,Ny
          do i = 1,Nx
             W(i,j,k) = (1.-ff(i,j))*Q1(i,j,k) & 
                  + ff(i,j)*0.25*(Q1(i,j,k+1)+2.*Q1(i,j,k)+Q1(i,j,k-1))
          end do
       end do
    end do
    do j = 1,Ny
       do i = 1,Nx
          W(i,j,1) = (1.-ff(i,j))*Q1(i,j,1) & 
               + ff(i,j)*0.25*(Q1(i,j,2)+2.*Q1(i,j,1)+Q1(i,j,Nz))
          W(i,j,Nz) = (1.-ff(i,j))*Q1(i,j,Nz) & 
               + ff(i,j)*0.25*(Q1(i,j,1)+2.*Q1(i,j,Nz)+Q1(i,j,Nz-1))
       end do
    end do
    Q1 = W

  END SUBROUTINE filterZ1

  SUBROUTINE initFilters
    
    integer :: i,j

    allocate(ffy(Ny),ffx(Nx),ff(Nx,Ny))

    ffx = 0.
    open(1,file='mix.filt.x')
    do i = 1,Nx        
       if (X(i).gt.Xph) ffx(i) = (X(i)-Xph)**2/(Lx-Xph)**2
       write(1,"(2F20.10)")X(i),ffx(i)
    end do
    close(1)

    ffy = 0.
    open(1,file='mix.filt.y')
    do j = 1,Ny
       if (ABS(Y(j)).gt.Yph) ffy(j) = (ABS(Y(j))-Yph)**2/(Ly/2.-Yph)**2
       write(1,"(2F20.10)")Y(j),ffy(j)
    end do
    close(1)

    do j = 1,Ny
       do i = 1,Nx
          ff(i,j) = MAX(ffx(i),ffy(j))
       end do
    end do
   
  END SUBROUTINE initFilters

END MODULE filter
