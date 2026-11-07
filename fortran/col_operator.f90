MODULE col_operator

  USE prms
  USE bspline
  USE quad
  USE col_fft

  IMPLICIT none

  ! matrices
  real, allocatable, dimension(:)       :: pivoty
  real, allocatable, dimension(:)       :: ppivoty
  real, allocatable, dimension(:, :)    :: LUy, PLUy
  real, allocatable, dimension(:, :)    :: Px, Py, Pz
  real, allocatable, dimension(:, :, :) :: By, dBy
  real, allocatable, dimension(:, :, :) :: Bx, Bz
  real, allocatable, dimension(:)       :: lamda, vector_lamda, pxrow, pzrow
  real, allocatable, dimension(:, :)    :: BBx, BBy, BBz, dBTy, dBTBPy
  real, allocatable, dimension(:, :)    :: BTx, BTy, BTz
  real, allocatable, dimension(:, :)    :: BTBPx, BTBPy, BTBPz

  ! number of bands
  integer                               :: BKu, BKl, PKu, PKl, pinfoy, infoy

  ! global NS variables
  real, allocatable, dimension(:, :)       :: global_tempx, global_tempy, global_tempz
  real, allocatable, dimension(:, :, :)    :: global_rho, global_T
  real, allocatable, dimension(:, :, :)    :: global_temp1, global_temp2
  real, allocatable, dimension(:, :, :, :) :: global_u

CONTAINS

  SUBROUTINE allocate_operator
    allocate(pivoty(Ny), ppivoty(Ny))
    allocate(Px(xNc, xNc), Py(Ny, Ny), Pz(zNc, zNc))
    allocate(Bx(xNc, xNc, 0:Mo), By(Ny, Ny, 0:Mo), Bz(zNc, zNc, 0:Mo))
    allocate(lamda(3), pxrow(xNc), pzrow(zNc), vector_lamda(Ny+2))
    allocate(BBx(xNc, xNc), BBy(Ny, Ny), BBz(zNc, zNc))
    allocate(BTx(xNc, xNc), BTy(Ny, Ny), BTz(zNc, zNc))
    allocate(BTBPx(xNc, xNc), BTBPy(Ny, Ny), BTBPz(zNc, zNc))

    allocate(global_tempx(xNc, 2), global_tempy(Ny, 2), global_tempz(zNc, 2))
    allocate(global_rho(xNc, Ny, zNc), global_T(xNc, Ny, zNc))
    allocate(global_temp1(xNc, Ny, zNc), global_temp2(xNc, Ny, zNc))
    allocate(global_u(xNc, Ny, zNc, 3))
  END SUBROUTINE allocate_operator

  SUBROUTINE deallocate_operator
    deallocate(LUy, pivoty, dBy)
    deallocate(Bx, By, Bz, Px, Py, Pz)
    deallocate(lamda, pxrow, pzrow, ppivoty, dBTy, vector_lamda)
    deallocate(BTx, BTy, BTz)
    deallocate(BBx, BBy, BBz, BTBPx, BTBPy, BTBPz)
    call destroy_fft

    deallocate(global_tempx, global_tempy, global_tempz)
    deallocate(global_temp1, global_temp2)
    deallocate(global_rho, global_u, global_T)
  END SUBROUTINE deallocate_operator

  SUBROUTINE kron(matXYZ, matXY, matX, matY, matZ)
    real, dimension(xNc, xNc) :: matX
    real, dimension(Ny, Ny) :: matY
    real, dimension(zNc, zNc) :: matZ
    real, dimension(xNc*Ny, xNc*Ny) :: matXY
    real, dimension(xNc*Ny*zNc, xNc*Ny*zNc) :: matXYZ

    real :: i, j
    integer :: m, row0, row1, col0, col1

    m = Ny

    do j=1, xNc
       do i=1, xNc
          row0 = (i-1) * m + 1
          col0 = (j-1) * m + 1
          row1 = row0 + m - 1
          col1 = col0 + m - 1

          matXY(row0:row1, col0:col1) = matX(i, j) * matY
       end do
    end do

    m = zNc

    do j=1, xNc*Ny
       do i=1, xNc*Ny
          row0 = (i-1) * m + 1
          col0 = (j-1) * m + 1
          row1 = row0 + m - 1
          col1 = col0 + m - 1

          matXYZ(row0:row1, col0:col1) = matXY(i, j) * matZ
       end do
    end do
  END SUBROUTINE kron

  SUBROUTINE generate_operator
    real, allocatable, dimension(:, :) :: xrow, zrow
    integer :: nodes
    integer :: i, j, k, deriv

    allocate(xrow(xNc, 0:Mo), zrow(zNc, 0:Mo))

    call allocate_operator

    ! x-dimension
    do deriv=0, Mo
       do j=1, xNc
            if (deriv .eq. 1) then
               xrow(j, deriv) = -xdBIJK(deriv, d, d, Xc(j))
            else
               xrow(j, deriv) = xdBIJK(deriv, d, d, Xc(j))
            end if
       end do
    end do

    Bx = 0.0
    do deriv=0, Mo
       do i=1, xNc
          Bx(i, :, deriv) = CSHIFT(xrow(:, deriv), (d+1)/2-i)
       end do
    end do
    xrow(:, 0) = CSHIFT(xrow(:, 0), (d+1)/2-1)

    ! y-dimension
    By = 0.0
    do deriv=0, Mo
       do j=1, Ny
           do i=1, Ny
              By(i, j, deriv) = ydBIJK(deriv, j, d, Yc(i))
           end do
       end do
    end do

    ! make banded matrix
    call getKuKl(Ny, By(:, :, 0), BKu, BKl)
    allocate(dBy(-BKl:BKu, Ny, 0:Mo))
    allocate(LUy(-2*BKl:BKu, Ny))

    do deriv=0, Mo
       call makebnd(Ny, By(:, :, deriv), BKu, BKl, dBy(:, :, deriv))
    end do

    LUy(-BKl:BKu, :) = dBy(:, :, 0)

    ! LU factorization
    call DGBTRF(Ny, Ny, BKl, BKu, LUy, 2*BKl+BKu+1, pivoty, infoy)

    ! z-dimension
    do deriv=0, Mo
       do j=1, zNc
            if (deriv .eq. 1) then
               zrow(j, deriv) = -zdBIJK(deriv, d, d, Zc(j))
            else
               zrow(j, deriv) = zdBIJK(deriv, d, d, Zc(j))
            end if
       end do
    end do

    Bz = 0.0
    do deriv=0, Mo
       do i=1, zNc
          Bz(i, :, deriv) = CSHIFT(zrow(:, deriv), (d+1)/2-i)
       end do
    end do
    zrow(:, 0)  = CSHIFT(zrow(:, 0), (d+1)/2-1)

    call init_fft(xrow(:, 0), zrow(:, 0))

    xrow = 0.0; zrow = 0.0

    ! Number of nodes give exact quadrature for degree d**2 B-spline
    nodes = d + 2

    ! Gauss weights and nodes
    call gauss_wts_pts(nodes)

    xrow = 0.0
    do j=1, 2*d+1
       do i=1, 2*d+1
          xrow(j, 2) = xrow(j, 2) + gauss_xbquad(nodes, Xk(i-1), Xk(i), d+1, j, 2, 2)
       end do
    end do

    Px = 0.0
    do i=1, xNc
       Px(i, :) = CSHIFT(xrow(:, 2), d-(i-1))
    end do

    Py = 0.0
    do j = 1, Ny
       do i = 1, Ny
          if (ABS(i-j) .le. d) then
             ! Loop over all knot intervals
             do k = 0, yNc-1
                Py(i, j) = Py(i, j) + gauss_ybquad(nodes, Yk(k), Yk(k+1), i, j, 2, 2)
             end do
          end if
       end do
    end do

    zrow = 0.0
    do j=1, 2*d+1
       do i=1, 2*d+1
          zrow(j, 2) = zrow(j, 2) + gauss_zbquad(nodes, Zk(i-1), Zk(i), d+1, j, 2, 2)
       end do
    end do

    Pz = 0.0
    do i=1, zNc
       Pz(i, :) = CSHIFT(zrow(:, 2), d-(i-1))
    end do

    ! scale double derivative by dx**2
    Px = (Px * (Lx/xNc)**2)
    Py = (Py * (Ly/Ny)**2)
    Pz = (Pz * (Lz/zNc)**2)

    ! make banded matrix
    BTx = TRANSPOSE(Bx(:,:,0))
    BTy = TRANSPOSE(By(:,:,0))
    BTz = TRANSPOSE(Bz(:,:,0))

    ! setting boundaries to have diagonal 1
    Py(1, :) = 0.0;  Py(Ny, :) = 0.0
    BTy(1, :) = 0.0; BTy(Ny, :) = 0.0
    BTy(1, 1) = 1.0; BTy(Ny, Ny) = 1.0

    allocate(dBTy(-BKl:BKu, Ny))
    call makebnd(Ny, BTy, BKu, BKl, dBTy)

    ! construct BB matrix
    BBx = MATMUL(BTx, Bx(:,:,0))
    BBy = MATMUL(BTy, By(:,:,0))
    BBz = MATMUL(BTz, Bz(:,:,0))

    call getKuKl(Ny, BBy, PKu, PKl)
    allocate(dBTBPy(-PKl:PKu, Ny))
    allocate(PLUy(-2*PKl:PKu, Ny))

    ! Deallocate array
    deallocate(xi, wts)
    deallocate(xrow, zrow)
  END SUBROUTINE generate_operator

  SUBROUTINE re_operator
    BTBPx = BBx + lamda(1) * Px
    BTBPy = BBy + lamda(2) * Py
    BTBPz = BBz + lamda(3) * Pz
    pxrow = BTBPx(:,1)
    pzrow = BTBPz(:,1)

    call makebnd(Ny, BTBPy, PKu, PKl, dBTBPy)
    PLUy(-PKl:PKu, :) = dBTBPy
    ! LU factorization
    call DGBTRF(Ny, Ny, PKl, PKu, PLUy, 2*PKl+PKu+1, ppivoty, pinfoy)

    if (allocated(pFFTx1)) then
        call pdestroy_fft
    end if
    call pinit_fft(pxrow, pzrow)
  END SUBROUTINE re_operator

  SUBROUTINE vector_re_operator
    integer :: i, j

    BTBPx = BBx + vector_lamda(1) * Px
    BTBPz = BBz + vector_lamda(2) * Pz

    BTBPy = 0.0
    do j=1, Ny
       do i=1, Ny
          BTBPy(i, j) = BBy(i, j) + vector_lamda(j+2) * Py(i, j)
       end do
    end do

    pxrow = BTBPx(:,1)
    pzrow = BTBPz(:,1)

    call makebnd(Ny, BTBPy, PKu, PKl, dBTBPy)
    PLUy(-PKl:PKu, :) = dBTBPy
    ! LU factorization
    call DGBTRF(Ny, Ny, PKl, PKu, PLUy, 2*PKl+PKu+1, ppivoty, pinfoy)

    if (allocated(pFFTx1)) then
        call pdestroy_fft
    end if
    call pinit_fft(pxrow, pzrow)
  END SUBROUTINE vector_re_operator

END MODULE col_operator
