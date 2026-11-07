MODULE bspline

  USE prms

  IMPLICIT none

  real, allocatable, dimension(:)      :: Xk, Yk, Zk
  real, allocatable, dimension(:)      :: Xc, Yc, Zc
  real, allocatable, dimension(:,:)    :: B
  real, allocatable, dimension(:,:,:)  :: dB
  integer                              :: Ku,Kl

  real, allocatable, dimension(:,:,:)  :: DM
  integer, allocatable, dimension(:,:) :: ipiv

  integer, parameter                   :: Mo=2

CONTAINS

  SUBROUTINE initbspline

    call initxknots
    call inityknots
    call initzknots
    call initxcolpts
    call initycolpts
    call initzcolpts

  END SUBROUTINE initbspline

  SUBROUTINE printsplines

    integer  :: j,i,io
    real     :: yy
    integer, parameter :: iden=100

    do io = d,d
       do i = 1,yNc+io
          do j = 1,iden*Ny
             yy = REAL(j-1)/REAL(iden*Ny-1)*Ly
             write((io+2)*100+i+d,*)yy,yBJK(i,io,yy)
          end do
          close((io+2)*100+i+d)
       end do
    end do

  END SUBROUTINE printsplines

  SUBROUTINE f2a(f,a)
    real, dimension(Ny) :: f,a

    integer   :: info

    a = f
    call DGBTRS('N',Ny,Kl,Ku,1,DM(:,:,0),2*Kl+Ku+1,ipiv(:,0),a,Ny,info)

  END SUBROUTINE f2a

  SUBROUTINE a2f(a,f)
    real, dimension(Ny) :: f,a

    call DGBMV('N',Ny,Ny,Kl,Ku,1.,dB(:,:,0),Ku+Kl+1,a,1,0.,f,1)

  END SUBROUTINE a2f

  SUBROUTINE initzknots

    real      :: Zkm, dz
    integer   :: j

    allocate(Zk(-d-1:zNc+2*d+1))

    Zk(-d-1:0) = 0.
    Zk(zNc:zNc+2*d+1) = Lz
    do j = 1,zNc-1
       Zk(j) = Lz*REAL(j)/REAL(zNc)
    end do

    dz = Zk(1) - Zk(0)
    do j =-d-1,0
       Zk(j) = j * dz
    end do
    do j =zNc+1,zNc+2*d+1
       Zk(j) = Zk(j-1) + dz
    end do

    Zkm = Zk(zNc)
    Zk(1:zNc) = Zk(1:zNc)*Lz/Zkm

    open(1,file='zknt.dat')
    do j = -d-1,zNc+2*d+1
       write(1,*)j,Zk(j)
    end do
    close(1)

  END SUBROUTINE initzknots

  SUBROUTINE inityknots

    real      :: Ykm
    integer   :: j

    allocate(Yk(-d-1:yNc+d+1))

    Yk(-d-1:0) = 0.
    Yk(yNc:yNc+d+1) = Ly
    do j = 1,yNc-1
       Yk(j) = Ly*REAL(j)/REAL(yNc)
    end do

    Ykm = Yk(yNc)
    Yk(1:yNc) = Yk(1:yNc)*Ly/Ykm

    open(1,file='yknt.dat')
    do j = -d-1,yNc+d+1
       write(1,*)j,Yk(j)
    end do
    close(1)

  END SUBROUTINE inityknots

  SUBROUTINE initxknots

    real      :: Xkm, dx
    integer   :: j

    allocate(Xk(-d-1:xNc+2*d+1))

    Xk(-d-1:0) = 0.
    Xk(xNc:xNc+2*d+1) = Lx
    do j = 1,xNc-1
       Xk(j) = Lx*REAL(j)/REAL(xNc)
    end do

    dx = Xk(1) - Xk(0)
    do j =-d-1,0
       Xk(j) = j * dx
    end do
    do j =xNc+1,xNc+2*d+1
       Xk(j) = Xk(j-1) + dx
    end do

    Xkm = Xk(xNc)
    Xk(1:xNc) = Xk(1:xNc)*Lx/Xkm

    open(1,file='xknt.dat')
    do j = -d-1,xNc+2*d+1
       write(1,*)j,Xk(j)
    end do
    close(1)

  END SUBROUTINE initxknots

  SUBROUTINE initxcolpts

    real     :: dx
    integer  :: j

    allocate(Xc(xNc))

    dx = Lx/(xNc)
    do j=1, xNc
       Xc(j) = (j-1) * dx
    end do
  END SUBROUTINE initxcolpts

  SUBROUTINE initycolpts

    real     :: d1,b1
    integer  :: j

    allocate(Yc(Ny))

    if (d.eq.2) then
       do j = 1,Ny
          Yc(j) = Ly*REAL(j-1)/REAL(Ny-1)
       end do
    else
       Yc(1) = 0.
       do j = 2,Ny-1
          call dBzero(j,Yc(j))
       end do
       Yc(Ny) = Ly
    end if
  END SUBROUTINE initycolpts

  SUBROUTINE initzcolpts

    real     :: dz
    integer  :: j

    allocate(Zc(zNc))

    dz = Lz/(zNc)
    do j=1, zNc
       Zc(j) = (j-1) * dz
    end do
  END SUBROUTINE initzcolpts

  SUBROUTINE dBzero(j,co)
    integer :: j
    real    :: co

    integer              :: Nsg
    real, allocatable, dimension(:) :: Bv
    integer :: i
    real    :: c1,c2,f1,f2,fo,yy

    Nsg = 20.*Ny
    allocate (Bv(Nsg))

    do i = 1,Nsg
       yy = Ly*REAL(i-1)/REAL(Nsg-1)
       Bv(i) = yBJK(j,d,yy)
    end do

    i = MAXLOC(Bv,1)

    c1 = Ly*REAL(i-1)/REAL(Nsg-1)
    c2 = Ly*REAL(i)/REAL(Nsg-1)

    do while (ABS(c2-c1).gt.TINY(1.))

       f1 = ydBIJK(1,j,d,c1)
       f2 = ydBIJK(1,j,d,c2)

       co = c1 - f1*(c2-c1)/(f2-f1)

       fo = ydBIJK(1,j,d,co)

       c2 = c1
       c1 = co
    end do

    deallocate(Bv)

  END SUBROUTINE dBzero

  SUBROUTINE generateBs

    integer  :: j,l,i

    real, dimension(Ny,Ny) :: BB

    BB = 0.
    do j = 1,Ny
       do l = 1,Ny
          BB(j,l) = yBJK(j,d,Yc(l))
       end do
    end do
    call getKuKl(Ny,BB,Ku,Kl)
    allocate(B(-Kl:Ku,Ny))
    call makebnd(Ny,BB,Ku,Kl,B)

    do j = 1,Ny
       write(67,"(10E14.6)")(B(l,j),l=-Kl,Ku)
    end do

    allocate(dB(-Kl:Ku,Ny,0:Mo))
    do i = 0,Mo
       do j = 1,Ny
          do l = 1,Ny
             BB(l,j) = ydBIJK(i,j,d,Yc(l))
          end do
       end do
       call getKuKl(Ny,BB,Ku,Kl)
       call makebnd(Ny,BB,Ku,Kl,dB(:,:,i))
    end do

  END SUBROUTINE generateBs

  SUBROUTINE makebnd(N,A,Ku,Kl,B)
    integer                    :: N
    real, dimension(N,N)       :: A
    integer                    :: Ku,Kl
    real, dimension(Ku+Kl+1,N) :: B

    integer                   :: i,j,k

    B = 0.
    do j = 1,N
       k = Ku + 1 - j
       do i = MAX(1,j-Ku), MIN(N,j+Kl)
          B(k+i,j) = A(i,j)
       end do
    end do

  END SUBROUTINE makebnd

  SUBROUTINE getKuKl(N,A,Ku,Kl)
    integer              :: N
    real, dimension(N,N) :: A
    integer              :: Ku,Kl

    integer              :: Kul,Kll
    integer              :: i,j

    Ku = 0; Kl = 0
    do i = 1,N
       Kll = 0
       do j = 1,i-1
          if (A(i,j).ne.0) Kll = Kll + 1
       end do
       Kul = 0
       do j = i+1,N
          if (A(i,j).ne.0) Kul = Kul + 1
       end do
       Ku = MAX(Kul,Ku)
       Kl = MAX(Kll,Kl)
    end do

  END SUBROUTINE getKuKl

  FUNCTION zdBIJK(i,j,k,zz) RESULT (dB)
    integer :: i,j,k,l
    real    :: zz,dB

    dB = 0.
    do l = j-i,j
       dB = dB + zalpha(j,l,i,k) * zBJK(l,k-i,zz)
    end do

  END FUNCTION zdBIJK

  RECURSIVE FUNCTION zalpha(j,l,i,k) RESULT (a)
    integer :: j,l,i,k
    real    :: a

    if (j.eq.l .and. i.eq.0) then
       a = 1.
    else if (i.eq.0) then
       a = 0.
    else if (Zk(l)-Zk(l-k+i-1).ne.0) then
       a = REAL(k-i+1)*(zalpha(j,l+1,i-1,k)-zalpha(j,l,i-1,k))/(Zk(l)-Zk(l-k+i-1))
    else
       a = 0.
    end if

  END FUNCTION zalpha

  FUNCTION ydBIJK(i,j,k,yy) RESULT (dB)
    integer :: i,j,k,l
    real    :: yy,dB

    dB = 0.
    do l = j-i,j
       dB = dB + yalpha(j,l,i,k) * yBJK(l,k-i,yy)
    end do

  END FUNCTION ydBIJK

  RECURSIVE FUNCTION yalpha(j,l,i,k) RESULT (a)
    integer :: j,l,i,k
    real    :: a

    if (j.eq.l .and. i.eq.0) then
       a = 1.
    else if (i.eq.0) then
       a = 0.
    else if (Yk(l)-Yk(l-k+i-1).ne.0) then
       a = REAL(k-i+1)*(yalpha(j,l+1,i-1,k)-yalpha(j,l,i-1,k))/(Yk(l)-Yk(l-k+i-1))
    else
       a = 0.
    end if

  END FUNCTION yalpha

  FUNCTION xdBIJK(i,j,k,xx) RESULT (dB)
    integer :: i,j,k,l
    real    :: xx,dB

    dB = 0.
    do l = j-i,j
       dB = dB + xalpha(j,l,i,k) * xBJK(l,k-i,xx)
    end do

  END FUNCTION xdBIJK

  RECURSIVE FUNCTION xalpha(j,l,i,k) RESULT (a)
    integer :: j,l,i,k
    real    :: a

    if (j.eq.l .and. i.eq.0) then
       a = 1.
    else if (i.eq.0) then
       a = 0.
    else if (Xk(l)-Xk(l-k+i-1).ne.0) then
       a = REAL(k-i+1)*(xalpha(j,l+1,i-1,k)-xalpha(j,l,i-1,k))/(Xk(l)-Xk(l-k+i-1))
    else
       a = 0.
    end if

  END FUNCTION xalpha

  FUNCTION zBJK(j,k,zz) RESULT (Bl)
    integer     :: j,k
    real        :: zz,Bl

    if (j.eq.1.and.zz.eq.0. .or. j.eq.zNc+k.and.zz.eq.Lz) then
       Bl = 1.
    else
       Bl = zBJKr(j,k,zz)
    end if

  END FUNCTION zBJK

  RECURSIVE FUNCTION zBJKr(j,k,zz) RESULT (Bl)
    integer     :: j,k
    real        :: zz,Bl

    real :: b1,b2

    if (k.eq.0) then
       if (zz.gt.Zk(j-1).and.zz.lt.Zk(j)) then
          Bl = 1.
       else if (zz.eq.Zk(j-1) .or. zz.eq.Zk(j)) then
          Bl = 0.5
       else
          Bl = 0.
       end if
    else
       b1 = zBJK(j-1,k-1,zz)
       b2 = zBJK(j,k-1,zz)
       Bl = 0.
       if (b1.ne.0.and.Zk(j-1).ne.Zk(j-k-1))  then
          Bl = (zz-Zk(j-k-1))/(Zk(j-1)-Zk(j-k-1))*b1
       end if
       if (b2.ne.0.and.Zk(j).ne.Zk(j-k))  then
          Bl = Bl + (Zk(j)-zz)/(Zk(j)-Zk(j-k))*b2
       end if
    end if

  END FUNCTION zBJKr

  FUNCTION yBJK(j,k,yy) RESULT (Bl)
    integer     :: j,k
    real        :: yy,Bl

    if (j.eq.1.and.yy.eq.0. .or. j.eq.yNc+k.and.yy.eq.Ly) then
       Bl = 1.
    else
       Bl = yBJKr(j,k,yy)
    end if

  END FUNCTION yBJK

  RECURSIVE FUNCTION yBJKr(j,k,yy) RESULT (Bl)
    integer     :: j,k
    real        :: yy,Bl

    real :: b1,b2

    if (k.eq.0) then
       if (yy.gt.Yk(j-1).and.yy.lt.Yk(j)) then
          Bl = 1.
       else if (yy.eq.Yk(j-1) .or. yy.eq.Yk(j)) then
          Bl = 0.5
       else
          Bl = 0.
       end if
    else
       b1 = yBJK(j-1,k-1,yy)
       b2 = yBJK(j,k-1,yy)
       Bl = 0.
       if (b1.ne.0.and.Yk(j-1).ne.Yk(j-k-1))  then
          Bl = (yy-Yk(j-k-1))/(Yk(j-1)-Yk(j-k-1))*b1
       end if
       if (b2.ne.0.and.Yk(j).ne.Yk(j-k))  then
          Bl = Bl + (Yk(j)-yy)/(Yk(j)-Yk(j-k))*b2
       end if
    end if

  END FUNCTION yBJKr

  FUNCTION xBJK(j,k,xx) RESULT (Bl)
    integer     :: j,k
    real        :: xx,Bl

    if (j.eq.1.and.xx.eq.0. .or. j.eq.xNc+k.and.xx.eq.Lx) then
       Bl = 1.
    else
       Bl = xBJKr(j,k,xx)
    end if

  END FUNCTION xBJK

  RECURSIVE FUNCTION xBJKr(j,k,xx) RESULT (Bl)
    integer     :: j,k
    real        :: xx,Bl

    real :: b1,b2

    if (k.eq.0) then
       if (xx.gt.Xk(j-1).and.xx.lt.Xk(j)) then
          Bl = 1.
       else if (xx.eq.Xk(j-1) .or. xx.eq.Xk(j)) then
          Bl = 0.5
       else
          Bl = 0.
       end if
    else
       b1 = xBJK(j-1,k-1,xx)
       b2 = xBJK(j,k-1,xx)
       Bl = 0.
       if (b1.ne.0.and.Xk(j-1).ne.Xk(j-k-1))  then
          Bl = (xx-Xk(j-k-1))/(Xk(j-1)-Xk(j-k-1))*b1
       end if
       if (b2.ne.0.and.Xk(j).ne.Xk(j-k))  then
          Bl = Bl + (Xk(j)-xx)/(Xk(j)-Xk(j-k))*b2
       end if
    end if

  END FUNCTION xBJKr

END MODULE bspline
