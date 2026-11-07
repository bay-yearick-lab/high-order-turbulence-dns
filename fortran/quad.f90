MODULE quad

  USE prms
  USE mms
  USE gausstable
  USE bspline

  IMPLICIT none

CONTAINS


  FUNCTION gauss_xbquad(nodes, a, b, Bi, Bj, derivi, derivj) RESULT (num)
    integer :: nodes, i
    integer :: Bi, Bj, derivi, derivj
    real :: a, b, num, gauss_pt

    num = 0.0
    do i = 1, nodes
       gauss_pt = 0.5 * (b - a) * xi(i) + 0.5 * (b + a)
       num = num + 0.5 * (b - a) * wts(i) * &
             xdBIJK(derivi, Bi, d, gauss_pt) * xdBIJK(derivj, Bj, d, gauss_pt)
    end do
  END FUNCTION gauss_xbquad

  FUNCTION gauss_ybquad(nodes, a, b, Bi, Bj, derivi, derivj) RESULT (num)
    integer :: nodes, i
    integer :: Bi, Bj, derivi, derivj
    real :: a, b, num, gauss_pt

    num = 0.0
    do i = 1, nodes
       gauss_pt = 0.5 * (b - a) * xi(i) + 0.5 * (b + a)
       num = num + 0.5 * (b - a) * wts(i) * &
             ydBIJK(derivi, Bi, d, gauss_pt) * ydBIJK(derivj, Bj, d, gauss_pt)
    end do
  END FUNCTION gauss_ybquad

  FUNCTION gauss_zbquad(nodes, a, b, Bi, Bj, derivi, derivj) RESULT (num)
    integer :: nodes, i
    integer :: Bi, Bj, derivi, derivj
    real :: a, b, num, gauss_pt

    num = 0.0
    do i = 1, nodes
       gauss_pt = 0.5 * (b - a) * xi(i) + 0.5 * (b + a)
       num = num + 0.5 * (b - a) * wts(i) * &
             zdBIJK(derivi, Bi, d, gauss_pt) * zdBIJK(derivj, Bj, d, gauss_pt)
    end do
  END FUNCTION gauss_zbquad

  FUNCTION gauss_xbbquad(nodes, a, b, Bi, Bj, Bl, derivi, derivj, derivl) RESULT (num)
    integer :: nodes, i
    integer :: Bi, Bj, Bl, derivi, derivj, derivl
    real :: a, b, num, gauss_pt

    num = 0.0
    do i = 1, nodes
       gauss_pt = 0.5 * (b - a) * xi(i) + 0.5 * (b + a)
       num = num + 0.5 * (b - a) * wts(i) * &
             xdBIJK(derivi, Bi, d, gauss_pt) * &
             xdBIJK(derivj, Bj, d, gauss_pt) * &
             xdBIJK(derivl, Bl, d, gauss_pt)
    end do
  END FUNCTION gauss_xbbquad

  FUNCTION gauss_ybbquad(nodes, a, b, Bi, Bj, Bl, derivi, derivj, derivl) RESULT (num)
    integer :: nodes, i
    integer :: Bi, Bj, Bl, derivi, derivj, derivl
    real :: a, b, num, gauss_pt

    num = 0.0
    do i = 1, nodes
       gauss_pt = 0.5 * (b - a) * xi(i) + 0.5 * (b + a)
       num = num + 0.5 * (b - a) * wts(i) * &
             ydBIJK(derivi, Bi, d, gauss_pt) * &
             ydBIJK(derivj, Bj, d, gauss_pt) * &
             ydBIJK(derivl, Bl, d, gauss_pt)
    end do
  END FUNCTION gauss_ybbquad

  FUNCTION gauss_zbbquad(nodes, a, b, Bi, Bj, Bl, derivi, derivj, derivl) RESULT (num)
    integer :: nodes, i
    integer :: Bi, Bj, Bl, derivi, derivj, derivl
    real :: a, b, num, gauss_pt

    num = 0.0
    do i = 1, nodes
       gauss_pt = 0.5 * (b - a) * xi(i) + 0.5 * (b + a)
       num = num + 0.5 * (b - a) * wts(i) * &
             zdBIJK(derivi, Bi, d, gauss_pt) * &
             zdBIJK(derivj, Bj, d, gauss_pt) * &
             zdBIJK(derivl, Bl, d, gauss_pt)
    end do
  END FUNCTION gauss_zbbquad

  FUNCTION gauss_xrhs(nodes, a, b, Bi, ptr) RESULT (num)
    integer :: nodes, i
    integer :: Bi
    real :: a, b, num, gauss_pt
    procedure (func1D), pointer :: ptr

    num = 0.0
    do i = 1, nodes
       gauss_pt = 0.5 * (b - a) * xi(i) + 0.5 * (b + a)
       num = num + 0.5 * (b - a) * wts(i) * &
             xBJK(Bi, d, gauss_pt) * ptr(gauss_pt, 0.0)
    end do
  END FUNCTION gauss_xrhs

  FUNCTION gauss_yrhs(nodes, a, b, Bj, ptr) RESULT (num)
    integer :: nodes, i
    integer :: Bj
    real :: a, b, num, gauss_pt
    procedure (func1D), pointer :: ptr

    num = 0.0
    do i = 1, nodes
       gauss_pt = 0.5 * (b - a) * xi(i) + 0.5 * (b + a)
       num = num + 0.5 * (b - a) * wts(i) * &
             yBJK(Bj, d, gauss_pt) * ptr(gauss_pt, 0.0)
    end do
  END FUNCTION gauss_yrhs

  FUNCTION gauss_zrhs(nodes, a, b, Bk, ptr) RESULT (num)
    integer :: nodes, i
    integer :: Bk
    real :: a, b, num, gauss_pt
    procedure (func1D), pointer :: ptr

    num = 0.0
    do i = 1, nodes
       gauss_pt = 0.5 * (b - a) * xi(i) + 0.5 * (b + a)
       num = num + 0.5 * (b - a) * wts(i) * &
             zBJK(Bk, d, gauss_pt) * ptr(gauss_pt, 0.0)
    end do
  END FUNCTION gauss_zrhs

  FUNCTION gauss_3Drhs(nodes, x1, x2, y1, y2, z1, z2, Bx, By, Bz, ptr) RESULT (num)
    integer :: nodes
    integer :: Bx, By, Bz
    real :: x1, x2, y1, y2, z1, z2
    procedure (func3D), pointer :: ptr

    integer :: i, j, k
    real :: num, integrand
    real :: x_pt, y_pt, z_pt

    num = 0.0
    do k = 1, nodes
       do j = 1, nodes
          do i = 1, nodes
             x_pt = 0.5 * (x2 - x1) * xi(i) + 0.5 * (x2 + x1)
             y_pt = 0.5 * (y2 - y1) * xi(j) + 0.5 * (y2 + y1)
             z_pt = 0.5 * (z2 - z1) * xi(k) + 0.5 * (z2 + z1)

             integrand = ptr(x_pt, y_pt, z_pt, 0.0)

             num = num + 0.125 * (z2 - z1) * (y2 - y1) * (x2 - x1) * &
                   wts(i) * wts(j) * wts(k) * integrand * &
                   xBJK(Bx, d, x_pt) * yBJK(By, d, y_pt) * zBJK(Bz, d, z_pt)
          end do
       end do
    end do
  END FUNCTION gauss_3Drhs

END MODULE quad
