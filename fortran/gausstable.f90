MODULE gausstable

  IMPLICIT none

  real, allocatable, dimension(:) :: xi
  real, allocatable, dimension(:) :: wts

  abstract interface
      function func(x, y, z, degree)
         real :: func
         real :: x
         real, optional :: y, z
         integer, optional :: degree
      end function func
  end interface

  procedure (func), pointer :: f_ptr => null()

CONTAINS


  SUBROUTINE gauss_wts_pts(nodes)
    integer :: nodes

    allocate(xi(nodes))
    allocate(wts(nodes))

    if (nodes .eq. 3) then
       wts(1) = 0.8888888888888888; xi(1) = 0.0000000000000000
       wts(2) = 0.5555555555555556; xi(2) = -0.7745966692414834
       wts(3) = 0.5555555555555556; xi(3) = 0.7745966692414834
    else if (nodes .eq. 5) then
       wts(1) = 0.5688888888888889; xi(1) = 0.0000000000000000
       wts(2) = 0.4786286704993665; xi(2) = -0.5384693101056831
       wts(3) = 0.4786286704993665; xi(3) = 0.5384693101056831
       wts(4) = 0.2369268850561891; xi(4) = -0.9061798459386640
       wts(5) = 0.2369268850561891; xi(5) = 0.9061798459386640
    else if (nodes .eq. 7) then
       wts(1) = 0.4179591836734694; xi(1) = 0.0000000000000000
       wts(2) = 0.3818300505051189; xi(2) = 0.4058451513773972
       wts(3) = 0.3818300505051189; xi(3) = -0.4058451513773972
       wts(4) = 0.2797053914892766; xi(4) = -0.7415311855993945
       wts(5) = 0.2797053914892766; xi(5) = 0.7415311855993945
       wts(6) = 0.1294849661688697; xi(6) = -0.9491079123427585
       wts(7) = 0.1294849661688697; xi(7) = 0.9491079123427585
    else if (nodes .eq. 9) then
       wts(1) = 0.3302393550012598; xi(1) = 0.0000000000000000
       wts(2) = 0.1806481606948574; xi(2) = -0.8360311073266358
       wts(3) = 0.1806481606948574; xi(3) = 0.8360311073266358
       wts(4) = 0.0812743883615744; xi(4) = -0.9681602395076261
       wts(5) = 0.0812743883615744; xi(5) = 0.9681602395076261
       wts(6) = 0.3123470770400029; xi(6) = -0.3242534234038089
       wts(7) = 0.3123470770400029; xi(7) = 0.3242534234038089
       wts(8) = 0.2606106964029354; xi(8) = -0.6133714327005904
       wts(9) = 0.2606106964029354; xi(9) = 0.6133714327005904
    else if (nodes .eq. 11) then
       wts(1) = 0.2729250867779006; xi(1) = 0.0000000000000000
       wts(2) = 0.2628045445102467; xi(2) = -0.2695431559523450
       wts(3) = 0.2628045445102467; xi(3) = 0.2695431559523450
       wts(4) = 0.2331937645919905; xi(4) = -0.5190961292068118
       wts(5) = 0.2331937645919905; xi(5) = 0.5190961292068118
       wts(6) = 0.1862902109277343; xi(6) = -0.7301520055740494
       wts(7) = 0.1862902109277343; xi(7) = 0.7301520055740494
       wts(8) = 0.1255803694649046; xi(8) = -0.8870625997680953
       wts(9) = 0.1255803694649046; xi(9) = 0.8870625997680953
       wts(10) = 0.0556685671161737; xi(10) = -0.9782286581460570
       wts(11) = 0.0556685671161737; xi(11) = 0.9782286581460570
    else if (nodes .eq. 13) then
       wts(1) = 0.2325515532308739; xi(1) = 0.0000000000000000
       wts(2) = 0.2262831802628972; xi(2) = -0.2304583159551348
       wts(3) = 0.2262831802628972; xi(3) = 0.2304583159551348
       wts(4) = 0.2078160475368885; xi(4) = -0.4484927510364469
       wts(5) = 0.2078160475368885; xi(5) = 0.4484927510364469
       wts(6) = 0.1781459807619457; xi(6) = -0.6423493394403402
       wts(7) = 0.1781459807619457; xi(7) = 0.6423493394403402
       wts(8) = 0.1388735102197872; xi(8) = -0.8015780907333099
       wts(9) = 0.1388735102197872; xi(9) = 0.8015780907333099
       wts(10) = 0.0921214998377285; xi(10) =-0.9175983992229779
       wts(11) = 0.0921214998377285; xi(11) =0.9175983992229779
       wts(12) = 0.0404840047653159; xi(12) =-0.9841830547185881
       wts(13) = 0.0404840047653159; xi(13) =0.9841830547185881
    else if (nodes .eq. 15) then
       wts(1) = 0.2025782419255613; xi(1) = 0.0000000000000000
       wts(2) = 0.1984314853271116; xi(2) = -0.2011940939974345
       wts(3) = 0.1984314853271116; xi(3) = 0.2011940939974345
       wts(4) = 0.1861610000155622; xi(4) = -0.3941513470775634
       wts(5) = 0.1861610000155622; xi(5) = 0.3941513470775634
       wts(6) = 0.1662692058169939; xi(6) = -0.5709721726085388
       wts(7) = 0.1662692058169939; xi(7) = 0.5709721726085388
       wts(8) = 0.1395706779261543; xi(8) = -0.7244177313601701
       wts(9) = 0.1395706779261543; xi(9) = 0.7244177313601701
       wts(10) = 0.1071592204671719; xi(10) = -0.8482065834104272
       wts(11) = 0.1071592204671719; xi(11) = 0.8482065834104272
       wts(12) = 0.0703660474881081; xi(12) = -0.9372733924007060
       wts(13) = 0.0703660474881081; xi(13) = 0.9372733924007060
       wts(14) = 0.0307532419961173; xi(14) = -0.9879925180204854
       wts(15) = 0.0307532419961173; xi(15) = 0.9879925180204854
    else if (nodes .eq. 21) then
       wts(1) = 0.1460811336496904; xi(1) = 0.0000000000000000
       wts(2) = 0.1445244039899700; xi(2) = -0.1455618541608951
       wts(3) = 0.1445244039899700; xi(3) = 0.1455618541608951
       wts(4) = 0.1398873947910731; xi(4) = -0.2880213168024011
       wts(5) = 0.1398873947910731; xi(5) = 0.2880213168024011
       wts(6) = 0.1322689386333375; xi(6) = -0.4243421202074388
       wts(7) = 0.1322689386333375; xi(7) = 0.4243421202074388
       wts(8) = 0.1218314160537285; xi(8) = -0.5516188358872198
       wts(9) = 0.1218314160537285; xi(9) = 0.5516188358872198
       wts(10) = 0.1087972991671484; xi(10) = -0.6671388041974123
       wts(11) = 0.1087972991671484; xi(11) = 0.6671388041974123
       wts(12) = 0.0934444234560339; xi(12) = -0.7684399634756779
       wts(13) = 0.0934444234560339; xi(13) = 0.7684399634756779
       wts(14) = 0.0761001136283793; xi(14) = -0.8533633645833173
       wts(15) = 0.0761001136283793; xi(15) = 0.8533633645833173
       wts(16) = 0.0571344254268572; xi(16) = -0.9200993341504008
       wts(17) = 0.0571344254268572; xi(17) = 0.9200993341504008
       wts(18) = 0.0369537897708525; xi(18) = -0.9672268385663063
       wts(19) = 0.0369537897708525; xi(19) = 0.9672268385663063
       wts(20) = 0.0160172282577743; xi(20) = -0.9937521706203895
       wts(21) = 0.0160172282577743; xi(21) = 0.9937521706203895
    else
       print*, 'Error in number of Gauss nodes.'
       stop
    end if
  END SUBROUTINE gauss_wts_pts

  FUNCTION f_test(x, y, z, d) RESULT (integrand)
    integer, optional :: d
    real :: x, integrand
    real, optional :: y, z

    if (present(z)) then
       integrand = (x**d + 10.0 * x**2) * &
                   (y**d + 10.0 * y**2) * &
                   (z**d + 10.0 * z**2)
    else if (present(y)) then
       integrand = (x**d + 10.0 * x**2) * &
                   (y**d + 10.0 * y**2)
    else
       integrand = (x**d + 10.0 * x**2)
    end if
  END FUNCTION f_test

  FUNCTION sol_test(x, y, z, d) RESULT (solution)
    integer, optional :: d
    real :: x, solution
    real, optional :: y, z

    if (present(z)) then
       solution = (1.0/(27.0 * (d + 1.0)**3)) * x * y * z * &
                  (3.0 * x**d + 10.0 * (d + 1.0) * x**2) *  &
                  (3.0 * y**d + 10.0 * (d + 1.0) * y**2) *  &
                  (3.0 * z**d + 10.0 * (d + 1.0) * z**2)
    else if (present(y)) then
       solution = (1.0/(9.0 * (d + 1.0)**2)) * x * y *     &
                  (3.0 * x**d + 10.0 * (d + 1.0) * x**2) * &
                  (3.0 * y**d + 10.0 * (d + 1.0) * y**2)
    else
       solution = (1.0/(d + 1.0)) * x**(d + 1) + (10.0/3.0) * x**3
    end if
  END FUNCTION sol_test

  FUNCTION gauss_3quad(f_ptr, n, a, b, c, d, e, f, degree) RESULT (num)

    integer :: n, i, j, k
    integer, optional :: degree
    real :: a, b, c, d, e, f, num, integrand
    procedure (func), pointer :: f_ptr

    num = 0.0
    do k = 1, n
       do j = 1, n
          do i = 1, n
             integrand = f_ptr(x=0.5 * (b - a) * xi(i) + 0.5 * (b + a), &
                               y=0.5 * (d - c) * xi(j) + 0.5 * (d + c), &
                               z=0.5 * (f - e) * xi(k) + 0.5 * (f + e), degree=degree)
             num = num + 0.125 * (f - e) * (d - c) * (b - a) * wts(i) * wts(j) * wts(k) * integrand
          end do
       end do
    end do
  END FUNCTION gauss_3quad

  FUNCTION gauss_2quad(f_ptr, n, a, b, c, d, degree) RESULT (num)
    integer :: n, i, j
    integer, optional :: degree
    real :: a, b, c, d, num, integrand
    procedure (func), pointer :: f_ptr

    num = 0.0
    do j = 1, n
       do i = 1, n
          integrand = f_ptr(x=0.5 * (b - a) * xi(i) + 0.5 * (b + a), &
                            y=0.5 * (d - c) * xi(j) + 0.5 * (d + c), degree=degree)
          num = num + 0.25 * (d - c) * (b - a) * wts(i) * wts(j) * integrand
       end do
    end do
  END FUNCTION gauss_2quad

  FUNCTION gauss_quad(f_ptr, n, a, b, degree) RESULT (num)
    integer :: n, i
    integer, optional :: degree
    real :: a, b, num
    procedure (func), pointer :: f_ptr

    num = 0.0
    do i = 1, n
       num = num + 0.5 * (b - a) * wts(i) * &
           f_ptr(x=0.5 * (b - a) * xi(i) + 0.5 * (b + a), degree=degree)
    end do
  END FUNCTION gauss_quad

  SUBROUTINE test_quad
    integer :: n, degree
    real :: num, sol
    real :: x1, x2, y1, y2, z1, z2

    n = 15
    degree = 10
    x1 = 1.0; x2 = 2.0; y1 = 0.0; y2 = 2.0; z1 = 0.5; z2 = 1.0

    write(*, '(A, I2)') 'Gauss points: ',n, 'Degree of integrand: ',  degree
    call gauss_wts_pts(n)

    ! Single intgral
    f_ptr => f_test
    num = gauss_quad(f_ptr, n, x1, x2, degree)

    f_ptr => sol_test
    sol = f_ptr(x2, degree=degree) - f_ptr(x1, degree=degree)
    print*, 'Single integral'
    print*, 'Integral value = ', num
    print*, 'Relative error = ', ABS(num - sol)/sol

    ! Double intgral
    f_ptr => f_test
    num = gauss_2quad(f_ptr, n, x1, x2, y1, y2, degree)

    f_ptr => sol_test
    sol = (f_ptr(x2, y=y2, degree=degree) - f_ptr(x1, y=y2, degree=degree)) - &
          (f_ptr(x2, y=y1, degree=degree) - f_ptr(x1, y=y1, degree=degree))
    print*, 'Double integral'
    print*, 'Integral value = ', num
    print*, 'Relative error = ', ABS(num - sol)/sol

    ! Triple integral
    f_ptr => f_test
    num = gauss_3quad(f_ptr, n, x1, x2, y1, y2, z1, z2, degree)

    f_ptr => sol_test
    sol = (f_ptr(x2, y2, z2, degree) - f_ptr(x1, y2, z2, degree)) - &
          (f_ptr(x2, y1, z2, degree) - f_ptr(x1, y1, z2, degree)) - &
          (f_ptr(x2, y2, z1, degree) - f_ptr(x1, y2, z1, degree)) - &
          (f_ptr(x2, y1, z1, degree) - f_ptr(x1, y1, z1, degree))
    print*, 'Triple integral'
    print*, 'Integral value = ', num
    print*, 'Relative error = ', ABS(num - sol)/sol

    nullify(f_ptr)
    deallocate(xi, wts)
    stop
  END SUBROUTINE test_quad

END MODULE gausstable
