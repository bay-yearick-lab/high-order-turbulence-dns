MODULE mms

  IMPLICIT none

  abstract interface
      function func1D(v, t)
         real :: func1D
         real :: v, t
      end function func1D

      function func3D(x, y, z, t)
         real :: func3D
         real :: x, y, z, t
      end function func3D
  end interface

  type pointer3D
      procedure (func1D), pointer, nopass :: ptr
  end type pointer3D

  type (pointer3D), dimension(3) :: fu, frho, frhou
  procedure (func1D), pointer :: S1_mms_ptr => null()
  procedure (func1D), pointer :: S2_mms_ptr => null()
  procedure (func3D), pointer :: u_mms_ptr => null()
  procedure (func3D), pointer :: rhou_mms_ptr => null()
  procedure (func3D), pointer :: rho_mms_ptr => null()
  procedure (func3D), pointer :: dudt_mms_ptr => null()
  procedure (func3D), pointer :: dudx_mms_ptr => null()

CONTAINS

  FUNCTION dudx_mms(x, y, z, t) RESULT (integrand)
    real :: x, y, z, t, integrand

    integrand = 6.0 * (COS(6.0 * x) - SIN(6.0 * x)) * EXP(-t)
  END FUNCTION dudx_mms

  FUNCTION dudt_mms(x, y, z, t) RESULT (integrand)
    real :: x, y, z, t, integrand

    integrand = - (COS(6.0 * x) + SIN(6.0 * x)) * EXP(-t)
  END FUNCTION dudt_mms

  FUNCTION S2_mms(x, t) RESULT (integrand)
    real :: x, t, integrand

    integrand = EXP(- 2.0 * t) * 6.0 * COS(12.0 * x)
  END FUNCTION S2_mms

  FUNCTION S1_mms(x, t) RESULT (integrand)
    real :: x, t, integrand

    integrand = 35.0 * EXP(-t) * (COS(6.0 * x) + SIN(6.0 * x))
  END FUNCTION S1_mms

  FUNCTION fcos4_mms(v, t) RESULT (integrand)
    real :: v, t, integrand

    integrand = COS(4.0 * v)
  END FUNCTION fcos4_mms

  FUNCTION fsin4_mms(v, t) RESULT (integrand)
    real :: v, t, integrand

    integrand = SIN(4.0 * v)
  END FUNCTION fsin4_mms

  FUNCTION fcos2_mms(v, t) RESULT (integrand)
    real :: v, t, integrand

    integrand = COS(2.0 * v)
  END FUNCTION fcos2_mms

  FUNCTION fsin2_mms(v, t) RESULT (integrand)
    real :: v, t, integrand

    integrand = SIN(2.0 * v)
  END FUNCTION fsin2_mms

  FUNCTION fcos2cos4_mms(v, t) RESULT (integrand)
    real :: v, t, integrand

    integrand = fcos2_mms(v, t) * fcos4_mms(v, t)
  END FUNCTION fcos2cos4_mms

  FUNCTION fsin2sin4_mms(v, t) RESULT (integrand)
    real :: v, t, integrand

    integrand = fsin2_mms(v, t) * fsin4_mms(v, t)
  END FUNCTION fsin2sin4_mms

  FUNCTION u_mms(x, y, z, t) RESULT (integrand)
    real :: x, y, z, t, integrand

    integrand = (fcos4_mms(x, t) * fsin4_mms(y, t) * fcos4_mms(z, t)) * EXP(-t)
  END FUNCTION u_mms

  FUNCTION rho_mms(x, y, z, t) RESULT (integrand)
    real :: x, y, z, t, integrand

    integrand = (fcos2_mms(x, t) * fsin2_mms(y, t) * fcos2_mms(z, t)) * EXP(-t)
  END FUNCTION rho_mms

  FUNCTION rhou_mms(x, y, z, t) RESULT (integrand)
    real :: x, y, z, t, integrand

    integrand = rho_mms(x, y, z, t) * u_mms(x, y, z, t)
  END FUNCTION rhou_mms
END MODULE mms
