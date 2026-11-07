MODULE col_FFTW3
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
END MODULE col_FFTW3

MODULE col_fft
  USE prms
  USE bspline
  USE col_FFTW3

  IMPLICIT NONE
  integer*8 :: FFT_xplan, iFFT_xplan, FFT_zplan, iFFT_zplan
  complex, allocatable, dimension(:) :: FFTx1, FFTx2, FFTx3
  complex, allocatable, dimension(:) :: FFTz1, FFTz2, FFTz3
  complex, allocatable, dimension(:) :: FFT_x_Matrix, FFT_z_Matrix

  integer*8 :: pFFT_xplan, piFFT_xplan, pFFT_zplan, piFFT_zplan
  complex, allocatable, dimension(:) :: pFFTx1, pFFTx2, pFFTx3
  complex, allocatable, dimension(:) :: pFFTz1, pFFTz2, pFFTz3
  complex, allocatable, dimension(:) :: pFFT_x_Matrix, pFFT_z_Matrix

CONTAINS

  SUBROUTINE init_FFT(xrow, zrow)
    real, dimension(xNc) :: xrow
    real, dimension(zNc) :: zrow
    integer :: plan

    ! z dimension
    allocate(FFT_z_Matrix(zNc))
    allocate(FFTz1(zNc), FFTz2(zNc), FFTz3(zNc))
    call dfftw_plan_dft_1d(plan, zNc, FFTz1, FFTz2, FFTW_FORWARD, FFTW_ESTIMATE)
    FFT_zplan = plan
    call dfftw_plan_dft_1d(plan, zNc, FFTz1, FFTz2, FFTW_BACKWARD, FFTW_ESTIMATE)
    iFFT_zplan = plan

    FFTz1 = zrow
    call dfftw_execute_dft(FFT_zplan, FFTz1, FFT_z_Matrix)

    ! x dimension
    allocate(FFT_x_Matrix(xNc))
    allocate(FFTx1(xNc), FFTx2(xNc), FFTx3(xNc))
    call dfftw_plan_dft_1d(plan, xNc, FFTx1, FFTx2, FFTW_FORWARD, FFTW_ESTIMATE)
    FFT_xplan = plan
    call dfftw_plan_dft_1d(plan, xNc, FFTx1, FFTx2, FFTW_BACKWARD, FFTW_ESTIMATE)
    iFFT_xplan = plan

    FFTx1 = xrow
    call dfftw_execute_dft(FFT_xplan, FFTx1, FFT_x_Matrix)
  END SUBROUTINE init_FFT

  SUBROUTINE destroy_FFT
    deallocate(FFTx1, FFTx2, FFTx3, FFT_x_Matrix)
    deallocate(FFTz1, FFTz2, FFTz3, FFT_z_Matrix)
  END SUBROUTINE destroy_FFT

  SUBROUTINE circulant_xsolve(rhs)
    real, dimension(xNc) :: rhs

    FFTx1 = rhs
    call dfftw_execute_dft(FFT_xplan, FFTx1, FFTx2)

    FFTx1 = FFTx2 / FFT_x_Matrix

    call dfftw_execute_dft(iFFT_xplan, FFTx1, FFTx2)

    rhs = REAL(FFTx2)/xNc
  END SUBROUTINE circulant_xsolve

  SUBROUTINE circulant_zsolve(rhs)
    real, dimension(zNc) :: rhs

    FFTz1 = rhs
    call dfftw_execute_dft(FFT_zplan, FFTz1, FFTz2)

    FFTz1 = FFTz2 / FFT_z_Matrix

    call dfftw_execute_dft(iFFT_zplan, FFTz1, FFTz2)

    rhs = REAL(FFTz2)/zNc
  END SUBROUTINE circulant_zsolve

  SUBROUTINE pinit_FFT(xrow, zrow)
    real, dimension(xNc) :: xrow
    real, dimension(zNc) :: zrow
    integer :: plan

    ! z dimension
    allocate(pFFT_z_Matrix(zNc))
    allocate(pFFTz1(zNc), pFFTz2(zNc), pFFTz3(zNc))
    call dfftw_plan_dft_1d(plan, zNc, pFFTz1, pFFTz2, FFTW_FORWARD, FFTW_ESTIMATE)
    pFFT_zplan = plan
    call dfftw_plan_dft_1d(plan, zNc, pFFTz1, pFFTz2, FFTW_BACKWARD, FFTW_ESTIMATE)
    piFFT_zplan = plan

    pFFTz1 = zrow
    call dfftw_execute_dft(pFFT_zplan, pFFTz1, pFFT_z_Matrix)

    ! x dimension
    allocate(pFFT_x_Matrix(xNc))
    allocate(pFFTx1(xNc), pFFTx2(xNc), pFFTx3(xNc))
    call dfftw_plan_dft_1d(plan, xNc, pFFTx1, pFFTx2, FFTW_FORWARD, FFTW_ESTIMATE)
    pFFT_xplan = plan
    call dfftw_plan_dft_1d(plan, xNc, pFFTx1, pFFTx2, FFTW_BACKWARD, FFTW_ESTIMATE)
    piFFT_xplan = plan

    pFFTx1 = xrow
    call dfftw_execute_dft(pFFT_xplan, pFFTx1, pFFT_x_Matrix)
  END SUBROUTINE pinit_FFT

  SUBROUTINE pdestroy_FFT
    deallocate(pFFTx1, pFFTx2, pFFTx3, pFFT_x_Matrix)
    deallocate(pFFTz1, pFFTz2, pFFTz3, pFFT_z_Matrix)
  END SUBROUTINE pdestroy_FFT

  SUBROUTINE pcirculant_xsolve(rhs)
    real, dimension(xNc) :: rhs

    pFFTx1 = rhs
    call dfftw_execute_dft(pFFT_xplan, pFFTx1, pFFTx2)

    pFFTx1 = pFFTx2 / pFFT_x_Matrix

    call dfftw_execute_dft(piFFT_xplan, pFFTx1, pFFTx2)

    rhs = REAL(pFFTx2)/xNc
  END SUBROUTINE pcirculant_xsolve

  SUBROUTINE pcirculant_zsolve(rhs)
    real, dimension(zNc) :: rhs

    pFFTz1 = rhs
    call dfftw_execute_dft(pFFT_zplan, pFFTz1, pFFTz2)

    pFFTz1 = pFFTz2 / pFFT_z_Matrix

    call dfftw_execute_dft(piFFT_zplan, pFFTz1, pFFTz2)

    rhs = REAL(pFFTz2)/zNc
  END SUBROUTINE pcirculant_zsolve

END MODULE col_fft
