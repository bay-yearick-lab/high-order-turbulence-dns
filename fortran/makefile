# nsbs: compressible Navier-Stokes solver with B-spline collocation and B-spline filtering
#
# Requirements: Intel Fortran (ifort/ifx) with MKL (ILP64 BLAS/LAPACK: DGBTRF, DGBTRS,
# DGBMV, DGEMV, DGEMM) and FFTW3 (Fortran 2003 interface, fftw3.f03).
# Reals are 8-byte and default integers are 8-byte (-r8 -i8), hence the ILP64 MKL link line.
#
#   make            build ./mixc
#   ./mixc          run with the parameters in mix.in

FORT    = ifort
NICE    = -r8 -i8
FFTW_INC ?= /usr/include
LIB     = ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a   \
          ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a \
          -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a \
                            ${MKLROOT}/lib/intel64/libmkl_intel_thread.a \
                            ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group \
          -liomp5 -lm -ldl -lfftw3
OPTS    = -O3 -heap-arrays -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -I$(FFTW_INC)
FC_CMD  = $(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)

OBJECTS = prms.o data.o rk.o rhs.o main.o ders.o filter.o ic.o out.o \
          bc.o mesh.o bspline.o gausstable.o quad.o mms.o \
          col_fft.o col_operator.o col_linalg.o test_collocation.o \
          gcv_objectives.o praxis.o

mixc: $(OBJECTS) makefile
	$(FORT) -o mixc $(PROF) $(OBJECTS) $(LIB)

main.o: main.f90 prms.o data.o rk.o ders.o ic.o out.o mesh.o filter.o \
	bspline.o gausstable.o quad.o mms.o col_operator.o test_collocation.o \
	col_fft.o col_linalg.o makefile
	$(FC_CMD) main.f90
rk.o: rk.f90 prms.o data.o bc.o rhs.o col_linalg.o makefile
	$(FC_CMD) rk.f90
rhs.o: rhs.f90 prms.o data.o ders.o bc.o makefile
	$(FC_CMD) rhs.f90
out.o: out.f90 prms.o ders.o mesh.o data.o makefile
	$(FC_CMD) out.f90
mesh.o: mesh.f90 prms.o bspline.o makefile
	$(FC_CMD) mesh.f90
ic.o: ic.f90 data.o prms.o mesh.o makefile
	$(FC_CMD) ic.f90
bspline.o: bspline.f90 prms.o makefile
	$(FC_CMD) bspline.f90
bc.o: bc.f90 prms.o mesh.o makefile
	$(FC_CMD) bc.f90
filter.o: filter.f90 prms.o mesh.o makefile
	$(FC_CMD) filter.f90
data.o: data.f90 prms.o makefile
	$(FC_CMD) data.f90
ders.o: ders.f90 prms.o mesh.o bspline.o makefile
	$(FC_CMD) ders.f90
prms.o: prms.f90 makefile
	$(FC_CMD) prms.f90
gausstable.o: gausstable.f90 makefile
	$(FC_CMD) gausstable.f90
quad.o: quad.f90 prms.o mms.o gausstable.o bspline.o makefile
	$(FC_CMD) quad.f90
mms.o: mms.f90 makefile
	$(FC_CMD) mms.f90
col_operator.o: col_operator.f90 prms.o bspline.o quad.o col_fft.o makefile
	$(FC_CMD) col_operator.f90
col_fft.o: col_fft.f90 bspline.o prms.o makefile
	$(FC_CMD) col_fft.f90
col_linalg.o: col_linalg.f90 prms.o col_operator.o col_fft.o makefile
	$(FC_CMD) col_linalg.f90
test_collocation.o: test_collocation.f90 prms.o bspline.o col_fft.o \
	col_operator.o col_linalg.o rk.o praxis.o makefile
	$(FC_CMD) test_collocation.f90
gcv_objectives.o: gcv_objectives.f90 col_linalg.o col_operator.o makefile
	$(FC_CMD) gcv_objectives.f90
praxis.o: third_party/praxis.f90 makefile
	$(FC_CMD) third_party/praxis.f90

clean:
	rm -f mixc *.o *.mod *~ fort.* mix.grid.* mix.filt.* col*.dat knt.dat mix.sig.*
