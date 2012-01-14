F90 = ifort -free -u -debug
INCLUDE = -I/home/dodek/src/lapack/include/intel64/lp64/
LDFLAGS = -L/opt/intel/composer_xe_2011_sp1.7.256/mkl/lib/intel64
LINK = -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

bairstow: bairstow.o polynomial.o
	${F90} bairstow.o polynomial.o -o bairstow ${LDFLAGS} ${LINK}

polynomial.o: polynomial.f90 
	${F90} -c polynomial.f90 -o polynomial.o 

bairstow.o: polynomial.o bairstow.f90 
	${F90} -c bairstow.f90 -o bairstow.o ${INCLUDE} 

clean:
	rm -f *.o *.mod bairstow
