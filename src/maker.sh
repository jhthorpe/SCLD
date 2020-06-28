gfortran -c P.f90
gfortran -c TS.f90 P.o
gfortran -o SCLD SCLD.f90 TS.o P.o
cp SCLD bin
