CC = gfortran
CCc = gfortran -x f95 -c
# ifort -assume byterecl


TEST: .libW.o .libP.o .test.o
	$(CC) -o TEST .test.o .libW.o .libP.o

.libW.o:
	$(CCc) -o .libW.o vtylib__WAVECAR.f90

.libP.o:
	$(CCc) -o .libP.o vtylib__POSCAR.f90

.test.o:
	$(CCc) -o .test.o wavecarTEST.f90



.PHONY : clean
clean:
	-rm -f TEST *.mod *.o .*.o