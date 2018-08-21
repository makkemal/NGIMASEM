CFLAGS = -Wall -O3  -I /usr/include/spooles -L/usr/lib -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE
FFLAGS = -Wall -O3 -fopenmp


CC=cc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.14.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)


LIBS = -lspooles -larpack -lpthread -lm -lc -lblas

ccx_2.14: $(OCCXMAIN) ccx_2.14.a 
	./date.pl; $(CC) $(CFLAGS) -c ccx_2.14.c; $(FC) -fopenmp -Wall -O3 -o $@ $(OCCXMAIN) ccx_2.14.a $(LIBS)

ccx_2.14.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

.PHONY: clean

clean:
	rm -f ccx_2.14 *.o *.a
