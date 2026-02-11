FC := gfortran
FFLAGS := -O3

all: driver1 driver2

driver1: driver1.f90 mgh.o set_precision.o
	$(FC) $(FFLAGS) -o $@ $^

driver2: driver2.f90 mgh.o mgh_wrapper.o set_precision.o
	$(FC) $(FFLAGS) -o $@ $^

mgh.o: mgh.f90 set_precision.o
	$(FC) $(FFLAGS) -c -o $@ $<

mgh_wrapper.o: mgh_wrapper.f90 mgh.o
	$(FC) $(FFLAGS) -c -o $@ $<

set_precision.o: set_precision.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f mgh.o mgh.mod
	rm -f mgh_wrapper.o
	rm -f set_precision.o set_precision.mod
	rm -f driver1
	rm -f driver2
