
FFLAGS=-O3  

PoreProcess : a.o b.o Modules.o secnds.o
	ftn $(FFLAGS) -o $@ PoreProcess.F90 a.o b.o Modules.o secnds.o $(LDFLAGS)

b.o : a.o
a.o : Modules.o

%.o: %.F90
	ftn $(FFLAGS) -o $@ -c $< $(LDFLAGS) 
%.o: %.f90
	ftn $(FFLAGS) -o $@ -c $< $(LDFLAGS)

clean :
	rm -f a.o b.o Modules.o secnds.o PoreProcess PoreProcess+pat
	rm -f *.mod buckingham_m.mod  constant_m.mod  f2dlocald_m.mod  flexiclinfo_m.mod  fluctuation_m.mod  fluid_solid_m.mod  mcsetting_m.mod  physical_m.mod  porefigure_m.mod  subbox_m.mod
	rm -rf PoreProcess.pl

PoreProcess+pat : PoreProcess
	pat_build $^

