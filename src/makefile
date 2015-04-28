#FFLAGS=-O3 -hnoomp
FFLAGS=-O3
#FFLAGS=-h pl=PoreProcess.pl -h profile_generate
#FFLAGS=-h pl=PoreProcess.pl -hnoomp -h profile_generate

PoreProcess : a.o b.o Modules.o secnds.o
	ftn -o $@ $(FFLAGS) PoreProcess.F90 a.o b.o Modules.o secnds.o

b.o : a.o
a.o : Modules.o

%.o: %.F90
	ftn -o $@ $(FFLAGS) -c $<

%.o: %.f90
	ftn -o $@ $(FFLAGS) -c $<

clean :
	rm -f a.o b.o Modules.o secnds.o PoreProcess PoreProcess+pat
	rm -f buckingham_m.mod  constant_m.mod  f2dlocald_m.mod  flexiclinfo_m.mod  fluctuation_m.mod  fluid_solid_m.mod  mcsetting_m.mod  physical_m.mod  porefigure_m.mod  subbox_m.mod
	rm -rf PoreProcess.pl

PoreProcess+pat : PoreProcess
	pat_build $^

