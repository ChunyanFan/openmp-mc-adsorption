#FFLAGS= -O3 -hnoomp
FFLAGS=-O3 
# -O vector3 #--> 50.776
# -O ipa5 #--> 51.27
#-hautothread -O vector3 -O scalar3 -O thread3 -homp -O ipa5 -O modinline 
#-O ipafrom=Modules.F90
#-hautothread 
#Following 2 FFLAGS are for REVEAL
#FFLAGS=-hnoomp -h profile_generate 
#FFLAGS= -O3 -hnoomp -h pl=PoreProcess.pl 
#FFLAGS= -hnoomp -h profile_generate -h pl=PoreProcess.pl 
#Following FFLAGS and 2 LDFLAGS (shared and then static libraries) are for MAP library- 
#FFLAGS= -G2 -O3
#LDFLAGS= -Wl,@/group/partner1032/ddeeptimahanti/src/allinea-profiler.ld

#LDFLAGS=-dynamic -L/group/partner1032/ddeeptimahanti/src -lmap-sampler-pmpi -lmap-sampler -Wl,--eh-frame-hdr
#LDFLAGS=-L/group/partner1032/ddeeptimahanti/src -lmap-sampler-pmpi -Wl,--undefined,allinea_init_sampler_now -lmap-sampler -lstdc++ -lgcc_eh -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -Wl,--eh-frame-hdr

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
	rm -f buckingham_m.mod  constant_m.mod  f2dlocald_m.mod  flexiclinfo_m.mod  fluctuation_m.mod  fluid_solid_m.mod  mcsetting_m.mod  physical_m.mod  porefigure_m.mod  subbox_m.mod
	rm -rf PoreProcess.pl

PoreProcess+pat : PoreProcess
	pat_build $^

