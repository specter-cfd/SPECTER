#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

MODS = boundary_mod.o

all: $(MODS)

mods: $(MODS)

$(MODS):
	$(COMPILE) $(@:.o=.f90)

edit:
	$(CPP) -D$(SOLVER)_SOL boundary_mod.fpp -o boundary_mod.f90

$(OBJS):
	$(COMPILE) $(@:.o=.f90)

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc*
	rm boundary_mod.f90
