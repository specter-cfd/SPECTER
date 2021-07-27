#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

MODS   = fprecision_mod.o fftp_mod.o cuda_bindings_mod.o fcgram_mod.o
FOBJS  =  fftp.o
CUOBJS = cuda_utils.o cuwrappers.o cuda_math.o

all: $(FOBJS) $(CUOBJS)

mods: $(MODS)

$(FOBJS): mods
	$(COMPILE) $(@:.o=.f90)

$(CUOBJS):
	$(CUCOMPILE) $(@:.o=.cu)

$(MODS):
	    $(COMPILE) $(@:.o=.f90)

edit:
	$(CPP) -DG$(PRECISION)_PRECISION fprecision_mod.fpp -o \
           fprecision_mod.f90
	$(CPP) -DIKIND_=$(IKIND) -DCSIZE_=$(CSIZE) -DNSTRIP_=$(NSTRIP) \
	   -DNSTREAMS_=$(NSTREAMS) -DG$(PRECISION)_PRECISION fftp_mod.fpp \
           -o fftp_mod.f90
	$(CPP) -DG$(PRECISION)_PRECISION fftp.fpp -o fftp.f90

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc*
	rm -f fprecision_mod.f90 fftp_mod.f90 fftp.f90
