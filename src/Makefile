#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

#****************************************************************************
# Edit library paths, compiler flags, and code options in Makefile.in
#****************************************************************************
include Makefile.in

#****************************************************************************
# Don't edit below this line
#****************************************************************************
CPP         = $(CPP_$(COMPILER))
FC          = $(FC_$(COMPILER))
FFLAGS      = $(FFLAGS_$(COMPILER))
FPSPEC      = $(FPSPEC_$(COMPILER))
TARGET_ARCH = $(TARGET_$(COMPILER))
MPIINC      = $(MPIINC_$(MPI))
MPILIB      = $(MPILIB_$(MPI))
MPILDFLAGS  = $(MPILD_$(MPI))

# OpenMP-MPI hybdridization
DO_HYBRIDyes = $($(COMPILER)_OMP)
DO_HYBRIDno  =

# CUDA FFT & RT library variables:
CUCC         = $(CUCC_$(CUCOMPILER))
CUFLAGS      = $(CUFLAGS_$(CUCOMPILER)) -DG$(PRECISION)_PRECISION
CULIBDIR     = $(CULIBDIR_$(CUCOMPILER))
CULIB        = $(CULIB_$(CUCOMPILER))
CUINC        = $(CUINC_$(CUCOMPILER))
CUGENCODE    = $(CUGENCODE_$(CUCOMPILER))
CUARCH       = $(CUARCH_$(CUCOMPILER))

GHOME        = $(PWD)
CONTRIB      = $(GHOME)/contrib
UTILS        = $(GHOME)/utils
IOLIB        = $(GHOME)/mpiio
BOUNDIR      = $(GHOME)/boundary

# FFTW library variables:
FFTWSINGLE  = f
FFTWDOUBLE  =
FFTWHYByes  = -lfftw3$(FFTW3$(PRECISION))_threads
FFTWHYBno   =
FFTWFLAGS   = $(FFTWHYB$(P_HYBRID)) -lfftw3$(FFTW$(PRECISION)) -lm 
FFTWLDFLAGS = $(FFTWFLAGS) $(DO_HYBRID$(P_HYBRID))
FFTWINC     = -I$(FFTWDIR)/include
FFTWLIBDIR  = -L$(FFTWDIR)/lib

CUDA_GHOSTno   =
CUDA_GHOSTyes  = -DDEF_GHOST_CUDA_ -DG_PPN_=$(PPN) $(CUBIND_NVIDIA)
CUDA_OBJSno    =
CUDA_OBJSyes   = $(CONTRIB)/nvidia/setaffinity_for_nvidia.o


# Common DEFs and Contrib Objects
GDEF           = $(CUDA_GHOST$(P_CUDA))
CONTRIBOBJS    = $(CUDA_OBJS$(P_CUDA))

# FFT Target variables derived from above FFT library quantities:
FFTTYPEno   = W
FFTTYPEyes  = CU
FFTPDIRno   = fftp
FFTPDIRyes  = fftp-cu
FFTTYPE     = $(FFTTYPE$(P_CUDA))
FFTTLIBDIR  = $(FFT$(FFTTYPE)LIBDIR)
FFTTLIB     = $(FFT$(FFTTYPE)LIB)
FFTTINC     = $(FFT$(FFTTYPE)INC)
FFTTLDFLAGS = $(FFT$(FFTTYPE)LDFLAGS)

FFTPDIR     = $(GHOME)/$(FFTPDIR$(P_CUDA))

FFTCULIBDIR = $(CULIBDIR)
FFTCULIB    = $(CULIB)
FFTCUINC    = $(CUINC) -I$(FFTPDIR)
FFTCULDFLAGS= $(CUFLAGS)

# Final flags and compiling commands
LDFLAGS     = $(MPILIB) $(FFTTLIBDIR) $(FFTTLIB) $(MPILDFLAGS) \
	      $(FFTTLDFLAGS)
INCLUDE     = $(MPIINC) $(FFTTINC) -I$(FFTPDIR) -I$(BOUNDIR) -I$(IOLIB) \
	      -I$(UTILS) -I$(GHOME)/include -I$(GHOME)
CPPFLAGS    = -D$(SOLVER)_SOL $(CUDA_GHOST$(P_CUDA))
PREPROCESS  = $(CPP) $(CPPFLAGS) $(GDEF)
COMPILE     = $(FC) $(FFLAGS) $(DO_HYBRID$(P_HYBRID)) $(INCLUDE) \
	      $(TARGET_ARCH) -c
SPECOMP     = $(FC) $(FPSPEC) $(DO_HYBRID$(P_HYBRID)) $(INCLUDE) \
	      $(TARGET_ARCH) -c
CCOMPILE    = $(CC) $(CFLAGS) $(INCLUDE) $(TARGET_ARCH) -c
CUCOMPILE   = $(CUCC) $(CUGENCODE) $(CUARCH) $(CUINC) -I$(FFTPDIR) \
	      $(CUFLAGS) -c
BIN         = ../bin

# Objects
HD_OBJS       = pseudospec_mod.o module_dns.o pseudospec_hd.o
BOUSS_OBJS    = $(HD_OBJS) pseudospec_phd.o
ROTBOUSS_OBJS = $(BOUSS_OBJS)
MHD_OBJS      = $(HD_OBJS) pseudospec_mhd.o
OBJS          = $($(SOLVER)_OBJS)
OBJSDIR       = pseudo
SUBDIRS       = mpiio fftp fftp-cu utils boundary $(CONTRIB)/nvidia
EXTRASDIR     = tools
GTIMER        = $(UTILS)/gtimer_mod.o
GUTILS        = $(UTILS)/gutils_mod.o 

all: edit $(CONTRIBOBJS) $(GTIMER) $(FFTPDIR)/fftp.o \
	$(GUTILS) $(IOLIB)/binary_io.o $(OBJS) $(BOUNDIR)/boundary_mod.o\
	$(SOLVER)



main: all

cancel: edit $(FFTPDIR)/fftp.o $(IOLIB)/binary_io.o $(OBJS)
	$(MAKE) $@$(SOLVER)
	$(MAKE)

boots: edit $(GTIMER) $(FFTPDIR)/fftp.o $(GUTILS) \
	$(IOLIB)/binary_io.o $(OBJS)
	$(PREPROCESS) -DG$(PRECISION)_PRECISION -I$(FFTPDIR) \
	   $(EXTRASDIR)/$(@).fpp -o $@.f90
	$(COMPILE) $@.f90 -o $@.o
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTPDIR)/*.o $(IOLIB)/*.o \
	   $(GTIMER) $(GUTILS) $(OBJS) $@.o $(LDFLAGS) -o $(BIN)/$@
	$(MAKE)

clean: 
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc* $(SOLVER).f90
	rm -f $(OBJSDIR)/pseudospec_mod.f90
	for dir in $(SUBDIRS); do \
            cd $${dir} && $(MAKE) clean && cd .. ; \
        done

dist:
	rm -f $(BIN)/{*MHD*,*CHANNEL*,*ROT*,*RAY*,boots,README*,*TESTS*} \
              *HD* *MHD* *ROT* *BOUSS* *CHANNEL* boots* *TESTS*
	rm -f $(BIN)/README_output.info
	$(MAKE) clean

edit:
	$(CPP) -DNX_=$(NX) -DNY_=$(NY) -DNZ_=$(NZ) -DORD_=$(ORD) \
    -DCX_=$(CX) -DCY_=$(CY) -DCZ_=$(CZ) -DOX_=$(OX) -DOY_=$(OY) -DOZ_=$(OZ)\
	   $(OBJSDIR)/pseudospec_mod.fpp -o $(OBJSDIR)/pseudospec_mod.f90
	cd $(FFTPDIR) && $(MAKE) edit CPP='$(CPP)' IKIND='$(IKIND)' \
	   CSIZE='$(CSIZE)' NSTRIP='$(NSTRIP)' PRECISION='$(PRECISION)' \
     NSTREAMS='$(NSTREAMS)' DO_HYBRID='DO_HYBRID$(P_HYBRID)'
	cd $(BOUNDIR) && $(MAKE) edit CPP='$(CPP)' SOLVER=$(SOLVER) \
	   DO_HYBRID='DO_HYBRID$(P_HYBRID)'
	cd ../
	echo "Output for this solver includes the following files:" >> \
	   $(BIN)/README_output.info
	echo "(note some output files may not be created, depending" >> \
	   $(BIN)/README_output.info
	echo "on specific solver and options being used)."  >> \
	   $(BIN)/README_output.info
	echo >> $(BIN)/README_output.info
	for file in $(OBJS:.o=.f90); do \
	 echo "From subroutines in" $$file "(see file for more details):"  >> \
	   $(BIN)/README_output.info ; \
	 perl -ne 'print if /Output files contain/ .. /Parameters/' \
	   $(OBJSDIR)/$$file >> $(BIN)/README_output.info ; \
	 echo >> $(BIN)/README_output.info ; \
	done
	for file in $(BOUNDIR)/*boundary.f90; do \
	 fname=$${file##*/}; \
	 if grep -q $$fname $(BOUNDIR)/boundary_mod.f90; then \
	  echo "From subroutines in" $$fname "(see file for more details):"  >> \
	    $(BIN)/README_output.info ; \
	  perl -ne 'print if /Output files contain/ .. /Parameters/' \
	    $$file >> $(BIN)/README_output.info ; \
	  echo >> $(BIN)/README_output.info ; \
	 fi \
	done
	perl -i -lpe 's/! Parameters//g' $(BIN)/README_output.info

$(SOLVER):
	$(PREPROCESS) specter.fpp -o $(SOLVER).f90
	$(COMPILE) $(SOLVER).f90
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(DO_HYBRID$(P_HYBRID)) \
           $(FFTPDIR)/*.o $(IOLIB)/*.o $(BOUNDIR)/*.o $(CONTRIBOBJS) $(GUTILS)\
           $(GTIMER) $(OBJS) $(SOLVER).o $(LDFLAGS) -o $(BIN)/$(SOLVER)

cancel$(SOLVER):
	$(PREPROCESS) $(EXTRASDIR)/$(@:$(SOLVER)=3D).fpp -o $@.f90
	$(COMPILE) $@.f90
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTP)/*.o $(IOLIB)/*.o $(OBJS) \
           $@.o $(LDFLAGS) -o $(BIN)/$@

$(FFTPDIR)/fftp.o:
	cd $(FFTPDIR) && $(MAKE) COMPILE='$(COMPILE)' CUCOMPILE='$(CUCOMPILE)'

$(BOUNDIR)/boundary_mod.o:
	cd $(BOUNDIR) && $(MAKE) COMPILE='$(COMPILE)'

$(IOLIB)/binary_io.o:
	cd $(IOLIB) && $(MAKE) COMPILE='$(COMPILE)'

$(UTILS)/gutils_mod.o:
	cd $(UTILS) && $(MAKE) COMPILE='$(COMPILE)' gutils_mod.o

$(UTILS)/gtimer_mod.o:
	cd $(UTILS) && $(MAKE) COMPILE='$(COMPILE)' gtimer_mod.o

$(CONTRIB)/nvidia/setaffinity_for_nvidia.o:
	cd $(CONTRIB)/nvidia && $(MAKE) CCOMPILE='$(CCOMPILE) -D_GNU_SOURCE'

$(OBJS):
	$(COMPILE) $(OBJSDIR)/$(@:.o=.f90)

tests:  edit $(CONTRIBOBJS) $(GTIMER) $(FFTPDIR)/fftp.o \
	   $(GUTILS) $(IOLIB)/binary_io.o $(OBJS)
	cd $(BOUNDIR) && $(MAKE) clean && $(MAKE) COMPILE='$(COMPILE)'
	$(PREPROCESS) -DTESTS_ specter.fpp -o TESTS.f90
	$(COMPILE) TESTS.f90
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(DO_HYBRID$(P_HYBRID)) \
           $(FFTPDIR)/*.o $(IOLIB)/*.o $(BOUNDIR)/*.o $(CONTRIBOBJS) $(GUTILS)\
           $(GTIMER) $(OBJS) TESTS.o $(LDFLAGS) -o $(BIN)/TESTS
