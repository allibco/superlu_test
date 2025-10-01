FC       = mpif90
FFLAGS   = -O2 -fopenmp
LDFLAGS  =

# SuperLU_DIST install root
SUPERLU_DIST = $(NCAR_ROOT_SUPERLU_DIST)

# Include and library paths
INCLUDES = -I$(NCAR_INC_SUPERLU_DIST)
LIBDIRS  = -L$(SUPERLU_DIST)/lib

# Libraries to link
LIBS = -lsuperlu_dist -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

# Targets
all: small_superlu.exe

small_superlu.exe: small_superlu.o superlu_mod.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBDIRS) $(LIBS)

%.o: %.F90
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

clean:
	rm -f *.o *.mod small_superlu.exe
