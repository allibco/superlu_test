# module load superlu-dist/9.1.0
#module load cray-mpich/8.1.29 

# Compiler (Cray MPI wrapper)
FC      = ftn
FFLAGS  = -g
LDFLAGS =

# SuperLU_DIST install root
SUPERLU_DIST = $(NCAR_ROOT_SUPERLU_DIST)

# oneAPI MKL root (from module load oneapi)
#MKLROOT 

# Include and library paths
INCLUDES = -I$(SUPERLU_DIST)/include
LIBDIRS  = -L$(SUPERLU_DIST)/lib

# Static MKL link line (LP64, sequential threading)
MKL_LIBS = -Wl,--start-group \
           $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_sequential.a \
           $(MKLROOT)/lib/intel64/libmkl_core.a \
           -Wl,--end-group -lpthread -lm -ldl

# SuperLU_DIST + MKL
LIBS = -lsuperlu_dist $(MKL_LIBS)


OBJS = superlu_mod.o small_superlu.o superlupara.o

all: small_superlu.exe

small_superlu.exe: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

superlu_mod.o: superlu_mod.f90 superlupara.o
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

superlupara.o: superlupara.f90 
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

small_superlu.o: small_superlu.F90 superlu_mod.o superlupara.o
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

clean:
	rm -f *.o *.mod small_superlu.exe


