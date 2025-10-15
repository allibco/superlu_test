# module load superlu-dist/9.1.0
#module load cray-mpich/8.1.29 

# Compiler (Cray MPI wrapper)
FC      = ftn
FFLAGS  = -g -traceback -check bounds -check uninit 
LDFLAGS =
CC = cc

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
LIBS = -lsuperlu_dist_fortran -lsuperlu_dist $(MKL_LIBS)

OBJS = superlu_mod.o small_superlu.o superlupara.o superlu_bindings.o f_pdvsmv.o

all: small_superlu.exe

small_superlu.exe: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

superlu_mod.o: superlu_mod.f90 superlupara.o
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

superlupara.o: superlupara.f90 
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

superlu_bindings.o: superlu_bindings.F90
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

small_superlu.o: small_superlu.F90 superlu_mod.o superlupara.o f_pdvsmv.o
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

f_pdgsmv.o: f_pdgsmv.c
	$(CC) -c $(INCLUDES) $< -o $@


clean:
	rm -f *.o *.mod small_superlu.exe


