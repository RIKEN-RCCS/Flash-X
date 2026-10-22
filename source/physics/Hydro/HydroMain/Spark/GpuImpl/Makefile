# Run generate.py first; execute this Makefile in its output directory.
# Supply RIKYU's OpenMP target compiler AND its architecture/offload flags.
FC = gfortran
FFLAGS = -O3 -fopenmp -ffree-line-length-none
TARGET_FLAGS =

all: periodic_batch
periodic_batch: hydro_types.f90 batch_kernels.f90 periodic_batch.F90
	$(FC) $(FFLAGS) $(TARGET_FLAGS) hydro_types.f90 batch_kernels.f90 periodic_batch.F90 -o $@
