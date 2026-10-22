# Batched GPU arithmetic baseline

This is a new OpenMP target execution layer. NewImpl is unchanged. Arithmetic
helpers are loaded directly from NewImpl/hydro_helpers.ini, update_helpers.ini,
diagnostic_helpers.ini, and driver_helpers.ini. The generator also extracts the reference face
solver's executable sequence, retaining HLLC star states, hybrid HLLE choice,
floors, finite checks, and original-state viscosity. Allocation helpers and
block-traversal helpers are not invoked on the device.

Initial scope: serial host orchestration, one GPU, uniform 1D periodic Cartesian
blocks, ideal gas gamma=1.4, no gravity/species/MHD, SSPRK3. This is a standalone
baseline, not yet a Flash-X adapter, multidimensional solver, or AMR/MPI driver.
The backend is spatial-index first, then variable, then block. Each collapsed
block/spatial loop covers the complete batch. Full minus/plus reconstruction
arrays remain; no stencil tiling, reconstruction/flux fusion, or tuned arithmetic
has been introduced. Small private arrays and original helper control flow remain.

## RIKYU build and execution

No compilation or numerical/device execution has been performed locally. The
RIKYU compiler, accelerator architecture, OpenMP runtime and module environment
must be supplied by the site. An OpenMP target Fortran compiler supporting
ASSOCIATE, IEEE intrinsics, reductions and persistent mappings is required.
The existing expander produces Fortran 2003 helper constructs; the example
program additionally uses NEWUNIT (Fortran 2008).

From this directory:

```sh
python3 generate.py --output build
cd build
make FC="YOUR_FORTRAN_COMPILER" FFLAGS="YOUR_OPTIMIZATION_AND_OPENMP_FLAGS" TARGET_FLAGS="YOUR_DEVICE_FLAGS"
export OMP_TARGET_OFFLOAD=MANDATORY
./periodic_batch 256 4 0.05 gpu_sod.csv
```

Do not use the Makefile's host defaults for performance measurements. Use the
site compiler's GPU target and architecture flags. The executable independently
checks omp_is_initial_device and rejects silent CPU fallback. GPU_ALLOW_HOST=1
permits deliberate host correctness development only. Argument order is global
cell count, block count, final time, output CSV. Cell count must be at least 8
and divisible by positive block count. CSV output replaces the selected file.
Periodicity introduces a reversed discontinuity at x=0 as well as the Sod
interface at x=0.5.

Solution arrays enter one target-data region for the full time integration.
Reference/reconstruction/flux/acoustic buffers are device allocations. Kernel
routines use map(alloc:...) to reuse those mappings. Only scalar reductions,
small inputs and status values need host interaction between kernels; the
solution returns at the end. Generated declaration names are explicitly private
per iteration. All launches are synchronous initially. Error reductions report
status, not the first failing cell. Failed steps are not retried.

Separate launches perform guards, acoustic preparation, shock marking, CFL
reduction, reference copy, flattening, reconstruction, face solve, periodic face
reconciliation, update, integrated-flux accumulation, and ideal-gas EOS recovery.
Directional scratch reuse remains unchanged because the first provider is 1D.
Integrated transfers use unit face area and RK stage weights; currently they
are retained on device but not exported as a flux register.

## Validation on RIKYU

First build NewImpl's standalone reference with the same compiler family and
precision (adjust its GNU-specific build flags if necessary). Generate reference
output using SSPRK3, the same cells, and the same end time:

```sh
# In NewImpl, after building its reference executable:
build/standalone/periodic_shock_tube 256 4 0.05 reference.csv 3
# In GpuImpl, after running the GPU executable:
python3 compare.py ../NewImpl/reference.csv build/gpu_sod.csv
```

Repeat with 1 and 4 blocks and compare both GPU profiles. Increase batch size and
cell count for throughput measurements after correctness passes. The executable
checks mass/momentum/energy conservation to 1e-9. compare.py checks finite,
positive states and profile agreement (default absolute tolerance 1e-9). GPU
floating-point contraction may require an explained tolerance change; establish
numerical agreement before performance interpretation. Invalid or extreme states
are not comprehensively tested by this first driver.

The executable reports elapsed wall time for resident stepping and for the
whole data region including initial/final transfers. Synchronous launches make
these timings inclusive of completion, but they are not individual kernel
profiling. Use RIKYU's profiler to verify mappings do not transfer full arrays
at each launch and to measure kernel time, memory traffic, register spills,
occupancy, and launch overhead. Match compiler optimization and runtime checks
when comparing CPU/reference timings. The NewImpl default builder enables
runtime checks, whereas this Makefile defaults to optimization without checks.
Repeat runs and report hardware/compiler/flags and problem sizes. No speedup or
performance penalty has yet been measured.
