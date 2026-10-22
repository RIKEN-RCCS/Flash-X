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

GPU expansion substitutes helper arguments directly instead of emitting
ASSOCIATE bindings. NVFORTRAN's unoptimized device compilation reported an
unsupported pghpf_ptr_assn_i8 runtime operation for those bindings. Array maps
and offsets now have named parameter arrays; spacing is initialized once per
host kernel call. Helper definitions in NewImpl are unchanged. Direct substitution
requires side-effect-free actual arguments and does not promise once-only
evaluation for expressions. Current device calls satisfy this restricted contract.
Run `python3 test_generation.py` for source-generation checks (no Fortran execution).

## RIKYU build and execution

No compilation or numerical/device execution has been performed locally. The
RIKYU compiler, accelerator architecture, OpenMP runtime and module environment
must be supplied by the site. An OpenMP target Fortran compiler supporting
IEEE intrinsics, reductions and persistent mappings is required.
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

## NVHPC diagnostic build on RIKYU

Regenerate first, then rebuild all generated sources with the same compiler:

```sh
python3 generate.py
cd build
nvfortran -O0 -mp=gpu -gpu=lineinfo -Minfo=mp hydro_types.f90 batch_kernels.f90 periodic_batch.F90 -o periodic_batch_debug
export OMP_TARGET_OFFLOAD=MANDATORY
compute-sanitizer --tool memcheck --log-file sanitizer.log ./periodic_batch_debug 32 4 0.001 debug.csv > debug.out 2>&1
```

This command requests line information without combining it with device debug
information from `-g`. Removing descriptor association addresses the reported
compilation barrier; a successful device run and sanitizer report are still
needed to establish whether the earlier illegal address is resolved.

## Optimized performance campaign on RIKYU

From GpuImpl, generate and build with the loaded NVHPC module:

```sh
python3 generate.py
cd build
make nvidia > build-opt.log 2>&1
cd ..
```

`make nvidia` always rebuilds `periodic_batch_opt` using `-O3 -mp=gpu`.
Supply `GPU_FLAGS` for site-specific architecture selection if required.
Keep build-opt.log with the measurements; the runner cannot reconstruct compiler
flags from an executable. No fast-math switches are added. Before timing, rerun
256 cells/4 blocks to t=0.05 with the optimized executable and compare its CSV to
the SSPRK3 NewImpl reference. Also compare its one-block output.

Run within an allocated GPU compute node, without Compute Sanitizer/profilers:

```sh
python3 benchmark.py --exe build/periodic_batch_opt
```

Defaults sweep 4096, 16384, and 65536 total cells with 1, 16, and 64 blocks,
100 timesteps, and 5 independent process repetitions per configuration. For a
short initial campaign:

```sh
python3 benchmark.py --cells 4096 --blocks 1 16 64 --steps 100 --repeats 3
```

The optional fifth solver argument is a maximum step count (0, the default,
means no step cap). The runner passes final time 1.0 and requires the requested
step count to complete before that time cap. Fixed step counts avoid the
quadratic increase in work that occurs when both resolution and the number of
CFL-limited steps grow in a fixed-physical-time sweep. These benchmark CSVs are
not necessarily at t=0.05 and must not be compared to the earlier reference
without matching their achieved time.

Campaign directories retain every profile/log, raw timings.csv and metadata.json
with GPU/compiler output and relevant environment. Only successful device runs
are accepted. The solver's conservation assertion remains enabled. Summary values
are median/min/max resident seconds per step. Cell updates/s counts full SSPRK3
timesteps, not individual stages and not guard-cell operations. These are
independent process repetitions, not an in-process warmup protocol.

At fixed total cells, changing block count measures decomposition and guard-work
effects. Increasing total cells measures workload scaling; it does not isolate
batch count because the block extent changes too. Compare cases with a fixed
cells-per-block ratio separately when studying batch growth. Start with the
small campaign; increase steps if durations are too short for stable timing.

Resident timing includes host orchestration, synchronous kernel launches,
reductions and stage execution, but excludes initial mapping and final solution
copy-back. Transfer-inclusive timing also includes those data-region operations
and the resident-time print. Both exclude program initialization, initial GPU
probe/context setup, final CSV writing and final conservation checks. Neither is
full process wall time or an isolated GPU kernel measurement.

The current NewImpl executable has no comparable internal timer. Do not report a
CPU/GPU speedup by comparing its full process wall time with GPU resident time.
For an initial controlled CPU comparison, the same optimized batch executable
can be run manually with OMP_TARGET_OFFLOAD=DISABLED and GPU_ALLOW_HOST=1;
record OMP_NUM_THREADS and CPU affinity. The campaign runner intentionally
rejects host execution. A comparison against the original Flash-X performance
requires its own matched numerical configuration and timing boundaries.
