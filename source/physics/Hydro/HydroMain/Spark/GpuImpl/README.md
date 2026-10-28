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

## Two-dimensional configurations

A separate 2D template, driver, and backend extend the same arithmetic to both
spatial directions. Existing 1D builds and command lines remain available.
Generate into a distinct directory and build on RIKYU:

```sh
python3 generate.py --dimension 2
cd build2d
make nvidia > build-opt.log 2>&1
# Or: make nvidia-debug > build-debug.log 2>&1
export OMP_TARGET_OFFLOAD=MANDATORY
./periodic_batch_opt 128 128 2 2 0.02 x_tube.csv 1
```

Arguments are global Nx, global Ny, blocks-x, blocks-y, final time, CSV path,
case number, and optional maximum steps (0 means no step cap). Defaults are
128x128 cells, 2x2 blocks, time 0.05, gpu2d.csv, case 1, and no step cap.
Both global extents must be at least eight and divisible by their block-grid
extents. The domain is [0,1)x[0,1), periodic in both axes; gamma=1.4 and SSPRK3
are unchanged. Cases are:

1. x-aligned Sod tube: high state where x<0.5, constant in y.
2. y-aligned tube: high state where y<0.5, constant in x.
3. diagonal periodic stripe: high state where modulo(x+y,1)<0.5.

The initial high/low states are the same as the 1D test, with both velocities
zero. Case 3 exercises both velocity components and multidimensional transport;
it is not compared to the ordinary 1D Sod exact profile.

Each stage computes x reconstruction/flux followed by y reconstruction/flux,
reusing minus/plus scratch sequentially. Both directional flux arrays survive
until the unsplit update. Corner guards copy directly from owning periodic
interiors, including when a block has fewer than three cells along an axis.
Positive x/y shared faces take their negative-face neighbor's flux. Integrated
fluxes include dy for x faces and dx for y faces (unit depth). Kernels cover
block,j,i with collapse(3). Shock neighborhoods and flattening use ndim=2;
CFL retains the original maximum directional rate and CFL=0.4. No alternate
multidimensional limiter, GPU arithmetic optimization, or AMR support is added.

### Validation before 2D timing

From GpuImpl, on the allocated GPU:

```sh
python3 validate_2d.py --exe build2d/periodic_batch_opt
```

This runs 64x64 cases 1/2/3 with block grids 1x1, 2x2, 4x1 and 1x4, compares
profiles across decompositions, checks x/y rotational agreement, and compares
96x64 diagonal runs with 1x1 and 3x2 blocks. Every run must pass the solver's
mass, both momentum components, and energy conservation check (1e-9). The
comparison tool additionally requires positive finite candidate states and
profile differences at most 1e-9. The output directory retains logs and profiles;
rerunning uses the same filenames and replaces them.

For an independent planar reference, run NewImpl's already-built CPU executable
with 64 cells and time 0.02, then pass that CSV to the validator:

```sh
# From NewImpl/build/standalone:
./periodic_shock_tube 64 1 0.02 reference64.csv 3
# From GpuImpl:
python3 validate_2d.py --reference-1d ../NewImpl/build/standalone/reference64.csv
```

Every row of the planar 2D solution is compared to the corresponding 1D profile;
the transverse velocity must also be below tolerance. This check assumes square
cells and equal global x/y resolution to match timestep selection. Without
--reference-1d the validator checks decomposition/rotation consistency and
conservation but has no independent numerical reference for the diagonal case.

Manual comparisons are also available:

```sh
python3 compare_2d.py validation2d/64x64-1x1-case3.csv validation2d/64x64-2x2-case3.csv
python3 compare_2d.py validation2d/64x64-1x1-case1.csv validation2d/64x64-1x1-case2.csv --transpose
```

CSV columns are x,y,density,velocity_x,velocity_y,pressure,internal_energy,
total_energy, in global y-then-x row order. The 1D compare.py/benchmark.py use
a different CLI/CSV schema and must not be used for these 2D executables.

For initial manual throughput runs, use the optional step cap:

```sh
# From build2d, after validation, without sanitizer:
./periodic_batch_opt 512 512 8 8 1.0 perf2d.csv 3 100
```

Record the actual step count, compiler flags, hardware, block dimensions and
resident/transfer-inclusive times; repeat configurations before drawing timing
conclusions. Large 2D grids increase memory and output volume. For local extents
nx,ny,nb, device buffers require approximately
8*nb*(34*(nx+6)*(ny+6)+20*(nx+1)*(ny+1)) bytes, excluding compiler temporaries
and runtime allocations. The host also retains backing allocations. Keep the
same caution about process startup and timing boundaries as in the 1D campaign.

Only source generation and Python comparison checks have been performed locally.
RIKYU compilation, sanitizer checks and numerical validation remain required for
this 2D extension.

## Macro definition file suffixes

The standalone storage definitions are `gpu_layout.ini-template` and
`gpu_layout_2d.ini-template`. They use NewImpl's custom macro syntax and are
loaded explicitly by generate.py. The suffix is intentional: selecting the AMR
unit also selects its GpuImpl parent, and Flash-X setup loads files ending in
`.ini` from selected directories with its native section-based macro parser.
Do not retain copies under the former `.ini` names in the source directory.
The AMR backend follows the same convention with `amr_layout.ini-template`.
