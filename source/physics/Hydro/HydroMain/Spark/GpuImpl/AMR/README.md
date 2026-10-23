# Flash-X AMR batch adapter

This optional backend connects the existing GPU arithmetic to Flash-X's
ordinary telescoping AMR advance in both 2D and 3D. Grid operations and reflux stay on the CPU.
The standalone 1D/2D GPU programs are unchanged.

## Supported initial configuration

- 2D or 3D Cartesian AMR, equal-sized full leaf blocks, variable-first host arrays.
- Telescoping RK2 and default WENO reconstruction; no rk3 subunit.
- Non-hybrid Gamma EOS; no gravity, MHD, species, mass scalars, or index reordering.
- Ordinary flux correction, not flux_correct_perLevel.
- One selected GPU per rank. Configure GPU visibility or OMP_DEFAULT_DEVICE in
  the job environment; the adapter does not assign devices to ranks automatically.

Unsupported compile-time combinations are rejected; runtime checks reject
non-Cartesian/non-telescoping/per-level execution and missing target devices.
These are implementation limits, not a general restriction on AMR.

## Generate and build on RIKYU

From the Flash-X root, before setup:

```sh
python3 source/physics/Hydro/HydroMain/Spark/GpuImpl/AMR/generate_amr.py
```

Add `gpuAmr=True` to your existing successful 2D or 3D Spark AMR/Gamma setup command.
Keep its normal machine, mesh, precision, and problem options. Do not select
`omp_offload`, `oacc`, `Milhoja`, `threadWithinBlock`, non-telescoping, or an
alternative reconstruction. The AMR Config requires the Gamma and AMR units,
defines HY_GPU_AMR and requests the OMP_OL library. This requests compiler/link
flags without selecting Spark's legacy per-block offload subunit.

Use the RIKYU site's nvfortran toolchain with OpenMP GPU compilation AND linking.
In that site's Makefile.h, verify FFLAGS_OMP_OL and LIB_OMP_OL include `-mp=gpu`
and any required GPU architecture flags. Regenerate setup after changing units.
Do not compile generated kernels without offload flags: the runtime device probe
will refuse CPU fallback. The new unit Makefile provides the module dependencies.

Use normal flash.par settings, including:

```text
flux_correct = .true.
flux_correct_perLevel = .false.
hy_gpuBatchSize = 32
```

Retain the usual gamma and hydrodynamic floors for the CPU comparison. Launch
Flash-X through the normal RIKYU allocation/launcher with
OMP_TARGET_OFFLOAD=MANDATORY. No local Flash-X compilation or solver execution
has been performed for this adapter; the generated units and full application
must be compiled and validated on RIKYU.

## Execution and ownership

Hydro.F90-mc still fills telescoping guards once on the host, zeros flux registers,
and later executes its existing flux communication and correction pass.
HY_GPU_AMR selects a replacement for only the ordinary per-block preparation/
advance/publication loop. The legacy per-level branch is not used.

The adapter iterates full leaf blocks without tiling. It calls the same host
shock detector used by Hydro_prepBlock, then packs eight hydrodynamic fields plus
the shock marker into contiguous block-major batches. Grid_tile_t descriptors
are copied by value; Grid data pointers are released immediately after packing
and reacquired for unpacking. The inspected Paramesh4 and AMReX descriptors store
block identity by value. Grid connectivity and descriptors never go to the GPU.

Batches flush on capacity, level/spacing changes, and end of iteration. Thus a
batch has one spacing and level. The final partial batch passes only its active
contiguous prefix. There is no level-dependent update schedule or level subcycling
in this adapter. Ranks with no leaf blocks launch no numerical batch but still
participate in the surrounding Hydro communication.

Host/device workspace is allocated once per Hydro invocation and reused for all
batches in that invocation. It is released before the CPU reflux pass, not retained
across timesteps/regridding. This avoids stale block associations but includes
allocation costs each Hydro call. No mutable global workspace is introduced.

Each batch transfers prepared state once, saves the complete reference state on
the device, and executes two RK stages. Stage one updates interior +/- NSTENCIL;
stage two updates the physical interior. Flattening and directional reconstruction
ranges include their stencil margins. Actual prepared guards must cover
NSTENCIL+3 cells for stage one; default WENO setup uses NSTENCIL=4 and NGUARD=8.
There are no GPU Grid calls, guard exchanges, or standalone periodic-face copies.

After each update the device reproduces the Eos.ini MODE_DENS_EI internal-energy
selection (eos_smalle/eos_eintSwitch), gamma-law pressure, and total-energy
reconciliation. The reference state remains immutable. Only pressure, gamma,
and energies are needed during subsequent numerical stages. After returning and
unpacking physical interiors, the ordinary CPU Eos_multiDim call populates all
other mapped EOS outputs such as temperature, GAME and optional entropy. This
is an additional final host EOS evaluation, not an intermediate state transfer.

Flux buffers contain sum(weight(stage)*flux_density), with stage-one reset and
only physical interior faces accumulated. They contain neither dt nor face area.
The host maps components to HY_MASS/HY_XMOM/.../HY_ENER and calls the same
Grid_putFluxData interface as Hydro_math.ini. Hydro's original CPU
Grid_communicateFluxes/getFluxCorrData/correctFluxes and correction EOS are retained.
Reconstruction workspace and instantaneous stage fluxes are not returned.

A numerical batch failure reports status, stage, level and batch size and aborts
before that batch's interiors are committed. First failing cell identification
is not yet implemented. Earlier batches and host shock diagnostics may already
have changed: this is an abort path, not a transactional timestep retry.

## Validation on RIKYU

Build a CPU comparison with identical mesh, precision, WENO helpers, RK2, gamma,
floors, and timestep configuration. Do not compare to the standalone SSPRK3
profiles. In particular, a CPU run using another reconstruction is not a direct
arithmetic-equivalence reference.

1. Run a fixed AMR hierarchy with at least one coarse-fine interface crossed by
   a shock, using a single MPI rank/GPU. Compare conserved totals after refluxing,
   solution profiles, and weighted flux buffers against the CPU run.
2. Repeat with batch sizes 1, 3 and 32, including a leaf count not divisible by
   batch size. Batch size 1 is a correctness case, not the intended performance mode.
3. Exercise multiple levels and all active interface orientations (x/y in 2D, x/y/z in 3D). Enable regridding
   only after the fixed hierarchy passes; then test MPI ranks, including empty ranks.
4. Run Compute Sanitizer on a small AMR application before timing. Verify all
   Grid/reflux operations remain host-side and large transfers occur per batch,
   not per stage.

Generation/interface checks can run without a Fortran compiler or GPU:

```sh
python3 source/physics/Hydro/HydroMain/Spark/GpuImpl/AMR/test_amr_generation.py
```

These checks do not establish device correctness or full Flash-X linkage. Legacy
Hydro scratch allocations are retained for compatibility with the CPU correction
pass; reducing their unused storage and persisting batch storage across Hydro
calls are subsequent optimizations. No CPU/GPU speedup is claimed for this
integration before its complete packing/transfer/EOS/reflux costs are measured.

## 3D extension and 2D compatibility

The same generated sources now compile with NDIM=2 or NDIM=3 from Simulation.h.
For a 3D run, regenerate the sources and use your normal **3D** Spark AMR/Gamma
setup with `gpuAmr=True`; do a fresh setup/build so incompatible 2D module files
are not reused. No extra generator dimension option is needed. Other restrictions
(RK2, Cartesian, WENO, ordinary reflux, host variable-first storage) still apply.

Compact cell batches have dimensions (i,j,k,variable,block), flattening scratch
(i,j,k,block), and directional face arrays (i,j,k,component,direction,block).
The backend definitions are in amr_layout.ini-template; the suffix keeps this
custom macro syntax out of Flash-X setup's standard INI parser. GPU loops collapse block/k/j/i;
2D has one unguarded k plane, while 3D includes the prepared z guards. Stage
bounds shrink in every active direction. Directional reconstruction remains
sequential so minus/plus scratch can be reused safely. Z flux divergence uses
the unchanged update helper. The adapter packs/unpacks all k planes and returns
z-face RK-weighted flux densities to the existing Grid_putFluxData/reflux path.

Start 3D validation with a small fixed AMR hierarchy and `hy_gpuBatchSize=1` or 2.
Repeat planar shock tests along x, y and z against matched CPU RK2 results, then
add a genuinely 3D initial condition. Check all momentum components and conserved
budgets after refluxing. Repeat with a partial batch, multiple refinement levels,
and finally regridding/MPI. Rerun the 2D AMR case too: its internal storage rank
has changed even though its physical equations and inactive dimension have not.

Memory grows with the fully guarded block volume. With real size R, batch capacity
B and dimension D, approximate device array storage is

```text
R * B * (33*Vguard + 5*D*(Vface_padded + Vface_interior))
```

Vguard includes guards only on active axes. The face volumes use a rectangular
allocation with one extra index in every active axis, even though a directional
face set needs the extra index only along its normal. Compiler temporaries and
runtime allocations are additional, and host backing arrays also consume memory.
Increase batch size only after checking available device and host memory.

Local tests perform source generation, setup-contract checks, and C preprocessing
for NDIM=2/3 (with synthetic headers), including balanced-loop and z-flux output
checks. These are not Fortran compilation or numerical/device validation.

## CPU testing on ghosha (gfortran + MPI)

For host correctness testing, add **both** variables to the same normal Spark
AMR/Gamma setup command:

```text
gpuAmr=True gpuAmrCpuTest=True
```

Perform a fresh setup and rebuild. This still selects the AMR batch adapter and
its generated numerical kernels. It defines HY_GPU_AMR_CPU_TEST and requests
OMP_HOST instead of OMP_OL. The ghosha site makefile now supplies
FFLAGS_OMP_HOST and LIB_OMP_HOST as `-fopenmp`. CFLAGS_OMP_HOST is empty
because ghosha's mpicc invokes Apple Clang; the batch OpenMP kernels are Fortran.
The MPI Fortran wrapper must invoke gfortran; inspect `mpif90 -show`.
Other sites can use this mode after providing equivalent host OpenMP flags for
their compiler. It does not request an OpenMP target/offload library.

Run from the configured object directory with host execution explicitly selected:

```sh
export OMP_TARGET_OFFLOAD=DISABLED
export OMP_NUM_THREADS=4
mpirun -n 1 ./flashx
```

Adjust thread/rank counts to the laptop's resources. The adapter prints a CPU-test
notice, accepts the initial device, and rejects accidental GPU execution in this
mode. This exercises the batched telescoping/update/EOS/flux-return path on CPU;
it does not validate GPU mappings or GPU performance. All normal AMR, geometry,
EOS and RK restrictions still apply. No environment-variable bypass is added to
the production configuration: without gpuAmrCpuTest=True the OMP_OL requirement
and mandatory target-device check remain unchanged.

The earlier request not to compile or run Flash-X locally is retained: local
verification covers setup selection, source generation and preprocessing in
both dimensions and both modes, not a compiled application run.

## Riemann failure diagnostics

HY_DEGENERATE (status 3) means the selected Riemann algebra encountered an
exactly zero denominator. On a face-kernel failure, the driver copies the current
stage/reconstruction arrays to the host and serially replays the same expanded
face calculation. It prints the first reproduced failing batch slot/cell/direction,
calculation name, repaired reconstructed left/right states, wave speeds, contact
speed and star pressure. These failure-only transfers do not occur on successful
steps. Capture standard output as well as the Flash-X log when reporting failures.
No HLLE fallback or arithmetic modification has been introduced. Host replay can
fail to reproduce an accelerator-only/compiler-dependent issue; this is explicitly
reported, rather than treating the replay as a successful timestep.

## HDF5-independent cell verification dumps

Optional rank-local CSV output is now available after the entire Hydro advance
and CPU reflux/correction-EOS pass. Regenerate setup and rebuild to register the
new runtime parameters and dump module. It is enabled for this AMR backend in
both CPU-testing and GPU builds, and requires no HDF5 calls or libraries itself.

For the fixed-level comparison, use the same problem, block dimensions, gamma,
CFL=0.4, precision and refinement level on both machines. For example:

```text
lrefine_min = 4
lrefine_max = 4
hy_gpuDumpInterval = 10
hy_gpuDumpPrefix = "laptop_verify"
```

Set the prefix to `rikyu_verify` for the GPU run. The level above is an example;
choose one that fits both machines. Interval zero (the default) disables all
dumps. A positive interval writes whenever the current Driver_getNStep index is
divisible by the interval. Each snapshot contains the post-Hydro state at the
passed timeEndAdv and dt; it is not an initial-condition dump or a complete
checkpoint. Example filenames are:

```text
laptop_verify_step0000000010_rank000000.csv
laptop_verify_step0000000010_rank000001.csv
```

Only owned leaf interiors are written, excluding guards and covered parents.
Columns include driver step/time/dt, refinement level, Grid i/j/k indices,
physical cell centers, widths, density, three velocities, pressure, and specific
internal/total energies. The inactive z coordinate is zero and dz is one in 2D
(unit depth). The console reports rank-local mass/momenta/energy; the comparison
script computes corresponding global totals by combining all rank files.
These are solution dumps, not flux-buffer dumps. Fixed-level tests do not exercise
coarse-fine reflux interfaces.

Files at the same prefix/step/rank path are replaced; other steps remain. Prefixes
may include a directory, but it must already exist on every rank. Disable output
when measuring performance, since dumping is included in the Hydro timer.

Copy the snapshots to one machine and compare exactly one step at a time:

```sh
python3 source/physics/Hydro/HydroMain/Spark/GpuImpl/AMR/compare_dumps.py \
  --reference 'cpu/laptop_verify_step0000000010_rank*.csv' \
  --candidate 'gpu/rikyu_verify_step0000000010_rank*.csv' \
  --report comparison-step10.json
```

Quoted patterns are expanded by the script. It matches records by refinement
level and rounded physical coordinates, independent of rank/file order and local
Grid indices. Defaults are 12 coordinate decimal places, cell-value absolute
tolerance 1e-10 and relative tolerance 1e-8. A cell quantity passes when its absolute
difference is at most atol+rtol*max(abs(reference),abs(candidate)). Time and dt must
agree within 1e-12 and step indices must match. Widths are also checked.

The JSON report includes maximum absolute differences by field, failed-cell
counts, global reference/candidate conserved totals and their differences.
Global total differences are reported, not independently thresholded. Missing
cells, coordinate-key collisions, mixed snapshots, invalid data or failed cell
tolerances produce a nonzero exit status. Use --coordinate-digits only when
appropriate for the domain scale/finest spacing. Use --time-atol/--atol/--rtol
with a numerical justification, not to hide a different timestep history or mesh.
The script accepts different MPI rank counts if the physical cell set agrees.

This output does not remove any separate HDF5 dependency selected by Flash-X's
normal IO configuration. That must still be disabled/replaced in the RIKYU
setup independently. These CSVs cannot be used to restart Flash-X.

Synthetic comparison tests (no solver execution):

```sh
python3 source/physics/Hydro/HydroMain/Spark/GpuImpl/AMR/test_compare_dumps.py
```
