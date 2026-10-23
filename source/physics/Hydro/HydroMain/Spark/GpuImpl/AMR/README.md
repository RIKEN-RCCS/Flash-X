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
