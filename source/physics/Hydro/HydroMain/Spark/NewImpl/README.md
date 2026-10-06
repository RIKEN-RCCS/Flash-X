# Single-stage hydrodynamic face-flux implementation

`hy_getFaceFlux.F90-mc` implements the active host mathematics in
`docs/designDocs/face_flux_solver_algorithm.pdf`: flattening, componentwise WENO-Z,
paired monotonicity, state floors, Davis speeds, HLLC contact and star states,
optional shock-triggered HLLE, original-state artificial viscosity, and output.
Magnetic fields, GLM, species, other mass scalars, and RK stage indexing are absent.
No EOS, boundary construction, or solution update is performed.

## Expand and compile

From this `NewImpl` directory:

```sh
python3 macro_expand.py --definitions hydro_layout.ini hydro_helpers.ini \
  --input hy_getFaceFlux.F90-mc --output build/hy_getFaceFlux.F90
cd build
gfortran -std=f2003 -Wall -Wextra -fcheck=all -c hy_getFaceFlux.F90
```

Run the macro and numerical tests from this directory:

```sh
python3 test_newimpl.py
```

Tests compile and execute the same solver against three backends: variable-first
arrays, spatial-first arrays, and spatially indexed derived records with allocated
component vectors. They use permuted, noncontiguous variable maps and two blocks.
Checks include dimensional reductions, physical fluxes, WENO linear reconstruction,
flattening, floors, a stationary contact, Sod HLLC/HLLE fluxes, hybrid viscosity,
viscosity after the tiny-state zero flux, and explicit error returns. Compilation
enables bounds checks and traps invalid arithmetic, division by zero, and overflow.
A missing gfortran skips compiled tests; the expander needs only Python 3.10+.

## Macro contract

The expander is a standalone adaptation of the repository's
`source/physics/Hydro/HydroMain/unsplit/macro_expand.py`, including a fix to avoid
substituting formal identifiers twice inside actual argument text.
Definitions use unsectioned INI entries with `args=` and `definition=`, matching
`hy_newMacros.ini`. Comments begin with `;` or `#`. Use `${argument}` placeholders.
When a definition contains `${...}`, only explicit placeholders are substituted;
plain identifiers are left alone. Definitions without placeholders retain legacy
bare-identifier substitution. This prevents arguments such as `block` from
replacing Fortran keywords.

Invocations are `@name@` or `@name(arguments)@`; nested macros and parenthesized
index expressions are supported. Undefined names, wrong argument counts, malformed
invocations, and recursive definitions are errors. Later files override earlier
ones unless `--reject-duplicates` is specified. CPP directive lines are preserved.
Computational helpers use `scope=block` as a logical scoping property. The file
expander targets Fortran 2003: it hoists uniquely named local declarations into
the enclosing procedure and emits a named single-pass `do` construct, rather than
a Fortran 2008 `block`. Once-evaluated `associate` argument bindings are retained,
and `return` becomes a helper-local `exit`. Fortran 95 is not supported because
`associate` and the supplied allocatable-component backends require Fortran 2003. Macro authors write ordinary short names:

```ini
example
args=q,u
scope=block
declarations=
  real(rk) :: kinetic
definition=
  kinetic=0.5_rk*q(QRHO)*sum(q(QUX:QUZ)**2)
  u(UENER)=kinetic+q(QEINT)
```

The `declarations=` section is optional. Local variable declarations must use `::`
and one complete declaration statement per line; scalar and array declarations
are supported, while initialized declarations are rejected. Intrinsic `use`
statements are allowed in that section. Local names must differ from formal
argument names. Renaming is case-insensitive and leaves component names, macro
names, quoted strings, and comments intact. Nested helpers receive independently
scoped aliases and temporaries. Spatial bounds remain explicit arguments.
Reserve generated `macro_scope_` names. The old `__macro_scope__` marker remains
supported for unscoped legacy definitions, but is no longer used in helpers.
The file writer wraps expanded `.f90` code to standard free-form line lengths.
Expansion is textual; reserve `@` for invocations outside CPP directives.

All spatial storage types, scalar references, allocation, and scratch release are
backend macros. Numerical kernels use small local vectors with fixed semantic
ordering; these do not prescribe the mesh layout. There are no mesh array slices
or assumed contiguous mesh-variable ranges.

For example:

```text
@cell_ref(data,inputMap(QRHO),block,i-i1,j-j1,k-k1)@
@face_ref(flux,outputMap(UENER),block,dir,i,j,k)@
```

To select the spatial-first or record backend, append its INI file after
`hydro_layout.ini` in `--definitions`. To integrate another storage system,
replace `field_types`, `cell_ref`, `face_ref`, `allocate_cell`, `allocate_flux`, and
`release_cell`/`release_flux`. The solver body is unchanged. The supplied backends are examples;
the per-cell allocation record backend demonstrates flexibility, not performance.

## Declaration insertion markers

Each procedure expanding scoped macros must contain these markers:

```fortran
subroutine example(...)
  !$macro_imports
  implicit none
  ! Existing argument and local declarations go here.
  !$macro_declarations
  ! Executable code and helper macro invocations go here.
end subroutine
```

The first marker receives required intrinsic `use` statements; place it before
`implicit none`. The second receives the helpers' uniquely named declarations;
place it after existing declarations and before executable statements. Markers
are removed during file expansion. Each contained procedure needs its own markers
if it expands scoped helpers. The solver and test callers already include them.
Local array bounds must be valid in the enclosing declaration section; bounds
that depend on generated executable argument aliases are not supported.
`expand_text` performs macro substitution for fragments; `expand_file` additionally
hoists declarations, lowers logical scopes, and wraps Fortran lines.

## Caller contract

The module has no mutable global variables. `use` imports only intrinsic kind and
IEEE procedures. Configuration, fields, block identity, bounds, maps, workspace,
and diagnostics are explicit arguments. The parameter and workspace objects are
caller-owned. The implementation is serial; local temporaries are per invocation.
Separate simultaneous calls need separate workspaces and disjoint output regions.

`inputMap(7)` identifies density, x/y/z velocity, pressure, gammaC, and **specific**
internal energy, in that order. Stencil gathering forms rho times specific energy
before reconstruction. `outputMap(5)` identifies mass, x/y/z momentum, and total
energy flux. All three velocity components are retained at every dimensionality.
`shockVariable` identifies a supplied shock marker when hybrid mode and
`shock_marker_present` are both enabled; otherwise it is not accessed.

Allocate input storage with enough variables for all maps and the optional marker,
and enough blocks for the requested ID. `availableLo/Hi` must describe genuinely
allocated, initialized cell data. The driver validates the declared neighborhood,
not the backend's actual storage extents or map upper bounds.

For requested cell bounds `lo:hi`, allocate workspace over `lo-1:hi+1` in active
directions (inactive extents can remain unchanged). Expand `@hy_allocateWorkspace(work,lo,hi,nblock)@` once and reuse the workspace.
Allocation and release helpers are macros expanded in caller code. Allocation
macros release previously allocated storage before allocating new storage. Input guard cells must cover `lo-3:hi+3` in active directions.
Allocate output over `lo:hi+1` in active directions and all five mapped components.
Face `(i,j,k,dir)` lies between `(i,j,k)-unit(dir)` and `(i,j,k)`; only its normal
axis extends through `hi(dir)+1`. Fluxes are positive toward the latter cell, per
unit area, and have no time-step or area scaling. Other output entries are untouched.

The workspace stores one flattening scalar and two seven-component reconstructed
states per workspace cell/block. Flattening is reused across directions. The two
reconstruction fields are overwritten for each direction. Five-point stencils,
WENO triples, primitive/conservative/physical-flux vectors, and two five-component
star-state vectors are local scratch. No full prepared-state or flux scratch field
is needed. Disabled flattening fills its work field with one.

`hy_getFaceFlux` is the only subroutine in the implementation module. All helpers,
including allocation/release, are defined in `hydro_helpers.ini` and invoked as
`@hy_helper(arguments)@`. The expander hoists local scratch declarations and generates a named single-pass
`do` plus `associate` aliases that evaluate argument expressions once. A helper's
early return is a named `exit` from its single-pass loop; status checks in the top-level
routine determine whether the full computation should stop. The helper's logical
scope imposes no spatial bounds. Loop helpers receive `lo/hi` explicitly;
single-cell and single-face helpers receive explicit indices.

Arguments that the macro assigns must be writable variables or array designators.
Input/output arguments must not alias each other in ways that change input values
before their last use. All macro arguments must have the kind, rank, and shape
expected by the original numerical interfaces. Macro callers import the public
types/constants from `newimpl_hydro` and supply the `rk` kind; IEEE checks import
only intrinsic procedures, hoisted into the enclosing procedure.

Each PDF step has a named kernel; Steps 8 and 9 are
`hy_computeHydroContact` and `hy_computeHydroStarState`. Conservative conversion is
shared by the Riemann calculation and original-state viscosity. Step 12 only emits
hydrodynamic fluxes because there are no advected fractions.

## Error and branch behavior

`HY_OK`, `HY_BAD_INPUT`, `HY_BAD_STATE`, and `HY_DEGENERATE` are public status codes.
The driver rejects invalid neighborhoods/configuration, nonfinite face states,
invalid sound speeds, and zero denominators. This makes the PDF's unchecked error
path explicit; no additional numerical fallback is silently selected.
`failedFace=[i,j,k,dir]` identifies face-loop failures, and is zero on success or
failures before the face loop. On failure, already emitted fluxes and the signal
speed reduction are partial; the caller must reject that computation.

Floors affect density and pressure only. Reconstructed internal-energy density,
gammaC, and velocity are not reset by an EOS. Tiny positive repaired density or
pressure yields zero Riemann flux, followed by the usual viscosity correction.
Hybrid HLLE faces skip unused HLLC intermediates. Supersonic HLLC faces similarly
return the relevant physical flux directly. Degenerate states can therefore behave
differently from a reference that evaluates unused star-state divisions first.
The default real kind is IEEE double precision, preserving the flattening eta
`1e-99`; arithmetic reordering may change rounding relative to the existing code.

This directory is standalone and has not been wired into Flash-X's setup/build
system or validated against the full Spark suite.

## Flash-X host compatibility wrapper

`hy_rk_getFaceFlux_wrapper.F90-mc` exports the original external
`hy_rk_getFaceFlux(stage,starState,flat3d,flx,fly,flz,limits,deltas,
scr_rope,scr_flux,scr_uPlus,scr_uMinus,loGC)` signature. It expands the same
computational helpers against raw Flash-X arrays using `hydro_layout_flashx.ini`.
It allocates no mesh arrays and makes no scratch copies. Small numerical vectors
remain automatic local scratch. The standalone `hy_getFaceFlux` likewise never
allocates its workspace; its allocation macros are called separately by callers.

The wrapper uses `limits(:,:,stage)` and the supplied `flat3d`; it does not
recompute flattening. `scr_rope` receives the original primitive state and rho*e
on the same NSTENCIL-expanded region as the legacy routine. `scr_uPlus/Minus`
receive reconstructed, density/pressure-floored values for every directional
reconstruction pass. `scr_flux` receives viscosity-corrected fluxes each pass.
As in the original routine, scratch contains the last direction's values in
regions overwritten by that direction; directional outputs are retained in
`flx/fly/flz`. Other entries have the original OUT-argument undefined semantics.
`deltas` is retained but does not enter this flux arithmetic.

The legacy signature does not pass configuration values. The adapter therefore
imports only the necessary configuration from `Hydro_data`; the computational
helpers continue to receive it explicitly. It imports immutable identifiers from
`newimpl_hydro`. Its temporaries use `kind(1.0)`, matching Flash-X default-real
arrays, including builds with default-real promotion. The standalone storage
backends retain their own double-precision kind.

From this directory, expand both files:

```sh
python3 macro_expand.py --definitions hydro_layout.ini hydro_helpers.ini \
  --input hy_getFaceFlux.F90-mc --output build/hy_getFaceFlux.F90
python3 macro_expand.py --definitions hydro_layout.ini hydro_helpers.ini hydro_layout_flashx.ini \
  --input hy_rk_getFaceFlux_wrapper.F90-mc --output build/hy_rk_getFaceFlux.F90
```

For a Flash-X test build, compile `newimpl_hydro` from the first generated file
before the wrapper; compile the wrapper with the build's generated Simulation.h,
constants.h, Spark.h, Hydro_data module, and normal compiler flags. Link the
wrapper **in place of** the old `hy_rk_getFaceFlux` implementation, avoiding two
objects exporting the same symbol. The existing interface declaration can remain.
The wrapper is a serial host adapter; accelerator variants and scheduling are not
implemented. This is a manual test-build integration, not a change to setup rules.

Unsupported MHD/species/mass-scalar configurations are rejected during
preprocessing. WENO5 requires at least three cells beyond the requested cell
region; the wrapper rejects NSTENCIL<3 and checks the declared input bounds.
Callers must provide sufficiently large output and scratch arrays as required by
the existing interface. Riemann errors call `Driver_abort` rather than returning
partial fluxes through a signature without a status argument.

The test runner compiles the wrapper with the repository's actual Spark.h,
mock generated headers, and mock Hydro_data/Driver_abort. It verifies the original
signature, stage-indexed bounds, supplied flattening, noncontiguous input IDs,
prepared and reconstructed scratch values, final-direction flux scratch, constant
physical fluxes, and nonuniform hybrid HLLE fluxes in 1D/2D/3D for both default-real
and promoted-double builds under Fortran 2003. The full Spark suite remains to be
run in an actual Flash-X build.

## Solution update

`hy_updateSolution.F90-mc` adds the `newimpl_update` module with a single top-level
routine. `update_helpers.ini` provides conservative loading, pre-update gravity
sources, directional flux divergence, arbitrary supplied coefficient combination,
density/energy recovery, optional potential extrapolation, and selective output.
It shares the existing expander and layout backends. No mesh-sized scratch is
allocated or required. Local five-component vectors are hoisted by the expander.

The update is Cartesian hydrodynamics without species or mass scalars. Its
`stateMap(6)` denotes density, x/y/z velocity, stored specific total energy
(ENER), and specific internal energy (EINT). Conservative input energy is rho*ENER;
EINT is never used to construct it. `fluxMap(5)` denotes mass, three momenta, and
energy. `gravityMap(3)` identifies acceleration components; `potentialMap(2)`
identifies current and previous potential fields. Inactive optional maps may be
zero. `referenceState` is immutable; `currentState` is updated in place and must
not alias the reference. Gravity and face fluxes are supplied immutable inputs.
All configuration is explicit; only immutable types/status identifiers and
intrinsic IEEE procedures are imported.

The caller provides valid storage for the requested cells, both input states,
and incident faces through `cellHi+unit(direction)`. Bounds are inclusive and
there is no need for neighboring state cells. Pressure, temperature, EOS
coefficients, other diagnostics, previous potential, and cells outside the
requested region remain unchanged. Gravity uses all three components of the
pre-update stage momentum even in reduced dimensions. Floors use the same
floored density throughout recovery; internal-energy flooring does not reset
stored total energy. Arbitrary coefficients are accepted, including a+b!=1 and
zero c*dt. Potential extrapolation requires nonzero dtOld only when enabled.
Status and failedCell identify errors; on failure previously updated cells remain
updated, so the caller must reject that partial computation.

Expand with:

```sh
python3 macro_expand.py --definitions hydro_layout.ini update_helpers.ini \
  --input hy_updateSolution.F90-mc --output build/hy_updateSolution.F90
python3 macro_expand.py --definitions hydro_layout.ini update_helpers.ini hydro_layout_flashx.ini \
  --input hy_rk_updateSoln_wrapper.F90-mc --output build/hy_rk_updateSoln.F90
```

Compile `newimpl_hydro` before `newimpl_update`. Append a chosen alternate storage
INI after the base definitions to select a different layout.

`hy_rk_updateSoln_wrapper.F90-mc` preserves the original external Flash-X signature,
including unused geometry arrays. It selects limits and coefficients by stage,
imports the legacy configuration only in the adapter, and rejects non-Cartesian
geometry. Gravity and potential follow the existing GRAVITY/GPOT_VAR preprocessing
flags. MHD/species/scalar builds are rejected. Compile and link it in place of the
old update routine, avoiding duplicate symbols. No EOS, gravity solve, abundance
postprocessing, or flux correction is added; those surrounding call paths remain
in `Hydro_advance`.

The default wrapper saved-state declaration uses loGC lower bounds. For a
non-telescoping build, append `update_flashx_nontelescoping.ini` **last** to the
wrapper definition list; this changes tmpState lower bounds to lo, matching
Spark/NonTelescoping/NI.ini. Select this according to the existing build, not from
the current stage number. The ABI matches the existing assumed-shape interface.

The existing test runner now also verifies update arithmetic across all three
layouts, noncontiguous maps, both blocks, 1D/2D/3D flux contributions, gravity work,
potential extrapolation, unrelated-field preservation, floor order, zero-stage
scale, invalid inputs, and telescoping interior-face conservation. Update wrapper
checks use the actual Spark.h with mocked build configuration in 1D/2D/3D,
default-real telescoping and promoted-double non-telescoping variants. All checks
compile with Fortran 2003, bounds checking, and floating-point traps. Full Flash-X
setup integration and the Spark suite remain separate validation steps.

## Shock detection and timestep calculation

`hy_shockDetect.F90-mc` and `hy_computeDt.F90-mc` provide `newimpl_shock` and
`newimpl_dt`. Each module has one top-level subroutine; computational helpers
are in `diagnostic_helpers.ini` and share the existing expander. Inputs and
configuration are explicit, with no mutable module globals in the cores.
The same three storage backends are supported. `line_ref` covers the supplied
one-dimensional coordinates, cell widths, and grid velocities.

The arithmetic follows `docs/designDocs/shock_detect_equations.pdf` and
`hydro_compute_dt_equations.pdf`, cross-checked against the active source files.
The shock PDF describes `Hydro_funcs.F90/shockDetect`, whose detection region is
explicit. The active `hy_rk_shockDetect` instead trims one cell from each active
axis of blkLimitsGC. The new core therefore accepts a sound/reset region
`cellLo:cellHi` and an independent detection region `detectLo:detectHi`.
The latter must allow one-cell neighbors within the former on each active axis.
To reproduce the PDF's full-array reset, pass the full allocated domain as the
sound/reset region. To reproduce the active entry point, use its guard limits
and the trimmed detection region, as the wrapper does.

Both cores use `map(6)` = density, ux, uy, uz, pressure, gammaC. Shock detection
writes the separately supplied shock-variable identifier and caller-owned sound
scratch through `flat_ref`; no mesh scratch allocation is performed. Nonpositive
shockVariable is a successful no-op. The sound denominator is max(rho,densityGuard),
with beta=0.5 and delta=0.1 supplied by the compatibility wrapper. Neighborhood
minima include diagonal cells (3/9/27 members). Velocity/pressure differences are
undivided, independent of spacing and geometry; both threshold comparisons are
inclusive. All marker values in the sound/reset region are cleared before
classification. Scratch and marker entries outside that region are untouched in
the core. Errors can leave a partially computed scratch/marker region.

The timestep core uses **no density floor**, includes velocities relative to the
supplied grid velocities, and takes the maximum directional rate rather than a
sum. Cartesian widths come from the lower interior index, matching the source
rather than introducing a per-cell Cartesian spacing formula. Cylindrical 3D
uses r*dz for the angular z width; spherical active angular directions use r*dy
and r*sin(theta)*dz. Unused angular directions are not evaluated. The direction
constants are DT_CARTESIAN, DT_CYLINDRICAL, DT_SPHERICAL. All active widths must be
positive; supplied arrays and maps must cover their requested indices.

`dtCheck`, `dtMinLoc(5)`, and `maxSignalSpeed` are caller-owned accumulators.
For positive CFL, the first positive maximum rate in i-fastest/j/k traversal
wins. A candidate replaces dtCheck/location only when strictly smaller, while
maxSignalSpeed accumulates across calls regardless of whether dtCheck changes.
Location is [i,j,k,level,rank], with level and rank explicit. Zero rate yields
huge(kind-real) and no limiting location. dt and signal-speed outputs are committed
only after a successful scan. Positive finite CFL and incoming dtCheck, finite
nonnegative incoming maximum speed, valid geometry, and valid acoustic states are
required. Invalid input returns status and failedCell instead of unchecked
square roots/divisions. This targets admissible-input arithmetic, not the source's
unguarded invalid-input behavior. CFL-weighted rate comparisons are simplified
to rate comparisons under the positive-CFL requirement.

Compatibility adapters are `hy_rk_shockDetect_wrapper.F90-mc` and
`Hydro_computeDt_wrapper.F90-mc`. They preserve the existing argument lists and
import legacy mutable configuration only at the adapter boundary. The shock
adapter uses caller-supplied Vc, zeros its whole array as the active routine does,
and preserves markers outside blkLimitsGC. Without SHOK_VAR it returns immediately.
The timestep adapter preserves explicit-shape coordinate-array arguments and the
pointer U argument; it retains the hydro-enabled early return, first-call flag,
hy_lChyp accumulation, level/rank, and unchanged optional extraInfo. Other legacy
geometry values use Cartesian factors, exactly as the source's default branch.
SPARK_GLM is rejected by the timestep adapter; magnetic waves remain outside this
hydrodynamic implementation. These are serial host adapters, not offload variants.

From this directory:

```sh
python3 macro_expand.py --definitions hydro_layout.ini diagnostic_helpers.ini \
  --input hy_shockDetect.F90-mc --output build/hy_shockDetect.F90
python3 macro_expand.py --definitions hydro_layout.ini diagnostic_helpers.ini \
  --input hy_computeDt.F90-mc --output build/hy_computeDt.F90
python3 macro_expand.py --definitions hydro_layout.ini diagnostic_helpers.ini hydro_layout_flashx.ini \
  --input hy_rk_shockDetect_wrapper.F90-mc --output build/hy_rk_shockDetect.F90
python3 macro_expand.py --definitions hydro_layout.ini diagnostic_helpers.ini hydro_layout_flashx.ini \
  --input Hydro_computeDt_wrapper.F90-mc --output build/Hydro_computeDt.F90
```

Compile newimpl_hydro first and newimpl_dt before the timestep adapter. Link each
adapter instead of the corresponding old symbol, with the build's normal headers
and modules. The tests cover all three layouts, reduced dimensions, pressure and
compression switches, threshold equality, diagonal minima, guarded sound speed,
Cartesian/cylindrical/spherical rates, grid motion, tie ordering, incoming
constraints, zero rate, and invalid states. Adapter checks run 1D/2D/3D in both
default and promoted-double precision with mocked Flash-X configuration. All use
Fortran 2003, bounds checks, and floating-point traps. Full Spark-suite validation
and build-system integration have not been performed.

## Combined standalone preparation and advance

See [STANDALONE.md](STANDALONE.md) for the Grid-provider callback contract,
hy_prepareAdvance driver, periodic Sod example, build commands, and end-to-end
tests. The driver uses caller-owned scratch, explicit SSPRK tables, periodic
exchange between stages, EOS recovery, and area/time-integrated flux buffers
without Flash-X dependencies.
