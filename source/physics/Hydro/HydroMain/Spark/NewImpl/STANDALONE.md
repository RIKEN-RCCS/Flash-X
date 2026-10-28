# Periodic uniform Cartesian driver

hy_prepareAdvance.F90-mc combines the preparation and advance roles of
Hydro_prepBlock and Hydro_advance. It uses the existing macro expander,
face-flux/update/shock kernels, and storage backends. It has no Flash-X header,
module, or mutable global dependency. The target is a static uniform Cartesian
block mesh, hydrodynamics without gravity, potential, species, or magnetic fields.

## Build and run

From this NewImpl directory:

~~~sh
python3 build_standalone.py
build/standalone/periodic_shock_tube 256 4 0.05 build/standalone/sod.csv 2
python3 test_standalone.py
python3 test_newimpl.py
~~~

Arguments are cell count, block count, final time, output CSV path, and SSPRK
order (2 or 3). Defaults are 256 cells, four blocks, final time 0.05, output
periodic_shock_tube.csv, and SSPRK2. Cell count must be divisible by block count.
The example uses an ideal gas with gamma=1.4, left rho=1 and p=1, right rho=0.125
and p=0.1, and zero initial velocity on [0,1).

Periodicity produces two discontinuities: the usual Sod interface at x=0.5 and
its reversed partner at the domain boundary. This is not an isolated shock tube
with outflow boundaries. The early-time plateau near x=0.5 can be compared with
the ordinary Sod solution before the two wave systems interact.

Use --layout spatial-first or --layout records with the build script to select
another backend; use separate --output directories when keeping several builds.
Compilation uses Fortran 2003, bounds checks, and floating-point traps. The CSV
contains x, density, x velocity, pressure, specific internal energy, and specific
total energy. The executable reports changes in conserved quantities and the
integrated periodic flux imbalance.

## Grid provider contract

hydro_grid_contract.F90 defines the metadata and explicit callback interfaces.
Supply a context extending hydro_grid_context_t and one hydro_block_t per local
block with a unique storage block identifier, interior/guard bounds, spacing,
face areas, volume, level, and rank. Supply these callbacks:

- fillGuards: exchange all variables, including shock markers, across neighbors
  and periodic boundaries. Include corner guards in multiple dimensions.
- synchronizeFlux: reconcile shared faces to one consistently oriented flux.
  Complete communication before returning.
- reduceLimits: global minimum timestep, maximum signal speed, and corresponding
  limiting location. The local timestep already includes the requested time cap.
- evaluateEos: recover pressure and gammaC from density and specific internal
  energy over each interior. Preserve stored total energy.

The driver needs no coordinates. The provider supplies geometry; this uniform
Cartesian version checks A_d*h_d/V=1. No moving mesh, curvilinear momentum sources,
AMR interpolation, refluxing, or level subcycling is implemented.

periodic_grid_example.F90-mc is a replaceable serial 1D provider. It maps guards
to periodic interior donors and makes each positive block-boundary flux equal
the next block's negative boundary flux. It supplies geometry, identity global
reductions, and ideal-gas EOS recovery. The driver accepts dimensions 1–3; a
production provider must implement the corresponding multidimensional exchange,
EOS traversal, and face reconciliation.

## Ownership and numerical inputs

The caller allocates separate solution, reference, stage, acoustic scratch, face
flux, integrated flux, and hydro_workspace objects once and reuses them. The
advance routine performs no mesh allocation. Writable objects must be disjoint.

- inputMap(7): density, ux, uy, uz, pressure, gammaC, specific internal energy.
- updateMap(6): density, ux, uy, uz, specific total energy, specific internal energy.
- fluxMap(5): mass, three momenta, total energy.
- shockVariable: marker identifier, or zero to disable detection.

Total and internal energy are distinct fields. Allocate three guards in every
active axis and include positive boundary faces in flux storage. The driver
checks metadata but cannot verify arbitrary backend allocation sizes or map
upper bounds. At least one local block is required; empty MPI ranks require
separate orchestration.

Provide CFL, timestep cap, floors, coefficients(3,nStages), and stage flux weights.
The demo uses Spark's SSPRK2/SSPRK3 tables. Other tables are accepted; their
stability and flux-weight consistency are the caller's responsibility.

## Sequence of one call

1. Evaluate the input EOS and exchange periodic guards.
2. Detect shocks, then exchange the updated markers.
3. Estimate the static Cartesian CFL limit and reduce timestep/signal speed.
4. Save the complete prepared state into reference and stage scratch.
5. For each stage: exchange guards; compute all blocks' fluxes; reconcile shared
   faces; accumulate integrated fluxes; update all interiors; evaluate the EOS.
6. Exchange final guards and copy the stage state back to the solution.

Unlike telescoping mode, guards are exchanged between stages and only interiors
are updated. Every block's fluxes are complete before any stage state changes.
Shock markers are computed once per timestep and retained through the stages.

integratedFlux stores sum(weight_s*dt*faceArea*flux_s), reset on stage one.
These are area- and time-integrated transfers, distinct from instantaneous flux
densities and suitable for conservation checks or a future flux-register adapter.

Floors can change the conserved budget; EOS recovery does not reset total energy
to internal plus kinetic energy. No positivity fallback or failed-step retry is
supplied. On failure, callbacks may have changed input diagnostics and scratch
may contain partial work. Check status and reject the call; the demo stops.
The main solution advance is committed after a successful final guard exchange.

## Validation

test_standalone.py runs SSPRK2 and SSPRK3 on one and four blocks for all three
layouts. It checks positive finite output, mass/momentum/energy conservation,
EOS consistency, the Sod left-star plateau, and agreement across layouts and
block decompositions. The executable checks integrated periodic flux balance.
Existing kernel/wrapper tests remain in test_newimpl.py. These are serial
standalone checks, not MPI or full Flash-X integration tests.
