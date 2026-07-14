# PLAN.md — modularize Multipole_new

The checklist the loop advances, and the concrete target it advances against.
Status: `[ ]` open · `[~]` done this run, awaiting review · `[x]` accepted by a
maintainer. The loop sets `[~]`; only a human sets `[x]`, after the acceptance
build + regression below passes. Correctness criteria are in SPEC.md.

## Target

- Unit:            `source/Grid/GridSolvers/Multipole_new`
- Data module:     `gr_mpoleData.F90`
- Interface decls: `gr_mpoleData.F90` — mirror every signature change here
- Acceptance build + regression (human-run):

      ./setup unitTest/Gravity/Poisson3 -auto -3d +newMpole && make

## Consistency checks (run each iteration; see SPEC.md § Verification)

    grep -rn "<subroutine_name>" --include=*.F90 source/Grid/GridSolvers/Multipole_new
    grep -rn "<deleted_file>.o"  --include=Makefile source/

## Phase 1 — Bring allocations into the data module
Skill: `dev/skills/move-allocatable-arrays-into-data-module/SKILL.md`

- [x] Radial-array allocations moved into `gr_mpoleData` (commit 730c9f6c)
- [x] Reconcile allocatables before removing `use` lines (commit f93c3468)
- [ ] `gr_mpoleDeallocateRadialArrays` — decide Case A vs B; if allocation-only,
      fold into the module and remove the standalone file + Makefile entry

## Phase 2 — Remove `use gr_mpoleData` from subroutines
Skill: `dev/skills/fortran-remove-use-statements/SKILL.md`

Start with the center-of-expansion family (the skill's own worked example), then
fan out. ~39 files in this unit still `use gr_mpoleData`.

- [ ] `gr_mpoleCen1Dspherical.F90`  — note: skill example; convert first
- [ ] `gr_mpoleCen2Dspherical.F90`
- [ ] `gr_mpoleCen2Dcylindrical.F90`
- [ ] `gr_mpoleCen3Dspherical.F90`
- [ ] `gr_mpoleCen3Dcartesian.F90`
- [ ] `gr_mpoleCenterOfExpansion.F90` — note: caller of the Cen* family; update
      its call sites in the same iteration it converts

## Pass log (newest first)
- (loop writes one line per pass here)
