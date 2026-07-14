# SPEC.md — modularize a Flash-X unit (remove `use`, own allocations)

The correctness specification for this loop: what "done" means and where the
verification boundary sits. It is the *oracle* — a change is done not when it
compiles or the agent says so, but when it meets the criteria below and a human
accepts it. The concrete target (unit path, data module, file list, acceptance
command) lives in PLAN.md; this file names no paths.

## Goal

Reduce a unit's dependence on module-level global state by:

1. Converting `use <module>, ONLY: ...` dependencies in a subroutine into
   explicit dummy arguments (skill: `fortran-remove-use-statements`), and
2. Bringing dynamic allocation of module-declared `allocatable, save` arrays into
   the data module itself (skill: `move-allocatable-arrays-into-data-module`).

## Invariants (what the transformation must preserve)

- **Behavior is unchanged.** This is a refactor, not a numerical change; the
  computed results must be identical.
- **Argument order matches the `use ... ONLY` order.** New dummy arguments are
  appended in the order the variables appeared in the removed `use` clause.
- **Intents are correct.** `intent(in)` read-only, `intent(out)` written-only,
  `intent(inout)` both. Allocatable dummies are `intent(inout)`.
- **All edges move together.** For every converted subroutine, the caller sites,
  the interface declaration, and any `localAPI` stub must match the new signature.
  A half-converted subroutine is a defect, not progress.
- **Same names throughout.** Keep the original module variable names as the dummy
  argument names.
- **Makefiles stay consistent.** If a standalone file is deleted, its `.o` entry
  is removed from the relevant Makefile.

## Boundaries

- Edit only files under the target unit (named in PLAN.md), plus the caller,
  interface, stub, and Makefile edges a single conversion requires.
- Never edit `dev/`, `sites/`, `lib/`, or any submodule.
- One subroutine per iteration; do not batch unrelated files.

## Verification boundary

The loop verifies *mechanically* each iteration; final acceptance is the human's.

### Automated checks the loop must run (grep-level consistency)

For each subroutine touched this iteration, a search (the concrete commands are
in PLAN.md) must show that:

- the removed `use` line is gone from the target file;
- no calling site still uses the old signature (argument count matches);
- the interface block signature matches the implementation signature;
- no stale reference to a deleted file/subroutine remains in `*.F90` or Makefiles.

Report the exact command output; never paraphrase or fabricate a result. An
iteration that cannot show these searches came back clean is INCOMPLETE.

### Human verification boundary (acceptance)

The loop never marks work accepted. A maintainer accepts a converted unit only
after running the build + regression named in PLAN.md and confirming the result
matches the reference within tolerance. Until then, converted subroutines are
*translated*, not *verified* (paper principle P1).
