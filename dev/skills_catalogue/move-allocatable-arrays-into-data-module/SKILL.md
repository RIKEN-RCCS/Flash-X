---
name: move-allocatable-arrays-into-data-module
description: Bring dynamic allocations in one subroutine into the data module 
---

# Move Allocatable Arrays Into a Fortran Data Module

This skill guides you through relocating allocation logic for module-declared allocatable arrays from standalone subroutines into the module itself. The subroutine's `use` statement is eliminated by passing the arrays and their dimension extents as arguments.

## When to Use This Skill

Use this skill when:
- A Fortran module declares `allocatable` arrays (with `save`) that are allocated elsewhere
- You want the allocation logic to live inside the module as a public procedure
- The project has a separate interface module and/or stub files

## Two Cases

### Case A: The subroutine does ONLY allocations

If the subroutine performs **only** `allocate` calls (with optional status checks and error handling), move the entire subroutine into the module. Delete the standalone file, any stubs, and interface declarations. Callers import the procedure directly from the module.

### Case B: The subroutine has allocations PLUS other operations

Extract just the allocation logic into a **new subroutine** inside the module. Replace the original allocation statements in the original subroutine with a call to the new module procedure. The original file remains but delegates to the module.

## General Steps

### Step 1: Identify the allocatable arrays and their extents

Read the data module and note:

- **Array declarations**: For each `allocatable, save` in the module, note its name and rank
- **Extent variables**: The scalar variables (usually integers) that determine the bounds used in `allocate(a(1:extent))` calls

### Step 2: Determine which case applies

Read the allocating subroutine:

```
grep -rn "subroutine.*Name\|call allocate" --include="*.F90" path/to/search
```

- If every statement is an `allocate(..., stat=status)`, `if (status > 0) call abort(...)`, or comments — it's **Case A**
- If there are calculations, loops, I/O, or other logic between allocations — it's **Case B**

### Step 3: Add the module procedure

In the data module, just before `end module`, add:

```fortran
  public :: allocating_subroutine_name

contains

subroutine allocating_subroutine_name (array1, array2, ..., extent1, extent2, ...)

  use error_handler_module, ONLY : abort_subroutine   ! if applicable

  implicit none

  ! Arguments: allocatable arrays as intent(inout)
  real,    allocatable, intent(inout) :: array1 (:)
  integer, allocatable, intent(inout) :: array2 (:)
  type_something, allocatable, intent(inout) :: array3 (:,:)

  ! Arguments: extent integers as intent(in)
  integer,              intent(in)    :: extent1
  integer,              intent(in)    :: extent2

  integer :: status

  ! ...all allocate calls ...
  allocate(array1(1:extent1), stat=status)
  if (status > 0) call abort_subroutine ('ERROR: ...')

  return
end subroutine allocating_subroutine_name
```

**Naming note**: The dummy argument names can be the same as the module variable names, or different. Using the same names keeps things recognizable.

### Step 4: Handle Case A — Remove the standalone file and update callers

1. **Delete** the standalone implementation file
2. **Delete** any stub files (often in a `localAPI`-style directory)
3. **Remove** entry from corresponding Makefiles
4. **Remove** the subroutine's interface block from the interface module file
5. **Update each caller**:
   - Remove the subroutine from the interface-module `use` statement
   - Add the subroutine and the needed module variables to the data-module `use` statement
   - Add the arguments to the `call`

### Step 5: Handle Case B — Consolidate and wrap

1. **Add** the new allocation-only subroutine to the module (Step 3)
2. **Modify the original subroutine** to call the module procedure:

```fortran
subroutine original_subroutine (original_args...)

  use data_module, ONLY : helper => new_allocation_subroutine

  implicit none

  ! ... existing code that does non-allocation work ...

  call helper (module_array1, module_array2, ..., extent1, extent2)

  ! ... existing code that does more non-allocation work ...

  return
end subroutine original_subroutine
```

3. If the original subroutine had `use data_module, ONLY: array1, array2, ...`, you may still need those for the call arguments
4. No Makefile changes are needed (original file still exists)
5. No stub or interface changes are needed unless the original subroutine signature changed

### Step 6: Update interface declarations (if applicable)

If the project uses a separate interface module (e.g., `*Interface.F90`):

- **Case A**: Remove the interface block entirely
- **Case B**: Update only if the original subroutine's signature changed

### Step 7: Update Makefiles (Case A only)

Find Makefiles that reference the deleted `.o` file:

```
grep -rn "deleted_filename.o" --include="Makefile" source/
```

Remove the entry.

### Step 8: Verify no stale references remain

Search for the old subroutine name:

```
grep -rn "subroutine_name" --include="*.F90" source/
grep -rn "subroutine_name" --include="Makefile" source/
```

All legitimate remaining references should be:
- Inside the data module itself (the procedure definition)
- In callers that now import from the data module and pass arguments

## Argument Intent Conventions

| Usage | Intent |
|---|---|
| Array is allocated inside the subroutine (initially unallocated) | `intent(inout)` or `intent(out)` — `inout` is safer for allocatables in practice |
| Extent variable is read-only | `intent(in)` |
| Non-allocatable scalar used for bounds | `intent(in)` |

## Common Pitfalls

- **Name collision**: When the original subroutine and the new module procedure have the same name, use `use module, ONLY: alias => original_name` in the wrapper
- **Allocatable intent**: Not all compilers support `intent(out)` on allocatable dummy arguments the same way; `intent(inout)` is the safer default
- **Missing extents**: If an extent is computed at runtime before the call, make sure it's available at the calling site (it typically is, since the original code already had it via `use`)

## Completion Checklist

- [ ] Module procedure added to the data module
- [ ] For Case A: standalone file deleted
- [ ] For Case A: stub file deleted
- [ ] For Case A: interface block removed
- [ ] For Case A: Makefile entry removed
- [ ] All callers updated to import from data module
- [ ] All `call` sites updated with arguments
- [ ] No stale references to the old interface
