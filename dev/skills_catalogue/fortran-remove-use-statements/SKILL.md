---
name: fortran-remove-use-statements
description: Convert Fortran subroutines to remove use statements by converting module variables into subroutine arguments. Handles all related files including calling sites, interface declarations, and stub implementations.
---

# Fortran Remove Use Statements

This skill guides you through converting Fortran subroutines to eliminate dependencies on module use statements by converting imported variables into subroutine arguments.

## When to Use This Skill

Use this skill when:
- A Fortran subroutine has a `use module_name, ONLY: ...` statement you want to remove
- The subroutine accesses module-level variables that should be passed as arguments
- You need to maintain the same functionality while changing the subroutine interface
- The conversion involves multiple files that need to be updated consistently

## Conversion Process

Follow these steps systematically to convert a subroutine:

### 1. Analyze the Target Subroutine

Read the target subroutine file and identify:
- The `use module_name, ONLY: ...` statement line(s)
- All variables imported from the module
- How each imported variable is used (read-only read/write, or written)
- Current subroutine signature and documentation

**Example:**
```fortran
subroutine gr_mpoleCen1Dspherical (idensvar)
  use gr_mpoleData, ONLY: gr_mpoleDrInnerZone, gr_mpoleTotalMass, ...
```

### 2. Read the Module Definition

Read the module file (e.g., `gr_mpoleData.F90`) to understand each imported variable:
- Data type (real, integer, logical, etc.)
- Intent (should be determined from usage context)
- Save attribute (usually global variables)
- Descriptions and meanings

### 3. Determine Variable Count and Intents

For each imported variable:
- **Input variables** (`intent(in)`): Variables that are only read within the subroutine
- **Output variables** (`intent(out)`): Variables that are written/calculated within the subroutine
- **Mixed usage**: Variables that are both read and written should be `intent(inout)`

Track the conversion systematically:
```
Variable Name           | Type         | Intent     | Usage Context
-----------------------|--------------|------------|----------------------------------
gr_mpoleDrInnerZone    | real         | out        | Calculated and assigned
gr_mpoleTotalMass      | real         | out        | Computed from local sums
gr_mpoleDomainRmin     | real         | in         | Only used for comparison
```

### 4. Update the Main Subroutine

Modify the target subroutine file:
- **Remove** the entire `use module_name, ONLY: ...` statement
- **Add** new parameters to the subroutine signature in the same order as the use statement
- **Add** parameter declarations with proper intent specifications
- **Update** the documentation (SYNOPSIS, ARGUMENTS sections) to describe each new parameter

**Example change:**
```fortran
subroutine gr_mpoleCen1Dspherical (idensvar, gr_mpoleDrInnerZone, gr_mpoleDrInnerZoneInv, ...)
  integer, intent(in)  :: idensvar
  real,    intent(out) :: gr_mpoleDrInnerZone    ! New output parameter
  real,    intent(out) :: gr_mpoleDrInnerZoneInv ! New output parameter
  ...
```

### 5. Find All Calling Sites

Search comprehensively for files that call this subroutine:
```bash
grep -r "subroutine_name" --include="*.F90" source/
```

Identify:
- Files containing subroutine calls (`call subroutine_name (...)`)
- Interface declaration files
- Stub implementation files
- Any other references

### 6. Update Calling Sites

For each file that calls the subroutine:
- **Add** the new variables to the `use module_name` statement if they come from that module
- **Update** the subroutine call to include all new parameters in the same order
- **Ensure** variable names match those defined in the module

**Example update:**
```fortran
use gr_mpoleData, ONLY: gr_mpoleGeometry, gr_mpoleDrInnerZone, ...
call gr_mpoleCen1Dspherical (idensvar, gr_mpoleDrInnerZone, gr_mpoleDrInnerZoneInv, ...)
```

### 7. Update Interface Declarations

Find and update interface declaration files:
- Locate the interface block for the subroutine
- **Add** all new parameters with their types and intents
- **Ensure** the interface matches the new subroutine signature exactly

**Example update:**
```fortran
interface
   subroutine gr_mpoleCen1Dspherical (idensvar, gr_mpoleDrInnerZone, ...)
     integer, intent (in) :: idensvar
     real,    intent (out) :: gr_mpoleDrInnerZone
     ...
   end subroutine gr_mpoleCen1Dspherical
end interface
```

### 8. Update Stub Implementations

If stub implementations exist (often in localAPI directories):
- **Update** the subroutine signature to match the main implementation
- **Add** parameter declarations with proper intents
- **Update** documentation to match the new signature
- **Keep** stub functionality minimal (typically just `return`)

### 9. Verify Consistency

After completing all changes:
- **Search** again for all references to ensure nothing was missed
- **Check** that parameter counts and types match across all files
- **Verify** documentation accuracy in all updated files
- **Confirm** that the same variable names are used throughout

### 10. Testing (when possible)

If the codebase has build/test capabilities:
- **Attempt** to build the modified files
- **Run** relevant tests to verify functionality is preserved
- **Report** any compilation errors or test failures for resolution

## Common Patterns and Gotchas

### Pattern 1: Multiple Use Statements
Some subroutines may have multiple `use` statements from different modules. Convert them systematically one at a time.

### Pattern 2: Different Variable Types
Handle different data types appropriately:
- **Logical variables**: Use `logical` type, typically `intent(in)` for configuration flags
- **Real/Integer**: Use appropriate type and precision matches the module definition
- **Arrays**: Preserve array dimensions and characteristics in parameter declarations

### Pattern 3: Upper/Lower Case Name Variations
Fortran is case-insensitive, but maintain consistency:
- Use the same variable names as defined in the module
- Preserve the same capitalization throughout the conversion

### Pattern 4: Documentation Inconsistencies
Always update documentation to reflect the new interface:
- Update SYNOPSIS to show the new signature
- Add detailed ARGUMENTS documentation for each new parameter
- Update any examples or usage notes

## Error Handling

If you encounter issues:
- **Missing variables**: Check the module file again for any overlooked imports
- **Type mismatches**: Verify parameter types match the module definitions exactly
- **Compilation errors**: Check for missing parameters, wrong intent declarations, or name typos
- **Interface mismatches**: Ensure all files use exactly the same signature

## Completion Criteria

The conversion is complete when:
1. The `use module_name` statement is removed from the target subroutine
2. All imported variables are declared as parameters with correct intents
3. All calling sites are updated with the new parameter list
4. Interface declarations match the new signature
5. Stub implementations are updated
6. Documentation accurately describes all parameters
7. No compilation errors exist (when verifiable)

## Notes

- Always use the same variable names as in the original module to avoid search/replace complexity
- Maintain the original order of parameters as they appeared in the use statement
- Focus on one subroutine conversion at a time to avoid confusion
- This process eliminates module dependencies while preserving functionality
- The converted subroutines are more modular and easier to test independently