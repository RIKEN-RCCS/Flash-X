# dev/ — agentic docs for Flash-X

A minimal, engine-independent way to embed AI-driven refactoring work in this
repository. The methodology follows this principle: keep the human-authored inputs (the
capability, the plan, the spec) as files in the repo, separate from whatever
agent runs them, so the same task is reproducible, auditable, and reusable.

    dev/
      skills/     reusable capabilities  — WHAT to do, once, well
      loops/      runnable transformations — apply skills across a work-list

## The two pieces

### `skills/` — reusable capabilities
A **skill** is a self-contained procedure for one kind of change, written
independently of any agent. Each lives in its own directory as `SKILL.md` with a
frontmatter `name`/`description` and a step-by-step process plus a completion
checklist. Because a skill names no agent, the same `SKILL.md` can be run by
CodeScribe, Claude Code, or a human. Current skills:

- `fortran-remove-use-statements` — turn `use <module>, ONLY: ...` dependencies
  into explicit subroutine arguments, updating callers, interfaces, and stubs.
- `move-allocatable-arrays-into-data-module` — relocate allocation of
  `allocatable, save` module arrays into the data module itself.

Add a skill by creating `skills/<kebab-name>/SKILL.md`. Keep it about the
*capability*, not any one target file.

### `loops/` — transformations that use skills
A **loop** applies one or more skills across a **plan** (a checklist of targets)
until a **spec** (the "done" definition) is met. A loop directory holds the
human-authored inputs; the loop engine reads them and drives the agents. See
`loops/README.md` for the file layout and how to run one. The worked example,
`loops/modularize-unit`, drives both skills over the `Multipole_new` grid solver.

## Vocabulary (from the paper)

| Term        | Here                                                              |
|-------------|-------------------------------------------------------------------|
| Skill       | `skills/<name>/SKILL.md` — reusable capability + checklist        |
| Plan        | a loop's `PLAN.md` — human-seeded work-list, the running record   |
| Spec        | a loop's `SPEC.md` — correctness definition + verification boundary|
| Loop        | `loops/<name>/` — a bounded executor→reviewer run over the plan   |
| Agent       | the engine that runs a loop (e.g. CodeScribe, Claude Code)        |
| Transformation | one run of a loop moving the code from one state toward the spec |

## Principles this layout encodes

- **Inputs are separate from the engine.** `SKILL.md` / `PLAN.md` / `SPEC.md`
  are plain files; swap the agent without rewriting the task.
- **Verification at the boundary.** A change is done when it meets `SPEC.md` and
  a human accepts it — not when it merely compiles or an agent claims success.
- **Humans own intent and acceptance; agents own mechanical breadth.** Loops mark
  work drafted (`[~]`); only a maintainer marks it accepted (`[x]`).
- **Reproducibility.** The plan, spec, and skill are the durable record of what a
  run was asked to do and how it would be checked.
