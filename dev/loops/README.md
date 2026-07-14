# dev/loops/ — transformations that apply skills

A **loop** applies one or more `dev/skills/` capabilities across a checklist of
targets until a correctness spec is met. It is the minimal implementation of a
workflow: a bounded, deterministic executor→reviewer cycle whose state lives in
files on disk.

The reference engine here is **CodeScribe** (`code-scribe loop`). It runs a
bounded number of loops (default 5); each loop is an *execution* agent that
advances the plan as far as it can in one session, followed by a *review* agent
that checks what actually happened and records what remains. Cross-loop state is
carried by the harness; durable artifacts land in `.codescribe/loop/`
(`run.toml`, `state.toml`, `execution.toml`, `review_output.toml`) for inspection
and crash-resume. A loop stops early when the execution agent reports
`STATUS: COMPLETE` or the reviewer finds no pending items and no blocker.

## Layout of a loop directory

    loops/<name>/
      loop.toml   the CodeScribe task file: a [tools].bash allowlist (a hard
                  guardrail) + a minimal [[chat.user]] seed that only points the
                  agent at the inputs below
      SPEC.md     the plan-agnostic oracle: what "done" means, the invariants,
                  and the verification boundary — names no paths
      PLAN.md     the concrete target (unit path, data module, acceptance
                  command, consistency-check commands) + the checklist the loop
                  updates ([ ] open · [~] drafted, awaiting review · [x] accepted)
      run.sh      the exact command, pinned to the repo root as --workdir

Division of labour: `PLAN.md` owns everything concrete (paths, commands, the
target list) and is the only file that changes during a run; `SPEC.md` is the
reusable correctness contract; the `dev/skills/` procedures are generic (how to
do it, gotchas, error handling). Acceptance (`[x]`) is always a human's.

## `loop.toml` in brief

CodeScribe reads it as a chat template:

    [tools]
    bash = ["grep", "find", "rg"]   # bare command names only — the hard guardrail

    [[chat.user]]
    content = '''
    Read PLAN.md and SPEC.md; load the named skill; advance one plan item...
    '''

The `bash` list bounds ONLY the execution phase; the review phase runs under its
own near-read-only allowlist. Multi-line content must use triple single quotes.
Keep a full `./setup`/`make` build OUT of the allowlist — that is the human's
verification boundary, not the loop's.

## Running a loop

    # from the repo root; model needs a backend prefix (anthropic-/openai-/...)
    dev/loops/modularize-unit/run.sh anthropic-claude-sonnet-4-5 5

Watch a run, edit `PLAN.md` between runs to re-steer, and review the `[~]` items
before promoting any to `[x]`.

## Adding a loop

1. `mkdir dev/loops/<name>` and write `SPEC.md` and `PLAN.md`.
2. Write `loop.toml`: a minimal `[tools].bash` allowlist + a `[[chat.user]]` seed
   pointing at those two files and the `dev/skills/` skill(s) it uses.
3. Copy `run.sh` and adjust the task-file path.

Validate the task file parses before running:

    python3 -c "from codescribe.lib import load_chat_template as f; \
      print(f('dev/loops/<name>/loop.toml', return_meta=True)[1])"

## Worked example

`modularize-unit/` drives `fortran-remove-use-statements` and
`move-allocatable-arrays-into-data-module` over
`source/Grid/GridSolvers/Multipole_new`, reducing that unit's dependence on
`gr_mpoleData` global state one subroutine at a time.
