#!/usr/bin/env bash
# Run the modularize-unit CodeScribe executor-reviewer loop.
#
# Usage:
#   dev/loops/modularize-unit/run.sh [model] [nloops]
#
# Must be run from the repo root (it pins --workdir to the repo root so the
# agent can reach source/ and dev/). The task file is read-only to the agent;
# PLAN.md is the editable running record. Durable run artifacts land in
# .codescribe/loop/ (run.toml, state.toml, execution.toml, review_output.toml)
# for inspection and crash-resume.
set -euo pipefail

# Model string needs a CodeScribe backend prefix: anthropic- / openai- / argo- /
# oaic- , or a local model path. Override via arg or the CS_MODEL env var.
MODEL="${1:-${CS_MODEL:-anthropic-claude-sonnet-4-5}}"
NLOOPS="${2:-5}"

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

exec code-scribe loop dev/loops/modularize-unit/loop.toml \
  --model "$MODEL" \
  --workdir "$REPO_ROOT" \
  --agent-loops "$NLOOPS" \
  --verbose --log
