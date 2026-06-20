# Agent Instructions

This repository is ProSig, a project for protein sequence signature discovery and function prediction from discovered signatures.

## Product Command Plan

- Use `setup-data` for the one-time or occasional workflow that downloads and
  caches external data for offline use.
- Use `build-library` for building the minimized GO graph, Leiden-based
  adjustable function clusters, and customizable motif library. Treat clustering
  as part of this workflow unless a later decision record explicitly separates
  it into its own user-facing command.
- Use `discover` for motif discovery from grouped function clusters and
  background sequences.
- Use `annotate` for scanning sequence(s), reporting motif hits, and predicting
  sequence function from those motif hits. Motif hits should remain visible as
  the evidence or reasoning behind predictions.

Avoid introducing top-level `cluster`, `scan`, or `predict` commands without a
documented design decision, because the current plan folds clustering into
`build-library` and folds scanning plus prediction into `annotate`.

## Branch Safety

- Before editing files, check the current branch with `git branch --show-current`.
- Do not edit files while on the `main` branch.
- If work is needed and the repository is on `main`, create or switch to a task branch first.
- If branch creation or switching is blocked by permissions, stop and ask for approval before editing.

## Neighbor Repository Access

- Agents may read from `../pclass`.
- Treat `../pclass` as reference material only unless the user explicitly asks for changes there.
- Use `../pclass` to inspect the current motif handling methods, naming conventions, data models, tests, and implementation tradeoffs before porting or redesigning similar behavior in ProSig.

## Work Style

- Prefer small, reviewable changes with focused specs or TODO updates.
- When implementing a function, module, or command, add or update focused tests in the same change.
- If tests cannot be added for an implementation change, document the reason and the remaining risk in the final response.
- Keep implementation notes in `docs/specs/` and task lists in `docs/todos/`.
- Keep large datasets, generated artifacts, and experiment outputs out of git unless explicitly approved.
- Preserve scientific assumptions and algorithm choices in specs or decision records before encoding them deeply in code.
