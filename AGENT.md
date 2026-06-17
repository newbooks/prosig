# Agent Instructions

This repository is ProSig, a project for protein sequence signature discovery and function prediction from discovered signatures.

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
- Keep implementation notes in `docs/specs/` and task lists in `docs/todos/`.
- Keep large datasets, generated artifacts, and experiment outputs out of git unless explicitly approved.
- Preserve scientific assumptions and algorithm choices in specs or decision records before encoding them deeply in code.
