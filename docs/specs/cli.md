# CLI Spec

## Goal

Provide a command and subcommand interface for ProSig workflows.

The CLI should be a thin user-facing layer over package modules. Command modules should parse arguments, validate user intent, and call implementation functions; they should not contain motif discovery, signature scanning, or prediction logic directly.

## CLI Framework

Use Typer for the command system.

Rationale:

- Typer has been used successfully in related prior work.
- It supports command and subcommand apps with `typer.Typer()`.
- Type hints can define CLI argument and option parsing.
- It keeps command definitions readable as the command tree grows.

## Package Layout

Recommended initial layout:

```text
src/prosig/
  cli/
    __init__.py
    app.py
    setup_data.py
    discover.py
    annotate.py
  discovery/
  prediction/
  signatures/
  io/
```

## Entry Point

Expose a `prosig` executable through `pyproject.toml`:

```toml
[project.scripts]
prosig = "prosig.cli.app:main"
```

Typer should be tracked as a runtime dependency:

```toml
dependencies = [
  "typer>=0.12",
]
```

## Initial Commands

Start with these top-level subcommands:

- `prosig setup-data`: download and cache external data for offline use.
- `prosig build-library`: build the minimized GO graph, adjustable function
  clusters, and motif library.
- `prosig discover`: discover sequence signatures from positive and background or negative sequence sets.
- `prosig annotate`: scan sequences against a signature library and predict
  protein function from signature hits.
- `prosig signatures`: inspect, validate, convert, or summarize signature libraries.

## Design Rules

- Keep CLI commands small and focused.
- Put reusable behavior in package modules outside `prosig.cli`.
- Use explicit names for input and output options.
- Prefer file paths and formats that can be used in reproducible scripts.
- Return non-zero exit codes for invalid inputs or failed workflows.
- Keep output formats stable once documented.

## Open Questions

- Which command should be implemented first after `setup-data`?
- What should be the first canonical signature library file format?
- Should command output default to human-readable text, structured JSON/TSV, or both?
