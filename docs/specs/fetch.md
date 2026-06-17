# Fetch Subcommand Spec

## Goal

Define the `prosig fetch` subcommand for retrieving external files needed by ProSig workflows.

## Scope

Use this spec to describe:

- What datasets or resources `fetch` should retrieve.
- Where fetched files should be written.
- Which files are required versus optional.
- How existing files should be handled.
- How downloads should be verified.
- What user-facing CLI options are required.

## Working Directory Convention

Fetched files should be treated as reproducible working artifacts and should default to the local ignored `work/` directory unless this spec says otherwise.

## Proposed Command Shape

```bash
prosig fetch [OPTIONS]
```

Potential options to specify:

- `--out-dir PATH`: directory where fetched files are written.
- `--force`: overwrite existing files.
- `--dry-run`: show what would be fetched without writing files.

## Tests

Implementation must add or update focused tests for:

- CLI help and option parsing.
- Expected output paths.
- Existing-file behavior.
- Download or fetch behavior, preferably with network access mocked.

## Open Questions

- Which external resources should be fetched first?
- Should `fetch` download from the network directly, copy from local paths, or support both?
- Should checksum or size validation be required?
- Should output be human-readable text, JSON, or both?
