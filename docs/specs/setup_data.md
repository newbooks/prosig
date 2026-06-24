# Setup Data Subcommand Spec

## Goal

Define the `prosig setup-data` subcommand for retrieving external files needed by ProSig workflows.

## Scope

Use this spec to describe:

- What datasets or resources `setup-data` should retrieve.
- Where fetched files should be written.
- Which files are required versus optional.
- How existing files should be handled.
- How downloads should be verified.
- What user-facing CLI options are required.

## Destination Convention

Fetched files are written to the current working directory by default.

Destination values in the download source table are filenames relative to the current working directory.

## Download Helper Policy

Shared download behavior belongs in `prosig.io.download`.

The shared `download_file()` helper must not decide whether overwriting is allowed. If `download_file()` is called, it always performs a fresh download and atomically replaces any existing destination file after the new download completes.

Force or skip behavior belongs in the calling command or workflow. `prosig setup-data` should check whether a destination already exists and decide whether to call `download_file()` based on user options such as `--force`.

Default existing-file behavior:

- If the destination exists and `--force` is not set, skip the file.
- Log that the file was skipped because it already exists.
- Include a message telling the user to use `--force` to overwrite existing files.
- If `--force` is set, call `download_file()` and atomically overwrite the destination after the new download completes.

Checksum validation is out of scope for the first implementation.

## Progress Logging

Downloads should report progress through the shared ProSig logger.

Log format:

```text
[LEVEL]: message
```

Example:

```text
[INFO]: Downloading GO Graph: https://current.geneontology.org/ontology/go-basic.obo -> go-basic.obo
```

The root `prosig` command should allow users to set the logger level.

Default progress log behavior:

- Log progress every 30 seconds by default.
- Include downloaded bytes, total bytes when known, and percentage when total bytes are known.
- Use human-readable KB units in progress messages.

Example:

```text
Downloaded 3,222 KB / 123,943 KB (2.5%)
```

If the total size is unknown, log only the downloaded amount:

```text
Downloaded 3,222 KB
```

## Threaded Download Policy

Use multi-threaded downloads when supported by the source server.

Default threaded download behavior:

- Use 16 download threads by default.
- Add a `setup-data` option to alter the thread count.
- Log the requested thread count when starting a download.
- Log whether the completed download used multi-threaded or single-threaded mode.
- Apply multi-threaded downloading to files larger than 50 MB when file size can be detected before downloading.
- If file size cannot be detected before downloading, use multi-threaded downloading by default.
- Fall back to single-threaded downloading if the server does not support ranged requests or if threaded downloading cannot be used safely.

## Proposed Command Shape

```bash
prosig setup-data [OPTIONS]
```

Potential options to specify:

- `--force`, `-f`: overwrite existing files.
- `--dry-run`: show what would be fetched without writing files.
- `--threads INTEGER`: number of download threads to use when threaded downloading is supported. Default: `16`.

## Tests

Implementation must add or update focused tests for:

- CLI help and option parsing.
- Expected output paths.
- Existing-file behavior.
- Download helper behavior with network access mocked.
- Atomic overwrite behavior when a destination already exists.
- Failed download behavior that does not corrupt an existing destination file.
- Skip behavior when a destination exists and `--force` is not set.
- Progress logging behavior.
- Thread count option parsing and propagation.
- Threaded versus single-threaded download selection based on detected file size and range support.

## Download Sources

Do not unpack compressed files when downloading. We will read and process these files in stream mode.

The `setup-data` command will only download from the internet.

| Description | Source URL | Destination |
| --- | --- | --- |
| GO Graph | https://current.geneontology.org/ontology/go-basic.obo | go-basic.obo |
| Swiss-Prot GO | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz | uniprot_sprot.dat.gz |
| PROSITE | https://ftp.expasy.org/databases/prosite/prosite.dat | prosite.dat |
