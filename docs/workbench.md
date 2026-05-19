# DotMatch Workbench

DotMatch Workbench is an optional local desktop app for reviewing and running
AssaySpec workflows. It is designed for labs that want a guided interface for
fixed-window known-target assays while keeping the command-line DotMatch engine
as the source of truth.

All sequencing data stays on the user's machine. The app runs the installed
`dotmatch` executable locally, writes ordinary AssaySpec outputs, and reads the
same `assay_manifest.json`, `assay_manifest.summary.tsv`, `sample_qc.tsv`,
audit, and autopsy files that workflow systems can already consume.

## Status

The Workbench is not part of the Bioconda recipe. It is a separate desktop app
app that detects `dotmatch` from `PATH`, a user-configured executable path,
or the `DOTMATCH_WORKBENCH_DOTMATCH` environment variable.

The current Bioconda target remains the command-line package. No hosted uploads,
accounts, telemetry, cloud storage, or external workflow adoption claims are
required for Workbench use.

## Local Model

Workbench sessions begin with a workspace directory. FASTQ, target, barcode,
AssaySpec, and output paths are accepted as paths relative to that workspace.
The backend canonicalizes paths before use and rejects absolute paths, parent
traversal, and symlink escapes.

Commands are run as explicit argument arrays. The app allows only the DotMatch
doctor checks and AssaySpec commands needed for local assay work:

```text
dotmatch --version
dotmatch dist ACGT AGGT
dotmatch assay infer ...
dotmatch assay check assay.toml
dotmatch assay plan assay.toml
dotmatch assay run assay.toml
dotmatch assay autopsy assay.toml --out-dir autopsy/
```

## Workflow

1. Choose a workspace and verify the local DotMatch executable.
2. Draft a count AssaySpec or run inference from local FASTQ and target files.
3. Review candidate windows, QC warnings, and generated TOML.
4. Check and plan the spec before running.
5. Run the assay locally, then review `assay_report.html`, sample QC, audit
   outputs, autopsy outputs, and workflow-friendly manifests.

The desktop app does not change target tables or silently promote draft specs.
Users remain responsible for reviewing inferred offsets and assay settings
before running production analyses.

## Threat Model

The first Workbench release focuses on local safety:

- workspace confinement for all file paths;
- command allowlisting for DotMatch-only actions;
- local-only data handling with no network service dependency;
- escaped display of paths, warnings, command output, and file values;
- command logs under `<workspace>/.dotmatch/workbench/`;
- no raw FASTQ copying into Workbench state.

The CLI remains the stable automation interface for servers, workflow managers,
and package-channel users.
