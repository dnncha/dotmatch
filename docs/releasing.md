# Releasing EditWitness

## Current distribution status

The delivered source and wheel are research-alpha artifacts. A configured workflow
is not a published release. `BUILD_STATUS.md` records what was actually checked.
Do not claim PyPI, Bioconda, a DOI, hosted documentation, or remote CI is complete
until those systems confirm it.

## Standalone GitHub publication

The ChatGPT GitHub connection used for this build can write to existing
repositories but cannot create a repository. If the code is supplied through an
isolated staging branch in `dnncha/dotmatch`, do **not merge that branch into
DotMatch**. EditWitness is an independent project.

With the source extracted locally and GitHub CLI already authenticated as
`dnncha`, first inspect the dry run:

```bash
python scripts/publish_github.py --public --dry-run
python scripts/publish_github.py --public
```

The script verifies the checked source inventory, refuses an existing target,
creates a fresh git history, and uses `gh repo create` to publish
`dnncha/editwitness`. It does not transfer DotMatch history, force-push, modify
existing repositories, or publish to PyPI. Choose `--private` instead of
`--public` to make the new repository private. Visibility is always explicit.

If creation succeeds but pushing fails, the script retains the temporary source
repository and prints its recovery location. Inspect the remote before retrying;
do not delete an existing repository to make a retry pass.

After reviewed source edits, regenerate the inventory with
`python scripts/release_manifest.py`. Do not do this merely to suppress an
unexpected integrity failure.

## Release checklist

1. Run tests, schema checks, source hygiene and strict typing. Build and install
   the wheel outside the source tree. Record actual environment and outcomes.
2. Check local and remote files for secrets or patient-derived data. Keep
   synthetic fixtures explicitly labeled. Enable private vulnerability reporting
   when the standalone repository is available.
3. Execute the Linux/macOS/Windows CI matrix. Resolve failures before tagging.
   Add repository topics and links only after the target exists.
4. Verify ownership and availability of the `editwitness` package namespace.
   Configure reviewed trusted publishing; do not store long-lived credentials in
   the repo. Publish an explicitly alpha version, never a misleading stable 1.0.
5. Check install instructions on a clean machine and archive versioned artifacts
   with SHA-256 checksums. Cite a real commit/version; register a DOI only if an
   actual archive integration has created one.
6. Treat empirical validation as a separate scientific gate. Passing packaging
   and CI cannot satisfy it.

The build workflow produces artifacts only. It does not contain a hidden
publishing step or assume package-index credentials exist.
