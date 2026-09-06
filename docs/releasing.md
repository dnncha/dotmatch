# Publishing a reviewed research alpha

The analysis package is offline. Publication is a separate, explicitly
requested authenticated action. Do not put credentials or private DNA into this
repository. Do not merge the historical staging branch into DotMatch.

## Prepare

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python -m mypy src/editwitness
python -m ruff check src/editwitness --select E4,E7,E9,F
python scripts/check_style.py
python scripts/generate_schemas.py --check
python scripts/release_manifest.py --check
```

Any intentional source change requires inventory regeneration **after review**.
Do not regenerate it automatically to hide unexplained differences. Record
actual results, not a previous release's successful checks.

## Publish source and create a GitHub prerelease

With Git and GitHub CLI installed and already authenticated as the intended
owner, run from the extracted release root:

```bash
python scripts/publish_github.py --public --release --dry-run
python scripts/publish_github.py --public --release
```

The default owner is `dnncha`. The script verifies the source inventory and user,
requires that the target repository not exist, creates `dnncha/editwitness`, and
pushes a fresh independent history. It does not clone or modify DotMatch. Without
`--release` it stops after source publication.

With `--release` it waits for the exact published commit's **push** run of
`ci.yml`, requires a successful conclusion, downloads the CI-built distribution
artifact, verifies SHA256SUMS and creates the alpha tag and prerelease. A PR run,
stale success, existing tag, mismatched artifact or failed CI is not sufficient.
The CI artifact contains the wheel and source distribution built only after the
matrix, coverage and typing checks succeed. The publisher does not publish to
PyPI and does not pretend a checksum is an authenticated signature.

## Partial completion and recovery

If source publication succeeds but CI fails, the public source may exist without
a release. Inspect the error and retained recovery directory. Never retry with a
force push. `--resume` is allowed only when the existing repository has the
requested visibility and **exactly** the reviewed source inventory; it does not
overwrite later work. For a real correction, make a reviewed follow-up commit,
update its inventory, push via normal Git, resolve all CI, and run the publisher
with `--resume --public --release` from that exact reviewed source.

This session supplied only read-only connector operations and no authenticated
CLI. Therefore remote publication could not be performed here. The local
publisher helper tests and dry run are not a live GitHub end-to-end test.

## Package-index release

PyPI namespace ownership, trusted publishing, package provenance attestations
and a clean index installation test are separate tasks. They are not configured
or performed by this script. Do not advertise an index command until the actual
index project and version have been verified.
