# Distribution Submission Dossier

This dossier turns the prepared release artifacts into a concrete public
distribution handoff. It is a pre-release checklist, not a public distribution claim.

## Required Checks

Run these before creating or pushing an immutable release tag:

```bash
make release-ready
make python-package-test
make bioconda-recipe-ready
make distribution-record-ready
make distribution-submission-ready
```

Run the post-release public-channel gate only after PyPI, Bioconda, GHCR,
BioContainers, and Zenodo are public:

```bash
make distribution-channels
```

## Submission Targets

- PyPI trusted publishing: push the immutable `v0.1.0` tag and publish the
  source distribution built and verified by
  `scripts/check_python_wheel.py --sdist-only --out-dir dist` plus the repaired manylinux/musllinux wheels from the `dotmatch-linux-repaired-wheels` artifact
  in `.github/workflows/release.yml`; keep raw `linux_x86_64` wheels out of
  PyPI.
- GHCR: verify `ghcr.io/dnncha/dotmatch:v0.1.0`, OCI labels,
  `dotmatch --version`, and `dotmatch dist ACGT AGGT` from the pushed image.
- Bioconda: replace `REPLACE_WITH_RELEASE_TARBALL_SHA256` in
  `packaging/bioconda/meta.yaml` with the immutable GitHub release-tarball hash,
  then submit the recipe to `bioconda-recipes`.
- BioContainers: after Bioconda acceptance, verify
  `quay.io/biocontainers/dotmatch` has a matching tag and runtime CLI smoke
  test.
- Zenodo: archive the immutable GitHub release, mint the release DOI, add the
  DOI to `CITATION.cff` and release metadata, and keep citation wording aligned
  with `docs/methods-and-citation.md`.

## Record Update

After every public channel is live, update `docs/distribution-release.json`:

- set the overall status to `released`;
- set every channel status to `verified`;
- add stable public HTTPS URLs in `public_url`;
- add evidence or review links in `evidence_url`;
- add `verified_date` in `YYYY-MM-DD` format.

Use unique channel IDs and replace all placeholders with non-example public
URLs before setting a channel to `verified`.

Do not set status to `released` until stable public HTTPS URLs exist for every
channel and `make distribution-channels` passes from a clean environment.
