# EditWitness — research-alpha source delivery

**Know what your CRISPR assay can—and cannot—see.**

This isolated branch delivers EditWitness 0.1.0a1, a separate Apache-2.0 bioinformatics package. It is NOT a change to DotMatch and must NOT be merged into DotMatch main.

The branch bootstrap verifies a SHA-256 checked source transport, materializes the package and its documentation, and runs the tests and strict type checks. After materialization, this README is replaced by the complete project README.

The software passed 572 local tests and a branch-aware covered run. This verifies an explicitly bounded software model; no biological validation or clinical suitability is claimed.

The current GitHub connection cannot create repositories. The complete source includes a conservative publisher: `python scripts/publish_github.py --public`. Run it with an authenticated GitHub CLI to create `dnncha/editwitness` with clean, independent history. It refuses existing repositories and never modifies DotMatch main. Nothing has been published on PyPI.

See STAGING.md for delivery details.
