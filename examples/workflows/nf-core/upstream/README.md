# DotMatch modules prepared for nf-core/modules PR

This directory contains the **exact tree** that can be copied into a fork of https://github.com/nf-core/modules for a PR.

Structure matches nf-core convention:
- modules/nf-core/dotmatch/crispr_count/main.nf
- modules/nf-core/dotmatch/crispr_count/meta.yml
- modules/nf-core/dotmatch/crispr_count/tests/main.nf.test
- (same for assay_run)

## How to use for PR

1. Clone nf-core/modules
2. cp -r <this-dir>/modules/nf-core/dotmatch/* modules/nf-core/dotmatch/
3. Update container hashes to the exact ones from the Bioconda release you are pinning (see packaging/bioconda or bioconda-recipes).
4. Run `nf-core modules lint dotmatch/crispr_count` and fix any issues.
5. Adapt tests if needed (the current tests reference this repo's fixtures; for upstream, nf-core often uses a test-datasets submodule or small embedded data).
6. Submit PR following nf-core contributing guide.
7. Once merged, record in this repo's docs/workflow-adoption.json with the PR URL as adoption_url and a usage example or the module page as evidence_url.

See the parent ../README.md for the full contribution guide and checklist.

These modules were prepped with:
- nf-core style container conditional
- cpus/memory/time resource labels
- when directive
- full ext.args support
- threads forwarding (auto CPU detection)
- versions.yml
- nf-test candidates using self-contained test data/ (copied fixtures for standalone PR)
- Standard nf-core meta.yml fields (authors, maintainers, license)
- Comprehensive meta.yml

After PR, update `make workflow-adoption-status` will pass and scientific claims can reference the integration.

This is the actionable step for "nf-core module upstreaming" in the roadmap.