# Ecosystem Status

This page records the public state of DotMatch distribution, registry, and
workflow integrations. It separates source readiness from upstream submission,
acceptance, release, and installability. States were checked against the linked
provider records on 2026-09-02.

| Surface | Version or revision | Local state | Public state | Installable from that surface? |
| --- | --- | --- | --- | --- |
| [PyPI](https://pypi.org/project/dotmatch/0.4.0/) | 0.4.0 | Tagged release workflow passed | Published with an sdist, a universal macOS wheel, and repaired manylinux/musllinux wheels for x86_64 and aarch64 | Yes: `python3 -m pip install dotmatch==0.4.0` |
| [Bioconda](https://anaconda.org/bioconda/dotmatch) | 0.2.2 | Recipe handoff for 0.4.0 is open | Released for linux-64, linux-aarch64, osx-64, and osx-arm64 at 0.2.2; [0.4.0 PR #68663](https://github.com/bioconda/bioconda-recipes/pull/68663) is the preferred update | Yes at 0.2.2; 0.4.0 awaits Bioconda merge |
| [Spack packages #6191](https://github.com/spack/spack-packages/pull/6191) | 0.2.2 | Package recipe reviewed upstream | Merged into `spack-packages` `develop`; [recipe present](https://github.com/spack/spack-packages/blob/develop/repos/spack_repo/builtin/packages/py_dotmatch/package.py) | Not recorded as released here; no dated Spack release containing the recipe was verified |
| [Galaxy IUC #8336](https://github.com/galaxyproject/tools-iuc/pull/8336) | `e936bbce3577492ff6c12c83d29534213dcb6ce6` | Four wrappers and fixtures are available; review comments addressed | Submitted, open; matching-mode select, profile 25.0, and stricter asserts pushed | No; IUC merge and ToolShed publication are separate gates |
| [nf-core/modules #12156](https://github.com/nf-core/modules/pull/12156) | `77e849b86cae10557a8b17f9c86fa87f6833ece2` | Scoped `crispr_count` module and tests are available | Submitted, open | No; merge and release have not occurred |
| [Snakemake wrappers #5825](https://github.com/snakemake/snakemake-wrappers/pull/5825) | `v9.17.1` | Local workflow and upstream wrapper fixtures are available | Accepted, merged, and present on public tags `v9.17.0` and `v9.17.1` | Yes from the released wrapper tag |
| [MultiQC #3629](https://github.com/MultiQC/MultiQC/pull/3629) | `166f94ce70f2bc1fdbc94f460a4c857511bf1416` | Upstream module is submitted; DotMatch includes the packaged-plugin search registration fix | Submitted, open; fixture PR merged and module CI green | The DotMatch plugin is installable from PyPI; no upstream MultiQC release includes the submitted module yet |
| [bio.tools](https://bio.tools/?q=dotmatch) | draft metadata | `docs/registries/biotools.yml` is ready for review | Exact tool API record returns 404; no accepted DotMatch record found | No |
| [WorkflowHub](https://workflowhub.eu/workflows) | none | Runnable repository workflow examples exist | No exact DotMatch workflow record found | No |

## Reproducible install checks

PyPI:

```bash
python3 -m venv .dotmatch-pypi-smoke
.dotmatch-pypi-smoke/bin/pip install dotmatch==0.4.0
.dotmatch-pypi-smoke/bin/dotmatch --version
.dotmatch-pypi-smoke/bin/dotmatch dist ACGT AGGT
```

Bioconda:

```bash
conda create -y -p .dotmatch-bioconda-smoke \
  -c conda-forge -c bioconda dotmatch=0.2.2
conda run -p .dotmatch-bioconda-smoke dotmatch --version
conda run -p .dotmatch-bioconda-smoke dotmatch dist ACGT AGGT
```

The MultiQC plugin fix can be checked from this source tree without a custom
search-pattern configuration:

```bash
python3 -m venv .dotmatch-multiqc-smoke
.dotmatch-multiqc-smoke/bin/pip install ".[multiqc]"
.dotmatch-multiqc-smoke/bin/multiqc \
  examples/workflows/multiqc/data --module dotmatch \
  -o .dotmatch-multiqc-report
```

An open pull request is not an accepted or released integration. The Snakemake
wrapper is released on `v9.17.1`. A merged repository recipe is not described as
available from a stable package-manager release until that release is verified.
`docs/workflow-adoption.json` remains the gate for accepted workflow-manager
records.
