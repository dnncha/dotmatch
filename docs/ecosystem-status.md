# Ecosystem Status

This page records the public state of DotMatch distribution, registry, and
workflow integrations. It separates source readiness from upstream submission,
acceptance, release, and installability. States were checked against the linked
provider records on 2026-08-28.

| Surface | Version or revision | Local state | Public state | Installable from that surface? |
| --- | --- | --- | --- | --- |
| [PyPI](https://pypi.org/project/dotmatch/0.2.2/) | 0.2.2 | Clean-wheel CLI smoke passes | Released | Yes: `python3 -m pip install dotmatch==0.2.2` |
| [Bioconda](https://anaconda.org/bioconda/dotmatch) | 0.2.2 | Clean Conda CLI smoke passes | Released for linux-64, linux-aarch64, osx-64, and osx-arm64 | Yes: `conda create -n dotmatch -c conda-forge -c bioconda dotmatch=0.2.2` |
| [Spack packages #6191](https://github.com/spack/spack-packages/pull/6191) | 0.2.2 | Package recipe reviewed upstream | Merged into `spack-packages` `develop`; [recipe present](https://github.com/spack/spack-packages/blob/develop/repos/spack_repo/builtin/packages/py_dotmatch/package.py) | Not recorded as released here; no dated Spack release containing the recipe was verified |
| [Galaxy IUC #8336](https://github.com/galaxyproject/tools-iuc/pull/8336) | `b547330d0b3f9fb38368eac3a94fd84098d51031` | Four wrappers and fixtures are available | Submitted, open; lint and containerized Planemo tests pass | No; IUC merge and ToolShed publication are separate gates |
| [nf-core/modules #12156](https://github.com/nf-core/modules/pull/12156) | `77e849b86cae10557a8b17f9c86fa87f6833ece2` | Scoped `crispr_count` module and tests are available | Submitted, open | No; merge and release have not occurred |
| [Snakemake wrappers #5825](https://github.com/snakemake/snakemake-wrappers/pull/5825) | `a72d2bb8bdeb97d8ac506cbc249925f969888b6f` | Local workflow and upstream wrapper fixtures are available; Black passes locally | Submitted, open; Code quality and Tests runs await upstream maintainer approval | No; acceptance and release have not occurred |
| [MultiQC #3629](https://github.com/MultiQC/MultiQC/pull/3629) | `166f94ce70f2bc1fdbc94f460a4c857511bf1416` | Upstream module is submitted; this source branch also fixes packaged-plugin search registration | Submitted, open; eight checks pass and one additional Python 3.9 check is cancelled | No upstream release; the packaged-plugin discovery fix is source-only until the next DotMatch release |
| [bio.tools](https://bio.tools/?q=dotmatch) | draft metadata | `docs/registries/biotools.yml` is ready for review | Exact tool API record returns 404; no accepted DotMatch record found | No |
| [WorkflowHub](https://workflowhub.eu/workflows) | none | Runnable repository workflow examples exist | No exact DotMatch workflow record found | No |

## Reproducible install checks

PyPI:

```bash
python3 -m venv /tmp/dotmatch-pypi-smoke
/tmp/dotmatch-pypi-smoke/bin/pip install dotmatch==0.2.2
/tmp/dotmatch-pypi-smoke/bin/dotmatch --version
/tmp/dotmatch-pypi-smoke/bin/dotmatch dist ACGT AGGT
```

Bioconda:

```bash
conda create -y -p /tmp/dotmatch-bioconda-smoke \
  -c conda-forge -c bioconda dotmatch=0.2.2
conda run -p /tmp/dotmatch-bioconda-smoke dotmatch --version
conda run -p /tmp/dotmatch-bioconda-smoke dotmatch dist ACGT AGGT
```

The MultiQC plugin fix can be checked from this source tree without a custom
search-pattern configuration:

```bash
python3 -m venv /tmp/dotmatch-multiqc-smoke
/tmp/dotmatch-multiqc-smoke/bin/pip install ".[multiqc]"
/tmp/dotmatch-multiqc-smoke/bin/multiqc \
  examples/workflows/multiqc/data --module dotmatch \
  -o /tmp/dotmatch-multiqc-report
```

An open pull request is not an accepted or released integration. A merged
repository recipe is not described as available from a stable package-manager
release until that release is verified. `docs/workflow-adoption.json` remains
the gate for accepted workflow-manager records.
