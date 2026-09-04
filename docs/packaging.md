# Packaging and installation

DotMatch is distributed through PyPI, GitHub releases, Bioconda, and
containers. Release 0.4.1 is the current release on PyPI, GitHub Releases, and
GHCR. Its version-specific Zenodo DOI is pending; the prior archived release is
0.4.0. Bioconda remains at 0.2.2, and its generated BioContainers images can lag
those channels, so check the provider records before pinning a version.

## PyPI

Install the current Python package with:

```bash
python3 -m pip install dotmatch
dotmatch --version
```

The wheel includes the Python package, the `dotmatch` command, the native
library, and the native command-line executable. Supported wheels are published
for Linux and macOS. Python 3.9 or newer is required.

To install an exact release:

```bash
python3 -m pip install dotmatch==<version>
```

The PyPI page is <https://pypi.org/project/dotmatch/>.

The release workflow is configured to build repaired `manylinux` and
`musllinux` wheels for `x86_64` and `aarch64`. Check the release record for
the architectures confirmed for a specific version before pinning it.

## Bioconda

When the required version is available in Bioconda:

```bash
conda create -n dotmatch -c conda-forge -c bioconda dotmatch
conda activate dotmatch
dotmatch --version
```

The recipe covers Linux, Intel macOS, and Apple Silicon (`osx-arm64`). Check the
[Anaconda package page](https://anaconda.org/bioconda/dotmatch) before pinning a
newly tagged version; Bioconda may still be building or reviewing it.

## GitHub releases

Each tagged release contains the source archive, Python distribution files,
checksums, and release notes:

<https://github.com/dnncha/dotmatch/releases>

Verify a downloaded artifact with the matching entry in `checksums.txt` before
installing it outside a package manager.

## Container image

Release images are published to GitHub Container Registry:

```bash
docker pull ghcr.io/dnncha/dotmatch:v<version>
docker run --rm ghcr.io/dnncha/dotmatch:v<version> dist ACGT AGGT
```

The release workflow is configured to publish `linux/amd64` and `linux/arm64`
image manifests and smoke-test both native CLI paths. Check the release record
before pinning a tag.

BioContainers images are generated after the corresponding Bioconda package is
published. Their tags include the Bioconda build number, so use the tag shown on
the package page rather than guessing it.

For the latest BioContainers release currently verified here, one available
Linux/Python 3.11 tag is:

```bash
docker pull quay.io/biocontainers/dotmatch:0.2.2--py311h13f8228_1
docker run --rm quay.io/biocontainers/dotmatch:0.2.2--py311h13f8228_1 dotmatch --version
```

See the [BioContainers package page](https://quay.io/repository/biocontainers/dotmatch)
for the other build tags.

## Build from source

The source build requires a C compiler, `make`, Python 3.9 or newer, and zlib.

```bash
git clone https://github.com/dnncha/dotmatch.git
cd dotmatch
make
python3 -m pip install .
dotmatch --version
```

For a development checkout:

```bash
python3 -m pip install -e .
make test
```

## What the package installs

The Python distribution provides:

- the `dotmatch` Python package;
- the `dotmatch` command;
- the native library and executable used by the Python command;
- the optional `assaycode` compatibility command for assay-specification work;
- entry points for supported workflow integrations.

The optional MultiQC integration in DotMatch 0.3.0 includes the module entry
point and the `before_config` search-pattern hook needed for direct file
discovery with MultiQC 1.35. See `ecosystem-status.md` for the separate state
of the upstream MultiQC submission.

The desktop Workbench is not included. It is maintained in the separate
[`dotmatch-community`](https://github.com/dnncha/dotmatch-community) repository.

## Check an installation

These commands confirm the package version, native library, and CLI path:

```bash
dotmatch --version
dotmatch dist ACGT AGGT
dotmatch leq 1 ACGT AGGT
python3 -c "import dotmatch; print(dotmatch.__version__); print(dotmatch.distance('ACGT', 'AGGT'))"
```

The expected distance is `1`, and `leq` should print `true`.

## Maintainer release checks

Maintainers use the repository checks before tagging a release:

```bash
make pretag-ready
make release-ready
```

The release workflow builds the source distribution and platform wheels,
repairs Linux wheels for `x86_64` and `aarch64`, checks their contents,
publishes through PyPI trusted publishing, and uploads the same artifacts to
the GitHub release. See [Release process](release-process.md) for the complete
maintainer sequence.

### Verify published artifacts

PyPI trusted publishing uploads a source distribution plus a macOS wheel and
repaired manylinux/musllinux wheels. The channel check rejects raw Linux wheels
(`linux_x86_64` or `linux_aarch64`) for the recorded architectures, creates a
clean virtual environment, and runs an exact install such as:

```bash
pip install dotmatch==<version>
```

Check the GitHub Container Registry image with:

```bash
docker run --rm ghcr.io/dnncha/dotmatch:v<version> --version
```

For a multi-architecture release, inspect the image index and confirm both
recorded platforms are present:

```bash
docker buildx imagetools inspect ghcr.io/dnncha/dotmatch:v<version>
```

Check the Bioconda package in a new prefix:

```bash
conda create -p <env> -c conda-forge -c bioconda dotmatch=<version>
```

BioContainers images for DotMatch are generated from the accepted Bioconda
recipe. Their tags include the build number:

```bash
docker run --rm quay.io/biocontainers/dotmatch:<version>--<build> dotmatch --version
```
