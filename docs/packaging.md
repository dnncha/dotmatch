# Packaging and installation

DotMatch is distributed through PyPI, GitHub releases, and containers. A
Bioconda recipe is maintained separately and can arrive later than a PyPI
release while its build is reviewed.

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
python3 -m pip install dotmatch==0.2.2
```

The PyPI page is <https://pypi.org/project/dotmatch/>.

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
docker pull ghcr.io/dnncha/dotmatch:0.2.2
docker run --rm ghcr.io/dnncha/dotmatch:0.2.2 dist ACGT AGGT
```

BioContainers images are generated after the corresponding Bioconda package is
published. Their tags include the Bioconda build number, so use the tag shown on
the package page rather than guessing it.

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
repairs Linux wheels, checks their contents, publishes through PyPI trusted
publishing, and uploads the same artifacts to the GitHub release. See
[Release process](release-process.md) for the complete maintainer sequence.

### Verify published artifacts

PyPI trusted publishing uploads a source distribution plus a macOS wheel and
repaired manylinux/musllinux wheels. The channel check rejects raw
`linux_x86_64` wheels, creates a clean virtual environment, and runs an exact
install such as:

```bash
pip install dotmatch==<version>
```

Check the GitHub Container Registry image with:

```bash
docker run --rm ghcr.io/dnncha/dotmatch:v<version> --version
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
