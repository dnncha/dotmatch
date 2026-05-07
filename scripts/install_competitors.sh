#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
PREFIX="${DOTMATCH_COMPETITOR_PREFIX:-$ROOT/build/competitor-env}"
GUIDE_ROOT="${DOTMATCH_GUIDE_COUNTER_ROOT:-$ROOT/build/guide-counter}"
MAMBA_ROOT="${DOTMATCH_MICROMAMBA_ROOT:-$ROOT/build/micromamba}"

bootstrap_micromamba() {
  if [ "${DOTMATCH_BOOTSTRAP_MICROMAMBA:-1}" = "0" ]; then
    return 1
  fi
  mkdir -p "$MAMBA_ROOT/bin"
  if [ -x "$MAMBA_ROOT/bin/micromamba" ]; then
    CONDA_EXE="$MAMBA_ROOT/bin/micromamba"
    return 0
  fi
  os=$(uname -s)
  arch=$(uname -m)
  case "$os:$arch" in
    Linux:x86_64) platform=linux-64 ;;
    Linux:aarch64|Linux:arm64) platform=linux-aarch64 ;;
    Darwin:arm64) platform=osx-arm64 ;;
    Darwin:x86_64) platform=osx-64 ;;
    *)
      echo "no bundled micromamba platform mapping for $os/$arch" >&2
      return 1
      ;;
  esac
  tmp="${MAMBA_ROOT}/micromamba.tar.bz2"
  url="https://micro.mamba.pm/api/micromamba/${platform}/latest"
  echo "Installing local micromamba from $url" >&2
  if command -v curl >/dev/null 2>&1; then
    curl -L "$url" -o "$tmp"
  elif command -v python3 >/dev/null 2>&1; then
    python3 - "$url" "$tmp" <<'PY'
import sys, urllib.request
urllib.request.urlretrieve(sys.argv[1], sys.argv[2])
PY
  else
    echo "curl or python3 is required to bootstrap micromamba" >&2
    return 1
  fi
  tar -xjf "$tmp" -C "$MAMBA_ROOT" bin/micromamba
  CONDA_EXE="$MAMBA_ROOT/bin/micromamba"
  return 0
}

if command -v mamba >/dev/null 2>&1; then
  CONDA_EXE=mamba
elif command -v conda >/dev/null 2>&1; then
  CONDA_EXE=conda
elif bootstrap_micromamba; then
  :
else
  echo "mamba/conda is required to install benchmark competitors, or allow local micromamba bootstrap" >&2
  exit 1
fi

if [ ! -x "$PREFIX/bin/mageck" ] || [ ! -x "$PREFIX/bin/cutadapt" ] || [ ! -x "$PREFIX/bin/bowtie2" ]; then
  "$CONDA_EXE" create -y -p "$PREFIX" \
    -c conda-forge -c bioconda \
    python=3.11 \
    mageck=0.5.9.5 \
    cutadapt=5.2 \
    bowtie2=2.5.4 \
    rust
fi

if [ ! -x "$GUIDE_ROOT/bin/guide-counter" ]; then
  PATH="$PREFIX/bin:$PATH" "$PREFIX/bin/cargo" install guide-counter --version 0.1.3 --root "$GUIDE_ROOT" --locked
fi

cat <<EOF
DOTMATCH_COMPETITOR_BIN=$PREFIX/bin
DOTMATCH_GUIDE_COUNTER_BIN=$GUIDE_ROOT/bin

To run the public CRISPR gauntlet with installed tools:

PATH="$GUIDE_ROOT/bin:$PREFIX/bin:\$PATH" python3 scripts/run_public_crispr_benchmark.py --small --run-mageck --run-cutadapt --run-bowtie2 --run-guide-counter
EOF
