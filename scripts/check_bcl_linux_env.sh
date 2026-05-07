#!/bin/sh
set -eu

echo "platform	$(uname -s)-$(uname -m)"
echo "cpu	$(sysctl -n machdep.cpu.brand_string 2>/dev/null || grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ //')"

for tool in bcl-convert bcl2fastq cuda-demux; do
  if command -v "$tool" >/dev/null 2>&1; then
    printf '%s\t%s\n' "$tool" "$(command -v "$tool")"
  else
    printf '%s\tnot_installed\n' "$tool"
  fi
done

cat <<'EOF'

Publication-grade raw-BCL comparator runs should be executed on Linux x86_64:

  make fetch-10x-bcl-demo
  DOTMATCH_BCL_THREADS=8 make bench-bcl-10x

For user-supplied real classic-BCL or CBCL runs:

  DOTMATCH_BCL_RUN_FOLDER=/data/run \
  DOTMATCH_BCL_SAMPLE_SHEET=/data/SampleSheet.csv \
  DOTMATCH_BCL_THREADS=8 \
  make bench-bcl-real

Install BCL Convert and/or bcl2fastq according to Illumina's license and
platform requirements. DotMatch does not redistribute proprietary Illumina
software. CUDA-Demux requires a CUDA-capable Linux host.
EOF
