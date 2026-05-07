#!/bin/sh
set -eu

for tool in bcl-convert bcl2fastq cuda-demux; do
  if command -v "$tool" >/dev/null 2>&1; then
    printf '%s\t%s\n' "$tool" "$(command -v "$tool")"
  else
    printf '%s\tnot_installed\n' "$tool"
  fi
done

cat <<'EOF'

DotMatch does not bundle proprietary Illumina tools. Install BCL Convert or
bcl2fastq separately, then rerun:

  python3 scripts/bench_bcl_demux.py --run-folder /path/to/run --sample-sheet /path/to/SampleSheet.csv --detect-competitors
  python3 scripts/generate_bcl_demux_report.py
EOF
