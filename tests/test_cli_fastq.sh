#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
DOTMATCH_BIN=${DOTMATCH_BIN:-"$ROOT/dotmatch"}
TMPDIR="${TMPDIR:-/tmp}/dotmatch-cli-$$"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

cat > "$TMPDIR/barcodes.tsv" <<'BARCODES'
bc0	ACGT
bc1	AGGT
bc2	ACGA
bc3	TTTT
BARCODES

cat > "$TMPDIR/reads.fastq" <<'FASTQ'
@r0
ACGTAAAA
+
IIIIIIII
@r1
TTTGAAAA
+
IIIIIIII
@r2
GGGGAAAA
+
IIIIIIII
@r3
AC
+
II
FASTQ

"$DOTMATCH_BIN" fastq-assign \
  --barcodes "$TMPDIR/barcodes.tsv" \
  --reads "$TMPDIR/reads.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 1 \
  --out "$TMPDIR/out.tsv"

cat > "$TMPDIR/expected.tsv" <<'EXPECTED'
read_id	observed_barcode	target_index	target_id	target_seq	best_distance	second_best_distance	match_count	status
r0	ACGT	0	bc0	ACGT	0	1	3	unique
r1	TTTG	3	bc3	TTTT	1	-1	1	unique
r2	GGGG	-1			-1	-1	0	none
r3		-1			-1	-1	0	invalid
EXPECTED

diff -u "$TMPDIR/expected.tsv" "$TMPDIR/out.tsv"

"$DOTMATCH_BIN" fastq-assign \
  --barcodes "$TMPDIR/barcodes.tsv" \
  --reads "$TMPDIR/reads.fastq" \
  --barcode-start 1 \
  --barcode-length 4 \
  --k 1 \
  --out "$TMPDIR/offset.tsv"

grep '^r1	TTGA	-1			-1	-1	0	none$' "$TMPDIR/offset.tsv" >/dev/null

gzip -c "$TMPDIR/reads.fastq" > "$TMPDIR/reads.fastq.gz"
"$DOTMATCH_BIN" fastq-assign \
  --barcodes "$TMPDIR/barcodes.tsv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 0 \
  --out "$TMPDIR/gz.tsv"

grep '^r0	ACGT	0	bc0	ACGT	0	-1	1	unique$' "$TMPDIR/gz.tsv" >/dev/null

cat > "$TMPDIR/mageck_seq_header.tsv" <<'TARGETS'
sgRNAID	Seq	gene
guide0	ACGT	GENE0
guide1	TTTT	GENE1
TARGETS

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/mageck_seq_header.tsv" \
  --reads "$TMPDIR/reads.fastq" \
  --sample-label sample \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --format mageck \
  --out "$TMPDIR/mageck_seq_counts.tsv" \
  --summary "$TMPDIR/mageck_seq_summary.json"

grep '^guide0	GENE0	1$' "$TMPDIR/mageck_seq_counts.tsv" >/dev/null
grep '"n_targets": 2' "$TMPDIR/mageck_seq_summary.json" >/dev/null

if "$DOTMATCH_BIN" fastq-assign \
  --barcodes "$TMPDIR/barcodes.tsv" \
  --reads "$TMPDIR/reads.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 2 \
  --out "$TMPDIR/bad.tsv" 2>/dev/null; then
  echo "k=2 should fail for fastq-assign first milestone" >&2
  exit 1
fi

cat > "$TMPDIR/bad.fastq" <<'BADFASTQ'
@bad
ACGT
+
BADFASTQ

if "$DOTMATCH_BIN" fastq-assign \
  --barcodes "$TMPDIR/barcodes.tsv" \
  --reads "$TMPDIR/bad.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 1 \
  --out "$TMPDIR/malformed.tsv" 2>/dev/null; then
  echo "malformed FASTQ should fail" >&2
  exit 1
fi

cat > "$TMPDIR/demux_reads.fastq" <<'DEMUXFASTQ'
@d0
ACGTAAAA
+
IIIIIIII
@d1
TTTGAAAA
+
IIIIIIII
@d2
AGGAAAAA
+
IIIIIIII
@d3
GGGGAAAA
+
IIIIIIII
@d4
AC
+
II
DEMUXFASTQ

mkdir "$TMPDIR/demux"
"$DOTMATCH_BIN" demux \
  --barcodes "$TMPDIR/barcodes.tsv" \
  --reads "$TMPDIR/demux_reads.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 1 \
  --metric hamming \
  --out-dir "$TMPDIR/demux" \
  --summary "$TMPDIR/demux_summary.json" \
  --assignments "$TMPDIR/demux_assignments.tsv" \
  --ambiguous-out "$TMPDIR/demux_ambiguous.fastq" \
  --unmatched-out "$TMPDIR/demux_unmatched.fastq"

grep '^@d0$' "$TMPDIR/demux/bc0.fastq" >/dev/null
grep '^ACGTAAAA$' "$TMPDIR/demux/bc0.fastq" >/dev/null
grep '^@d1$' "$TMPDIR/demux/bc3.fastq" >/dev/null
grep '^@d2$' "$TMPDIR/demux_ambiguous.fastq" >/dev/null
grep '^@d3$' "$TMPDIR/demux_unmatched.fastq" >/dev/null
grep '^@d4$' "$TMPDIR/demux_unmatched.fastq" >/dev/null
grep '"assigned_unique": 2' "$TMPDIR/demux_summary.json" >/dev/null
grep '"ambiguous": 1' "$TMPDIR/demux_summary.json" >/dev/null
grep '"unmatched": 1' "$TMPDIR/demux_summary.json" >/dev/null
grep '"invalid": 1' "$TMPDIR/demux_summary.json" >/dev/null
grep '^d2	AGGA	1	bc1	AGGT	1	-1	2	ambiguous$' "$TMPDIR/demux_assignments.tsv" >/dev/null

python3 - "$TMPDIR/bcl_run" <<'PY'
import gzip
import struct
import sys
from pathlib import Path

root = Path(sys.argv[1])
base = root / "Data" / "Intensities" / "BaseCalls" / "L001"
for cycle in range(1, 9):
    (base / f"C{cycle}.1").mkdir(parents=True, exist_ok=True)

(root / "RunInfo.xml").write_text("""<?xml version=\"1.0\"?>
<RunInfo>
  <Run Id=\"tiny\" Number=\"1\">
    <Flowcell>TEST</Flowcell>
    <Reads>
      <Read Number=\"1\" NumCycles=\"4\" IsIndexedRead=\"N\"/>
      <Read Number=\"2\" NumCycles=\"4\" IsIndexedRead=\"Y\"/>
    </Reads>
    <FlowcellLayout LaneCount=\"1\" SurfaceCount=\"1\" SwathCount=\"1\" TileCount=\"1\"/>
  </Run>
</RunInfo>
""")

(root / "SampleSheet.csv").write_text("""[Header]
IEMFileVersion,4
[Data]
Sample_ID,Sample_Name,index
s1,Sample One,ACGT
s2,Sample Two,AGGT
s3,Sample Three,ACGA
""")

reads = ["TTTT", "CCCC", "GGGG", "AAAA"]
indexes = ["ACGT", "AGGA", "GGGG", "AGGT"]
base_code = {"A": 0, "C": 1, "G": 2, "T": 3}

def write_bcl(path, bases):
    with gzip.open(path, "wb") as fh:
        fh.write(struct.pack("<I", len(bases)))
        for base in bases:
            fh.write(bytes([(40 << 2) | base_code[base]]))

for pos in range(4):
    write_bcl(base / f"C{pos + 1}.1" / "s_1_1101.bcl.gz", [seq[pos] for seq in reads])
for pos in range(4):
    write_bcl(base / f"C{pos + 5}.1" / "s_1_1101.bcl.gz", [seq[pos] for seq in indexes])

with (base / "s_1_1101.filter").open("wb") as fh:
    fh.write(struct.pack("<II", 0, 4))
    fh.write(bytes([1, 1, 1, 0]))
PY

"$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_run" \
  --sample-sheet "$TMPDIR/bcl_run/SampleSheet.csv" \
  --out-dir "$TMPDIR/bcl_out" \
  --barcode-mismatches 1 \
  --summary "$TMPDIR/bcl_summary.json"

gzip -cd "$TMPDIR/bcl_out/s1_S1_L001_R1_001.fastq.gz" | grep '^TTTT$' >/dev/null
gzip -cd "$TMPDIR/bcl_out/Undetermined_S0_L001_R1_001.fastq.gz" | grep '^CCCC$' >/dev/null
gzip -cd "$TMPDIR/bcl_out/Undetermined_S0_L001_R1_001.fastq.gz" | grep '^GGGG$' >/dev/null
grep '^s1,1,1$' "$TMPDIR/bcl_out/Demultiplex_Stats.csv" >/dev/null
grep '^Undetermined,2,2$' "$TMPDIR/bcl_out/Demultiplex_Stats.csv" >/dev/null
grep '^AGGA,1$' "$TMPDIR/bcl_out/Top_Unknown_Barcodes.csv" >/dev/null
grep '^GGGG,1$' "$TMPDIR/bcl_out/Top_Unknown_Barcodes.csv" >/dev/null
grep '"assigned_reads": 1' "$TMPDIR/bcl_summary.json" >/dev/null
grep '"undetermined_reads": 2' "$TMPDIR/bcl_summary.json" >/dev/null
grep '"filtered_clusters": 1' "$TMPDIR/bcl_summary.json" >/dev/null

"$DOTMATCH_BIN" bcl-validate \
  --dotmatch-out "$TMPDIR/bcl_out" \
  --truth-out "$TMPDIR/bcl_out" | grep '"mismatched_fastq_files": 0' >/dev/null

cat > "$TMPDIR/bcl_run/SampleSheet.aliases.csv" <<'SHEET'
[Header]
IEMFileVersion,4
[Data]
Sample_ID,Sample_Name,index
s1,Sample One,ACGT
s1,Sample One,AGGA
s2,Sample Two,ACGA
SHEET

"$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_run" \
  --sample-sheet "$TMPDIR/bcl_run/SampleSheet.aliases.csv" \
  --out-dir "$TMPDIR/bcl_alias_out" \
  --barcode-mismatches 0 \
  --summary "$TMPDIR/bcl_alias_summary.json"

gzip -cd "$TMPDIR/bcl_alias_out/s1_S1_L001_R1_001.fastq.gz" | grep '^TTTT$' >/dev/null
gzip -cd "$TMPDIR/bcl_alias_out/s1_S1_L001_R1_001.fastq.gz" | grep '^CCCC$' >/dev/null
test ! -e "$TMPDIR/bcl_alias_out/s1_S2_L001_R1_001.fastq.gz"
grep '^s1,2,2$' "$TMPDIR/bcl_alias_out/Demultiplex_Stats.csv" >/dev/null
grep '"assigned_reads": 2' "$TMPDIR/bcl_alias_summary.json" >/dev/null

python3 - "$TMPDIR/bcl_pe_run" <<'PY'
import gzip
import struct
import sys
from pathlib import Path

root = Path(sys.argv[1])
base = root / "Data" / "Intensities" / "BaseCalls" / "L001"
for cycle in range(1, 7):
    (base / f"C{cycle}.1").mkdir(parents=True, exist_ok=True)

(root / "RunInfo.xml").write_text("""<?xml version=\"1.0\"?>
<RunInfo>
  <Run Id=\"pe\" Number=\"1\">
    <Flowcell>TEST</Flowcell>
    <Reads>
      <Read Number=\"1\" NumCycles=\"2\" IsIndexedRead=\"N\"/>
      <Read Number=\"2\" NumCycles=\"2\" IsIndexedRead=\"Y\"/>
      <Read Number=\"3\" NumCycles=\"2\" IsIndexedRead=\"N\"/>
    </Reads>
    <FlowcellLayout LaneCount=\"1\" SurfaceCount=\"1\" SwathCount=\"1\" TileCount=\"1\"/>
  </Run>
</RunInfo>
""")
(root / "SampleSheet.csv").write_text("""[Header]
IEMFileVersion,4
[Data]
Sample_ID,Sample_Name,index
s1,Sample One,AA
""")

base_code = {"A": 0, "C": 1, "G": 2, "T": 3}
cycles = ["AG", "CT", "AC", "AC", "TA", "GC"]

def write_bcl(path, bases):
    with gzip.open(path, "wb") as fh:
        fh.write(struct.pack("<I", len(bases)))
        for base in bases:
            fh.write(bytes([(40 << 2) | base_code[base]]))

for i, bases in enumerate(cycles, start=1):
    write_bcl(base / f"C{i}.1" / "s_1_1101.bcl.gz", list(bases))
with (base / "s_1_1101.filter").open("wb") as fh:
    fh.write(struct.pack("<II", 0, 2))
    fh.write(bytes([1, 1]))
PY

"$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_pe_run" \
  --sample-sheet "$TMPDIR/bcl_pe_run/SampleSheet.csv" \
  --out-dir "$TMPDIR/bcl_pe_out" \
  --barcode-mismatches 0 \
  --emit-index-fastqs \
  --summary "$TMPDIR/bcl_pe_summary.json"

gzip -cd "$TMPDIR/bcl_pe_out/s1_S1_L001_R1_001.fastq.gz" | grep '^AC$' >/dev/null
gzip -cd "$TMPDIR/bcl_pe_out/s1_S1_L001_R2_001.fastq.gz" | grep '^TG$' >/dev/null
gzip -cd "$TMPDIR/bcl_pe_out/s1_S1_L001_I1_001.fastq.gz" | grep '^AA$' >/dev/null
gzip -cd "$TMPDIR/bcl_pe_out/Undetermined_S0_L001_R1_001.fastq.gz" | grep '^GT$' >/dev/null
gzip -cd "$TMPDIR/bcl_pe_out/Undetermined_S0_L001_R2_001.fastq.gz" | grep '^AC$' >/dev/null
gzip -cd "$TMPDIR/bcl_pe_out/Undetermined_S0_L001_I1_001.fastq.gz" | grep '^CC$' >/dev/null
grep '^s1,1,1,1$' "$TMPDIR/bcl_pe_out/Demultiplex_Stats.csv" >/dev/null
grep '^Undetermined,1,1,1$' "$TMPDIR/bcl_pe_out/Demultiplex_Stats.csv" >/dev/null

"$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_pe_run" \
  --sample-sheet "$TMPDIR/bcl_pe_run/SampleSheet.csv" \
  --out-dir "$TMPDIR/bcl_pe_threads_out" \
  --barcode-mismatches 0 \
  --emit-index-fastqs \
  --threads 2 \
  --summary "$TMPDIR/bcl_pe_threads_summary.json"

for fq in s1_S1_L001_R1_001.fastq.gz s1_S1_L001_R2_001.fastq.gz s1_S1_L001_I1_001.fastq.gz Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz Undetermined_S0_L001_I1_001.fastq.gz; do
  gzip -cd "$TMPDIR/bcl_pe_out/$fq" > "$TMPDIR/serial.fastq"
  gzip -cd "$TMPDIR/bcl_pe_threads_out/$fq" > "$TMPDIR/threaded.fastq"
  diff -u "$TMPDIR/serial.fastq" "$TMPDIR/threaded.fastq"
done
grep '"effective_threads": 2' "$TMPDIR/bcl_pe_threads_summary.json" >/dev/null

cat > "$TMPDIR/targets.csv" <<'TARGETS'
id,gRNA.sequence,Gene
bc0,ACGT,G0
bc1,AGGT,G1
bc2,ACGA,G2
bc3,TTTT,G3
TARGETS

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label sample1 \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --out "$TMPDIR/counts.tsv" \
  --assignments "$TMPDIR/assignments.tsv" \
  --summary "$TMPDIR/summary.json" \
  --report "$TMPDIR/report.html" \
  --sample-qc "$TMPDIR/sample_qc.tsv" \
  --target-counts-long "$TMPDIR/target_counts.long.tsv" \
  --ambiguous report \
  --ambiguous-out "$TMPDIR/ambiguous.tsv" \
  --unmatched-out "$TMPDIR/unmatched.tsv"

grep '^bc0	ACGT	G0	1	1	0	0	0	0	1$' "$TMPDIR/counts.tsv" >/dev/null
grep '^bc3	TTTT	G3	0	0	1	0	0	0	1$' "$TMPDIR/counts.tsv" >/dev/null
grep '"assigned_unique": 2' "$TMPDIR/summary.json" >/dev/null
grep '"metric": "levenshtein"' "$TMPDIR/summary.json" >/dev/null
grep '"library_covered_targets": 2' "$TMPDIR/summary.json" >/dev/null
grep '^sample1	r2	GGGG	-1			-1	-1	0	none	none$' "$TMPDIR/unmatched.tsv" >/dev/null
grep '^sample1	' "$TMPDIR/sample_qc.tsv" | grep '	4	3	2	1	1	1	0	0	' >/dev/null
grep '^sample1	bc0	G0	ACGT	1	0	0	0	0	1	1$' "$TMPDIR/target_counts.long.tsv" >/dev/null
grep '<title>DotMatch Report</title>' "$TMPDIR/report.html" >/dev/null
grep 'sample1' "$TMPDIR/report.html" >/dev/null
grep 'Assignment rate' "$TMPDIR/report.html" >/dev/null
grep 'Library coverage' "$TMPDIR/report.html" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label sample_lev \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric levenshtein \
  --format mageck \
  --out "$TMPDIR/counts_lev.tsv" \
  --summary "$TMPDIR/summary_lev.json"

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label sample_lev \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric levenshtein \
  --format mageck \
  --threads 2 \
  --out "$TMPDIR/counts_lev_threads.tsv" \
  --summary "$TMPDIR/summary_lev_threads.json"

diff -u "$TMPDIR/counts_lev.tsv" "$TMPDIR/counts_lev_threads.tsv"
grep '"read_threads": 2' "$TMPDIR/summary_lev_threads.json" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label sample_hamming \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --format mageck \
  --out "$TMPDIR/counts_hamming.tsv" \
  --summary "$TMPDIR/summary_hamming.json"

grep '^bc0	G0	1$' "$TMPDIR/counts_hamming.tsv" >/dev/null
grep '^bc3	G3	1$' "$TMPDIR/counts_hamming.tsv" >/dev/null
grep '"count_engine": "hamming_lookup_direct_single_offset"' "$TMPDIR/summary_hamming.json" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label sample_hamming \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --format mageck \
  --threads 2 \
  --out "$TMPDIR/counts_hamming_threads.tsv" \
  --summary "$TMPDIR/summary_hamming_threads.json"

diff -u "$TMPDIR/counts_hamming.tsv" "$TMPDIR/counts_hamming_threads.tsv"
grep '"read_threads": 2' "$TMPDIR/summary_hamming_threads.json" >/dev/null

python3 - "$TMPDIR/long_header.fastq.gz" <<'PY'
import gzip
import sys

with gzip.open(sys.argv[1], "wt") as fh:
    fh.write("@" + ("h" * 9000) + "\n")
    fh.write("ACGTAAAA\n")
    fh.write("+\n")
    fh.write("IIIIIIII\n")
PY

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/long_header.fastq.gz" \
  --sample-label long_header \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --format mageck \
  --out "$TMPDIR/counts_long_header.tsv" \
  --summary "$TMPDIR/summary_long_header.json"

grep '^bc0	G0	1$' "$TMPDIR/counts_long_header.tsv" >/dev/null
grep '"total_reads": 1' "$TMPDIR/summary_long_header.json" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/long_header.fastq.gz" \
  --sample-label long_header_lev \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric levenshtein \
  --format mageck \
  --out "$TMPDIR/counts_long_header_lev.tsv" \
  --summary "$TMPDIR/summary_long_header_lev.json"

grep '^bc0	G0	1$' "$TMPDIR/counts_long_header_lev.tsv" >/dev/null
grep '"total_reads": 1' "$TMPDIR/summary_long_header_lev.json" >/dev/null

"$DOTMATCH_BIN" inspect-unmatched \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --top 10 \
  --out "$TMPDIR/top_unmatched.tsv"

grep '^sequence	count	length	nearest_target	nearest_distance	nearest_edit_class	possible_reason	reverse_complement	revcomp_nearest_target	revcomp_nearest_distance	offset_hint	adapter_hint$' "$TMPDIR/top_unmatched.tsv" >/dev/null
grep '^GGGG	1	4	bc1	2	other	near_known_target_above_k	CCCC	bc0	3		$' "$TMPDIR/top_unmatched.tsv" >/dev/null

"$DOTMATCH_BIN" inspect-unmatched \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --adapter AAAA \
  --top 10 \
  --out "$TMPDIR/top_unmatched_adapter.tsv"

grep '^GGGG	1	4	bc1	2	other	adapter_or_primer_candidate	CCCC	bc0	3		AAAA$' "$TMPDIR/top_unmatched_adapter.tsv" >/dev/null

cat > "$TMPDIR/low_quality.fastq" <<'LOWQUALITY'
@lowq
GGGG
+
!!!!
LOWQUALITY
cat > "$TMPDIR/lowq_targets.tsv" <<'LOWQTARGETS'
bc0	ACGT	G0
LOWQTARGETS

"$DOTMATCH_BIN" inspect-unmatched \
  --targets "$TMPDIR/lowq_targets.tsv" \
  --reads "$TMPDIR/low_quality.fastq" \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --low-quality-threshold 20 \
  --top 10 \
  --out "$TMPDIR/top_unmatched_lowq.tsv"

grep '^GGGG	1	4	bc0	3	other	low_quality_candidate	CCCC	bc0	3		$' "$TMPDIR/top_unmatched_lowq.tsv" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label sample1 \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --out "$TMPDIR/counts_hamming.tsv" \
  --summary "$TMPDIR/summary_hamming.json"

grep '"metric": "hamming"' "$TMPDIR/summary_hamming.json" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label radius \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --ambiguity-policy radius \
  --out "$TMPDIR/counts_radius.tsv" \
  --summary "$TMPDIR/summary_radius.json" \
  --ambiguous report \
  --ambiguous-out "$TMPDIR/ambiguous_radius.tsv"

grep '^bc0	ACGT	G0	1	0	0	0	0	0	0$' "$TMPDIR/counts_radius.tsv" >/dev/null
grep '^bc3	TTTT	G3	0	0	1	0	0	0	1$' "$TMPDIR/counts_radius.tsv" >/dev/null
grep '"ambiguity_policy": "radius"' "$TMPDIR/summary_radius.json" >/dev/null
grep '^radius	r0	ACGT	0	bc0	ACGT	0	1	3	ambiguous	ambiguous$' "$TMPDIR/ambiguous_radius.tsv" >/dev/null

cat > "$TMPDIR/short.fastq" <<'SHORTFASTQ'
@short_del
ACG
+
III
SHORTFASTQ

cat > "$TMPDIR/one_target.tsv" <<'ONETARGET'
bc0	ACGT	G0
ONETARGET

cat > "$TMPDIR/shifted.fastq" <<'SHIFTED'
@shifted
NNACGT
+
IIIIII
SHIFTED

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/one_target.tsv" \
  --reads "$TMPDIR/shifted.fastq" \
  --sample-label shifted \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --auto-offset 2 \
  --auto-offset-sample 1 \
  --out "$TMPDIR/counts_shifted.tsv" \
  --summary "$TMPDIR/summary_shifted.json"

grep '^bc0	ACGT	G0	0	1	0	0	0	0	1$' "$TMPDIR/counts_shifted.tsv" >/dev/null
grep '"selected_target_start": 2' "$TMPDIR/summary_shifted.json" >/dev/null

cat > "$TMPDIR/multi_targets.tsv" <<'MULTITARGETS'
dm0	ACGT	G0
dm1	TTTT	G1
MULTITARGETS

cat > "$TMPDIR/multi_offset.fastq" <<'MULTIFASTQ'
@same_target
ACGTNNACGT
+
IIIIIIIIII
@diff_targets
ACGTNNTTTT
+
IIIIIIIIII
@exact_plus_worse
ACGTNNTTTG
+
IIIIIIIIII
@one_mismatch
NNACGANN
+
IIIIIIII
@one_n
NNACGNNN
+
IIIIIIII
MULTIFASTQ

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/multi_targets.tsv" \
  --reads "$TMPDIR/multi_offset.fastq" \
  --sample-label multi \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --hamming-index precompute \
  --auto-offset 6 \
  --auto-offset-sample 3 \
  --offset-mode multi \
  --offset-min-fraction 0.0 \
  --out "$TMPDIR/counts_multi_offset.tsv" \
  --summary "$TMPDIR/summary_multi_offset.json" \
  --assignments "$TMPDIR/assignments_multi_offset.tsv" \
  --ambiguous report \
  --ambiguous-out "$TMPDIR/ambiguous_multi_offset.tsv"

grep '^dm0	ACGT	G0	0	2	2	0	0	0	4$' "$TMPDIR/counts_multi_offset.tsv" >/dev/null
grep '^dm1	TTTT	G1	0	0	0	0	0	0	0$' "$TMPDIR/counts_multi_offset.tsv" >/dev/null
grep '^multi	same_target	ACGT	0	dm0	ACGT	0	-1	1	unique	exact$' "$TMPDIR/assignments_multi_offset.tsv" >/dev/null
grep '^multi	exact_plus_worse	ACGT	0	dm0	ACGT	0	1	2	unique	exact$' "$TMPDIR/assignments_multi_offset.tsv" >/dev/null
grep '^multi	diff_targets	ACGT	0	dm0	ACGT	0	-1	2	ambiguous	ambiguous$' "$TMPDIR/ambiguous_multi_offset.tsv" >/dev/null
grep '"offset_mode": "multi"' "$TMPDIR/summary_multi_offset.json" >/dev/null
grep '"selected_target_starts": \[0, 1, 2, 3, 4, 5, 6\]' "$TMPDIR/summary_multi_offset.json" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/multi_targets.tsv" \
  --reads "$TMPDIR/multi_offset.fastq" \
  --sample-label multi_fast \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --hamming-index precompute \
  --auto-offset 6 \
  --auto-offset-sample 3 \
  --offset-mode multi \
  --offset-min-fraction 0.0 \
  --out "$TMPDIR/counts_multi_offset_fast.tsv" \
  --summary "$TMPDIR/summary_multi_offset_fast.json"

grep '^dm0	ACGT	G0	0	2	2	0	0	0	4$' "$TMPDIR/counts_multi_offset_fast.tsv" >/dev/null
grep '^dm1	TTTT	G1	0	0	0	0	0	0	0$' "$TMPDIR/counts_multi_offset_fast.tsv" >/dev/null
grep '"offset_detection_strategy": "fused"' "$TMPDIR/summary_multi_offset_fast.json" >/dev/null
grep '"count_engine": "hamming_lookup_direct"' "$TMPDIR/summary_multi_offset_fast.json" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/multi_targets.tsv" \
  --reads "$TMPDIR/multi_offset.fastq" \
  --sample-label multi_radius \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --hamming-index precompute \
  --auto-offset 6 \
  --auto-offset-sample 3 \
  --offset-mode multi \
  --offset-min-fraction 0.0 \
  --ambiguity-policy radius \
  --out "$TMPDIR/counts_multi_offset_radius.tsv" \
  --summary "$TMPDIR/summary_multi_offset_radius.json" \
  --ambiguous report \
  --ambiguous-out "$TMPDIR/ambiguous_multi_offset_radius.tsv"

grep '^dm0	ACGT	G0	0	1	2	0	0	0	3$' "$TMPDIR/counts_multi_offset_radius.tsv" >/dev/null
grep '^dm1	TTTT	G1	0	0	0	0	0	0	0$' "$TMPDIR/counts_multi_offset_radius.tsv" >/dev/null
grep '^multi_radius	same_target	ACGT	0	dm0	ACGT	0	-1	1	unique	exact$' "$TMPDIR/ambiguous_multi_offset_radius.tsv" && exit 1 || true
grep '^multi_radius	diff_targets	ACGT	0	dm0	ACGT	0	-1	2	ambiguous	ambiguous$' "$TMPDIR/ambiguous_multi_offset_radius.tsv" >/dev/null
grep '^multi_radius	exact_plus_worse	ACGT	0	dm0	ACGT	0	1	2	ambiguous	ambiguous$' "$TMPDIR/ambiguous_multi_offset_radius.tsv" >/dev/null
grep '"ambiguity_policy": "radius"' "$TMPDIR/summary_multi_offset_radius.json" >/dev/null

"$DOTMATCH_BIN" validate \
  --targets "$TMPDIR/multi_targets.tsv" \
  --reads "$TMPDIR/multi_offset.fastq" \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --auto-offset 6 \
  --auto-offset-sample 3 \
  --offset-mode multi \
  --offset-min-fraction 0.0 \
  --oracle scan \
  --sample 3 > "$TMPDIR/validate_multi_offset.json"

grep '"mismatches": 0' "$TMPDIR/validate_multi_offset.json" >/dev/null
grep '"offset_mode": "multi"' "$TMPDIR/validate_multi_offset.json" >/dev/null
grep '"selected_target_starts": \[0, 1, 2, 3, 4, 5, 6\]' "$TMPDIR/validate_multi_offset.json" >/dev/null

cat > "$TMPDIR/no_offset_hits.fastq" <<'NOHITSFASTQ'
@nohits
GGGGGG
+
IIIIII
NOHITSFASTQ

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/multi_targets.tsv" \
  --reads "$TMPDIR/no_offset_hits.fastq" \
  --sample-label nohits \
  --target-start 1 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --auto-offset 2 \
  --auto-offset-sample 1 \
  --offset-mode multi \
  --out "$TMPDIR/counts_multi_fallback.tsv" \
  --summary "$TMPDIR/summary_multi_fallback.json"

grep '"selected_target_start": 1' "$TMPDIR/summary_multi_fallback.json" >/dev/null
grep '"selected_target_starts": \[1\]' "$TMPDIR/summary_multi_fallback.json" >/dev/null

"$DOTMATCH_BIN" inspect-unmatched \
  --targets "$TMPDIR/one_target.tsv" \
  --reads "$TMPDIR/shifted.fastq" \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --offset-window 2 \
  --top 10 \
  --out "$TMPDIR/top_unmatched_shifted.tsv"

grep '^NNAC	1	4	bc0	4	other	offset_shift_candidate	GTNN	bc0	4	2	$' "$TMPDIR/top_unmatched_shifted.tsv" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/one_target.tsv" \
  --reads "$TMPDIR/short.fastq" \
  --sample-label indel \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric levenshtein \
  --indel-window 1 \
  --out "$TMPDIR/counts_indel.tsv"

grep '^bc0	ACGT	G0	0	0	0	0	1	0	1$' "$TMPDIR/counts_indel.tsv" >/dev/null

cat > "$TMPDIR/mixed.tsv" <<'MIXED'
short	ACG
long	ACGT
MIXED

if "$DOTMATCH_BIN" count \
  --targets "$TMPDIR/mixed.tsv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label bad \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --out "$TMPDIR/bad_hamming.tsv" 2>/dev/null; then
  echo "hamming metric should reject mixed target lengths" >&2
  exit 1
fi

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label plasmid,esc \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --format mageck \
  --out "$TMPDIR/mageck.tsv"

cat > "$TMPDIR/expected_mageck.tsv" <<'MAGECK'
sgRNA	Gene	plasmid	esc
bc0	G0	1	1
bc1	G1	0	0
bc2	G2	0	0
bc3	G3	1	1
MAGECK
diff -u "$TMPDIR/expected_mageck.tsv" "$TMPDIR/mageck.tsv"

cat > "$TMPDIR/samples.tsv" <<SAMPLES
sample_id	fastq
plasmid	$TMPDIR/reads.fastq
esc	$TMPDIR/reads.fastq.gz
SAMPLES

"$DOTMATCH_BIN" crispr-count \
  --library "$TMPDIR/targets.csv" \
  --samples "$TMPDIR/samples.tsv" \
  --guide-start 0 \
  --guide-length 4 \
  --k 1 \
  --metric levenshtein \
  --threads 2 \
  --out "$TMPDIR/crispr_mageck.tsv" \
  --summary "$TMPDIR/crispr_qc.json"

diff -u "$TMPDIR/expected_mageck.tsv" "$TMPDIR/crispr_mageck.tsv"
grep '"k1_rescued_reads": 1' "$TMPDIR/crispr_qc.json" >/dev/null
grep '"percent_rescued_by_k1": 25.000000' "$TMPDIR/crispr_qc.json" >/dev/null

"$DOTMATCH_BIN" audit \
  --targets "$TMPDIR/targets.csv" \
  --k 1 \
  --out-dir "$TMPDIR/audit"

grep '^targets	4$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '^safe_at_k1	no$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '^risk_pairs_for_k1	3$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '^ambiguous_query_variants_k1	14$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '"audit_mode": "exact"' "$TMPDIR/audit/audit_summary.json" >/dev/null
grep '"k": 1' "$TMPDIR/audit/audit_summary.json" >/dev/null
grep '"safe_at_k1": false' "$TMPDIR/audit/audit_summary.json" >/dev/null
grep '"risk_pairs_for_k1": 3' "$TMPDIR/audit/audit_summary.json" >/dev/null
grep '^bc0	bc1	ACGT	AGGT	1	yes	yes	$' "$TMPDIR/audit/collision_pairs.tsv" >/dev/null
grep '^bc0	ACGT	bc1	1	no	no	2$' "$TMPDIR/audit/target_safety.tsv" >/dev/null
grep '^ACG	2$' "$TMPDIR/audit/ambiguous_variants.tsv" >/dev/null

"$DOTMATCH_BIN" audit \
  --targets "$TMPDIR/targets.csv" \
  --k 1 \
  --audit-mode fast \
  --out-dir "$TMPDIR/audit_fast"

grep '^audit_mode	fast$' "$TMPDIR/audit_fast/audit_summary.tsv" >/dev/null
grep '^targets	4$' "$TMPDIR/audit_fast/audit_summary.tsv" >/dev/null
grep '^safe_at_k1	no$' "$TMPDIR/audit_fast/audit_summary.tsv" >/dev/null
grep '^risk_pairs_for_k1	3$' "$TMPDIR/audit_fast/audit_summary.tsv" >/dev/null
grep '^ambiguous_query_variants_k1	14$' "$TMPDIR/audit_fast/audit_summary.tsv" >/dev/null
grep '"audit_mode": "fast"' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
grep '"safe_at_k1": false' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
grep '"safe_at_k2": null' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
grep '"risk_pairs_for_k2": null' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
grep '^ACG	2$' "$TMPDIR/audit_fast/ambiguous_variants.tsv" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label sample1 \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --out "$TMPDIR/counts_report.tsv" \
  --report "$TMPDIR/report_rich.html" \
  --report-audit-dir "$TMPDIR/audit" \
  --report-unmatched "$TMPDIR/top_unmatched.tsv"

grep 'Library Audit' "$TMPDIR/report_rich.html" >/dev/null
grep 'Top Unmatched' "$TMPDIR/report_rich.html" >/dev/null
grep 'ambiguous_query_variants_k1' "$TMPDIR/report_rich.html" >/dev/null
grep 'near_known_target_above_k' "$TMPDIR/report_rich.html" >/dev/null

"$DOTMATCH_BIN" validate \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --oracle scan \
  --sample 10 | grep '"mismatches": 0' >/dev/null

if [ -x "$ROOT/build/dotmatch_edlib_validate" ]; then
  "$DOTMATCH_BIN" validate \
    --targets "$TMPDIR/targets.csv" \
    --reads "$TMPDIR/reads.fastq.gz" \
    --target-start 0 \
    --target-length 4 \
    --k 1 \
    --indel-window 1 \
    --oracle edlib \
    --threads 2 \
    --sample 10 > "$TMPDIR/edlib_validate.json"
  grep '"mismatches": 0' "$TMPDIR/edlib_validate.json" >/dev/null
  grep '"oracle_strategy": "bounded_edlib_candidates"' "$TMPDIR/edlib_validate.json" >/dev/null
  grep '"edlib_alignments": 9' "$TMPDIR/edlib_validate.json" >/dev/null

  cat > "$TMPDIR/sgrnaid_seq_targets.tsv" <<'TARGETS'
sgRNAID	Seq	gene
guide0	ACGT	G0
guide1	AGGT	G1
TARGETS
  "$DOTMATCH_BIN" validate \
    --targets "$TMPDIR/sgrnaid_seq_targets.tsv" \
    --reads "$TMPDIR/reads.fastq.gz" \
    --target-start 0 \
    --target-length 4 \
    --k 1 \
    --oracle edlib \
    --sample 1 > "$TMPDIR/edlib_sgrnaid_header.json"
  grep '"mismatches": 0' "$TMPDIR/edlib_sgrnaid_header.json" >/dev/null
  grep '"bounded_windows": 1' "$TMPDIR/edlib_sgrnaid_header.json" >/dev/null
  grep '"fallback_windows": 0' "$TMPDIR/edlib_sgrnaid_header.json" >/dev/null
else
  if "$DOTMATCH_BIN" validate \
    --targets "$TMPDIR/targets.csv" \
    --reads "$TMPDIR/reads.fastq.gz" \
    --target-start 0 \
    --target-length 4 \
    --k 1 \
    --oracle edlib 2>/dev/null; then
    echo "edlib oracle should require edlib-tools helper" >&2
    exit 1
  fi
fi
