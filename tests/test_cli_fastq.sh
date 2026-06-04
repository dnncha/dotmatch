#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
DOTMATCH_BIN=${DOTMATCH_BIN:-"$ROOT/dotmatch"}
TMPDIR="${TMPDIR:-/tmp}/dotmatch-cli-$$"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

EXPECTED_VERSION=$(python3 - <<'PY'
import re
from pathlib import Path
text = Path("pyproject.toml").read_text(encoding="utf-8")
match = re.search(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
if match is None:
    raise SystemExit("pyproject.toml version not found")
print(match.group(1))
PY
)
if [ "$("$DOTMATCH_BIN" --version)" != "dotmatch $EXPECTED_VERSION" ]; then
  echo "dotmatch --version must report pyproject version" >&2
  exit 1
fi
for help_arg in --help -h help; do
  "$DOTMATCH_BIN" "$help_arg" > "$TMPDIR/help-$help_arg.txt" 2> "$TMPDIR/help-$help_arg.err"
  if [ -s "$TMPDIR/help-$help_arg.err" ]; then
    echo "dotmatch $help_arg must write help to stdout without stderr" >&2
    exit 1
  fi
done
cp "$TMPDIR/help---help.txt" "$TMPDIR/help.txt"
grep "DotMatch $EXPECTED_VERSION" "$TMPDIR/help.txt" >/dev/null
grep "Counting and demultiplexing:" "$TMPDIR/help.txt" >/dev/null
grep "Assignment outcomes:" "$TMPDIR/help.txt" >/dev/null
grep "citation" "$TMPDIR/help.txt" >/dev/null
if grep "Packaging note:" "$TMPDIR/help.txt" >/dev/null; then
  echo "native help must not include packaging policy notes" >&2
  exit 1
fi
"$DOTMATCH_BIN" citation > "$TMPDIR/citation.txt"
grep "Software release: v$EXPECTED_VERSION" "$TMPDIR/citation.txt" >/dev/null
grep "O'Toole D. DotMatch: Streaming Exact One-Edit Barcode and Guide Assignment Without Exhaustive Scanning. Software release v$EXPECTED_VERSION." "$TMPDIR/citation.txt" >/dev/null
grep "Citation metadata: CITATION.cff" "$TMPDIR/citation.txt" >/dev/null
grep "DOI: 10.5281/zenodo.20541629" "$TMPDIR/citation.txt" >/dev/null
grep "DOI URL: https://doi.org/10.5281/zenodo.20541629" "$TMPDIR/citation.txt" >/dev/null

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

cat > "$TMPDIR/assign_reads.tsv" <<'ASSIGNREADS'
r0	ACGT
r1	ACGC
r2	TTTT
ASSIGNREADS

"$DOTMATCH_BIN" assign 1 "$TMPDIR/barcodes.tsv" "$TMPDIR/assign_reads.tsv" > "$TMPDIR/assign_radius.tsv"
grep '^assign	r0	ACGT	0	ACGT	0	ambiguous	3	1$' "$TMPDIR/assign_radius.tsv" >/dev/null
grep '^assign	r1	ACGC	0	ACGT	1	ambiguous	2	-1$' "$TMPDIR/assign_radius.tsv" >/dev/null

if "$DOTMATCH_BIN" assign 1x "$TMPDIR/barcodes.tsv" "$TMPDIR/assign_reads.tsv" > "$TMPDIR/assign_bad_k.tsv" 2> "$TMPDIR/assign_bad_k.err"; then
  echo "assign should reject malformed edit-distance radius values" >&2
  exit 1
fi
grep 'assign K barcodes.txt reads.txt' "$TMPDIR/assign_bad_k.err" >/dev/null
test ! -s "$TMPDIR/assign_bad_k.tsv"

"$DOTMATCH_BIN" assign 1 "$TMPDIR/barcodes.tsv" "$TMPDIR/assign_reads.tsv" --ambiguity-policy best > "$TMPDIR/assign_best.tsv"
grep '^assign	r0	ACGT	0	ACGT	0	unique	3	1$' "$TMPDIR/assign_best.tsv" >/dev/null

if "$DOTMATCH_BIN" leq 1x ACGT AGGT > "$TMPDIR/leq_bad_k.txt" 2> "$TMPDIR/leq_bad_k.err"; then
  echo "leq should reject malformed edit-distance radius values" >&2
  exit 1
fi
grep 'leq K SEQ1 SEQ2' "$TMPDIR/leq_bad_k.err" >/dev/null
test ! -s "$TMPDIR/leq_bad_k.txt"

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

cat > "$TMPDIR/expected_radius.tsv" <<'EXPECTEDRADIUS'
read_id	observed_barcode	target_index	target_id	target_seq	best_distance	second_best_distance	match_count	status
r0	ACGT	0	bc0	ACGT	0	1	3	ambiguous
r1	TTTG	3	bc3	TTTT	1	-1	1	unique
r2	GGGG	-1			-1	-1	0	none
r3		-1			-1	-1	0	invalid
EXPECTEDRADIUS

diff -u "$TMPDIR/expected_radius.tsv" "$TMPDIR/out.tsv"

cat > "$TMPDIR/duplicate_barcode_ids.tsv" <<'DUPBARCODES'
dup	ACGT
dup	TTTT
DUPBARCODES

if "$DOTMATCH_BIN" fastq-assign \
  --barcodes "$TMPDIR/duplicate_barcode_ids.tsv" \
  --reads "$TMPDIR/reads.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 1 \
  --out "$TMPDIR/duplicate_barcode_assignments.tsv" \
  2> "$TMPDIR/duplicate_barcode_ids.err"; then
  echo "fastq-assign should reject duplicate barcode IDs" >&2
  exit 1
fi
grep 'barcode IDs must be unique; duplicate ID: "dup"' "$TMPDIR/duplicate_barcode_ids.err" >/dev/null
test ! -e "$TMPDIR/duplicate_barcode_assignments.tsv"

cat > "$TMPDIR/empty_barcode_id.tsv" <<'EMPTYBARCODEID'
	ACGT
EMPTYBARCODEID

if "$DOTMATCH_BIN" fastq-assign \
  --barcodes "$TMPDIR/empty_barcode_id.tsv" \
  --reads "$TMPDIR/reads.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 1 \
  --out "$TMPDIR/empty_barcode_assignments.tsv" \
  2> "$TMPDIR/empty_barcode_id.err"; then
  echo "fastq-assign should reject empty barcode IDs" >&2
  exit 1
fi
grep 'record ID and sequence must be non-empty' "$TMPDIR/empty_barcode_id.err" >/dev/null
test ! -e "$TMPDIR/empty_barcode_assignments.tsv"

"$DOTMATCH_BIN" fastq-assign \
  --barcodes "$TMPDIR/barcodes.tsv" \
  --reads "$TMPDIR/reads.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 1 \
  --ambiguity-policy best \
  --out "$TMPDIR/out_best.tsv"

diff -u "$TMPDIR/expected.tsv" "$TMPDIR/out_best.tsv"

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

cat > "$TMPDIR/pair_left.tsv" <<'PAIRLEFT'
L0	ACGT
L1	TTTT
L2	AGGA
PAIRLEFT

cat > "$TMPDIR/pair_right.tsv" <<'PAIRRIGHT'
R0	GGAA
R1	CCCC
PAIRRIGHT

cat > "$TMPDIR/pair_reads.fastq" <<'PAIRFASTQ'
@p0
ACGTGGAA
+
IIIIIIII
@p1
TTTTCCCC
+
IIIIIIII
@p2
ACGCGGAA
+
IIIIIIII
@p3
GGGGGGAA
+
IIIIIIII
@p4
ACGTAAAA
+
IIIIIIII
@p5
AGGTGGAA
+
IIIIIIII
@p6
AC
+
II
PAIRFASTQ

"$DOTMATCH_BIN" pair-count \
  --left-targets "$TMPDIR/pair_left.tsv" \
  --right-targets "$TMPDIR/pair_right.tsv" \
  --reads "$TMPDIR/pair_reads.fastq" \
  --left-start 0 \
  --left-length 4 \
  --right-start 4 \
  --right-length 4 \
  --k 1 \
  --metric hamming \
  --out "$TMPDIR/pair_counts.tsv" \
  --summary "$TMPDIR/pair_summary.json" \
  --assignments "$TMPDIR/pair_assignments.tsv"

grep '^L0	R0	2$' "$TMPDIR/pair_counts.tsv" >/dev/null
grep '^L1	R1	1$' "$TMPDIR/pair_counts.tsv" >/dev/null
grep '"assigned_pairs": 3' "$TMPDIR/pair_summary.json" >/dev/null
grep '"pair_ambiguous": 1' "$TMPDIR/pair_summary.json" >/dev/null
grep '"left_unmatched": 1' "$TMPDIR/pair_summary.json" >/dev/null
grep '"right_unmatched": 1' "$TMPDIR/pair_summary.json" >/dev/null
grep '"invalid": 1' "$TMPDIR/pair_summary.json" >/dev/null
grep '^p5	AGGT	0	L0	ambiguous	1	GGAA	0	R0	unique	0	ambiguous$' "$TMPDIR/pair_assignments.tsv" >/dev/null
grep '^p6		-1		invalid	-1		-1		invalid	-1	invalid$' "$TMPDIR/pair_assignments.tsv" >/dev/null

cat > "$TMPDIR/pair_left_duplicate.tsv" <<'TARGETS'
Ldup	ACGT
Ldup	TTTT
TARGETS

if "$DOTMATCH_BIN" pair-count \
  --left-targets "$TMPDIR/pair_left_duplicate.tsv" \
  --right-targets "$TMPDIR/pair_right.tsv" \
  --reads "$TMPDIR/pair_reads.fastq" \
  --left-start 0 \
  --left-length 4 \
  --right-start 4 \
  --right-length 4 \
  --k 1 \
  --metric hamming \
  --out "$TMPDIR/pair_duplicate_left.tsv" 2> "$TMPDIR/pair_duplicate_left.err"; then
  echo "pair-count accepted duplicate left target IDs" >&2
  exit 1
fi
grep 'left target IDs must be unique; duplicate ID: "Ldup"' "$TMPDIR/pair_duplicate_left.err" >/dev/null
test ! -s "$TMPDIR/pair_duplicate_left.tsv"

cat > "$TMPDIR/pair_right_duplicate.tsv" <<'TARGETS'
Rdup	GGAA
Rdup	CCCC
TARGETS

if "$DOTMATCH_BIN" pair-count \
  --left-targets "$TMPDIR/pair_left.tsv" \
  --right-targets "$TMPDIR/pair_right_duplicate.tsv" \
  --reads "$TMPDIR/pair_reads.fastq" \
  --left-start 0 \
  --left-length 4 \
  --right-start 4 \
  --right-length 4 \
  --k 1 \
  --metric hamming \
  --out "$TMPDIR/pair_duplicate_right.tsv" 2> "$TMPDIR/pair_duplicate_right.err"; then
  echo "pair-count accepted duplicate right target IDs" >&2
  exit 1
fi
grep 'right target IDs must be unique; duplicate ID: "Rdup"' "$TMPDIR/pair_duplicate_right.err" >/dev/null
test ! -s "$TMPDIR/pair_duplicate_right.tsv"

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

cat > "$TMPDIR/barcode_header.tsv" <<'BARCODEHEADER'
barcode_id	barcode_seq
s1	ACGT
s2	TTTT
BARCODEHEADER

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/barcode_header.tsv" \
  --reads "$TMPDIR/reads.fastq" \
  --target-start 0 \
  --target-length 4 \
  --k 0 \
  --metric hamming \
  --out "$TMPDIR/barcode_header_counts.tsv" \
  --summary "$TMPDIR/barcode_header_summary.json"

grep '^s1	ACGT		0	1	0	0	0	0	1$' "$TMPDIR/barcode_header_counts.tsv" >/dev/null
grep '"n_targets": 2' "$TMPDIR/barcode_header_summary.json" >/dev/null

cat > "$TMPDIR/k3_targets.tsv" <<'K3TARGETS'
only	ACGT
K3TARGETS

cat > "$TMPDIR/k3_reads.fastq" <<'K3FASTQ'
@exact
ACGT
+
IIII
@two_mismatch
AGGA
+
IIII
@three_mismatch
TGGA
+
IIII
@four_mismatch
TGCA
+
IIII
K3FASTQ

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/k3_targets.tsv" \
  --reads "$TMPDIR/k3_reads.fastq" \
  --sample-label k3 \
  --target-start 0 \
  --target-length 4 \
  --k 2 \
  --metric hamming \
  --ambiguity-policy best \
  --out "$TMPDIR/k2_hamming_counts.tsv" \
  --summary "$TMPDIR/k2_hamming_summary.json"

grep '^only	ACGT		0	1	0	0	0	1	2$' "$TMPDIR/k2_hamming_counts.tsv" >/dev/null
grep '"k": 2' "$TMPDIR/k2_hamming_summary.json" >/dev/null
grep '"metric": "hamming"' "$TMPDIR/k2_hamming_summary.json" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/k3_targets.tsv" \
  --reads "$TMPDIR/k3_reads.fastq" \
  --sample-label k3 \
  --target-start 0 \
  --target-length 4 \
  --k 3 \
  --metric hamming \
  --ambiguity-policy best \
  --out "$TMPDIR/k3_hamming_counts.tsv" \
  --summary "$TMPDIR/k3_hamming_summary.json"

grep '^only	ACGT		0	1	0	0	0	2	3$' "$TMPDIR/k3_hamming_counts.tsv" >/dev/null
grep '"k": 3' "$TMPDIR/k3_hamming_summary.json" >/dev/null
grep '"metric": "hamming"' "$TMPDIR/k3_hamming_summary.json" >/dev/null

if "$DOTMATCH_BIN" count \
  --targets "$TMPDIR/k3_targets.tsv" \
  --reads "$TMPDIR/k3_reads.fastq" \
  --sample-label k3 \
  --target-start 0 \
  --target-length 4 \
  --k 3 \
  --metric levenshtein \
  --out "$TMPDIR/bad_levenshtein_k3.tsv" 2>/dev/null; then
  echo "Levenshtein k=3 should fail for count" >&2
  exit 1
fi

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
  --ambiguity-policy best \
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
grep '"alphabet_policy": "literal-byte; A/C/G/T/N/IUPAC symbols are ordinary byte symbols; no wildcard expansion"' "$TMPDIR/demux_summary.json" >/dev/null
grep '"ambiguous": 1' "$TMPDIR/demux_summary.json" >/dev/null
grep '"unmatched": 1' "$TMPDIR/demux_summary.json" >/dev/null
grep '"invalid": 1' "$TMPDIR/demux_summary.json" >/dev/null
grep '^d2	AGGA	1	bc1	AGGT	1	-1	2	ambiguous$' "$TMPDIR/demux_assignments.tsv" >/dev/null

mkdir "$TMPDIR/demux_radius_default"
"$DOTMATCH_BIN" demux \
  --barcodes "$TMPDIR/barcodes.tsv" \
  --reads "$TMPDIR/demux_reads.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 1 \
  --metric hamming \
  --out-dir "$TMPDIR/demux_radius_default" \
  --summary "$TMPDIR/demux_radius_default_summary.json" \
  --assignments "$TMPDIR/demux_radius_default_assignments.tsv" \
  --ambiguous-out "$TMPDIR/demux_radius_default_ambiguous.fastq"

grep '"ambiguity_policy": "radius"' "$TMPDIR/demux_radius_default_summary.json" >/dev/null
grep '^@d0$' "$TMPDIR/demux_radius_default_ambiguous.fastq" >/dev/null
grep '^d0	ACGT	0	bc0	ACGT	0	1	3	ambiguous$' "$TMPDIR/demux_radius_default_assignments.tsv" >/dev/null

cat > "$TMPDIR/demux_colliding_ids.tsv" <<'DEMUXCOLLIDE'
a/b	ACGT
a:b	TTTT
DEMUXCOLLIDE

if "$DOTMATCH_BIN" demux \
  --barcodes "$TMPDIR/demux_colliding_ids.tsv" \
  --reads "$TMPDIR/demux_reads.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 0 \
  --metric hamming \
  --out-dir "$TMPDIR/demux_colliding_ids" \
  2> "$TMPDIR/demux_colliding_ids.err"; then
  echo "demux should reject barcode IDs that collide after filename sanitization" >&2
  exit 1
fi
grep 'barcode IDs produce the same output filename after sanitization' "$TMPDIR/demux_colliding_ids.err" >/dev/null
test ! -d "$TMPDIR/demux_colliding_ids"

cat > "$TMPDIR/demux_duplicate_ids.tsv" <<'DEMUXDUP'
dup	ACGT
dup	TTTT
DEMUXDUP

if "$DOTMATCH_BIN" demux \
  --barcodes "$TMPDIR/demux_duplicate_ids.tsv" \
  --reads "$TMPDIR/demux_reads.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 0 \
  --metric hamming \
  --out-dir "$TMPDIR/demux_duplicate_ids" \
  2> "$TMPDIR/demux_duplicate_ids.err"; then
  echo "demux should reject duplicate barcode IDs" >&2
  exit 1
fi
grep 'barcode IDs must be unique; duplicate ID: "dup"' "$TMPDIR/demux_duplicate_ids.err" >/dev/null
test ! -d "$TMPDIR/demux_duplicate_ids"

cat > "$TMPDIR/demux_empty_id.tsv" <<'DEMUXEMPTY'
barcode_id	barcode_seq
	ACGT
DEMUXEMPTY

if "$DOTMATCH_BIN" demux \
  --barcodes "$TMPDIR/demux_empty_id.tsv" \
  --reads "$TMPDIR/demux_reads.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 0 \
  --metric hamming \
  --out-dir "$TMPDIR/demux_empty_id" \
  2> "$TMPDIR/demux_empty_id.err"; then
  echo "demux should reject empty barcode IDs" >&2
  exit 1
fi
grep 'target ID and sequence must be non-empty' "$TMPDIR/demux_empty_id.err" >/dev/null
test ! -d "$TMPDIR/demux_empty_id"

cat > "$TMPDIR/demux_quality.fastq" <<'DEMUXQUAL'
@dq_exact
ACGTAAAA
+
IIIIIIII
@dq_low
ACGCAAAA
+
III!IIII
@dq_high
ACGCAAAA
+
IIIIIIII
DEMUXQUAL

cat > "$TMPDIR/demux_quality_barcodes.tsv" <<'DEMUXQUALBC'
bc0	ACGT
DEMUXQUALBC

mkdir "$TMPDIR/demux_quality"
"$DOTMATCH_BIN" demux \
  --barcodes "$TMPDIR/demux_quality_barcodes.tsv" \
  --reads "$TMPDIR/demux_quality.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 1 \
  --metric hamming \
  --max-correction-qual 20 \
  --out-dir "$TMPDIR/demux_quality" \
  --summary "$TMPDIR/demux_quality_summary.json" \
  --assignments "$TMPDIR/demux_quality_assignments.tsv" \
  --unmatched-out "$TMPDIR/demux_quality_unmatched.fastq"

grep '^@dq_exact$' "$TMPDIR/demux_quality/bc0.fastq" >/dev/null
grep '^@dq_low$' "$TMPDIR/demux_quality/bc0.fastq" >/dev/null
grep '^@dq_high$' "$TMPDIR/demux_quality_unmatched.fastq" >/dev/null
grep '"max_correction_qual": 20' "$TMPDIR/demux_quality_summary.json" >/dev/null
grep '"unmatched": 1' "$TMPDIR/demux_quality_summary.json" >/dev/null
grep '^dq_high	ACGC	-1			-1	-1	0	none$' "$TMPDIR/demux_quality_assignments.tsv" >/dev/null

cat > "$TMPDIR/demux_k2.fastq" <<'DEMUXK2'
@dk_exact
ACGTAAAA
+
IIIIIIII
@dk_two_sub
AATTAAAA
+
IIIIIIII
@dk_none
TTTTAAAA
+
IIIIIIII
DEMUXK2

mkdir "$TMPDIR/demux_k2"
"$DOTMATCH_BIN" demux \
  --barcodes "$TMPDIR/demux_quality_barcodes.tsv" \
  --reads "$TMPDIR/demux_k2.fastq" \
  --barcode-start 0 \
  --barcode-length 4 \
  --k 2 \
  --metric levenshtein \
  --out-dir "$TMPDIR/demux_k2" \
  --summary "$TMPDIR/demux_k2_summary.json" \
  --assignments "$TMPDIR/demux_k2_assignments.tsv" \
  --unmatched-out "$TMPDIR/demux_k2_unmatched.fastq"

grep '^@dk_exact$' "$TMPDIR/demux_k2/bc0.fastq" >/dev/null
grep '^@dk_two_sub$' "$TMPDIR/demux_k2/bc0.fastq" >/dev/null
grep '^@dk_none$' "$TMPDIR/demux_k2_unmatched.fastq" >/dev/null
grep '"k": 2' "$TMPDIR/demux_k2_summary.json" >/dev/null
grep '"assigned_corrected": 1' "$TMPDIR/demux_k2_summary.json" >/dev/null
grep '^dk_two_sub	AATT	0	bc0	ACGT	2	-1	1	unique$' "$TMPDIR/demux_k2_assignments.tsv" >/dev/null

cat > "$TMPDIR/variable_barcodes.tsv" <<'VARBC'
long	ACGA
short	TTT
prefix_short	GG
prefix_long	GGGG
VARBC

cat > "$TMPDIR/variable_reads.fastq" <<'VARFASTQ'
@v0
ACGAAAAA
+
IIIIIIII
@v1
TTTAAAAA
+
IIIIIIII
@v2
GGGGAAAA
+
IIIIIIII
@v3
CCCCAAAA
+
IIIIIIII
VARFASTQ

mkdir "$TMPDIR/demux_variable"
"$DOTMATCH_BIN" demux \
  --barcodes "$TMPDIR/variable_barcodes.tsv" \
  --reads "$TMPDIR/variable_reads.fastq" \
  --barcode-start 0 \
  --barcode-length auto \
  --k 0 \
  --metric hamming \
  --out-dir "$TMPDIR/demux_variable" \
  --summary "$TMPDIR/demux_variable_summary.json" \
  --assignments "$TMPDIR/demux_variable_assignments.tsv" \
  --ambiguous-out "$TMPDIR/demux_variable_ambiguous.fastq" \
  --unmatched-out "$TMPDIR/demux_variable_unmatched.fastq"

grep '^@v0$' "$TMPDIR/demux_variable/long.fastq" >/dev/null
grep '^@v1$' "$TMPDIR/demux_variable/short.fastq" >/dev/null
grep '^@v2$' "$TMPDIR/demux_variable_ambiguous.fastq" >/dev/null
grep '^@v3$' "$TMPDIR/demux_variable_unmatched.fastq" >/dev/null
grep '"barcode_length_mode": "auto"' "$TMPDIR/demux_variable_summary.json" >/dev/null
grep '"assigned_unique": 2' "$TMPDIR/demux_variable_summary.json" >/dev/null
grep '"ambiguous": 1' "$TMPDIR/demux_variable_summary.json" >/dev/null
grep '^v2	GG	2	prefix_short	GG	0	-1	2	ambiguous$' "$TMPDIR/demux_variable_assignments.tsv" >/dev/null

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

if "$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_run" \
  --sample-sheet "$TMPDIR/bcl_run/SampleSheet.csv" \
  --out-dir "$TMPDIR/bcl_bad_lane_out" \
  --lanes 2 \
  --summary "$TMPDIR/bcl_bad_lane_summary.json" 2>/dev/null; then
  echo "bcl-demux should reject unsupported non-lane-1 requests" >&2
  exit 1
fi

cp -R "$TMPDIR/bcl_run" "$TMPDIR/bcl_bad_read_number_run"
python3 - "$TMPDIR/bcl_bad_read_number_run/RunInfo.xml" <<'PY'
import sys
from pathlib import Path

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
path.write_text(text.replace('Read Number="1"', 'Read Number="not_a_read"', 1), encoding="utf-8")
PY

if "$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_bad_read_number_run" \
  --sample-sheet "$TMPDIR/bcl_bad_read_number_run/SampleSheet.csv" \
  --out-dir "$TMPDIR/bcl_bad_read_number_out" \
  --barcode-mismatches 0 \
  2> "$TMPDIR/bcl_bad_read_number.err"; then
  echo "bcl-demux should reject malformed RunInfo read numbers" >&2
  exit 1
fi
grep 'failed to parse RunInfo.xml' "$TMPDIR/bcl_bad_read_number.err" >/dev/null
test ! -d "$TMPDIR/bcl_bad_read_number_out"

cp -R "$TMPDIR/bcl_run" "$TMPDIR/bcl_bad_indexed_flag_run"
python3 - "$TMPDIR/bcl_bad_indexed_flag_run/RunInfo.xml" <<'PY'
import sys
from pathlib import Path

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
path.write_text(text.replace('IsIndexedRead="Y"', 'IsIndexedRead="maybe"', 1), encoding="utf-8")
PY

if "$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_bad_indexed_flag_run" \
  --sample-sheet "$TMPDIR/bcl_bad_indexed_flag_run/SampleSheet.csv" \
  --out-dir "$TMPDIR/bcl_bad_indexed_flag_out" \
  --barcode-mismatches 0 \
  2> "$TMPDIR/bcl_bad_indexed_flag.err"; then
  echo "bcl-demux should reject malformed RunInfo IsIndexedRead values" >&2
  exit 1
fi
grep 'failed to parse RunInfo.xml' "$TMPDIR/bcl_bad_indexed_flag.err" >/dev/null
test ! -d "$TMPDIR/bcl_bad_indexed_flag_out"

if "$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_run" \
  --sample-sheet "$TMPDIR/bcl_run/SampleSheet.csv" \
  --out-dir "$TMPDIR/bcl_bad_mismatch_out" \
  --barcode-mismatches x \
  2> "$TMPDIR/bcl_bad_mismatch.err"; then
  echo "bcl-demux should reject nonnumeric --barcode-mismatches values" >&2
  exit 1
fi
grep 'bcl-demux --run-folder RUN --sample-sheet SampleSheet.csv' "$TMPDIR/bcl_bad_mismatch.err" >/dev/null
test ! -d "$TMPDIR/bcl_bad_mismatch_out"

if "$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_run" \
  --sample-sheet "$TMPDIR/bcl_run/SampleSheet.csv" \
  --out-dir "$TMPDIR/bcl_bad_dual_mismatch_out" \
  --barcode-mismatches 1,x \
  2> "$TMPDIR/bcl_bad_dual_mismatch.err"; then
  echo "bcl-demux should reject malformed dual-index --barcode-mismatches values" >&2
  exit 1
fi
grep 'bcl-demux --run-folder RUN --sample-sheet SampleSheet.csv' "$TMPDIR/bcl_bad_dual_mismatch.err" >/dev/null
test ! -d "$TMPDIR/bcl_bad_dual_mismatch_out"

cat > "$TMPDIR/bcl_run/SampleSheet.empty-id.csv" <<'SHEET'
[Header]
IEMFileVersion,4
[Data]
Sample_ID,Sample_Name,index
,Missing ID,ACGT
SHEET

if "$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_run" \
  --sample-sheet "$TMPDIR/bcl_run/SampleSheet.empty-id.csv" \
  --out-dir "$TMPDIR/bcl_empty_id_out" \
  --barcode-mismatches 0 \
  2> "$TMPDIR/bcl_empty_id.err"; then
  echo "bcl-demux should reject empty Sample_ID values" >&2
  exit 1
fi
grep 'BCL sample sheet Sample_ID and index must be non-empty' "$TMPDIR/bcl_empty_id.err" >/dev/null
test ! -d "$TMPDIR/bcl_empty_id_out"

cat > "$TMPDIR/bcl_run/SampleSheet.empty-index.csv" <<'SHEET'
[Header]
IEMFileVersion,4
[Data]
Sample_ID,Sample_Name,index
s1,Missing Index,
SHEET

if "$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_run" \
  --sample-sheet "$TMPDIR/bcl_run/SampleSheet.empty-index.csv" \
  --out-dir "$TMPDIR/bcl_empty_index_out" \
  --barcode-mismatches 0 \
  2> "$TMPDIR/bcl_empty_index.err"; then
  echo "bcl-demux should reject empty index values" >&2
  exit 1
fi
grep 'BCL sample sheet Sample_ID and index must be non-empty' "$TMPDIR/bcl_empty_index.err" >/dev/null
test ! -d "$TMPDIR/bcl_empty_index_out"

cat > "$TMPDIR/bcl_run/SampleSheet.bad-lane.csv" <<'SHEET'
[Header]
IEMFileVersion,4
[Data]
Sample_ID,Sample_Name,index,Lane
s1,Bad Lane,ACGT,not_a_lane
SHEET

if "$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_run" \
  --sample-sheet "$TMPDIR/bcl_run/SampleSheet.bad-lane.csv" \
  --out-dir "$TMPDIR/bcl_bad_sample_lane_out" \
  --barcode-mismatches 0 \
  2> "$TMPDIR/bcl_bad_sample_lane.err"; then
  echo "bcl-demux should reject nonnumeric sample-sheet Lane values" >&2
  exit 1
fi
grep 'BCL sample sheet Lane must be a positive integer' "$TMPDIR/bcl_bad_sample_lane.err" >/dev/null
test ! -d "$TMPDIR/bcl_bad_sample_lane_out"

cat > "$TMPDIR/bcl_run/SampleSheet.unsupported-lane.csv" <<'SHEET'
[Header]
IEMFileVersion,4
[Data]
Sample_ID,Sample_Name,index,Lane
s1,Unsupported Lane,ACGT,2
SHEET

if "$DOTMATCH_BIN" bcl-demux \
  --run-folder "$TMPDIR/bcl_run" \
  --sample-sheet "$TMPDIR/bcl_run/SampleSheet.unsupported-lane.csv" \
  --out-dir "$TMPDIR/bcl_unsupported_sample_lane_out" \
  --barcode-mismatches 0 \
  2> "$TMPDIR/bcl_unsupported_sample_lane.err"; then
  echo "bcl-demux should reject unsupported sample-sheet Lane values" >&2
  exit 1
fi
grep 'classic BCL demux currently supports sample-sheet Lane 1 only' "$TMPDIR/bcl_unsupported_sample_lane.err" >/dev/null
test ! -d "$TMPDIR/bcl_unsupported_sample_lane_out"

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
  --ambiguity-policy best \
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
grep '"alphabet_policy": "literal-byte; A/C/G/T/N/IUPAC symbols are ordinary byte symbols; no wildcard expansion"' "$TMPDIR/summary.json" >/dev/null
grep '"metric": "levenshtein"' "$TMPDIR/summary.json" >/dev/null
grep '"library_covered_targets": 2' "$TMPDIR/summary.json" >/dev/null
grep '^sample1	r2	GGGG	-1			-1	-1	0	none	none$' "$TMPDIR/unmatched.tsv" >/dev/null
grep '^sample1	' "$TMPDIR/sample_qc.tsv" | grep '	4	3	2	1	1	1	0	0	' >/dev/null
awk -F '\t' 'NR == 2 && ($14 != "0.66666667" || $15 != "0.33333333" || $16 != "0.33333333" || $18 != "0.33333333") { exit 1 }' "$TMPDIR/sample_qc.tsv"
awk -F '\t' 'NR == 2 && ($21 < 0 || $21 > 1) { exit 1 }' "$TMPDIR/sample_qc.tsv"
grep '^sample1	bc0	G0	ACGT	1	0	0	0	0	1	1$' "$TMPDIR/target_counts.long.tsv" >/dev/null
grep '<title>DotMatch Report</title>' "$TMPDIR/report.html" >/dev/null
grep 'sample1' "$TMPDIR/report.html" >/dev/null
grep 'Run Status' "$TMPDIR/report.html" >/dev/null
grep 'Target Assignment QC' "$TMPDIR/report.html" >/dev/null
grep 'Inputs and Configuration' "$TMPDIR/report.html" >/dev/null
grep 'Assignment rate' "$TMPDIR/report.html" >/dev/null
grep 'Valid windows' "$TMPDIR/report.html" >/dev/null
grep '66.67%' "$TMPDIR/report.html" >/dev/null
grep 'Library coverage' "$TMPDIR/report.html" >/dev/null
REPORT_MODE=$(python3 - "$TMPDIR/report.html" <<'PY'
import os
import stat
import sys
print(oct(stat.S_IMODE(os.stat(sys.argv[1]).st_mode)))
PY
)
if [ "$REPORT_MODE" != "0o600" ]; then
  echo "HTML report mode should be 0600, got $REPORT_MODE" >&2
  exit 1
fi

if "$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label bad_numeric \
  --target-start -1 \
  --target-length 4 \
  --k 1 \
  --out "$TMPDIR/negative_target_start.tsv" 2> "$TMPDIR/negative_target_start.err"; then
  echo "count should reject negative unsigned numeric arguments" >&2
  exit 1
fi
grep 'count --targets targets.tsv|targets.csv --reads reads.fastq' "$TMPDIR/negative_target_start.err" >/dev/null
test ! -e "$TMPDIR/negative_target_start.tsv"

if "$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label bad_numeric \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --auto-offset 1 \
  --offset-min-fraction NaN \
  --out "$TMPDIR/nan_offset_fraction.tsv" 2> "$TMPDIR/nan_offset_fraction.err"; then
  echo "count should reject non-finite numeric arguments" >&2
  exit 1
fi
grep 'count --targets targets.tsv|targets.csv --reads reads.fastq' "$TMPDIR/nan_offset_fraction.err" >/dev/null
test ! -e "$TMPDIR/nan_offset_fraction.tsv"

if "$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label duplicate,duplicate \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --out "$TMPDIR/duplicate_sample_counts.tsv" \
  2> "$TMPDIR/duplicate_sample_labels.err"; then
  echo "count should reject duplicate sample labels" >&2
  exit 1
fi
grep 'duplicate sample label: "duplicate"' "$TMPDIR/duplicate_sample_labels.err" >/dev/null
test ! -e "$TMPDIR/duplicate_sample_counts.tsv"

cat > "$TMPDIR/duplicate_target_ids.tsv" <<'DUPTARGETS'
target_id	target_seq	gene
dup	ACGT	G0
dup	TTTT	G1
DUPTARGETS

if "$DOTMATCH_BIN" count \
  --targets "$TMPDIR/duplicate_target_ids.tsv" \
  --reads "$TMPDIR/reads.fastq" \
  --sample-label sample1 \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --out "$TMPDIR/duplicate_target_counts.tsv" \
  2> "$TMPDIR/duplicate_target_ids.err"; then
  echo "count should reject duplicate target IDs" >&2
  exit 1
fi
grep 'target IDs must be unique; duplicate ID: "dup"' "$TMPDIR/duplicate_target_ids.err" >/dev/null
test ! -e "$TMPDIR/duplicate_target_counts.tsv"

cat > "$TMPDIR/empty_target_id.tsv" <<'EMPTYTARGET'
target_id	target_seq	gene
	ACGT	G0
EMPTYTARGET

if "$DOTMATCH_BIN" count \
  --targets "$TMPDIR/empty_target_id.tsv" \
  --reads "$TMPDIR/reads.fastq" \
  --sample-label sample1 \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --out "$TMPDIR/empty_target_counts.tsv" \
  2> "$TMPDIR/empty_target_id.err"; then
  echo "count should reject empty target IDs" >&2
  exit 1
fi
grep 'target ID and sequence must be non-empty' "$TMPDIR/empty_target_id.err" >/dev/null
test ! -e "$TMPDIR/empty_target_counts.tsv"

cat > "$TMPDIR/quality_reads.fastq" <<'QUALFASTQ'
@q_exact
ACGTAAAA
+
IIIIIIII
@q_low
ACGCAAAA
+
III!IIII
@q_high
ACGCAAAA
+
IIIIIIII
QUALFASTQ

cat > "$TMPDIR/quality_targets.tsv" <<'QUALTARGETS'
bc0	ACGT	G0
QUALTARGETS

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/quality_targets.tsv" \
  --reads "$TMPDIR/quality_reads.fastq" \
  --sample-label quality \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --max-correction-qual 20 \
  --out "$TMPDIR/quality_counts.tsv" \
  --summary "$TMPDIR/quality_summary.json" \
  --unmatched-out "$TMPDIR/quality_unmatched.tsv"

grep '^bc0	ACGT	G0	0	1	1	0	0	0	2$' "$TMPDIR/quality_counts.tsv" >/dev/null
grep '"max_correction_qual": 20' "$TMPDIR/quality_summary.json" >/dev/null
grep '"unmatched": 1' "$TMPDIR/quality_summary.json" >/dev/null
grep '^quality	q_high	ACGC	-1			-1	-1	0	none	quality_rejected$' "$TMPDIR/quality_unmatched.tsv" >/dev/null

cat > "$TMPDIR/k2_reads.fastq" <<'K2FASTQ'
@k2_exact
ACGTAAAA
+
IIIIIIII
@k2_two_sub
AATTAAAA
+
IIIIIIII
@k2_none
TTTTAAAA
+
IIIIIIII
K2FASTQ

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/quality_targets.tsv" \
  --reads "$TMPDIR/k2_reads.fastq" \
  --sample-label k2 \
  --target-start 0 \
  --target-length 4 \
  --k 2 \
  --metric levenshtein \
  --out "$TMPDIR/k2_counts.tsv" \
  --summary "$TMPDIR/k2_summary.json" \
  --unmatched-out "$TMPDIR/k2_unmatched.tsv"

grep '^bc0	ACGT	G0	0	1	0	0	0	1	2$' "$TMPDIR/k2_counts.tsv" >/dev/null
grep '"k": 2' "$TMPDIR/k2_summary.json" >/dev/null
grep '"assigned_corrected": 1' "$TMPDIR/k2_summary.json" >/dev/null
grep '^k2	k2_none	TTTT	-1			-1	-1	0	none	none$' "$TMPDIR/k2_unmatched.tsv" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label sample_lev \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric levenshtein \
  --ambiguity-policy best \
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
  --ambiguity-policy best \
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
  --ambiguity-policy best \
  --format mageck \
  --out "$TMPDIR/counts_hamming.tsv" \
  --summary "$TMPDIR/summary_hamming.json"

grep '^bc0	G0	1$' "$TMPDIR/counts_hamming.tsv" >/dev/null
grep '^bc3	G3	1$' "$TMPDIR/counts_hamming.tsv" >/dev/null
grep '"count_engine": "hamming_lookup_direct_single_offset"' "$TMPDIR/summary_hamming.json" >/dev/null

cat > "$TMPDIR/gc_library.tsv" <<'GCLIB'
guide	bases	gene
g_exact	ACGT	GENE_ESS
g_control	TTTT	CTRL_SAFE
g_other	GGCC	GENE_OTHER
GCLIB

cat > "$TMPDIR/gc_essential.tsv" <<'GCESS'
GENE_ESS
GCESS

cat > "$TMPDIR/gc_nonessential.tsv" <<'GCNON'
GENE_OTHER
GCNON

cat > "$TMPDIR/gc_control.tsv" <<'GCCTRL'
g_control
GCCTRL

cat > "$TMPDIR/gc_sample.fastq" <<'GCFASTQ'
@gc0
NNACGTAA
+
IIIIIIII
@gc1
NNTTTTAA
+
IIIIIIII
@gc2
NNACGAAA
+
IIIIIIII
@gc3
NNGGCCAA
+
IIIIIIII
@gc4
NNAAAAAA
+
IIIIIIII
GCFASTQ

"$DOTMATCH_BIN" guide-counter count \
  --input "$TMPDIR/gc_sample.fastq" \
  --library "$TMPDIR/gc_library.tsv" \
  --essential-genes "$TMPDIR/gc_essential.tsv" \
  --nonessential-genes "$TMPDIR/gc_nonessential.tsv" \
  --control-guides "$TMPDIR/gc_control.tsv" \
  --control-pattern ctrl \
  --offset-sample-size 5 \
  --offset-min-fraction 0.1 \
  --output "$TMPDIR/gc_out"

cat > "$TMPDIR/gc_expected_counts.tsv" <<'GCCOUNTS'
guide	gene	gc_sample
g_exact	GENE_ESS	2
g_control	CTRL_SAFE	1
g_other	GENE_OTHER	1
GCCOUNTS
diff -u "$TMPDIR/gc_expected_counts.tsv" "$TMPDIR/gc_out.counts.txt"
grep '^guide	gene	guide_type	gc_sample$' "$TMPDIR/gc_out.extended-counts.txt" >/dev/null
grep '^g_exact	GENE_ESS	Essential	2$' "$TMPDIR/gc_out.extended-counts.txt" >/dev/null
grep '^g_control	CTRL_SAFE	Control	1$' "$TMPDIR/gc_out.extended-counts.txt" >/dev/null
grep '^g_other	GENE_OTHER	Nonessential	1$' "$TMPDIR/gc_out.extended-counts.txt" >/dev/null
grep '^file	label	total_guides	total_reads	mapped_reads	frac_mapped	mean_reads_per_guide	mean_reads_essential	mean_reads_nonessential	mean_reads_control	mean_reads_other	zero_read_guides$' "$TMPDIR/gc_out.stats.txt" >/dev/null
grep "^$TMPDIR/gc_sample.fastq	gc_sample	3	5	4	0.8000	1.33	2.00	1.00	1.00	0.00	0$" "$TMPDIR/gc_out.stats.txt" >/dev/null

if "$DOTMATCH_BIN" guide-counter count \
  --input "$TMPDIR/gc_sample.fastq" "$TMPDIR/gc_sample.fastq" \
  --samples duplicate duplicate \
  --library "$TMPDIR/gc_library.tsv" \
  --output "$TMPDIR/gc_duplicate" \
  2> "$TMPDIR/gc_duplicate.err"; then
  echo "guide-counter compatibility mode should reject duplicate sample labels" >&2
  exit 1
fi
grep 'duplicate sample label: "duplicate"' "$TMPDIR/gc_duplicate.err" >/dev/null
test ! -e "$TMPDIR/gc_duplicate.counts.txt"

cat > "$TMPDIR/gc_duplicate_guides.tsv" <<'GCDUPGUIDES'
guide	bases	gene
gdup	ACGT	GENE_ESS
gdup	TTTT	CTRL_SAFE
GCDUPGUIDES

if "$DOTMATCH_BIN" guide-counter count \
  --input "$TMPDIR/gc_sample.fastq" \
  --library "$TMPDIR/gc_duplicate_guides.tsv" \
  --output "$TMPDIR/gc_duplicate_guides" \
  2> "$TMPDIR/gc_duplicate_guides.err"; then
  echo "guide-counter compatibility mode should reject duplicate guide IDs" >&2
  exit 1
fi
grep 'guide IDs must be unique; duplicate ID: "gdup"' "$TMPDIR/gc_duplicate_guides.err" >/dev/null
test ! -e "$TMPDIR/gc_duplicate_guides.counts.txt"

cat > "$TMPDIR/gc_empty_guide.tsv" <<'GCEMPTYGUIDE'
guide	bases	gene
	ACGT	GENE_ESS
GCEMPTYGUIDE

if "$DOTMATCH_BIN" guide-counter count \
  --input "$TMPDIR/gc_sample.fastq" \
  --library "$TMPDIR/gc_empty_guide.tsv" \
  --output "$TMPDIR/gc_empty_guide" \
  2> "$TMPDIR/gc_empty_guide.err"; then
  echo "guide-counter compatibility mode should reject empty guide IDs" >&2
  exit 1
fi
grep 'target ID and sequence must be non-empty' "$TMPDIR/gc_empty_guide.err" >/dev/null
test ! -e "$TMPDIR/gc_empty_guide.counts.txt"

"$DOTMATCH_BIN" guide-counter-count \
  --input "$TMPDIR/gc_sample.fastq" \
  --samples explicit_sample \
  --library "$TMPDIR/gc_library.tsv" \
  --exact-match \
  --offset-sample-size 5 \
  --offset-min-fraction 0.1 \
  --output "$TMPDIR/gc_exact"

grep '^guide	gene	explicit_sample$' "$TMPDIR/gc_exact.counts.txt" >/dev/null
grep '^g_exact	GENE_ESS	1$' "$TMPDIR/gc_exact.counts.txt" >/dev/null
grep '^g_control	CTRL_SAFE	1$' "$TMPDIR/gc_exact.counts.txt" >/dev/null
grep '^g_other	GENE_OTHER	1$' "$TMPDIR/gc_exact.counts.txt" >/dev/null
grep "^$TMPDIR/gc_sample.fastq	explicit_sample	3	5	3	0.6000	1.00" "$TMPDIR/gc_exact.stats.txt" >/dev/null

"$DOTMATCH_BIN" guide-count \
  --input "$TMPDIR/gc_sample.fastq" \
  --samples short_alias \
  --library "$TMPDIR/gc_library.tsv" \
  --exact-match \
  --offset-sample-size 5 \
  --offset-min-fraction 0.1 \
  --output "$TMPDIR/gc_short_alias"

grep '^guide	gene	short_alias$' "$TMPDIR/gc_short_alias.counts.txt" >/dev/null
grep '^g_exact	GENE_ESS	1$' "$TMPDIR/gc_short_alias.counts.txt" >/dev/null
grep '^g_control	CTRL_SAFE	1$' "$TMPDIR/gc_short_alias.counts.txt" >/dev/null
grep '^g_other	GENE_OTHER	1$' "$TMPDIR/gc_short_alias.counts.txt" >/dev/null

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/targets.csv" \
  --reads "$TMPDIR/reads.fastq.gz" \
  --sample-label sample_hamming \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --metric hamming \
  --ambiguity-policy best \
  --format mageck \
  --threads 2 \
  --out "$TMPDIR/counts_hamming_threads.tsv" \
  --summary "$TMPDIR/summary_hamming_threads.json"

diff -u "$TMPDIR/counts_hamming.tsv" "$TMPDIR/counts_hamming_threads.tsv"
grep '"read_threads": 2' "$TMPDIR/summary_hamming_threads.json" >/dev/null

python3 - "$TMPDIR/sparse_targets.tsv" "$TMPDIR/sparse_reads.fastq" <<'PY'
import random
import sys

targets_path, reads_path = sys.argv[1:3]
alphabet = "ACGT"
rng = random.Random(7)
targets = []

while len(targets) < 2048:
    seq = "".join(rng.choice(alphabet) for _ in range(12))
    if any(sum(a != b for a, b in zip(seq, other)) < 3 for other in targets):
        continue
    targets.append(seq)

with open(targets_path, "w", encoding="utf-8") as out:
    for i, seq in enumerate(targets):
        out.write(f"sparse_{i}\t{seq}\tG{i % 17}\n")

touched = [targets[3], targets[701], targets[1777]]
with open(reads_path, "w", encoding="utf-8") as out:
    for i in range(4096):
        seq = touched[i % len(touched)]
        if i % 11 == 0:
            repl = "A" if seq[0] != "A" else "C"
            seq = repl + seq[1:]
        elif i % 17 == 0:
            seq = "N" + seq[1:]
        elif i % 29 == 0:
            seq = "TTTTTTTTTTTT"
        out.write(f"@sparse_{i}\n{seq}AAAA\n+\n{'I' * 16}\n")
PY

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/sparse_targets.tsv" \
  --reads "$TMPDIR/sparse_reads.fastq" \
  --sample-label sparse_hamming \
  --target-start 0 \
  --target-length 12 \
  --k 1 \
  --metric hamming \
  --ambiguity-policy best \
  --out "$TMPDIR/sparse_counts_hamming.tsv" \
  --summary "$TMPDIR/sparse_summary_hamming.json"

"$DOTMATCH_BIN" count \
  --targets "$TMPDIR/sparse_targets.tsv" \
  --reads "$TMPDIR/sparse_reads.fastq" \
  --sample-label sparse_hamming \
  --target-start 0 \
  --target-length 12 \
  --k 1 \
  --metric hamming \
  --ambiguity-policy best \
  --threads 4 \
  --out "$TMPDIR/sparse_counts_hamming_threads.tsv" \
  --summary "$TMPDIR/sparse_summary_hamming_threads.json"

diff -u "$TMPDIR/sparse_counts_hamming.tsv" "$TMPDIR/sparse_counts_hamming_threads.tsv"
grep '^sparse_3	' "$TMPDIR/sparse_counts_hamming_threads.tsv" | grep '	[1-9][0-9]*$' >/dev/null
grep '^sparse_701	' "$TMPDIR/sparse_counts_hamming_threads.tsv" | grep '	[1-9][0-9]*$' >/dev/null
grep '^sparse_1777	' "$TMPDIR/sparse_counts_hamming_threads.tsv" | grep '	[1-9][0-9]*$' >/dev/null
grep '"read_threads": 4' "$TMPDIR/sparse_summary_hamming_threads.json" >/dev/null

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
  --ambiguity-policy best \
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
  --ambiguity-policy best \
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

if "$DOTMATCH_BIN" count \
  --targets "$TMPDIR/one_target.tsv" \
  --reads "$TMPDIR/shifted.fastq" \
  --sample-label shifted \
  --target-start 0 \
  --target-length 4 \
  --k 1 \
  --auto-offset 1025 \
  --out "$TMPDIR/counts_offset_too_large.tsv" 2> "$TMPDIR/counts_offset_too_large.err"; then
  echo "oversized --auto-offset must fail" >&2
  exit 1
fi
grep -- '^--auto-offset must be <= 1024$' "$TMPDIR/counts_offset_too_large.err" >/dev/null

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
  --ambiguity-policy best \
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
  --ambiguity-policy best \
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
  --ambiguity-policy best \
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

cat > "$TMPDIR/duplicate_samples.tsv" <<SAMPLES
sample_id	fastq
plasmid	$TMPDIR/reads.fastq
plasmid	$TMPDIR/reads.fastq.gz
SAMPLES

if "$DOTMATCH_BIN" crispr-count \
  --library "$TMPDIR/targets.csv" \
  --samples "$TMPDIR/duplicate_samples.tsv" \
  --guide-start 0 \
  --guide-length 4 \
  --k 1 \
  --metric levenshtein \
  --out "$TMPDIR/duplicate_crispr.tsv" \
  2> "$TMPDIR/duplicate_crispr.err"; then
  echo "crispr-count should reject duplicate sample labels from sample sheets" >&2
  exit 1
fi
grep 'duplicate sample label: "plasmid"' "$TMPDIR/duplicate_crispr.err" >/dev/null
test ! -e "$TMPDIR/duplicate_crispr.tsv"

cat > "$TMPDIR/duplicate_guides.csv" <<'GUIDEDUP'
sgRNA,gRNA.sequence,Gene
gdup,ACGT,G0
gdup,TTTT,G1
GUIDEDUP

if "$DOTMATCH_BIN" crispr-count \
  --library "$TMPDIR/duplicate_guides.csv" \
  --samples "$TMPDIR/samples.tsv" \
  --guide-start 0 \
  --guide-length 4 \
  --k 1 \
  --metric levenshtein \
  --out "$TMPDIR/duplicate_guides.tsv" \
  2> "$TMPDIR/duplicate_guides.err"; then
  echo "crispr-count should reject duplicate guide IDs" >&2
  exit 1
fi
grep 'guide IDs must be unique; duplicate ID: "gdup"' "$TMPDIR/duplicate_guides.err" >/dev/null
test ! -e "$TMPDIR/duplicate_guides.tsv"

cat > "$TMPDIR/empty_guides.csv" <<'GUIDEEMPTY'
sgRNA,gRNA.sequence,Gene
,ACGT,G0
GUIDEEMPTY

if "$DOTMATCH_BIN" crispr-count \
  --library "$TMPDIR/empty_guides.csv" \
  --samples "$TMPDIR/samples.tsv" \
  --guide-start 0 \
  --guide-length 4 \
  --k 1 \
  --metric levenshtein \
  --out "$TMPDIR/empty_guides.tsv" \
  2> "$TMPDIR/empty_guides.err"; then
  echo "crispr-count should reject empty guide IDs" >&2
  exit 1
fi
grep 'target ID and sequence must be non-empty' "$TMPDIR/empty_guides.err" >/dev/null
test ! -e "$TMPDIR/empty_guides.tsv"

"$DOTMATCH_BIN" crispr-count \
  --library "$TMPDIR/targets.csv" \
  --samples "$TMPDIR/samples.tsv" \
  --guide-start 0 \
  --guide-length 4 \
  --k 1 \
  --metric levenshtein \
  --ambiguity-policy best \
  --threads 2 \
  --out "$TMPDIR/crispr_mageck.tsv" \
  --summary "$TMPDIR/crispr_qc.json"

diff -u "$TMPDIR/expected_mageck.tsv" "$TMPDIR/crispr_mageck.tsv"
grep '"k1_rescued_reads": 1' "$TMPDIR/crispr_qc.json" >/dev/null
grep '"percent_rescued_by_k1": 25.000000' "$TMPDIR/crispr_qc.json" >/dev/null

"$DOTMATCH_BIN" crispr-count --help > "$TMPDIR/crispr_count_help.txt"
grep '^DotMatch .* crispr-count$' "$TMPDIR/crispr_count_help.txt" >/dev/null
grep 'MAGeCK-ready' "$TMPDIR/crispr_count_help.txt" >/dev/null
grep 'Hamming supports k=0..3' "$TMPDIR/crispr_count_help.txt" >/dev/null

"$DOTMATCH_BIN" audit --help > "$TMPDIR/audit_help.txt"
grep '^DotMatch .* audit$' "$TMPDIR/audit_help.txt" >/dev/null
grep 'safe_at_hamming_k3' "$TMPDIR/audit_help.txt" >/dev/null
grep '2k+1 substitutions' "$TMPDIR/audit_help.txt" >/dev/null

cat > "$TMPDIR/samples_reordered.csv" <<SAMPLES
fastq_path,sample
$TMPDIR/reads.fastq,plasmid
$TMPDIR/reads.fastq.gz,esc
SAMPLES

cat > "$TMPDIR/library_aliases.csv" <<'LIBRARYALIASES'
sgRNA,sgRNA_sequence,gene_symbol
bc0,ACGT,G0
bc1,AGGT,G1
bc2,ACGA,G2
bc3,TTTT,G3
LIBRARYALIASES

"$DOTMATCH_BIN" crispr-count \
  --library "$TMPDIR/library_aliases.csv" \
  --samples "$TMPDIR/samples_reordered.csv" \
  --guide-start 0 \
  --guide-length 4 \
  --k 1 \
  --metric hamming \
  --ambiguity-policy best \
  --out "$TMPDIR/crispr_mageck_reordered.tsv"

diff -u "$TMPDIR/expected_mageck.tsv" "$TMPDIR/crispr_mageck_reordered.tsv"

"$DOTMATCH_BIN" audit \
  --targets "$TMPDIR/targets.csv" \
  --k 1 \
  --out-dir "$TMPDIR/audit"

grep '^targets	4$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '^safe_at_k1	no$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '^safe_at_hamming_k2	no$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '^safe_at_hamming_k3	no$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '^risk_pairs_for_k1	3$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '^risk_pairs_for_hamming_k2	6$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '^risk_pairs_for_hamming_k3	6$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '^ambiguous_query_variants_k1	14$' "$TMPDIR/audit/audit_summary.tsv" >/dev/null
grep '"audit_mode": "exact"' "$TMPDIR/audit/audit_summary.json" >/dev/null
grep '"k": 1' "$TMPDIR/audit/audit_summary.json" >/dev/null
grep '"safe_at_k1": false' "$TMPDIR/audit/audit_summary.json" >/dev/null
grep '"safe_at_hamming_k2": false' "$TMPDIR/audit/audit_summary.json" >/dev/null
grep '"safe_at_hamming_k3": false' "$TMPDIR/audit/audit_summary.json" >/dev/null
grep '"risk_pairs_for_k1": 3' "$TMPDIR/audit/audit_summary.json" >/dev/null
grep '"risk_pairs_for_hamming_k2": 6' "$TMPDIR/audit/audit_summary.json" >/dev/null
grep '"risk_pairs_for_hamming_k3": 6' "$TMPDIR/audit/audit_summary.json" >/dev/null
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
grep '^safe_at_hamming_k2	not_computed$' "$TMPDIR/audit_fast/audit_summary.tsv" >/dev/null
grep '^safe_at_hamming_k3	not_computed$' "$TMPDIR/audit_fast/audit_summary.tsv" >/dev/null
grep '^risk_pairs_for_k1	3$' "$TMPDIR/audit_fast/audit_summary.tsv" >/dev/null
grep '^risk_pairs_for_hamming_k2	not_computed$' "$TMPDIR/audit_fast/audit_summary.tsv" >/dev/null
grep '^risk_pairs_for_hamming_k3	not_computed$' "$TMPDIR/audit_fast/audit_summary.tsv" >/dev/null
grep '^ambiguous_query_variants_k1	14$' "$TMPDIR/audit_fast/audit_summary.tsv" >/dev/null
grep '"audit_mode": "fast"' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
grep '"safe_at_k1": false' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
grep '"safe_at_hamming_k2": null' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
grep '"safe_at_hamming_k3": null' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
grep '"safe_at_k2": null' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
grep '"risk_pairs_for_k2": null' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
grep '"risk_pairs_for_hamming_k2": null' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
grep '"risk_pairs_for_hamming_k3": null' "$TMPDIR/audit_fast/audit_summary.json" >/dev/null
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
