#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _DARWIN_C_SOURCE
#define _DARWIN_C_SOURCE 1
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include "qdalign.h"
#include "qdmetal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <pthread.h>
#include <errno.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <fcntl.h>
#include <dirent.h>
#include <regex.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <unistd.h>
#include <zlib.h>

#ifndef DOTMATCH_VERSION
#define DOTMATCH_VERSION "0.1.8"
#endif

#define MAX_AUTO_OFFSET 1024
#define MAX_BCL_CYCLE_CLUSTERS 100000000
#define MAX_BCL_READ_CYCLES 1024
#define MAX_BCL_TOTAL_CYCLES 1024
#define MAX_BCL_SAMPLE_ROWS 100000

typedef struct seq_record {
    char *id;
    char *seq;
    char *gene;
    size_t len;
} seq_record;

typedef struct seq_table {
    seq_record *records;
    size_t count;
    size_t cap;
} seq_table;

static double seconds_now(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

static size_t get_cpu_count(void) {
    long n = 1;
    n = sysconf(_SC_NPROCESSORS_ONLN);
    if (n < 1) n = 1;
    if ((size_t)n > 1024) n = 1024; /* cap for thread pools in this tool */
    return (size_t)n;
}

static void usage(const char *argv0) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s --help\n", argv0);
    fprintf(stderr, "  %s --version\n", argv0);
    fprintf(stderr, "  %s citation\n", argv0);
    fprintf(stderr, "  %s dist SEQ1 SEQ2\n", argv0);
    fprintf(stderr, "  %s leq K SEQ1 SEQ2\n", argv0);
    fprintf(stderr, "  %s assign K barcodes.txt reads.txt [--ambiguity-policy radius|best]\n", argv0);
    fprintf(stderr, "  %s match K targets.txt reads.txt [--ambiguity-policy radius|best]\n", argv0);
    fprintf(stderr, "  %s fastq-assign --barcodes barcodes.tsv --reads reads.fastq[.gz] --barcode-start N --barcode-length L --k 0|1 [--ambiguity-policy radius|best] --out assignments.tsv\n", argv0);
    fprintf(stderr, "  %s pair-count --left-targets left.tsv --right-targets right.tsv (--reads reads.fastq[.gz] | --left-reads left.fastq[.gz] --right-reads right.fastq[.gz]) --left-start N --left-length L --right-start N --right-length L --k 0|1|2 --metric hamming|levenshtein [--ambiguity-policy radius|best] --out pair_counts.tsv [--summary summary.json]\n", argv0);
    fprintf(stderr, "  %s demux --barcodes barcodes.tsv|barcodes.csv --reads reads.fastq[.gz] --barcode-start N --barcode-length L|auto --k 0|1|2 --metric hamming|levenshtein [--ambiguity-policy radius|best] [--max-correction-qual Q] --out-dir demux_dir [--summary qc.json]\n", argv0);
    fprintf(stderr, "  %s bcl-demux --run-folder RUN --sample-sheet SampleSheet.csv --out-dir demux_dir --barcode-mismatches 0|1|1,1 [--threads N] (0=auto) [--gzip-level 0..9] [--emit-index-fastqs] [--summary summary.json]\n", argv0);
    fprintf(stderr, "  %s bcl-validate --dotmatch-out DIR --truth-out DIR\n", argv0);
    fprintf(stderr, "  %s count --targets targets.tsv|targets.csv --reads reads.fastq[.gz] [--reads more.fastq.gz] --sample-label labels --target-start N --target-length L --k 0|1|2|3 --metric hamming|levenshtein [--hamming-index auto|query|precompute] [--max-correction-qual Q] [--ambiguity-policy radius|best] --offset-mode best|multi --out counts.tsv [--format dotmatch|mageck]\n", argv0);
    fprintf(stderr, "  %s crispr-count --library guides.tsv|guides.csv --samples samples.tsv|samples.csv --guide-start N --guide-length L --k 0|1|2|3 [--ambiguity-policy radius|best] [--backend auto|cpu|gpu-metal-experimental] --out counts.mageck.tsv [--summary qc.json] [--sample-qc sample_qc.tsv]\n", argv0);
    fprintf(stderr, "  %s guide-counter count --input reads.fastq[.gz]... --library guides.tsv|guides.csv --output prefix [--samples labels...] [--exact-match]\n", argv0);
    fprintf(stderr, "  %s inspect-unmatched --targets targets.tsv|targets.csv --reads reads.fastq[.gz] --target-start N --target-length L --k 0|1 --top N --out top_unmatched.tsv [--low-quality-threshold Q]\n", argv0);
    fprintf(stderr, "  %s audit --targets targets.tsv|targets.csv --k 0|1|2|3 --out-dir audit_dir [--audit-mode auto|exact|fast]\n", argv0);
    fprintf(stderr, "  %s validate --targets targets.tsv|targets.csv --reads reads.fastq[.gz] --target-start N --target-length L --k 0|1 [--metric hamming|levenshtein] [--indel-window 0|1] [--offset-mode best|multi] [--threads N] (0=auto) --oracle scan|edlib\n", argv0);
}

static void help_manual(FILE *out, const char *argv0) {
    fprintf(out, "DotMatch %s\n", DOTMATCH_VERSION);
    fprintf(out, "\n");
    fprintf(out, "Deterministic known-target short-DNA assignment for fixed read windows.\n");
    fprintf(out, "DotMatch is for cases where the expected guides, barcodes, primers, or\n");
    fprintf(out, "panel targets are already known. It is not a genome aligner, basecaller,\n");
    fprintf(out, "variant caller, adapter trimmer, cell/UMI quantifier, or screen statistics tool.\n");
    fprintf(out, "\n");
    fprintf(out, "Usage:\n");
    fprintf(out, "  %s --help\n", argv0);
    fprintf(out, "  %s --version\n", argv0);
    fprintf(out, "  %s citation\n", argv0);
    fprintf(out, "  %s <command> [options]\n", argv0);
    fprintf(out, "\n");
    fprintf(out, "Core commands:\n");
    fprintf(out, "  dist SEQ1 SEQ2\n");
    fprintf(out, "      Print the global edit distance between two short DNA strings.\n");
    fprintf(out, "  leq K SEQ1 SEQ2\n");
    fprintf(out, "      Print true when the edit distance is <= K, otherwise false.\n");
    fprintf(out, "  assign K targets.tsv reads.tsv [--ambiguity-policy radius|best]\n");
    fprintf(out, "      Assign tabular read sequences to known targets.\n");
    fprintf(out, "  fastq-assign --barcodes barcodes.tsv --reads reads.fastq[.gz] \\\n");
    fprintf(out, "      --barcode-start N --barcode-length L --k 0|1 --out assignments.tsv\n");
    fprintf(out, "      Write per-read FASTQ barcode assignments.\n");
    fprintf(out, "\n");
    fprintf(out, "Counting and demultiplexing:\n");
    fprintf(out, "  count --targets targets.tsv|targets.csv --reads reads.fastq[.gz] \\\n");
    fprintf(out, "      --sample-label sample --target-start N --target-length L \\\n");
    fprintf(out, "      --k 0|1|2|3 --metric hamming|levenshtein --out counts.tsv\n");
    fprintf(out, "      Count fixed-window target assignments. Add --format mageck for MAGeCK-style counts.\n");
    fprintf(out, "  crispr-count --library guides.tsv|guides.csv --samples samples.tsv|samples.csv \\\n");
    fprintf(out, "      --guide-start N --guide-length L --k 0|1|2|3 --out counts.mageck.tsv\n");
    fprintf(out, "      MAGeCK-ready guide counts. Sample sheets may use sample_id/sample and fastq/fastq_path columns.\n");
    fprintf(out, "  guide-counter count --input reads.fastq[.gz]... --library guides.tsv|guides.csv \\\n");
    fprintf(out, "      --output prefix [--samples labels...] [--exact-match]\n");
    fprintf(out, "      GuideCounter-compatible count, extended-counts, and stats outputs.\n");
    fprintf(out, "  demux --barcodes barcodes.tsv|barcodes.csv --reads reads.fastq[.gz] \\\n");
    fprintf(out, "      --barcode-start N --barcode-length L|auto --k 0|1|2 --out-dir demux_dir\n");
    fprintf(out, "      Split reads by fixed-position inline barcodes.\n");
    fprintf(out, "  pair-count --left-targets left.tsv --right-targets right.tsv \\\n");
    fprintf(out, "      (--reads reads.fastq[.gz] | --left-reads left.fastq[.gz] --right-reads right.fastq[.gz]) \\\n");
    fprintf(out, "      --left-start N --left-length L --right-start N --right-length L --out pair_counts.tsv\n");
    fprintf(out, "      Count independent fixed-window targets from one read or synchronized FASTQ mates.\n");
    fprintf(out, "\n");
    fprintf(out, "Diagnostics and validation:\n");
    fprintf(out, "  audit --targets targets.tsv|targets.csv --k K --out-dir audit_dir\n");
    fprintf(out, "      Report nearby target pairs that make correction ambiguous or unsafe.\n");
    fprintf(out, "      Exact audit includes Hamming k=2/k=3 safety fields; fast audit reports them as not_computed.\n");
    fprintf(out, "  inspect-unmatched --targets targets.tsv|targets.csv --reads reads.fastq[.gz] \\\n");
    fprintf(out, "      --target-start N --target-length L --k 0|1 --top N --out top_unmatched.tsv\n");
    fprintf(out, "      Summarize the most frequent unmatched read windows.\n");
    fprintf(out, "  validate --targets targets.tsv|targets.csv --reads reads.fastq[.gz] \\\n");
    fprintf(out, "      --target-start N --target-length L --k 0|1 --oracle scan|edlib\n");
    fprintf(out, "      Compare indexed assignment with a validation oracle. Installed packages should use --oracle scan.\n");
    fprintf(out, "  bcl-demux --run-folder RUN --sample-sheet SampleSheet.csv --out-dir demux_dir \\\n");
    fprintf(out, "      --barcode-mismatches 0|1|1,1 [--threads N] (0=auto) [--gzip-level 0..9]\n");
    fprintf(out, "      Classic per-cycle BCL parser milestone; not a production CBCL/NovaSeq replacement.\n");
    fprintf(out, "\n");
    fprintf(out, "Assignment outcomes:\n");
    fprintf(out, "  unique       exactly one target is compatible with the read window\n");
    fprintf(out, "  ambiguous    more than one target is compatible; not silently forced\n");
    fprintf(out, "  none         no target is close enough\n");
    fprintf(out, "  invalid      the requested read window cannot be extracted\n");
    fprintf(out, "\n");
    fprintf(out, "Defaults and conventions:\n");
    fprintf(out, "  --ambiguity-policy radius is the conservative default for assignment commands.\n");
    fprintf(out, "  --ambiguity-policy best is available only when best-distance compatibility is intended.\n");
    fprintf(out, "  --threads N (default 0 = auto-detect from CPU count) enables parallel read processing for count/demux/validate/bcl.\n");
    fprintf(out, "  --metric levenshtein supports k=0,1,2 for count/demux fixed windows; hamming supports k=0,1,2,3 fixed-length comparisons.\n");
    fprintf(out, "  N and IUPAC ambiguity symbols are treated as literal bytes, not wildcards.\n");
    fprintf(out, "\n");
    fprintf(out, "Examples:\n");
    fprintf(out, "  %s dist ACGT AGGT\n", argv0);
    fprintf(out, "  %s leq 1 ACGT AGGT\n", argv0);
    fprintf(out, "  %s citation\n", argv0);
    fprintf(out, "  %s count --targets guides.tsv --reads sample.fastq.gz --sample-label sample \\\n", argv0);
    fprintf(out, "      --target-start 23 --target-length 20 --k 1 --metric hamming --out counts.tsv\n");
    fprintf(out, "  %s demux --barcodes barcodes.tsv --reads pooled.fastq.gz \\\n", argv0);
    fprintf(out, "      --barcode-start 0 --barcode-length auto --k 0 --out-dir demuxed\n");
}

static void print_citation(FILE *out) {
    fprintf(out, "DotMatch citation\n");
    fprintf(out, "Software release: v%s\n", DOTMATCH_VERSION);
    fprintf(out, "Preferred citation:\n");
    fprintf(out, "O'Toole D. DotMatch: deterministic known-target short-DNA assignment for sequencing workflows. Software release v%s. https://github.com/dnncha/dotmatch\n", DOTMATCH_VERSION);
    fprintf(out, "Citation metadata: CITATION.cff\n");
    fprintf(out, "DOI: 10.5281/zenodo.20541628\n");
    fprintf(out, "DOI URL: https://doi.org/10.5281/zenodo.20541628\n");
}

static int help_requested(int argc, char **argv) {
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) return 1;
    }
    return 0;
}

static void count_help_manual(FILE *out, const char *argv0, int crispr_mode) {
    const char *cmd = crispr_mode ? "crispr-count" : "count";
    fprintf(out, "DotMatch %s %s\n\n", DOTMATCH_VERSION, cmd);
    if (crispr_mode) {
        fprintf(out, "Usage:\n");
        fprintf(out, "  %s crispr-count --library guides.tsv|guides.csv --samples samples.tsv|samples.csv \\\n", argv0);
        fprintf(out, "      --guide-start N --guide-length L --k 0|1|2|3 --metric hamming|levenshtein \\\n");
        fprintf(out, "      --out counts.mageck.tsv [--summary summary.json] [--sample-qc sample_qc.tsv]\n\n");
        fprintf(out, "CRISPR guide counting wrapper. Outputs are MAGeCK-ready by default and include\n");
        fprintf(out, "guide-level counts plus optional run and sample QC summaries.\n\n");
        fprintf(out, "Required inputs:\n");
        fprintf(out, "  --library PATH       Guide table. Common columns include id/sgRNA, sequence/guide_seq,\n");
        fprintf(out, "                       and gene/Gene; CSV and TSV are accepted.\n");
        fprintf(out, "  --samples PATH       Sample sheet with sample_id/sample and fastq/fastq_path columns.\n");
        fprintf(out, "  --guide-start N      Zero-based start of the guide window in each read.\n");
        fprintf(out, "  --guide-length L     Length of the guide window to extract.\n");
        fprintf(out, "  --out PATH           Count matrix output.\n");
    } else {
        fprintf(out, "Usage:\n");
        fprintf(out, "  %s count --targets targets.tsv|targets.csv --reads reads.fastq[.gz] [--reads more.fastq.gz] \\\n", argv0);
        fprintf(out, "      --sample-label labels --target-start N --target-length L --k 0|1|2|3 \\\n");
        fprintf(out, "      --metric hamming|levenshtein --out counts.tsv [--format dotmatch|mageck]\n\n");
        fprintf(out, "General fixed-window known-target counting for guides, barcodes, primers, or panels.\n\n");
        fprintf(out, "Required inputs:\n");
        fprintf(out, "  --targets PATH       Target table with id and sequence columns; CSV and TSV accepted.\n");
        fprintf(out, "  --reads PATH         FASTQ or FASTQ.GZ input. Repeat for multiple samples.\n");
        fprintf(out, "  --sample-label LIST  Comma-separated sample names matching --reads order.\n");
        fprintf(out, "  --target-start N     Zero-based start of the target window in each read.\n");
        fprintf(out, "  --target-length L    Length of the target window to extract.\n");
        fprintf(out, "  --out PATH           Count matrix output.\n");
    }
    fprintf(out, "\nAssignment options:\n");
    fprintf(out, "  --k 0|1|2|3              Maximum correction radius. Levenshtein supports k=0..2;\n");
    fprintf(out, "                           Hamming supports k=0..3 for fixed-length windows.\n");
    fprintf(out, "  --metric NAME            hamming for substitutions only; levenshtein for short indels.\n");
    fprintf(out, "  --ambiguity-policy NAME  radius keeps all targets within k; best keeps only nearest targets.\n");
    fprintf(out, "  --ambiguous NAME         discard, include, or separate ambiguous counts.\n");
    fprintf(out, "  --indel-window N         Levenshtein-only extraction slack for insertion/deletion rescue.\n");
    fprintf(out, "  --hamming-index NAME     auto, query, or precompute for Hamming counting.\n");
    fprintf(out, "  --backend NAME           auto, cpu, or gpu-metal-experimental (Darwin Metal Hamming k<=1).\n");
    fprintf(out, "  --progress               Emit read-processing progress to stderr (default on interactive stderr).\n");
    fprintf(out, "  --no-progress            Disable progress reporting.\n");
    fprintf(out, "  --progress-interval N    Report every N reads (default 250000).\n");
    fprintf(out, "  --metal-validate         With gpu-metal-experimental, shadow-count on CPU and require agreement.\n");
    fprintf(out, "                           Also enabled by DOTMATCH_METAL_VALIDATE=1.\n");
    fprintf(out, "\nOutputs:\n");
    fprintf(out, "  --summary PATH       JSON run summary with assignment rates and command metadata.\n");
    fprintf(out, "  --sample-qc PATH     Per-sample QC TSV; crispr-count defaults to sample_qc.tsv beside --out.\n");
    fprintf(out, "  --format mageck      Emit sgRNA/Gene columns for MAGeCK-style downstream analysis.\n");
    fprintf(out, "\nSafety:\n");
    fprintf(out, "  Before production Hamming k=2 or k=3 runs, use:\n");
    fprintf(out, "    %s audit --targets guides.tsv --k 3 --audit-mode exact --out-dir audit\n", argv0);
    fprintf(out, "  Proceed only when safe_at_hamming_k2 or safe_at_hamming_k3 is true for the radius you use.\n");
}

static void audit_help_manual(FILE *out, const char *argv0) {
    fprintf(out, "DotMatch %s audit\n\n", DOTMATCH_VERSION);
    fprintf(out, "Usage:\n");
    fprintf(out, "  %s audit --targets targets.tsv|targets.csv --k 0|1|2|3 --out-dir audit_dir \\\n", argv0);
    fprintf(out, "      [--audit-mode auto|exact|fast]\n\n");
    fprintf(out, "Audit a target library before correction-based assignment. The audit reports\n");
    fprintf(out, "duplicate targets and nearby target pairs that can make reads ambiguous.\n\n");
    fprintf(out, "Options:\n");
    fprintf(out, "  --targets PATH       Guide, barcode, primer, or target table; CSV and TSV accepted.\n");
    fprintf(out, "  --k 0|1|2|3          Radius to audit. Use k=3 to expose all Hamming k2/k3 fields.\n");
    fprintf(out, "  --out-dir DIR        Directory for audit_summary.tsv/json and risk-pair tables.\n");
    fprintf(out, "  --audit-mode exact   Pairwise exact audit. Required for Hamming k=2/k=3 safety fields.\n");
    fprintf(out, "  --audit-mode fast    Scalable duplicate/k1-style audit; Hamming k2/k3 fields are not_computed.\n");
    fprintf(out, "  --audit-mode auto    Exact for small libraries, fast for larger libraries.\n\n");
    fprintf(out, "Important fields:\n");
    fprintf(out, "  min_hamming_distance       Smallest fixed-length distance between any target pair.\n");
    fprintf(out, "  safe_at_hamming_k2         true only when Hamming radius 2 correction is unambiguous.\n");
    fprintf(out, "  safe_at_hamming_k3         true only when Hamming radius 3 correction is unambiguous.\n");
    fprintf(out, "  risk_pairs_for_hamming_k2  Number of target pairs too close for Hamming k=2.\n");
    fprintf(out, "  risk_pairs_for_hamming_k3  Number of target pairs too close for Hamming k=3.\n\n");
    fprintf(out, "Rule of thumb:\n");
    fprintf(out, "  Radius k is safe under radius-policy correction only when target pairs are\n");
    fprintf(out, "  separated by at least 2k+1 substitutions. If exact audit says unsafe,\n");
    fprintf(out, "  lower k, switch policy intentionally, or remove/merge conflicting targets.\n");
}

static char *xstrndup(const char *s, size_t n) {
    char *out = (char *)malloc(n + 1);
    if (out == NULL) return NULL;
    memcpy(out, s, n);
    out[n] = '\0';
    return out;
}

static void trim_line(char *line) {
    size_t n = strlen(line);
    while (n > 0 && (line[n - 1] == '\n' || line[n - 1] == '\r')) {
        line[--n] = '\0';
    }
}

static size_t trim_line_len(char *line, size_t n) {
    while (n > 0 && (line[n - 1] == '\n' || line[n - 1] == '\r')) {
        line[--n] = '\0';
    }
    return n;
}

static FILE *open_output_file(const char *path) {
    int fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
    if (fd < 0) return NULL;
    FILE *out = fdopen(fd, "w");
    if (out == NULL) {
        close(fd);
        return NULL;
    }
    return out;
}

static void uppercase_ascii(char *s) {
    for (; *s != '\0'; ++s) {
        if (*s >= 'a' && *s <= 'z') *s = (char)(*s - 'a' + 'A');
    }
}

static const char *status_name(int status) {
    switch (status) {
        case QDALN_MATCH_NONE:
            return "none";
        case QDALN_MATCH_UNIQUE:
            return "unique";
        case QDALN_MATCH_AMBIGUOUS:
            return "ambiguous";
        default:
            return "invalid";
    }
}

static int parse_size_value(const char *s, size_t *out) {
    if (s == NULL || s[0] == '-' || s[0] == '\0') return -1;
    char *end = NULL;
    errno = 0;
    unsigned long v = strtoul(s, &end, 10);
    if (errno == ERANGE || end == s || *end != '\0' || v > SIZE_MAX) return -1;
    *out = (size_t)v;
    return 0;
}

static int parse_int_value(const char *s, int *out) {
    if (s == NULL || s[0] == '\0') return -1;
    char *end = NULL;
    errno = 0;
    long v = strtol(s, &end, 10);
    if (errno == ERANGE || end == s || *end != '\0' || v < INT_MIN || v > INT_MAX) return -1;
    *out = (int)v;
    return 0;
}

static int parse_double_value(const char *s, double *out) {
    if (s == NULL || s[0] == '\0') return -1;
    char *end = NULL;
    errno = 0;
    double v = strtod(s, &end);
    if (errno == ERANGE || end == s || *end != '\0' || !isfinite(v)) return -1;
    *out = v;
    return 0;
}

static int offset_count_for_range(size_t range, size_t *out) {
    if (range > MAX_AUTO_OFFSET || range > (SIZE_MAX - 1) / 2) return -1;
    *out = range * 2 + 1;
    return 0;
}

static int checked_mul_size(size_t a, size_t b, size_t *out) {
    if (a != 0 && b > SIZE_MAX / a) return -1;
    *out = a * b;
    return 0;
}

static size_t alloc_count_or_one(size_t count) {
    return count == 0 ? 1 : count;
}

static void free_table(seq_table *table) {
    for (size_t i = 0; i < table->count; ++i) {
        free(table->records[i].id);
        free(table->records[i].seq);
        free(table->records[i].gene);
    }
    free(table->records);
    table->records = NULL;
    table->count = 0;
    table->cap = 0;
}

static int push_record_gene(seq_table *table, const char *id, size_t id_len, const char *seq, size_t seq_len,
                            const char *gene, size_t gene_len);

static int push_record(seq_table *table, const char *id, size_t id_len, const char *seq, size_t seq_len) {
    return push_record_gene(table, id, id_len, seq, seq_len, "", 0);
}

static int push_record_gene(seq_table *table, const char *id, size_t id_len, const char *seq, size_t seq_len,
                            const char *gene, size_t gene_len) {
    if (table->count == table->cap) {
        size_t next_cap = table->cap == 0 ? 16 : table->cap * 2;
        seq_record *next = (seq_record *)realloc(table->records, next_cap * sizeof(seq_record));
        if (next == NULL) return -1;
        table->records = next;
        table->cap = next_cap;
    }

    seq_record *r = &table->records[table->count];
    r->id = xstrndup(id, id_len);
    r->seq = xstrndup(seq, seq_len);
    r->gene = xstrndup(gene, gene_len);
    if (r->id == NULL || r->seq == NULL || r->gene == NULL) {
        free(r->id);
        free(r->seq);
        free(r->gene);
        return -1;
    }
    uppercase_ascii(r->seq);
    r->len = seq_len;
    ++table->count;
    return 0;
}

static int read_table(const char *path, seq_table *table) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return -1;

    char buf[8192];
    size_t row = 0;
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        trim_line(buf);
        if (buf[0] == '\0') continue;

        char *tab = strchr(buf, '\t');
        const char *id = NULL;
        size_t id_len = 0;
        const char *seq = NULL;
        if (tab != NULL) {
            *tab = '\0';
            id = buf;
            id_len = strlen(id);
            seq = tab + 1;
            if (id_len == 0 || seq[0] == '\0') {
                fprintf(stderr, "%s:%zu: record ID and sequence must be non-empty\n", path, row + 1);
                fclose(fp);
                return -1;
            }
        } else {
            char id_buf[32];
            int n = snprintf(id_buf, sizeof(id_buf), "%zu", row);
            if (n < 0 || (size_t)n >= sizeof(id_buf)) {
                fclose(fp);
                return -1;
            }
            if (push_record(table, id_buf, (size_t)n, buf, strlen(buf)) != 0) {
                fclose(fp);
                return -1;
            }
            ++row;
            continue;
        }

        if (push_record(table, id, id_len, seq, strlen(seq)) != 0) {
            fclose(fp);
            return -1;
        }
        ++row;
    }

    if (ferror(fp)) {
        fclose(fp);
        return -1;
    }
    fclose(fp);
    return 0;
}

static size_t split_fields(char *line, char delim, char **fields, size_t max_fields) {
    size_t n = 0;
    char *p = line;
    while (n < max_fields) {
        fields[n++] = p;
        char *next = strchr(p, delim);
        if (next == NULL) break;
        *next = '\0';
        p = next + 1;
    }
    return n;
}

static int field_eq(const char *a, const char *b) {
    while (*a != '\0' && *b != '\0') {
        char ca = *a;
        char cb = *b;
        if (ca >= 'A' && ca <= 'Z') ca = (char)(ca - 'A' + 'a');
        if (cb >= 'A' && cb <= 'Z') cb = (char)(cb - 'A' + 'a');
        if (ca != cb) return 0;
        ++a;
        ++b;
    }
    return *a == '\0' && *b == '\0';
}

static int find_column(char **fields, size_t n, const char *a, const char *b, const char *c) {
    for (size_t i = 0; i < n; ++i) {
        if (field_eq(fields[i], a) || (b != NULL && field_eq(fields[i], b)) || (c != NULL && field_eq(fields[i], c))) {
            return (int)i;
        }
    }
    return -1;
}

#include "target_table.h"

static int run_batch(const char *argv0, int argc, char **argv, const char *mode) {
    if (argc != 5 && argc != 7) {
        usage(argv0);
        return 2;
    }

    int k = 0;
    if (parse_int_value(argv[2], &k) != 0 || k < 0) {
        usage(argv0);
        return 2;
    }

    int radius_policy = 1;
    if (argc == 7) {
        if (strcmp(argv[5], "--ambiguity-policy") != 0) {
            usage(argv0);
            return 2;
        }
        if (strcmp(argv[6], "radius") == 0) {
            radius_policy = 1;
        } else if (strcmp(argv[6], "best") == 0) {
            radius_policy = 0;
        } else {
            usage(argv0);
            return 2;
        }
    }

    seq_table targets = {0};
    seq_table reads = {0};
    int rc = 1;

    if (read_table(argv[3], &targets) != 0 || read_table(argv[4], &reads) != 0) {
        fprintf(stderr, "failed to read input files\n");
        goto done;
    }

    const char **read_ptrs = (const char **)malloc(reads.count * sizeof(char *));
    const char **target_ptrs = (const char **)malloc(targets.count * sizeof(char *));
    size_t *read_lens = (size_t *)malloc(reads.count * sizeof(size_t));
    size_t *target_lens = (size_t *)malloc(targets.count * sizeof(size_t));
    qdaln_match_result *results = (qdaln_match_result *)malloc(reads.count * sizeof(qdaln_match_result));
    if ((reads.count != 0 && (read_ptrs == NULL || read_lens == NULL || results == NULL)) ||
        (targets.count != 0 && (target_ptrs == NULL || target_lens == NULL))) {
        fprintf(stderr, "out of memory\n");
        free(read_ptrs);
        free(target_ptrs);
        free(read_lens);
        free(target_lens);
        free(results);
        goto done;
    }

    for (size_t i = 0; i < reads.count; ++i) {
        read_ptrs[i] = reads.records[i].seq;
        read_lens[i] = reads.records[i].len;
    }
    for (size_t i = 0; i < targets.count; ++i) {
        target_ptrs[i] = targets.records[i].seq;
        target_lens[i] = targets.records[i].len;
    }

    if (qdaln_match_many(read_ptrs, read_lens, reads.count, target_ptrs, target_lens, targets.count, k, results) != 0) {
        fprintf(stderr, "batch match failed\n");
        free(read_ptrs);
        free(target_ptrs);
        free(read_lens);
        free(target_lens);
        free(results);
        goto done;
    }

    printf("mode\tread_id\tread_seq\ttarget_index\ttarget_seq\tdistance\tstatus\tmatch_count\tsecond_best_distance\n");
    for (size_t i = 0; i < reads.count; ++i) {
        qdaln_match_result r = results[i];
        if (radius_policy && r.status == QDALN_MATCH_UNIQUE && r.match_count > 1) {
            r.status = QDALN_MATCH_AMBIGUOUS;
        }
        const char *target_seq = r.target_index >= 0 ? targets.records[r.target_index].seq : "";
        printf("%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\n",
               mode, reads.records[i].id, reads.records[i].seq, r.target_index,
               target_seq, r.best_distance, status_name(r.status), r.match_count,
               r.second_best_distance);
    }

    free(read_ptrs);
    free(target_ptrs);
    free(read_lens);
    free(target_lens);
    free(results);
    rc = 0;

done:
    free_table(&targets);
    free_table(&reads);
    return rc;
}

typedef struct fastq_reader {
    FILE *fp;
    gzFile gz;
    unsigned char *gz_buf;
    size_t gz_pos;
    size_t gz_len;
    size_t gz_cap;
    int gz_eof;
    int is_gz;
} fastq_reader;

static int ends_with(const char *s, const char *suffix) {
    size_t n = strlen(s);
    size_t m = strlen(suffix);
    return n >= m && strcmp(s + n - m, suffix) == 0;
}

static int fastq_reader_open(fastq_reader *reader, const char *path) {
    memset(reader, 0, sizeof(*reader));
    reader->is_gz = ends_with(path, ".gz");
    if (reader->is_gz) {
        reader->gz = gzopen(path, "rb");
        if (reader->gz == NULL) return -1;
        gzbuffer(reader->gz, 1 << 20);
        reader->gz_cap = 1 << 20;
        reader->gz_buf = (unsigned char *)malloc(reader->gz_cap);
        if (reader->gz_buf == NULL) {
            gzclose(reader->gz);
            memset(reader, 0, sizeof(*reader));
            return -1;
        }
        return 0;
    }
    reader->fp = fopen(path, "r");
    return reader->fp == NULL ? -1 : 0;
}

static void fastq_reader_close(fastq_reader *reader) {
    if (reader->is_gz && reader->gz != NULL) gzclose(reader->gz);
    if (!reader->is_gz && reader->fp != NULL) fclose(reader->fp);
    free(reader->gz_buf);
    memset(reader, 0, sizeof(*reader));
}

static int fastq_getline_len(fastq_reader *reader, char *buf, size_t cap, size_t *len_out) {
    if (reader->is_gz) {
        if (cap == 0) return -1;
        size_t out = 0;
        for (;;) {
            if (reader->gz_pos == reader->gz_len) {
                if (reader->gz_eof) {
                    if (out == 0) return 0;
                    buf[out] = '\0';
                    if (len_out != NULL) *len_out = out;
                    return 1;
                }
                int n = gzread(reader->gz, reader->gz_buf, (unsigned int)reader->gz_cap);
                if (n < 0) return -1;
                if (n == 0) {
                    reader->gz_eof = 1;
                    continue;
                }
                reader->gz_pos = 0;
                reader->gz_len = (size_t)n;
            }

            size_t avail = reader->gz_len - reader->gz_pos;
            unsigned char *src = reader->gz_buf + reader->gz_pos;
            unsigned char *nl = (unsigned char *)memchr(src, '\n', avail);
            size_t take = nl == NULL ? avail : (size_t)(nl - src) + 1;
            if (out + take >= cap) return -1;
            memcpy(buf + out, src, take);
            out += take;
            reader->gz_pos += take;
            if (nl != NULL) {
                buf[out] = '\0';
                if (len_out != NULL) *len_out = out;
                return 1;
            }
        }
        return 1;
    }
    if (fgets(buf, (int)cap, reader->fp) == NULL) {
        return ferror(reader->fp) ? -1 : 0;
    }
    if (len_out != NULL) *len_out = strlen(buf);
    return 1;
}

static int fastq_skip_line_len(fastq_reader *reader, int *first_char_out, size_t *len_out) {
    int first = -1;
    size_t len = 0;
    unsigned char last = 0;
    int have_last = 0;

    if (reader->is_gz) {
        for (;;) {
            if (reader->gz_pos == reader->gz_len) {
                if (reader->gz_eof) {
                    if (len == 0) return 0;
                    if (have_last && last == '\r') --len;
                    if (first_char_out != NULL) *first_char_out = first;
                    if (len_out != NULL) *len_out = len;
                    return 1;
                }
                int n = gzread(reader->gz, reader->gz_buf, (unsigned int)reader->gz_cap);
                if (n < 0) return -1;
                if (n == 0) {
                    reader->gz_eof = 1;
                    continue;
                }
                reader->gz_pos = 0;
                reader->gz_len = (size_t)n;
            }

            size_t avail = reader->gz_len - reader->gz_pos;
            unsigned char *src = reader->gz_buf + reader->gz_pos;
            if (first < 0 && avail > 0) first = src[0];
            unsigned char *nl = (unsigned char *)memchr(src, '\n', avail);
            size_t take = nl == NULL ? avail : (size_t)(nl - src);
            if (take > 0) {
                last = src[take - 1];
                have_last = 1;
            }
            len += take;
            reader->gz_pos += take + (nl == NULL ? 0 : 1);
            if (nl != NULL) {
                if (have_last && last == '\r') --len;
                if (first_char_out != NULL) *first_char_out = first;
                if (len_out != NULL) *len_out = len;
                return 1;
            }
        }
    }

    int c = 0;
    while ((c = fgetc(reader->fp)) != EOF) {
        if (first < 0) first = c;
        if (c == '\n') {
            if (have_last && last == '\r') --len;
            if (first_char_out != NULL) *first_char_out = first;
            if (len_out != NULL) *len_out = len;
            return 1;
        }
        last = (unsigned char)c;
        have_last = 1;
        ++len;
    }
    if (ferror(reader->fp)) return -1;
    if (len == 0) return 0;
    if (have_last && last == '\r') --len;
    if (first_char_out != NULL) *first_char_out = first;
    if (len_out != NULL) *len_out = len;
    return 1;
}

static int fastq_read_record_len(fastq_reader *reader, char *header, char *seq, char *plus, char *qual,
                                 size_t cap, size_t *seq_len_out) {
    size_t header_len = 0;
    size_t seq_len = 0;
    size_t plus_len = 0;
    size_t qual_len = 0;
    int got = fastq_getline_len(reader, header, cap, &header_len);
    if (got <= 0) return got;
    if (fastq_getline_len(reader, seq, cap, &seq_len) != 1 ||
        fastq_getline_len(reader, plus, cap, &plus_len) != 1 ||
        fastq_getline_len(reader, qual, cap, &qual_len) != 1) {
        return -1;
    }
    header_len = trim_line_len(header, header_len);
    seq_len = trim_line_len(seq, seq_len);
    plus_len = trim_line_len(plus, plus_len);
    qual_len = trim_line_len(qual, qual_len);
    (void)header_len;
    (void)plus_len;
    if (header[0] != '@' || plus[0] != '+') return -1;
    if (seq_len != qual_len) return -1;
    if (seq_len_out != NULL) *seq_len_out = seq_len;
    return 1;
}

static int fastq_read_sequence_record_len(fastq_reader *reader, char *seq, size_t cap, size_t *seq_len_out) {
    int header_first = 0;
    size_t header_len = 0;
    int got = fastq_skip_line_len(reader, &header_first, &header_len);
    if (got <= 0) return got;
    if (header_first != '@') return -1;

    size_t seq_len = 0;
    if (fastq_getline_len(reader, seq, cap, &seq_len) != 1) return -1;
    seq_len = trim_line_len(seq, seq_len);

    int plus_first = 0;
    size_t plus_len = 0;
    size_t qual_len = 0;
    if (fastq_skip_line_len(reader, &plus_first, &plus_len) != 1 ||
        fastq_skip_line_len(reader, NULL, &qual_len) != 1) {
        return -1;
    }
    (void)header_len;
    (void)plus_len;
    if (plus_first != '+') return -1;
    if (seq_len != qual_len) return -1;
    if (seq_len_out != NULL) *seq_len_out = seq_len;
    return 1;
}

static void fastq_read_id(const char *header, char *out, size_t out_cap) {
    const char *start = header[0] == '@' ? header + 1 : header;
    size_t n = 0;
    while (start[n] != '\0' && start[n] != ' ' && start[n] != '\t') ++n;
    if (n >= out_cap) n = out_cap - 1;
    memcpy(out, start, n);
    out[n] = '\0';
}

/*
 * Paired FASTQ files commonly use either a shared Illumina identifier followed
 * by a read-number field, or a terminal /1 and /2 suffix. Compare the stable
 * identifier in both forms while retaining fastq_read_id behavior for
 * single-read commands.
 */
static void fastq_pair_read_id(const char *header, char *out, size_t out_cap) {
    fastq_read_id(header, out, out_cap);
    size_t n = strlen(out);
    if (n >= 2 && out[n - 2] == '/' && (out[n - 1] == '1' || out[n - 1] == '2')) {
        out[n - 2] = '\0';
    }
}

static void print_fastq_row(FILE *out, const seq_table *targets, const char *read_id,
                            const char *observed, qdaln_match_result r) {
    const char *target_id = "";
    const char *target_seq = "";
    if (r.target_index >= 0) {
        target_id = targets->records[r.target_index].id;
        target_seq = targets->records[r.target_index].seq;
    }
    fprintf(out, "%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%s\n",
            read_id, observed, r.target_index, target_id, target_seq, r.best_distance,
            r.second_best_distance, r.match_count, status_name(r.status));
}

typedef struct string_list {
    char **items;
    size_t count;
    size_t cap;
} string_list;

static void free_string_list(string_list *list) {
    for (size_t i = 0; i < list->count; ++i) free(list->items[i]);
    free(list->items);
    list->items = NULL;
    list->count = 0;
    list->cap = 0;
}

static int push_string(string_list *list, const char *s) {
    if (list->count == list->cap) {
        size_t next_cap = list->cap == 0 ? 4 : list->cap * 2;
        char **next = (char **)realloc(list->items, next_cap * sizeof(char *));
        if (next == NULL) return -1;
        list->items = next;
        list->cap = next_cap;
    }
    list->items[list->count] = xstrndup(s, strlen(s));
    if (list->items[list->count] == NULL) return -1;
    ++list->count;
    return 0;
}

static int split_string_list(string_list *list, const char *s, char delim) {
    const char *start = s;
    for (;;) {
        const char *p = strchr(start, delim);
        size_t n = p == NULL ? strlen(start) : (size_t)(p - start);
        char *tmp = xstrndup(start, n);
        if (tmp == NULL) return -1;
        int rc = push_string(list, tmp);
        free(tmp);
        if (rc != 0) return -1;
        if (p == NULL) break;
        start = p + 1;
    }
    return 0;
}

static int string_list_has_duplicates(const string_list *list, const char **first, const char **second) {
    for (size_t i = 0; i < list->count; ++i) {
        for (size_t j = i + 1; j < list->count; ++j) {
            if (strcmp(list->items[i], list->items[j]) == 0) {
                if (first != NULL) *first = list->items[i];
                if (second != NULL) *second = list->items[j];
                return 1;
            }
        }
    }
    return 0;
}

static int validate_unique_sample_labels(const string_list *labels, const char *option_name) {
    const char *first = NULL;
    const char *second = NULL;
    if (!string_list_has_duplicates(labels, &first, &second)) return 0;
    fprintf(stderr, "%s values must be unique; duplicate sample label: \"%s\"\n", option_name, first);
    (void)second;
    return -1;
}

typedef struct {
    const char *id;
    size_t index;
} seq_id_entry;

static int compare_seq_id_entry(const void *a, const void *b) {
    const seq_id_entry *ea = (const seq_id_entry *)a;
    const seq_id_entry *eb = (const seq_id_entry *)b;
    int cmp = strcmp(ea->id, eb->id);
    if (cmp != 0) return cmp;
    if (ea->index < eb->index) return -1;
    if (ea->index > eb->index) return 1;
    return 0;
}

static int validate_unique_seq_ids(const seq_table *targets, const char *record_kind) {
    if (targets->count < 2) return 0;
    seq_id_entry *entries = (seq_id_entry *)malloc(targets->count * sizeof(seq_id_entry));
    if (entries == NULL) return -1;
    for (size_t i = 0; i < targets->count; ++i) {
        entries[i].id = targets->records[i].id;
        entries[i].index = i;
    }
    qsort(entries, targets->count, sizeof(seq_id_entry), compare_seq_id_entry);
    for (size_t i = 1; i < targets->count; ++i) {
        if (strcmp(entries[i - 1].id, entries[i].id) == 0) {
            fprintf(stderr, "%s IDs must be unique; duplicate ID: \"%s\"\n", record_kind, entries[i].id);
            free(entries);
            return -2;
        }
    }
    free(entries);
    return 0;
}

static int read_samples_file(const char *path, string_list *labels, string_list *reads) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return -1;
    char buf[8192];
    int sample_col = 0;
    int reads_col = 1;
    int first_data = 1;
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        trim_line(buf);
        if (buf[0] == '\0' || buf[0] == '#') continue;

        char delim = strchr(buf, ',') != NULL && strchr(buf, '\t') == NULL ? ',' : '\t';
        char *fields[32];
        size_t nf = split_fields(buf, delim, fields, sizeof(fields) / sizeof(fields[0]));
        if (first_data) {
            int maybe_sample = find_column(fields, nf, "sample_id", "sample", "label");
            if (maybe_sample < 0) maybe_sample = find_column(fields, nf, "sample_name", "name", "id");
            int maybe_reads = find_column(fields, nf, "fastq", "fastq_path", "reads");
            if (maybe_reads < 0) maybe_reads = find_column(fields, nf, "read", "path", "file");
            if (maybe_sample >= 0 && maybe_reads >= 0) {
                sample_col = maybe_sample;
                reads_col = maybe_reads;
                first_data = 0;
                continue;
            }
        }
        first_data = 0;

        if ((size_t)sample_col >= nf || (size_t)reads_col >= nf) {
            fclose(fp);
            return -1;
        }
        char *sample = fields[sample_col];
        char *path_field = fields[reads_col];
        if (sample[0] == '\0' || path_field[0] == '\0') {
            fclose(fp);
            return -1;
        }
        if (push_string(labels, sample) != 0 || push_string(reads, path_field) != 0) {
            fclose(fp);
            return -1;
        }
    }
    if (ferror(fp)) {
        fclose(fp);
        return -1;
    }
    fclose(fp);
    return 0;
}

static const char *path_basename(const char *path) {
    const char *slash = strrchr(path, '/');
    return slash == NULL ? path : slash + 1;
}

static int one_delete_matches(const char *longer, size_t longer_len, const char *shorter, size_t shorter_len) {
    if (longer_len != shorter_len + 1) return 0;
    size_t i = 0;
    size_t j = 0;
    int edits = 0;
    while (i < longer_len && j < shorter_len) {
        if (longer[i] == shorter[j]) {
            ++i;
            ++j;
        } else {
            ++edits;
            if (edits > 1) return 0;
            ++i;
        }
    }
    return 1;
}

static int correction_kind(const char *observed, size_t observed_len, const char *target, size_t target_len, int d) {
    if (d == 0) return 0;
    if (d != 1) return 4;
    if (observed_len == target_len) return 1;
    if (one_delete_matches(observed, observed_len, target, target_len)) return 2;
    if (one_delete_matches(target, target_len, observed, observed_len)) return 3;
    return 4;
}

static const char *correction_name(int kind) {
    switch (kind) {
        case 0:
            return "exact";
        case 1:
            return "substitution";
        case 2:
            return "insertion";
        case 3:
            return "deletion";
        default:
            return "other";
    }
}

typedef struct count_stats {
    unsigned long long total;
    unsigned long long unique;
    unsigned long long exact;
    unsigned long long corrected;
    unsigned long long ambiguous;
    unsigned long long unmatched;
    unsigned long long invalid;
    unsigned long long candidates_considered;
    unsigned long long candidates_verified;
} count_stats;

typedef struct count_progress {
    int enabled;
    const char *sample_label;
    size_t interval_reads;
    unsigned long long reads_done;
    double start_seconds;
    double last_report_seconds;
    pthread_mutex_t lock;
} count_progress;

typedef struct sample_qc_metrics {
    double assignment_rate;
    double ambiguous_rate;
    double no_match_rate;
    double invalid_rate;
    double coverage_fraction;
    double zero_count_fraction;
    double gini_index;
    double top_1pct_fraction;
} sample_qc_metrics;

static void count_progress_init(count_progress *progress, const char *sample_label, size_t interval_reads) {
    if (progress == NULL) return;
    progress->enabled = 1;
    progress->sample_label = sample_label;
    progress->interval_reads = interval_reads > 0 ? interval_reads : 250000;
    progress->reads_done = 0;
    progress->start_seconds = seconds_now();
    progress->last_report_seconds = progress->start_seconds;
    pthread_mutex_init(&progress->lock, NULL);
}

static void count_progress_fini(count_progress *progress) {
    if (progress == NULL) return;
    pthread_mutex_destroy(&progress->lock);
}

static void count_progress_tick(count_progress *progress) {
    if (progress == NULL || !progress->enabled) return;
    pthread_mutex_lock(&progress->lock);
    ++progress->reads_done;
    unsigned long long total = progress->reads_done;
    double now = seconds_now();
    if (total % progress->interval_reads == 0 || now - progress->last_report_seconds >= 5.0) {
        double elapsed = now - progress->start_seconds;
        double rate = elapsed > 0.0 ? (double)total / elapsed : 0.0;
        fprintf(stderr, "dotmatch: %s: %llu reads (%.0f reads/s)\n", progress->sample_label, total, rate);
        progress->last_report_seconds = now;
    }
    pthread_mutex_unlock(&progress->lock);
}

static void count_progress_finish(count_progress *progress) {
    if (progress == NULL || !progress->enabled || progress->reads_done == 0) return;
    pthread_mutex_lock(&progress->lock);
    double elapsed = seconds_now() - progress->start_seconds;
    fprintf(stderr, "dotmatch: %s: finished %llu reads in %.1fs\n", progress->sample_label, progress->reads_done,
            elapsed);
    pthread_mutex_unlock(&progress->lock);
}

static int derive_output_sibling_path(const char *out_path, const char *filename, char *buf, size_t cap) {
    if (out_path == NULL || filename == NULL || buf == NULL || cap == 0) return -1;
    const char *slash = strrchr(out_path, '/');
    if (slash == NULL) {
        int n = snprintf(buf, cap, "%s", filename);
        return n < 0 || (size_t)n >= cap ? -1 : 0;
    }
    size_t dir_len = (size_t)(slash - out_path) + 1;
    int n = snprintf(buf, cap, "%.*s%s", (int)dir_len, out_path, filename);
    return n < 0 || (size_t)n >= cap ? -1 : 0;
}

typedef enum count_metric {
    COUNT_METRIC_LEVENSHTEIN = 0,
    COUNT_METRIC_HAMMING = 1
} count_metric;

static const char *metric_name(count_metric metric) {
    return metric == COUNT_METRIC_HAMMING ? "hamming" : "levenshtein";
}

typedef enum hamming_index_strategy {
    HAMMING_INDEX_QUERY = 0,
    HAMMING_INDEX_PRECOMPUTE = 1,
    HAMMING_INDEX_AUTO = 2
} hamming_index_strategy;

typedef enum count_backend_mode {
    COUNT_BACKEND_AUTO = 0,
    COUNT_BACKEND_CPU = 1,
    COUNT_BACKEND_METAL = 2
} count_backend_mode;

typedef enum offset_mode {
    OFFSET_MODE_BEST = 0,
    OFFSET_MODE_MULTI = 1
} offset_mode;

static const char *offset_mode_name(offset_mode mode) {
    return mode == OFFSET_MODE_MULTI ? "multi" : "best";
}

typedef struct offset_list {
    size_t *items;
    size_t count;
    size_t cap;
} offset_list;

static void free_offset_list(offset_list *list) {
    if (list == NULL) return;
    free(list->items);
    list->items = NULL;
    list->count = 0;
    list->cap = 0;
}

static int offset_list_contains(const offset_list *list, size_t offset) {
    for (size_t i = 0; i < list->count; ++i) {
        if (list->items[i] == offset) return 1;
    }
    return 0;
}

static int push_offset_unique(offset_list *list, size_t offset) {
    if (offset_list_contains(list, offset)) return 0;
    if (list->count == list->cap) {
        size_t next_cap = list->cap == 0 ? 8 : list->cap * 2;
        size_t *next = (size_t *)realloc(list->items, next_cap * sizeof(size_t));
        if (next == NULL) return -1;
        list->items = next;
        list->cap = next_cap;
    }
    list->items[list->count++] = offset;
    return 0;
}

static size_t first_selected_offset(const offset_list *list, size_t fallback) {
    return list != NULL && list->count != 0 ? list->items[0] : fallback;
}

typedef struct hamming_lookup_entry {
    uint64_t code;
    int target_index;
    int match_count;
} hamming_lookup_entry;

typedef struct hamming_seed_entry {
    uint64_t code;
    int target_index;
    int next;
    unsigned char seed_id;
} hamming_seed_entry;

typedef struct hamming_lookup {
    hamming_lookup_entry *exact;
    hamming_lookup_entry *mismatch;
    hamming_seed_entry *seeds;
    int *seed_heads;
    uint64_t *target_codes;
    size_t exact_cap;
    size_t mismatch_cap;
    size_t seed_hash_cap;
    size_t n_seeds;
    size_t target_len;
    size_t seed0_len;
    int ready;
    int seed_ready;
} hamming_lookup;

typedef struct levenshtein1_lookup {
    hamming_lookup_entry *exact;
    hamming_lookup_entry *substitution;
    hamming_lookup_entry *target_deletion;
    hamming_lookup_entry *target_insertion;
    size_t exact_cap;
    size_t substitution_cap;
    size_t target_deletion_cap;
    size_t target_insertion_cap;
    size_t target_len;
    int ready;
} levenshtein1_lookup;

static const char *hamming_lookup_kind(const hamming_lookup *lookup) {
    if (lookup == NULL || !lookup->ready) return "query";
    if (lookup->mismatch != NULL && lookup->mismatch_cap != 0) return "precompute";
    if (lookup->seed_ready) return "seed";
    return "exact";
}

static size_t next_pow2_local(size_t n) {
    size_t p = 1;
    while (p < n && p <= (SIZE_MAX >> 1)) p <<= 1;
    return p < n ? n : p;
}

static inline size_t code_hash_local(uint64_t code, size_t len, size_t cap) {
    uint64_t x = code ^ ((uint64_t)len * 0x9e3779b97f4a7c15ULL);
    x *= 0x9e3779b97f4a7c15ULL;
    x ^= x >> 32;
    return (size_t)x & (cap - 1);
}

static inline size_t seed_hash_local(uint64_t code, size_t len, unsigned char seed_id, size_t cap) {
    return code_hash_local(code ^ ((uint64_t)seed_id * 0x517cc1b727220a95ULL), len + seed_id * 37U, cap);
}

static uint64_t code_low_mask_local(size_t len) {
    if (len == 0) return 0;
    if (len >= 32) return UINT64_MAX;
    return (1ULL << (2 * len)) - 1ULL;
}

static uint64_t code_segment_local(uint64_t code, size_t start, size_t len) {
    return (code >> (2 * start)) & code_low_mask_local(len);
}

static int hamming_code_distance_local(uint64_t a, uint64_t b, size_t len) {
    uint64_t diff = a ^ b;
    diff |= diff >> 1;
    diff &= code_low_mask_local(len);
    diff &= 0x5555555555555555ULL;
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(diff);
#else
    int d = 0;
    while (diff != 0) {
        d += (int)(diff & 1ULL);
        diff >>= 2;
    }
    return d;
#endif
}

static int dna2_code_local(const char *s, size_t len, uint64_t *code_out) {
    if (s == NULL && len != 0) return 0;
    if (len > 32) return 0;
    uint64_t code = 0;
    for (size_t i = 0; i < len; ++i) {
        uint64_t v;
        switch (s[i]) {
            case 'A':
                v = 0;
                break;
            case 'C':
                v = 1;
                break;
            case 'G':
                v = 2;
                break;
            case 'T':
                v = 3;
                break;
            default:
                return 0;
        }
        code |= v << (2 * i);
    }
    *code_out = code;
    return 1;
}

static int dna2_code_local_fold(const char *s, size_t len, uint64_t *code_out) {
    if (s == NULL && len != 0) return 0;
    if (len > 32) return 0;
    uint64_t code = 0;
    for (size_t i = 0; i < len; ++i) {
        uint64_t v;
        switch (s[i]) {
            case 'A':
            case 'a':
                v = 0;
                break;
            case 'C':
            case 'c':
                v = 1;
                break;
            case 'G':
            case 'g':
                v = 2;
                break;
            case 'T':
            case 't':
                v = 3;
                break;
            default:
                return 0;
        }
        code |= v << (2 * i);
    }
    *code_out = code;
    return 1;
}

static int dna2_base_fold_value(char c, uint64_t *value_out) {
    switch (c) {
        case 'A':
        case 'a':
            *value_out = 0;
            return 1;
        case 'C':
        case 'c':
            *value_out = 1;
            return 1;
        case 'G':
        case 'g':
            *value_out = 2;
            return 1;
        case 'T':
        case 't':
            *value_out = 3;
            return 1;
        default:
            *value_out = 0;
            return 0;
    }
}

static void copy_upper_ascii_window(char *dst, size_t dst_cap, const char *src, size_t len) {
    if (dst_cap == 0) return;
    if (len >= dst_cap) len = dst_cap - 1;
    for (size_t i = 0; i < len; ++i) {
        unsigned char c = (unsigned char)src[i];
        dst[i] = (char)(c >= 'a' && c <= 'z' ? c - 32 : c);
    }
    dst[len] = '\0';
}

static char uppercase_ascii_char(char c) {
    unsigned char uc = (unsigned char)c;
    return (char)(uc >= 'a' && uc <= 'z' ? uc - 32 : uc);
}

static int phred33_quality(char c) {
    unsigned char uc = (unsigned char)c;
    return uc < 33 ? 0 : (int)uc - 33;
}

static int window_matches_observed_folded(const char *seq, size_t seq_len, size_t start,
                                          const char *observed, size_t observed_len) {
    if (start > seq_len || observed_len > seq_len - start) return 0;
    for (size_t i = 0; i < observed_len; ++i) {
        if (uppercase_ascii_char(seq[start + i]) != observed[i]) return 0;
    }
    return 1;
}

static int observed_base_qualities_within_edit(const char *observed, size_t observed_len,
                                               const char *target, size_t target_len,
                                               const char *qual, size_t qual_start,
                                               int max_correction_qual) {
    if (observed_len == target_len) {
        int saw_mismatch = 0;
        for (size_t i = 0; i < observed_len; ++i) {
            if (observed[i] == target[i]) continue;
            saw_mismatch = 1;
            if (phred33_quality(qual[qual_start + i]) > max_correction_qual) return 0;
        }
        return saw_mismatch;
    }

    if (observed_len == target_len + 1) {
        for (size_t drop = 0; drop < observed_len; ++drop) {
            size_t ti = 0;
            int matches = 1;
            for (size_t oi = 0; oi < observed_len; ++oi) {
                if (oi == drop) continue;
                if (observed[oi] != target[ti++]) {
                    matches = 0;
                    break;
                }
            }
            if (matches && phred33_quality(qual[qual_start + drop]) <= max_correction_qual) return 1;
        }
        return 0;
    }

    if (observed_len + 1 == target_len) {
        return 1;
    }

    return 0;
}

static int hamming_lookup_insert(hamming_lookup_entry *table, size_t cap, uint64_t code, int target_index) {
    size_t slot = code_hash_local(code, 0, cap);
    for (;;) {
        hamming_lookup_entry *entry = &table[slot];
        if (entry->target_index < 0) {
            entry->code = code;
            entry->target_index = target_index;
            entry->match_count = 1;
            return 0;
        }
        if (entry->code == code) {
            if (entry->target_index != target_index) {
                if (target_index < entry->target_index) entry->target_index = target_index;
                ++entry->match_count;
            }
            return 0;
        }
        slot = (slot + 1) & (cap - 1);
    }
}

static const hamming_lookup_entry *hamming_lookup_find(const hamming_lookup_entry *table, size_t cap, uint64_t code) {
    if (table == NULL || cap == 0) return NULL;
    size_t slot = code_hash_local(code, 0, cap);
    for (;;) {
        const hamming_lookup_entry *entry = &table[slot];
        if (entry->target_index < 0) return NULL;
        if (entry->code == code) return entry;
        slot = (slot + 1) & (cap - 1);
    }
}

static int levenshtein1_lookup_insert(hamming_lookup_entry *table, size_t cap, uint64_t code, size_t len,
                                      int target_index) {
    size_t slot = code_hash_local(code, len, cap);
    for (;;) {
        hamming_lookup_entry *entry = &table[slot];
        if (entry->target_index < 0) {
            entry->code = code;
            entry->target_index = target_index;
            entry->match_count = 1;
            return 0;
        }
        if (entry->code == code) {
            if (entry->target_index != target_index) {
                if (target_index < entry->target_index) entry->target_index = target_index;
                ++entry->match_count;
            }
            return 0;
        }
        slot = (slot + 1) & (cap - 1);
    }
}

static const hamming_lookup_entry *levenshtein1_lookup_find(const hamming_lookup_entry *table, size_t cap,
                                                           uint64_t code, size_t len) {
    if (table == NULL || cap == 0) return NULL;
    size_t slot = code_hash_local(code, len, cap);
    for (;;) {
        const hamming_lookup_entry *entry = &table[slot];
        if (entry->target_index < 0) return NULL;
        if (entry->code == code) return entry;
        slot = (slot + 1) & (cap - 1);
    }
}

static hamming_lookup_entry *alloc_hamming_table(size_t cap) {
    hamming_lookup_entry *table = (hamming_lookup_entry *)malloc(cap * sizeof(hamming_lookup_entry));
    if (table == NULL) return NULL;
    for (size_t i = 0; i < cap; ++i) {
        table[i].code = 0;
        table[i].target_index = -1;
        table[i].match_count = 0;
    }
    return table;
}

static int hamming_seed_insert(hamming_lookup *lookup, unsigned char seed_id, uint64_t code, int target_index) {
    if (lookup->n_seeds > (size_t)INT32_MAX) return -1;
    size_t seed_len = seed_id == 0 ? lookup->seed0_len : lookup->target_len - lookup->seed0_len;
    size_t slot = seed_hash_local(code, seed_len, seed_id, lookup->seed_hash_cap);
    size_t e = lookup->n_seeds++;
    lookup->seeds[e].code = code;
    lookup->seeds[e].target_index = target_index;
    lookup->seeds[e].seed_id = seed_id;
    lookup->seeds[e].next = lookup->seed_heads[slot];
    lookup->seed_heads[slot] = (int)e;
    return 0;
}

static void free_hamming_lookup(hamming_lookup *lookup) {
    if (lookup == NULL) return;
    free(lookup->exact);
    free(lookup->mismatch);
    free(lookup->seeds);
    free(lookup->seed_heads);
    free(lookup->target_codes);
    lookup->exact = NULL;
    lookup->mismatch = NULL;
    lookup->seeds = NULL;
    lookup->seed_heads = NULL;
    lookup->target_codes = NULL;
    lookup->exact_cap = 0;
    lookup->mismatch_cap = 0;
    lookup->seed_hash_cap = 0;
    lookup->n_seeds = 0;
    lookup->target_len = 0;
    lookup->seed0_len = 0;
    lookup->ready = 0;
    lookup->seed_ready = 0;
}

static void free_levenshtein1_lookup(levenshtein1_lookup *lookup) {
    if (lookup == NULL) return;
    free(lookup->exact);
    free(lookup->substitution);
    free(lookup->target_deletion);
    free(lookup->target_insertion);
    memset(lookup, 0, sizeof(*lookup));
}

static uint64_t code_remove_base_local(uint64_t code, size_t pos, size_t len) {
    uint64_t low = code & code_low_mask_local(pos);
    uint64_t high = code >> (2 * (pos + 1));
    (void)len;
    return low | (high << (2 * pos));
}

static uint64_t code_insert_base_local(uint64_t code, size_t pos, size_t len, uint64_t base) {
    uint64_t low = code & code_low_mask_local(pos);
    uint64_t high = (code >> (2 * pos)) & code_low_mask_local(len - pos);
    return low | (base << (2 * pos)) | (high << (2 * (pos + 1)));
}

static int build_levenshtein1_lookup(const seq_table *targets, size_t target_len, levenshtein1_lookup *lookup) {
    memset(lookup, 0, sizeof(*lookup));
    if (target_len == 0 || target_len > 31) return 0;
    for (size_t i = 0; i < targets->count; ++i) {
        uint64_t code = 0;
        if (targets->records[i].len != target_len ||
            !dna2_code_local(targets->records[i].seq, target_len, &code)) {
            return 0;
        }
    }

    lookup->exact_cap = next_pow2_local(targets->count * 2 + 16);
    lookup->substitution_cap = next_pow2_local(targets->count * target_len * 3 + 16);
    lookup->target_deletion_cap = next_pow2_local(targets->count * target_len + 16);
    lookup->target_insertion_cap = next_pow2_local(targets->count * (target_len + 1) * 4 + 16);
    lookup->exact = alloc_hamming_table(lookup->exact_cap);
    lookup->substitution = alloc_hamming_table(lookup->substitution_cap);
    lookup->target_deletion = alloc_hamming_table(lookup->target_deletion_cap);
    lookup->target_insertion = alloc_hamming_table(lookup->target_insertion_cap);
    if (lookup->exact == NULL || lookup->substitution == NULL || lookup->target_deletion == NULL ||
        lookup->target_insertion == NULL) {
        free_levenshtein1_lookup(lookup);
        return -1;
    }
    lookup->target_len = target_len;

    for (size_t i = 0; i < targets->count; ++i) {
        uint64_t code = 0;
        if (!dna2_code_local(targets->records[i].seq, target_len, &code)) {
            free_levenshtein1_lookup(lookup);
            return 0;
        }
        if (levenshtein1_lookup_insert(lookup->exact, lookup->exact_cap, code, target_len, (int)i) != 0) {
            free_levenshtein1_lookup(lookup);
            return -1;
        }
        for (size_t pos = 0; pos < target_len; ++pos) {
            uint64_t shift = (uint64_t)2 * pos;
            uint64_t old_base = (code >> shift) & 3ULL;
            uint64_t mask = 3ULL << shift;
            for (uint64_t b = 0; b < 4; ++b) {
                if (b == old_base) continue;
                uint64_t mutated = (code & ~mask) | (b << shift);
                if (levenshtein1_lookup_insert(lookup->substitution, lookup->substitution_cap, mutated,
                                               target_len, (int)i) != 0) {
                    free_levenshtein1_lookup(lookup);
                    return -1;
                }
            }
            uint64_t deleted = code_remove_base_local(code, pos, target_len);
            if (levenshtein1_lookup_insert(lookup->target_deletion, lookup->target_deletion_cap, deleted,
                                           target_len - 1, (int)i) != 0) {
                free_levenshtein1_lookup(lookup);
                return -1;
            }
        }
        for (size_t pos = 0; pos <= target_len; ++pos) {
            for (uint64_t b = 0; b < 4; ++b) {
                uint64_t inserted = code_insert_base_local(code, pos, target_len, b);
                if (levenshtein1_lookup_insert(lookup->target_insertion, lookup->target_insertion_cap, inserted,
                                               target_len + 1, (int)i) != 0) {
                    free_levenshtein1_lookup(lookup);
                    return -1;
                }
            }
        }
    }
    lookup->ready = 1;
    return 0;
}

static int build_hamming_lookup(const seq_table *targets, size_t target_len, hamming_lookup *lookup) {
    memset(lookup, 0, sizeof(*lookup));
    if (target_len == 0 || target_len > 32) return 0;
    for (size_t i = 0; i < targets->count; ++i) {
        if (targets->records[i].len != target_len) return 0;
        uint64_t code = 0;
        if (!dna2_code_local(targets->records[i].seq, target_len, &code)) return 0;
    }

    size_t exact_need = targets->count * 2 + 16;
    /* Tuned *3 (not *4) since exactly 3 mutations per position; smaller mismatch table
       improves cache residency for hamming_lookup precompute k=1 guide/barcode counting. */
    size_t mismatch_need = targets->count * target_len * 3 + 16;
    lookup->exact_cap = next_pow2_local(exact_need);
    lookup->mismatch_cap = next_pow2_local(mismatch_need);
    lookup->exact = alloc_hamming_table(lookup->exact_cap);
    lookup->mismatch = alloc_hamming_table(lookup->mismatch_cap);
    if (lookup->exact == NULL || lookup->mismatch == NULL) {
        free_hamming_lookup(lookup);
        return -1;
    }
    lookup->target_len = target_len;

    for (size_t i = 0; i < targets->count; ++i) {
        uint64_t code = 0;
        if (!dna2_code_local(targets->records[i].seq, target_len, &code)) {
            free_hamming_lookup(lookup);
            return 0;
        }
        if (hamming_lookup_insert(lookup->exact, lookup->exact_cap, code, (int)i) != 0) {
            free_hamming_lookup(lookup);
            return -1;
        }
        for (size_t pos = 0; pos < target_len; ++pos) {
            uint64_t shift = (uint64_t)2 * pos;
            uint64_t old_base = (code >> shift) & 3ULL;
            uint64_t mask = 3ULL << shift;
            for (uint64_t b = 0; b < 4; ++b) {
                if (b == old_base) continue;
                uint64_t mutated = (code & ~mask) | (b << shift);
                if (hamming_lookup_insert(lookup->mismatch, lookup->mismatch_cap, mutated, (int)i) != 0) {
                    free_hamming_lookup(lookup);
                    return -1;
                }
            }
        }
    }
    lookup->ready = 1;
    return 0;
}

static int build_hamming_exact_lookup(const seq_table *targets, size_t target_len, hamming_lookup *lookup) {
    memset(lookup, 0, sizeof(*lookup));
    if (target_len == 0 || target_len > 32) return 0;
    for (size_t i = 0; i < targets->count; ++i) {
        if (targets->records[i].len != target_len) return 0;
        uint64_t code = 0;
        if (!dna2_code_local(targets->records[i].seq, target_len, &code)) return 0;
    }

    size_t exact_need = targets->count * 2 + 16;
    lookup->exact_cap = next_pow2_local(exact_need);
    lookup->exact = alloc_hamming_table(lookup->exact_cap);
    if (lookup->exact == NULL) {
        free_hamming_lookup(lookup);
        return -1;
    }
    lookup->target_len = target_len;

    for (size_t i = 0; i < targets->count; ++i) {
        uint64_t code = 0;
        if (!dna2_code_local(targets->records[i].seq, target_len, &code)) {
            free_hamming_lookup(lookup);
            return 0;
        }
        if (hamming_lookup_insert(lookup->exact, lookup->exact_cap, code, (int)i) != 0) {
            free_hamming_lookup(lookup);
            return -1;
        }
    }
    lookup->ready = 1;
    return 0;
}

static int build_hamming_seed_lookup(const seq_table *targets, size_t target_len, hamming_lookup *lookup) {
    memset(lookup, 0, sizeof(*lookup));
    if (target_len < 2 || target_len > 32) return build_hamming_lookup(targets, target_len, lookup);
    for (size_t i = 0; i < targets->count; ++i) {
        if (targets->records[i].len != target_len) return 0;
        uint64_t code = 0;
        if (!dna2_code_local(targets->records[i].seq, target_len, &code)) return 0;
    }

    size_t exact_need = targets->count * 2 + 16;
    size_t seed_need = targets->count * 2 + 16;
    lookup->exact_cap = next_pow2_local(exact_need);
    lookup->seed_hash_cap = next_pow2_local(seed_need * 2 + 1);
    lookup->exact = alloc_hamming_table(lookup->exact_cap);
    lookup->seeds = (hamming_seed_entry *)malloc(seed_need * sizeof(hamming_seed_entry));
    lookup->seed_heads = (int *)malloc(lookup->seed_hash_cap * sizeof(int));
    lookup->target_codes = (uint64_t *)malloc((targets->count == 0 ? 1 : targets->count) * sizeof(uint64_t));
    if (lookup->exact == NULL || lookup->seeds == NULL || lookup->seed_heads == NULL || lookup->target_codes == NULL) {
        free_hamming_lookup(lookup);
        return -1;
    }
    for (size_t i = 0; i < lookup->seed_hash_cap; ++i) lookup->seed_heads[i] = -1;
    lookup->target_len = target_len;
    lookup->seed0_len = target_len / 2;

    for (size_t i = 0; i < targets->count; ++i) {
        uint64_t code = 0;
        if (!dna2_code_local(targets->records[i].seq, target_len, &code)) {
            free_hamming_lookup(lookup);
            return 0;
        }
        lookup->target_codes[i] = code;
        if (hamming_lookup_insert(lookup->exact, lookup->exact_cap, code, (int)i) != 0) {
            free_hamming_lookup(lookup);
            return -1;
        }
        uint64_t seed0 = code_segment_local(code, 0, lookup->seed0_len);
        uint64_t seed1 = code_segment_local(code, lookup->seed0_len, target_len - lookup->seed0_len);
        if (hamming_seed_insert(lookup, 0, seed0, (int)i) != 0 ||
            hamming_seed_insert(lookup, 1, seed1, (int)i) != 0) {
            free_hamming_lookup(lookup);
            return -1;
        }
    }
    lookup->ready = 1;
    lookup->seed_ready = 1;
    return 0;
}

static int cmp_ull_desc(const void *a, const void *b) {
    unsigned long long aa = *(const unsigned long long *)a;
    unsigned long long bb = *(const unsigned long long *)b;
    return aa < bb ? 1 : (aa > bb ? -1 : 0);
}

static int cmp_ull_asc(const void *a, const void *b) {
    unsigned long long aa = *(const unsigned long long *)a;
    unsigned long long bb = *(const unsigned long long *)b;
    return aa > bb ? 1 : (aa < bb ? -1 : 0);
}

static double gini_from_counts(const unsigned long long *values, size_t n) {
    if (n == 0) return 0.0;
    unsigned long long *tmp = (unsigned long long *)malloc(n * sizeof(unsigned long long));
    if (tmp == NULL) return 0.0;
    unsigned long long sum = 0;
    for (size_t i = 0; i < n; ++i) {
        tmp[i] = values[i];
        sum += values[i];
    }
    if (sum == 0) {
        free(tmp);
        return 0.0;
    }
    qsort(tmp, n, sizeof(unsigned long long), cmp_ull_asc);
    long double weighted = 0.0;
    for (size_t i = 0; i < n; ++i) weighted += (long double)(i + 1) * (long double)tmp[i];
    free(tmp);
    long double gini = (2.0L * weighted / ((long double)n * (long double)sum)) -
                       (((long double)n + 1.0L) / (long double)n);
    if (gini < 0.0L) return 0.0;
    if (gini > 1.0L) return 1.0;
    return (double)gini;
}

static double top_fraction_from_counts(const unsigned long long *values, size_t n, double fraction) {
    if (n == 0) return 0.0;
    unsigned long long *tmp = (unsigned long long *)malloc(n * sizeof(unsigned long long));
    if (tmp == NULL) return 0.0;
    unsigned long long sum = 0;
    for (size_t i = 0; i < n; ++i) {
        tmp[i] = values[i];
        sum += values[i];
    }
    if (sum == 0) {
        free(tmp);
        return 0.0;
    }
    qsort(tmp, n, sizeof(unsigned long long), cmp_ull_desc);
    size_t top_n = (size_t)((double)n * fraction);
    if (top_n == 0) top_n = 1;
    if (top_n > n) top_n = n;
    unsigned long long top_sum = 0;
    for (size_t i = 0; i < top_n; ++i) top_sum += tmp[i];
    free(tmp);
    return (double)top_sum / (double)sum;
}

static int compute_sample_qc_metrics(const seq_table *targets, const unsigned long long *counts, size_t sample_index,
                                     const count_stats *stats, sample_qc_metrics *metrics_out) {
    if (targets == NULL || counts == NULL || stats == NULL || metrics_out == NULL) return -1;
    unsigned long long *target_totals =
            (unsigned long long *)calloc(targets->count == 0 ? 1 : targets->count, sizeof(unsigned long long));
    if (target_totals == NULL) return -1;
    unsigned long long observed_targets = 0;
    for (size_t t = 0; t < targets->count; ++t) {
        for (size_t kind = 0; kind < 5; ++kind) {
            target_totals[t] += counts[((sample_index * targets->count + t) * 5) + kind];
        }
        if (target_totals[t] != 0) ++observed_targets;
    }
    unsigned long long valid = stats->total >= stats->invalid ? stats->total - stats->invalid : 0;
    double valid_denom = valid == 0 ? 1.0 : (double)valid;
    metrics_out->assignment_rate = (double)stats->unique / valid_denom;
    metrics_out->ambiguous_rate = (double)stats->ambiguous / valid_denom;
    metrics_out->no_match_rate = (double)stats->unmatched / valid_denom;
    metrics_out->invalid_rate = stats->total == 0 ? 0.0 : (double)stats->invalid / (double)stats->total;
    metrics_out->coverage_fraction =
            targets->count == 0 ? 0.0 : (double)observed_targets / (double)targets->count;
    metrics_out->zero_count_fraction =
            targets->count == 0 ? 0.0 : (double)(targets->count - observed_targets) / (double)targets->count;
    metrics_out->gini_index = gini_from_counts(target_totals, targets->count);
    metrics_out->top_1pct_fraction = top_fraction_from_counts(target_totals, targets->count, 0.01);
    free(target_totals);
    return 0;
}

static void emit_sample_qc_review_warnings(const string_list *labels, const sample_qc_metrics *metrics, size_t count) {
    if (labels == NULL || metrics == NULL || count == 0) return;
    int any = 0;
    for (size_t sample = 0; sample < count; ++sample) {
        const sample_qc_metrics *m = &metrics[sample];
        int sample_warn = 0;
        if (m->assignment_rate < 0.80) sample_warn = 1;
        if (m->ambiguous_rate > 0.05) sample_warn = 1;
        if (m->no_match_rate > 0.15) sample_warn = 1;
        if (m->invalid_rate > 0.02) sample_warn = 1;
        if (m->coverage_fraction < 0.90) sample_warn = 1;
        if (m->zero_count_fraction > 0.10) sample_warn = 1;
        if (m->gini_index > 0.50) sample_warn = 1;
        if (m->top_1pct_fraction > 0.30) sample_warn = 1;
        if (!sample_warn) continue;
        any = 1;
        fprintf(stderr, "dotmatch: QC review recommended for sample %s:", labels->items[sample]);
        if (m->assignment_rate < 0.80) fprintf(stderr, " assignment_rate=%.1f%%", 100.0 * m->assignment_rate);
        if (m->ambiguous_rate > 0.05) fprintf(stderr, " ambiguous_rate=%.1f%%", 100.0 * m->ambiguous_rate);
        if (m->no_match_rate > 0.15) fprintf(stderr, " no_match_rate=%.1f%%", 100.0 * m->no_match_rate);
        if (m->invalid_rate > 0.02) fprintf(stderr, " invalid_rate=%.1f%%", 100.0 * m->invalid_rate);
        if (m->coverage_fraction < 0.90) fprintf(stderr, " coverage=%.1f%%", 100.0 * m->coverage_fraction);
        if (m->zero_count_fraction > 0.10) fprintf(stderr, " zero_count_guides=%.1f%%", 100.0 * m->zero_count_fraction);
        if (m->gini_index > 0.50) fprintf(stderr, " gini=%.2f", m->gini_index);
        if (m->top_1pct_fraction > 0.30) fprintf(stderr, " top_1pct_fraction=%.1f%%", 100.0 * m->top_1pct_fraction);
        fprintf(stderr, "\n");
    }
    if (any) {
        fprintf(stderr,
                "dotmatch: review sample_qc.tsv and summary.json before downstream MAGeCK/BAGEL analysis; "
                "thresholds are conservative diagnostics, not biological pass/fail rules\n");
    }
}

typedef enum ambiguity_policy {
    AMBIGUITY_POLICY_BEST = 0,
    AMBIGUITY_POLICY_RADIUS = 1
} ambiguity_policy;

static const char *ambiguity_policy_name(ambiguity_policy policy) {
    return policy == AMBIGUITY_POLICY_RADIUS ? "radius" : "best";
}

static int apply_ambiguity_policy(qdaln_match_result *result, ambiguity_policy policy) {
    if (policy == AMBIGUITY_POLICY_RADIUS && result->status == QDALN_MATCH_UNIQUE && result->match_count > 1) {
        result->status = QDALN_MATCH_AMBIGUOUS;
    }
    return 0;
}

static void html_escape(FILE *out, const char *s);

static void write_tsv_preview_table(FILE *out, const char *title, const char *path, size_t max_rows) {
    FILE *in = fopen(path, "r");
    if (in == NULL) return;
    fprintf(out, "<h2>");
    html_escape(out, title);
    fprintf(out, "</h2><table>\n");
    char line[16384];
    size_t row = 0;
    while (row <= max_rows && fgets(line, sizeof(line), in) != NULL) {
        trim_line(line);
        fprintf(out, "<tr>");
        char *fields[128];
        size_t n = split_fields(line, '\t', fields, 128);
        for (size_t i = 0; i < n; ++i) {
            fprintf(out, row == 0 ? "<th>" : "<td>");
            html_escape(out, fields[i]);
            fprintf(out, row == 0 ? "</th>" : "</td>");
        }
        fprintf(out, "</tr>\n");
        ++row;
    }
    fprintf(out, "</table>\n");
    fclose(in);
}

static int write_count_html_report(const char *path, const seq_table *targets, const string_list *reads,
                                   const string_list *labels, const unsigned long long *counts,
                                   const count_stats *stats_by_sample, const offset_list *selected_offsets,
                                   int k, count_metric metric, ambiguity_policy policy, size_t target_len,
                                   const char *audit_dir, const char *unmatched_report_path) {
    FILE *out = open_output_file(path);
    if (out == NULL) return -1;

    int needs_review = 0;
    for (size_t sample = 0; sample < reads->count; ++sample) {
        const count_stats *s = &stats_by_sample[sample];
        unsigned long long valid = s->total >= s->invalid ? s->total - s->invalid : 0;
        double denom = valid == 0 ? 1.0 : (double)valid;
        if ((double)s->ambiguous / denom > 0.01 || (double)s->unmatched / denom > 0.10) needs_review = 1;
    }

    fprintf(out,
            "<!doctype html>\n<html><head><meta charset=\"utf-8\"><title>DotMatch Report</title>"
            "<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:0;line-height:1.45;color:#17202a;background:#f7f9fb}"
            "main{max-width:1160px;margin:0 auto;padding:32px}h1{font-size:32px;margin:0 0 8px}h2{margin-top:28px;border-bottom:1px solid #d8dee4;padding-bottom:6px}"
            ".lede{color:#57606a;margin:0 0 20px}.metric{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}"
            ".metric div,.note{background:#fff;border:1px solid #d8dee4;padding:12px;border-radius:8px}.label{font-size:12px;color:#57606a;text-transform:uppercase;letter-spacing:.04em}"
            ".value{font-size:20px;font-weight:650;margin-top:4px}table{border-collapse:collapse;width:100%%;margin:16px 0;background:#fff}"
            "th,td{border:1px solid #d8dee4;padding:7px 9px;text-align:right;vertical-align:top}th:first-child,td:first-child{text-align:left}th{background:#eef2f7}"
            ".warn{color:#9a6700}.ok{color:#1a7f37}.bad{color:#cf222e}</style></head><body><main>\n");
    fprintf(out, "<h1>DotMatch Report</h1>\n");
    fprintf(out, "<p class=\"lede\">Known-target assignment report for %zu target%s and %zu sample%s. Ambiguous reads are not silently counted.</p>\n",
            targets->count, targets->count == 1 ? "" : "s", reads->count, reads->count == 1 ? "" : "s");
    fprintf(out, "<h2>Run Status</h2><div class=\"metric\">"
                 "<div><span class=\"label\">Status</span><div class=\"value %s\">%s</div></div>"
                 "<div><span class=\"label\">Targets</span><div class=\"value\">%zu</div></div>"
                 "<div><span class=\"label\">Samples</span><div class=\"value\">%zu</div></div>"
                 "<div><span class=\"label\">Target length</span><div class=\"value\">%zu</div></div></div>\n",
            needs_review ? "warn" : "ok", needs_review ? "Needs Review" : "Ready",
            targets->count, reads->count, target_len);
    fprintf(out, "<h2>Inputs and Configuration</h2><div class=\"metric\"><div><span class=\"label\">k</span><div class=\"value\">%d</div></div>"
                 "<div><span class=\"label\">Metric</span><div class=\"value\">%s</div></div>"
                 "<div><span class=\"label\">Ambiguity policy</span><div class=\"value\">%s</div></div>"
                 "<div><span class=\"label\">Assignment</span><div class=\"value\">Known target</div></div></div>\n",
            k, metric_name(metric), ambiguity_policy_name(policy));

    fprintf(out, "<h2>Target Assignment QC</h2><table><tr><th>Sample</th><th>Total reads</th><th>Valid windows</th><th>Assignment rate</th>"
                 "<th>Exact rate</th><th>Rescue rate</th><th>Ambiguous rate</th><th>No-match rate</th>"
                 "<th>Library coverage</th><th>Candidates verified</th></tr>\n");
    for (size_t sample = 0; sample < reads->count; ++sample) {
        const count_stats *s = &stats_by_sample[sample];
        unsigned long long covered = 0;
        for (size_t t = 0; t < targets->count; ++t) {
            unsigned long long total = 0;
            for (size_t kind = 0; kind < 5; ++kind) total += counts[((sample * targets->count + t) * 5) + kind];
            if (total != 0) ++covered;
        }
        unsigned long long valid = s->total >= s->invalid ? s->total - s->invalid : 0;
        double denom = valid == 0 ? 1.0 : (double)valid;
        fprintf(out, "<tr><td>");
        html_escape(out, labels->items[sample]);
        fprintf(out, "</td><td>%llu</td><td>%llu</td><td>%.2f%%</td><td>%.2f%%</td><td>%.2f%%</td><td>%.2f%%</td>"
                     "<td>%.2f%%</td><td>%.2f%%</td><td>%llu</td></tr>\n",
                s->total, valid, 100.0 * (double)s->unique / denom, 100.0 * (double)s->exact / denom,
                100.0 * (double)s->corrected / denom, 100.0 * (double)s->ambiguous / denom,
                100.0 * (double)s->unmatched / denom,
                targets->count == 0 ? 0.0 : 100.0 * (double)covered / (double)targets->count,
                s->candidates_verified);
    }
    fprintf(out, "</table>\n");

    fprintf(out, "<h2>Warnings</h2><ul>\n");
    if (!needs_review) {
        fprintf(out, "<li class=\"ok\">No high ambiguous or no-match warning thresholds were crossed.</li>\n");
    }
    for (size_t sample = 0; sample < reads->count; ++sample) {
        const count_stats *s = &stats_by_sample[sample];
        unsigned long long valid = s->total >= s->invalid ? s->total - s->invalid : 0;
        double denom = valid == 0 ? 1.0 : (double)valid;
        if ((double)s->ambiguous / denom > 0.01) {
            fprintf(out, "<li class=\"warn\">Sample ");
            html_escape(out, labels->items[sample]);
            fprintf(out, " has ambiguous assignments above 1%% of valid extracted windows.</li>\n");
        }
        if ((double)s->unmatched / denom > 0.10) {
            fprintf(out, "<li class=\"warn\">Sample ");
            html_escape(out, labels->items[sample]);
            fprintf(out, " has no-match reads above 10%% of valid extracted windows.</li>\n");
        }
    }
    fprintf(out, "<li class=\"ok\">Ambiguous reads are not silently counted.</li></ul>\n");

    fprintf(out, "<h2>Input Files</h2><table><tr><th>Sample</th><th>FASTQ</th><th>Selected start(s)</th></tr>\n");
    for (size_t sample = 0; sample < reads->count; ++sample) {
        fprintf(out, "<tr><td>");
        html_escape(out, labels->items[sample]);
        fprintf(out, "</td><td>");
        html_escape(out, reads->items[sample]);
        fprintf(out, "</td><td>");
        for (size_t i = 0; i < selected_offsets[sample].count; ++i) {
            if (i != 0) fprintf(out, ", ");
            fprintf(out, "%zu", selected_offsets[sample].items[i]);
        }
        fprintf(out, "</td></tr>\n");
    }
    fprintf(out, "</table>\n");

    if (audit_dir != NULL) {
        char audit_path[4096];
        int n = snprintf(audit_path, sizeof(audit_path), "%s/%s", audit_dir, "audit_summary.tsv");
        if (n >= 0 && (size_t)n < sizeof(audit_path)) {
            write_tsv_preview_table(out, "Library Audit", audit_path, 40);
        }
    }
    if (unmatched_report_path != NULL) {
        write_tsv_preview_table(out, "Top Unmatched", unmatched_report_path, 25);
    }

    fprintf(out, "</main></body></html>\n");
    fclose(out);
    return 0;
}

typedef struct count_dirty_slot {
    size_t slot;
    unsigned long long count;
} count_dirty_slot;

typedef struct count_dirty_slots {
    count_dirty_slot *items;
    size_t count;
    size_t cap;
    size_t *table;
    size_t table_cap;
} count_dirty_slots;

typedef struct count_sample_job {
    const qdaln_index *index;
    const hamming_lookup *hlookup;
    const levenshtein1_lookup *levlookup;
    const seq_table *targets;
    const char **target_ptrs;
    const size_t *target_lens;
    const char *reads_path;
    const char *sample_label;
    size_t sample_index;
    offset_list *selected_offsets;
    size_t target_len;
    int k;
    count_metric metric;
    size_t indel_window;
    unsigned long long *counts;
    count_stats *stats;
    FILE *assignments;
    FILE *ambiguous_out;
    FILE *unmatched_out;
    const char *ambiguous_policy;
    ambiguity_policy assignment_policy;
    int direct_hamming_counts;
    int metal_hamming_counts;
    const uint64_t *metal_target_codes;
    int fused_offset_detection;
    size_t target_start;
    size_t auto_offset;
    size_t auto_offset_sample;
    offset_mode offsets_mode;
    double offset_min_fraction;
    size_t read_threads;
    int max_correction_qual;
    int rc;
    count_dirty_slots *dirty_slots;
    count_progress *progress;
} count_sample_job;

static void write_assignment_like_row(FILE *out, const seq_table *targets, const char *sample, const char *read_id,
                                      const char *observed, qdaln_match_result r, const char *correction) {
    const char *target_id = "";
    const char *target_seq = "";
    if (r.target_index >= 0) {
        target_id = targets->records[r.target_index].id;
        target_seq = targets->records[r.target_index].seq;
    }
    fprintf(out, "%s\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n",
            sample, read_id, observed, r.target_index, target_id, target_seq, r.best_distance,
            r.second_best_distance, r.match_count, status_name(r.status), correction);
}

static int find_observed_quality_window(const char *seq, size_t seq_len, const offset_list *offsets,
                                        size_t fallback_offset, size_t target_len, count_metric metric,
                                        size_t indel_window, int k, const char *observed,
                                        size_t observed_len, size_t *start_out) {
    size_t min_len = target_len;
    size_t max_len = target_len;
    if (metric == COUNT_METRIC_LEVENSHTEIN && indel_window != 0 && k == 1) {
        min_len = target_len > indel_window ? target_len - indel_window : 0;
        max_len = target_len + indel_window;
    }
    if (observed_len < min_len || observed_len > max_len) return 0;

    size_t n_offsets = offsets == NULL || offsets->count == 0 ? 1 : offsets->count;
    for (size_t i = 0; i < n_offsets; ++i) {
        size_t offset = offsets == NULL || offsets->count == 0 ? fallback_offset : offsets->items[i];
        if (window_matches_observed_folded(seq, seq_len, offset, observed, observed_len)) {
            *start_out = offset;
            return 1;
        }
    }
    return 0;
}

static int quality_allows_unique_correction(const char *seq, size_t seq_len, const char *qual,
                                            const offset_list *offsets, size_t fallback_offset,
                                            size_t target_len, count_metric metric, size_t indel_window,
                                            int k, const char *observed, const seq_record *target,
                                            qdaln_match_result result, int max_correction_qual) {
    if (max_correction_qual < 0 || qual == NULL) return 1;
    if (result.status != QDALN_MATCH_UNIQUE || result.best_distance <= 0 || result.target_index < 0) return 1;

    size_t observed_len = strlen(observed);
    size_t qual_start = 0;
    if (!find_observed_quality_window(seq, seq_len, offsets, fallback_offset, target_len, metric, indel_window,
                                      k, observed, observed_len, &qual_start)) {
        return 0;
    }
    return observed_base_qualities_within_edit(observed, observed_len, target->seq, target->len, qual,
                                               qual_start, max_correction_qual);
}

static void html_escape(FILE *out, const char *s) {
    for (; s != NULL && *s != '\0'; ++s) {
        switch (*s) {
            case '&':
                fputs("&amp;", out);
                break;
            case '<':
                fputs("&lt;", out);
                break;
            case '>':
                fputs("&gt;", out);
                break;
            case '"':
                fputs("&quot;", out);
                break;
            default:
                fputc(*s, out);
                break;
        }
    }
}

static int build_target_arrays(const seq_table *targets, const char ***target_ptrs_out, size_t **target_lens_out) {
    const char **target_ptrs = (const char **)malloc(targets->count * sizeof(char *));
    size_t *target_lens = (size_t *)malloc(targets->count * sizeof(size_t));
    if (targets->count != 0 && (target_ptrs == NULL || target_lens == NULL)) {
        free(target_ptrs);
        free(target_lens);
        return -1;
    }
    for (size_t i = 0; i < targets->count; ++i) {
        target_ptrs[i] = targets->records[i].seq;
        target_lens[i] = targets->records[i].len;
    }
    *target_ptrs_out = target_ptrs;
    *target_lens_out = target_lens;
    return 0;
}

static int all_targets_have_length(const seq_table *targets, size_t len) {
    for (size_t i = 0; i < targets->count; ++i) {
        if (targets->records[i].len != len) return 0;
    }
    return 1;
}

static int cmp_size_asc(const void *a, const void *b) {
    size_t aa = *(const size_t *)a;
    size_t bb = *(const size_t *)b;
    return aa > bb ? 1 : (aa < bb ? -1 : 0);
}

static int collect_target_lengths(const seq_table *targets, size_t **lengths_out, size_t *count_out) {
    size_t *lengths = (size_t *)malloc((targets->count == 0 ? 1 : targets->count) * sizeof(size_t));
    if (lengths == NULL) return -1;
    size_t count = 0;
    for (size_t i = 0; i < targets->count; ++i) {
        size_t len = targets->records[i].len;
        int seen = 0;
        for (size_t j = 0; j < count; ++j) {
            if (lengths[j] == len) {
                seen = 1;
                break;
            }
        }
        if (!seen) lengths[count++] = len;
    }
    qsort(lengths, count, sizeof(size_t), cmp_size_asc);
    *lengths_out = lengths;
    *count_out = count;
    return 0;
}

typedef struct match_merge_hit {
    int target_index;
    int distance;
} match_merge_hit;

typedef struct match_merge_state {
    match_merge_hit inline_hits[16];
    match_merge_hit *hits;
    size_t count;
    size_t cap;
    int next_synthetic_target;
    int saw_none;
} match_merge_state;

static void merge_state_init(match_merge_state *state) {
    state->hits = state->inline_hits;
    state->count = 0;
    state->cap = sizeof(state->inline_hits) / sizeof(state->inline_hits[0]);
    state->next_synthetic_target = -2;
    state->saw_none = 0;
}

static void merge_state_free(match_merge_state *state) {
    if (state->hits != state->inline_hits) free(state->hits);
    merge_state_init(state);
}

static void copy_merge_observed(char *dst, size_t dst_cap, const char *src) {
    if (dst_cap == 0) return;
    size_t n = 0;
    while (n + 1 < dst_cap && src[n] != '\0') ++n;
    memcpy(dst, src, n);
    dst[n] = '\0';
}

static int merge_state_grow(match_merge_state *state) {
    size_t next_cap = state->cap * 2;
    match_merge_hit *next = (match_merge_hit *)malloc(next_cap * sizeof(match_merge_hit));
    if (next == NULL) return -1;
    memcpy(next, state->hits, state->count * sizeof(match_merge_hit));
    if (state->hits != state->inline_hits) free(state->hits);
    state->hits = next;
    state->cap = next_cap;
    return 0;
}

static int merge_state_add_hit(match_merge_state *state, int target_index, int distance) {
    if (target_index >= 0) {
        for (size_t i = 0; i < state->count; ++i) {
            if (state->hits[i].target_index == target_index) {
                if (distance < state->hits[i].distance) state->hits[i].distance = distance;
                return 0;
            }
        }
    }
    if (state->count == state->cap && merge_state_grow(state) != 0) return -1;
    state->hits[state->count].target_index = target_index;
    state->hits[state->count].distance = distance;
    ++state->count;
    return 0;
}

static int merge_state_best_distance(const match_merge_state *state) {
    int best = -1;
    for (size_t i = 0; i < state->count; ++i) {
        int d = state->hits[i].distance;
        if (best < 0 || d < best) best = d;
    }
    return best;
}

static int merge_state_add_result(match_merge_state *state, char *best_observed, size_t best_observed_cap,
                                  const char *observed, qdaln_match_result r) {
    if (r.status == QDALN_MATCH_INVALID) return 0;
    if (r.match_count == 0) {
        if (state->count == 0 && !state->saw_none && best_observed_cap != 0) {
            copy_merge_observed(best_observed, best_observed_cap, observed);
        }
        state->saw_none = 1;
        return 0;
    }

    int prior_best = merge_state_best_distance(state);
    if ((prior_best < 0 || r.best_distance < prior_best) && best_observed_cap != 0) {
        copy_merge_observed(best_observed, best_observed_cap, observed);
    }

    if (r.target_index >= 0 && merge_state_add_hit(state, r.target_index, r.best_distance) != 0) return -1;
    int remaining = r.match_count - (r.target_index >= 0 ? 1 : 0);
    if (remaining <= 0) return 0;

    int synthetic_distance = r.best_distance;
    if (r.status != QDALN_MATCH_AMBIGUOUS && r.second_best_distance >= 0) {
        synthetic_distance = r.second_best_distance;
    }
    for (int i = 0; i < remaining; ++i) {
        if (merge_state_add_hit(state, state->next_synthetic_target--, synthetic_distance) != 0) return -1;
    }
    return 0;
}

static void merge_state_finish(const match_merge_state *state, qdaln_match_result *result) {
    *result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_INVALID};
    if (state->count == 0) {
        result->status = state->saw_none ? QDALN_MATCH_NONE : QDALN_MATCH_INVALID;
        return;
    }

    int best_ties = 0;
    for (size_t i = 0; i < state->count; ++i) {
        int d = state->hits[i].distance;
        ++result->match_count;
        if (result->best_distance < 0 || d < result->best_distance) {
            result->second_best_distance = result->best_distance;
            result->best_distance = d;
            result->target_index = state->hits[i].target_index >= 0 ? state->hits[i].target_index : -1;
            best_ties = 1;
        } else if (d == result->best_distance) {
            if (state->hits[i].target_index >= 0 &&
                (result->target_index < 0 || state->hits[i].target_index < result->target_index)) {
                result->target_index = state->hits[i].target_index;
            }
            ++best_ties;
        } else if (result->second_best_distance < 0 || d < result->second_best_distance) {
            result->second_best_distance = d;
        }
    }
    result->status = best_ties > 1 ? QDALN_MATCH_AMBIGUOUS : QDALN_MATCH_UNIQUE;
}

static void merge_summary_result(qdaln_match_result *best, char *best_observed, size_t best_observed_cap,
                                 const char *observed, qdaln_match_result r) {
    if (r.status == QDALN_MATCH_INVALID) return;
    if (r.match_count == 0) {
        if (best->status == QDALN_MATCH_INVALID) {
            *best = r;
            copy_merge_observed(best_observed, best_observed_cap, observed);
        }
        return;
    }
    if (best->match_count == 0 || best->best_distance < 0 || r.best_distance < best->best_distance) {
        *best = r;
        copy_merge_observed(best_observed, best_observed_cap, observed);
        return;
    }
    if (r.best_distance == best->best_distance) {
        if (r.target_index != best->target_index || r.status == QDALN_MATCH_AMBIGUOUS ||
            best->status == QDALN_MATCH_AMBIGUOUS) {
            if (r.target_index >= 0 && (best->target_index < 0 || r.target_index < best->target_index)) {
                best->target_index = r.target_index;
            }
            best->status = QDALN_MATCH_AMBIGUOUS;
            best->match_count += r.match_count;
        }
    } else if (best->second_best_distance < 0 || r.best_distance < best->second_best_distance) {
        best->second_best_distance = r.best_distance;
        best->match_count += r.match_count;
    }
}

static int assign_count_window(const qdaln_index *index, const char *seq, size_t seq_len, size_t target_start,
                               size_t target_len, int k, count_metric metric, size_t indel_window,
                               qdaln_match_result *result, qdaln_index_stats *stats, char *observed,
                               size_t observed_cap, int best_exact_shortcut) {
    *result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_INVALID};
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }
    if (observed_cap != 0) observed[0] = '\0';

    match_merge_state merge;
    merge_state_init(&merge);
    int rc = 0;
    size_t min_len = target_len;
    size_t max_len = target_len;
    if (metric == COUNT_METRIC_LEVENSHTEIN && indel_window != 0 && k == 1) {
        min_len = target_len > indel_window ? target_len - indel_window : 0;
        max_len = target_len + indel_window;
    }

    for (size_t len = min_len; len <= max_len; ++len) {
        if (len >= observed_cap) continue;
        if (target_start > seq_len || len > seq_len - target_start) continue;

        if (metric == COUNT_METRIC_HAMMING && k == 0) {
            qdaln_match_result r;
            qdaln_index_stats s = {0, 0};
            int exact_rc = qdaln_index_lookup_exact_ascii_stats(index, seq + target_start, len, &r, &s);
            if (exact_rc != 0) {
                rc = -1;
                goto done;
            }
            if (stats != NULL) {
                stats->candidates_considered += s.candidates_considered;
                stats->candidates_verified += s.candidates_verified;
            }
            char candidate[8192];
            if (len >= sizeof(candidate)) continue;
            memcpy(candidate, seq + target_start, len);
            candidate[len] = '\0';
            uppercase_ascii(candidate);
            if (merge_state_add_result(&merge, observed, observed_cap, candidate, r) != 0) {
                rc = -1;
                goto done;
            }
            continue;
        }

        char candidate[8192];
        if (len >= sizeof(candidate)) continue;
        memcpy(candidate, seq + target_start, len);
        candidate[len] = '\0';
        uppercase_ascii(candidate);

        const char *read_ptr = candidate;
        size_t read_len = len;
        qdaln_match_result r;
        qdaln_index_stats s = {0, 0};
        if (best_exact_shortcut && metric == COUNT_METRIC_HAMMING && k == 1) {
            int exact_rc = qdaln_index_lookup_exact_stats(index, read_ptr, read_len, &r, &s);
            if (exact_rc != 0) {
                rc = -1;
                goto done;
            }
            if (stats != NULL) {
                stats->candidates_considered += s.candidates_considered;
                stats->candidates_verified += s.candidates_verified;
            }
            if (r.status == QDALN_MATCH_UNIQUE || r.status == QDALN_MATCH_AMBIGUOUS || r.status == QDALN_MATCH_INVALID) {
                if (merge_state_add_result(&merge, observed, observed_cap, candidate, r) != 0) {
                    rc = -1;
                    goto done;
                }
                continue;
            }
            s.candidates_considered = 0;
            s.candidates_verified = 0;
        }
        int assign_rc = metric == COUNT_METRIC_HAMMING
                ? qdaln_index_assign_hamming_stats(index, &read_ptr, &read_len, 1, k, &r, &s)
                : best_exact_shortcut
                        ? qdaln_index_assign_status_stats(index, &read_ptr, &read_len, 1, k, &r, &s)
                        : qdaln_index_assign_stats(index, &read_ptr, &read_len, 1, k, &r, &s);
        if (assign_rc != 0) {
            rc = -1;
            goto done;
        }
        if (stats != NULL) {
            stats->candidates_considered += s.candidates_considered;
            stats->candidates_verified += s.candidates_verified;
        }
        if (merge_state_add_result(&merge, observed, observed_cap, candidate, r) != 0) {
            rc = -1;
            goto done;
        }
    }

done:
    merge_state_finish(&merge, result);
    merge_state_free(&merge);
    return rc;
}

static int assign_count_length_set(const qdaln_index *index, const char *seq, size_t seq_len, size_t target_start,
                                   const size_t *lengths, size_t n_lengths, int k, count_metric metric,
                                   size_t indel_window, qdaln_match_result *result, qdaln_index_stats *stats,
                                   char *observed, size_t observed_cap, int best_exact_shortcut) {
    *result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_INVALID};
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }
    if (observed_cap != 0) observed[0] = '\0';
    for (size_t i = 0; i < n_lengths; ++i) {
        qdaln_match_result r = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
        qdaln_index_stats s = {0, 0};
        char candidate[8192];
        if (assign_count_window(index, seq, seq_len, target_start, lengths[i], k, metric, indel_window,
                                &r, &s, candidate, sizeof(candidate), best_exact_shortcut) != 0) {
            return -1;
        }
        if (stats != NULL) {
            stats->candidates_considered += s.candidates_considered;
            stats->candidates_verified += s.candidates_verified;
        }
        merge_summary_result(result, observed, observed_cap, candidate, r);
    }
    return 0;
}

static int assign_count_offsets(const qdaln_index *index, const char *seq, size_t seq_len,
                                const offset_list *offsets, size_t fallback_offset, size_t target_len,
                                int k, count_metric metric, size_t indel_window,
                                qdaln_match_result *result, qdaln_index_stats *stats,
                                char *observed, size_t observed_cap, int best_exact_shortcut) {
    *result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_INVALID};
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }
    if (observed_cap != 0) observed[0] = '\0';

    size_t n_offsets = offsets == NULL || offsets->count == 0 ? 1 : offsets->count;
    if (best_exact_shortcut && k == 1 && target_len < observed_cap) {
        match_merge_state exact_merge;
        merge_state_init(&exact_merge);
        qdaln_match_result exact_result = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
        qdaln_index_stats exact_stats_total = {0, 0};
        char exact_observed[8192];
        int exact_rc_total = 0;
        for (size_t i = 0; i < n_offsets; ++i) {
            size_t offset = offsets == NULL || offsets->count == 0 ? fallback_offset : offsets->items[i];
            if (offset > seq_len || target_len > seq_len - offset || target_len >= sizeof(exact_observed)) continue;
            qdaln_match_result exact_one;
            qdaln_index_stats exact_stats = {0, 0};
            if (qdaln_index_lookup_exact_ascii_stats(index, seq + offset, target_len, &exact_one, &exact_stats) != 0) {
                merge_state_free(&exact_merge);
                return -1;
            }
            exact_stats_total.candidates_considered += exact_stats.candidates_considered;
            exact_stats_total.candidates_verified += exact_stats.candidates_verified;
            if (exact_one.match_count == 0) continue;
            memcpy(exact_observed, seq + offset, target_len);
            exact_observed[target_len] = '\0';
            uppercase_ascii(exact_observed);
            if (merge_state_add_result(&exact_merge, observed, observed_cap, exact_observed, exact_one) != 0) {
                exact_rc_total = -1;
                break;
            }
        }
        merge_state_finish(&exact_merge, &exact_result);
        merge_state_free(&exact_merge);
        if (exact_rc_total != 0) return -1;
        if (stats != NULL) {
            stats->candidates_considered += exact_stats_total.candidates_considered;
            stats->candidates_verified += exact_stats_total.candidates_verified;
        }
        if (exact_result.match_count > 0) {
            *result = exact_result;
            return 0;
        }
    }
    match_merge_state merge;
    merge_state_init(&merge);
    int rc = 0;
    for (size_t i = 0; i < n_offsets; ++i) {
        size_t offset = offsets == NULL || offsets->count == 0 ? fallback_offset : offsets->items[i];
        qdaln_match_result candidate = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
        qdaln_index_stats local_stats = {0, 0};
        char local_observed[8192];
        if (assign_count_window(index, seq, seq_len, offset, target_len, k, metric, indel_window,
                                &candidate, &local_stats, local_observed, sizeof(local_observed),
                                best_exact_shortcut) != 0) {
            rc = -1;
            break;
        }
        if (stats != NULL) {
            stats->candidates_considered += local_stats.candidates_considered;
            stats->candidates_verified += local_stats.candidates_verified;
        }
        if (merge_state_add_result(&merge, observed, observed_cap, local_observed, candidate) != 0) {
            rc = -1;
            break;
        }
    }

    merge_state_finish(&merge, result);
    merge_state_free(&merge);
    return rc;
}

static int hamming_distance_within_k_cli(const char *a, size_t a_len, const char *b, size_t b_len, int k) {
    if (a_len != b_len) return -1;
    int d = 0;
    for (size_t i = 0; i < a_len; ++i) {
        if (a[i] != b[i] && ++d > k) return -1;
    }
    return d;
}

static int hamming_distance_cli(const char *a, size_t a_len, const char *b, size_t b_len) {
    if (a_len != b_len) return -1;
    int d = 0;
    for (size_t i = 0; i < a_len; ++i) {
        if (a[i] != b[i]) ++d;
    }
    return d;
}

static int scan_assign_metric(const char *read, size_t read_len, const char *const *targets,
                              const size_t *target_lens, size_t n_targets, int k,
                              count_metric metric, qdaln_match_result *result) {
    *result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_NONE};
    int best_ties = 0;
    for (size_t i = 0; i < n_targets; ++i) {
        int d = metric == COUNT_METRIC_HAMMING
                ? hamming_distance_within_k_cli(read, read_len, targets[i], target_lens[i], k)
                : qdaln_edit_distance_leq(read, read_len, targets[i], target_lens[i], k) > 0
                        ? qdaln_edit_distance(read, read_len, targets[i], target_lens[i])
                        : -1;
        if (d < 0 || d > k) continue;
        ++result->match_count;
        if (result->best_distance < 0 || d < result->best_distance) {
            result->second_best_distance = result->best_distance;
            result->best_distance = d;
            result->target_index = (int)i;
            best_ties = 1;
        } else if (d == result->best_distance) {
            if (result->target_index < 0 || (int)i < result->target_index) result->target_index = (int)i;
            ++best_ties;
        } else if (result->second_best_distance < 0 || d < result->second_best_distance) {
            result->second_best_distance = d;
        }
    }
    if (result->match_count == 0) {
        result->status = QDALN_MATCH_NONE;
    } else if (best_ties > 1) {
        result->status = QDALN_MATCH_AMBIGUOUS;
    } else {
        result->status = QDALN_MATCH_UNIQUE;
    }
    return 0;
}

static int scan_count_window(const char *const *targets, const size_t *target_lens, size_t n_targets,
                             const char *seq, size_t seq_len, size_t target_start, size_t target_len,
                             int k, count_metric metric, size_t indel_window,
                             qdaln_match_result *result, char *observed, size_t observed_cap) {
    *result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_INVALID};
    if (observed_cap != 0) observed[0] = '\0';

    match_merge_state merge;
    merge_state_init(&merge);
    int rc = 0;
    size_t min_len = target_len;
    size_t max_len = target_len;
    if (metric == COUNT_METRIC_LEVENSHTEIN && indel_window != 0 && k == 1) {
        min_len = target_len > indel_window ? target_len - indel_window : 0;
        max_len = target_len + indel_window;
    }

    for (size_t len = min_len; len <= max_len; ++len) {
        if (len >= observed_cap) continue;
        if (target_start > seq_len || len > seq_len - target_start) continue;
        char candidate[8192];
        if (len >= sizeof(candidate)) continue;
        memcpy(candidate, seq + target_start, len);
        candidate[len] = '\0';
        uppercase_ascii(candidate);
        qdaln_match_result r;
        if (scan_assign_metric(candidate, len, targets, target_lens, n_targets, k, metric, &r) != 0) {
            rc = -1;
            break;
        }
        if (merge_state_add_result(&merge, observed, observed_cap, candidate, r) != 0) {
            rc = -1;
            break;
        }
    }
    merge_state_finish(&merge, result);
    merge_state_free(&merge);
    return rc;
}

static int scan_count_offsets(const char *const *targets, const size_t *target_lens, size_t n_targets,
                              const char *seq, size_t seq_len, const offset_list *offsets,
                              size_t fallback_offset, size_t target_len, int k, count_metric metric,
                              size_t indel_window, qdaln_match_result *result,
                              char *observed, size_t observed_cap) {
    *result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_INVALID};
    if (observed_cap != 0) observed[0] = '\0';
    size_t n_offsets = offsets == NULL || offsets->count == 0 ? 1 : offsets->count;
    match_merge_state merge;
    merge_state_init(&merge);
    int rc = 0;
    for (size_t i = 0; i < n_offsets; ++i) {
        size_t offset = offsets == NULL || offsets->count == 0 ? fallback_offset : offsets->items[i];
        qdaln_match_result candidate = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
        char local_observed[8192];
        if (scan_count_window(targets, target_lens, n_targets, seq, seq_len, offset, target_len, k, metric,
                              indel_window, &candidate, local_observed, sizeof(local_observed)) != 0) {
            rc = -1;
            break;
        }
        if (merge_state_add_result(&merge, observed, observed_cap, local_observed, candidate) != 0) {
            rc = -1;
            break;
        }
    }
    merge_state_finish(&merge, result);
    merge_state_free(&merge);
    return rc;
}

static void hamming_lookup_result_from_entry(const hamming_lookup_entry *entry, int distance,
                                             qdaln_match_result *result) {
    result->target_index = entry->target_index;
    result->best_distance = distance;
    result->second_best_distance = -1;
    result->match_count = entry->match_count;
    result->status = entry->match_count > 1 ? QDALN_MATCH_AMBIGUOUS : QDALN_MATCH_UNIQUE;
}

static void levenshtein1_lookup_result_from_entry(const hamming_lookup_entry *entry, int distance,
                                                  qdaln_match_result *result) {
    result->target_index = entry->target_index;
    result->best_distance = distance;
    result->second_best_distance = -1;
    result->match_count = entry->match_count;
    result->status = entry->match_count > 1 ? QDALN_MATCH_AMBIGUOUS : QDALN_MATCH_UNIQUE;
}

static int assign_levenshtein1_lookup_offset(const levenshtein1_lookup *lookup, const char *seq, size_t seq_len,
                                             size_t offset, qdaln_match_result *result, qdaln_index_stats *stats,
                                             char *observed, size_t observed_cap) {
    if (lookup == NULL || !lookup->ready || lookup->target_len == 0 || lookup->target_len > 31) return 0;
    *result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_INVALID};
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }
    if (observed_cap != 0) observed[0] = '\0';

    size_t target_len = lookup->target_len;
    if (offset > seq_len || target_len - 1 > seq_len - offset) return 1;

    uint64_t code = 0;
    const hamming_lookup_entry *entry = NULL;
    qdaln_match_result local = {-1, -1, -1, 0, QDALN_MATCH_NONE};
    int have_len_l = 0;
    uint64_t code_len_l = 0;
    if (target_len < observed_cap && target_len <= seq_len - offset &&
        dna2_code_local_fold(seq + offset, target_len, &code_len_l)) {
        have_len_l = 1;
        code = code_len_l;
        entry = levenshtein1_lookup_find(lookup->exact, lookup->exact_cap, code, target_len);
        if (entry != NULL) {
            levenshtein1_lookup_result_from_entry(entry, 0, &local);
            if (stats != NULL) {
                stats->candidates_considered += (size_t)entry->match_count;
                stats->candidates_verified += (size_t)entry->match_count;
            }
            copy_upper_ascii_window(observed, observed_cap, seq + offset, target_len);
            *result = local;
            return 1;
        }
    }

    match_merge_state merge;
    merge_state_init(&merge);
    int rc = 0;
    char candidate_observed[128];

    if (have_len_l) {
        entry = levenshtein1_lookup_find(lookup->substitution, lookup->substitution_cap, code_len_l, target_len);
        if (entry != NULL) {
            levenshtein1_lookup_result_from_entry(entry, 1, &local);
            if (stats != NULL) {
                stats->candidates_considered += (size_t)entry->match_count;
                stats->candidates_verified += (size_t)entry->match_count;
            }
            copy_upper_ascii_window(candidate_observed, sizeof(candidate_observed), seq + offset, target_len);
            rc = merge_state_add_result(&merge, observed, observed_cap, candidate_observed, local);
        }
    }

    if (rc == 0 && target_len - 1 < observed_cap && dna2_code_local_fold(seq + offset, target_len - 1, &code)) {
        entry = levenshtein1_lookup_find(lookup->target_deletion, lookup->target_deletion_cap, code, target_len - 1);
        if (entry != NULL) {
            levenshtein1_lookup_result_from_entry(entry, 1, &local);
            if (stats != NULL) {
                stats->candidates_considered += (size_t)entry->match_count;
                stats->candidates_verified += (size_t)entry->match_count;
            }
            copy_upper_ascii_window(candidate_observed, sizeof(candidate_observed), seq + offset, target_len - 1);
            rc = merge_state_add_result(&merge, observed, observed_cap, candidate_observed, local);
        }
    }

    if (rc == 0 && target_len + 1 < observed_cap && target_len + 1 <= seq_len - offset &&
        dna2_code_local_fold(seq + offset, target_len + 1, &code)) {
        entry = levenshtein1_lookup_find(lookup->target_insertion, lookup->target_insertion_cap, code, target_len + 1);
        if (entry != NULL) {
            levenshtein1_lookup_result_from_entry(entry, 1, &local);
            if (stats != NULL) {
                stats->candidates_considered += (size_t)entry->match_count;
                stats->candidates_verified += (size_t)entry->match_count;
            }
            copy_upper_ascii_window(candidate_observed, sizeof(candidate_observed), seq + offset, target_len + 1);
            rc = merge_state_add_result(&merge, observed, observed_cap, candidate_observed, local);
        }
    }

    if (rc != 0) {
        merge_state_free(&merge);
        return -1;
    }
    merge_state_finish(&merge, result);
    merge_state_free(&merge);
    if (result->match_count > 0) {
        return 1;
    }

    if (target_len <= seq_len - offset && target_len < observed_cap) {
        copy_upper_ascii_window(observed, observed_cap, seq + offset, target_len);
    } else if (target_len - 1 < observed_cap) {
        copy_upper_ascii_window(observed, observed_cap, seq + offset, target_len - 1);
    }
    *result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_NONE};
    return 1;
}

static int hamming_lookup_counts_eligible(int count_only, int max_correction_qual, count_metric metric,
                                          size_t indel_window, int k, size_t target_len,
                                          hamming_index_strategy hamming_strategy) {
    if (!count_only || max_correction_qual >= 0) return 0;
    if (metric != COUNT_METRIC_HAMMING || indel_window != 0) return 0;
    if (k != 0 && k != 1) return 0;
    if (target_len > 32) return 0;
    return hamming_strategy == HAMMING_INDEX_PRECOMPUTE || hamming_strategy == HAMMING_INDEX_AUTO;
}

static int levenshtein1_lookup_counts_eligible(int count_only, int max_correction_qual, count_metric metric,
                                               size_t indel_window, int k, size_t target_len,
                                               size_t max_selected_offsets, FILE *assignments,
                                               FILE *ambiguous_out, FILE *unmatched_out,
                                               ambiguity_policy assignment_policy) {
    if (!count_only || max_correction_qual >= 0) return 0;
    if (assignments != NULL || ambiguous_out != NULL || unmatched_out != NULL) return 0;
    if (assignment_policy != AMBIGUITY_POLICY_BEST) return 0;
    if (metric != COUNT_METRIC_LEVENSHTEIN || indel_window != 1 || k != 1) return 0;
    if (max_selected_offsets > 1 || target_len == 0 || target_len > 31) return 0;
    return 1;
}

static int hamming_direct_worker_eligible(int lookup_eligible, ambiguity_policy assignment_policy, int k) {
    if (!lookup_eligible) return 0;
    if (k == 0) return 1;
    return assignment_policy == AMBIGUITY_POLICY_BEST || assignment_policy == AMBIGUITY_POLICY_RADIUS;
}

static const char *count_backend_mode_name(count_backend_mode mode) {
    switch (mode) {
        case COUNT_BACKEND_CPU:
            return "cpu";
        case COUNT_BACKEND_METAL:
            return "gpu-metal-experimental";
        case COUNT_BACKEND_AUTO:
        default:
            return "auto";
    }
}

static int parse_count_backend_mode(const char *value, count_backend_mode *mode_out) {
    if (strcmp(value, "auto") == 0) {
        *mode_out = COUNT_BACKEND_AUTO;
        return 0;
    }
    if (strcmp(value, "cpu") == 0) {
        *mode_out = COUNT_BACKEND_CPU;
        return 0;
    }
    if (strcmp(value, "gpu-metal-experimental") == 0 || strcmp(value, "metal") == 0) {
        *mode_out = COUNT_BACKEND_METAL;
        return 0;
    }
    return -1;
}

static int metal_hamming_count_eligible(count_backend_mode backend, int hamming_lookup_eligible,
                                        ambiguity_policy assignment_policy, int k, size_t max_selected_offsets,
                                        int fused_offset_detection) {
    if (backend != COUNT_BACKEND_METAL) return 0;
    if (!hamming_lookup_eligible || !qdmetal_available()) return 0;
    if (fused_offset_detection) return 0;
    if (max_selected_offsets > 1) return 0;
    if (k == 1 && assignment_policy != AMBIGUITY_POLICY_BEST) return 0;
    return 1;
}

static int build_packed_target_codes(const seq_table *targets, size_t target_len, uint64_t **codes_out) {
    if (targets == NULL || codes_out == NULL || target_len == 0 || target_len > 32) return 0;
    uint64_t *codes = (uint64_t *)calloc(targets->count == 0 ? 1 : targets->count, sizeof(uint64_t));
    if (codes == NULL) return -1;
    for (size_t i = 0; i < targets->count; ++i) {
        if (targets->records[i].len != target_len || !dna2_code_local(targets->records[i].seq, target_len, &codes[i])) {
            free(codes);
            return 0;
        }
    }
    *codes_out = codes;
    return 1;
}

static const char *metal_count_engine_name(size_t n_targets) {
    return n_targets >= 1024 ? "hamming_metal_seed_index" : "hamming_metal_brute_force";
}

static int direct_hamming_merge_lookup_entry(match_merge_state *merge, const hamming_lookup_entry *entry,
                                             int distance) {
    if (entry == NULL) return 0;
    qdaln_match_result r;
    hamming_lookup_result_from_entry(entry, distance, &r);
    return merge_state_add_result(merge, NULL, 0, "", r);
}

static int direct_hamming_collect_seed_hits(count_sample_job *job, match_merge_state *merge, unsigned char seed_id,
                                            uint64_t seed_code, uint64_t read_code) {
    const hamming_lookup *lookup = job->hlookup;
    size_t seed_len = seed_id == 0 ? lookup->seed0_len : lookup->target_len - lookup->seed0_len;
    size_t slot = seed_hash_local(seed_code, seed_len, seed_id, lookup->seed_hash_cap);
    for (int e = lookup->seed_heads[slot]; e >= 0; e = lookup->seeds[e].next) {
        const hamming_seed_entry *entry = &lookup->seeds[e];
        if (entry->seed_id != seed_id || entry->code != seed_code || entry->target_index < 0) continue;
        job->stats->candidates_considered += 1;
        job->stats->candidates_verified += 1;
        if (hamming_code_distance_local(read_code, lookup->target_codes[entry->target_index], lookup->target_len) > 1) {
            continue;
        }
        qdaln_match_result r = {entry->target_index, 1, -1, 1, QDALN_MATCH_UNIQUE};
        if (merge_state_add_result(merge, NULL, 0, "", r) != 0) return -1;
    }
    return 0;
}

static int assign_hamming_lookup_offsets(const hamming_lookup *lookup, const char *seq, size_t seq_len,
                                         const offset_list *offsets, size_t fallback_offset, int k,
                                         qdaln_match_result *result, qdaln_index_stats *stats,
                                         char *observed, size_t observed_cap, int exact_merge) {
    if (lookup == NULL || !lookup->ready || lookup->target_len >= observed_cap || (k != 0 && k != 1)) return 0;
    *result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_INVALID};
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }
    if (observed_cap != 0) observed[0] = '\0';

    size_t n_offsets = offsets == NULL || offsets->count == 0 ? 1 : offsets->count;
    match_merge_state merge;
    merge_state_init(&merge);
    qdaln_match_result fast_result = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
    int fast_saw_window = 0;
    int rc = 1;
    for (size_t i = 0; i < n_offsets; ++i) {
        size_t offset = offsets == NULL || offsets->count == 0 ? fallback_offset : offsets->items[i];
        if (offset > seq_len || lookup->target_len > seq_len - offset) continue;
        if (!exact_merge) fast_saw_window = 1;
        char candidate[8192];
        uint64_t code = 0;
        const char *window = seq + offset;
        if (!dna2_code_local_fold(window, lookup->target_len, &code)) {
            rc = 0;
            break;
        }

        const hamming_lookup_entry *entry = hamming_lookup_find(lookup->exact, lookup->exact_cap, code);
        qdaln_match_result r = {-1, -1, -1, 0, QDALN_MATCH_NONE};
        if (entry != NULL) {
            hamming_lookup_result_from_entry(entry, 0, &r);
            if (stats != NULL) {
                stats->candidates_considered += (size_t)entry->match_count;
                stats->candidates_verified += (size_t)entry->match_count;
            }
            copy_upper_ascii_window(candidate, sizeof(candidate), window, lookup->target_len);
            if (exact_merge ? merge_state_add_result(&merge, observed, observed_cap, candidate, r) != 0
                            : (merge_summary_result(&fast_result, observed, observed_cap, candidate, r), 0)) {
                rc = -1;
                break;
            }
            if (!(exact_merge && k == 1)) continue;
        }
        if (k == 1) {
            if (!exact_merge && fast_result.best_distance == 0) continue;
            entry = hamming_lookup_find(lookup->mismatch, lookup->mismatch_cap, code);
            if (entry != NULL) {
                hamming_lookup_result_from_entry(entry, 1, &r);
                if (stats != NULL) {
                    stats->candidates_considered += (size_t)entry->match_count;
                    stats->candidates_verified += (size_t)entry->match_count;
                }
                copy_upper_ascii_window(candidate, sizeof(candidate), window, lookup->target_len);
                if (exact_merge ? merge_state_add_result(&merge, observed, observed_cap, candidate, r) != 0
                                : (merge_summary_result(&fast_result, observed, observed_cap, candidate, r), 0)) {
                    rc = -1;
                    break;
                }
            } else {
                if (!exact_merge) continue;
                copy_upper_ascii_window(candidate, sizeof(candidate), window, lookup->target_len);
                if (merge_state_add_result(&merge, observed, observed_cap, candidate, r) != 0) {
                    rc = -1;
                    break;
                }
            }
        } else {
            if (!exact_merge) continue;
            copy_upper_ascii_window(candidate, sizeof(candidate), window, lookup->target_len);
            if (merge_state_add_result(&merge, observed, observed_cap, candidate, r) != 0) {
                rc = -1;
                break;
            }
        }
    }

    if (exact_merge) merge_state_finish(&merge, result);
    else {
        if (fast_result.status == QDALN_MATCH_INVALID && fast_saw_window) {
            fast_result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_NONE};
        }
        *result = fast_result;
    }
    merge_state_free(&merge);
    return rc;
}

typedef struct seq_buffer {
    char **items;
    size_t *lens;
    size_t count;
    size_t cap;
    char *fixed_items;
    size_t fixed_len;
    size_t fixed_cap;
    int fixed_active;
} seq_buffer;

static int seq_buffer_ptr_in_fixed(const seq_buffer *buffer, const char *ptr) {
    if (buffer->fixed_items == NULL || buffer->fixed_cap == 0) return 0;
    uintptr_t p = (uintptr_t)ptr;
    uintptr_t start = (uintptr_t)buffer->fixed_items;
    uintptr_t end = start + buffer->fixed_cap * (buffer->fixed_len + 1);
    return p >= start && p < end;
}

static void free_seq_buffer(seq_buffer *buffer) {
    if (buffer == NULL) return;
    for (size_t i = 0; i < buffer->count; ++i) {
        if (!seq_buffer_ptr_in_fixed(buffer, buffer->items[i])) free(buffer->items[i]);
    }
    free(buffer->fixed_items);
    free(buffer->items);
    free(buffer->lens);
    buffer->items = NULL;
    buffer->lens = NULL;
    buffer->fixed_items = NULL;
    buffer->count = 0;
    buffer->cap = 0;
    buffer->fixed_len = 0;
    buffer->fixed_cap = 0;
    buffer->fixed_active = 0;
}

/* reset_seq_buffer reclaims per-sequence mallocs (for variable-length batches) but
 * retains the items/lens arrays and fixed block (when active for uniform lengths)
 * so subsequent batches reuse the same allocations without per-batch malloc/free.
 * This is used in the high-throughput batched reader paths for large FASTQ.
 */
static void reset_seq_buffer(seq_buffer *buffer) {
    if (buffer == NULL) return;
    for (size_t i = 0; i < buffer->count; ++i) {
        if (!seq_buffer_ptr_in_fixed(buffer, buffer->items[i])) free(buffer->items[i]);
    }
    if (!buffer->fixed_active) {
        free(buffer->fixed_items);
        buffer->fixed_items = NULL;
        buffer->fixed_len = 0;
        buffer->fixed_cap = 0;
    }
    /* keep items/lens (and fixed block + active if still set) and cap for reuse */
    buffer->count = 0;
}

static int grow_seq_buffer(seq_buffer *buffer) {
    size_t old_cap = buffer->cap;
    size_t next_cap = buffer->cap == 0 ? 1024 : buffer->cap * 2;
    char **next_items = (char **)realloc(buffer->items, next_cap * sizeof(char *));
    if (next_items == NULL) return -1;
    buffer->items = next_items;
    size_t *next_lens = (size_t *)realloc(buffer->lens, next_cap * sizeof(size_t));
    if (next_lens == NULL) return -1;
    buffer->lens = next_lens;
    buffer->cap = next_cap;
    if (buffer->fixed_active) {
        char *next_fixed = (char *)realloc(buffer->fixed_items, next_cap * (buffer->fixed_len + 1));
        if (next_fixed == NULL) {
            buffer->cap = old_cap;
            return -1;
        }
        buffer->fixed_items = next_fixed;
        buffer->fixed_cap = next_cap;
        for (size_t i = 0; i < buffer->count; ++i) {
            buffer->items[i] = buffer->fixed_items + i * (buffer->fixed_len + 1);
        }
    }
    return 0;
}

static int reserve_seq_buffer(seq_buffer *buffer, size_t requested_cap) {
    if (requested_cap <= buffer->cap) return 0;
    char **next_items = (char **)realloc(buffer->items, requested_cap * sizeof(char *));
    if (next_items == NULL) return -1;
    buffer->items = next_items;
    size_t *next_lens = (size_t *)realloc(buffer->lens, requested_cap * sizeof(size_t));
    if (next_lens == NULL) return -1;
    buffer->lens = next_lens;
    buffer->cap = requested_cap;
    if (buffer->fixed_active) {
        char *next_fixed = (char *)realloc(buffer->fixed_items, requested_cap * (buffer->fixed_len + 1));
        if (next_fixed == NULL) return -1;
        buffer->fixed_items = next_fixed;
        buffer->fixed_cap = requested_cap;
        for (size_t i = 0; i < buffer->count; ++i) {
            buffer->items[i] = buffer->fixed_items + i * (buffer->fixed_len + 1);
        }
    }
    return 0;
}

static int push_seq_buffer(seq_buffer *buffer, const char *seq, size_t len) {
    if (buffer->count == buffer->cap && grow_seq_buffer(buffer) != 0) return -1;
    if (buffer->count == 0 && buffer->fixed_items == NULL && len <= 8191) {
        buffer->fixed_items = (char *)malloc(buffer->cap * (len + 1));
        if (buffer->fixed_items != NULL) {
            buffer->fixed_len = len;
            buffer->fixed_cap = buffer->cap;
            buffer->fixed_active = 1;
        }
    }
    if (buffer->fixed_active && len == buffer->fixed_len) {
        char *dst = buffer->fixed_items + buffer->count * (buffer->fixed_len + 1);
        memcpy(dst, seq, len);
        dst[len] = '\0';
        buffer->items[buffer->count] = dst;
    } else {
        if (buffer->fixed_active && len != buffer->fixed_len) buffer->fixed_active = 0;
        buffer->items[buffer->count] = xstrndup(seq, len);
        if (buffer->items[buffer->count] == NULL) return -1;
    }
    buffer->lens[buffer->count] = len;
    ++buffer->count;
    return 0;
}

static void direct_hamming_record_hit(int target_index, int match_count, int *best_target, int *ambiguous) {
    if (match_count > 1) *ambiguous = 1;
    if (*best_target < 0) {
        *best_target = target_index;
    } else if (target_index != *best_target) {
        if (target_index >= 0 && target_index < *best_target) *best_target = target_index;
        *ambiguous = 1;
    }
}

static void merge_count_stats(count_stats *dst, const count_stats *src) {
    dst->total += src->total;
    dst->unique += src->unique;
    dst->exact += src->exact;
    dst->corrected += src->corrected;
    dst->ambiguous += src->ambiguous;
    dst->unmatched += src->unmatched;
    dst->invalid += src->invalid;
    dst->candidates_considered += src->candidates_considered;
    dst->candidates_verified += src->candidates_verified;
}

static void free_count_dirty_slots(count_dirty_slots *dirty) {
    if (dirty == NULL) return;
    free(dirty->items);
    free(dirty->table);
    dirty->items = NULL;
    dirty->table = NULL;
    dirty->count = 0;
    dirty->cap = 0;
    dirty->table_cap = 0;
}

static size_t count_dirty_slot_hash(size_t slot) {
    uint64_t x = (uint64_t)slot;
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return (size_t)x;
}

static int count_dirty_slots_rehash(count_dirty_slots *dirty, size_t min_cap) {
    size_t cap = 32;
    while (cap < min_cap) cap *= 2;
    size_t *table = (size_t *)calloc(cap, sizeof(size_t));
    if (table == NULL) return -1;

    for (size_t i = 0; i < dirty->count; ++i) {
        size_t mask = cap - 1;
        size_t pos = count_dirty_slot_hash(dirty->items[i].slot) & mask;
        while (table[pos] != 0) pos = (pos + 1) & mask;
        table[pos] = i + 1;
    }

    free(dirty->table);
    dirty->table = table;
    dirty->table_cap = cap;
    return 0;
}

static int mark_count_dirty_slot(count_dirty_slots *dirty, size_t slot) {
    if (dirty == NULL) return 0;

    if (dirty->table_cap == 0 || ((dirty->count + 1) * 2) > dirty->table_cap) {
        if (count_dirty_slots_rehash(dirty, (dirty->count + 1) * 4) != 0) return -1;
    }

    size_t mask = dirty->table_cap - 1;
    size_t pos = count_dirty_slot_hash(slot) & mask;
    while (dirty->table[pos] != 0) {
        count_dirty_slot *item = &dirty->items[dirty->table[pos] - 1];
        if (item->slot == slot) {
            ++item->count;
            return 0;
        }
        pos = (pos + 1) & mask;
    }

    if (dirty->count == dirty->cap) {
        size_t next_cap = dirty->cap == 0 ? 16 : dirty->cap * 2;
        count_dirty_slot *next = (count_dirty_slot *)realloc(dirty->items, next_cap * sizeof(count_dirty_slot));
        if (next == NULL) return -1;
        dirty->items = next;
        dirty->cap = next_cap;
    }
    dirty->items[dirty->count].slot = slot;
    dirty->items[dirty->count].count = 1;
    dirty->table[pos] = dirty->count + 1;
    ++dirty->count;
    return 0;
}

static int increment_count_slot(count_sample_job *job, size_t slot) {
    if (job->dirty_slots != NULL) return mark_count_dirty_slot(job->dirty_slots, slot);
    ++job->counts[slot];
    return 0;
}

static int direct_hamming_apply_match_result(count_sample_job *job, qdaln_match_result result, int saw_window) {
    apply_ambiguity_policy(&result, job->assignment_policy);
    if (result.status == QDALN_MATCH_UNIQUE && result.target_index >= 0) {
        int kind = result.best_distance == 0 ? 0 : 1;
        if (increment_count_slot(job, ((job->sample_index * job->targets->count + (size_t)result.target_index) * 5) + (size_t)kind) != 0) {
            return -1;
        }
        ++job->stats->unique;
        if (result.best_distance == 0) ++job->stats->exact;
        else ++job->stats->corrected;
    } else if (result.status == QDALN_MATCH_AMBIGUOUS) {
        ++job->stats->ambiguous;
    } else if (result.status == QDALN_MATCH_NONE) {
        ++job->stats->unmatched;
    } else if (saw_window) {
        ++job->stats->unmatched;
    } else {
        ++job->stats->invalid;
    }
    return 0;
}

static void direct_hamming_visit_seed(const count_sample_job *job, unsigned char seed_id, uint64_t seed_code,
                                      uint64_t read_code, int *best_target, int *ambiguous) {
    const hamming_lookup *lookup = job->hlookup;
    size_t seed_len = seed_id == 0 ? lookup->seed0_len : lookup->target_len - lookup->seed0_len;
    size_t slot = seed_hash_local(seed_code, seed_len, seed_id, lookup->seed_hash_cap);
    for (int e = lookup->seed_heads[slot]; e >= 0; ) {
        int nexte = lookup->seeds[e].next;
        const hamming_seed_entry *entry = &lookup->seeds[e];
        if (entry->seed_id != seed_id || entry->code != seed_code) {
            e = nexte;
            continue;
        }
        int target_index = entry->target_index;
        if (target_index < 0) {
            e = nexte;
            continue;
        }
        job->stats->candidates_considered += 1;
        job->stats->candidates_verified += 1;
        if (hamming_code_distance_local(read_code, lookup->target_codes[target_index], lookup->target_len) > 1) {
            e = nexte;
            continue;
        }
        direct_hamming_record_hit(target_index, 1, best_target, ambiguous);
        e = nexte;
    }
}

static size_t selected_offset_at(const offset_list *offsets, size_t fallback_offset, size_t i) {
    return offsets == NULL || offsets->count == 0 ? fallback_offset : offsets->items[i];
}

static int selected_offsets_are_sorted(const offset_list *offsets) {
    if (offsets == NULL || offsets->count < 2) return 1;
    for (size_t i = 1; i < offsets->count; ++i) {
        if (offsets->items[i] < offsets->items[i - 1]) return 0;
    }
    return 1;
}

static void fill_direct_hamming_codes(const char *seq, size_t seq_len, const offset_list *offsets,
                                      size_t fallback_offset, size_t target_len, uint64_t *codes,
                                      unsigned char *valid, unsigned char *invalid_counts,
                                      unsigned char *bad_positions, size_t n_offsets,
                                      int *saw_window, int *saw_non_acgt_window) {
    *saw_window = 0;
    *saw_non_acgt_window = 0;
    memset(valid, 0, n_offsets);
    memset(invalid_counts, 0, n_offsets);
    memset(bad_positions, 0, n_offsets);
    if (target_len == 0 || target_len > 32) return;

    if (!selected_offsets_are_sorted(offsets)) {
        for (size_t i = 0; i < n_offsets; ++i) {
            size_t offset = selected_offset_at(offsets, fallback_offset, i);
            if (offset > seq_len || target_len > seq_len - offset) continue;
            *saw_window = 1;
            uint64_t code = 0;
            unsigned char n_bad = 0;
            unsigned char bad_pos = 0;
            for (size_t j = 0; j < target_len; ++j) {
                uint64_t value = 0;
                if (!dna2_base_fold_value(seq[offset + j], &value)) {
                    if (n_bad < 255) ++n_bad;
                    bad_pos = (unsigned char)j;
                }
                code |= value << (2 * j);
            }
            if (n_bad != 0) {
                *saw_non_acgt_window = 1;
                codes[i] = code;
                invalid_counts[i] = n_bad;
                bad_positions[i] = bad_pos;
                continue;
            }
            valid[i] = 1;
            codes[i] = code;
        }
        return;
    }

    uint64_t code = 0;
    size_t invalid_count = 0;
    size_t current_offset = 0;
    int have_window = 0;
    for (size_t i = 0; i < n_offsets; ++i) {
        size_t offset = selected_offset_at(offsets, fallback_offset, i);
        if (offset > seq_len || target_len > seq_len - offset) break;

        if (!have_window) {
            code = 0;
            invalid_count = 0;
            size_t last_bad = 0;
            for (size_t j = 0; j < target_len; ++j) {
                uint64_t value = 0;
                if (!dna2_base_fold_value(seq[offset + j], &value)) {
                    ++invalid_count;
                    last_bad = j;
                }
                code |= value << (2 * j);
            }
            current_offset = offset;
            have_window = 1;
            if (invalid_count == 1) bad_positions[i] = (unsigned char)last_bad;
        } else {
            while (current_offset < offset) {
                uint64_t outgoing = 0;
                if (!dna2_base_fold_value(seq[current_offset], &outgoing) && invalid_count != 0) --invalid_count;
                (void)outgoing;
                code >>= 2;
                uint64_t incoming = 0;
                if (!dna2_base_fold_value(seq[current_offset + target_len], &incoming)) ++invalid_count;
                code |= incoming << (2 * (target_len - 1));
                ++current_offset;
            }
            if (invalid_count == 1) {
                for (size_t j = 0; j < target_len; ++j) {
                    uint64_t value = 0;
                    if (!dna2_base_fold_value(seq[current_offset + j], &value)) {
                        bad_positions[i] = (unsigned char)j;
                        break;
                    }
                }
            }
        }

        *saw_window = 1;
        if (invalid_count == 0) {
            valid[i] = 1;
            codes[i] = code;
        } else {
            *saw_non_acgt_window = 1;
            codes[i] = code;
            invalid_counts[i] = invalid_count > 255 ? 255 : (unsigned char)invalid_count;
        }
    }
}

static int direct_hamming_count_seq(count_sample_job *job, const char *seq, size_t seq_len) {
    if (job->hlookup == NULL || !job->hlookup->ready) return 0;
    ++job->stats->total;
    count_progress_tick(job->progress);

    size_t n_offsets = job->selected_offsets == NULL || job->selected_offsets->count == 0
            ? 1 : job->selected_offsets->count;
    /* The single-offset shortcut returns an exact hit before collecting its
     * one-mismatch neighbours. That is valid for best-distance semantics, but
     * radius semantics must retain every compatible target so an exact hit
     * plus a distance-one hit remains ambiguous. Route that case through the
     * merge path below. */
    if (n_offsets == 1 && !(job->k == 1 && job->assignment_policy == AMBIGUITY_POLICY_RADIUS)) {
        size_t offset = job->selected_offsets == NULL || job->selected_offsets->count == 0
                ? job->target_start : job->selected_offsets->items[0];
        if (offset > seq_len || job->hlookup->target_len > seq_len - offset) {
            ++job->stats->invalid;
            return 0;
        }

        uint64_t code = 0;
        unsigned char invalid_count = 0;
        unsigned char bad_position = 0;
        for (size_t j = 0; j < job->hlookup->target_len; ++j) {
            uint64_t value = 0;
            if (!dna2_base_fold_value(seq[offset + j], &value)) {
                if (invalid_count < 255) ++invalid_count;
                bad_position = (unsigned char)j;
            }
            code |= value << (2 * j);
        }

        int exact_target = -1;
        int exact_ambiguous = 0;
        if (invalid_count == 0) {
            const hamming_lookup_entry *entry =
                    hamming_lookup_find(job->hlookup->exact, job->hlookup->exact_cap, code);
            if (entry != NULL) {
                job->stats->candidates_considered += (unsigned long long)entry->match_count;
                job->stats->candidates_verified += (unsigned long long)entry->match_count;
                if (job->k == 0) {
                    if (entry->match_count > 1) {
                        ++job->stats->ambiguous;
                    } else {
                        if (increment_count_slot(job, ((job->sample_index * job->targets->count +
                                                        (size_t)entry->target_index) * 5) + 0) != 0) {
                            return -1;
                        }
                        ++job->stats->unique;
                        ++job->stats->exact;
                    }
                    return 0;
                }
                direct_hamming_record_hit(entry->target_index, entry->match_count, &exact_target,
                                          &exact_ambiguous);
            }
        }

        if (job->k == 0) {
            ++job->stats->unmatched;
            return 0;
        }

        if (exact_target >= 0) {
            if (exact_ambiguous) {
                ++job->stats->ambiguous;
            } else {
                if (increment_count_slot(job, ((job->sample_index * job->targets->count + (size_t)exact_target) * 5) + 0) != 0) {
                    return -1;
                }
                ++job->stats->unique;
                ++job->stats->exact;
            }
            return 0;
        }

        int mismatch_target = -1;
        int mismatch_ambiguous = 0;
        if (job->k == 1) {
            if (invalid_count == 0) {
                if (job->hlookup->seed_ready) {
                    uint64_t seed0 = code_segment_local(code, 0, job->hlookup->seed0_len);
                    uint64_t seed1 = code_segment_local(code, job->hlookup->seed0_len,
                                                        job->hlookup->target_len - job->hlookup->seed0_len);
                    direct_hamming_visit_seed(job, 0, seed0, code, &mismatch_target, &mismatch_ambiguous);
                    direct_hamming_visit_seed(job, 1, seed1, code, &mismatch_target, &mismatch_ambiguous);
                } else {
                    const hamming_lookup_entry *entry =
                            hamming_lookup_find(job->hlookup->mismatch, job->hlookup->mismatch_cap, code);
                    if (entry != NULL) {
                        job->stats->candidates_considered += (unsigned long long)entry->match_count;
                        job->stats->candidates_verified += (unsigned long long)entry->match_count;
                        direct_hamming_record_hit(entry->target_index, entry->match_count, &mismatch_target,
                                                  &mismatch_ambiguous);
                    }
                }
            } else if (invalid_count == 1) {
                uint64_t shift = (uint64_t)2 * bad_position;
                uint64_t mask = 3ULL << shift;
                for (uint64_t b = 0; b < 4; ++b) {
                    uint64_t patched = (code & ~mask) | (b << shift);
                    const hamming_lookup_entry *entry =
                            hamming_lookup_find(job->hlookup->exact, job->hlookup->exact_cap, patched);
                    if (entry == NULL) continue;
                    job->stats->candidates_considered += (unsigned long long)entry->match_count;
                    job->stats->candidates_verified += (unsigned long long)entry->match_count;
                    direct_hamming_record_hit(entry->target_index, entry->match_count, &mismatch_target,
                                              &mismatch_ambiguous);
                }
            }
        }

        if (mismatch_target >= 0) {
            if (mismatch_ambiguous) {
                ++job->stats->ambiguous;
            } else {
                if (increment_count_slot(job, ((job->sample_index * job->targets->count + (size_t)mismatch_target) * 5) + 1) != 0) {
                    return -1;
                }
                ++job->stats->unique;
                ++job->stats->corrected;
            }
        } else {
            ++job->stats->unmatched;
        }
        return 0;
    }

    uint64_t inline_codes[64];
    unsigned char inline_valid[64];
    unsigned char inline_invalid_counts[64];
    unsigned char inline_bad_positions[64];
    uint64_t *codes = inline_codes;
    unsigned char *valid = inline_valid;
    unsigned char *invalid_counts = inline_invalid_counts;
    unsigned char *bad_positions = inline_bad_positions;
    if (n_offsets > sizeof(inline_codes) / sizeof(inline_codes[0])) {
        codes = (uint64_t *)malloc(n_offsets * sizeof(uint64_t));
        valid = (unsigned char *)malloc(n_offsets);
        invalid_counts = (unsigned char *)malloc(n_offsets);
        bad_positions = (unsigned char *)malloc(n_offsets);
        if (codes == NULL || valid == NULL || invalid_counts == NULL || bad_positions == NULL) {
            free(codes);
            free(valid);
            free(invalid_counts);
            free(bad_positions);
            ++job->stats->invalid;
            return -1;
        }
    }

    int saw_window = 0;
    int saw_non_acgt_window = 0;
    fill_direct_hamming_codes(seq, seq_len, job->selected_offsets, job->target_start,
                              job->hlookup->target_len, codes, valid, invalid_counts,
                              bad_positions, n_offsets,
                              &saw_window, &saw_non_acgt_window);
    (void)saw_non_acgt_window;

    if (job->k == 1 && job->assignment_policy == AMBIGUITY_POLICY_RADIUS) {
        match_merge_state merge;
        merge_state_init(&merge);
        int rc = 0;
        for (size_t i = 0; i < n_offsets && rc == 0; ++i) {
            if (valid[i]) {
                const hamming_lookup_entry *entry =
                        hamming_lookup_find(job->hlookup->exact, job->hlookup->exact_cap, codes[i]);
                if (entry != NULL) {
                    job->stats->candidates_considered += (unsigned long long)entry->match_count;
                    job->stats->candidates_verified += (unsigned long long)entry->match_count;
                    rc = direct_hamming_merge_lookup_entry(&merge, entry, 0);
                }
                if (rc == 0) {
                    if (job->hlookup->seed_ready) {
                        uint64_t seed0 = code_segment_local(codes[i], 0, job->hlookup->seed0_len);
                        uint64_t seed1 = code_segment_local(codes[i], job->hlookup->seed0_len,
                                                            job->hlookup->target_len - job->hlookup->seed0_len);
                        if (direct_hamming_collect_seed_hits(job, &merge, 0, seed0, codes[i]) != 0 ||
                            direct_hamming_collect_seed_hits(job, &merge, 1, seed1, codes[i]) != 0) {
                            rc = -1;
                        }
                    } else {
                        entry = hamming_lookup_find(job->hlookup->mismatch, job->hlookup->mismatch_cap, codes[i]);
                        if (entry != NULL) {
                            job->stats->candidates_considered += (unsigned long long)entry->match_count;
                            job->stats->candidates_verified += (unsigned long long)entry->match_count;
                            rc = direct_hamming_merge_lookup_entry(&merge, entry, 1);
                        }
                    }
                }
            } else if (invalid_counts[i] == 1) {
                uint64_t shift = (uint64_t)2 * bad_positions[i];
                uint64_t mask = 3ULL << shift;
                for (uint64_t b = 0; b < 4 && rc == 0; ++b) {
                    uint64_t patched = (codes[i] & ~mask) | (b << shift);
                    const hamming_lookup_entry *entry =
                            hamming_lookup_find(job->hlookup->exact, job->hlookup->exact_cap, patched);
                    if (entry == NULL) continue;
                    job->stats->candidates_considered += (unsigned long long)entry->match_count;
                    job->stats->candidates_verified += (unsigned long long)entry->match_count;
                    rc = direct_hamming_merge_lookup_entry(&merge, entry, 1);
                }
            }
        }
        qdaln_match_result result = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
        if (rc == 0) {
            merge_state_finish(&merge, &result);
            if (result.status == QDALN_MATCH_INVALID && saw_window) {
                result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_NONE};
            }
            rc = direct_hamming_apply_match_result(job, result, saw_window);
        }
        merge_state_free(&merge);
        if (codes != inline_codes) {
            free(codes);
            free(valid);
            free(invalid_counts);
            free(bad_positions);
        }
        return rc;
    }

    int exact_target = -1;
    int exact_ambiguous = 0;
    for (size_t i = 0; i < n_offsets; ++i) {
        if (!valid[i]) continue;
        const hamming_lookup_entry *entry =
                hamming_lookup_find(job->hlookup->exact, job->hlookup->exact_cap, codes[i]);
        if (entry == NULL) continue;
        job->stats->candidates_considered += (unsigned long long)entry->match_count;
        job->stats->candidates_verified += (unsigned long long)entry->match_count;
        direct_hamming_record_hit(entry->target_index, entry->match_count, &exact_target, &exact_ambiguous);
    }

    if (exact_target >= 0) {
        if (exact_ambiguous) {
            ++job->stats->ambiguous;
        } else {
            if (increment_count_slot(job, ((job->sample_index * job->targets->count + (size_t)exact_target) * 5) + 0) != 0) {
                if (codes != inline_codes) {
                    free(codes);
                    free(valid);
                    free(invalid_counts);
                    free(bad_positions);
                }
                return -1;
            }
            ++job->stats->unique;
            ++job->stats->exact;
        }
        if (codes != inline_codes) {
            free(codes);
            free(valid);
            free(invalid_counts);
            free(bad_positions);
        }
        return 0;
    }

    int mismatch_target = -1;
    int mismatch_ambiguous = 0;
    if (job->k == 1) {
        for (size_t i = 0; i < n_offsets; ++i) {
            if (valid[i]) {
                if (job->hlookup->seed_ready) {
                    uint64_t seed0 = code_segment_local(codes[i], 0, job->hlookup->seed0_len);
                    uint64_t seed1 = code_segment_local(codes[i], job->hlookup->seed0_len,
                                                        job->hlookup->target_len - job->hlookup->seed0_len);
                    direct_hamming_visit_seed(job, 0, seed0, codes[i], &mismatch_target, &mismatch_ambiguous);
                    direct_hamming_visit_seed(job, 1, seed1, codes[i], &mismatch_target, &mismatch_ambiguous);
                } else {
                    const hamming_lookup_entry *entry =
                            hamming_lookup_find(job->hlookup->mismatch, job->hlookup->mismatch_cap, codes[i]);
                    if (entry == NULL) continue;
                    job->stats->candidates_considered += (unsigned long long)entry->match_count;
                    job->stats->candidates_verified += (unsigned long long)entry->match_count;
                    direct_hamming_record_hit(entry->target_index, entry->match_count, &mismatch_target,
                                              &mismatch_ambiguous);
                }
            } else if (invalid_counts[i] == 1) {
                uint64_t shift = (uint64_t)2 * bad_positions[i];
                uint64_t mask = 3ULL << shift;
                for (uint64_t b = 0; b < 4; ++b) {
                    uint64_t patched = (codes[i] & ~mask) | (b << shift);
                    const hamming_lookup_entry *entry =
                            hamming_lookup_find(job->hlookup->exact, job->hlookup->exact_cap, patched);
                    if (entry == NULL) continue;
                    job->stats->candidates_considered += (unsigned long long)entry->match_count;
                    job->stats->candidates_verified += (unsigned long long)entry->match_count;
                    direct_hamming_record_hit(entry->target_index, entry->match_count, &mismatch_target,
                                              &mismatch_ambiguous);
                }
            }
        }
    }

    if (mismatch_target >= 0) {
        if (mismatch_ambiguous) {
            ++job->stats->ambiguous;
        } else {
            if (increment_count_slot(job, ((job->sample_index * job->targets->count + (size_t)mismatch_target) * 5) + 1) != 0) {
                if (codes != inline_codes) {
                    free(codes);
                    free(valid);
                    free(invalid_counts);
                    free(bad_positions);
                }
                return -1;
            }
            ++job->stats->unique;
            ++job->stats->corrected;
        }
    } else if (saw_window) {
        ++job->stats->unmatched;
    } else {
        ++job->stats->invalid;
    }

    if (codes != inline_codes) {
        free(codes);
        free(valid);
        free(invalid_counts);
        free(bad_positions);
    }
    return 0;
}

typedef struct direct_hamming_batch_job {
    count_sample_job job;
    char **items;
    size_t *lens;
    size_t start;
    size_t end;
    count_dirty_slots dirty_slots;
    count_stats local_stats;
    int rc;
} direct_hamming_batch_job;

static void *direct_hamming_batch_worker(void *arg) {
    direct_hamming_batch_job *batch = (direct_hamming_batch_job *)arg;
    batch->job.sample_index = 0;
    batch->job.stats = &batch->local_stats;
    batch->job.dirty_slots = &batch->dirty_slots;
    batch->rc = 0;
    for (size_t i = batch->start; i < batch->end; ++i) {
        if (direct_hamming_count_seq(&batch->job, batch->items[i], batch->lens[i]) != 0) {
            batch->rc = 1;
            break;
        }
    }
    return NULL;
}

static int process_exact_hamming_window_buffer(count_sample_job *job, const seq_buffer *buffer) {
    if (buffer->count == 0) return 0;
    if (job->index == NULL) return 1;

    size_t offset = job->selected_offsets == NULL || job->selected_offsets->count == 0
            ? job->target_start
            : job->selected_offsets->items[0];
    size_t alloc_count = alloc_count_or_one(buffer->count);
    const char **windows = (const char **)malloc(alloc_count * sizeof(char *));
    size_t *window_lens = (size_t *)malloc(alloc_count * sizeof(size_t));
    qdaln_match_result *results = (qdaln_match_result *)malloc(alloc_count * sizeof(qdaln_match_result));
    if (windows == NULL || window_lens == NULL || results == NULL) {
        free(windows);
        free(window_lens);
        free(results);
        return 1;
    }

    size_t n_windows = 0;
    for (size_t i = 0; i < buffer->count; ++i) {
        ++job->stats->total;
        count_progress_tick(job->progress);
        if (offset > buffer->lens[i] || job->target_len > buffer->lens[i] - offset) {
            ++job->stats->invalid;
            continue;
        }
        windows[n_windows] = buffer->items[i] + offset;
        window_lens[n_windows] = job->target_len;
        ++n_windows;
    }

    int rc = 0;
    if (n_windows != 0) {
        qdaln_index_stats batch_stats = {0, 0};
        if (qdaln_index_lookup_exact_ascii_many_stats(job->index, windows, window_lens,
                                                      n_windows, results, &batch_stats) != 0) {
            rc = 1;
        } else {
            job->stats->candidates_considered += (unsigned long long)batch_stats.candidates_considered;
            job->stats->candidates_verified += (unsigned long long)batch_stats.candidates_verified;
            for (size_t i = 0; i < n_windows; ++i) {
                if (direct_hamming_apply_match_result(job, results[i], 1) != 0) {
                    rc = 1;
                    break;
                }
            }
        }
    }

    free(windows);
    free(window_lens);
    free(results);
    return rc;
}

static int process_direct_hamming_buffer(count_sample_job *job, const seq_buffer *buffer) {
    if (buffer->count == 0) return 0;
    if (job->k == 0 && job->index != NULL && job->hlookup != NULL && job->hlookup->ready &&
        job->hlookup->target_len == job->target_len && (job->selected_offsets == NULL || job->selected_offsets->count <= 1) &&
        job->read_threads <= 1) {
        return process_exact_hamming_window_buffer(job, buffer);
    }
    size_t read_threads = job->read_threads;
    if (read_threads <= 1 || buffer->count < 1024) {
        for (size_t i = 0; i < buffer->count; ++i) {
            if (direct_hamming_count_seq(job, buffer->items[i], buffer->lens[i]) != 0) return 1;
        }
        return 0;
    }
    if (read_threads > buffer->count) read_threads = buffer->count;

    size_t target_slots = job->targets->count * 5;
    pthread_t *thread_ids = (pthread_t *)calloc(read_threads, sizeof(pthread_t));
    direct_hamming_batch_job *jobs = (direct_hamming_batch_job *)calloc(read_threads, sizeof(direct_hamming_batch_job));
    if (thread_ids == NULL || jobs == NULL) {
        free(thread_ids);
        free(jobs);
        return 1;
    }

    size_t launched = 0;
    int rc = 0;
    for (size_t t = 0; t < read_threads; ++t) {
        size_t start = (buffer->count * t) / read_threads;
        size_t end = (buffer->count * (t + 1)) / read_threads;
        jobs[t].job = *job;
        jobs[t].items = buffer->items;
        jobs[t].lens = buffer->lens;
        jobs[t].start = start;
        jobs[t].end = end;
        if (pthread_create(&thread_ids[t], NULL, direct_hamming_batch_worker, &jobs[t]) != 0) {
            rc = 1;
            break;
        }
        ++launched;
    }

    for (size_t t = 0; t < launched; ++t) {
        pthread_join(thread_ids[t], NULL);
        if (jobs[t].rc != 0) rc = 1;
    }
    if (rc == 0) {
        size_t dst_offset = job->sample_index * target_slots;
        for (size_t t = 0; t < launched; ++t) {
            merge_count_stats(job->stats, &jobs[t].local_stats);
            for (size_t i = 0; i < jobs[t].dirty_slots.count; ++i) {
                size_t slot = jobs[t].dirty_slots.items[i].slot;
                job->counts[dst_offset + slot] += jobs[t].dirty_slots.items[i].count;
            }
        }
    }

    for (size_t t = 0; t < read_threads; ++t) free_count_dirty_slots(&jobs[t].dirty_slots);
    free(thread_ids);
    free(jobs);
    return rc;
}

static void score_offsets_for_seq(const hamming_lookup *lookup, const char *seq, size_t seq_len,
                                  size_t target_start, size_t target_len, size_t range,
                                  unsigned long long *scores) {
    if (lookup == NULL || !lookup->ready || lookup->target_len != target_len) return;
    size_t n_offsets = 0;
    if (offset_count_for_range(range, &n_offsets) != 0) return;
    for (size_t oi = 0; oi < n_offsets; ++oi) {
        long delta = (long)oi - (long)range;
        if (delta < 0 && target_start < (size_t)(-delta)) continue;
        size_t offset = delta < 0 ? target_start - (size_t)(-delta) : target_start + (size_t)delta;
        if (offset > seq_len || target_len > seq_len - offset) continue;
        uint64_t code = 0;
        if (!dna2_code_local_fold(seq + offset, target_len, &code)) continue;
        const hamming_lookup_entry *entry = hamming_lookup_find(lookup->exact, lookup->exact_cap, code);
        if (entry != NULL && entry->match_count == 1) ++scores[oi];
    }
}

static int select_offsets_from_scores(size_t target_start, size_t range, const unsigned long long *scores,
                                      size_t checked, offset_mode mode, double min_fraction,
                                      offset_list *selected_offsets);

static int count_sample_worker_direct_hamming(count_sample_job *job) {
    fastq_reader reader = {0};
    if (fastq_reader_open(&reader, job->reads_path) != 0) {
        fprintf(stderr, "failed to open FASTQ input\n");
        return 1;
    }

    char seq[8192];
    int got = 0;

    if (job->fused_offset_detection) {
        size_t n_offsets = 0;
        if (offset_count_for_range(job->auto_offset, &n_offsets) != 0) {
            fastq_reader_close(&reader);
            return 1;
        }
        unsigned long long *scores = (unsigned long long *)calloc(n_offsets, sizeof(unsigned long long));
        if (scores == NULL) {
            fastq_reader_close(&reader);
            return 1;
        }
        seq_buffer buffered = {0};
        if (reserve_seq_buffer(&buffered, job->auto_offset_sample) != 0) {
            free_seq_buffer(&buffered);
            free(scores);
            fastq_reader_close(&reader);
            return 1;
        }
        size_t checked = 0;
        size_t seq_len = 0;
        while (checked < job->auto_offset_sample &&
               (got = fastq_read_sequence_record_len(&reader, seq, sizeof(seq), &seq_len)) == 1) {
            score_offsets_for_seq(job->hlookup, seq, seq_len, job->target_start, job->target_len,
                                  job->auto_offset, scores);
            if (push_seq_buffer(&buffered, seq, seq_len) != 0) {
                free_seq_buffer(&buffered);
                free(scores);
                fastq_reader_close(&reader);
                return 1;
            }
            ++checked;
        }
        if (got < 0 ||
            select_offsets_from_scores(job->target_start, job->auto_offset, scores, checked, job->offsets_mode,
                                       job->offset_min_fraction, job->selected_offsets) != 0) {
            free_seq_buffer(&buffered);
            free(scores);
            fastq_reader_close(&reader);
            return 1;
        }
        free(scores);
        if (process_direct_hamming_buffer(job, &buffered) != 0) {
            free_seq_buffer(&buffered);
            fastq_reader_close(&reader);
            return 1;
        }
        free_seq_buffer(&buffered);
    }

    size_t seq_len = 0;
    if (job->read_threads <= 1) {
        while ((got = fastq_read_sequence_record_len(&reader, seq, sizeof(seq), &seq_len)) == 1) {
            if (direct_hamming_count_seq(job, seq, seq_len) != 0) {
                fastq_reader_close(&reader);
                return 1;
            }
        }
    } else {
        const size_t batch_reads = 1048576;
        seq_buffer batch = {0};
        if (reserve_seq_buffer(&batch, batch_reads) != 0) {
            free_seq_buffer(&batch);
            fastq_reader_close(&reader);
            return 1;
        }
        while ((got = fastq_read_sequence_record_len(&reader, seq, sizeof(seq), &seq_len)) == 1) {
            if (push_seq_buffer(&batch, seq, seq_len) != 0) {
                free_seq_buffer(&batch);
                fastq_reader_close(&reader);
                return 1;
            }
            if (batch.count == batch_reads) {
                if (process_direct_hamming_buffer(job, &batch) != 0) {
                    free_seq_buffer(&batch);
                    fastq_reader_close(&reader);
                    return 1;
                }
                reset_seq_buffer(&batch);
                if (reserve_seq_buffer(&batch, batch_reads) != 0) {
                    free_seq_buffer(&batch);
                    fastq_reader_close(&reader);
                    return 1;
                }
            }
        }
        if (got >= 0 && batch.count != 0 && process_direct_hamming_buffer(job, &batch) != 0) {
            free_seq_buffer(&batch);
            fastq_reader_close(&reader);
            return 1;
        }
        free_seq_buffer(&batch);
    }
    fastq_reader_close(&reader);
    if (got < 0) {
        fprintf(stderr, "malformed FASTQ input\n");
        return 1;
    }
    return 0;
}

static int pack_read_window_code(const char *seq, size_t seq_len, size_t offset, size_t target_len, uint64_t *code_out) {
    if (offset > seq_len || target_len > seq_len - offset) return 0;
    return dna2_code_local(seq + offset, target_len, code_out);
}

static int apply_metal_match_to_counts(count_sample_job *job, const qdmetal_match_result *metal_result) {
    qdaln_match_result result = {metal_result->target_index, metal_result->best_distance,
                                 metal_result->second_best_distance, metal_result->match_count,
                                 metal_result->status};
    return direct_hamming_apply_match_result(job, result, 1);
}

static int count_sample_worker_metal_hamming(count_sample_job *job) {
    if (job->metal_target_codes == NULL) return 1;
    fastq_reader reader = {0};
    if (fastq_reader_open(&reader, job->reads_path) != 0) {
        fprintf(stderr, "failed to open FASTQ input\n");
        return 1;
    }

    size_t window_offset = job->selected_offsets == NULL || job->selected_offsets->count == 0
            ? job->target_start
            : job->selected_offsets->items[0];
    const size_t batch_cap = 262144;
    uint64_t *read_codes = (uint64_t *)malloc(batch_cap * sizeof(uint64_t));
    qdmetal_match_result *metal_results =
            (qdmetal_match_result *)malloc(batch_cap * sizeof(qdmetal_match_result));
    if (read_codes == NULL || metal_results == NULL) {
        free(read_codes);
        free(metal_results);
        fastq_reader_close(&reader);
        return 1;
    }

    char seq[8192];
    size_t seq_len = 0;
    size_t batch_count = 0;
    int got = 0;
    int rc = 0;

    while (rc == 0 && (got = fastq_read_sequence_record_len(&reader, seq, sizeof(seq), &seq_len)) == 1) {
        ++job->stats->total;
        count_progress_tick(job->progress);
        uint64_t code = 0;
        if (!pack_read_window_code(seq, seq_len, window_offset, job->target_len, &code)) {
            if (window_offset > seq_len || job->target_len > seq_len - window_offset) ++job->stats->invalid;
            else ++job->stats->unmatched;
            continue;
        }
        read_codes[batch_count++] = code;
        if (batch_count == batch_cap) {
            qdmetal_assign_stats mstats = {0, 0, NULL, NULL};
            if (qdmetal_hamming_assign(read_codes, batch_count, job->metal_target_codes, job->targets->count,
                                       job->target_len, job->k, metal_results, &mstats) != 0) {
                rc = 1;
                break;
            }
            job->stats->candidates_considered += (unsigned long long)mstats.candidates_considered;
            job->stats->candidates_verified += (unsigned long long)mstats.candidates_verified;
            for (size_t i = 0; i < batch_count; ++i) {
                if (apply_metal_match_to_counts(job, &metal_results[i]) != 0) {
                    rc = 1;
                    break;
                }
            }
            batch_count = 0;
        }
    }

    if (rc == 0 && got >= 0 && batch_count != 0) {
        qdmetal_assign_stats mstats = {0, 0, NULL, NULL};
        if (qdmetal_hamming_assign(read_codes, batch_count, job->metal_target_codes, job->targets->count,
                                   job->target_len, job->k, metal_results, &mstats) != 0) {
            rc = 1;
        } else {
            job->stats->candidates_considered += (unsigned long long)mstats.candidates_considered;
            job->stats->candidates_verified += (unsigned long long)mstats.candidates_verified;
            for (size_t i = 0; i < batch_count; ++i) {
                if (apply_metal_match_to_counts(job, &metal_results[i]) != 0) {
                    rc = 1;
                    break;
                }
            }
        }
    }

    free(read_codes);
    free(metal_results);
    fastq_reader_close(&reader);
    if (got < 0) {
        fprintf(stderr, "malformed FASTQ input\n");
        return 1;
    }
    return rc;
}

static int count_sample_sequence(count_sample_job *job, const char *seq, size_t seq_len, const char *qual,
                                 const char *read_id) {
    char observed[8192];
    qdaln_match_result result = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
    qdaln_index_stats istats = {0, 0};
    observed[0] = '\0';
    ++job->stats->total;
    count_progress_tick(job->progress);
    int best_exact_shortcut = job->k == 1 &&
            job->assignment_policy == AMBIGUITY_POLICY_BEST && job->assignments == NULL &&
            job->ambiguous_out == NULL && job->unmatched_out == NULL;
    int handled = 0;
    if (job->metric == COUNT_METRIC_LEVENSHTEIN && job->k == 1 && job->indel_window == 1 &&
        job->levlookup != NULL && job->levlookup->ready &&
        (job->selected_offsets == NULL || job->selected_offsets->count <= 1)) {
        size_t offset = job->selected_offsets == NULL || job->selected_offsets->count == 0
                ? job->target_start : job->selected_offsets->items[0];
        int lookup_rc = assign_levenshtein1_lookup_offset(job->levlookup, seq, seq_len, offset, &result,
                                                          &istats, observed, sizeof(observed));
        if (lookup_rc < 0) return -1;
        handled = lookup_rc;
    }
    if (job->metric == COUNT_METRIC_HAMMING && job->indel_window == 0 && job->hlookup != NULL && job->hlookup->ready) {
        int exact_merge = job->assignment_policy == AMBIGUITY_POLICY_RADIUS ||
                          job->assignments != NULL || job->ambiguous_out != NULL ||
                          job->unmatched_out != NULL;
        int lookup_rc = assign_hamming_lookup_offsets(job->hlookup, seq, seq_len, job->selected_offsets, 0,
                                                      job->k, &result, &istats, observed, sizeof(observed),
                                                      exact_merge);
        if (lookup_rc < 0) return -1;
        handled = lookup_rc;
    }
    if (!handled &&
        assign_count_offsets(job->index, seq, seq_len, job->selected_offsets, 0, job->target_len, job->k,
                             job->metric, job->indel_window, &result, &istats, observed,
                             sizeof(observed), best_exact_shortcut) != 0) {
        return -1;
    }
    apply_ambiguity_policy(&result, job->assignment_policy);
    if (result.status != QDALN_MATCH_INVALID) {
        job->stats->candidates_considered += (unsigned long long)istats.candidates_considered;
        job->stats->candidates_verified += (unsigned long long)istats.candidates_verified;
    }

    const char *correction = "invalid";
    int quality_rejected = 0;
    if (result.status == QDALN_MATCH_UNIQUE && result.target_index >= 0 && result.best_distance > 0) {
        seq_record *target = &job->targets->records[result.target_index];
        if (!quality_allows_unique_correction(seq, seq_len, qual, job->selected_offsets, 0, job->target_len,
                                              job->metric, job->indel_window, job->k, observed, target, result,
                                              job->max_correction_qual)) {
            result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_NONE};
            quality_rejected = 1;
        }
    }
    if (result.status == QDALN_MATCH_UNIQUE && result.target_index >= 0) {
        seq_record *target = &job->targets->records[result.target_index];
        int kind = correction_kind(observed, strlen(observed), target->seq, target->len, result.best_distance);
        size_t slot = ((job->sample_index * job->targets->count + (size_t)result.target_index) * 5) + (size_t)kind;
        if (increment_count_slot(job, slot) != 0) return -1;
        ++job->stats->unique;
        if (result.best_distance == 0) ++job->stats->exact;
        else ++job->stats->corrected;
        correction = correction_name(kind);
    } else if (result.status == QDALN_MATCH_AMBIGUOUS) {
        ++job->stats->ambiguous;
        correction = "ambiguous";
    } else if (result.status == QDALN_MATCH_NONE) {
        ++job->stats->unmatched;
        correction = quality_rejected ? "quality_rejected" : "none";
    } else {
        ++job->stats->invalid;
    }

    const char *id = read_id == NULL ? "" : read_id;
    if (job->assignments != NULL &&
        (result.status != QDALN_MATCH_AMBIGUOUS || strcmp(job->ambiguous_policy, "report") == 0)) {
        write_assignment_like_row(job->assignments, job->targets, job->sample_label, id, observed, result,
                                  correction);
    }
    if (job->ambiguous_out != NULL && result.status == QDALN_MATCH_AMBIGUOUS) {
        write_assignment_like_row(job->ambiguous_out, job->targets, job->sample_label, id, observed, result,
                                  correction);
    }
    if (job->unmatched_out != NULL &&
        (result.status == QDALN_MATCH_NONE || result.status == QDALN_MATCH_INVALID)) {
        write_assignment_like_row(job->unmatched_out, job->targets, job->sample_label, id, observed, result,
                                  correction);
    }
    return 0;
}

typedef struct count_batch_job {
    count_sample_job job;
    char **items;
    size_t *lens;
    size_t start;
    size_t end;
    count_dirty_slots dirty_slots;
    count_stats local_stats;
    int rc;
} count_batch_job;

static void *count_batch_worker(void *arg) {
    count_batch_job *batch = (count_batch_job *)arg;
    batch->job.sample_index = 0;
    batch->job.stats = &batch->local_stats;
    batch->job.dirty_slots = &batch->dirty_slots;
    batch->job.assignments = NULL;
    batch->job.ambiguous_out = NULL;
    batch->job.unmatched_out = NULL;
    batch->rc = 0;
    for (size_t i = batch->start; i < batch->end; ++i) {
        if (count_sample_sequence(&batch->job, batch->items[i], batch->lens[i], NULL, NULL) != 0) {
            batch->rc = 1;
            break;
        }
    }
    return NULL;
}

static int process_count_buffer(count_sample_job *job, const seq_buffer *buffer) {
    if (buffer->count == 0) return 0;
    if (job->k == 0 && job->metric == COUNT_METRIC_HAMMING && job->indel_window == 0 && job->index != NULL &&
        job->hlookup != NULL && job->hlookup->ready && job->hlookup->target_len == job->target_len &&
        job->assignments == NULL && job->ambiguous_out == NULL && job->unmatched_out == NULL &&
        job->max_correction_qual < 0 && job->assignment_policy == AMBIGUITY_POLICY_BEST &&
        job->read_threads <= 1) {
        return process_exact_hamming_window_buffer(job, buffer);
    }
    size_t read_threads = job->read_threads;
    if (read_threads <= 1 || buffer->count < 1024) {
        for (size_t i = 0; i < buffer->count; ++i) {
            if (count_sample_sequence(job, buffer->items[i], buffer->lens[i], NULL, NULL) != 0) return 1;
        }
        return 0;
    }
    if (read_threads > buffer->count) read_threads = buffer->count;

    size_t target_slots = job->targets->count * 5;
    pthread_t *thread_ids = (pthread_t *)calloc(read_threads, sizeof(pthread_t));
    count_batch_job *jobs = (count_batch_job *)calloc(read_threads, sizeof(count_batch_job));
    if (thread_ids == NULL || jobs == NULL) {
        free(thread_ids);
        free(jobs);
        return 1;
    }

    size_t launched = 0;
    int rc = 0;
    for (size_t t = 0; t < read_threads; ++t) {
        size_t start = (buffer->count * t) / read_threads;
        size_t end = (buffer->count * (t + 1)) / read_threads;
        jobs[t].job = *job;
        jobs[t].items = buffer->items;
        jobs[t].lens = buffer->lens;
        jobs[t].start = start;
        jobs[t].end = end;
        if (pthread_create(&thread_ids[t], NULL, count_batch_worker, &jobs[t]) != 0) {
            rc = 1;
            break;
        }
        ++launched;
    }

    for (size_t t = 0; t < launched; ++t) {
        pthread_join(thread_ids[t], NULL);
        if (jobs[t].rc != 0) rc = 1;
    }
    if (rc == 0) {
        size_t dst_offset = job->sample_index * target_slots;
        for (size_t t = 0; t < launched; ++t) {
            merge_count_stats(job->stats, &jobs[t].local_stats);
            for (size_t i = 0; i < jobs[t].dirty_slots.count; ++i) {
                size_t slot = jobs[t].dirty_slots.items[i].slot;
                job->counts[dst_offset + slot] += jobs[t].dirty_slots.items[i].count;
            }
        }
    }

    for (size_t t = 0; t < read_threads; ++t) free_count_dirty_slots(&jobs[t].dirty_slots);
    free(thread_ids);
    free(jobs);
    return rc;
}

static void *count_sample_worker(void *arg);

typedef struct count_samples_args {
    const qdaln_index *index;
    const hamming_lookup *hlookup;
    const levenshtein1_lookup *levlookup;
    const seq_table *targets;
    const char **target_ptrs;
    const size_t *target_lens;
    const string_list *reads;
    const string_list *labels;
    offset_list *selected_offsets;
    size_t target_len;
    int k;
    count_metric metric;
    size_t indel_window;
    unsigned long long *counts;
    count_stats *stats_by_sample;
    FILE *assignments;
    FILE *ambiguous_out;
    FILE *unmatched_out;
    const char *ambiguous_policy;
    ambiguity_policy assignment_policy;
    int direct_hamming_counts;
    int metal_hamming_counts;
    const uint64_t *metal_target_codes;
    int fused_offset_detection;
    size_t target_start;
    size_t auto_offset;
    size_t auto_offset_sample;
    offset_mode offsets_mode;
    double offset_min_fraction;
    size_t effective_read_threads;
    int max_correction_qual;
    size_t sample_threads;
    count_progress *progress_by_sample;
} count_samples_args;

static int run_count_samples_phase(const count_samples_args *args) {
    if (args == NULL || args->reads == NULL || args->labels == NULL || args->counts == NULL ||
        args->stats_by_sample == NULL) {
        return -1;
    }
    size_t sample_threads = args->sample_threads;
    if (sample_threads <= 1 || args->reads->count <= 1) {
        for (size_t sample = 0; sample < args->reads->count; ++sample) {
            count_sample_job job = {
                args->index, args->hlookup, args->levlookup, args->targets, args->target_ptrs, args->target_lens,
                args->reads->items[sample], args->labels->items[sample], sample, &args->selected_offsets[sample],
                args->target_len, args->k, args->metric, args->indel_window, args->counts, &args->stats_by_sample[sample],
                args->assignments, args->ambiguous_out, args->unmatched_out, args->ambiguous_policy,
                args->assignment_policy, args->direct_hamming_counts, args->metal_hamming_counts,
                args->metal_target_codes, args->fused_offset_detection, args->target_start, args->auto_offset,
                args->auto_offset_sample, args->offsets_mode, args->offset_min_fraction, args->effective_read_threads,
                args->max_correction_qual, 1, NULL,
                args->progress_by_sample != NULL ? &args->progress_by_sample[sample] : NULL};
            count_sample_worker(&job);
            if (job.rc != 0) return job.rc;
        }
        return 0;
    }

    pthread_t *thread_ids = (pthread_t *)calloc(sample_threads, sizeof(pthread_t));
    count_sample_job *jobs = (count_sample_job *)calloc(args->reads->count, sizeof(count_sample_job));
    if (thread_ids == NULL || jobs == NULL) {
        free(thread_ids);
        free(jobs);
        fprintf(stderr, "out of memory\n");
        return -1;
    }
    int rc = 0;
    size_t next_sample = 0;
    while (next_sample < args->reads->count && rc == 0) {
        size_t batch = args->reads->count - next_sample;
        if (batch > sample_threads) batch = sample_threads;
        for (size_t i = 0; i < batch; ++i) {
            size_t sample = next_sample + i;
            jobs[sample] = (count_sample_job){
                args->index, args->hlookup, args->levlookup, args->targets, args->target_ptrs, args->target_lens,
                args->reads->items[sample], args->labels->items[sample], sample, &args->selected_offsets[sample],
                args->target_len, args->k, args->metric, args->indel_window, args->counts, &args->stats_by_sample[sample],
                NULL, NULL, NULL, args->ambiguous_policy, args->assignment_policy, args->direct_hamming_counts,
                args->metal_hamming_counts, args->metal_target_codes, args->fused_offset_detection, args->target_start,
                args->auto_offset, args->auto_offset_sample, args->offsets_mode, args->offset_min_fraction, 1,
                args->max_correction_qual, 1, NULL,
                args->progress_by_sample != NULL ? &args->progress_by_sample[sample] : NULL};
            if (pthread_create(&thread_ids[i], NULL, count_sample_worker, &jobs[sample]) != 0) {
                fprintf(stderr, "failed to create worker thread\n");
                batch = i;
                rc = -1;
                break;
            }
        }
        for (size_t i = 0; i < batch; ++i) {
            pthread_join(thread_ids[i], NULL);
            if (rc == 0 && jobs[next_sample + i].rc != 0) rc = jobs[next_sample + i].rc;
        }
        next_sample += batch;
    }
    free(thread_ids);
    free(jobs);
    return rc;
}

static unsigned long long count_matrix_cell_total(const unsigned long long *counts, size_t n_targets, size_t sample,
                                                  size_t target) {
    unsigned long long total = 0;
    for (size_t kind = 0; kind < 5; ++kind) total += counts[((sample * n_targets + target) * 5) + kind];
    return total;
}

static int validate_metal_counts_against_cpu(const unsigned long long *metal_counts, const unsigned long long *cpu_counts,
                                             size_t n_samples, size_t n_targets, size_t *diff_guides_out,
                                             long long *delta_reads_out, char *example_guide, size_t example_cap) {
    if (metal_counts == NULL || cpu_counts == NULL || diff_guides_out == NULL || delta_reads_out == NULL) return -1;
    *diff_guides_out = 0;
    *delta_reads_out = 0;
    if (example_guide != NULL && example_cap > 0) example_guide[0] = '\0';
    for (size_t t = 0; t < n_targets; ++t) {
        int guide_diff = 0;
        for (size_t s = 0; s < n_samples; ++s) {
            unsigned long long metal_total = count_matrix_cell_total(metal_counts, n_targets, s, t);
            unsigned long long cpu_total = count_matrix_cell_total(cpu_counts, n_targets, s, t);
            if (metal_total != cpu_total) {
                guide_diff = 1;
                *delta_reads_out += (long long)metal_total - (long long)cpu_total;
                if (example_guide != NULL && example_guide[0] == '\0' && example_cap > 0) {
                    snprintf(example_guide, example_cap, "sample_index=%zu target_index=%zu metal=%llu cpu=%llu", s, t,
                             metal_total, cpu_total);
                }
            }
        }
        if (guide_diff) ++*diff_guides_out;
    }
    return *diff_guides_out == 0 ? 0 : 1;
}

static int env_truthy(const char *value) {
    if (value == NULL || value[0] == '\0') return 0;
    if (strcmp(value, "0") == 0) return 0;
    if (strcmp(value, "false") == 0 || strcmp(value, "FALSE") == 0) return 0;
    if (strcmp(value, "no") == 0 || strcmp(value, "NO") == 0) return 0;
    return 1;
}

static void *count_sample_worker(void *arg) {
    count_sample_job *job = (count_sample_job *)arg;
    if (job->metal_hamming_counts) {
        job->rc = count_sample_worker_metal_hamming(job);
        return NULL;
    }
    if (job->direct_hamming_counts) {
        job->rc = count_sample_worker_direct_hamming(job);
        return NULL;
    }

    fastq_reader reader = {0};
    if (fastq_reader_open(&reader, job->reads_path) != 0) {
        fprintf(stderr, "failed to open FASTQ input\n");
        job->rc = 1;
        return NULL;
    }

    char header[8192];
    char seq[8192];
    char plus[8192];
    char qual[8192];
    char read_id[8192];
    int got = 0;
    int need_read_id = job->assignments != NULL || job->ambiguous_out != NULL || job->unmatched_out != NULL;
    int need_quality = job->max_correction_qual >= 0;
    int need_full_record = need_read_id || need_quality;
    size_t seq_len = 0;
    if (job->read_threads > 1 && !need_full_record) {
        const size_t batch_reads = 1048576;
        seq_buffer batch = {0};
        if (reserve_seq_buffer(&batch, batch_reads) != 0) {
            free_seq_buffer(&batch);
            fastq_reader_close(&reader);
            job->rc = 1;
            return NULL;
        }
        while ((got = fastq_read_sequence_record_len(&reader, seq, sizeof(seq), &seq_len)) == 1) {
            if (push_seq_buffer(&batch, seq, seq_len) != 0) {
                free_seq_buffer(&batch);
                fastq_reader_close(&reader);
                job->rc = 1;
                return NULL;
            }
            if (batch.count == batch_reads) {
                if (process_count_buffer(job, &batch) != 0) {
                    free_seq_buffer(&batch);
                    fastq_reader_close(&reader);
                    fprintf(stderr, "FASTQ assignment failed\n");
                    job->rc = 1;
                    return NULL;
                }
                reset_seq_buffer(&batch);
                if (reserve_seq_buffer(&batch, batch_reads) != 0) {
                    free_seq_buffer(&batch);
                    fastq_reader_close(&reader);
                    job->rc = 1;
                    return NULL;
                }
            }
        }
        if (got >= 0 && batch.count != 0 && process_count_buffer(job, &batch) != 0) {
            free_seq_buffer(&batch);
            fastq_reader_close(&reader);
            fprintf(stderr, "FASTQ assignment failed\n");
            job->rc = 1;
            return NULL;
        }
        free_seq_buffer(&batch);
    } else {
        while ((got = need_full_record
                ? fastq_read_record_len(&reader, header, seq, plus, qual, sizeof(header), &seq_len)
                : fastq_read_sequence_record_len(&reader, seq, sizeof(seq), &seq_len)) == 1) {
            read_id[0] = '\0';
            if (need_read_id) fastq_read_id(header, read_id, sizeof(read_id));
            if (count_sample_sequence(job, seq, seq_len, job->max_correction_qual >= 0 ? qual : NULL,
                                      read_id) != 0) {
                fastq_reader_close(&reader);
                fprintf(stderr, "FASTQ assignment failed\n");
                job->rc = 1;
                return NULL;
            }
        }
    }
    fastq_reader_close(&reader);
    if (got < 0) {
        fprintf(stderr, "malformed FASTQ input\n");
        job->rc = 1;
        return NULL;
    }
    job->rc = 0;
    return NULL;
}

static int select_offsets_from_scores(size_t target_start, size_t range, const unsigned long long *scores,
                                      size_t checked, offset_mode mode, double min_fraction,
                                      offset_list *selected_offsets) {
    free_offset_list(selected_offsets);
    if (range == 0) return push_offset_unique(selected_offsets, target_start);

    size_t n_offsets = 0;
    if (offset_count_for_range(range, &n_offsets) != 0) return -1;

    size_t best_i = range;
    for (size_t oi = 0; oi < n_offsets; ++oi) {
        size_t best_dist = best_i > range ? best_i - range : range - best_i;
        size_t this_dist = oi > range ? oi - range : range - oi;
        if (scores[oi] > scores[best_i] || (scores[oi] == scores[best_i] && this_dist < best_dist)) {
            best_i = oi;
        }
    }

    int rc = 0;
    if (mode == OFFSET_MODE_MULTI && checked != 0) {
        for (size_t oi = 0; oi < n_offsets; ++oi) {
            double fraction = (double)scores[oi] / (double)checked;
            if (fraction + 1e-12 < min_fraction) continue;
            long delta = (long)oi - (long)range;
            if (delta < 0 && target_start < (size_t)(-delta)) continue;
            size_t offset = delta < 0 ? target_start - (size_t)(-delta) : target_start + (size_t)delta;
            if (push_offset_unique(selected_offsets, offset) != 0) {
                rc = -1;
                break;
            }
        }
    }
    if (rc == 0 && selected_offsets->count == 0 && scores[best_i] != 0) {
        long best_delta = (long)best_i - (long)range;
        size_t offset = best_delta < 0 ? target_start - (size_t)(-best_delta) : target_start + (size_t)best_delta;
        rc = push_offset_unique(selected_offsets, offset);
    }
    if (rc == 0 && selected_offsets->count == 0) {
        rc = push_offset_unique(selected_offsets, target_start);
    }
    return rc;
}

static int detect_offsets(const qdaln_index *index, const hamming_lookup *exact_lookup, const char *reads_path, size_t target_start,
                          size_t target_len, size_t range, size_t sample_limit, offset_mode mode,
                          double min_fraction, offset_list *selected_offsets) {
    if (range == 0) {
        free_offset_list(selected_offsets);
        return push_offset_unique(selected_offsets, target_start);
    }

    size_t n_offsets = 0;
    if (offset_count_for_range(range, &n_offsets) != 0) return -1;
    unsigned long long *scores = (unsigned long long *)calloc(n_offsets, sizeof(unsigned long long));
    if (scores == NULL) return -1;

    fastq_reader reader = {0};
    if (fastq_reader_open(&reader, reads_path) != 0) {
        free(scores);
        return -1;
    }

    char header[8192];
    char seq[8192];
    char plus[8192];
    char qual[8192];
    size_t checked = 0;
    int got = 0;
    size_t seq_len = 0;
    while (checked < sample_limit &&
           (got = fastq_read_record_len(&reader, header, seq, plus, qual, sizeof(header), &seq_len)) == 1) {
        if (exact_lookup != NULL && exact_lookup->ready && exact_lookup->target_len == target_len) {
            score_offsets_for_seq(exact_lookup, seq, seq_len, target_start, target_len, range, scores);
        } else {
            for (size_t oi = 0; oi < n_offsets; ++oi) {
                long delta = (long)oi - (long)range;
                if (delta < 0 && target_start < (size_t)(-delta)) continue;
                size_t offset = delta < 0 ? target_start - (size_t)(-delta) : target_start + (size_t)delta;
                if (offset > seq_len || target_len > seq_len - offset) continue;
                qdaln_match_result r;
                qdaln_index_stats s;
                if (qdaln_index_lookup_exact_ascii_stats(index, seq + offset, target_len, &r, &s) != 0) {
                    fastq_reader_close(&reader);
                    free(scores);
                    return -1;
                }
                if (r.status == QDALN_MATCH_UNIQUE) ++scores[oi];
            }
        }
        ++checked;
    }

    fastq_reader_close(&reader);
    if (got < 0) {
        free(scores);
        return -1;
    }

    int rc = select_offsets_from_scores(target_start, range, scores, checked, mode, min_fraction, selected_offsets);
    free(scores);
    return rc;
}

typedef struct offset_detect_job {
    const qdaln_index *index;
    const hamming_lookup *exact_lookup;
    const char *reads_path;
    size_t target_start;
    size_t target_len;
    size_t range;
    size_t sample_limit;
    offset_mode mode;
    double min_fraction;
    offset_list *selected_offsets;
    int rc;
} offset_detect_job;

static void *detect_offsets_worker(void *arg) {
    offset_detect_job *job = (offset_detect_job *)arg;
    job->rc = detect_offsets(job->index, job->exact_lookup, job->reads_path, job->target_start, job->target_len,
                             job->range, job->sample_limit, job->mode, job->min_fraction, job->selected_offsets);
    return NULL;
}

static int detect_offsets_for_samples(const qdaln_index *index, const hamming_lookup *exact_lookup,
                                      const string_list *reads, size_t target_start, size_t target_len,
                                      size_t range, size_t sample_limit, offset_mode mode, double min_fraction,
                                      offset_list *selected_offsets, size_t threads) {
    if (reads->count == 0) return 0;
    if (threads > reads->count) threads = reads->count;
    if (threads <= 1 || reads->count <= 1) {
        for (size_t sample = 0; sample < reads->count; ++sample) {
            if (detect_offsets(index, exact_lookup, reads->items[sample], target_start, target_len, range,
                               sample_limit, mode, min_fraction, &selected_offsets[sample]) != 0) {
                return -1;
            }
        }
        return 0;
    }

    pthread_t *thread_ids = (pthread_t *)calloc(threads, sizeof(pthread_t));
    offset_detect_job *jobs = (offset_detect_job *)calloc(reads->count, sizeof(offset_detect_job));
    if (thread_ids == NULL || jobs == NULL) {
        free(thread_ids);
        free(jobs);
        return -1;
    }

    int rc = 0;
    size_t next_sample = 0;
    while (next_sample < reads->count && rc == 0) {
        size_t batch = reads->count - next_sample;
        if (batch > threads) batch = threads;
        for (size_t i = 0; i < batch; ++i) {
            size_t sample = next_sample + i;
            jobs[sample] = (offset_detect_job){index, exact_lookup, reads->items[sample], target_start, target_len,
                                               range, sample_limit, mode, min_fraction,
                                               &selected_offsets[sample], 0};
            if (pthread_create(&thread_ids[i], NULL, detect_offsets_worker, &jobs[sample]) != 0) {
                batch = i;
                rc = -1;
                break;
            }
        }
        for (size_t i = 0; i < batch; ++i) {
            pthread_join(thread_ids[i], NULL);
            if (jobs[next_sample + i].rc != 0) rc = -1;
        }
        next_sample += batch;
    }

    free(thread_ids);
    free(jobs);
    return rc;
}

static int run_count(const char *argv0, int argc, char **argv) {
    const char *targets_path = NULL;
    const char *samples_path = NULL;
    const char *out_path = NULL;
    const char *assignments_path = NULL;
    const char *summary_path = NULL;
    const char *report_path = NULL;
    const char *report_audit_dir = NULL;
    const char *report_unmatched_path = NULL;
    const char *ambiguous_path = NULL;
    const char *unmatched_path = NULL;
    const char *sample_qc_path = NULL;
    const char *target_counts_long_path = NULL;
    const int crispr_mode = strcmp(argv[1], "crispr-count") == 0;
    const char *format = crispr_mode ? "mageck" : "dotmatch";
    const char *ambiguous_policy = "discard";
    ambiguity_policy assignment_policy = AMBIGUITY_POLICY_RADIUS;
    count_metric metric = COUNT_METRIC_LEVENSHTEIN;
    hamming_index_strategy hamming_strategy = HAMMING_INDEX_AUTO;
    size_t target_start = 0;
    size_t target_len = 0;
    size_t indel_window = 0;
    size_t auto_offset = 0;
    size_t auto_offset_sample = 1000;
    offset_mode offsets_mode = OFFSET_MODE_BEST;
    double offset_min_fraction = 0.005;
    size_t threads = 0;
    int max_correction_qual = -1;
    int k = -1;
    count_backend_mode backend_mode = COUNT_BACKEND_AUTO;
    string_list reads = {0};
    string_list labels = {0};
    int show_progress = isatty(STDERR_FILENO);
    size_t progress_interval_reads = 250000;
    char default_sample_qc_path[PATH_MAX];
    count_progress *progress_by_sample = NULL;
    int metal_validate = env_truthy(getenv("DOTMATCH_METAL_VALIDATE"));
    const char *metal_validation_status = NULL;

    int i = 2;
    while (i < argc) {
        const char *arg = argv[i++];
        if ((strcmp(arg, "--targets") == 0 || strcmp(arg, "--library") == 0) && i < argc) {
            targets_path = argv[i++];
        } else if (strcmp(arg, "--samples") == 0 && i < argc) {
            samples_path = argv[i++];
        } else if (strcmp(arg, "--reads") == 0 && i < argc) {
            if (push_string(&reads, argv[i++]) != 0) {
                fprintf(stderr, "out of memory\n");
                goto fail_args;
            }
        } else if (strcmp(arg, "--sample-label") == 0 && i < argc) {
            if (split_string_list(&labels, argv[i++], ',') != 0) {
                fprintf(stderr, "out of memory\n");
                goto fail_args;
            }
        } else if ((strcmp(arg, "--target-start") == 0 || strcmp(arg, "--guide-start") == 0) && i < argc) {
            if (parse_size_value(argv[i++], &target_start) != 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if ((strcmp(arg, "--target-length") == 0 || strcmp(arg, "--guide-length") == 0) && i < argc) {
            if (parse_size_value(argv[i++], &target_len) != 0 || target_len == 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--k") == 0 && i < argc) {
            if (parse_int_value(argv[i++], &k) != 0 || k < 0 || k > 3) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--metric") == 0 && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "hamming") == 0) {
                metric = COUNT_METRIC_HAMMING;
            } else if (strcmp(value, "levenshtein") == 0) {
                metric = COUNT_METRIC_LEVENSHTEIN;
            } else {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--hamming-index") == 0 && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "auto") == 0) {
                hamming_strategy = HAMMING_INDEX_AUTO;
            } else if (strcmp(value, "query") == 0) {
                hamming_strategy = HAMMING_INDEX_QUERY;
            } else if (strcmp(value, "precompute") == 0) {
                hamming_strategy = HAMMING_INDEX_PRECOMPUTE;
            } else {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--indel-window") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &indel_window) != 0 || indel_window > 1) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--auto-offset") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &auto_offset) != 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--auto-offset-sample") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &auto_offset_sample) != 0 || auto_offset_sample == 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--offset-mode") == 0 && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "best") == 0) {
                offsets_mode = OFFSET_MODE_BEST;
            } else if (strcmp(value, "multi") == 0) {
                offsets_mode = OFFSET_MODE_MULTI;
            } else {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--offset-min-fraction") == 0 && i < argc) {
            if (parse_double_value(argv[i++], &offset_min_fraction) != 0 ||
                offset_min_fraction < 0.0 || offset_min_fraction > 1.0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--backend") == 0 && i < argc) {
            if (parse_count_backend_mode(argv[i++], &backend_mode) != 0) {
                fprintf(stderr, "--backend must be auto, cpu, or gpu-metal-experimental\n");
                goto fail_args;
            }
        } else if (strcmp(arg, "--progress") == 0) {
            show_progress = 1;
        } else if (strcmp(arg, "--no-progress") == 0) {
            show_progress = 0;
        } else if (strcmp(arg, "--progress-interval") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &progress_interval_reads) != 0 || progress_interval_reads == 0) {
                fprintf(stderr, "--progress-interval must be a positive integer\n");
                goto fail_args;
            }
        } else if (strcmp(arg, "--metal-validate") == 0) {
            metal_validate = 1;
        } else if (strcmp(arg, "--threads") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &threads) != 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--max-correction-qual") == 0 && i < argc) {
            if (parse_int_value(argv[i++], &max_correction_qual) != 0 ||
                max_correction_qual < 0 || max_correction_qual > 93) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--out") == 0 && i < argc) {
            out_path = argv[i++];
        } else if (strcmp(arg, "--assignments") == 0 && i < argc) {
            assignments_path = argv[i++];
        } else if (strcmp(arg, "--summary") == 0 && i < argc) {
            summary_path = argv[i++];
        } else if (strcmp(arg, "--report") == 0 && i < argc) {
            report_path = argv[i++];
        } else if (strcmp(arg, "--report-audit-dir") == 0 && i < argc) {
            report_audit_dir = argv[i++];
        } else if (strcmp(arg, "--report-unmatched") == 0 && i < argc) {
            report_unmatched_path = argv[i++];
        } else if (strcmp(arg, "--qc") == 0 && i < argc) {
            sample_qc_path = argv[i++];
        } else if (strcmp(arg, "--sample-qc") == 0 && i < argc) {
            sample_qc_path = argv[i++];
        } else if (strcmp(arg, "--target-counts-long") == 0 && i < argc) {
            target_counts_long_path = argv[i++];
        } else if (strcmp(arg, "--ambiguous-out") == 0 && i < argc) {
            ambiguous_path = argv[i++];
        } else if (strcmp(arg, "--unmatched-out") == 0 && i < argc) {
            unmatched_path = argv[i++];
        } else if (strcmp(arg, "--ambiguous") == 0 && i < argc) {
            ambiguous_policy = argv[i++];
            if (strcmp(ambiguous_policy, "discard") != 0 && strcmp(ambiguous_policy, "report") != 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--ambiguity-policy") == 0 && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "best") == 0) {
                assignment_policy = AMBIGUITY_POLICY_BEST;
            } else if (strcmp(value, "radius") == 0) {
                assignment_policy = AMBIGUITY_POLICY_RADIUS;
            } else {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(arg, "--format") == 0 && i < argc) {
            format = argv[i++];
            if (strcmp(format, "dotmatch") != 0 && strcmp(format, "mageck") != 0) {
                usage(argv0);
                goto fail_args;
            }
        } else {
            usage(argv0);
            goto fail_args;
        }
    }

    if (samples_path != NULL && read_samples_file(samples_path, &labels, &reads) != 0) {
        fprintf(stderr, "failed to read samples file\n");
        goto fail_args;
    }
    if (targets_path == NULL || reads.count == 0 || out_path == NULL || target_len == 0 || k < 0) {
        usage(argv0);
        goto fail_args;
    }
    if (crispr_mode && sample_qc_path == NULL &&
        derive_output_sibling_path(out_path, "sample_qc.tsv", default_sample_qc_path,
                                   sizeof(default_sample_qc_path)) == 0) {
        sample_qc_path = default_sample_qc_path;
    }
    if (auto_offset > MAX_AUTO_OFFSET) {
        fprintf(stderr, "--auto-offset must be <= %d\n", MAX_AUTO_OFFSET);
        goto fail_args;
    }
    if (metric == COUNT_METRIC_HAMMING && indel_window != 0) {
        fprintf(stderr, "--indel-window is only valid with --metric levenshtein\n");
        goto fail_args;
    }
    if (metric == COUNT_METRIC_LEVENSHTEIN && k > 2) {
        fprintf(stderr, "--metric levenshtein supports --k up to 2\n");
        goto fail_args;
    }
    if (indel_window != 0 && k != 1) {
        fprintf(stderr, "--indel-window requires --k 1\n");
        goto fail_args;
    }
    if (labels.count == 0) {
        for (size_t i = 0; i < reads.count; ++i) {
            if (push_string(&labels, path_basename(reads.items[i])) != 0) {
                fprintf(stderr, "out of memory\n");
                goto fail_args;
            }
        }
    }
    if (labels.count != reads.count) {
        fprintf(stderr, "--sample-label count must match --reads count\n");
        goto fail_args;
    }
    if (validate_unique_sample_labels(&labels, crispr_mode ? "--samples sample" : "--sample-label") != 0) {
        goto fail_args;
    }
    if (threads > 1 && (assignments_path != NULL || ambiguous_path != NULL || unmatched_path != NULL)) {
        fprintf(stderr, "--threads > 1 is not supported with row-level diagnostic outputs\n");
        goto fail_args;
    }
    if (threads == 0) {
        size_t auto_t = get_cpu_count();
        if (auto_t > 1 && (assignments_path != NULL || ambiguous_path != NULL || unmatched_path != NULL)) {
            threads = 1; /* fall back for safe ordered diagnostic outputs */
        } else {
            threads = auto_t;
        }
    }
    int count_only = assignments_path == NULL && ambiguous_path == NULL && unmatched_path == NULL;

    seq_table targets = {0};
    qdaln_index *index = NULL;
    hamming_lookup hlookup = {0};
    hamming_lookup offset_lookup = {0};
    levenshtein1_lookup levlookup = {0};
    const char **target_ptrs = NULL;
    size_t *target_lens = NULL;
    unsigned char *ambiguous_nearby = NULL;
    unsigned long long *counts = NULL;
    count_stats *stats_by_sample = NULL;
    offset_list *selected_offsets = NULL;
    FILE *out = NULL;
    FILE *assignments = NULL;
    FILE *ambiguous_out = NULL;
    FILE *unmatched_out = NULL;
    int rc = 1;
    double run_start_seconds = seconds_now();
    double target_index_seconds = 0.0;
    double offset_detection_seconds = 0.0;
    double hamming_precompute_seconds = 0.0;
    double counting_seconds = 0.0;
    const char *offset_detection_strategy = auto_offset == 0 ? "none" : "prepass";
    const char *count_engine = "generic_indexed";
    size_t effective_read_threads = threads;
    const char *backend_effective = "cpu";
    int metal_hamming_counts = 0;
    uint64_t *metal_target_codes = NULL;

    double phase_start_seconds = seconds_now();
    if (read_target_table(targets_path, &targets) != 0) {
        fprintf(stderr, "failed to read targets\n");
        goto done;
    }
    int target_id_check = validate_unique_seq_ids(&targets, crispr_mode ? "guide" : "target");
    if (target_id_check != 0) {
        if (target_id_check == -1) fprintf(stderr, "out of memory\n");
        goto done;
    }
    if (metric == COUNT_METRIC_HAMMING && !all_targets_have_length(&targets, target_len)) {
        fprintf(stderr, "--metric hamming requires every target to have --target-length bases\n");
        goto done;
    }
    target_index_seconds = seconds_now() - phase_start_seconds;

    int hamming_lookup_eligible = hamming_lookup_counts_eligible(
            count_only, max_correction_qual, metric, indel_window, k, target_len, hamming_strategy);
    int direct_hamming_counts = hamming_direct_worker_eligible(hamming_lookup_eligible, assignment_policy, k);
    int may_use_metal = backend_mode == COUNT_BACKEND_METAL && hamming_lookup_eligible && qdmetal_available() &&
            (k == 0 || assignment_policy == AMBIGUITY_POLICY_BEST);
    if (direct_hamming_counts) {
        phase_start_seconds = seconds_now();
        int use_mismatch_precompute_now = k == 1 &&
                (hamming_strategy == HAMMING_INDEX_PRECOMPUTE ||
                 (hamming_strategy == HAMMING_INDEX_AUTO && auto_offset != 0 &&
                  offsets_mode == OFFSET_MODE_MULTI));
        int lookup_rc = 0;
        if (k == 0) {
            lookup_rc = build_hamming_exact_lookup(&targets, target_len, &hlookup);
        } else if (use_mismatch_precompute_now) {
            lookup_rc = build_hamming_lookup(&targets, target_len, &hlookup);
        } else {
            lookup_rc = build_hamming_seed_lookup(&targets, target_len, &hlookup);
        }
        if (lookup_rc != 0) {
            fprintf(stderr, "failed to build Hamming lookup\n");
            goto done;
        }
        hamming_precompute_seconds = seconds_now() - phase_start_seconds;
        if (hlookup.ready) {
            count_engine = "hamming_lookup_direct";
        } else {
            direct_hamming_counts = 0;
        }
    }

    int need_general_index = !direct_hamming_counts || strcmp(format, "dotmatch") == 0;
    if (need_general_index) {
        phase_start_seconds = seconds_now();
        if (build_target_arrays(&targets, &target_ptrs, &target_lens) != 0) {
            fprintf(stderr, "out of memory\n");
            goto done;
        }
        index = qdaln_index_build(target_ptrs, target_lens, targets.count);
        if (index == NULL) {
            fprintf(stderr, "failed to build target index\n");
            goto done;
        }
        target_index_seconds += seconds_now() - phase_start_seconds;
    }

    size_t sample_target_slots = 0;
    size_t total_slots = 0;
    if (checked_mul_size(reads.count, targets.count, &sample_target_slots) != 0 ||
        checked_mul_size(sample_target_slots, 5, &total_slots) != 0) {
        fprintf(stderr, "count matrix is too large\n");
        goto done;
    }
    counts = (unsigned long long *)calloc(alloc_count_or_one(total_slots), sizeof(unsigned long long));
    stats_by_sample = (count_stats *)calloc(alloc_count_or_one(reads.count), sizeof(count_stats));
    ambiguous_nearby = (unsigned char *)calloc(alloc_count_or_one(targets.count), sizeof(unsigned char));
    selected_offsets = (offset_list *)calloc(alloc_count_or_one(reads.count), sizeof(offset_list));
    if (counts == NULL || stats_by_sample == NULL || ambiguous_nearby == NULL || selected_offsets == NULL) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    for (size_t sample = 0; sample < reads.count; ++sample) {
        if (push_offset_unique(&selected_offsets[sample], target_start) != 0) {
            fprintf(stderr, "out of memory\n");
            goto done;
        }
    }

    if (strcmp(format, "dotmatch") == 0) {
        for (size_t i = 0; i < targets.count; ++i) {
            qdaln_match_result r;
            qdaln_index_stats s;
            const char *seq_ptr = targets.records[i].seq;
            size_t seq_len = targets.records[i].len;
            int assign_rc = metric == COUNT_METRIC_HAMMING
                    ? qdaln_index_assign_hamming_stats(index, &seq_ptr, &seq_len, 1, k, &r, &s)
                    : qdaln_index_assign_stats(index, &seq_ptr, &seq_len, 1, k, &r, &s);
            if (assign_rc != 0) {
                fprintf(stderr, "target ambiguity check failed\n");
                goto done;
            }
            ambiguous_nearby[i] = r.match_count > 1 ? 1 : 0;
        }
    }

    if (assignments_path != NULL) {
        assignments = open_output_file(assignments_path);
        if (assignments == NULL) {
            fprintf(stderr, "failed to open assignments output\n");
            goto done;
        }
        fprintf(assignments, "sample\tread_id\tobserved_seq\ttarget_index\ttarget_id\ttarget_seq\tbest_distance\tsecond_best_distance\tmatch_count\tstatus\tcorrection\n");
    }
    if (ambiguous_path != NULL) {
        ambiguous_out = open_output_file(ambiguous_path);
        if (ambiguous_out == NULL) {
            fprintf(stderr, "failed to open ambiguous output\n");
            goto done;
        }
        fprintf(ambiguous_out, "sample\tread_id\tobserved_seq\ttarget_index\ttarget_id\ttarget_seq\tbest_distance\tsecond_best_distance\tmatch_count\tstatus\tcorrection\n");
    }
    if (unmatched_path != NULL) {
        unmatched_out = open_output_file(unmatched_path);
        if (unmatched_out == NULL) {
            fprintf(stderr, "failed to open unmatched output\n");
            goto done;
        }
        fprintf(unmatched_out, "sample\tread_id\tobserved_seq\ttarget_index\ttarget_id\ttarget_seq\tbest_distance\tsecond_best_distance\tmatch_count\tstatus\tcorrection\n");
    }

    int metal_blocks_fused_offset = may_use_metal && offsets_mode != OFFSET_MODE_MULTI;
    int fused_offset_detection = direct_hamming_counts && auto_offset != 0 && !metal_blocks_fused_offset;
    if (fused_offset_detection) {
        offset_detection_strategy = "fused";
    }

    if (auto_offset != 0 && !fused_offset_detection) {
        phase_start_seconds = seconds_now();
        const hamming_lookup *offset_lookup_ptr = hlookup.ready ? &hlookup : NULL;
        if (metric == COUNT_METRIC_HAMMING && offset_lookup_ptr == NULL) {
            int lookup_rc = build_hamming_exact_lookup(&targets, target_len, &offset_lookup);
            if (lookup_rc != 0) {
                fprintf(stderr, "failed to build offset detection lookup\n");
                goto done;
            }
            if (offset_lookup.ready) offset_lookup_ptr = &offset_lookup;
        }
        if (detect_offsets_for_samples(index, offset_lookup_ptr, &reads, target_start, target_len, auto_offset,
                                       auto_offset_sample, offsets_mode, offset_min_fraction, selected_offsets,
                                       threads) != 0) {
            fprintf(stderr, "automatic offset detection failed\n");
            goto done;
        }
        offset_detection_seconds = seconds_now() - phase_start_seconds;
        free_hamming_lookup(&offset_lookup);
    }

    size_t max_selected_offsets = 0;
    for (size_t sample = 0; sample < reads.count; ++sample) {
        if (selected_offsets[sample].count > max_selected_offsets) max_selected_offsets = selected_offsets[sample].count;
    }
    int direct_levenshtein_counts = levenshtein1_lookup_counts_eligible(
            count_only, max_correction_qual, metric, indel_window, k, target_len, max_selected_offsets,
            assignments, ambiguous_out, unmatched_out, assignment_policy);
    if (direct_levenshtein_counts) {
        phase_start_seconds = seconds_now();
        int lookup_rc = build_levenshtein1_lookup(&targets, target_len, &levlookup);
        if (lookup_rc != 0) {
            fprintf(stderr, "failed to build Levenshtein k=1 lookup\n");
            goto done;
        }
        if (levlookup.ready) {
            count_engine = "levenshtein_k1_lookup_direct";
            target_index_seconds += seconds_now() - phase_start_seconds;
        } else {
            direct_levenshtein_counts = 0;
        }
    }
    if (metal_hamming_count_eligible(backend_mode, hamming_lookup_eligible, assignment_policy, k,
                                     max_selected_offsets, fused_offset_detection)) {
        int pack_rc = build_packed_target_codes(&targets, target_len, &metal_target_codes);
        if (pack_rc == 1) {
            metal_hamming_counts = 1;
            direct_hamming_counts = 0;
            count_engine = metal_count_engine_name(targets.count);
            backend_effective = "gpu-metal-experimental";
        } else if (pack_rc < 0) {
            goto done;
        } else if (backend_mode == COUNT_BACKEND_METAL) {
            fprintf(stderr,
                    "Metal backend requires packable A/C/G/T targets, count-only output, single offset, and "
                    "best-distance policy for k=1\n");
            goto done;
        }
    } else if (backend_mode == COUNT_BACKEND_METAL) {
        fprintf(stderr,
                "Metal backend unavailable for this workload; requires Darwin Metal, --metric hamming, k 0|1, "
                "count-only output, and single offset");
        if (k == 1) fprintf(stderr, " with --ambiguity-policy best");
        fprintf(stderr, "\n");
        goto done;
    } else if (direct_hamming_counts && !fused_offset_detection && max_selected_offsets <= 1) {
        count_engine = "hamming_lookup_direct_single_offset";
    }
    int use_precomputed_hamming = metric == COUNT_METRIC_HAMMING && k == 1 &&
            (hamming_strategy == HAMMING_INDEX_PRECOMPUTE ||
             (hamming_strategy == HAMMING_INDEX_AUTO && max_selected_offsets > 1));
    if (use_precomputed_hamming && (!hlookup.ready || hlookup.mismatch == NULL)) {
        phase_start_seconds = seconds_now();
        free_hamming_lookup(&hlookup);
        int lookup_rc = build_hamming_lookup(&targets, target_len, &hlookup);
        if (lookup_rc != 0) {
            fprintf(stderr, "failed to build Hamming lookup\n");
            goto done;
        }
        hamming_precompute_seconds = seconds_now() - phase_start_seconds;
    }

    if (show_progress) {
        progress_by_sample = (count_progress *)calloc(reads.count, sizeof(count_progress));
        if (progress_by_sample == NULL) {
            fprintf(stderr, "out of memory\n");
            goto done;
        }
        for (size_t sample = 0; sample < reads.count; ++sample) {
            count_progress_init(&progress_by_sample[sample], labels.items[sample], progress_interval_reads);
        }
    }

    phase_start_seconds = seconds_now();
    size_t sample_threads = threads;
    if ((direct_hamming_counts || metal_hamming_counts) && reads.count == 1 && threads > 1) {
        effective_read_threads = threads;
        sample_threads = 1;
    } else if (count_only && reads.count == 1 && threads > 1) {
        effective_read_threads = threads;
        sample_threads = 1;
    } else if (sample_threads > reads.count) {
        sample_threads = reads.count;
    }
    count_samples_args sample_args = {
        index, &hlookup, &levlookup, &targets, target_ptrs, target_lens, &reads, &labels, selected_offsets, target_len, k, metric,
        indel_window, counts, stats_by_sample, assignments, ambiguous_out, unmatched_out, ambiguous_policy,
        assignment_policy, direct_hamming_counts, metal_hamming_counts, metal_target_codes, fused_offset_detection,
        target_start, auto_offset, auto_offset_sample, offsets_mode, offset_min_fraction, effective_read_threads,
        max_correction_qual, sample_threads, progress_by_sample};
    if (run_count_samples_phase(&sample_args) != 0) goto done;

    if (metal_validate && metal_hamming_counts) {
        size_t count_slots = total_slots == 0 ? 1 : total_slots;
        size_t sample_slots = reads.count == 0 ? 1 : reads.count;
        unsigned long long *metal_counts_snapshot =
                (unsigned long long *)calloc(count_slots, sizeof(unsigned long long));
        count_stats *metal_stats_snapshot = (count_stats *)calloc(sample_slots, sizeof(count_stats));
        if (metal_counts_snapshot == NULL || metal_stats_snapshot == NULL) {
            free(metal_counts_snapshot);
            free(metal_stats_snapshot);
            fprintf(stderr, "out of memory\n");
            goto done;
        }
        memcpy(metal_counts_snapshot, counts, total_slots * sizeof(unsigned long long));
        memcpy(metal_stats_snapshot, stats_by_sample, reads.count * sizeof(count_stats));
        memset(counts, 0, count_slots * sizeof(unsigned long long));
        memset(stats_by_sample, 0, sample_slots * sizeof(count_stats));

        int cpu_direct = hamming_direct_worker_eligible(hamming_lookup_eligible, assignment_policy, k) && hlookup.ready;
        if (!cpu_direct) {
            fprintf(stderr, "Metal validation requires a CPU Hamming direct lookup for this workload\n");
            free(metal_counts_snapshot);
            free(metal_stats_snapshot);
            goto done;
        }
        count_samples_args cpu_args = sample_args;
        cpu_args.direct_hamming_counts = 1;
        cpu_args.metal_hamming_counts = 0;
        cpu_args.metal_target_codes = NULL;
        cpu_args.assignments = NULL;
        cpu_args.ambiguous_out = NULL;
        cpu_args.unmatched_out = NULL;
        cpu_args.progress_by_sample = NULL;
        if (run_count_samples_phase(&cpu_args) != 0) {
            free(metal_counts_snapshot);
            free(metal_stats_snapshot);
            goto done;
        }

        size_t diff_guides = 0;
        long long delta_reads = 0;
        char example[256];
        example[0] = '\0';
        if (validate_metal_counts_against_cpu(metal_counts_snapshot, counts, reads.count, targets.count, &diff_guides,
                                              &delta_reads, example, sizeof(example)) != 0) {
            metal_validation_status = "failed";
            fprintf(stderr,
                    "dotmatch: Metal validation failed: %zu guides differ across samples (net read delta %lld); %s\n",
                    diff_guides, delta_reads, example[0] == '\0' ? "no example" : example);
            fprintf(stderr, "dotmatch: rerun with --backend cpu for authoritative counts\n");
            free(metal_counts_snapshot);
            free(metal_stats_snapshot);
            goto done;
        }
        metal_validation_status = "passed";
        fprintf(stderr, "dotmatch: Metal validation passed (%zu guides, CPU authority check)\n", targets.count);
        memcpy(counts, metal_counts_snapshot, total_slots * sizeof(unsigned long long));
        memcpy(stats_by_sample, metal_stats_snapshot, reads.count * sizeof(count_stats));
        free(metal_counts_snapshot);
        free(metal_stats_snapshot);
    } else if (metal_validate && backend_mode == COUNT_BACKEND_METAL) {
        fprintf(stderr, "dotmatch: --metal-validate requires an active Metal counting backend for this workload\n");
        goto done;
    }

    counting_seconds = seconds_now() - phase_start_seconds;
    if (show_progress && progress_by_sample != NULL) {
        for (size_t sample = 0; sample < reads.count; ++sample) {
            count_progress_finish(&progress_by_sample[sample]);
        }
    }

    out = open_output_file(out_path);
    if (out == NULL) {
        fprintf(stderr, "failed to open count output\n");
        goto done;
    }
    if (strcmp(format, "mageck") == 0) {
        fprintf(out, "sgRNA\tGene");
        for (size_t sample = 0; sample < reads.count; ++sample) fprintf(out, "\t%s", labels.items[sample]);
        fprintf(out, "\n");
        for (size_t t = 0; t < targets.count; ++t) {
            fprintf(out, "%s\t%s", targets.records[t].id, targets.records[t].gene);
            for (size_t sample = 0; sample < reads.count; ++sample) {
                unsigned long long total = 0;
                for (size_t kind = 0; kind < 5; ++kind) total += counts[((sample * targets.count + t) * 5) + kind];
                fprintf(out, "\t%llu", total);
            }
            fprintf(out, "\n");
        }
    } else {
        fprintf(out, "target_id\ttarget_seq\tgene\tambiguous_nearby");
        for (size_t sample = 0; sample < reads.count; ++sample) {
            fprintf(out, "\t%s_count_exact\t%s_count_corrected_substitution\t%s_count_corrected_insertion\t%s_count_corrected_deletion\t%s_count_corrected_other\t%s_count_total",
                    labels.items[sample], labels.items[sample], labels.items[sample], labels.items[sample], labels.items[sample], labels.items[sample]);
        }
        fprintf(out, "\n");
        for (size_t t = 0; t < targets.count; ++t) {
            fprintf(out, "%s\t%s\t%s\t%d", targets.records[t].id, targets.records[t].seq, targets.records[t].gene, (int)ambiguous_nearby[t]);
            for (size_t sample = 0; sample < reads.count; ++sample) {
                unsigned long long exact = counts[((sample * targets.count + t) * 5) + 0];
                unsigned long long sub = counts[((sample * targets.count + t) * 5) + 1];
                unsigned long long ins = counts[((sample * targets.count + t) * 5) + 2];
                unsigned long long del = counts[((sample * targets.count + t) * 5) + 3];
                unsigned long long other = counts[((sample * targets.count + t) * 5) + 4];
                fprintf(out, "\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu", exact, sub, ins, del, other, exact + sub + ins + del + other);
            }
            fprintf(out, "\n");
        }
    }

    if (target_counts_long_path != NULL) {
        FILE *long_out = open_output_file(target_counts_long_path);
        if (long_out == NULL) {
            fprintf(stderr, "failed to open long target-count output\n");
            goto done;
        }
        fprintf(long_out, "sample_id\ttarget_id\tgroup\tsequence\texact_count\tk1_sub_count\tk1_ins_count\tk1_del_count\tother_count\ttotal_count\tambiguous_nearby\n");
        for (size_t sample = 0; sample < reads.count; ++sample) {
            for (size_t t = 0; t < targets.count; ++t) {
                unsigned long long exact = counts[((sample * targets.count + t) * 5) + 0];
                unsigned long long sub = counts[((sample * targets.count + t) * 5) + 1];
                unsigned long long ins = counts[((sample * targets.count + t) * 5) + 2];
                unsigned long long del = counts[((sample * targets.count + t) * 5) + 3];
                unsigned long long other = counts[((sample * targets.count + t) * 5) + 4];
                fprintf(long_out, "%s\t%s\t%s\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%d\n",
                        labels.items[sample], targets.records[t].id, targets.records[t].gene, targets.records[t].seq,
                        exact, sub, ins, del, other, exact + sub + ins + del + other, (int)ambiguous_nearby[t]);
            }
        }
        fclose(long_out);
    }

    if (sample_qc_path != NULL) {
        FILE *qc = open_output_file(sample_qc_path);
        if (qc == NULL) {
            fprintf(stderr, "failed to open sample QC output\n");
            goto done;
        }
        fprintf(qc, "sample_id\tfastq\ttotal_reads\tvalid_extracted_reads\tassigned_reads\texact_reads\tk1_rescued_reads\tk1_sub_reads\tk1_ins_reads\tk1_del_reads\tambiguous_reads\tno_match_reads\tinvalid_reads\tassignment_rate\texact_rate\trescue_rate\tambiguous_rate\tno_match_rate\ttargets_observed\tzero_count_targets\tgini_index\ttop_1pct_read_fraction\tcandidates_verified\n");
        for (size_t sample = 0; sample < reads.count; ++sample) {
            unsigned long long *target_totals = (unsigned long long *)calloc(targets.count == 0 ? 1 : targets.count, sizeof(unsigned long long));
            if (target_totals == NULL) {
                fclose(qc);
                fprintf(stderr, "out of memory\n");
                goto done;
            }
            unsigned long long sub = 0;
            unsigned long long ins = 0;
            unsigned long long del = 0;
            unsigned long long observed_targets = 0;
            for (size_t t = 0; t < targets.count; ++t) {
                unsigned long long exact = counts[((sample * targets.count + t) * 5) + 0];
                sub += counts[((sample * targets.count + t) * 5) + 1];
                ins += counts[((sample * targets.count + t) * 5) + 2];
                del += counts[((sample * targets.count + t) * 5) + 3];
                target_totals[t] = exact + counts[((sample * targets.count + t) * 5) + 1] +
                                   counts[((sample * targets.count + t) * 5) + 2] +
                                   counts[((sample * targets.count + t) * 5) + 3] +
                                   counts[((sample * targets.count + t) * 5) + 4];
                if (target_totals[t] != 0) ++observed_targets;
            }
            count_stats *s = &stats_by_sample[sample];
            unsigned long long valid = s->total >= s->invalid ? s->total - s->invalid : 0;
            double valid_denom = valid == 0 ? 1.0 : (double)valid;
            fprintf(qc, "%s\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%llu\t%llu\t%.8f\t%.8f\t%llu\n",
                    labels.items[sample], reads.items[sample], s->total, valid, s->unique, s->exact, s->corrected,
                    sub, ins, del, s->ambiguous, s->unmatched, s->invalid,
                    (double)s->unique / valid_denom, (double)s->exact / valid_denom, (double)s->corrected / valid_denom,
                    (double)s->ambiguous / valid_denom, (double)s->unmatched / valid_denom,
                    observed_targets, (unsigned long long)(targets.count - observed_targets),
                    gini_from_counts(target_totals, targets.count),
                    top_fraction_from_counts(target_totals, targets.count, 0.01),
                    s->candidates_verified);
            free(target_totals);
        }
        fclose(qc);
    }

    if (summary_path != NULL) {
        FILE *summary = open_output_file(summary_path);
        if (summary == NULL) {
            fprintf(stderr, "failed to open summary output\n");
            goto done;
        }
        double total_before_summary_seconds = seconds_now() - run_start_seconds;
        fprintf(summary,
                "{\n  \"k\": %d,\n  \"metric\": \"%s\",\n  \"ambiguity_policy\": \"%s\",\n  \"alphabet_policy\": \"%s\",\n  \"max_correction_qual\": ",
                k, metric_name(metric), ambiguity_policy_name(assignment_policy), qdaln_alphabet_policy());
        if (max_correction_qual >= 0) {
            fprintf(summary, "%d", max_correction_qual);
        } else {
            fprintf(summary, "null");
        }
        fprintf(summary,
                ",\n  \"indel_window\": %zu,\n  \"target_start\": %zu,\n  \"auto_offset\": %zu,\n  \"offset_mode\": \"%s\",\n  \"offset_min_fraction\": %.8f,\n  \"offset_detection_strategy\": \"%s\",\n  \"backend_requested\": \"%s\",\n  \"backend_effective\": \"%s\",\n  \"metal_device\": ",
                indel_window, target_start, auto_offset, offset_mode_name(offsets_mode), offset_min_fraction,
                offset_detection_strategy, count_backend_mode_name(backend_mode), backend_effective);
        if (metal_hamming_counts && qdmetal_device_name() != NULL) {
            fprintf(summary, "\"%s\"", qdmetal_device_name());
        } else {
            fprintf(summary, "null");
        }
        fprintf(summary, ",\n  \"metal_validation\": ");
        if (metal_validation_status != NULL) {
            fprintf(summary, "\"%s\"", metal_validation_status);
        } else {
            fprintf(summary, "null");
        }
        fprintf(summary,
                ",\n  \"count_engine\": \"%s\",\n  \"hamming_index\": \"%s\",\n  \"target_length\": %zu,\n  \"n_targets\": %zu,\n  \"read_threads\": %zu,\n  \"phase_seconds\": {\"target_index\": %.6f, \"offset_detection\": %.6f, \"hamming_precompute\": %.6f, \"counting\": %.6f, \"total_before_summary\": %.6f},\n  \"samples\": [\n",
                count_engine, hamming_lookup_kind(&hlookup), target_len, targets.count, effective_read_threads,
                target_index_seconds, offset_detection_seconds, hamming_precompute_seconds, counting_seconds,
                total_before_summary_seconds);
        for (size_t sample = 0; sample < reads.count; ++sample) {
            count_stats *s = &stats_by_sample[sample];
            unsigned long long covered = 0;
            unsigned long long top_count = 0;
            size_t top_target = 0;
            for (size_t t = 0; t < targets.count; ++t) {
                unsigned long long total = 0;
                for (size_t kind = 0; kind < 5; ++kind) total += counts[((sample * targets.count + t) * 5) + kind];
                if (total != 0) ++covered;
                if (total > top_count) {
                    top_count = total;
                    top_target = t;
                }
            }
            double rescued_percent = s->total == 0 ? 0.0 : 100.0 * (double)s->corrected / (double)s->total;
            double ambiguous_percent = s->total == 0 ? 0.0 : 100.0 * (double)s->ambiguous / (double)s->total;
            double unmatched_percent = s->total == 0 ? 0.0 : 100.0 * (double)s->unmatched / (double)s->total;
            fprintf(summary,
                    "    {\"sample\": \"%s\", \"selected_target_start\": %zu, \"selected_target_starts\": [",
                    labels.items[sample], first_selected_offset(&selected_offsets[sample], target_start));
            for (size_t oi = 0; oi < selected_offsets[sample].count; ++oi) {
                if (oi != 0) fprintf(summary, ", ");
                fprintf(summary, "%zu", selected_offsets[sample].items[oi]);
            }
            fprintf(summary,
                    "], \"total_reads\": %llu, \"assigned_unique\": %llu, \"assigned_exact\": %llu, \"assigned_corrected\": %llu, \"k1_rescued_reads\": %llu, \"percent_rescued_by_k1\": %.6f, \"ambiguous\": %llu, \"percent_ambiguous\": %.6f, \"unmatched\": %llu, \"percent_unmatched\": %.6f, \"invalid\": %llu, \"library_covered_targets\": %llu, \"library_coverage_fraction\": %.6f, \"top_target_id\": \"%s\", \"top_target_count\": %llu, \"candidates_considered\": %llu, \"candidates_verified\": %llu}%s\n",
                    s->total, s->unique, s->exact, s->corrected,
                    s->corrected, rescued_percent, s->ambiguous, ambiguous_percent, s->unmatched, unmatched_percent,
                    s->invalid,
                    covered, targets.count == 0 ? 0.0 : (double)covered / (double)targets.count,
                    targets.count == 0 ? "" : targets.records[top_target].id, top_count, s->candidates_considered,
                    s->candidates_verified, sample + 1 == reads.count ? "" : ",");
        }
        fprintf(summary, "  ]\n}\n");
        fclose(summary);
    }
    if (report_path != NULL) {
        if (write_count_html_report(report_path, &targets, &reads, &labels, counts, stats_by_sample, selected_offsets,
                                    k, metric, assignment_policy, target_len, report_audit_dir,
                                    report_unmatched_path) != 0) {
            fprintf(stderr, "failed to write HTML report\n");
            goto done;
        }
    }
    {
        sample_qc_metrics *qc_metrics = (sample_qc_metrics *)calloc(reads.count, sizeof(sample_qc_metrics));
        if (qc_metrics != NULL) {
            int metrics_ok = 1;
            for (size_t sample = 0; sample < reads.count; ++sample) {
                if (compute_sample_qc_metrics(&targets, counts, sample, &stats_by_sample[sample],
                                              &qc_metrics[sample]) != 0) {
                    metrics_ok = 0;
                    break;
                }
            }
            if (metrics_ok) {
                int enough_reads = 0;
                for (size_t sample = 0; sample < reads.count; ++sample) {
                    if (stats_by_sample[sample].total >= 1000) enough_reads = 1;
                }
                if (enough_reads || crispr_mode) emit_sample_qc_review_warnings(&labels, qc_metrics, reads.count);
            }
            free(qc_metrics);
        }
    }
    rc = 0;

done:
    if (progress_by_sample != NULL) {
        for (size_t sample = 0; sample < reads.count; ++sample) count_progress_fini(&progress_by_sample[sample]);
        free(progress_by_sample);
    }
    if (out != NULL) fclose(out);
    if (assignments != NULL) fclose(assignments);
    if (ambiguous_out != NULL) fclose(ambiguous_out);
    if (unmatched_out != NULL) fclose(unmatched_out);
    qdaln_index_free(index);
    free_hamming_lookup(&hlookup);
    free_hamming_lookup(&offset_lookup);
    free_levenshtein1_lookup(&levlookup);
    free(metal_target_codes);
    free(target_ptrs);
    free(target_lens);
    free(ambiguous_nearby);
    free(counts);
    free(stats_by_sample);
    if (selected_offsets != NULL) {
        for (size_t sample = 0; sample < reads.count; ++sample) free_offset_list(&selected_offsets[sample]);
    }
    free(selected_offsets);
    free_table(&targets);
    free_string_list(&reads);
    free_string_list(&labels);
    return rc;

fail_args:
    free_string_list(&reads);
    free_string_list(&labels);
    return 2;
}

static int string_list_contains_exact(const string_list *list, const char *s) {
    if (list == NULL || s == NULL) return 0;
    for (size_t i = 0; i < list->count; ++i) {
        if (strcmp(list->items[i], s) == 0) return 1;
    }
    return 0;
}

static int read_first_column_values(const char *path, string_list *values) {
    if (path == NULL) return 0;
    fastq_reader reader = {0};
    if (fastq_reader_open(&reader, path) != 0) return -1;
    char buf[8192];
    size_t len = 0;
    int rc = 0;
    while ((rc = fastq_getline_len(&reader, buf, sizeof(buf), &len)) > 0) {
        (void)len;
        trim_line(buf);
        if (buf[0] == '\0' || buf[0] == '#') continue;
        char *tab = strchr(buf, '\t');
        if (tab != NULL) *tab = '\0';
        if (push_string(values, buf) != 0) {
            fastq_reader_close(&reader);
            return -1;
        }
    }
    fastq_reader_close(&reader);
    return rc < 0 ? -1 : 0;
}

static int guide_counter_push_sample_name(string_list *labels, const char *path, size_t idx) {
    const char *base = path_basename(path);
    char fallback[32];
    if (base == NULL || base[0] == '\0') {
        int n = snprintf(fallback, sizeof(fallback), "s%zu", idx + 1);
        if (n < 0 || (size_t)n >= sizeof(fallback)) return -1;
        return push_string(labels, fallback);
    }
    char *name = xstrndup(base, strlen(base));
    if (name == NULL) return -1;
    if (ends_with(name, ".gz")) name[strlen(name) - 3] = '\0';
    if (ends_with(name, ".fastq")) {
        name[strlen(name) - 6] = '\0';
    } else if (ends_with(name, ".fq")) {
        name[strlen(name) - 3] = '\0';
    }
    int rc = push_string(labels, name);
    free(name);
    return rc;
}

static const char *guide_counter_type_for_target(const seq_record *target, const string_list *essential_genes,
                                                 const string_list *nonessential_genes,
                                                 const string_list *control_guides,
                                                 regex_t *control_re) {
    if (string_list_contains_exact(essential_genes, target->gene)) return "Essential";
    if (string_list_contains_exact(nonessential_genes, target->gene)) return "Nonessential";
    if (string_list_contains_exact(control_guides, target->id)) return "Control";
    if (control_re != NULL &&
        (regexec(control_re, target->id, 0, NULL, 0) == 0 ||
         regexec(control_re, target->gene, 0, NULL, 0) == 0)) {
        return "Control";
    }
    return "Other";
}

static int parse_ull_value(const char *s, unsigned long long *out) {
    if (s == NULL || s[0] == '-' || s[0] == '\0') return -1;
    char *end = NULL;
    errno = 0;
    unsigned long long v = strtoull(s, &end, 10);
    if (errno == ERANGE || end == s || *end != '\0') return -1;
    *out = v;
    return 0;
}

static double round_positive_dp(double value, int places) {
    double factor = 1.0;
    for (int i = 0; i < places; ++i) factor *= 10.0;
    unsigned long long scaled = (unsigned long long)(value * factor + 0.5);
    return (double)scaled / factor;
}

static int guide_counter_write_outputs(const char *output_prefix, const char *tmp_counts_path,
                                       const char *tmp_qc_path, const seq_table *targets,
                                       const string_list *reads, const string_list *labels,
                                       const string_list *essential_genes,
                                       const string_list *nonessential_genes,
                                       const string_list *control_guides, regex_t *control_re) {
    char counts_path[4096];
    char extended_path[4096];
    char stats_path[4096];
    int n = snprintf(counts_path, sizeof(counts_path), "%s.counts.txt", output_prefix);
    if (n < 0 || (size_t)n >= sizeof(counts_path)) return -1;
    n = snprintf(extended_path, sizeof(extended_path), "%s.extended-counts.txt", output_prefix);
    if (n < 0 || (size_t)n >= sizeof(extended_path)) return -1;
    n = snprintf(stats_path, sizeof(stats_path), "%s.stats.txt", output_prefix);
    if (n < 0 || (size_t)n >= sizeof(stats_path)) return -1;

    size_t matrix_slots = 0;
    if (checked_mul_size(targets->count, labels->count, &matrix_slots) != 0) return -1;
    unsigned long long *matrix = (unsigned long long *)calloc(alloc_count_or_one(matrix_slots),
                                                             sizeof(unsigned long long));
    const char **types = (const char **)calloc(alloc_count_or_one(targets->count), sizeof(const char *));
    if (matrix == NULL || types == NULL) {
        free(matrix);
        free(types);
        return -1;
    }
    for (size_t t = 0; t < targets->count; ++t) {
        types[t] = guide_counter_type_for_target(&targets->records[t], essential_genes, nonessential_genes,
                                                 control_guides, control_re);
    }

    FILE *in = fopen(tmp_counts_path, "r");
    FILE *counts = open_output_file(counts_path);
    FILE *extended = open_output_file(extended_path);
    if (in == NULL || counts == NULL || extended == NULL) {
        if (in != NULL) fclose(in);
        if (counts != NULL) fclose(counts);
        if (extended != NULL) fclose(extended);
        free(matrix);
        free(types);
        return -1;
    }

    fprintf(counts, "guide\tgene");
    fprintf(extended, "guide\tgene\tguide_type");
    for (size_t sample = 0; sample < labels->count; ++sample) {
        fprintf(counts, "\t%s", labels->items[sample]);
        fprintf(extended, "\t%s", labels->items[sample]);
    }
    fprintf(counts, "\n");
    fprintf(extended, "\n");

    char buf[65536];
    size_t row = 0;
    int first = 1;
    while (fgets(buf, sizeof(buf), in) != NULL) {
        trim_line(buf);
        if (first) {
            first = 0;
            continue;
        }
        char *fields[1024];
        size_t nf = split_fields(buf, '\t', fields, sizeof(fields) / sizeof(fields[0]));
        if (nf < 2 + labels->count || row >= targets->count) {
            fclose(in);
            fclose(counts);
            fclose(extended);
            free(matrix);
            free(types);
            return -1;
        }
        fprintf(counts, "%s\t%s", targets->records[row].id, targets->records[row].gene);
        fprintf(extended, "%s\t%s\t%s", targets->records[row].id, targets->records[row].gene, types[row]);
        for (size_t sample = 0; sample < labels->count; ++sample) {
            unsigned long long value = 0;
            if (parse_ull_value(fields[2 + sample], &value) != 0) {
                fclose(in);
                fclose(counts);
                fclose(extended);
                free(matrix);
                free(types);
                return -1;
            }
            matrix[row * labels->count + sample] = value;
            fprintf(counts, "\t%llu", value);
            fprintf(extended, "\t%llu", value);
        }
        fprintf(counts, "\n");
        fprintf(extended, "\n");
        ++row;
    }
    int matrix_ok = !ferror(in) && row == targets->count;
    fclose(in);
    fclose(counts);
    fclose(extended);
    if (!matrix_ok) {
        free(matrix);
        free(types);
        return -1;
    }

    unsigned long long *total_reads = (unsigned long long *)calloc(alloc_count_or_one(labels->count),
                                                                   sizeof(unsigned long long));
    if (total_reads == NULL) {
        free(matrix);
        free(types);
        return -1;
    }
    FILE *qc = fopen(tmp_qc_path, "r");
    if (qc == NULL) {
        free(total_reads);
        free(matrix);
        free(types);
        return -1;
    }
    int total_reads_col = -1;
    size_t qc_row = 0;
    first = 1;
    while (fgets(buf, sizeof(buf), qc) != NULL) {
        trim_line(buf);
        char *fields[64];
        size_t nf = split_fields(buf, '\t', fields, sizeof(fields) / sizeof(fields[0]));
        if (first) {
            total_reads_col = find_column(fields, nf, "total_reads", NULL, NULL);
            first = 0;
            continue;
        }
        if (total_reads_col < 0 || (size_t)total_reads_col >= nf || qc_row >= labels->count ||
            parse_ull_value(fields[total_reads_col], &total_reads[qc_row]) != 0) {
            fclose(qc);
            free(total_reads);
            free(matrix);
            free(types);
            return -1;
        }
        ++qc_row;
    }
    int qc_ok = !ferror(qc) && qc_row == labels->count;
    fclose(qc);
    if (!qc_ok) {
        free(total_reads);
        free(matrix);
        free(types);
        return -1;
    }

    FILE *stats = open_output_file(stats_path);
    if (stats == NULL) {
        free(total_reads);
        free(matrix);
        free(types);
        return -1;
    }
    fprintf(stats, "file\tlabel\ttotal_guides\ttotal_reads\tmapped_reads\tfrac_mapped\tmean_reads_per_guide\tmean_reads_essential\tmean_reads_nonessential\tmean_reads_control\tmean_reads_other\tzero_read_guides\n");
    for (size_t sample = 0; sample < labels->count; ++sample) {
        unsigned long long mapped = 0;
        unsigned long long zero = 0;
        double essential_sum = 0.0;
        double nonessential_sum = 0.0;
        double control_sum = 0.0;
        double other_sum = 0.0;
        size_t essential_count = 0;
        size_t nonessential_count = 0;
        size_t control_count = 0;
        size_t other_count = 0;
        for (size_t t = 0; t < targets->count; ++t) {
            unsigned long long value = matrix[t * labels->count + sample];
            mapped += value;
            if (value == 0) ++zero;
            if (strcmp(types[t], "Essential") == 0) {
                essential_sum += (double)value;
                ++essential_count;
            } else if (strcmp(types[t], "Nonessential") == 0) {
                nonessential_sum += (double)value;
                ++nonessential_count;
            } else if (strcmp(types[t], "Control") == 0) {
                control_sum += (double)value;
                ++control_count;
            } else {
                other_sum += (double)value;
                ++other_count;
            }
        }
        double total = (double)total_reads[sample];
        double frac = total == 0.0 ? 0.0 : round_positive_dp((double)mapped / total, 4);
        double mean_all = targets->count == 0 ? 0.0 : round_positive_dp((double)mapped / (double)targets->count, 2);
        double mean_essential = essential_count == 0 ? 0.0 : round_positive_dp(essential_sum / (double)essential_count, 2);
        double mean_nonessential = nonessential_count == 0 ? 0.0 : round_positive_dp(nonessential_sum / (double)nonessential_count, 2);
        double mean_control = control_count == 0 ? 0.0 : round_positive_dp(control_sum / (double)control_count, 2);
        double mean_other = other_count == 0 ? 0.0 : round_positive_dp(other_sum / (double)other_count, 2);
        fprintf(stats, "%s\t%s\t%zu\t%llu\t%llu\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%llu\n",
                reads->items[sample], labels->items[sample], targets->count, total_reads[sample], mapped,
                frac, mean_all, mean_essential, mean_nonessential, mean_control, mean_other, zero);
    }
    fclose(stats);
    free(total_reads);
    free(matrix);
    free(types);
    return 0;
}

static int push_count_arg(string_list *args, const char *s) {
    return push_string(args, s);
}

static int run_guide_counter_compatible(const char *argv0, int argc, char **argv) {
    int start = 2;
    if (strcmp(argv[1], "guide-counter") == 0) {
        if (argc < 3 || strcmp(argv[2], "count") != 0) {
            usage(argv0);
            return 2;
        }
        start = 3;
    }

    const char *library_path = NULL;
    const char *output_prefix = NULL;
    const char *essential_path = NULL;
    const char *nonessential_path = NULL;
    const char *control_guides_path = NULL;
    const char *control_pattern = NULL;
    size_t offset_sample_size = 100000;
    double offset_min_fraction = 0.0025;
    int exact_match = 0;
    string_list reads = {0};
    string_list labels = {0};

    int i = start;
    while (i < argc) {
        const char *arg = argv[i++];
        if ((strcmp(arg, "--input") == 0 || strcmp(arg, "-i") == 0) && i < argc) {
            while (i < argc && argv[i][0] != '-') {
                if (push_string(&reads, argv[i++]) != 0) goto oom;
            }
        } else if ((strcmp(arg, "--samples") == 0 || strcmp(arg, "-s") == 0) && i < argc) {
            while (i < argc && argv[i][0] != '-') {
                if (push_string(&labels, argv[i++]) != 0) goto oom;
            }
        } else if ((strcmp(arg, "--library") == 0 || strcmp(arg, "-l") == 0) && i < argc) {
            library_path = argv[i++];
        } else if ((strcmp(arg, "--output") == 0 || strcmp(arg, "-o") == 0) && i < argc) {
            output_prefix = argv[i++];
        } else if ((strcmp(arg, "--essential-genes") == 0 || strcmp(arg, "-e") == 0) && i < argc) {
            essential_path = argv[i++];
        } else if ((strcmp(arg, "--nonessential-genes") == 0 || strcmp(arg, "-n") == 0) && i < argc) {
            nonessential_path = argv[i++];
        } else if ((strcmp(arg, "--control-guides") == 0 || strcmp(arg, "-c") == 0) && i < argc) {
            control_guides_path = argv[i++];
        } else if ((strcmp(arg, "--control-pattern") == 0 || strcmp(arg, "-C") == 0) && i < argc) {
            control_pattern = argv[i++];
        } else if ((strcmp(arg, "--offset-sample-size") == 0 || strcmp(arg, "-N") == 0) && i < argc) {
            if (parse_size_value(argv[i++], &offset_sample_size) != 0 || offset_sample_size == 0) goto bad_args;
        } else if ((strcmp(arg, "--offset-min-fraction") == 0 || strcmp(arg, "-f") == 0) && i < argc) {
            if (parse_double_value(argv[i++], &offset_min_fraction) != 0 ||
                offset_min_fraction < 0.0 || offset_min_fraction > 1.0) {
                goto bad_args;
            }
        } else if (strcmp(arg, "--exact-match") == 0 || strcmp(arg, "-x") == 0) {
            exact_match = 1;
        } else if (strcmp(arg, "--help") == 0 || strcmp(arg, "-h") == 0) {
            usage(argv0);
            free_string_list(&reads);
            free_string_list(&labels);
            return 0;
        } else {
            goto bad_args;
        }
    }

    if (library_path == NULL || output_prefix == NULL || reads.count == 0) goto bad_args;
    if (labels.count == 0) {
        for (size_t sample = 0; sample < reads.count; ++sample) {
            if (guide_counter_push_sample_name(&labels, reads.items[sample], sample) != 0) goto oom;
        }
    }
    if (labels.count != reads.count) {
        fprintf(stderr, "--samples count must match --input count\n");
        free_string_list(&reads);
        free_string_list(&labels);
        return 2;
    }
    if (validate_unique_sample_labels(&labels, "--samples") != 0) {
        free_string_list(&reads);
        free_string_list(&labels);
        return 2;
    }

    seq_table targets = {0};
    string_list essential_genes = {0};
    string_list nonessential_genes = {0};
    string_list control_guides = {0};
    regex_t control_re;
    int have_control_re = 0;
    string_list count_args = {0};
    int rc = 1;
    char target_len_s[32];
    char k_s[8];
    char auto_offset_s[16];
    char offset_sample_s[32];
    char offset_min_s[64];
    char label_csv[8192];
    char tmp_counts_path[4096] = "";
    char tmp_qc_path[4096] = "";

    if (read_target_table(library_path, &targets) != 0 || targets.count == 0) {
        fprintf(stderr, "failed to read guide library\n");
        goto done;
    }
    int guide_id_check = validate_unique_seq_ids(&targets, "guide");
    if (guide_id_check != 0) {
        if (guide_id_check == -1) fprintf(stderr, "out of memory\n");
        goto done;
    }
    size_t guide_len = targets.records[0].len;
    for (size_t t = 1; t < targets.count; ++t) {
        if (targets.records[t].len != guide_len) {
            fprintf(stderr, "GuideCounter compatibility requires one guide length\n");
            goto done;
        }
    }
    if (read_first_column_values(essential_path, &essential_genes) != 0 ||
        read_first_column_values(nonessential_path, &nonessential_genes) != 0 ||
        read_first_column_values(control_guides_path, &control_guides) != 0) {
        fprintf(stderr, "failed to read guide annotation files\n");
        goto done;
    }
    if (control_pattern != NULL) {
        if (regcomp(&control_re, control_pattern, REG_EXTENDED | REG_ICASE | REG_NOSUB) != 0) {
            fprintf(stderr, "failed to compile --control-pattern\n");
            goto done;
        }
        have_control_re = 1;
    }

    label_csv[0] = '\0';
    for (size_t sample = 0; sample < labels.count; ++sample) {
        size_t used = strlen(label_csv);
        int n = snprintf(label_csv + used, sizeof(label_csv) - used, "%s%s",
                         sample == 0 ? "" : ",", labels.items[sample]);
        if (n < 0 || (size_t)n >= sizeof(label_csv) - used) {
            fprintf(stderr, "too many sample labels for GuideCounter compatibility wrapper\n");
            goto done;
        }
    }
    int n = snprintf(tmp_counts_path, sizeof(tmp_counts_path), "%s.dotmatch-counts.tmp", output_prefix);
    if (n < 0 || (size_t)n >= sizeof(tmp_counts_path)) goto done;
    n = snprintf(tmp_qc_path, sizeof(tmp_qc_path), "%s.dotmatch-qc.tmp", output_prefix);
    if (n < 0 || (size_t)n >= sizeof(tmp_qc_path)) goto done;
    snprintf(target_len_s, sizeof(target_len_s), "%zu", guide_len);
    snprintf(k_s, sizeof(k_s), "%d", exact_match ? 0 : 1);
    snprintf(auto_offset_s, sizeof(auto_offset_s), "%d", 499);
    snprintf(offset_sample_s, sizeof(offset_sample_s), "%zu", offset_sample_size);
    snprintf(offset_min_s, sizeof(offset_min_s), "%.8g", offset_min_fraction);

    if (push_count_arg(&count_args, argv0) != 0 ||
        push_count_arg(&count_args, "count") != 0 ||
        push_count_arg(&count_args, "--targets") != 0 ||
        push_count_arg(&count_args, library_path) != 0) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    for (size_t sample = 0; sample < reads.count; ++sample) {
        if (push_count_arg(&count_args, "--reads") != 0 ||
            push_count_arg(&count_args, reads.items[sample]) != 0) {
            fprintf(stderr, "out of memory\n");
            goto done;
        }
    }
    const char *fixed_args[] = {
        "--sample-label", label_csv,
        "--target-start", "0",
        "--target-length", target_len_s,
        "--k", k_s,
        "--metric", "hamming",
        "--ambiguity-policy", "best",
        "--format", "mageck",
        "--auto-offset", auto_offset_s,
        "--auto-offset-sample", offset_sample_s,
        "--offset-mode", "multi",
        "--offset-min-fraction", offset_min_s,
        "--out", tmp_counts_path,
        "--sample-qc", tmp_qc_path
    };
    for (size_t ai = 0; ai < sizeof(fixed_args) / sizeof(fixed_args[0]); ++ai) {
        if (push_count_arg(&count_args, fixed_args[ai]) != 0) {
            fprintf(stderr, "out of memory\n");
            goto done;
        }
    }

    rc = run_count(argv0, (int)count_args.count, count_args.items);
    if (rc != 0) goto done;
    if (guide_counter_write_outputs(output_prefix, tmp_counts_path, tmp_qc_path, &targets, &reads, &labels,
                                    &essential_genes, &nonessential_genes, &control_guides,
                                    have_control_re ? &control_re : NULL) != 0) {
        fprintf(stderr, "failed to write GuideCounter-compatible outputs\n");
        rc = 1;
        goto done;
    }
    rc = 0;

done:
    if (have_control_re) regfree(&control_re);
    unlink(tmp_counts_path);
    unlink(tmp_qc_path);
    free_string_list(&count_args);
    free_string_list(&essential_genes);
    free_string_list(&nonessential_genes);
    free_string_list(&control_guides);
    free_table(&targets);
    free_string_list(&reads);
    free_string_list(&labels);
    return rc;

oom:
    fprintf(stderr, "out of memory\n");
    free_string_list(&reads);
    free_string_list(&labels);
    return 1;

bad_args:
    usage(argv0);
    free_string_list(&reads);
    free_string_list(&labels);
    return 2;
}

static int run_fastq_assign(const char *argv0, int argc, char **argv) {
    const char *barcodes_path = NULL;
    const char *reads_path = NULL;
    const char *out_path = NULL;
    ambiguity_policy assignment_policy = AMBIGUITY_POLICY_RADIUS;
    size_t barcode_start = 0;
    size_t barcode_len = 0;
    int k = -1;

    int i = 2;
    while (i < argc) {
        const char *arg = argv[i++];
        if (strcmp(arg, "--barcodes") == 0 && i < argc) {
            barcodes_path = argv[i++];
        } else if (strcmp(arg, "--reads") == 0 && i < argc) {
            reads_path = argv[i++];
        } else if (strcmp(arg, "--barcode-start") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &barcode_start) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--barcode-length") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &barcode_len) != 0 || barcode_len == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--k") == 0 && i < argc) {
            if (parse_int_value(argv[i++], &k) != 0 || (k != 0 && k != 1)) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--ambiguity-policy") == 0 && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "best") == 0) {
                assignment_policy = AMBIGUITY_POLICY_BEST;
            } else if (strcmp(value, "radius") == 0) {
                assignment_policy = AMBIGUITY_POLICY_RADIUS;
            } else {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--out") == 0 && i < argc) {
            out_path = argv[i++];
        } else {
            usage(argv0);
            return 2;
        }
    }

    if (barcodes_path == NULL || reads_path == NULL || out_path == NULL || barcode_len == 0 || k < 0) {
        usage(argv0);
        return 2;
    }

    seq_table targets = {0};
    fastq_reader reader = {0};
    FILE *out = NULL;
    qdaln_index *index = NULL;
    int rc = 1;

    if (read_table(barcodes_path, &targets) != 0) {
        fprintf(stderr, "failed to read barcode file\n");
        goto done;
    }
    int barcode_id_check = validate_unique_seq_ids(&targets, "barcode");
    if (barcode_id_check != 0) {
        if (barcode_id_check == -1) fprintf(stderr, "out of memory\n");
        goto done;
    }

    const char **target_ptrs = (const char **)malloc(targets.count * sizeof(char *));
    size_t *target_lens = (size_t *)malloc(targets.count * sizeof(size_t));
    if (targets.count != 0 && (target_ptrs == NULL || target_lens == NULL)) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    for (size_t i = 0; i < targets.count; ++i) {
        target_ptrs[i] = targets.records[i].seq;
        target_lens[i] = targets.records[i].len;
    }
    index = qdaln_index_build(target_ptrs, target_lens, targets.count);
    free(target_ptrs);
    free(target_lens);
    if (index == NULL) {
        fprintf(stderr, "failed to build barcode index\n");
        goto done;
    }

    if (fastq_reader_open(&reader, reads_path) != 0) {
        fprintf(stderr, "failed to open FASTQ input\n");
        goto done;
    }
    out = open_output_file(out_path);
    if (out == NULL) {
        fprintf(stderr, "failed to open output file\n");
        goto done;
    }

    fprintf(out, "read_id\tobserved_barcode\ttarget_index\ttarget_id\ttarget_seq\tbest_distance\tsecond_best_distance\tmatch_count\tstatus\n");

    char header[8192];
    char seq[8192];
    char plus[8192];
    char qual[8192];
    char read_id[8192];
    char observed[8192];
    int got = 0;
    size_t seq_len = 0;
    while ((got = fastq_read_record_len(&reader, header, seq, plus, qual, sizeof(header), &seq_len)) == 1) {
        fastq_read_id(header, read_id, sizeof(read_id));
        qdaln_match_result result = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
        observed[0] = '\0';
        if (barcode_start <= seq_len && barcode_len <= seq_len - barcode_start && barcode_len < sizeof(observed)) {
            memcpy(observed, seq + barcode_start, barcode_len);
            observed[barcode_len] = '\0';
            const char *read_ptr = observed;
            size_t read_len = barcode_len;
            qdaln_index_stats stats;
            if (qdaln_index_assign_stats(index, &read_ptr, &read_len, 1, k, &result, &stats) != 0) {
                fprintf(stderr, "FASTQ assignment failed\n");
                goto done;
            }
            apply_ambiguity_policy(&result, assignment_policy);
        }
        print_fastq_row(out, &targets, read_id, observed, result);
    }
    if (got < 0) {
        fprintf(stderr, "malformed FASTQ input\n");
        goto done;
    }
    rc = 0;

done:
    if (out != NULL) fclose(out);
    fastq_reader_close(&reader);
    qdaln_index_free(index);
    free_table(&targets);
    return rc;
}

static void sanitize_filename(const char *in, char *out, size_t out_cap) {
    size_t j = 0;
    if (out_cap == 0) return;
    for (size_t i = 0; in[i] != '\0' && j + 1 < out_cap; ++i) {
        char c = in[i];
        if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ||
            c == '_' || c == '-' || c == '.') {
            out[j++] = c;
        } else {
            out[j++] = '_';
        }
    }
    if (j == 0 && out_cap > 1) out[j++] = '_';
    out[j] = '\0';
}

typedef struct sanitized_name_entry {
    char *name;
    const char *id;
} sanitized_name_entry;

static int compare_sanitized_name_entry(const void *a, const void *b) {
    const sanitized_name_entry *ea = (const sanitized_name_entry *)a;
    const sanitized_name_entry *eb = (const sanitized_name_entry *)b;
    return strcmp(ea->name, eb->name);
}

static int validate_unique_sanitized_filenames(const seq_table *targets) {
    sanitized_name_entry *entries = (sanitized_name_entry *)calloc(targets->count == 0 ? 1 : targets->count,
                                                                   sizeof(sanitized_name_entry));
    if (entries == NULL) return -1;

    int rc = 0;
    for (size_t i = 0; i < targets->count; ++i) {
        char safe_id[512];
        sanitize_filename(targets->records[i].id, safe_id, sizeof(safe_id));
        entries[i].name = xstrndup(safe_id, strlen(safe_id));
        entries[i].id = targets->records[i].id;
        if (entries[i].name == NULL) {
            rc = -1;
            goto done;
        }
    }

    qsort(entries, targets->count, sizeof(entries[0]), compare_sanitized_name_entry);
    for (size_t i = 1; i < targets->count; ++i) {
        if (strcmp(entries[i - 1].name, entries[i].name) == 0) {
            fprintf(stderr,
                    "barcode IDs produce the same output filename after sanitization: \"%s\" and \"%s\" -> %s.fastq\n",
                    entries[i - 1].id, entries[i].id, entries[i].name);
            rc = -2;
            goto done;
        }
    }

done:
    for (size_t i = 0; i < targets->count; ++i) free(entries[i].name);
    free(entries);
    return rc;
}

static int ensure_dir(const char *path) {
    if (mkdir(path, 0777) == 0) return 0;
    if (errno == EEXIST) {
        struct stat st;
        return stat(path, &st) == 0 && S_ISDIR(st.st_mode) ? 0 : -1;
    }
    return -1;
}

static int path_join(char *out, size_t out_cap, const char *dir, const char *name) {
    int n = snprintf(out, out_cap, "%s/%s", dir, name);
    return n < 0 || (size_t)n >= out_cap ? -1 : 0;
}

static size_t uf_find(size_t *parent, size_t x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static void uf_union(size_t *parent, size_t a, size_t b) {
    size_t ra = uf_find(parent, a);
    size_t rb = uf_find(parent, b);
    if (ra == rb) return;
    if (ra < rb) parent[rb] = ra;
    else parent[ra] = rb;
}

static int string_list_contains(const string_list *list, const char *s) {
    for (size_t i = 0; i < list->count; ++i) {
        if (strcmp(list->items[i], s) == 0) return 1;
    }
    return 0;
}

static int push_unique_string(string_list *list, const char *s) {
    if (string_list_contains(list, s)) return 0;
    return push_string(list, s);
}

static int add_k1_variants_for_target(string_list *variants, const char *seq, size_t len) {
    static const char dna[] = "ACGT";
    if (push_unique_string(variants, seq) != 0) return -1;
    char buf[8192];
    if (len + 2 > sizeof(buf)) return -1;

    for (size_t pos = 0; pos < len; ++pos) {
        for (size_t bi = 0; bi < 4; ++bi) {
            if (seq[pos] == dna[bi]) continue;
            memcpy(buf, seq, len);
            buf[pos] = dna[bi];
            buf[len] = '\0';
            if (push_unique_string(variants, buf) != 0) return -1;
        }
    }

    if (len > 0) {
        for (size_t pos = 0; pos < len; ++pos) {
            memcpy(buf, seq, pos);
            memcpy(buf + pos, seq + pos + 1, len - pos - 1);
            buf[len - 1] = '\0';
            if (push_unique_string(variants, buf) != 0) return -1;
        }
    }

    for (size_t pos = 0; pos <= len; ++pos) {
        for (size_t bi = 0; bi < 4; ++bi) {
            memcpy(buf, seq, pos);
            buf[pos] = dna[bi];
            memcpy(buf + pos + 1, seq + pos, len - pos);
            buf[len + 1] = '\0';
            if (push_unique_string(variants, buf) != 0) return -1;
        }
    }
    return 0;
}

typedef struct variant_record {
    char *key;
    size_t target;
} variant_record;

typedef struct variant_record_list {
    variant_record *items;
    size_t count;
    size_t cap;
} variant_record_list;

static void free_variant_record_list(variant_record_list *list) {
    for (size_t i = 0; i < list->count; ++i) free(list->items[i].key);
    free(list->items);
    list->items = NULL;
    list->count = 0;
    list->cap = 0;
}

static int push_variant_record(variant_record_list *list, const char *key, size_t target) {
    if (list->count == list->cap) {
        size_t next_cap = list->cap == 0 ? 1024 : list->cap * 2;
        variant_record *next = (variant_record *)realloc(list->items, next_cap * sizeof(variant_record));
        if (next == NULL) return -1;
        list->items = next;
        list->cap = next_cap;
    }
    list->items[list->count].key = xstrndup(key, strlen(key));
    if (list->items[list->count].key == NULL) return -1;
    list->items[list->count].target = target;
    ++list->count;
    return 0;
}

static int cmp_variant_record(const void *a, const void *b) {
    const variant_record *aa = (const variant_record *)a;
    const variant_record *bb = (const variant_record *)b;
    int c = strcmp(aa->key, bb->key);
    if (c != 0) return c;
    return aa->target > bb->target ? 1 : (aa->target < bb->target ? -1 : 0);
}

typedef struct pair_record {
    size_t a;
    size_t b;
} pair_record;

typedef struct pair_record_list {
    pair_record *items;
    size_t count;
    size_t cap;
} pair_record_list;

static void free_pair_record_list(pair_record_list *list) {
    free(list->items);
    list->items = NULL;
    list->count = 0;
    list->cap = 0;
}

static int cmp_pair_record(const void *a, const void *b) {
    const pair_record *aa = (const pair_record *)a;
    const pair_record *bb = (const pair_record *)b;
    if (aa->a != bb->a) return aa->a > bb->a ? 1 : -1;
    return aa->b > bb->b ? 1 : (aa->b < bb->b ? -1 : 0);
}

static int push_pair_record(pair_record_list *list, size_t a, size_t b) {
    if (a > b) {
        size_t tmp = a;
        a = b;
        b = tmp;
    }
    if (list->count == list->cap) {
        size_t next_cap = list->cap == 0 ? 1024 : list->cap * 2;
        pair_record *next = (pair_record *)realloc(list->items, next_cap * sizeof(pair_record));
        if (next == NULL) return -1;
        list->items = next;
        list->cap = next_cap;
    }
    list->items[list->count++] = (pair_record){a, b};
    return 0;
}

typedef struct seq_ref {
    const char *seq;
    size_t len;
} seq_ref;

static int cmp_seq_ref(const void *a, const void *b) {
    const seq_ref *aa = (const seq_ref *)a;
    const seq_ref *bb = (const seq_ref *)b;
    size_t min_len = aa->len < bb->len ? aa->len : bb->len;
    int c = memcmp(aa->seq, bb->seq, min_len);
    if (c != 0) return c;
    return aa->len > bb->len ? 1 : (aa->len < bb->len ? -1 : 0);
}

static size_t count_unique_target_sequences(const seq_table *targets) {
    if (targets->count == 0) return 0;
    seq_ref *refs = (seq_ref *)malloc(targets->count * sizeof(seq_ref));
    if (refs == NULL) return 0;
    for (size_t i = 0; i < targets->count; ++i) {
        refs[i].seq = targets->records[i].seq;
        refs[i].len = targets->records[i].len;
    }
    qsort(refs, targets->count, sizeof(seq_ref), cmp_seq_ref);
    size_t unique = 1;
    for (size_t i = 1; i < targets->count; ++i) {
        if (refs[i].len != refs[i - 1].len || memcmp(refs[i].seq, refs[i - 1].seq, refs[i].len) != 0) {
            ++unique;
        }
    }
    free(refs);
    return unique;
}

static int write_audit_summary_json(const char *out_dir, const char *audit_mode, int k,
                                    size_t n_targets, size_t unique_sequences,
                                    const char *min_edit_distance_json, const char *min_hamming_distance_json,
                                    int safe_at_k0, int safe_at_k1, const char *safe_at_k2_json,
                                    const char *safe_at_hamming_k2_json, const char *safe_at_hamming_k3_json,
                                    unsigned long long pairs_d0, unsigned long long pairs_d1,
                                    unsigned long long pairs_d2, unsigned long long pairs_within_k,
                                    unsigned long long risk_pairs_k1, const char *risk_pairs_k2_json,
                                    const char *risk_pairs_hamming_k2_json, const char *risk_pairs_hamming_k3_json,
                                    unsigned long long ambiguous_query_variants_k1, int recommended_k) {
    char path[4096];
    if (path_join(path, sizeof(path), out_dir, "audit_summary.json") != 0) return -1;
    FILE *out = open_output_file(path);
    if (out == NULL) return -1;
    fprintf(out,
            "{\n"
            "  \"audit_mode\": \"%s\",\n"
            "  \"k\": %d,\n"
            "  \"targets\": %zu,\n"
            "  \"unique_sequences\": %zu,\n"
            "  \"duplicate_sequences\": %zu,\n"
            "  \"min_edit_distance\": %s,\n"
            "  \"min_hamming_distance\": %s,\n"
            "  \"safe_at_k0\": %s,\n"
            "  \"safe_at_k1\": %s,\n"
            "  \"safe_at_k2\": %s,\n"
            "  \"safe_at_hamming_k2\": %s,\n"
            "  \"safe_at_hamming_k3\": %s,\n"
            "  \"pairs_distance_0\": %llu,\n"
            "  \"pairs_distance_1\": %llu,\n"
            "  \"pairs_distance_2\": %llu,\n"
            "  \"pairs_within_requested_k\": %llu,\n"
            "  \"risk_pairs_for_k1\": %llu,\n"
            "  \"risk_pairs_for_k2\": %s,\n"
            "  \"risk_pairs_for_hamming_k2\": %s,\n"
            "  \"risk_pairs_for_hamming_k3\": %s,\n"
            "  \"ambiguous_query_variants_k1\": %llu,\n"
            "  \"recommended_k\": %d\n"
            "}\n",
            audit_mode, k, n_targets, unique_sequences, n_targets - unique_sequences,
            min_edit_distance_json, min_hamming_distance_json,
            safe_at_k0 ? "true" : "false", safe_at_k1 ? "true" : "false", safe_at_k2_json,
            safe_at_hamming_k2_json, safe_at_hamming_k3_json,
            pairs_d0, pairs_d1, pairs_d2, pairs_within_k, risk_pairs_k1, risk_pairs_k2_json,
            risk_pairs_hamming_k2_json, risk_pairs_hamming_k3_json,
            ambiguous_query_variants_k1, recommended_k);
    if (fclose(out) != 0) return -1;
    return 0;
}

static int audit_fast_outputs(const seq_table *targets, const char *out_dir, int k) {
    int rc = -1;
    int min_dist = -1;
    unsigned long long pairs_d0 = 0;
    unsigned long long pairs_d1 = 0;
    unsigned long long pairs_d2 = 0;
    unsigned long long pairs_within_k = 0;
    unsigned long long risk_pairs_k1 = 0;
    unsigned long long ambiguous_query_variants_k1 = 0;
    int *nearest_dist = NULL;
    size_t *nearest_idx = NULL;
    unsigned long long *near_k1 = NULL;
    size_t *parent = NULL;
    variant_record_list variants = {0};
    pair_record_list candidate_pairs = {0};
    pair_record_list unique_pairs = {0};
    FILE *pairs = NULL;
    FILE *clusters = NULL;
    FILE *safety = NULL;
    FILE *summary = NULL;
    FILE *variants_out = NULL;
    char path[4096];

    nearest_dist = (int *)malloc((targets->count == 0 ? 1 : targets->count) * sizeof(int));
    nearest_idx = (size_t *)malloc((targets->count == 0 ? 1 : targets->count) * sizeof(size_t));
    near_k1 = (unsigned long long *)calloc(targets->count == 0 ? 1 : targets->count, sizeof(unsigned long long));
    parent = (size_t *)malloc((targets->count == 0 ? 1 : targets->count) * sizeof(size_t));
    if (nearest_dist == NULL || nearest_idx == NULL || near_k1 == NULL || parent == NULL) goto done;
    for (size_t i = 0; i < targets->count; ++i) {
        nearest_dist[i] = -1;
        nearest_idx[i] = (size_t)-1;
        parent[i] = i;
    }

    for (size_t i = 0; i < targets->count; ++i) {
        string_list local = {0};
        if (add_k1_variants_for_target(&local, targets->records[i].seq, targets->records[i].len) != 0) {
            free_string_list(&local);
            goto done;
        }
        for (size_t v = 0; v < local.count; ++v) {
            if (push_variant_record(&variants, local.items[v], i) != 0) {
                free_string_list(&local);
                goto done;
            }
        }
        free_string_list(&local);
    }
    qsort(variants.items, variants.count, sizeof(variant_record), cmp_variant_record);

    if (path_join(path, sizeof(path), out_dir, "ambiguous_variants.tsv") != 0) goto done;
    variants_out = open_output_file(path);
    if (variants_out == NULL) goto done;
    fprintf(variants_out, "query_variant\ttargets_within_k1\n");

    for (size_t start = 0; start < variants.count;) {
        size_t end = start + 1;
        while (end < variants.count && strcmp(variants.items[start].key, variants.items[end].key) == 0) ++end;
        size_t unique_targets = 0;
        size_t last_target = (size_t)-1;
        for (size_t i = start; i < end; ++i) {
            if (variants.items[i].target != last_target) {
                variants.items[start + unique_targets].target = variants.items[i].target;
                last_target = variants.items[i].target;
                ++unique_targets;
            }
        }
        if (unique_targets >= 2) {
            ++ambiguous_query_variants_k1;
            fprintf(variants_out, "%s\t%zu\n", variants.items[start].key, unique_targets);
            for (size_t i = 0; i < unique_targets; ++i) {
                for (size_t j = i + 1; j < unique_targets; ++j) {
                    if (push_pair_record(&candidate_pairs, variants.items[start + i].target,
                                         variants.items[start + j].target) != 0) {
                        goto done;
                    }
                }
            }
        }
        start = end;
    }
    fclose(variants_out);
    variants_out = NULL;

    qsort(candidate_pairs.items, candidate_pairs.count, sizeof(pair_record), cmp_pair_record);
    for (size_t i = 0; i < candidate_pairs.count; ++i) {
        if (i > 0 && candidate_pairs.items[i].a == candidate_pairs.items[i - 1].a &&
            candidate_pairs.items[i].b == candidate_pairs.items[i - 1].b) {
            continue;
        }
        if (push_pair_record(&unique_pairs, candidate_pairs.items[i].a, candidate_pairs.items[i].b) != 0) goto done;
    }

    if (path_join(path, sizeof(path), out_dir, "collision_pairs.tsv") != 0) goto done;
    pairs = open_output_file(path);
    if (pairs == NULL) goto done;
    fprintf(pairs, "target_a\ttarget_b\tsequence_a\tsequence_b\tdistance\trisk_at_k1\trisk_at_k2\texample_ambiguous_query\n");

    for (size_t p = 0; p < unique_pairs.count; ++p) {
        size_t i = unique_pairs.items[p].a;
        size_t j = unique_pairs.items[p].b;
        int d = qdaln_edit_distance(targets->records[i].seq, targets->records[i].len,
                                    targets->records[j].seq, targets->records[j].len);
        if (d < 0) goto done;
        if (min_dist < 0 || d < min_dist) min_dist = d;
        if (nearest_dist[i] < 0 || d < nearest_dist[i]) {
            nearest_dist[i] = d;
            nearest_idx[i] = j;
        }
        if (nearest_dist[j] < 0 || d < nearest_dist[j]) {
            nearest_dist[j] = d;
            nearest_idx[j] = i;
        }
        if (d == 0) ++pairs_d0;
        if (d == 1) ++pairs_d1;
        if (d == 2) ++pairs_d2;
        if (d <= k) ++pairs_within_k;
        if (d <= 2) {
            ++risk_pairs_k1;
            ++near_k1[i];
            ++near_k1[j];
            uf_union(parent, i, j);
        }
        fprintf(pairs, "%s\t%s\t%s\t%s\t%d\t%s\tnot_computed\t\n",
                targets->records[i].id, targets->records[j].id, targets->records[i].seq, targets->records[j].seq,
                d, d <= 2 ? "yes" : "no");
    }
    fclose(pairs);
    pairs = NULL;

    if (path_join(path, sizeof(path), out_dir, "target_safety.tsv") != 0) goto done;
    safety = open_output_file(path);
    if (safety == NULL) goto done;
    fprintf(safety, "target_id\tsequence\tnearest_target\tnearest_distance\tsafe_at_k1\tsafe_at_k2\tnum_nearby_k1_risk_targets\n");
    for (size_t i = 0; i < targets->count; ++i) {
        const char *near_id = nearest_idx[i] == (size_t)-1 ? "" : targets->records[nearest_idx[i]].id;
        int nd = nearest_dist[i];
        fprintf(safety, "%s\t%s\t%s\t%d\t%s\tnot_computed\t%llu\n",
                targets->records[i].id, targets->records[i].seq, near_id, nd,
                (nd < 0 || nd >= 3) ? "yes" : "no", near_k1[i]);
    }
    fclose(safety);
    safety = NULL;

    if (path_join(path, sizeof(path), out_dir, "collision_clusters.tsv") != 0) goto done;
    clusters = open_output_file(path);
    if (clusters == NULL) goto done;
    fprintf(clusters, "cluster_id\ttarget_id\tsequence\n");
    for (size_t i = 0; i < targets->count; ++i) {
        if (near_k1[i] == 0) continue;
        fprintf(clusters, "%zu\t%s\t%s\n", uf_find(parent, i), targets->records[i].id, targets->records[i].seq);
    }
    fclose(clusters);
    clusters = NULL;

    size_t unique_sequences = count_unique_target_sequences(targets);
    if (path_join(path, sizeof(path), out_dir, "audit_summary.tsv") != 0) goto done;
    summary = open_output_file(path);
    if (summary == NULL) goto done;
    fprintf(summary, "metric\tvalue\n");
    fprintf(summary, "audit_mode\tfast\n");
    fprintf(summary, "targets\t%zu\n", targets->count);
    fprintf(summary, "unique_sequences\t%zu\n", unique_sequences);
    fprintf(summary, "duplicate_sequences\t%zu\n", targets->count - unique_sequences);
    fprintf(summary, "min_edit_distance\t%s\n", min_dist < 0 ? ">=3" : (min_dist == 0 ? "0" : (min_dist == 1 ? "1" : "2")));
    fprintf(summary, "safe_at_k0\t%s\n", pairs_d0 == 0 ? "yes" : "no");
    fprintf(summary, "safe_at_k1\t%s\n", risk_pairs_k1 == 0 ? "yes" : "no");
    fprintf(summary, "safe_at_k2\tnot_computed\n");
    fprintf(summary, "safe_at_hamming_k2\tnot_computed\n");
    fprintf(summary, "safe_at_hamming_k3\tnot_computed\n");
    fprintf(summary, "pairs_distance_0\t%llu\n", pairs_d0);
    fprintf(summary, "pairs_distance_1\t%llu\n", pairs_d1);
    fprintf(summary, "pairs_distance_2\t%llu\n", pairs_d2);
    fprintf(summary, "pairs_within_requested_k\t%llu\n", pairs_within_k);
    fprintf(summary, "risk_pairs_for_k1\t%llu\n", risk_pairs_k1);
    fprintf(summary, "risk_pairs_for_k2\tnot_computed\n");
    fprintf(summary, "risk_pairs_for_hamming_k2\tnot_computed\n");
    fprintf(summary, "risk_pairs_for_hamming_k3\tnot_computed\n");
    fprintf(summary, "ambiguous_query_variants_k1\t%llu\n", ambiguous_query_variants_k1);
    fprintf(summary, "recommended_k\t%d\n", pairs_d0 == 0 && (k == 0 || risk_pairs_k1 == 0) ? k : 0);
    fclose(summary);
    summary = NULL;
    if (write_audit_summary_json(out_dir, "fast", k, targets->count, unique_sequences,
                                 min_dist < 0 ? "\">=3\"" : (min_dist == 0 ? "0" : (min_dist == 1 ? "1" : "2")),
                                 "null",
                                 pairs_d0 == 0, risk_pairs_k1 == 0, "null",
                                 "null", "null",
                                 pairs_d0, pairs_d1, pairs_d2, pairs_within_k, risk_pairs_k1, "null",
                                 "null", "null",
                                 ambiguous_query_variants_k1,
                                 pairs_d0 == 0 && (k == 0 || risk_pairs_k1 == 0) ? k : 0) != 0) {
        goto done;
    }
    rc = 0;

done:
    if (pairs != NULL) fclose(pairs);
    if (clusters != NULL) fclose(clusters);
    if (safety != NULL) fclose(safety);
    if (summary != NULL) fclose(summary);
    if (variants_out != NULL) fclose(variants_out);
    free(nearest_dist);
    free(nearest_idx);
    free(near_k1);
    free(parent);
    free_variant_record_list(&variants);
    free_pair_record_list(&candidate_pairs);
    free_pair_record_list(&unique_pairs);
    return rc;
}

static int run_audit(const char *argv0, int argc, char **argv) {
    const char *targets_path = NULL;
    const char *out_dir = NULL;
    const char *audit_mode = "auto";
    int k = 1;

    int i = 2;
    while (i < argc) {
        const char *arg = argv[i++];
        if ((strcmp(arg, "--targets") == 0 || strcmp(arg, "--library") == 0) && i < argc) {
            targets_path = argv[i++];
        } else if (strcmp(arg, "--k") == 0 && i < argc) {
            if (parse_int_value(argv[i++], &k) != 0 || k < 0 || k > 3) {
                usage(argv0);
                return 2;
            }
        } else if ((strcmp(arg, "--out-dir") == 0 || strcmp(arg, "--out") == 0) && i < argc) {
            out_dir = argv[i++];
        } else if (strcmp(arg, "--audit-mode") == 0 && i < argc) {
            audit_mode = argv[i++];
            if (strcmp(audit_mode, "auto") != 0 && strcmp(audit_mode, "exact") != 0 && strcmp(audit_mode, "fast") != 0) {
                usage(argv0);
                return 2;
            }
        } else {
            usage(argv0);
            return 2;
        }
    }
    if (targets_path == NULL || out_dir == NULL) {
        usage(argv0);
        return 2;
    }

    seq_table targets = {0};
    int rc = 1;
    int min_dist = -1;
    unsigned long long pairs_d0 = 0;
    unsigned long long pairs_d1 = 0;
    unsigned long long pairs_d2 = 0;
    unsigned long long pairs_within_k = 0;
    unsigned long long risk_pairs_k1 = 0;
    unsigned long long risk_pairs_k2 = 0;
    unsigned long long risk_pairs_hamming_k2 = 0;
    unsigned long long risk_pairs_hamming_k3 = 0;
    int min_hamming_dist = -1;
    int *nearest_dist = NULL;
    size_t *nearest_idx = NULL;
    unsigned long long *near_k1 = NULL;
    size_t *parent = NULL;
    FILE *pairs = NULL;
    FILE *clusters = NULL;
    FILE *safety = NULL;
    FILE *summary = NULL;
    FILE *variants_out = NULL;
    string_list k1_variants = {0};
    unsigned long long ambiguous_query_variants_k1 = 0;
    char path[4096];

    if (read_target_table(targets_path, &targets) != 0) {
        fprintf(stderr, "failed to read targets\n");
        goto done;
    }
    if (ensure_dir(out_dir) != 0) {
        fprintf(stderr, "failed to create audit output directory\n");
        goto done;
    }
    int use_fast = strcmp(audit_mode, "fast") == 0 || (strcmp(audit_mode, "auto") == 0 && targets.count > 2000);
    if (use_fast) {
        rc = audit_fast_outputs(&targets, out_dir, k) == 0 ? 0 : 1;
        if (rc == 0) printf("%s\n", out_dir);
        goto done;
    }
    nearest_dist = (int *)malloc((targets.count == 0 ? 1 : targets.count) * sizeof(int));
    nearest_idx = (size_t *)malloc((targets.count == 0 ? 1 : targets.count) * sizeof(size_t));
    near_k1 = (unsigned long long *)calloc(targets.count == 0 ? 1 : targets.count, sizeof(unsigned long long));
    parent = (size_t *)malloc((targets.count == 0 ? 1 : targets.count) * sizeof(size_t));
    if (nearest_dist == NULL || nearest_idx == NULL || near_k1 == NULL || parent == NULL) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    for (size_t i = 0; i < targets.count; ++i) {
        nearest_dist[i] = -1;
        nearest_idx[i] = (size_t)-1;
        parent[i] = i;
    }

    if (path_join(path, sizeof(path), out_dir, "collision_pairs.tsv") != 0) goto done;
    pairs = open_output_file(path);
    if (pairs == NULL) goto done;
    fprintf(pairs, "target_a\ttarget_b\tsequence_a\tsequence_b\tdistance\trisk_at_k1\trisk_at_k2\texample_ambiguous_query\n");

    for (size_t i = 0; i < targets.count; ++i) {
        for (size_t j = i + 1; j < targets.count; ++j) {
            int d = qdaln_edit_distance(targets.records[i].seq, targets.records[i].len,
                                        targets.records[j].seq, targets.records[j].len);
            if (d < 0) goto done;
            if (min_dist < 0 || d < min_dist) min_dist = d;
            if (nearest_dist[i] < 0 || d < nearest_dist[i]) {
                nearest_dist[i] = d;
                nearest_idx[i] = j;
            }
            if (nearest_dist[j] < 0 || d < nearest_dist[j]) {
                nearest_dist[j] = d;
                nearest_idx[j] = i;
            }
            if (d == 0) ++pairs_d0;
            if (d == 1) ++pairs_d1;
            if (d == 2) ++pairs_d2;
            if (d <= k) ++pairs_within_k;
            if (d <= 2) {
                ++risk_pairs_k1;
                ++near_k1[i];
                ++near_k1[j];
                uf_union(parent, i, j);
            }
            if (d <= 4) ++risk_pairs_k2;
            int hd = hamming_distance_cli(targets.records[i].seq, targets.records[i].len,
                                          targets.records[j].seq, targets.records[j].len);
            if (hd >= 0) {
                if (min_hamming_dist < 0 || hd < min_hamming_dist) min_hamming_dist = hd;
                if (hd <= 4) ++risk_pairs_hamming_k2;
                if (hd <= 6) ++risk_pairs_hamming_k3;
            }
            if (d <= 2 || d <= 2 * k) {
                const char *example = d == 0 ? targets.records[i].seq : "";
                fprintf(pairs, "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n",
                        targets.records[i].id, targets.records[j].id, targets.records[i].seq, targets.records[j].seq,
                        d, d <= 2 ? "yes" : "no", d <= 4 ? "yes" : "no", example);
            }
        }
    }
    fclose(pairs);
    pairs = NULL;

    if (path_join(path, sizeof(path), out_dir, "target_safety.tsv") != 0) goto done;
    safety = open_output_file(path);
    if (safety == NULL) goto done;
    fprintf(safety, "target_id\tsequence\tnearest_target\tnearest_distance\tsafe_at_k1\tsafe_at_k2\tnum_nearby_k1_risk_targets\n");
    for (size_t i = 0; i < targets.count; ++i) {
        const char *near_id = nearest_idx[i] == (size_t)-1 ? "" : targets.records[nearest_idx[i]].id;
        int nd = nearest_dist[i];
        fprintf(safety, "%s\t%s\t%s\t%d\t%s\t%s\t%llu\n",
                targets.records[i].id, targets.records[i].seq, near_id, nd,
                (nd < 0 || nd >= 3) ? "yes" : "no",
                (nd < 0 || nd >= 5) ? "yes" : "no",
                near_k1[i]);
    }
    fclose(safety);
    safety = NULL;

    if (path_join(path, sizeof(path), out_dir, "collision_clusters.tsv") != 0) goto done;
    clusters = open_output_file(path);
    if (clusters == NULL) goto done;
    fprintf(clusters, "cluster_id\ttarget_id\tsequence\n");
    for (size_t i = 0; i < targets.count; ++i) {
        if (near_k1[i] == 0) continue;
        fprintf(clusters, "%zu\t%s\t%s\n", uf_find(parent, i), targets.records[i].id, targets.records[i].seq);
    }
    fclose(clusters);
    clusters = NULL;

    size_t unique_sequences = 0;
    for (size_t i = 0; i < targets.count; ++i) {
        int seen = 0;
        for (size_t j = 0; j < i; ++j) {
            if (targets.records[i].len == targets.records[j].len &&
                memcmp(targets.records[i].seq, targets.records[j].seq, targets.records[i].len) == 0) {
                seen = 1;
                break;
            }
        }
        if (!seen) ++unique_sequences;
    }

    if (path_join(path, sizeof(path), out_dir, "ambiguous_variants.tsv") != 0) goto done;
    variants_out = open_output_file(path);
    if (variants_out == NULL) goto done;
    fprintf(variants_out, "query_variant\ttargets_within_k1\n");
    for (size_t i = 0; i < targets.count; ++i) {
        if (add_k1_variants_for_target(&k1_variants, targets.records[i].seq, targets.records[i].len) != 0) {
            fprintf(stderr, "failed to enumerate k=1 variants\n");
            goto done;
        }
    }
    for (size_t vi = 0; vi < k1_variants.count; ++vi) {
        unsigned long long within = 0;
        size_t q_len = strlen(k1_variants.items[vi]);
        for (size_t ti = 0; ti < targets.count; ++ti) {
            int ok = qdaln_edit_distance_leq(k1_variants.items[vi], q_len, targets.records[ti].seq,
                                             targets.records[ti].len, 1);
            if (ok < 0) goto done;
            if (ok) ++within;
        }
        if (within >= 2) {
            ++ambiguous_query_variants_k1;
            fprintf(variants_out, "%s\t%llu\n", k1_variants.items[vi], within);
        }
    }
    fclose(variants_out);
    variants_out = NULL;

    if (path_join(path, sizeof(path), out_dir, "audit_summary.tsv") != 0) goto done;
    summary = open_output_file(path);
    if (summary == NULL) goto done;
    fprintf(summary, "metric\tvalue\n");
    fprintf(summary, "audit_mode\texact\n");
    fprintf(summary, "targets\t%zu\n", targets.count);
    fprintf(summary, "unique_sequences\t%zu\n", unique_sequences);
    fprintf(summary, "duplicate_sequences\t%zu\n", targets.count - unique_sequences);
    fprintf(summary, "min_edit_distance\t%d\n", min_dist);
    if (min_hamming_dist >= 0) {
        fprintf(summary, "min_hamming_distance\t%d\n", min_hamming_dist);
    } else {
        fprintf(summary, "min_hamming_distance\tnot_computed\n");
    }
    fprintf(summary, "safe_at_k0\t%s\n", pairs_d0 == 0 ? "yes" : "no");
    fprintf(summary, "safe_at_k1\t%s\n", risk_pairs_k1 == 0 ? "yes" : "no");
    fprintf(summary, "safe_at_k2\t%s\n", risk_pairs_k2 == 0 ? "yes" : "no");
    fprintf(summary, "safe_at_hamming_k2\t%s\n", risk_pairs_hamming_k2 == 0 ? "yes" : "no");
    fprintf(summary, "safe_at_hamming_k3\t%s\n", risk_pairs_hamming_k3 == 0 ? "yes" : "no");
    fprintf(summary, "pairs_distance_0\t%llu\n", pairs_d0);
    fprintf(summary, "pairs_distance_1\t%llu\n", pairs_d1);
    fprintf(summary, "pairs_distance_2\t%llu\n", pairs_d2);
    fprintf(summary, "pairs_within_requested_k\t%llu\n", pairs_within_k);
    fprintf(summary, "risk_pairs_for_k1\t%llu\n", risk_pairs_k1);
    fprintf(summary, "risk_pairs_for_k2\t%llu\n", risk_pairs_k2);
    fprintf(summary, "risk_pairs_for_hamming_k2\t%llu\n", risk_pairs_hamming_k2);
    fprintf(summary, "risk_pairs_for_hamming_k3\t%llu\n", risk_pairs_hamming_k3);
    fprintf(summary, "ambiguous_query_variants_k1\t%llu\n", ambiguous_query_variants_k1);
    fprintf(summary, "recommended_k\t%d\n", pairs_d0 == 0 && (k == 0 || risk_pairs_k1 == 0) ? k : 0);
    fclose(summary);
    summary = NULL;
    char min_dist_json[32];
    char min_hamming_dist_json[32];
    char risk_pairs_k2_json[32];
    char risk_pairs_hamming_k2_json[32];
    char risk_pairs_hamming_k3_json[32];
    snprintf(min_dist_json, sizeof(min_dist_json), "%d", min_dist);
    snprintf(min_hamming_dist_json, sizeof(min_hamming_dist_json), "%d", min_hamming_dist);
    snprintf(risk_pairs_k2_json, sizeof(risk_pairs_k2_json), "%llu", risk_pairs_k2);
    snprintf(risk_pairs_hamming_k2_json, sizeof(risk_pairs_hamming_k2_json), "%llu", risk_pairs_hamming_k2);
    snprintf(risk_pairs_hamming_k3_json, sizeof(risk_pairs_hamming_k3_json), "%llu", risk_pairs_hamming_k3);
    if (write_audit_summary_json(out_dir, "exact", k, targets.count, unique_sequences,
                                 min_dist_json, min_hamming_dist >= 0 ? min_hamming_dist_json : "null",
                                 pairs_d0 == 0, risk_pairs_k1 == 0,
                                 risk_pairs_k2 == 0 ? "true" : "false",
                                 risk_pairs_hamming_k2 == 0 ? "true" : "false",
                                 risk_pairs_hamming_k3 == 0 ? "true" : "false",
                                 pairs_d0, pairs_d1, pairs_d2, pairs_within_k, risk_pairs_k1,
                                 risk_pairs_k2_json, risk_pairs_hamming_k2_json, risk_pairs_hamming_k3_json,
                                 ambiguous_query_variants_k1,
                                 pairs_d0 == 0 && (k == 0 || risk_pairs_k1 == 0) ? k : 0) != 0) {
        goto done;
    }

    printf("%s\n", out_dir);
    rc = 0;

done:
    if (pairs != NULL) fclose(pairs);
    if (clusters != NULL) fclose(clusters);
    if (safety != NULL) fclose(safety);
    if (summary != NULL) fclose(summary);
    if (variants_out != NULL) fclose(variants_out);
    free_string_list(&k1_variants);
    free(nearest_dist);
    free(nearest_idx);
    free(near_k1);
    free(parent);
    free_table(&targets);
    return rc;
}

static void write_fastq_record(FILE *out, const char *header, const char *seq, const char *plus, const char *qual) {
    fprintf(out, "%s\n%s\n%s\n%s\n", header, seq, plus, qual);
}

typedef struct unmatched_entry {
    char *seq;
    unsigned long long count;
    int offset_hint;
    unsigned long long low_quality_count;
    char *adapter_hint;
} unmatched_entry;

typedef struct unmatched_table {
    unmatched_entry *entries;
    size_t count;
    size_t cap;
} unmatched_table;

static void free_unmatched_table(unmatched_table *table) {
    for (size_t i = 0; i < table->count; ++i) {
        free(table->entries[i].seq);
        free(table->entries[i].adapter_hint);
    }
    free(table->entries);
    table->entries = NULL;
    table->count = 0;
    table->cap = 0;
}

static int add_unmatched_observation(unmatched_table *table, const char *seq, int offset_hint,
                                     int low_quality, const char *adapter_hint) {
    for (size_t i = 0; i < table->count; ++i) {
        if (strcmp(table->entries[i].seq, seq) == 0) {
            ++table->entries[i].count;
            if (table->entries[i].offset_hint == 0 && offset_hint != 0) table->entries[i].offset_hint = offset_hint;
            if (low_quality) ++table->entries[i].low_quality_count;
            if ((table->entries[i].adapter_hint == NULL || table->entries[i].adapter_hint[0] == '\0') &&
                adapter_hint != NULL && adapter_hint[0] != '\0') {
                free(table->entries[i].adapter_hint);
                table->entries[i].adapter_hint = xstrndup(adapter_hint, strlen(adapter_hint));
                if (table->entries[i].adapter_hint == NULL) return -1;
            }
            return 0;
        }
    }
    if (table->count == table->cap) {
        size_t next_cap = table->cap == 0 ? 16 : table->cap * 2;
        unmatched_entry *next = (unmatched_entry *)realloc(table->entries, next_cap * sizeof(unmatched_entry));
        if (next == NULL) return -1;
        table->entries = next;
        table->cap = next_cap;
    }
    table->entries[table->count].seq = xstrndup(seq, strlen(seq));
    if (table->entries[table->count].seq == NULL) return -1;
    table->entries[table->count].count = 1;
    table->entries[table->count].offset_hint = offset_hint;
    table->entries[table->count].low_quality_count = low_quality ? 1 : 0;
    table->entries[table->count].adapter_hint = adapter_hint == NULL ? xstrndup("", 0) : xstrndup(adapter_hint, strlen(adapter_hint));
    if (table->entries[table->count].adapter_hint == NULL) {
        free(table->entries[table->count].seq);
        return -1;
    }
    ++table->count;
    return 0;
}

static int cmp_unmatched_entry_desc(const void *a, const void *b) {
    const unmatched_entry *aa = (const unmatched_entry *)a;
    const unmatched_entry *bb = (const unmatched_entry *)b;
    if (aa->count != bb->count) return aa->count < bb->count ? 1 : -1;
    return strcmp(aa->seq, bb->seq);
}

static int contains_base_n(const char *seq) {
    for (; *seq != '\0'; ++seq) {
        if (*seq == 'N' || *seq == 'n') return 1;
    }
    return 0;
}

static char complement_base(char c) {
    switch (c) {
        case 'A':
        case 'a':
            return 'T';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        case 'T':
        case 't':
            return 'A';
        default:
            return 'N';
    }
}

static int reverse_complement_seq(const char *seq, char *out, size_t out_cap) {
    size_t len = strlen(seq);
    if (len + 1 > out_cap) return -1;
    for (size_t i = 0; i < len; ++i) out[i] = complement_base(seq[len - 1 - i]);
    out[len] = '\0';
    return 0;
}

static int nearest_target_for_query(const seq_table *targets, const char *query, int *nearest_index, int *nearest_dist) {
    *nearest_index = -1;
    *nearest_dist = -1;
    size_t q_len = strlen(query);
    for (size_t i = 0; i < targets->count; ++i) {
        int d = qdaln_edit_distance(query, q_len, targets->records[i].seq, targets->records[i].len);
        if (d < 0) return -1;
        if (*nearest_dist < 0 || d < *nearest_dist) {
            *nearest_dist = d;
            *nearest_index = (int)i;
        }
    }
    return 0;
}

static int find_offset_hint(const qdaln_index *index, const char *seq, size_t seq_len, size_t target_start,
                            size_t target_len, int k, size_t offset_window) {
    if (offset_window == 0) return 0;
    char observed[8192];
    if (target_len >= sizeof(observed)) return 0;
    for (size_t step = 1; step <= offset_window; ++step) {
        for (int sign = 1; sign >= -1; sign -= 2) {
            if (sign < 0 && target_start < step) continue;
            size_t offset = sign > 0 ? target_start + step : target_start - step;
            if (offset > seq_len || target_len > seq_len - offset) continue;
            qdaln_match_result r;
            qdaln_index_stats stats;
            if (k == 0) {
                if (qdaln_index_lookup_exact_ascii_stats(index, seq + offset, target_len, &r, &stats) != 0) {
                    return 0;
                }
            } else {
                memcpy(observed, seq + offset, target_len);
                observed[target_len] = '\0';
                uppercase_ascii(observed);
                const char *read_ptr = observed;
                size_t read_len = target_len;
                if (qdaln_index_assign_stats(index, &read_ptr, &read_len, 1, k, &r, &stats) != 0) return 0;
            }
            if (r.status == QDALN_MATCH_UNIQUE) return sign > 0 ? (int)step : -(int)step;
        }
    }
    return 0;
}

static int window_has_low_quality(const char *qual, size_t target_start, size_t target_len, int threshold) {
    if (threshold < 0) return 0;
    size_t qual_len = strlen(qual);
    if (target_start > qual_len || target_len > qual_len - target_start) return 0;
    for (size_t i = 0; i < target_len; ++i) {
        int phred = (int)((unsigned char)qual[target_start + i]) - 33;
        if (phred < threshold) return 1;
    }
    return 0;
}

static int run_inspect_unmatched(const char *argv0, int argc, char **argv) {
    const char *targets_path = NULL;
    const char *reads_path = NULL;
    const char *out_path = NULL;
    size_t target_start = 0;
    size_t target_len = 0;
    size_t top_n = 100;
    size_t offset_window = 0;
    char adapter[1024] = "";
    int low_quality_threshold = -1;
    int k = -1;

    int i = 2;
    while (i < argc) {
        const char *arg = argv[i++];
        if ((strcmp(arg, "--targets") == 0 || strcmp(arg, "--library") == 0) && i < argc) {
            targets_path = argv[i++];
        } else if (strcmp(arg, "--reads") == 0 && i < argc) {
            reads_path = argv[i++];
        } else if ((strcmp(arg, "--target-start") == 0 || strcmp(arg, "--guide-start") == 0) && i < argc) {
            if (parse_size_value(argv[i++], &target_start) != 0) {
                usage(argv0);
                return 2;
            }
        } else if ((strcmp(arg, "--target-length") == 0 || strcmp(arg, "--guide-length") == 0) && i < argc) {
            if (parse_size_value(argv[i++], &target_len) != 0 || target_len == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--k") == 0 && i < argc) {
            if (parse_int_value(argv[i++], &k) != 0 || (k != 0 && k != 1)) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--top") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &top_n) != 0 || top_n == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--offset-window") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &offset_window) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--adapter") == 0 && i < argc) {
            const char *value = argv[i++];
            size_t n = strlen(value);
            if (n == 0 || n >= sizeof(adapter)) {
                usage(argv0);
                return 2;
            }
            memcpy(adapter, value, n + 1);
            uppercase_ascii(adapter);
        } else if (strcmp(arg, "--low-quality-threshold") == 0 && i < argc) {
            if (parse_int_value(argv[i++], &low_quality_threshold) != 0 ||
                low_quality_threshold < 0 || low_quality_threshold > 93) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--out") == 0 && i < argc) {
            out_path = argv[i++];
        } else {
            usage(argv0);
            return 2;
        }
    }
    if (targets_path == NULL || reads_path == NULL || out_path == NULL || target_len == 0 || k < 0) {
        usage(argv0);
        return 2;
    }

    seq_table targets = {0};
    const char **target_ptrs = NULL;
    size_t *target_lens = NULL;
    qdaln_index *index = NULL;
    fastq_reader reader = {0};
    unmatched_table unmatched = {0};
    FILE *out = NULL;
    int rc = 1;

    if (read_target_table(targets_path, &targets) != 0) {
        fprintf(stderr, "failed to read targets\n");
        goto done;
    }
    if (build_target_arrays(&targets, &target_ptrs, &target_lens) != 0) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    index = qdaln_index_build(target_ptrs, target_lens, targets.count);
    if (index == NULL) {
        fprintf(stderr, "failed to build target index\n");
        goto done;
    }
    if (fastq_reader_open(&reader, reads_path) != 0) {
        fprintf(stderr, "failed to open FASTQ input\n");
        goto done;
    }

    char header[8192];
    char seq[8192];
    char plus[8192];
    char qual[8192];
    char observed[8192];
    int got = 0;
    size_t seq_len = 0;
    while ((got = fastq_read_record_len(&reader, header, seq, plus, qual, sizeof(header), &seq_len)) == 1) {
        qdaln_match_result r = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
        observed[0] = '\0';
        if (target_start <= seq_len && target_len <= seq_len - target_start && target_len < sizeof(observed)) {
            memcpy(observed, seq + target_start, target_len);
            observed[target_len] = '\0';
            uppercase_ascii(observed);
            qdaln_index_stats stats;
            int assign_rc = 0;
            if (k == 0) {
                assign_rc = qdaln_index_lookup_exact_ascii_stats(index, seq + target_start, target_len, &r, &stats);
            } else {
                const char *read_ptr = observed;
                size_t read_len = target_len;
                assign_rc = qdaln_index_assign_stats(index, &read_ptr, &read_len, 1, k, &r, &stats);
            }
            if (assign_rc != 0) {
                fprintf(stderr, "assignment failed\n");
                goto done;
            }
        } else {
            strncpy(observed, "<invalid>", sizeof(observed) - 1);
            observed[sizeof(observed) - 1] = '\0';
        }
        if (r.status == QDALN_MATCH_NONE || r.status == QDALN_MATCH_INVALID) {
            int offset_hint = strcmp(observed, "<invalid>") == 0 ? 0 :
                    find_offset_hint(index, seq, seq_len, target_start, target_len, k, offset_window);
            char seq_upper[8192];
            seq_upper[0] = '\0';
            const char *adapter_hint = "";
            int low_quality = strcmp(observed, "<invalid>") == 0 ? 0 :
                    window_has_low_quality(qual, target_start, target_len, low_quality_threshold);
            if (adapter[0] != '\0') {
                strncpy(seq_upper, seq, sizeof(seq_upper) - 1);
                seq_upper[sizeof(seq_upper) - 1] = '\0';
                uppercase_ascii(seq_upper);
                if (strstr(seq_upper, adapter) != NULL || strstr(observed, adapter) != NULL) adapter_hint = adapter;
            }
            if (add_unmatched_observation(&unmatched, observed, offset_hint, low_quality, adapter_hint) != 0) {
                fprintf(stderr, "out of memory\n");
                goto done;
            }
        }
    }
    if (got < 0) {
        fprintf(stderr, "malformed FASTQ input\n");
        goto done;
    }

    qsort(unmatched.entries, unmatched.count, sizeof(unmatched_entry), cmp_unmatched_entry_desc);
    out = open_output_file(out_path);
    if (out == NULL) {
        fprintf(stderr, "failed to open unmatched inspection output\n");
        goto done;
    }
    fprintf(out, "sequence\tcount\tlength\tnearest_target\tnearest_distance\tnearest_edit_class\tpossible_reason\treverse_complement\trevcomp_nearest_target\trevcomp_nearest_distance\toffset_hint\tadapter_hint\n");
    size_t limit = unmatched.count < top_n ? unmatched.count : top_n;
    for (size_t i = 0; i < limit; ++i) {
        int nearest_index = -1;
        int nearest_dist = -1;
        int rc_nearest_index = -1;
        int rc_nearest_dist = -1;
        const char *nearest_id = "";
        const char *rc_nearest_id = "";
        const char *edit_class = "invalid";
        const char *reason = "wrong_length";
        char rc_seq[8192] = "";
        if (strcmp(unmatched.entries[i].seq, "<invalid>") != 0 &&
            nearest_target_for_query(&targets, unmatched.entries[i].seq, &nearest_index, &nearest_dist) == 0) {
            if (nearest_index >= 0) {
                nearest_id = targets.records[nearest_index].id;
                int kind = correction_kind(unmatched.entries[i].seq, strlen(unmatched.entries[i].seq),
                                           targets.records[nearest_index].seq, targets.records[nearest_index].len,
                                           nearest_dist);
                edit_class = correction_name(kind);
            }
            if (contains_base_n(unmatched.entries[i].seq)) reason = "contains_N";
            else if (nearest_dist > k) reason = "near_known_target_above_k";
            else reason = "unknown";
            if (unmatched.entries[i].low_quality_count != 0) reason = "low_quality_candidate";
            if (unmatched.entries[i].adapter_hint != NULL && unmatched.entries[i].adapter_hint[0] != '\0') {
                reason = "adapter_or_primer_candidate";
            }
            if (unmatched.entries[i].offset_hint != 0) reason = "offset_shift_candidate";
            if (reverse_complement_seq(unmatched.entries[i].seq, rc_seq, sizeof(rc_seq)) == 0 &&
                nearest_target_for_query(&targets, rc_seq, &rc_nearest_index, &rc_nearest_dist) == 0 &&
                rc_nearest_index >= 0) {
                rc_nearest_id = targets.records[rc_nearest_index].id;
                if (unmatched.entries[i].offset_hint == 0 &&
                    rc_nearest_dist <= k && (nearest_dist < 0 || rc_nearest_dist < nearest_dist)) {
                    reason = "reverse_complement_candidate";
                }
            }
        }
        fprintf(out, "%s\t%llu\t%zu\t%s\t%d\t%s\t%s\t%s\t%s\t%d\t",
                unmatched.entries[i].seq, unmatched.entries[i].count, strlen(unmatched.entries[i].seq),
                nearest_id, nearest_dist, edit_class, reason, rc_seq, rc_nearest_id, rc_nearest_dist);
        if (unmatched.entries[i].offset_hint != 0) fprintf(out, "%d", unmatched.entries[i].offset_hint);
        fprintf(out, "\t%s\n", unmatched.entries[i].adapter_hint == NULL ? "" : unmatched.entries[i].adapter_hint);
    }
    rc = 0;

done:
    if (out != NULL) fclose(out);
    fastq_reader_close(&reader);
    qdaln_index_free(index);
    free(target_ptrs);
    free(target_lens);
    free_unmatched_table(&unmatched);
    free_table(&targets);
    return rc;
}

typedef struct pair_count_stats {
    unsigned long long total_reads;
    unsigned long long assigned_pairs;
    unsigned long long pair_ambiguous;
    unsigned long long left_unmatched;
    unsigned long long right_unmatched;
    unsigned long long invalid;
    unsigned long long left_invalid;
    unsigned long long right_invalid;
    unsigned long long candidates_considered;
    unsigned long long candidates_verified;
} pair_count_stats;

static const char *pair_status_name(qdaln_match_result left, qdaln_match_result right) {
    if (left.status == QDALN_MATCH_INVALID || right.status == QDALN_MATCH_INVALID) return "invalid";
    if (left.status == QDALN_MATCH_AMBIGUOUS || right.status == QDALN_MATCH_AMBIGUOUS) return "ambiguous";
    if (left.status == QDALN_MATCH_NONE || right.status == QDALN_MATCH_NONE) return "none";
    if (left.status == QDALN_MATCH_UNIQUE && right.status == QDALN_MATCH_UNIQUE) return "unique";
    return "invalid";
}

static void print_pair_assignment_row(FILE *out, const char *read_id, const seq_table *left_targets,
                                      const seq_table *right_targets, const char *left_observed,
                                      qdaln_match_result left, const char *right_observed,
                                      qdaln_match_result right) {
    const char *left_id = left.target_index >= 0 ? left_targets->records[left.target_index].id : "";
    const char *right_id = right.target_index >= 0 ? right_targets->records[right.target_index].id : "";
    fprintf(out, "%s\t%s\t%d\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%d\t%s\n",
            read_id, left_observed, left.target_index, left_id, status_name(left.status), left.best_distance,
            right_observed, right.target_index, right_id, status_name(right.status), right.best_distance,
            pair_status_name(left, right));
}

static int run_pair_count(const char *argv0, int argc, char **argv) {
    const char *left_path = NULL;
    const char *right_path = NULL;
    const char *reads_path = NULL;
    const char *left_reads_path = NULL;
    const char *right_reads_path = NULL;
    const char *out_path = NULL;
    const char *summary_path = NULL;
    const char *assignments_path = NULL;
    size_t left_start = 0;
    size_t right_start = 0;
    size_t left_len = 0;
    size_t right_len = 0;
    int k = -1;
    count_metric metric = COUNT_METRIC_LEVENSHTEIN;
    ambiguity_policy assignment_policy = AMBIGUITY_POLICY_RADIUS;

    int i = 2;
    while (i < argc) {
        const char *arg = argv[i++];
        if (strcmp(arg, "--left-targets") == 0 && i < argc) {
            left_path = argv[i++];
        } else if (strcmp(arg, "--right-targets") == 0 && i < argc) {
            right_path = argv[i++];
        } else if (strcmp(arg, "--reads") == 0 && i < argc) {
            reads_path = argv[i++];
        } else if (strcmp(arg, "--left-reads") == 0 && i < argc) {
            left_reads_path = argv[i++];
        } else if (strcmp(arg, "--right-reads") == 0 && i < argc) {
            right_reads_path = argv[i++];
        } else if (strcmp(arg, "--left-start") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &left_start) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--left-length") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &left_len) != 0 || left_len == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--right-start") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &right_start) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--right-length") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &right_len) != 0 || right_len == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--k") == 0 && i < argc) {
            if (parse_int_value(argv[i++], &k) != 0 || k < 0 || k > 2) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--metric") == 0 && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "hamming") == 0) {
                metric = COUNT_METRIC_HAMMING;
            } else if (strcmp(value, "levenshtein") == 0) {
                metric = COUNT_METRIC_LEVENSHTEIN;
            } else {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--ambiguity-policy") == 0 && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "radius") == 0) {
                assignment_policy = AMBIGUITY_POLICY_RADIUS;
            } else if (strcmp(value, "best") == 0) {
                assignment_policy = AMBIGUITY_POLICY_BEST;
            } else {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--out") == 0 && i < argc) {
            out_path = argv[i++];
        } else if (strcmp(arg, "--summary") == 0 && i < argc) {
            summary_path = argv[i++];
        } else if (strcmp(arg, "--assignments") == 0 && i < argc) {
            assignments_path = argv[i++];
        } else {
            usage(argv0);
            return 2;
        }
    }

    int paired_fastq = reads_path == NULL;
    if (left_path == NULL || right_path == NULL || out_path == NULL || left_len == 0 || right_len == 0 || k < 0 ||
        (reads_path != NULL && (left_reads_path != NULL || right_reads_path != NULL)) ||
        (paired_fastq && (left_reads_path == NULL || right_reads_path == NULL))) {
        fprintf(stderr, "pair-count requires --reads or both --left-reads and --right-reads\n");
        usage(argv0);
        return 2;
    }
    if (metric == COUNT_METRIC_HAMMING && k > 1) {
        fprintf(stderr, "--k 2 is only valid with --metric levenshtein\n");
        return 2;
    }

    seq_table left_targets = {0};
    seq_table right_targets = {0};
    const char **left_ptrs = NULL;
    const char **right_ptrs = NULL;
    size_t *left_lens = NULL;
    size_t *right_lens = NULL;
    qdaln_index *left_index = NULL;
    qdaln_index *right_index = NULL;
    fastq_reader left_reader = {0};
    fastq_reader right_reader = {0};
    FILE *out = NULL;
    FILE *summary = NULL;
    FILE *assignments = NULL;
    unsigned long long *pair_counts = NULL;
    pair_count_stats stats = {0};
    int rc = 1;

    if (read_target_table(left_path, &left_targets) != 0 || read_target_table(right_path, &right_targets) != 0) {
        fprintf(stderr, "failed to read pair target tables\n");
        goto done;
    }
    int left_id_check = validate_unique_seq_ids(&left_targets, "left target");
    int right_id_check = validate_unique_seq_ids(&right_targets, "right target");
    if (left_id_check != 0 || right_id_check != 0) {
        if (left_id_check == -1 || right_id_check == -1) fprintf(stderr, "out of memory\n");
        goto done;
    }
    if (metric == COUNT_METRIC_HAMMING &&
        (!all_targets_have_length(&left_targets, left_len) || !all_targets_have_length(&right_targets, right_len))) {
        fprintf(stderr, "--metric hamming requires targets to match their configured window lengths\n");
        goto done;
    }
    if (build_target_arrays(&left_targets, &left_ptrs, &left_lens) != 0 ||
        build_target_arrays(&right_targets, &right_ptrs, &right_lens) != 0) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    left_index = qdaln_index_build(left_ptrs, left_lens, left_targets.count);
    right_index = qdaln_index_build(right_ptrs, right_lens, right_targets.count);
    if (left_index == NULL || right_index == NULL) {
        fprintf(stderr, "failed to build pair target indexes\n");
        goto done;
    }
    size_t pair_count_slots = 0;
    if (checked_mul_size(left_targets.count, right_targets.count, &pair_count_slots) != 0) {
        fprintf(stderr, "pair count matrix is too large\n");
        goto done;
    }
    pair_counts = (unsigned long long *)calloc(alloc_count_or_one(pair_count_slots), sizeof(unsigned long long));
    if (pair_counts == NULL) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    if (fastq_reader_open(&left_reader, paired_fastq ? left_reads_path : reads_path) != 0 ||
        (paired_fastq && fastq_reader_open(&right_reader, right_reads_path) != 0)) {
        fprintf(stderr, "failed to open FASTQ input\n");
        goto done;
    }
    if (assignments_path != NULL) {
        assignments = open_output_file(assignments_path);
        if (assignments == NULL) {
            fprintf(stderr, "failed to open assignments output\n");
            goto done;
        }
        fprintf(assignments, "read_id\tleft_observed\tleft_index\tleft_id\tleft_status\tleft_distance\tright_observed\tright_index\tright_id\tright_status\tright_distance\tpair_status\n");
    }

    char left_header[8192];
    char left_seq[8192];
    char left_plus[8192];
    char left_qual[8192];
    char right_header[8192];
    char right_seq[8192];
    char right_plus[8192];
    char right_qual[8192];
    char read_id[8192];
    char right_read_id[8192];
    char left_observed[8192];
    char right_observed[8192];
    for (;;) {
        size_t left_seq_len = 0;
        size_t right_seq_len = 0;
        int left_got = fastq_read_record_len(&left_reader, left_header, left_seq, left_plus, left_qual,
                                             sizeof(left_header), &left_seq_len);
        if (left_got < 0) {
            fprintf(stderr, "malformed FASTQ input\n");
            goto done;
        }
        const char *right_seq_ptr = left_seq;
        if (paired_fastq) {
            int right_got = fastq_read_record_len(&right_reader, right_header, right_seq, right_plus, right_qual,
                                                  sizeof(right_header), &right_seq_len);
            if (right_got < 0) {
                fprintf(stderr, "malformed paired FASTQ input\n");
                goto done;
            }
            if (left_got != right_got) {
                fprintf(stderr, "paired FASTQ inputs have different record counts\n");
                goto done;
            }
            if (left_got == 0) break;
            fastq_pair_read_id(left_header, read_id, sizeof(read_id));
            fastq_pair_read_id(right_header, right_read_id, sizeof(right_read_id));
            if (read_id[0] == '\0' || right_read_id[0] == '\0') {
                fprintf(stderr, "paired FASTQ records require non-empty read IDs\n");
                goto done;
            }
            if (strcmp(read_id, right_read_id) != 0) {
                fprintf(stderr, "paired FASTQ read IDs do not match: %s != %s\n", read_id, right_read_id);
                goto done;
            }
            right_seq_ptr = right_seq;
        } else {
            if (left_got == 0) break;
            right_seq_len = left_seq_len;
            fastq_read_id(left_header, read_id, sizeof(read_id));
        }
        qdaln_match_result left = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
        qdaln_match_result right = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
        qdaln_index_stats left_stats = {0, 0};
        qdaln_index_stats right_stats = {0, 0};
        left_observed[0] = '\0';
        right_observed[0] = '\0';
        ++stats.total_reads;

        if (assign_count_window(left_index, left_seq, left_seq_len, left_start, left_len, k, metric, 0,
                                &left, &left_stats, left_observed, sizeof(left_observed), 0) != 0 ||
            assign_count_window(right_index, right_seq_ptr, right_seq_len, right_start, right_len, k, metric, 0,
                                &right, &right_stats, right_observed, sizeof(right_observed), 0) != 0) {
            fprintf(stderr, "FASTQ pair assignment failed\n");
            goto done;
        }
        apply_ambiguity_policy(&left, assignment_policy);
        apply_ambiguity_policy(&right, assignment_policy);
        stats.candidates_considered += left_stats.candidates_considered + right_stats.candidates_considered;
        stats.candidates_verified += left_stats.candidates_verified + right_stats.candidates_verified;

        const char *pair_status = pair_status_name(left, right);
        if (strcmp(pair_status, "unique") == 0) {
            size_t slot = (size_t)left.target_index * right_targets.count + (size_t)right.target_index;
            ++pair_counts[slot];
            ++stats.assigned_pairs;
        } else if (strcmp(pair_status, "invalid") == 0) {
            ++stats.invalid;
            if (left.status == QDALN_MATCH_INVALID) ++stats.left_invalid;
            if (right.status == QDALN_MATCH_INVALID) ++stats.right_invalid;
        } else {
            if (left.status == QDALN_MATCH_AMBIGUOUS || right.status == QDALN_MATCH_AMBIGUOUS) ++stats.pair_ambiguous;
            if (left.status == QDALN_MATCH_NONE) ++stats.left_unmatched;
            if (right.status == QDALN_MATCH_NONE) ++stats.right_unmatched;
        }
        if (assignments != NULL) {
            print_pair_assignment_row(assignments, read_id, &left_targets, &right_targets,
                                      left_observed, left, right_observed, right);
        }
    }

    out = open_output_file(out_path);
    if (out == NULL) {
        fprintf(stderr, "failed to open pair-count output\n");
        goto done;
    }
    fprintf(out, "left_id\tright_id\tcount\n");
    for (size_t li = 0; li < left_targets.count; ++li) {
        for (size_t ri = 0; ri < right_targets.count; ++ri) {
            unsigned long long count = pair_counts[li * right_targets.count + ri];
            if (count == 0) continue;
            fprintf(out, "%s\t%s\t%llu\n", left_targets.records[li].id, right_targets.records[ri].id, count);
        }
    }

    if (summary_path != NULL) {
        summary = open_output_file(summary_path);
        if (summary == NULL) {
            fprintf(stderr, "failed to open pair-count summary\n");
            goto done;
        }
        fprintf(summary,
                "{\n  \"workflow\": \"pair-count\",\n  \"input_mode\": \"%s\",\n  \"input_sync\": \"%s\",\n  \"k\": %d,\n  \"metric\": \"%s\",\n  \"ambiguity_policy\": \"%s\",\n  \"alphabet_policy\": \"%s\",\n  \"left_start\": %zu,\n  \"left_length\": %zu,\n  \"right_start\": %zu,\n  \"right_length\": %zu,\n  \"n_left_targets\": %zu,\n  \"n_right_targets\": %zu,\n  \"total_reads\": %llu,\n  \"total_pairs\": %llu,\n  \"assigned_pairs\": %llu,\n  \"pair_ambiguous\": %llu,\n  \"left_unmatched\": %llu,\n  \"right_unmatched\": %llu,\n  \"invalid\": %llu,\n  \"left_invalid\": %llu,\n  \"right_invalid\": %llu,\n  \"candidates_considered\": %llu,\n  \"candidates_verified\": %llu\n}\n",
                paired_fastq ? "paired-fastq" : "single-read", paired_fastq ? "canonical-read-id" : "not-applicable",
                k, metric_name(metric), ambiguity_policy_name(assignment_policy), qdaln_alphabet_policy(), left_start, left_len, right_start, right_len,
                left_targets.count, right_targets.count, stats.total_reads, stats.total_reads, stats.assigned_pairs,
                stats.pair_ambiguous, stats.left_unmatched, stats.right_unmatched, stats.invalid,
                stats.left_invalid, stats.right_invalid, stats.candidates_considered, stats.candidates_verified);
    }

    rc = 0;

done:
    if (out != NULL) fclose(out);
    if (summary != NULL) fclose(summary);
    if (assignments != NULL) fclose(assignments);
    fastq_reader_close(&left_reader);
    fastq_reader_close(&right_reader);
    qdaln_index_free(left_index);
    qdaln_index_free(right_index);
    free(left_ptrs);
    free(right_ptrs);
    free(left_lens);
    free(right_lens);
    free(pair_counts);
    free_table(&left_targets);
    free_table(&right_targets);
    return rc;
}

static FILE *open_demux_target_file(FILE **files, const seq_table *targets, size_t target_index, const char *out_dir) {
    if (files[target_index] != NULL) return files[target_index];
    char safe_id[512];
    sanitize_filename(targets->records[target_index].id, safe_id, sizeof(safe_id));
    char name[600];
    int n = snprintf(name, sizeof(name), "%s.fastq", safe_id);
    if (n < 0 || (size_t)n >= sizeof(name)) return NULL;
    char path[4096];
    if (path_join(path, sizeof(path), out_dir, name) != 0) return NULL;
    files[target_index] = open_output_file(path);
    return files[target_index];
}

static int run_demux(const char *argv0, int argc, char **argv) {
    const char *barcodes_path = NULL;
    const char *reads_path = NULL;
    const char *out_dir = NULL;
    const char *summary_path = NULL;
    const char *assignments_path = NULL;
    const char *ambiguous_path = NULL;
    const char *unmatched_path = NULL;
    count_metric metric = COUNT_METRIC_LEVENSHTEIN;
    size_t barcode_start = 0;
    size_t barcode_len = 0;
    int auto_barcode_len = 0;
    size_t indel_window = 0;
    int max_correction_qual = -1;
    int k = -1;
    ambiguity_policy assignment_policy = AMBIGUITY_POLICY_RADIUS;

    int i = 2;
    while (i < argc) {
        const char *arg = argv[i++];
        if ((strcmp(arg, "--barcodes") == 0 || strcmp(arg, "--targets") == 0) && i < argc) {
            barcodes_path = argv[i++];
        } else if (strcmp(arg, "--reads") == 0 && i < argc) {
            reads_path = argv[i++];
        } else if ((strcmp(arg, "--barcode-start") == 0 || strcmp(arg, "--target-start") == 0) && i < argc) {
            if (parse_size_value(argv[i++], &barcode_start) != 0) {
                usage(argv0);
                return 2;
            }
        } else if ((strcmp(arg, "--barcode-length") == 0 || strcmp(arg, "--target-length") == 0) && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "auto") == 0) {
                auto_barcode_len = 1;
                barcode_len = 0;
            } else if (parse_size_value(value, &barcode_len) != 0 || barcode_len == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--k") == 0 && i < argc) {
            if (parse_int_value(argv[i++], &k) != 0 || k < 0 || k > 2) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--metric") == 0 && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "hamming") == 0) {
                metric = COUNT_METRIC_HAMMING;
            } else if (strcmp(value, "levenshtein") == 0) {
                metric = COUNT_METRIC_LEVENSHTEIN;
            } else {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--ambiguity-policy") == 0 && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "radius") == 0) {
                assignment_policy = AMBIGUITY_POLICY_RADIUS;
            } else if (strcmp(value, "best") == 0) {
                assignment_policy = AMBIGUITY_POLICY_BEST;
            } else {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--indel-window") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &indel_window) != 0 || indel_window > 1) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--max-correction-qual") == 0 && i < argc) {
            if (parse_int_value(argv[i++], &max_correction_qual) != 0 ||
                max_correction_qual < 0 || max_correction_qual > 93) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--out-dir") == 0 && i < argc) {
            out_dir = argv[i++];
        } else if (strcmp(arg, "--summary") == 0 && i < argc) {
            summary_path = argv[i++];
        } else if (strcmp(arg, "--qc") == 0 && i < argc) {
            summary_path = argv[i++];
        } else if (strcmp(arg, "--assignments") == 0 && i < argc) {
            assignments_path = argv[i++];
        } else if (strcmp(arg, "--ambiguous-out") == 0 && i < argc) {
            ambiguous_path = argv[i++];
        } else if (strcmp(arg, "--unmatched-out") == 0 && i < argc) {
            unmatched_path = argv[i++];
        } else {
            usage(argv0);
            return 2;
        }
    }

    if (barcodes_path == NULL || reads_path == NULL || out_dir == NULL || (barcode_len == 0 && !auto_barcode_len) || k < 0) {
        usage(argv0);
        return 2;
    }
    if (metric == COUNT_METRIC_HAMMING && indel_window != 0) {
        fprintf(stderr, "--indel-window is only valid with --metric levenshtein\n");
        return 2;
    }
    if (metric == COUNT_METRIC_HAMMING && k > 1) {
        fprintf(stderr, "--k 2 is only valid with --metric levenshtein\n");
        return 2;
    }
    if (indel_window != 0 && k != 1) {
        fprintf(stderr, "--indel-window requires --k 1\n");
        return 2;
    }

    seq_table targets = {0};
    fastq_reader reader = {0};
    qdaln_index *index = NULL;
    hamming_lookup hlookup = {0};
    levenshtein1_lookup levlookup = {0};
    const char **target_ptrs = NULL;
    size_t *target_lens = NULL;
    size_t *auto_barcode_lens = NULL;
    size_t auto_barcode_lens_count = 0;
    size_t fixed_barcode_lens[1] = {0};
    const size_t *barcode_lens = NULL;
    size_t barcode_lens_count = 0;
    FILE **target_files = NULL;
    FILE *assignments = NULL;
    FILE *ambiguous_out = NULL;
    FILE *unmatched_out = NULL;
    count_stats stats = {0};
    unsigned long long *target_counts = NULL;
    const char *assignment_engine = "generic_indexed";
    int rc = 1;

    if (read_target_table(barcodes_path, &targets) != 0) {
        fprintf(stderr, "failed to read barcodes\n");
        goto done;
    }
    int barcode_id_check = validate_unique_seq_ids(&targets, "barcode");
    if (barcode_id_check != 0) {
        if (barcode_id_check == -1) fprintf(stderr, "out of memory\n");
        goto done;
    }
    if (!auto_barcode_len && metric == COUNT_METRIC_HAMMING && !all_targets_have_length(&targets, barcode_len)) {
        fprintf(stderr, "--metric hamming requires every barcode to have --barcode-length bases\n");
        goto done;
    }
    if (auto_barcode_len) {
        if (collect_target_lengths(&targets, &auto_barcode_lens, &auto_barcode_lens_count) != 0) {
            fprintf(stderr, "out of memory\n");
            goto done;
        }
        barcode_lens = auto_barcode_lens;
        barcode_lens_count = auto_barcode_lens_count;
    } else {
        fixed_barcode_lens[0] = barcode_len;
        barcode_lens = fixed_barcode_lens;
        barcode_lens_count = 1;
    }
    if (build_target_arrays(&targets, &target_ptrs, &target_lens) != 0) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    index = qdaln_index_build(target_ptrs, target_lens, targets.count);
    if (index == NULL) {
        fprintf(stderr, "failed to build barcode index\n");
        goto done;
    }
    if (!auto_barcode_len && metric == COUNT_METRIC_HAMMING && indel_window == 0 && (k == 0 || k == 1) &&
        barcode_len <= 32) {
        int lookup_rc = k == 0 ? build_hamming_exact_lookup(&targets, barcode_len, &hlookup)
                               : build_hamming_lookup(&targets, barcode_len, &hlookup);
        if (lookup_rc < 0) {
            fprintf(stderr, "failed to build barcode Hamming lookup\n");
            goto done;
        }
        if (hlookup.ready) assignment_engine = k == 0 ? "hamming_exact_lookup_direct" : "hamming_k1_lookup_direct";
    }
    if (!auto_barcode_len && metric == COUNT_METRIC_LEVENSHTEIN && indel_window == 1 && k == 1 &&
        assignment_policy == AMBIGUITY_POLICY_BEST && barcode_len <= 31) {
        int lookup_rc = build_levenshtein1_lookup(&targets, barcode_len, &levlookup);
        if (lookup_rc < 0) {
            fprintf(stderr, "failed to build barcode Levenshtein lookup\n");
            goto done;
        }
        if (levlookup.ready) assignment_engine = "levenshtein_k1_lookup_direct";
    }
    int filename_check = validate_unique_sanitized_filenames(&targets);
    if (filename_check != 0) {
        if (filename_check == -1) fprintf(stderr, "out of memory\n");
        goto done;
    }
    if (ensure_dir(out_dir) != 0) {
        fprintf(stderr, "failed to create output directory\n");
        goto done;
    }
    target_files = (FILE **)calloc(targets.count == 0 ? 1 : targets.count, sizeof(FILE *));
    target_counts = (unsigned long long *)calloc(targets.count == 0 ? 1 : targets.count, sizeof(unsigned long long));
    if (target_files == NULL || target_counts == NULL) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    if (fastq_reader_open(&reader, reads_path) != 0) {
        fprintf(stderr, "failed to open FASTQ input\n");
        goto done;
    }
    if (assignments_path != NULL) {
        assignments = open_output_file(assignments_path);
        if (assignments == NULL) {
            fprintf(stderr, "failed to open assignments output\n");
            goto done;
        }
        fprintf(assignments, "read_id\tobserved_barcode\ttarget_index\ttarget_id\ttarget_seq\tbest_distance\tsecond_best_distance\tmatch_count\tstatus\n");
    }
    if (ambiguous_path != NULL) {
        ambiguous_out = open_output_file(ambiguous_path);
        if (ambiguous_out == NULL) {
            fprintf(stderr, "failed to open ambiguous FASTQ output\n");
            goto done;
        }
    }
    if (unmatched_path != NULL) {
        unmatched_out = open_output_file(unmatched_path);
        if (unmatched_out == NULL) {
            fprintf(stderr, "failed to open unmatched FASTQ output\n");
            goto done;
        }
    }

    char header[8192];
    char seq[8192];
    char plus[8192];
    char qual[8192];
    char read_id[8192];
    char observed[8192];
    int got = 0;
    size_t seq_len = 0;
    while ((got = fastq_read_record_len(&reader, header, seq, plus, qual, sizeof(header), &seq_len)) == 1) {
        fastq_read_id(header, read_id, sizeof(read_id));
        qdaln_match_result result = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
        qdaln_index_stats istats = {0, 0};
        observed[0] = '\0';
        ++stats.total;

        int handled = 0;
        if (hlookup.ready) {
            int exact_merge = assignment_policy == AMBIGUITY_POLICY_RADIUS ||
                    assignments != NULL || ambiguous_out != NULL || unmatched_out != NULL;
            int lookup_rc = assign_hamming_lookup_offsets(&hlookup, seq, seq_len, NULL, barcode_start, k,
                                                          &result, &istats, observed, sizeof(observed),
                                                          exact_merge);
            if (lookup_rc < 0) {
                fprintf(stderr, "FASTQ assignment failed\n");
                goto done;
            }
            handled = lookup_rc;
        }
        if (!handled && levlookup.ready) {
            int lookup_rc = assign_levenshtein1_lookup_offset(&levlookup, seq, seq_len, barcode_start, &result,
                                                              &istats, observed, sizeof(observed));
            if (lookup_rc < 0) {
                fprintf(stderr, "FASTQ assignment failed\n");
                goto done;
            }
            handled = lookup_rc;
        }
        if (!handled &&
            assign_count_length_set(index, seq, seq_len, barcode_start, barcode_lens, barcode_lens_count, k, metric,
                                    indel_window, &result, &istats, observed, sizeof(observed), 0) != 0) {
            fprintf(stderr, "FASTQ assignment failed\n");
            goto done;
        }
        apply_ambiguity_policy(&result, assignment_policy);
        if (result.status == QDALN_MATCH_UNIQUE && result.target_index >= 0 && result.best_distance > 0) {
            seq_record *target = &targets.records[result.target_index];
            offset_list barcode_offset = {0};
            if (push_offset_unique(&barcode_offset, barcode_start) != 0) {
                free_offset_list(&barcode_offset);
                fprintf(stderr, "out of memory\n");
                goto done;
            }
            if (!quality_allows_unique_correction(seq, seq_len, qual, &barcode_offset, barcode_start,
                                                  target->len, metric, indel_window, k, observed, target, result,
                                                  max_correction_qual)) {
                result = (qdaln_match_result){-1, -1, -1, 0, QDALN_MATCH_NONE};
            }
            free_offset_list(&barcode_offset);
        }
        if (result.status != QDALN_MATCH_INVALID) {
            stats.candidates_considered += (unsigned long long)istats.candidates_considered;
            stats.candidates_verified += (unsigned long long)istats.candidates_verified;
        }

        if (assignments != NULL) print_fastq_row(assignments, &targets, read_id, observed, result);

        if (result.status == QDALN_MATCH_UNIQUE && result.target_index >= 0) {
            FILE *target_out = open_demux_target_file(target_files, &targets, (size_t)result.target_index, out_dir);
            if (target_out == NULL) {
                fprintf(stderr, "failed to open per-barcode FASTQ output\n");
                goto done;
            }
            write_fastq_record(target_out, header, seq, plus, qual);
            ++target_counts[result.target_index];
            ++stats.unique;
            if (result.best_distance == 0) ++stats.exact;
            else ++stats.corrected;
        } else if (result.status == QDALN_MATCH_AMBIGUOUS) {
            ++stats.ambiguous;
            if (ambiguous_out != NULL) write_fastq_record(ambiguous_out, header, seq, plus, qual);
        } else if (result.status == QDALN_MATCH_NONE) {
            ++stats.unmatched;
            if (unmatched_out != NULL) write_fastq_record(unmatched_out, header, seq, plus, qual);
        } else {
            ++stats.invalid;
            if (unmatched_out != NULL) write_fastq_record(unmatched_out, header, seq, plus, qual);
        }
    }
    if (got < 0) {
        fprintf(stderr, "malformed FASTQ input\n");
        goto done;
    }

    if (summary_path != NULL) {
        FILE *summary = open_output_file(summary_path);
        if (summary == NULL) {
            fprintf(stderr, "failed to open summary output\n");
            goto done;
        }
        unsigned long long top_count = 0;
        size_t top_target = 0;
        size_t nonempty = 0;
        for (size_t i = 0; i < targets.count; ++i) {
            if (target_counts[i] != 0) ++nonempty;
            if (target_counts[i] > top_count) {
                top_count = target_counts[i];
                top_target = i;
            }
        }
        fprintf(summary,
                "{\n  \"workflow\": \"demux\",\n  \"k\": %d,\n  \"metric\": \"%s\",\n  \"ambiguity_policy\": \"%s\",\n  \"assignment_engine\": \"%s\",\n  \"alphabet_policy\": \"%s\",\n  \"max_correction_qual\": ",
                k, metric_name(metric), ambiguity_policy_name(assignment_policy), assignment_engine,
                qdaln_alphabet_policy());
        if (max_correction_qual >= 0) {
            fprintf(summary, "%d", max_correction_qual);
        } else {
            fprintf(summary, "null");
        }
        fprintf(summary,
                ",\n  \"indel_window\": %zu,\n  \"barcode_start\": %zu,\n  \"barcode_length\": %zu,\n  \"barcode_length_mode\": \"%s\",\n  \"barcode_lengths\": [",
                indel_window, barcode_start, barcode_len, auto_barcode_len ? "auto" : "fixed");
        for (size_t i = 0; i < barcode_lens_count; ++i) {
            fprintf(summary, "%s%zu", i == 0 ? "" : ", ", barcode_lens[i]);
        }
        fprintf(summary,
                "],\n  \"n_barcodes\": %zu,\n  \"total_reads\": %llu,\n  \"assigned_unique\": %llu,\n  \"assigned_exact\": %llu,\n  \"assigned_corrected\": %llu,\n  \"ambiguous\": %llu,\n  \"unmatched\": %llu,\n  \"invalid\": %llu,\n  \"nonempty_outputs\": %zu,\n  \"top_barcode_id\": \"%s\",\n  \"top_barcode_count\": %llu,\n  \"candidates_considered\": %llu,\n  \"candidates_verified\": %llu\n}\n",
                targets.count, stats.total,
                stats.unique, stats.exact, stats.corrected, stats.ambiguous, stats.unmatched, stats.invalid,
                nonempty, targets.count == 0 ? "" : targets.records[top_target].id, top_count,
                stats.candidates_considered, stats.candidates_verified);
        fclose(summary);
    }

    rc = 0;

done:
    if (target_files != NULL) {
        for (size_t i = 0; i < targets.count; ++i) {
            if (target_files[i] != NULL) fclose(target_files[i]);
        }
    }
    if (assignments != NULL) fclose(assignments);
    if (ambiguous_out != NULL) fclose(ambiguous_out);
    if (unmatched_out != NULL) fclose(unmatched_out);
    fastq_reader_close(&reader);
    qdaln_index_free(index);
    free_hamming_lookup(&hlookup);
    free_levenshtein1_lookup(&levlookup);
    free(target_ptrs);
    free(target_lens);
    free(auto_barcode_lens);
    free(target_files);
    free(target_counts);
    free_table(&targets);
    return rc;
}

typedef struct bcl_sample {
    char *id;
    char *name;
    char *index1;
    char *index2;
    int lane;
    size_t output_index;
    int is_alias;
    unsigned long long assigned;
} bcl_sample;

typedef struct bcl_sample_table {
    bcl_sample *items;
    size_t count;
    size_t cap;
} bcl_sample_table;

typedef struct bcl_read_info {
    int number;
    size_t cycles;
    int indexed;
    size_t start_cycle;
} bcl_read_info;

typedef struct bcl_run_info {
    bcl_read_info reads[16];
    size_t read_count;
    size_t total_cycles;
} bcl_run_info;

typedef struct bcl_unknown_barcode {
    char *index;
    unsigned long long count;
} bcl_unknown_barcode;

typedef struct bcl_unknown_table {
    bcl_unknown_barcode *items;
    size_t count;
    size_t cap;
} bcl_unknown_table;

typedef struct text_buffer {
    char *data;
    size_t len;
    size_t cap;
} text_buffer;

static void free_bcl_unknowns(bcl_unknown_table *table) {
    for (size_t i = 0; i < table->count; ++i) free(table->items[i].index);
    free(table->items);
    table->items = NULL;
    table->count = 0;
    table->cap = 0;
}

static int add_bcl_unknown_count(bcl_unknown_table *table, const char *index, unsigned long long count) {
    if (count == 0) return 0;
    for (size_t i = 0; i < table->count; ++i) {
        if (strcmp(table->items[i].index, index) == 0) {
            table->items[i].count += count;
            return 0;
        }
    }
    if (table->count == table->cap) {
        size_t next_cap = table->cap == 0 ? 64 : table->cap * 2;
        bcl_unknown_barcode *next = (bcl_unknown_barcode *)realloc(table->items, next_cap * sizeof(bcl_unknown_barcode));
        if (next == NULL) return -1;
        table->items = next;
        table->cap = next_cap;
    }
    table->items[table->count].index = xstrndup(index, strlen(index));
    if (table->items[table->count].index == NULL) return -1;
    table->items[table->count].count = count;
    ++table->count;
    return 0;
}

static int add_bcl_unknown(bcl_unknown_table *table, const char *index) {
    return add_bcl_unknown_count(table, index, 1);
}

static int cmp_bcl_unknown_desc(const void *a, const void *b) {
    const bcl_unknown_barcode *aa = (const bcl_unknown_barcode *)a;
    const bcl_unknown_barcode *bb = (const bcl_unknown_barcode *)b;
    if (aa->count < bb->count) return 1;
    if (aa->count > bb->count) return -1;
    return strcmp(aa->index, bb->index);
}

static int merge_bcl_unknowns(bcl_unknown_table *dst, const bcl_unknown_table *src) {
    for (size_t i = 0; i < src->count; ++i) {
        if (add_bcl_unknown_count(dst, src->items[i].index, src->items[i].count) != 0) return -1;
    }
    return 0;
}

static void free_text_buffer(text_buffer *buf) {
    free(buf->data);
    buf->data = NULL;
    buf->len = 0;
    buf->cap = 0;
}

static int text_buffer_reserve(text_buffer *buf, size_t extra) {
    if (extra > SIZE_MAX - buf->len) return -1;
    size_t need = buf->len + extra;
    if (need <= buf->cap) return 0;
    size_t next = buf->cap == 0 ? 65536 : buf->cap;
    while (next < need) {
        if (next > SIZE_MAX / 2) {
            next = need;
            break;
        }
        next *= 2;
    }
    char *p = (char *)realloc(buf->data, next);
    if (p == NULL) return -1;
    buf->data = p;
    buf->cap = next;
    return 0;
}

static int text_buffer_append(text_buffer *buf, const char *s, size_t n) {
    if (text_buffer_reserve(buf, n) != 0) return -1;
    memcpy(buf->data + buf->len, s, n);
    buf->len += n;
    return 0;
}

static int gzwrite_all(gzFile gz, const char *data, size_t len) {
    size_t written = 0;
    while (written < len) {
        size_t remaining = len - written;
        unsigned int chunk = remaining > (size_t)UINT_MAX ? UINT_MAX : (unsigned int)remaining;
        int rc = gzwrite(gz, data + written, chunk);
        if (rc <= 0 || (unsigned int)rc != chunk) return -1;
        written += (size_t)rc;
    }
    return 0;
}

static void free_bcl_samples(bcl_sample_table *table) {
    for (size_t i = 0; i < table->count; ++i) {
        free(table->items[i].id);
        free(table->items[i].name);
        free(table->items[i].index1);
        free(table->items[i].index2);
    }
    free(table->items);
    table->items = NULL;
    table->count = 0;
    table->cap = 0;
}

static int push_bcl_sample(bcl_sample_table *table, const char *id, const char *name,
                           const char *index1, const char *index2, int lane) {
    size_t output_index = table->count;
    int is_alias = 0;
    for (size_t i = 0; i < table->count; ++i) {
        bcl_sample *existing = &table->items[i];
        if (strcmp(existing->id, id) == 0 && existing->lane == lane) {
            output_index = existing->output_index;
            is_alias = 1;
            break;
        }
    }
    if (table->count == table->cap) {
        size_t next_cap = table->cap == 0 ? 16 : table->cap * 2;
        bcl_sample *next = (bcl_sample *)realloc(table->items, next_cap * sizeof(bcl_sample));
        if (next == NULL) return -1;
        table->items = next;
        table->cap = next_cap;
    }
    bcl_sample *s = &table->items[table->count];
    memset(s, 0, sizeof(*s));
    s->id = xstrndup(id, strlen(id));
    s->name = xstrndup(name != NULL && name[0] != '\0' ? name : id, strlen(name != NULL && name[0] != '\0' ? name : id));
    s->index1 = xstrndup(index1, strlen(index1));
    s->index2 = xstrndup(index2 != NULL ? index2 : "", strlen(index2 != NULL ? index2 : ""));
    s->lane = lane;
    s->output_index = output_index;
    s->is_alias = is_alias;
    if (s->id == NULL || s->name == NULL || s->index1 == NULL || s->index2 == NULL) return -1;
    uppercase_ascii(s->index1);
    uppercase_ascii(s->index2);
    ++table->count;
    return 0;
}

static char *read_text_file(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return NULL;
    if (fseek(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        return NULL;
    }
    long n = ftell(fp);
    if (n < 0) {
        fclose(fp);
        return NULL;
    }
    rewind(fp);
    char *text = (char *)malloc((size_t)n + 1);
    if (text == NULL) {
        fclose(fp);
        return NULL;
    }
    if (fread(text, 1, (size_t)n, fp) != (size_t)n) {
        free(text);
        fclose(fp);
        return NULL;
    }
    text[n] = '\0';
    fclose(fp);
    return text;
}

static int xml_attr_value(const char *tag, const char *name, char *out, size_t out_cap) {
    char pattern[64];
    snprintf(pattern, sizeof(pattern), "%s=\"", name);
    const char *p = strstr(tag, pattern);
    if (p == NULL) return -1;
    p += strlen(pattern);
    const char *end = strchr(p, '"');
    if (end == NULL) return -1;
    size_t n = (size_t)(end - p);
    if (n >= out_cap) n = out_cap - 1;
    memcpy(out, p, n);
    out[n] = '\0';
    return 0;
}

static int parse_run_info(const char *run_folder, bcl_run_info *info) {
    char path[4096];
    if (path_join(path, sizeof(path), run_folder, "RunInfo.xml") != 0) return -1;
    char *xml = read_text_file(path);
    if (xml == NULL) return -1;
    memset(info, 0, sizeof(*info));
    const char *p = xml;
    while ((p = strstr(p, "<Read")) != NULL) {
        const char *end = strchr(p, '>');
        if (end == NULL) {
            free(xml);
            return -1;
        }
        char tag[1024];
        size_t tag_len = (size_t)(end - p + 1);
        if (tag_len >= sizeof(tag)) tag_len = sizeof(tag) - 1;
        memcpy(tag, p, tag_len);
        tag[tag_len] = '\0';
        char number[32] = "";
        char cycles[32] = "";
        char indexed[32] = "";
        if (xml_attr_value(tag, "Number", number, sizeof(number)) == 0 &&
            xml_attr_value(tag, "NumCycles", cycles, sizeof(cycles)) == 0 &&
            xml_attr_value(tag, "IsIndexedRead", indexed, sizeof(indexed)) == 0) {
            if (info->read_count >= 16) {
                free(xml);
                return -1;
            }
            int parsed_number = 0;
            size_t parsed_cycles = 0;
            if (parse_int_value(number, &parsed_number) != 0 || parsed_number < 1 ||
                parse_size_value(cycles, &parsed_cycles) != 0 || parsed_cycles == 0 ||
                parsed_cycles > MAX_BCL_READ_CYCLES ||
                info->total_cycles > MAX_BCL_TOTAL_CYCLES - parsed_cycles) {
                free(xml);
                return -1;
            }
            if (!(indexed[0] == 'Y' || indexed[0] == 'y' || indexed[0] == 'N' || indexed[0] == 'n') ||
                indexed[1] != '\0') {
                free(xml);
                return -1;
            }
            bcl_read_info *r = &info->reads[info->read_count++];
            r->number = parsed_number;
            r->cycles = parsed_cycles;
            r->indexed = indexed[0] == 'Y' || indexed[0] == 'y';
            r->start_cycle = info->total_cycles + 1;
            info->total_cycles += r->cycles;
        }
        p = end + 1;
    }
    free(xml);
    return info->read_count == 0 ? -1 : 0;
}

static int read_bcl_sample_sheet(const char *path, bcl_sample_table *samples) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return -1;
    char buf[16384];
    int in_data = 0;
    int have_header = 0;
    int id_col = -1;
    int name_col = -1;
    int index_col = -1;
    int index2_col = -1;
    int lane_col = -1;
    size_t row = 0;
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        ++row;
        trim_line(buf);
        if (buf[0] == '\0') continue;
        if (buf[0] == '[') {
            in_data = field_eq(buf, "[Data]") || field_eq(buf, "[BCLConvert_Data]");
            have_header = 0;
            continue;
        }
        if (!in_data) continue;
        char *fields[64];
        size_t nf = split_fields(buf, ',', fields, 64);
        if (!have_header) {
            id_col = find_column(fields, nf, "Sample_ID", "SampleID", "Sample_ID");
            if (id_col < 0) id_col = find_column(fields, nf, "Sample_Project", "SampleName", "Sample_Name");
            name_col = find_column(fields, nf, "Sample_Name", "SampleName", "sample_name");
            index_col = find_column(fields, nf, "index", "Index", "Index1");
            index2_col = find_column(fields, nf, "index2", "Index2", "index_2");
            lane_col = find_column(fields, nf, "Lane", "lane", NULL);
            have_header = 1;
            if (id_col < 0 || index_col < 0) {
                fclose(fp);
                return -1;
            }
            continue;
        }
        if ((size_t)id_col >= nf || (size_t)index_col >= nf || fields[id_col][0] == '\0' || fields[index_col][0] == '\0') {
            fprintf(stderr, "%s:%zu: BCL sample sheet Sample_ID and index must be non-empty\n", path, row);
            fclose(fp);
            return -1;
        }
        const char *name = (name_col >= 0 && (size_t)name_col < nf) ? fields[name_col] : fields[id_col];
        const char *index2 = (index2_col >= 0 && (size_t)index2_col < nf) ? fields[index2_col] : "";
        int lane = 0;
        if (lane_col >= 0 && (size_t)lane_col < nf && fields[lane_col][0] != '\0') {
            if (parse_int_value(fields[lane_col], &lane) != 0 || lane < 1) {
                fprintf(stderr, "%s:%zu: BCL sample sheet Lane must be a positive integer\n", path, row);
                fclose(fp);
                return -1;
            }
            if (lane != 1) {
                fprintf(stderr, "%s:%zu: classic BCL demux currently supports sample-sheet Lane 1 only\n", path, row);
                fclose(fp);
                return -1;
            }
        }
        if (push_bcl_sample(samples, fields[id_col], name, fields[index_col], index2, lane) != 0) {
            fclose(fp);
            return -1;
        }
    }
    fclose(fp);
    return samples->count == 0 ? -1 : 0;
}

static int hamming_distance_limit(const char *a, const char *b, int limit) {
    size_t na = strlen(a);
    size_t nb = strlen(b);
    if (na != nb) return limit + 1;
    int d = 0;
    for (size_t i = 0; i < na; ++i) {
        if (a[i] != b[i] && ++d > limit) return d;
    }
    return d;
}

static int assign_bcl_sample(const bcl_sample_table *samples, int lane, const char *index1, const char *index2,
                             int k1, int k2, int *match_count_out) {
    int best = -1;
    int best_d = 1000000;
    size_t best_output = (size_t)-1;
    int matches = 0;
    for (size_t i = 0; i < samples->count; ++i) {
        const bcl_sample *s = &samples->items[i];
        if (s->lane != 0 && s->lane != lane) continue;
        int d1 = hamming_distance_limit(index1, s->index1, k1);
        if (d1 > k1) continue;
        int d2 = 0;
        if (s->index2[0] != '\0' || index2[0] != '\0') {
            d2 = hamming_distance_limit(index2, s->index2, k2);
            if (d2 > k2) continue;
        }
        int d = d1 + d2;
        if (d < best_d) {
            best_d = d;
            best = (int)i;
            best_output = s->output_index;
            matches = 1;
        } else if (d == best_d && s->output_index != best_output) {
            ++matches;
        }
    }
    *match_count_out = matches;
    return matches == 1 ? best : -1;
}

static int read_u32_le(const unsigned char *p) {
    return (int)((unsigned int)p[0] | ((unsigned int)p[1] << 8) | ((unsigned int)p[2] << 16) | ((unsigned int)p[3] << 24));
}

static int path_exists(const char *path) {
    struct stat st;
    return stat(path, &st) == 0;
}

static int build_bcl_path(char *out, size_t out_cap, const char *basecalls, int lane, size_t cycle, const char *tile) {
    int n = snprintf(out, out_cap, "%s/L%03d/C%zu.1/s_%d_%s.bcl.gz", basecalls, lane, cycle, lane, tile);
    if (n < 0 || (size_t)n >= out_cap) return -1;
    if (path_exists(out)) return 0;
    n = snprintf(out, out_cap, "%s/L%03d/C%zu.1/s_%d_%s.bcl", basecalls, lane, cycle, lane, tile);
    if (n < 0 || (size_t)n >= out_cap) return -1;
    return path_exists(out) ? 0 : -1;
}

static int read_bcl_cycle(const char *path, unsigned char **bytes_out, size_t *count_out) {
    gzFile gz = gzopen(path, "rb");
    if (gz == NULL) return -1;
    unsigned char header[4];
    if (gzread(gz, header, 4) != 4) {
        gzclose(gz);
        return -1;
    }
    int n = read_u32_le(header);
    if (n < 0 || n > MAX_BCL_CYCLE_CLUSTERS) {
        gzclose(gz);
        return -1;
    }
    unsigned char *bytes = (unsigned char *)malloc((size_t)n == 0 ? 1 : (size_t)n);
    if (bytes == NULL) {
        gzclose(gz);
        return -1;
    }
    if (n != 0 && gzread(gz, bytes, (unsigned int)n) != n) {
        free(bytes);
        gzclose(gz);
        return -1;
    }
    gzclose(gz);
    *bytes_out = bytes;
    *count_out = (size_t)n;
    return 0;
}

static int read_filter_file(const char *basecalls, int lane, const char *tile, unsigned char **pf_out, size_t *count_out) {
    char path[4096];
    int n = snprintf(path, sizeof(path), "%s/L%03d/s_%d_%s.filter", basecalls, lane, lane, tile);
    if (n < 0 || (size_t)n >= sizeof(path)) return -1;
    if (!path_exists(path)) {
        *pf_out = NULL;
        *count_out = 0;
        return 0;
    }
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    unsigned char header[8];
    if (fread(header, 1, 8, fp) != 8) {
        fclose(fp);
        return -1;
    }
    int count = read_u32_le(header + 4);
    if (count < 0 || count > MAX_BCL_CYCLE_CLUSTERS) {
        fclose(fp);
        return -1;
    }
    unsigned char *pf = (unsigned char *)malloc((size_t)count == 0 ? 1 : (size_t)count);
    if (pf == NULL) {
        fclose(fp);
        return -1;
    }
    if (count != 0 && fread(pf, 1, (size_t)count, fp) != (size_t)count) {
        free(pf);
        fclose(fp);
        return -1;
    }
    fclose(fp);
    *pf_out = pf;
    *count_out = (size_t)count;
    return 0;
}

static char bcl_base(unsigned char b) {
    if (b == 0) return 'N';
    switch (b & 3u) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        default: return 'T';
    }
}

static char bcl_qual(unsigned char b) {
    if (b == 0) return '#';
    unsigned int q = b >> 2;
    if (q > 93) q = 93;
    return (char)(q + 33);
}

static int collect_tiles(const char *basecalls, int lane, char ***tiles_out, size_t *tile_count_out) {
    char path[4096];
    int n = snprintf(path, sizeof(path), "%s/L%03d/C1.1", basecalls, lane);
    if (n < 0 || (size_t)n >= sizeof(path)) return -1;
    DIR *dir = opendir(path);
    if (dir == NULL) return -1;
    string_list tiles = {0};
    struct dirent *ent;
    char prefix[32];
    n = snprintf(prefix, sizeof(prefix), "s_%d_", lane);
    if (n < 0 || (size_t)n >= sizeof(prefix)) {
        closedir(dir);
        return -1;
    }
    while ((ent = readdir(dir)) != NULL) {
        if (strncmp(ent->d_name, prefix, strlen(prefix)) != 0) continue;
        char *start = ent->d_name + strlen(prefix);
        char *bcl = strstr(start, ".bcl");
        if (bcl == NULL) continue;
        size_t len = (size_t)(bcl - start);
        char tile[128];
        if (len == 0 || len >= sizeof(tile)) continue;
        memcpy(tile, start, len);
        tile[len] = '\0';
        if (push_string(&tiles, tile) != 0) {
            closedir(dir);
            free_string_list(&tiles);
            return -1;
        }
    }
    closedir(dir);
    *tiles_out = tiles.items;
    *tile_count_out = tiles.count;
    return tiles.count == 0 ? -1 : 0;
}

static gzFile open_bcl_fastq(const char *out_dir, const char *sample_id, size_t sample_number, int lane,
                             char read_kind, int read_number, int gzip_level) {
    char safe_id[512];
    sanitize_filename(sample_id, safe_id, sizeof(safe_id));
    char name[700];
    int n = snprintf(name, sizeof(name), "%s_S%zu_L%03d_%c%d_001.fastq.gz", safe_id, sample_number,
                     lane, read_kind, read_number);
    if (n < 0 || (size_t)n >= sizeof(name)) return NULL;
    char path[4096];
    if (path_join(path, sizeof(path), out_dir, name) != 0) return NULL;
    char mode[8];
    snprintf(mode, sizeof(mode), "wb%d", gzip_level);
    gzFile gz = gzopen(path, mode);
    if (gz != NULL) gzbuffer(gz, 1024 * 1024);
    return gz;
}

static int bcl_output_enabled(const bcl_read_info *read, int emit_index_fastqs) {
    return !read->indexed || emit_index_fastqs;
}

static char bcl_output_kind(const bcl_read_info *read) {
    return read->indexed ? 'I' : 'R';
}

static int bcl_output_number(const bcl_run_info *run, size_t read_i) {
    int n = 0;
    for (size_t i = 0; i <= read_i && i < run->read_count; ++i) {
        if (run->reads[i].indexed == run->reads[read_i].indexed) ++n;
    }
    return n == 0 ? 1 : n;
}

static gzFile *bcl_output_slot(gzFile *files, size_t sample_i, size_t read_i, size_t read_count) {
    return &files[sample_i * read_count + read_i];
}

static int append_fastq_record(text_buffer *out, const char *header, const char *seq, const char *qual) {
    char buf[20000];
    int n = snprintf(buf, sizeof(buf), "%s\n%s\n+\n%s\n", header, seq, qual);
    if (n < 0 || (size_t)n >= sizeof(buf)) return -1;
    return text_buffer_append(out, buf, (size_t)n);
}

typedef struct bcl_block_result {
    text_buffer *sample_buffers;
    text_buffer *undetermined_buffers;
    unsigned long long *sample_assigned;
    unsigned long long passed_clusters;
    unsigned long long filtered_clusters;
    unsigned long long undetermined_reads;
    bcl_unknown_table unknowns;
    int error;
} bcl_block_result;

typedef struct bcl_block_job {
    const bcl_run_info *run;
    const bcl_sample_table *samples;
    unsigned char **cycles;
    const unsigned char *pf;
    const char *tile;
    size_t start;
    size_t end;
    int k1;
    int k2;
    int emit_index_fastqs;
    bcl_block_result *result;
} bcl_block_job;

static void free_bcl_block_result(bcl_block_result *result, size_t sample_count, size_t read_count) {
    if (result->sample_buffers != NULL) {
        size_t n = 0;
        if (checked_mul_size(sample_count, read_count, &n) != 0) n = 0;
        for (size_t i = 0; i < n; ++i) free_text_buffer(&result->sample_buffers[i]);
    }
    if (result->undetermined_buffers != NULL) {
        for (size_t i = 0; i < read_count; ++i) free_text_buffer(&result->undetermined_buffers[i]);
    }
    free(result->sample_buffers);
    free(result->undetermined_buffers);
    free(result->sample_assigned);
    free_bcl_unknowns(&result->unknowns);
    memset(result, 0, sizeof(*result));
}

static int init_bcl_block_result(bcl_block_result *result, size_t sample_count, size_t read_count) {
    memset(result, 0, sizeof(*result));
    size_t sample_read_slots = 0;
    if (checked_mul_size(sample_count, read_count, &sample_read_slots) != 0) return -1;
    result->sample_buffers = (text_buffer *)calloc(alloc_count_or_one(sample_read_slots), sizeof(text_buffer));
    result->undetermined_buffers = (text_buffer *)calloc(alloc_count_or_one(read_count), sizeof(text_buffer));
    result->sample_assigned = (unsigned long long *)calloc(alloc_count_or_one(sample_count), sizeof(unsigned long long));
    if (result->sample_buffers == NULL || result->undetermined_buffers == NULL || result->sample_assigned == NULL) {
        free_bcl_block_result(result, sample_count, read_count);
        return -1;
    }
    return 0;
}

static int process_bcl_block(const bcl_block_job *job) {
    const bcl_run_info *run = job->run;
    const bcl_sample_table *samples = job->samples;
    bcl_block_result *result = job->result;
    for (size_t cluster = job->start; cluster < job->end; ++cluster) {
        if (job->pf != NULL && job->pf[cluster] == 0) {
            ++result->filtered_clusters;
            continue;
        }
        ++result->passed_clusters;
        char index1[MAX_BCL_READ_CYCLES + 1] = "";
        char index2[MAX_BCL_READ_CYCLES + 1] = "";
        char seqs[16][MAX_BCL_READ_CYCLES + 1];
        char quals[16][MAX_BCL_READ_CYCLES + 1];
        memset(seqs, 0, sizeof(seqs));
        memset(quals, 0, sizeof(quals));
        int indexed_seen = 0;
        for (size_t r = 0; r < run->read_count; ++r) {
            const bcl_read_info *ri = &run->reads[r];
            char *seq_out = seqs[r];
            char *qual_out = quals[r];
            size_t cap = sizeof(seqs[r]);
            for (size_t j = 0; j < ri->cycles && j + 1 < cap; ++j) {
                unsigned char b = job->cycles[ri->start_cycle - 1 + j][cluster];
                seq_out[j] = bcl_base(b);
                qual_out[j] = bcl_qual(b);
            }
            seq_out[ri->cycles] = '\0';
            qual_out[ri->cycles] = '\0';
            if (ri->indexed) {
                if (indexed_seen == 0) {
                    int n = snprintf(index1, sizeof(index1), "%s", seq_out);
                    if (n < 0 || (size_t)n >= sizeof(index1)) return -1;
                } else if (indexed_seen == 1) {
                    int n = snprintf(index2, sizeof(index2), "%s", seq_out);
                    if (n < 0 || (size_t)n >= sizeof(index2)) return -1;
                }
                ++indexed_seen;
            }
        }
        int match_count = 0;
        int sample_index = assign_bcl_sample(samples, 1, index1, index2, job->k1, job->k2, &match_count);
        if (sample_index >= 0) {
            const bcl_sample *s = &samples->items[sample_index];
            size_t out_i = s->output_index;
            for (size_t r = 0; r < run->read_count; ++r) {
                const bcl_read_info *ri = &run->reads[r];
                if (!bcl_output_enabled(ri, job->emit_index_fastqs)) continue;
                char header[4096];
                int n = snprintf(header, sizeof(header), "@DOTMATCH:1:%s:%zu %d:N:0:%s%s%s",
                                 job->tile, cluster + 1, bcl_output_number(run, r),
                                 index1, index2[0] ? "+" : "", index2);
                if (n < 0 || (size_t)n >= sizeof(header)) return -1;
                text_buffer *buf = &result->sample_buffers[out_i * run->read_count + r];
                if (append_fastq_record(buf, header, seqs[r], quals[r]) != 0) return -1;
            }
            ++result->sample_assigned[out_i];
        } else {
            char unknown_index[(MAX_BCL_READ_CYCLES * 2) + 2];
            int n = snprintf(unknown_index, sizeof(unknown_index), "%s%s%s", index1, index2[0] ? "+" : "", index2);
            if (n < 0 || (size_t)n >= sizeof(unknown_index)) return -1;
            if (add_bcl_unknown(&result->unknowns, unknown_index) != 0) return -1;
            for (size_t r = 0; r < run->read_count; ++r) {
                const bcl_read_info *ri = &run->reads[r];
                if (!bcl_output_enabled(ri, job->emit_index_fastqs)) continue;
                char header[4096];
                n = snprintf(header, sizeof(header), "@DOTMATCH:1:%s:%zu %d:N:0:%s%s%s",
                             job->tile, cluster + 1, bcl_output_number(run, r),
                             index1, index2[0] ? "+" : "", index2);
                if (n < 0 || (size_t)n >= sizeof(header)) return -1;
                if (append_fastq_record(&result->undetermined_buffers[r], header, seqs[r], quals[r]) != 0) return -1;
            }
            ++result->undetermined_reads;
        }
    }
    return 0;
}

static void *bcl_block_worker(void *arg) {
    bcl_block_job *job = (bcl_block_job *)arg;
    job->result->error = process_bcl_block(job);
    return NULL;
}

static int write_bcl_block_result(const bcl_block_result *result, gzFile *sample_fastqs, gzFile *undetermined_fastqs,
                                  size_t sample_count, size_t read_count) {
    for (size_t i = 0; i < sample_count; ++i) {
        for (size_t r = 0; r < read_count; ++r) {
            const text_buffer *buf = &result->sample_buffers[i * read_count + r];
            if (buf->len == 0) continue;
            gzFile gz = *bcl_output_slot(sample_fastqs, i, r, read_count);
            if (gz != NULL && gzwrite_all(gz, buf->data, buf->len) != 0) return -1;
        }
    }
    for (size_t r = 0; r < read_count; ++r) {
        const text_buffer *buf = &result->undetermined_buffers[r];
        if (buf->len == 0) continue;
        if (undetermined_fastqs[r] != NULL && gzwrite_all(undetermined_fastqs[r], buf->data, buf->len) != 0) return -1;
    }
    return 0;
}

static int parse_mismatches(const char *s, int *k1, int *k2) {
    char *comma = strchr(s, ',');
    if (comma == NULL) {
        int k = 0;
        if (parse_int_value(s, &k) != 0) return -1;
        if (k < 0 || k > 1) return -1;
        *k1 = k;
        *k2 = k;
        return 0;
    }
    char left[16];
    size_t n = (size_t)(comma - s);
    if (n >= sizeof(left)) return -1;
    memcpy(left, s, n);
    left[n] = '\0';
    int a = 0;
    int b = 0;
    if (parse_int_value(left, &a) != 0 || parse_int_value(comma + 1, &b) != 0) return -1;
    if (a < 0 || a > 1 || b < 0 || b > 1) return -1;
    *k1 = a;
    *k2 = b;
    return 0;
}

static int run_bcl_demux(const char *argv0, int argc, char **argv) {
    const char *run_folder = NULL;
    const char *sample_sheet = NULL;
    const char *out_dir = NULL;
    const char *summary_path = NULL;
    const char *mismatches = "1";
    const char *lanes = "1";
    int k1 = 1;
    int k2 = 1;
    int emit_index_fastqs = 0;
    size_t requested_threads = 0;
    int gzip_level = 1;

    int i = 2;
    while (i < argc) {
        const char *arg = argv[i++];
        if (strcmp(arg, "--run-folder") == 0 && i < argc) {
            run_folder = argv[i++];
        } else if (strcmp(arg, "--sample-sheet") == 0 && i < argc) {
            sample_sheet = argv[i++];
        } else if (strcmp(arg, "--out-dir") == 0 && i < argc) {
            out_dir = argv[i++];
        } else if (strcmp(arg, "--summary") == 0 && i < argc) {
            summary_path = argv[i++];
        } else if (strcmp(arg, "--barcode-mismatches") == 0 && i < argc) {
            mismatches = argv[i++];
        } else if (strcmp(arg, "--emit-index-fastqs") == 0) {
            emit_index_fastqs = 1;
        } else if (strcmp(arg, "--threads") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &requested_threads) != 0) {
                fprintf(stderr, "invalid --threads value\n");
                return 2;
            }
        } else if (strcmp(arg, "--gzip-level") == 0 && i < argc) {
            if (parse_int_value(argv[i++], &gzip_level) != 0 || gzip_level < 0 || gzip_level > 9) {
                fprintf(stderr, "invalid --gzip-level value\n");
                return 2;
            }
        } else if (strcmp(arg, "--lanes") == 0 && i < argc) {
            lanes = argv[i++];
        } else if (strcmp(arg, "--interop-dir") == 0 && i < argc) {
            i++;
        } else {
            usage(argv0);
            return 2;
        }
    }
    if (run_folder == NULL || sample_sheet == NULL || out_dir == NULL || parse_mismatches(mismatches, &k1, &k2) != 0) {
        usage(argv0);
        return 2;
    }
    if (strcmp(lanes, "1") != 0 && strcmp(lanes, "001") != 0) {
        fprintf(stderr, "classic BCL demux currently supports lane 1 only; rerun with --lanes 1 or split the run externally\n");
        return 2;
    }

    char basecalls[4096];
    int basecalls_n = snprintf(basecalls, sizeof(basecalls), "%s/Data/Intensities/BaseCalls", run_folder);
    if (basecalls_n < 0 || (size_t)basecalls_n >= sizeof(basecalls)) {
        fprintf(stderr, "run folder path is too long\n");
        return 2;
    }
    bcl_run_info run = {0};
    bcl_sample_table samples = {0};
    gzFile *sample_fastqs = NULL;
    gzFile *undetermined_fastqs = NULL;
    char **tiles = NULL;
    size_t tile_count = 0;
    unsigned long long total_clusters = 0;
    unsigned long long passed_clusters = 0;
    unsigned long long filtered_clusters = 0;
    unsigned long long undetermined_reads = 0;
    size_t effective_threads = 1;
    bcl_unknown_table unknowns = {0};
    int rc = 1;

    if (parse_run_info(run_folder, &run) != 0) {
        fprintf(stderr, "failed to parse RunInfo.xml\n");
        goto done;
    }
    if (read_bcl_sample_sheet(sample_sheet, &samples) != 0) {
        fprintf(stderr, "failed to parse sample sheet\n");
        goto done;
    }
    if (samples.count > MAX_BCL_SAMPLE_ROWS) {
        fprintf(stderr, "sample sheet exceeds supported BCL sample count\n");
        goto done;
    }
    if (requested_threads == 0) {
        requested_threads = get_cpu_count();
    }
    if (ensure_dir(out_dir) != 0) {
        fprintf(stderr, "failed to create BCL output directory\n");
        goto done;
    }
    if (collect_tiles(basecalls, 1, &tiles, &tile_count) != 0) {
        fprintf(stderr, "failed to find classic BCL tiles; CBCL is not supported in this milestone\n");
        goto done;
    }

    size_t output_file_count = 0;
    if (checked_mul_size(samples.count == 0 ? 1 : samples.count,
                         run.read_count == 0 ? 1 : run.read_count,
                         &output_file_count) != 0) {
        fprintf(stderr, "BCL output file count overflow\n");
        goto done;
    }
    sample_fastqs = (gzFile *)calloc(output_file_count, sizeof(gzFile));
    undetermined_fastqs = (gzFile *)calloc(run.read_count == 0 ? 1 : run.read_count, sizeof(gzFile));
    if (sample_fastqs == NULL || undetermined_fastqs == NULL) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    for (size_t i = 0; i < samples.count; ++i) {
        if (samples.items[i].is_alias) continue;
        for (size_t r = 0; r < run.read_count; ++r) {
            bcl_read_info *ri = &run.reads[r];
            if (!bcl_output_enabled(ri, emit_index_fastqs)) continue;
            gzFile *slot = bcl_output_slot(sample_fastqs, i, r, run.read_count);
            *slot = open_bcl_fastq(out_dir, samples.items[i].id, i + 1, 1, bcl_output_kind(ri),
                                   bcl_output_number(&run, r), gzip_level);
            if (*slot == NULL) {
                fprintf(stderr, "failed to open sample FASTQ\n");
                goto done;
            }
        }
    }
    for (size_t r = 0; r < run.read_count; ++r) {
        bcl_read_info *ri = &run.reads[r];
        if (!bcl_output_enabled(ri, emit_index_fastqs)) continue;
        undetermined_fastqs[r] = open_bcl_fastq(out_dir, "Undetermined", 0, 1, bcl_output_kind(ri),
                                                bcl_output_number(&run, r), gzip_level);
        if (undetermined_fastqs[r] == NULL) {
            fprintf(stderr, "failed to open undetermined FASTQ\n");
            goto done;
        }
    }

    for (size_t tile_i = 0; tile_i < tile_count; ++tile_i) {
        unsigned char **cycles = (unsigned char **)calloc(run.total_cycles == 0 ? 1 : run.total_cycles, sizeof(unsigned char *));
        size_t cluster_count = 0;
        if (cycles == NULL) {
            fprintf(stderr, "out of memory\n");
            goto done;
        }
        for (size_t c = 1; c <= run.total_cycles; ++c) {
            char bcl_path[4096];
            size_t n = 0;
            if (build_bcl_path(bcl_path, sizeof(bcl_path), basecalls, 1, c, tiles[tile_i]) != 0 ||
                read_bcl_cycle(bcl_path, &cycles[c - 1], &n) != 0) {
                fprintf(stderr, "failed to read BCL cycle\n");
                for (size_t j = 0; j < run.total_cycles; ++j) free(cycles[j]);
                free(cycles);
                goto done;
            }
            if (c == 1) cluster_count = n;
            else if (n != cluster_count) {
                fprintf(stderr, "BCL cycle cluster counts do not match\n");
                for (size_t j = 0; j < run.total_cycles; ++j) free(cycles[j]);
                free(cycles);
                goto done;
            }
        }
        unsigned char *pf = NULL;
        size_t pf_count = 0;
        if (read_filter_file(basecalls, 1, tiles[tile_i], &pf, &pf_count) != 0) {
            fprintf(stderr, "failed to read filter file\n");
            for (size_t j = 0; j < run.total_cycles; ++j) free(cycles[j]);
            free(cycles);
            goto done;
        }
        if (pf != NULL && pf_count != cluster_count) {
            fprintf(stderr, "filter cluster count does not match BCL\n");
            free(pf);
            for (size_t j = 0; j < run.total_cycles; ++j) free(cycles[j]);
            free(cycles);
            goto done;
        }

        total_clusters += cluster_count;
        size_t threads = requested_threads;
        if (threads > cluster_count) threads = cluster_count == 0 ? 1 : cluster_count;
        if (threads > effective_threads) effective_threads = threads;
        const size_t block_size = 8192;
        for (size_t block_start = 0; block_start < cluster_count;) {
            size_t batch = 0;
            pthread_t *thread_ids = NULL;
            bcl_block_job *jobs = (bcl_block_job *)calloc(threads, sizeof(bcl_block_job));
            bcl_block_result *results = (bcl_block_result *)calloc(threads, sizeof(bcl_block_result));
            if (jobs == NULL || results == NULL) {
                free(jobs);
                free(results);
                fprintf(stderr, "out of memory\n");
                goto done;
            }
            if (threads > 1) {
                thread_ids = (pthread_t *)calloc(threads, sizeof(pthread_t));
                if (thread_ids == NULL) {
                    free(jobs);
                    free(results);
                    fprintf(stderr, "out of memory\n");
                    goto done;
                }
            }
            while (batch < threads && block_start < cluster_count) {
                size_t block_end = block_start + block_size;
                if (block_end > cluster_count) block_end = cluster_count;
                if (init_bcl_block_result(&results[batch], samples.count, run.read_count) != 0) {
                    fprintf(stderr, "out of memory\n");
                    goto done;
                }
                jobs[batch].run = &run;
                jobs[batch].samples = &samples;
                jobs[batch].cycles = cycles;
                jobs[batch].pf = pf;
                jobs[batch].tile = tiles[tile_i];
                jobs[batch].start = block_start;
                jobs[batch].end = block_end;
                jobs[batch].k1 = k1;
                jobs[batch].k2 = k2;
                jobs[batch].emit_index_fastqs = emit_index_fastqs;
                jobs[batch].result = &results[batch];
                if (threads > 1) {
                    if (pthread_create(&thread_ids[batch], NULL, bcl_block_worker, &jobs[batch]) != 0) {
                        fprintf(stderr, "failed to create BCL worker\n");
                        goto done;
                    }
                } else {
                    results[batch].error = process_bcl_block(&jobs[batch]);
                }
                ++batch;
                block_start = block_end;
            }
            if (threads > 1) {
                for (size_t i = 0; i < batch; ++i) pthread_join(thread_ids[i], NULL);
            }
            for (size_t i = 0; i < batch; ++i) {
                if (results[i].error != 0) {
                    fprintf(stderr, "failed to format BCL block\n");
                    goto done;
                }
                passed_clusters += results[i].passed_clusters;
                filtered_clusters += results[i].filtered_clusters;
                undetermined_reads += results[i].undetermined_reads;
                for (size_t s = 0; s < samples.count; ++s) samples.items[s].assigned += results[i].sample_assigned[s];
                if (merge_bcl_unknowns(&unknowns, &results[i].unknowns) != 0) {
                    fprintf(stderr, "out of memory\n");
                    goto done;
                }
                if (write_bcl_block_result(&results[i], sample_fastqs, undetermined_fastqs, samples.count, run.read_count) != 0) {
                    fprintf(stderr, "failed to write BCL block\n");
                    goto done;
                }
                free_bcl_block_result(&results[i], samples.count, run.read_count);
            }
            free(thread_ids);
            free(jobs);
            free(results);
        }
        free(pf);
        for (size_t j = 0; j < run.total_cycles; ++j) free(cycles[j]);
        free(cycles);
    }

    char stats_path[4096];
    if (path_join(stats_path, sizeof(stats_path), out_dir, "Demultiplex_Stats.csv") != 0) goto done;
    FILE *stats = open_output_file(stats_path);
    if (stats == NULL) goto done;
    fprintf(stats, "sample_id,assigned_reads");
    int non_index_read_count = 0;
    for (size_t r = 0; r < run.read_count; ++r) {
        if (!run.reads[r].indexed) fprintf(stats, ",read%d_records", ++non_index_read_count);
    }
    fprintf(stats, "\n");
    unsigned long long assigned_reads = 0;
    for (size_t i = 0; i < samples.count; ++i) {
        if (samples.items[i].is_alias) continue;
        fprintf(stats, "%s,%llu", samples.items[i].id, samples.items[i].assigned);
        for (int r = 0; r < non_index_read_count; ++r) fprintf(stats, ",%llu", samples.items[i].assigned);
        fprintf(stats, "\n");
        assigned_reads += samples.items[i].assigned;
    }
    fprintf(stats, "Undetermined,%llu", undetermined_reads);
    for (int r = 0; r < non_index_read_count; ++r) fprintf(stats, ",%llu", undetermined_reads);
    fprintf(stats, "\n");
    fclose(stats);

    if (unknowns.count > 0) {
        qsort(unknowns.items, unknowns.count, sizeof(unknowns.items[0]), cmp_bcl_unknown_desc);
        char unknown_path[4096];
        if (path_join(unknown_path, sizeof(unknown_path), out_dir, "Top_Unknown_Barcodes.csv") != 0) goto done;
        FILE *unknown = open_output_file(unknown_path);
        if (unknown == NULL) goto done;
        fprintf(unknown, "index,count\n");
        size_t n = unknowns.count < 100 ? unknowns.count : 100;
        for (size_t i = 0; i < n; ++i) fprintf(unknown, "%s,%llu\n", unknowns.items[i].index, unknowns.items[i].count);
        fclose(unknown);
    }

    char normalized_path[4096];
    if (path_join(normalized_path, sizeof(normalized_path), out_dir, "SampleSheet.normalized.csv") != 0) goto done;
    FILE *normalized = open_output_file(normalized_path);
    if (normalized != NULL) {
        fprintf(normalized, "sample_id,sample_name,lane,index,index2\n");
        for (size_t i = 0; i < samples.count; ++i) {
            fprintf(normalized, "%s,%s,%d,%s,%s\n", samples.items[i].id, samples.items[i].name,
                    samples.items[i].lane, samples.items[i].index1, samples.items[i].index2);
        }
        fclose(normalized);
    }

    if (summary_path != NULL) {
        FILE *summary = open_output_file(summary_path);
        if (summary == NULL) goto done;
        fprintf(summary,
                "{\n  \"workflow\": \"bcl-demux\",\n  \"format\": \"classic_bcl\",\n  \"lanes\": 1,\n  \"tiles\": %zu,\n  \"total_clusters\": %llu,\n  \"passed_filter_clusters\": %llu,\n  \"filtered_clusters\": %llu,\n  \"assigned_reads\": %llu,\n  \"undetermined_reads\": %llu,\n  \"barcode_mismatches_index1\": %d,\n  \"barcode_mismatches_index2\": %d,\n  \"requested_threads\": %zu,\n  \"effective_threads\": %zu,\n  \"gzip_level\": %d,\n  \"emit_index_fastqs\": %s\n}\n",
                tile_count, total_clusters, passed_clusters, filtered_clusters, assigned_reads, undetermined_reads,
                k1, k2, requested_threads, effective_threads, gzip_level, emit_index_fastqs ? "true" : "false");
        fclose(summary);
    }

    rc = 0;

done:
    if (sample_fastqs != NULL) {
        for (size_t i = 0; i < samples.count; ++i) {
            if (samples.items[i].is_alias) continue;
            for (size_t r = 0; r < run.read_count; ++r) {
                gzFile *slot = bcl_output_slot(sample_fastqs, i, r, run.read_count);
                if (*slot != NULL) gzclose(*slot);
            }
        }
    }
    if (undetermined_fastqs != NULL) {
        for (size_t r = 0; r < run.read_count; ++r) {
            if (undetermined_fastqs[r] != NULL) gzclose(undetermined_fastqs[r]);
        }
    }
    if (tiles != NULL) {
        for (size_t i = 0; i < tile_count; ++i) free(tiles[i]);
        free(tiles);
    }
    free(sample_fastqs);
    free(undetermined_fastqs);
    free_bcl_unknowns(&unknowns);
    free_bcl_samples(&samples);
    return rc;
}

static int compare_gzip_fastq_files(const char *a_path, const char *b_path, unsigned long long *records_out) {
    gzFile a = gzopen(a_path, "rb");
    gzFile b = gzopen(b_path, "rb");
    if (a == NULL || b == NULL) {
        if (a != NULL) gzclose(a);
        if (b != NULL) gzclose(b);
        return -1;
    }
    char abuf[8192];
    char bbuf[8192];
    unsigned long long lines = 0;
    int mismatch = 0;
    for (;;) {
        char *ag = gzgets(a, abuf, sizeof(abuf));
        char *bg = gzgets(b, bbuf, sizeof(bbuf));
        if (ag == NULL || bg == NULL) {
            if (ag != bg) mismatch = 1;
            break;
        }
        if (strcmp(abuf, bbuf) != 0) mismatch = 1;
        ++lines;
    }
    gzclose(a);
    gzclose(b);
    *records_out = lines / 4;
    return mismatch ? 1 : 0;
}

static int run_bcl_validate(const char *argv0, int argc, char **argv) {
    const char *dotmatch_out = NULL;
    const char *truth_out = NULL;
    int i = 2;
    while (i < argc) {
        const char *arg = argv[i++];
        if (strcmp(arg, "--dotmatch-out") == 0 && i < argc) {
            dotmatch_out = argv[i++];
        } else if (strcmp(arg, "--truth-out") == 0 && i < argc) {
            truth_out = argv[i++];
        } else {
            usage(argv0);
            return 2;
        }
    }
    if (dotmatch_out == NULL || truth_out == NULL) {
        usage(argv0);
        return 2;
    }
    DIR *dir = opendir(truth_out);
    if (dir == NULL) {
        fprintf(stderr, "failed to open truth output directory\n");
        return 1;
    }
    unsigned long long compared_files = 0;
    unsigned long long compared_records = 0;
    unsigned long long missing_files = 0;
    unsigned long long mismatched_files = 0;
    struct dirent *ent;
    while ((ent = readdir(dir)) != NULL) {
        if (!ends_with(ent->d_name, ".fastq.gz")) continue;
        char truth_path[4096];
        char dotmatch_path[4096];
        if (path_join(truth_path, sizeof(truth_path), truth_out, ent->d_name) != 0 ||
            path_join(dotmatch_path, sizeof(dotmatch_path), dotmatch_out, ent->d_name) != 0) {
            closedir(dir);
            fprintf(stderr, "BCL validation path is too long\n");
            return 1;
        }
        if (!path_exists(dotmatch_path)) {
            ++missing_files;
            continue;
        }
        unsigned long long records = 0;
        int cmp = compare_gzip_fastq_files(dotmatch_path, truth_path, &records);
        if (cmp != 0) ++mismatched_files;
        compared_records += records;
        ++compared_files;
    }
    closedir(dir);
    printf("{\n  \"compared_fastq_files\": %llu,\n  \"compared_records\": %llu,\n  \"missing_fastq_files\": %llu,\n  \"mismatched_fastq_files\": %llu\n}\n",
           compared_files, compared_records, missing_files, mismatched_files);
    return missing_files == 0 && mismatched_files == 0 ? 0 : 1;
}

static int run_edlib_validate_helper(const char *targets_path, const char *reads_path,
        size_t target_start, size_t target_len, int k, size_t indel_window, size_t sample_limit,
        size_t auto_offset, size_t auto_offset_sample, offset_mode offsets_mode, double offset_min_fraction,
        size_t threads) {
    const char *helper_path = "./build/dotmatch_edlib_validate";
    if (access(helper_path, X_OK) != 0) {
        fprintf(stderr, "edlib oracle validation requires build/dotmatch_edlib_validate; run `make edlib-tools`\n");
        return 2;
    }

    char target_start_buf[32];
    char target_len_buf[32];
    char k_buf[32];
    char indel_window_buf[32];
    char sample_buf[32];
    char auto_offset_buf[32];
    char auto_offset_sample_buf[32];
    char threads_buf[32];
    char offset_min_fraction_buf[64];
    snprintf(target_start_buf, sizeof(target_start_buf), "%zu", target_start);
    snprintf(target_len_buf, sizeof(target_len_buf), "%zu", target_len);
    snprintf(k_buf, sizeof(k_buf), "%d", k);
    snprintf(indel_window_buf, sizeof(indel_window_buf), "%zu", indel_window);
    snprintf(sample_buf, sizeof(sample_buf), "%zu", sample_limit);
    snprintf(auto_offset_buf, sizeof(auto_offset_buf), "%zu", auto_offset);
    snprintf(auto_offset_sample_buf, sizeof(auto_offset_sample_buf), "%zu", auto_offset_sample);
    snprintf(threads_buf, sizeof(threads_buf), "%zu", threads);
    snprintf(offset_min_fraction_buf, sizeof(offset_min_fraction_buf), "%.8f", offset_min_fraction);

    pid_t pid = fork();
    if (pid < 0) {
        perror("fork");
        return 1;
    }
    if (pid == 0) {
        execl(helper_path, helper_path,
              "--targets", targets_path,
              "--reads", reads_path,
              "--target-start", target_start_buf,
              "--target-length", target_len_buf,
              "--k", k_buf,
              "--indel-window", indel_window_buf,
              "--auto-offset", auto_offset_buf,
              "--auto-offset-sample", auto_offset_sample_buf,
              "--offset-mode", offset_mode_name(offsets_mode),
              "--offset-min-fraction", offset_min_fraction_buf,
              "--sample", sample_buf,
              "--threads", threads_buf,
              (char *)NULL);
        perror("execl");
        _exit(127);
    }

    int status = 0;
    if (waitpid(pid, &status, 0) < 0) {
        perror("waitpid");
        return 1;
    }
    if (WIFEXITED(status)) {
        return WEXITSTATUS(status);
    }
    return 1;
}

static int run_validate(const char *argv0, int argc, char **argv) {
    const char *targets_path = NULL;
    const char *reads_path = NULL;
    const char *oracle = "scan";
    size_t target_start = 0;
    size_t target_len = 0;
    size_t indel_window = 0;
    size_t sample_limit = 100000;
    size_t auto_offset = 0;
    size_t auto_offset_sample = 1000;
    size_t threads = 0;
    offset_mode offsets_mode = OFFSET_MODE_BEST;
    double offset_min_fraction = 0.005;
    count_metric metric = COUNT_METRIC_LEVENSHTEIN;
    int k = -1;

    int i = 2;
    while (i < argc) {
        const char *arg = argv[i++];
        if (strcmp(arg, "--targets") == 0 && i < argc) {
            targets_path = argv[i++];
        } else if (strcmp(arg, "--reads") == 0 && i < argc) {
            reads_path = argv[i++];
        } else if (strcmp(arg, "--target-start") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &target_start) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--target-length") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &target_len) != 0 || target_len == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--k") == 0 && i < argc) {
            if (parse_int_value(argv[i++], &k) != 0 || (k != 0 && k != 1)) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--indel-window") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &indel_window) != 0 || indel_window > 1) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--metric") == 0 && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "hamming") == 0) {
                metric = COUNT_METRIC_HAMMING;
            } else if (strcmp(value, "levenshtein") == 0) {
                metric = COUNT_METRIC_LEVENSHTEIN;
            } else {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--auto-offset") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &auto_offset) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--auto-offset-sample") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &auto_offset_sample) != 0 || auto_offset_sample == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--offset-mode") == 0 && i < argc) {
            const char *value = argv[i++];
            if (strcmp(value, "best") == 0) {
                offsets_mode = OFFSET_MODE_BEST;
            } else if (strcmp(value, "multi") == 0) {
                offsets_mode = OFFSET_MODE_MULTI;
            } else {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--offset-min-fraction") == 0 && i < argc) {
            if (parse_double_value(argv[i++], &offset_min_fraction) != 0 ||
                offset_min_fraction < 0.0 || offset_min_fraction > 1.0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--oracle") == 0 && i < argc) {
            oracle = argv[i++];
        } else if (strcmp(arg, "--sample") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &sample_limit) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(arg, "--threads") == 0 && i < argc) {
            if (parse_size_value(argv[i++], &threads) != 0) {
                usage(argv0);
                return 2;
            }
        } else {
            usage(argv0);
            return 2;
        }
    }

    if (targets_path == NULL || reads_path == NULL || target_len == 0 || k < 0) {
        usage(argv0);
        return 2;
    }
    if (auto_offset > MAX_AUTO_OFFSET) {
        fprintf(stderr, "--auto-offset must be <= %d\n", MAX_AUTO_OFFSET);
        return 2;
    }
    if (metric == COUNT_METRIC_HAMMING && indel_window != 0) {
        fprintf(stderr, "--indel-window is only valid with --metric levenshtein\n");
        return 2;
    }
    if (threads == 0) {
        threads = get_cpu_count();
    }
    if (strcmp(oracle, "edlib") == 0) {
        if (metric != COUNT_METRIC_LEVENSHTEIN) {
            fprintf(stderr, "--oracle edlib is only valid with --metric levenshtein\n");
            return 2;
        }
        int status = run_edlib_validate_helper(targets_path, reads_path, target_start, target_len, k, indel_window,
                                               sample_limit, auto_offset, auto_offset_sample, offsets_mode,
                                               offset_min_fraction, threads);
        return status == 0 ? 0 : status;
    }
    if (strcmp(oracle, "scan") != 0) {
        usage(argv0);
        return 2;
    }

    seq_table targets = {0};
    fastq_reader reader = {0};
    qdaln_index *index = NULL;
    const char **target_ptrs = NULL;
    size_t *target_lens = NULL;
    int rc = 1;
    size_t checked = 0;
    size_t mismatches = 0;
    offset_list offsets = {0};

    if (read_target_table(targets_path, &targets) != 0) {
        fprintf(stderr, "failed to read targets\n");
        goto done;
    }
    if (build_target_arrays(&targets, &target_ptrs, &target_lens) != 0) {
        fprintf(stderr, "out of memory\n");
        goto done;
    }
    index = qdaln_index_build(target_ptrs, target_lens, targets.count);
    if (index == NULL) {
        fprintf(stderr, "failed to build target index\n");
        goto done;
    }
    if (fastq_reader_open(&reader, reads_path) != 0) {
        fprintf(stderr, "failed to open FASTQ input\n");
        goto done;
    }
    if (detect_offsets(index, NULL, reads_path, target_start, target_len, auto_offset, auto_offset_sample,
                       offsets_mode, offset_min_fraction, &offsets) != 0) {
        fprintf(stderr, "automatic offset detection failed\n");
        goto done;
    }

    char header[8192];
    char seq[8192];
    char plus[8192];
    char qual[8192];
    char observed[8192];
    char scan_observed[8192];
    int got = 0;
    size_t seq_len = 0;
    while ((sample_limit == 0 || checked < sample_limit) &&
           (got = fastq_read_record_len(&reader, header, seq, plus, qual, sizeof(header), &seq_len)) == 1) {
        qdaln_match_result indexed;
        qdaln_match_result scan;
        qdaln_index_stats stats;
        if (assign_count_offsets(index, seq, seq_len, &offsets, target_start, target_len, k, metric,
                                 indel_window, &indexed, &stats, observed, sizeof(observed), 0) != 0 ||
            scan_count_offsets(target_ptrs, target_lens, targets.count, seq, seq_len, &offsets, target_start,
                               target_len, k, metric, indel_window, &scan, scan_observed,
                               sizeof(scan_observed)) != 0) {
            fprintf(stderr, "validation assignment failed\n");
            goto done;
        }
        if (indexed.target_index != scan.target_index ||
            indexed.best_distance != scan.best_distance ||
            indexed.second_best_distance != scan.second_best_distance ||
            indexed.match_count != scan.match_count ||
            indexed.status != scan.status) {
            ++mismatches;
        }
        ++checked;
    }
    if (got < 0) {
        fprintf(stderr, "malformed FASTQ input\n");
        goto done;
    }
    printf("{\n  \"oracle\": \"native_scan\",\n  \"checked_reads\": %zu,\n  \"mismatches\": %zu,\n  \"k\": %d,\n  \"metric\": \"%s\",\n  \"target_start\": %zu,\n  \"target_length\": %zu,\n  \"offset_mode\": \"%s\",\n  \"selected_target_starts\": [",
           checked, mismatches, k, metric_name(metric), target_start, target_len, offset_mode_name(offsets_mode));
    for (size_t i = 0; i < offsets.count; ++i) {
        if (i != 0) printf(", ");
        printf("%zu", offsets.items[i]);
    }
    printf("]\n}\n");
    rc = mismatches == 0 ? 0 : 1;

done:
    fastq_reader_close(&reader);
    qdaln_index_free(index);
    free_offset_list(&offsets);
    free(target_ptrs);
    free(target_lens);
    free_table(&targets);
    return rc;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 2;
    }

    if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "help") == 0) {
        if (argc != 2) {
            usage(argv[0]);
            return 2;
        }
        help_manual(stdout, argv[0]);
        return 0;
    }

    if (strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "version") == 0) {
        if (argc != 2) {
            usage(argv[0]);
            return 2;
        }
        printf("dotmatch %s\n", DOTMATCH_VERSION);
        return 0;
    }

    if (strcmp(argv[1], "citation") == 0 || strcmp(argv[1], "cite") == 0) {
        if (argc != 2) {
            usage(argv[0]);
            return 2;
        }
        print_citation(stdout);
        return 0;
    }

    if (strcmp(argv[1], "dist") == 0) {
        if (argc != 4) {
            usage(argv[0]);
            return 2;
        }
        int d = qdaln_edit_distance(argv[2], strlen(argv[2]), argv[3], strlen(argv[3]));
        if (d < 0) return 1;
        printf("%d\n", d);
        return 0;
    }

    if (strcmp(argv[1], "leq") == 0) {
        if (argc != 5) {
            usage(argv[0]);
            return 2;
        }
        int k = 0;
        if (parse_int_value(argv[2], &k) != 0 || k < 0) {
            usage(argv[0]);
            return 2;
        }
        int ok = qdaln_edit_distance_leq(argv[3], strlen(argv[3]), argv[4], strlen(argv[4]), k);
        if (ok < 0) return 1;
        printf("%s\n", ok ? "true" : "false");
        return 0;
    }

    if (strcmp(argv[1], "assign") == 0 || strcmp(argv[1], "match") == 0) {
        return run_batch(argv[0], argc, argv, argv[1]);
    }

    if (strcmp(argv[1], "fastq-assign") == 0) {
        return run_fastq_assign(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "pair-count") == 0) {
        return run_pair_count(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "demux") == 0) {
        return run_demux(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "bcl-demux") == 0) {
        return run_bcl_demux(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "bcl-validate") == 0) {
        return run_bcl_validate(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "count") == 0 || strcmp(argv[1], "crispr-count") == 0) {
        if (help_requested(argc, argv)) {
            count_help_manual(stdout, argv[0], strcmp(argv[1], "crispr-count") == 0);
            return 0;
        }
        return run_count(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "guide-counter") == 0 || strcmp(argv[1], "guide-counter-count") == 0 ||
        strcmp(argv[1], "guide-count") == 0) {
        return run_guide_counter_compatible(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "inspect-unmatched") == 0) {
        return run_inspect_unmatched(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "audit") == 0 || strcmp(argv[1], "audit-targets") == 0) {
        if (help_requested(argc, argv)) {
            audit_help_manual(stdout, argv[0]);
            return 0;
        }
        return run_audit(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "validate") == 0) {
        return run_validate(argv[0], argc, argv);
    }

    usage(argv[0]);
    return 2;
}
