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
#define DOTMATCH_VERSION "0.1.9"
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
    fprintf(stderr, "  %s pair-count --left-targets left.tsv --right-targets right.tsv --reads reads.fastq[.gz] --left-start N --left-length L --right-start N --right-length L --k 0|1|2 --metric hamming|levenshtein [--ambiguity-policy radius|best] --out pair_counts.tsv [--summary summary.json]\n", argv0);
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
    fprintf(out, "  pair-count --left-targets left.tsv --right-targets right.tsv --reads reads.fastq[.gz] \\\n");
    fprintf(out, "      --left-start N --left-length L --right-start N --right-length L --out pair_counts.tsv\n");
    fprintf(out, "      Count pairs of independent fixed-window targets.\n");
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

static int read_target_table(const char *path, seq_table *table) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return -1;

    char buf[16384];
    int id_col = 0;
    int seq_col = 1;
    int gene_col = 2;
    int have_header = 0;
    int first_data = 1;
    size_t row = 0;
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        trim_line(buf);
        if (buf[0] == '\0' || buf[0] == '#') continue;

        char delim = strchr(buf, ',') != NULL && strchr(buf, '\t') == NULL ? ',' : '\t';
        char *fields[16];
        size_t nf = split_fields(buf, delim, fields, 16);
        if (first_data) {
            int maybe_id = find_column(fields, nf, "id", "target_id", "barcode_id");
            if (maybe_id < 0) maybe_id = find_column(fields, nf, "guide", "sgRNA", "sgrna");
            if (maybe_id < 0) maybe_id = find_column(fields, nf, "sgRNAID", "sgrnaid", "guide_id");
            if (maybe_id < 0) maybe_id = find_column(fields, nf, "sgRNA_ID", "sgrna_id", NULL);
            int maybe_seq = find_column(fields, nf, "gRNA.sequence", "target_seq", "sequence");
            if (maybe_seq < 0) maybe_seq = find_column(fields, nf, "bases", NULL, NULL);
            if (maybe_seq < 0) maybe_seq = find_column(fields, nf, "Seq", "seq", "barcode_seq");
            if (maybe_seq < 0) maybe_seq = find_column(fields, nf, "guide_seq", "sgRNA.sequence", "sgrna_sequence");
            if (maybe_seq < 0) maybe_seq = find_column(fields, nf, "sgRNA_seq", "guide_sequence", "GuideSequence");
            int maybe_gene = find_column(fields, nf, "Gene", "gene", "gene_symbol");
            if (maybe_gene < 0) maybe_gene = find_column(fields, nf, "gene.symbol", "target_gene", NULL);
            if (maybe_id >= 0 && maybe_seq >= 0) {
                id_col = maybe_id;
                seq_col = maybe_seq;
                gene_col = maybe_gene;
                have_header = 1;
                first_data = 0;
                continue;
            }
        }
        first_data = 0;

        const char *id = NULL;
        const char *seq = NULL;
        const char *gene = "";
        char id_buf[32];
        if (nf == 1) {
            int n = snprintf(id_buf, sizeof(id_buf), "%zu", row);
            if (n < 0 || (size_t)n >= sizeof(id_buf)) {
                fclose(fp);
                return -1;
            }
            id = id_buf;
            seq = fields[0];
        } else {
            if ((size_t)id_col >= nf || (size_t)seq_col >= nf) {
                fclose(fp);
                return -1;
            }
            id = fields[id_col];
            seq = fields[seq_col];
            if (have_header && gene_col >= 0 && (size_t)gene_col < nf) gene = fields[gene_col];
            if (!have_header && nf > 2) gene = fields[2];
        }
        if (id[0] == '\0' || seq[0] == '\0') {
            fprintf(stderr, "%s:%zu: target ID and sequence must be non-empty\n", path, row + 1);
            fclose(fp);
            return -1;
        }
        if (push_record_gene(table, id, strlen(id), seq, strlen(seq), gene, strlen(gene)) != 0) {
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
    return table->count == 0 ? -1 : 0;
}

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
        if (list->items[i] == offset) r