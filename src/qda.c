#include "qdalign.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>
#include <stdint.h>
#include <limits.h>
#include <fcntl.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <unistd.h>
#include <zlib.h>

#ifndef DOTMATCH_VERSION
#define DOTMATCH_VERSION "0.1.1"
#endif

#define MAX_AUTO_OFFSET 1024

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

static void usage(const char *argv0) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s --version\n", argv0);
    fprintf(stderr, "  %s dist SEQ1 SEQ2\n", argv0);
    fprintf(stderr, "  %s leq K SEQ1 SEQ2\n", argv0);
    fprintf(stderr, "  %s assign K barcodes.txt reads.txt\n", argv0);
    fprintf(stderr, "  %s match K targets.txt reads.txt\n", argv0);
    fprintf(stderr, "  %s fastq-assign --barcodes barcodes.tsv --reads reads.fastq[.gz] --barcode-start N --barcode-length L --k 0|1 --out assignments.tsv\n", argv0);
    fprintf(stderr, "  %s pair-count --left-targets left.tsv --right-targets right.tsv --reads reads.fastq[.gz] --left-start N --left-length L --right-start N --right-length L --k 0|1|2 --metric hamming|levenshtein --out pair_counts.tsv [--summary summary.json]\n", argv0);
    fprintf(stderr, "  %s demux --barcodes barcodes.tsv|barcodes.csv --reads reads.fastq[.gz] --barcode-start N --barcode-length L|auto --k 0|1|2 --metric hamming|levenshtein [--max-correction-qual Q] --out-dir demux_dir [--summary qc.json]\n", argv0);
    fprintf(stderr, "  %s bcl-demux --run-folder RUN --sample-sheet SampleSheet.csv --out-dir demux_dir --barcode-mismatches 0|1|1,1 [--threads N] [--gzip-level 0..9] [--emit-index-fastqs] [--summary summary.json]\n", argv0);
    fprintf(stderr, "  %s bcl-validate --dotmatch-out DIR --truth-out DIR\n", argv0);
    fprintf(stderr, "  %s count --targets targets.tsv|targets.csv --reads reads.fastq[.gz] [--reads more.fastq.gz] --sample-label labels --target-start N --target-length L --k 0|1|2 --metric hamming|levenshtein [--hamming-index auto|query|precompute] [--max-correction-qual Q] --ambiguity-policy best|radius --offset-mode best|multi --out counts.tsv [--format dotmatch|mageck]\n", argv0);
    fprintf(stderr, "  %s crispr-count --library guides.tsv|guides.csv --samples samples.tsv --guide-start N --guide-length L --k 0|1|2 --out counts.tsv [--summary qc.json]\n", argv0);
    fprintf(stderr, "  %s inspect-unmatched --targets targets.tsv|targets.csv --reads reads.fastq[.gz] --target-start N --target-length L --k 0|1 --top N --out top_unmatched.tsv [--low-quality-threshold Q]\n", argv0);
    fprintf(stderr, "  %s audit --targets targets.tsv|targets.csv --k 1 --out-dir audit_dir [--audit-mode auto|exact|fast]\n", argv0);
    fprintf(stderr, "  %s validate --targets targets.tsv|targets.csv --reads reads.fastq[.gz] --target-start N --target-length L --k 0|1 [--metric hamming|levenshtein] [--indel-window 0|1] [--offset-mode best|multi] [--threads N] --oracle scan|edlib\n", argv0);
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
    char *end = NULL;
    unsigned long v = strtoul(s, &end, 10);
    if (end == s || *end != '\0') return -1;
    *out = (size_t)v;
    return 0;
}

static int parse_int_value(const char *s, int *out) {
    char *end = NULL;
    long v = strtol(s, &end, 10);
    if (end == s || *end != '\0') return -1;
    *out = (int)v;
    return 0;
}

static int parse_double_value(const char *s, double *out) {
    char *end = NULL;
    double v = strtod(s, &end);
    if (end == s || *end != '\0') return -1;
    *out = v;
    return 0;
}

static int offset_count_for_range(size_t range, size_t *out) {
    if (range > MAX_AUTO_OFFSET || range > (SIZE_MAX - 1) / 2) return -1;
    *out = range * 2 + 1;
    return 0;
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
            int maybe_id = find_column(fields, nf, "id", "target_id", "sgRNA");
            if (maybe_id < 0) maybe_id = find_column(fields, nf, "sgRNAID", "sgrnaid", "guide_id");
            int maybe_seq = find_column(fields, nf, "gRNA.sequence", "target_seq", "sequence");
            if (maybe_seq < 0) maybe_seq = find_column(fields, nf, "Seq", "seq", "barcode_seq");
            if (maybe_seq < 0) maybe_seq = find_column(fields, nf, "guide_seq", "sgRNA.sequence", "sgrna_sequence");
            int maybe_gene = find_column(fields, nf, "Gene", "gene", NULL);
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
        if (seq[0] == '\0') {
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
    if (argc != 5) {
        usage(argv0);
        return 2;
    }

    int k = 0;
    if (sscanf(argv[2], "%d", &k) != 1 || k < 0) {
        usage(argv0);
        return 2;
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

static int parse_fastq_record(FILE *fp, char **id_out, char **seq_out, char **qual_out) {
    size_t cap = 0;
    char *line = NULL;
    ssize_t n = getline(&line, &cap, fp);
    if (n < 0) {
        free(line);
        return 0;
    }
    trim_line(line);
    if (line[0] != '@') {
        free(line);
        return -1;
    }
    char *id = xstrndup(line + 1, strlen(line + 1));
    free(line);
    line = NULL;
    cap = 0;
    n = getline(&line, &cap, fp);
    if (n < 0) {
        free(id);
        free(line);
        return -1;
    }
    trim_line(line);
    char *seq = xstrndup(line, strlen(line));
    uppercase_ascii(seq);
    free(line);
    line = NULL;
    cap = 0;
    n = getline(&line, &cap, fp);
    if (n < 0 || line[0] != '+') {
        free(id);
        free(seq);
        free(line);
        return -1;
    }
    free(line);
    line = NULL;
    cap = 0;
    n = getline(&line, &cap, fp);
    if (n < 0) {
        free(id);
        free(seq);
        free(line);
        return -1;
    }
    trim_line(line);
    char *qual = xstrndup(line, strlen(line));
    free(line);
    if (qual == NULL) {
        free(id);
        free(seq);
        return -1;
    }
    *id_out = id;
    *seq_out = seq;
    *qual_out = qual;
    return 1;
}

static int parse_gz_fastq_record(gzFile fp, char **id_out, char **seq_out, char **qual_out) {
    char buf[16384];
    if (gzgets(fp, buf, sizeof(buf)) == NULL) return 0;
    trim_line(buf);
    if (buf[0] != '@') return -1;
    char *id = xstrndup(buf + 1, strlen(buf + 1));
    if (id == NULL) return -1;
    if (gzgets(fp, buf, sizeof(buf)) == NULL) {
        free(id);
        return -1;
    }
    trim_line(buf);
    char *seq = xstrndup(buf, strlen(buf));
    if (seq == NULL) {
        free(id);
        return -1;
    }
    uppercase_ascii(seq);
    if (gzgets(fp, buf, sizeof(buf)) == NULL || buf[0] != '+') {
        free(id);
        free(seq);
        return -1;
    }
    if (gzgets(fp, buf, sizeof(buf)) == NULL) {
        free(id);
        free(seq);
        return -1;
    }
    trim_line(buf);
    char *qual = xstrndup(buf, strlen(buf));
    if (qual == NULL) {
        free(id);
        free(seq);
        return -1;
    }
    *id_out = id;
    *seq_out = seq;
    *qual_out = qual;
    return 1;
}

static int load_fastq_records(const char *path, seq_table *records) {
    const size_t len = strlen(path);
    if (len >= 3 && strcmp(path + len - 3, ".gz") == 0) {
        gzFile fp = gzopen(path, "rb");
        if (fp == NULL) return -1;
        int rc = 0;
        while (1) {
            char *id = NULL;
            char *seq = NULL;
            char *qual = NULL;
            int ok = parse_gz_fastq_record(fp, &id, &seq, &qual);
            free(qual);
            if (ok == 0) break;
            if (ok < 0 || push_record(records, id, strlen(id), seq, strlen(seq)) != 0) {
                free(id);
                free(seq);
                rc = -1;
                break;
            }
            free(id);
            free(seq);
        }
        gzclose(fp);
        return rc;
    }

    FILE *fp = fopen(path, "r");
    if (fp == NULL) return -1;
    int rc = 0;
    while (1) {
        char *id = NULL;
        char *seq = NULL;
        char *qual = NULL;
        int ok = parse_fastq_record(fp, &id, &seq, &qual);
        free(qual);
        if (ok == 0) break;
        if (ok < 0 || push_record(records, id, strlen(id), seq, strlen(seq)) != 0) {
            free(id);
            free(seq);
            rc = -1;
            break;
        }
        free(id);
        free(seq);
    }
    fclose(fp);
    return rc;
}

static int parse_text_fastq_record(FILE *fp, char **id_out, char **seq_out, char **qual_out) {
    return parse_fastq_record(fp, id_out, seq_out, qual_out);
}

static int parse_text_or_gz_fastq_record(const char *path, FILE **fp_out, gzFile *gz_out, int *is_gz_out,
                                         char **id_out, char **seq_out, char **qual_out) {
    if (*fp_out == NULL && *gz_out == NULL) {
        size_t len = strlen(path);
        if (len >= 3 && strcmp(path + len - 3, ".gz") == 0) {
            *gz_out = gzopen(path, "rb");
            if (*gz_out == NULL) return -1;
            *is_gz_out = 1;
        } else {
            *fp_out = fopen(path, "r");
            if (*fp_out == NULL) return -1;
            *is_gz_out = 0;
        }
    }
    if (*is_gz_out) return parse_gz_fastq_record(*gz_out, id_out, seq_out, qual_out);
    return parse_text_fastq_record(*fp_out, id_out, seq_out, qual_out);
}

static int slice_window(const char *seq, size_t seq_len, size_t start, size_t length, char **out_seq) {
    if (start > seq_len || length > seq_len - start) return -1;
    *out_seq = xstrndup(seq + start, length);
    return *out_seq == NULL ? -1 : 0;
}

static int all_targets_same_length(const seq_table *targets, size_t *length_out) {
    if (targets->count == 0) return -1;
    size_t len = targets->records[0].len;
    for (size_t i = 1; i < targets->count; ++i) {
        if (targets->records[i].len != len) return 0;
    }
    *length_out = len;
    return 1;
}

static int compare_lengths_desc(const void *a, const void *b) {
    const size_t la = *(const size_t *)a;
    const size_t lb = *(const size_t *)b;
    if (la < lb) return 1;
    if (la > lb) return -1;
    return 0;
}

static size_t collect_unique_target_lengths(const seq_table *targets, size_t **lengths_out) {
    if (targets->count == 0) {
        *lengths_out = NULL;
        return 0;
    }
    size_t *lengths = (size_t *)malloc(targets->count * sizeof(size_t));
    if (lengths == NULL) {
        *lengths_out = NULL;
        return 0;
    }
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
    qsort(lengths, count, sizeof(size_t), compare_lengths_desc);
    *lengths_out = lengths;
    return count;
}

static int is_prefix_of(const char *shorter, size_t shorter_len, const char *longer, size_t longer_len) {
    if (shorter_len > longer_len) return 0;
    return memcmp(shorter, longer, shorter_len) == 0;
}

static void write_assign_result(FILE *out, const char *mode, const char *read_id, const char *read_seq,
                                const seq_table *targets, qdaln_match_result r) {
    const char *target_seq = r.target_index >= 0 ? targets->records[r.target_index].seq : "";
    fprintf(out, "%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\n",
            mode, read_id, read_seq, r.target_index, target_seq, r.best_distance,
            status_name(r.status), r.match_count, r.second_best_distance);
}

static int run_fastq_assign(const char *argv0, int argc, char **argv) {
    (void)argv0;
    const char *barcodes_path = NULL;
    const char *reads_path = NULL;
    const char *out_path = NULL;
    size_t barcode_start = 0;
    size_t barcode_length = 0;
    int barcode_length_auto = 0;
    int k = 0;
    qdaln_metric metric = QDALN_METRIC_HAMMING;
    int max_correction_qual = -1;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--barcodes") == 0 && i + 1 < argc) {
            barcodes_path = argv[++i];
        } else if (strcmp(argv[i], "--reads") == 0 && i + 1 < argc) {
            reads_path = argv[++i];
        } else if (strcmp(argv[i], "--barcode-start") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &barcode_start) != 0) return 2;
        } else if (strcmp(argv[i], "--barcode-length") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "auto") == 0) {
                barcode_length_auto = 1;
            } else {
                if (parse_size_value(value, &barcode_length) != 0) return 2;
                barcode_length_auto = 0;
            }
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0) return 2;
        } else if (strcmp(argv[i], "--metric") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "hamming") == 0) {
                metric = QDALN_METRIC_HAMMING;
            } else if (strcmp(value, "levenshtein") == 0) {
                metric = QDALN_METRIC_LEVENSHTEIN;
            } else {
                return 2;
            }
        } else if (strcmp(argv[i], "--max-correction-qual") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &max_correction_qual) != 0) return 2;
        } else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else {
            return 2;
        }
    }

    if (barcodes_path == NULL || reads_path == NULL || out_path == NULL || k < 0) return 2;
    if (metric == QDALN_METRIC_HAMMING && k > 1) return 2;
    if (max_correction_qual > 93) return 2;

    seq_table targets = {0};
    if (read_target_table(barcodes_path, &targets) != 0) {
        free_table(&targets);
        return 1;
    }
    size_t fixed_length = 0;
    int same_length = all_targets_same_length(&targets, &fixed_length);
    if (!barcode_length_auto) {
        if (barcode_length == 0 && same_length > 0) barcode_length = fixed_length;
        if (barcode_length == 0) {
            free_table(&targets);
            return 2;
        }
    } else if (same_length > 0) {
        barcode_length = fixed_length;
        barcode_length_auto = 0;
    }

    const char **target_ptrs = (const char **)malloc(targets.count * sizeof(char *));
    size_t *target_lens = (size_t *)malloc(targets.count * sizeof(size_t));
    if ((targets.count != 0) && (target_ptrs == NULL || target_lens == NULL)) {
        free(target_ptrs);
        free(target_lens);
        free_table(&targets);
        return 1;
    }
    for (size_t i = 0; i < targets.count; ++i) {
        target_ptrs[i] = targets.records[i].seq;
        target_lens[i] = targets.records[i].len;
    }

    FILE *out = open_output_file(out_path);
    if (out == NULL) {
        free(target_ptrs);
        free(target_lens);
        free_table(&targets);
        return 1;
    }

    fprintf(out, "mode\tread_id\tread_seq\ttarget_index\ttarget_seq\tdistance\tstatus\tmatch_count\tsecond_best_distance\n");

    FILE *text_fp = NULL;
    gzFile gz_fp = NULL;
    int is_gz = 0;
    int rc = 0;
    size_t *candidate_lengths = NULL;
    size_t candidate_length_count = 0;
    if (barcode_length_auto) {
        candidate_length_count = collect_unique_target_lengths(&targets, &candidate_lengths);
        if (candidate_length_count == 0) rc = 1;
    }
    while (rc == 0) {
        char *id = NULL;
        char *seq = NULL;
        char *qual = NULL;
        int ok = parse_text_or_gz_fastq_record(reads_path, &text_fp, &gz_fp, &is_gz, &id, &seq, &qual);
        if (ok == 0) break;
        if (ok < 0) {
            rc = 1;
            break;
        }
        char *observed = NULL;
        qdaln_match_result r = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
        if (!barcode_length_auto) {
            if (slice_window(seq, strlen(seq), barcode_start, barcode_length, &observed) != 0) {
                write_assign_result(out, "fastq-assign", id, seq, &targets, r);
            } else {
                qdaln_match_result tmp = {0};
                if (qdaln_match_read_metric(observed, barcode_length, target_ptrs, target_lens, targets.count,
                                            k, QDALN_POLICY_BEST, metric, &tmp) != 0) {
                    rc = 1;
                } else {
                    r = tmp;
                    if (qual != NULL && max_correction_qual >= 0 && r.status == QDALN_MATCH_UNIQUE &&
                        r.target_index >= 0 && r.best_distance == 1) {
                        int qok = qdaln_read_correction_quality_ok(
                            observed, barcode_length,
                            targets.records[r.target_index].seq,
                            targets.records[r.target_index].len,
                            qual + barcode_start, barcode_length,
                            max_correction_qual,
                            metric);
                        if (qok == 0) {
                            r.target_index = -1;
                            r.best_distance = -1;
                            r.second_best_distance = -1;
                            r.match_count = 0;
                            r.status = QDALN_MATCH_NONE;
                        } else if (qok < 0) {
                            rc = 1;
                        }
                    }
                    if (rc == 0) write_assign_result(out, "fastq-assign", id, observed, &targets, r);
                }
            }
        } else {
            qdaln_match_result best = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
            int saw_valid = 0;
            int ambiguous_prefix = 0;
            size_t ambiguous_count = 0;
            for (size_t li = 0; li < candidate_length_count && rc == 0; ++li) {
                size_t current_length = candidate_lengths[li];
                char *candidate = NULL;
                if (slice_window(seq, strlen(seq), barcode_start, current_length, &candidate) != 0) continue;
                qdaln_match_result tmp = {0};
                if (qdaln_match_read_metric(candidate, current_length, target_ptrs, target_lens, targets.count,
                                            k, QDALN_POLICY_BEST, metric, &tmp) != 0) {
                    free(candidate);
                    rc = 1;
                    break;
                }
                if (qual != NULL && max_correction_qual >= 0 && tmp.status == QDALN_MATCH_UNIQUE &&
                    tmp.target_index >= 0 && tmp.best_distance == 1) {
                    int qok = qdaln_read_correction_quality_ok(
                        candidate, current_length,
                        targets.records[tmp.target_index].seq,
                        targets.records[tmp.target_index].len,
                        qual + barcode_start, current_length,
                        max_correction_qual,
                        metric);
                    if (qok == 0) {
                        tmp.target_index = -1;
                        tmp.best_distance = -1;
                        tmp.second_best_distance = -1;
                        tmp.match_count = 0;
                        tmp.status = QDALN_MATCH_NONE;
                    } else if (qok < 0) {
                        free(candidate);
                        rc = 1;
                        break;
                    }
                }
                if (tmp.status == QDALN_MATCH_UNIQUE && tmp.target_index >= 0) {
                    if (!saw_valid) {
                        best = tmp;
                        observed = candidate;
                        saw_valid = 1;
                    } else {
                        const seq_record *best_target = &targets.records[best.target_index];
                        const seq_record *next_target = &targets.records[tmp.target_index];
                        if (best_target->len != next_target->len &&
                            is_prefix_of(best_target->seq, best_target->len, next_target->seq, next_target->len)) {
                            ambiguous_prefix = 1;
                            ++ambiguous_count;
                        } else if (best_target->len != next_target->len &&
                                   is_prefix_of(next_target->seq, next_target->len, best_target->seq, best_target->len)) {
                            ambiguous_prefix = 1;
                            ++ambiguous_count;
                        }
                        free(candidate);
                    }
                } else {
                    free(candidate);
                }
            }
            if (rc == 0) {
                if (ambiguous_prefix) {
                    qdaln_match_result amb = {-1, -1, -1, (int)(ambiguous_count + 1), QDALN_MATCH_AMBIGUOUS};
                    write_assign_result(out, "fastq-assign", id, observed ? observed : seq, &targets, amb);
                } else if (saw_valid) {
                    write_assign_result(out, "fastq-assign", id, observed, &targets, best);
                } else {
                    write_assign_result(out, "fastq-assign", id, seq, &targets, r);
                }
            }
        }
        free(id);
        free(seq);
        free(qual);
        free(observed);
    }

    if (text_fp != NULL) fclose(text_fp);
    if (gz_fp != NULL) gzclose(gz_fp);
    free(candidate_lengths);
    fclose(out);
    free(target_ptrs);
    free(target_lens);
    free_table(&targets);
    return rc;
}

static int compute_edit_kind(const char *observed, size_t observed_len, const char *target, size_t target_len) {
    if (observed_len == target_len && memcmp(observed, target, observed_len) == 0) return 0;
    if (observed_len == target_len) return 1;
    if (observed_len == target_len + 1) return 2;
    if (observed_len + 1 == target_len) return 3;
    return 4;
}

static int target_ambiguity_flags(const seq_table *targets, int k, int *flags) {
    for (size_t i = 0; i < targets->count; ++i) flags[i] = 0;
    if (k < 1) {
        for (size_t i = 0; i < targets->count; ++i) {
            for (size_t j = i + 1; j < targets->count; ++j) {
                if (targets->records[i].len == targets->records[j].len &&
                    memcmp(targets->records[i].seq, targets->records[j].seq, targets->records[i].len) == 0) {
                    flags[i] = 1;
                    flags[j] = 1;
                }
            }
        }
        return 0;
    }
    for (size_t i = 0; i < targets->count; ++i) {
        for (size_t j = i + 1; j < targets->count; ++j) {
            int d = qdaln_edit_distance(targets->records[i].seq, targets->records[i].len,
                                        targets->records[j].seq, targets->records[j].len);
            if (d < 0) return -1;
            if (d <= k) {
                flags[i] = 1;
                flags[j] = 1;
            }
        }
    }
    return 0;
}

static int parse_sample_labels(int argc, char **argv, size_t reads_count, char ***labels_out) {
    char **labels = (char **)calloc(reads_count, sizeof(char *));
    if (labels == NULL) return -1;
    const char *csv = NULL;
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--sample-label") == 0 && i + 1 < argc) {
            csv = argv[++i];
            break;
        }
    }
    if (csv == NULL) {
        for (size_t i = 0; i < reads_count; ++i) {
            char buf[32];
            int n = snprintf(buf, sizeof(buf), "sample_%zu", i + 1);
            labels[i] = xstrndup(buf, (size_t)n);
            if (labels[i] == NULL) {
                for (size_t j = 0; j < i; ++j) free(labels[j]);
                free(labels);
                return -1;
            }
        }
        *labels_out = labels;
        return 0;
    }
    const char *start = csv;
    size_t idx = 0;
    while (*start != '\0' && idx < reads_count) {
        const char *comma = strchr(start, ',');
        size_t len = comma == NULL ? strlen(start) : (size_t)(comma - start);
        labels[idx] = xstrndup(start, len);
        if (labels[idx] == NULL) {
            for (size_t j = 0; j < idx; ++j) free(labels[j]);
            free(labels);
            return -1;
        }
        ++idx;
        if (comma == NULL) break;
        start = comma + 1;
    }
    if (idx != reads_count) {
        for (size_t j = 0; j < idx; ++j) free(labels[j]);
        free(labels);
        return -1;
    }
    *labels_out = labels;
    return 0;
}

static void free_sample_labels(char **labels, size_t count) {
    if (labels == NULL) return;
    for (size_t i = 0; i < count; ++i) free(labels[i]);
    free(labels);
}

static int write_count_table(FILE *out, const seq_table *targets, const int *exact, const int *sub, const int *ins,
                             const int *del, const int *other, const int *amb_flags) {
    fprintf(out, "target_id\ttarget_seq\tgene\tcount_exact\tcount_corrected_substitution\t"
                 "count_corrected_insertion\tcount_corrected_deletion\tcount_corrected_other\tcount_total\tambiguous_nearby\n");
    for (size_t i = 0; i < targets->count; ++i) {
        int total = exact[i] + sub[i] + ins[i] + del[i] + other[i];
        fprintf(out, "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                targets->records[i].id, targets->records[i].seq, targets->records[i].gene,
                exact[i], sub[i], ins[i], del[i], other[i], total, amb_flags[i]);
    }
    return 0;
}

static int write_mageck_table(FILE *out, const seq_table *targets, char **sample_labels, size_t sample_count,
                              const int *counts) {
    fprintf(out, "sgRNA\tGene");
    for (size_t s = 0; s < sample_count; ++s) fprintf(out, "\t%s", sample_labels[s]);
    fprintf(out, "\n");
    for (size_t t = 0; t < targets->count; ++t) {
        fprintf(out, "%s\t%s", targets->records[t].id, targets->records[t].gene[0] ? targets->records[t].gene : "NA");
        for (size_t s = 0; s < sample_count; ++s) {
            fprintf(out, "\t%d", counts[s * targets->count + t]);
        }
        fprintf(out, "\n");
    }
    return 0;
}

static int parse_reads_option(int argc, char **argv, const char ***reads_out, size_t *count_out) {
    size_t cap = 4;
    size_t count = 0;
    const char **reads = (const char **)malloc(cap * sizeof(char *));
    if (reads == NULL) return -1;
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--reads") == 0 && i + 1 < argc) {
            if (count == cap) {
                cap *= 2;
                const char **next = (const char **)realloc((void *)reads, cap * sizeof(char *));
                if (next == NULL) {
                    free((void *)reads);
                    return -1;
                }
                reads = next;
            }
            reads[count++] = argv[++i];
        }
    }
    if (count == 0) {
        free((void *)reads);
        return -1;
    }
    *reads_out = reads;
    *count_out = count;
    return 0;
}

static int parse_target_window(const char *seq, size_t seq_len, size_t start, size_t length,
                               int auto_offset_enabled, size_t offset_window,
                               const qdaln_index *index, qdaln_metric metric,
                               int max_correction_qual,
                               const char *qual, int *best_offset_out, char **window_out) {
    if (!auto_offset_enabled) {
        *best_offset_out = (int)start;
        return slice_window(seq, seq_len, start, length, window_out);
    }
    int best_matches = -1;
    size_t best_start = start;
    char *best_window = NULL;
    for (size_t offset = 0; offset <= offset_window * 2; ++offset) {
        size_t candidate_start = start + offset;
        if (offset > offset_window) {
            size_t delta = offset - offset_window;
            if (delta > start) continue;
            candidate_start = start - delta;
        }
        char *candidate = NULL;
        if (slice_window(seq, seq_len, candidate_start, length, &candidate) != 0) continue;
        qdaln_match_result tmp = {0};
        if (qdaln_index_match_read(index, candidate, length, 0, QDALN_POLICY_BEST, metric, &tmp) != 0) {
            free(candidate);
            free(best_window);
            return -1;
        }
        if (qual != NULL && max_correction_qual >= 0 && tmp.status == QDALN_MATCH_UNIQUE &&
            tmp.target_index >= 0 && tmp.best_distance == 1) {
            int qok = qdaln_read_correction_quality_ok(
                candidate, length,
                index->targets[tmp.target_index],
                index->target_lens[tmp.target_index],
                qual + candidate_start, length,
                max_correction_qual,
                metric);
            if (qok == 0) {
                tmp.target_index = -1;
                tmp.best_distance = -1;
                tmp.second_best_distance = -1;
                tmp.match_count = 0;
                tmp.status = QDALN_MATCH_NONE;
            } else if (qok < 0) {
                free(candidate);
                free(best_window);
                return -1;
            }
        }
        int score = tmp.status == QDALN_MATCH_UNIQUE ? 1 : 0;
        if (score > best_matches) {
            best_matches = score;
            best_start = candidate_start;
            free(best_window);
            best_window = candidate;
        } else {
            free(candidate);
        }
    }
    if (best_window == NULL) return -1;
    *best_offset_out = (int)best_start;
    *window_out = best_window;
    return 0;
}

static int target_total_count(const int *exact, const int *sub, const int *ins, const int *del, const int *other,
                              size_t index) {
    return exact[index] + sub[index] + ins[index] + del[index] + other[index];
}

static double gini_index(const int *values, size_t n) {
    long long total = 0;
    for (size_t i = 0; i < n; ++i) total += values[i];
    if (n == 0 || total == 0) return 0.0;
    int *sorted = (int *)malloc(n * sizeof(int));
    if (sorted == NULL) return 0.0;
    memcpy(sorted, values, n * sizeof(int));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (sorted[i] > sorted[j]) {
                int tmp = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = tmp;
            }
        }
    }
    long long weighted = 0;
    for (size_t i = 0; i < n; ++i) {
        weighted += (long long)(i + 1) * (long long)sorted[i];
    }
    free(sorted);
    return (2.0 * (double)weighted) / ((double)n * (double)total) - ((double)n + 1.0) / (double)n;
}

static int cmp_int_desc(const void *a, const void *b) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;
    if (ia < ib) return 1;
    if (ia > ib) return -1;
    return 0;
}

static double top_one_percent_dominance(const int *values, size_t n) {
    long long total = 0;
    for (size_t i = 0; i < n; ++i) total += values[i];
    if (n == 0 || total == 0) return 0.0;
    int *sorted = (int *)malloc(n * sizeof(int));
    if (sorted == NULL) return 0.0;
    memcpy(sorted, values, n * sizeof(int));
    qsort(sorted, n, sizeof(int), cmp_int_desc);
    size_t top_n = n / 100;
    if (top_n == 0) top_n = 1;
    long long top_total = 0;
    for (size_t i = 0; i < top_n; ++i) top_total += sorted[i];
    free(sorted);
    return (double)top_total / (double)total;
}

static void count_provenance_fields(const qdaln_match_result *r, const char *observed, size_t observed_len,
                                    const seq_table *targets, int *exact, int *sub, int *ins, int *del, int *other) {
    if (r->status != QDALN_MATCH_UNIQUE || r->target_index < 0) return;
    const seq_record *target = &targets->records[r->target_index];
    int kind = compute_edit_kind(observed, observed_len, target->seq, target->len);
    switch (kind) {
        case 0:
            exact[r->target_index] += 1;
            break;
        case 1:
            sub[r->target_index] += 1;
            break;
        case 2:
            ins[r->target_index] += 1;
            break;
        case 3:
            del[r->target_index] += 1;
            break;
        default:
            other[r->target_index] += 1;
            break;
    }
}

static void print_json_string(FILE *out, const char *s) {
    fputc('"', out);
    for (; *s != '\0'; ++s) {
        unsigned char c = (unsigned char)*s;
        switch (c) {
            case '\\':
                fputs("\\\\", out);
                break;
            case '"':
                fputs("\\\"", out);
                break;
            case '\n':
                fputs("\\n", out);
                break;
            case '\r':
                fputs("\\r", out);
                break;
            case '\t':
                fputs("\\t", out);
                break;
            default:
                if (c < 0x20) {
                    fprintf(out, "\\u%04x", c);
                } else {
                    fputc(c, out);
                }
        }
    }
    fputc('"', out);
}

static int read_json_fastq_list(const char *path, seq_table *table) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return -1;
    char line[4096];
    while (fgets(line, sizeof(line), fp) != NULL) {
        trim_line(line);
        if (line[0] == '\0') continue;
        const char *fastq = strstr(line, "\"fastq\"");
        const char *sample = strstr(line, "\"sample_id\"");
        if (fastq == NULL || sample == NULL) continue;
        const char *sample_q = strchr(sample + 11, '"');
        if (sample_q == NULL) continue;
        const char *sample_end = strchr(sample_q + 1, '"');
        if (sample_end == NULL) continue;
        const char *fastq_q = strchr(fastq + 7, '"');
        if (fastq_q == NULL) continue;
        const char *fastq_end = strchr(fastq_q + 1, '"');
        if (fastq_end == NULL) continue;
        if (push_record(table, sample_q + 1, (size_t)(sample_end - sample_q - 1),
                        fastq_q + 1, (size_t)(fastq_end - fastq_q - 1)) != 0) {
            fclose(fp);
            return -1;
        }
    }
    fclose(fp);
    return table->count == 0 ? -1 : 0;
}

static int read_samples_table(const char *path, seq_table *table) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return -1;
    char buf[8192];
    int first_data = 1;
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        trim_line(buf);
        if (buf[0] == '\0' || buf[0] == '#') continue;
        char *fields[8];
        size_t nf = split_fields(buf, '\t', fields, 8);
        if (first_data) {
            if (nf >= 2 && (field_eq(fields[0], "sample_id") || field_eq(fields[0], "sample")) &&
                (field_eq(fields[1], "fastq") || field_eq(fields[1], "reads"))) {
                first_data = 0;
                continue;
            }
        }
        first_data = 0;
        if (nf < 2) {
            fclose(fp);
            return -1;
        }
        if (push_record(table, fields[0], strlen(fields[0]), fields[1], strlen(fields[1])) != 0) {
            fclose(fp);
            return -1;
        }
    }
    fclose(fp);
    return table->count == 0 ? -1 : 0;
}

static int load_sample_manifest(const char *path, seq_table *table) {
    const char *ext = strrchr(path, '.');
    if (ext != NULL && strcmp(ext, ".json") == 0) return read_json_fastq_list(path, table);
    return read_samples_table(path, table);
}

static int write_summary_json(FILE *out, long long total_reads, long long exact_reads, long long rescued_reads,
                              long long ambiguous_reads, long long none_reads, long long invalid_reads,
                              long long candidates_considered, long long candidates_verified,
                              long long offset_rescues, const char *sample_label) {
    fprintf(out, "{\n");
    fprintf(out, "  \"sample_label\": ");
    print_json_string(out, sample_label);
    fprintf(out, ",\n  \"total_reads\": %lld,\n", total_reads);
    fprintf(out, "  \"assigned_exact\": %lld,\n", exact_reads);
    fprintf(out, "  \"assigned_rescued\": %lld,\n", rescued_reads);
    fprintf(out, "  \"ambiguous_reads\": %lld,\n", ambiguous_reads);
    fprintf(out, "  \"none_reads\": %lld,\n", none_reads);
    fprintf(out, "  \"invalid_reads\": %lld,\n", invalid_reads);
    fprintf(out, "  \"candidates_considered\": %lld,\n", candidates_considered);
    fprintf(out, "  \"candidates_verified\": %lld,\n", candidates_verified);
    fprintf(out, "  \"auto_offset_rescues\": %lld\n", offset_rescues);
    fprintf(out, "}\n");
    return 0;
}

static int write_sample_qc(FILE *out, const char *sample_label, long long total_reads, long long unique_reads,
                           long long exact_reads, long long rescued_reads, long long ambiguous_reads,
                           long long none_reads, size_t target_count, const int *totals,
                           long long candidates_considered, long long candidates_verified,
                           long long offset_rescues) {
    size_t zero_targets = 0;
    int covered = 0;
    for (size_t i = 0; i < target_count; ++i) {
        if (totals[i] == 0) {
            ++zero_targets;
        } else {
            ++covered;
        }
    }
    double assignment_rate = total_reads ? (double)unique_reads / (double)total_reads : 0.0;
    double exact_rate = total_reads ? (double)exact_reads / (double)total_reads : 0.0;
    double rescued_rate = total_reads ? (double)rescued_reads / (double)total_reads : 0.0;
    double ambiguous_rate = total_reads ? (double)ambiguous_reads / (double)total_reads : 0.0;
    double none_rate = total_reads ? (double)none_reads / (double)total_reads : 0.0;
    double coverage_rate = target_count ? (double)covered / (double)target_count : 0.0;
    double gini = gini_index(totals, target_count);
    double dominance = top_one_percent_dominance(totals, target_count);
    fprintf(out, "%s\t%lld\t%lld\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%zu\t%.6f\t%.6f\t%lld\t%lld\t%lld\n",
            sample_label, total_reads, unique_reads, assignment_rate, exact_rate, rescued_rate,
            ambiguous_rate, none_rate, coverage_rate, zero_targets, gini, dominance,
            candidates_considered, candidates_verified, offset_rescues);
    return 0;
}

static int write_target_counts_long(FILE *out, const char *sample_label, const seq_table *targets,
                                    const int *exact, const int *sub, const int *ins, const int *del,
                                    const int *other, const int *amb_flags) {
    for (size_t i = 0; i < targets->count; ++i) {
        int total = exact[i] + sub[i] + ins[i] + del[i] + other[i];
        fprintf(out, "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                sample_label,
                targets->records[i].id,
                targets->records[i].seq,
                targets->records[i].gene,
                exact[i], sub[i], ins[i], del[i], other[i], total, amb_flags[i]);
    }
    return 0;
}

static int write_report_header(FILE *out) {
    fprintf(out,
            "<!doctype html>\n<html><head><meta charset=\"utf-8\"><title>DotMatch Report</title>"
            "<style>body{font-family:system-ui, sans-serif;max-width:960px;margin:2rem auto;padding:0 1rem;color:#1b1f23;}"
            "table{border-collapse:collapse;width:100%;margin:1rem 0;}th,td{border:1px solid #d0d7de;padding:0.5rem;text-align:left;}"
            "th{background:#f6f8fa;} .warn{color:#b42318;font-weight:600;} .mono{font-family:ui-monospace, monospace;}"
            "</style></head><body>\n<h1>DotMatch Count Report</h1>\n");
    return 0;
}

static int write_report_footer(FILE *out) {
    fprintf(out, "</body></html>\n");
    return 0;
}

static int ensure_dir_exists(const char *path) {
    struct stat st;
    if (stat(path, &st) == 0) {
        return S_ISDIR(st.st_mode) ? 0 : -1;
    }
    return mkdir(path, 0700);
}

static int write_text_file(const char *path, const char *text) {
    FILE *out = open_output_file(path);
    if (out == NULL) return -1;
    fputs(text, out);
    fclose(out);
    return 0;
}

static int write_audit_summary(const char *audit_dir, const seq_table *targets, int k) {
    if (ensure_dir_exists(audit_dir) != 0) return -1;
    char path[PATH_MAX];
    snprintf(path, sizeof(path), "%s/audit_summary.tsv", audit_dir);
    FILE *out = open_output_file(path);
    if (out == NULL) return -1;
    fprintf(out, "target_id\ttarget_seq\tmin_distance\n");
    for (size_t i = 0; i < targets->count; ++i) {
        int min_d = -1;
        for (size_t j = 0; j < targets->count; ++j) {
            if (i == j) continue;
            int d = qdaln_edit_distance(targets->records[i].seq, targets->records[i].len,
                                        targets->records[j].seq, targets->records[j].len);
            if (d < 0) {
                fclose(out);
                return -1;
            }
            if (min_d < 0 || d < min_d) min_d = d;
        }
        fprintf(out, "%s\t%s\t%d\n", targets->records[i].id, targets->records[i].seq, min_d);
    }
    fclose(out);
    snprintf(path, sizeof(path), "%s/audit_summary.json", audit_dir);
    char json_buf[256];
    snprintf(json_buf, sizeof(json_buf), "{\n  \"k\": %d,\n  \"target_count\": %zu\n}\n", k, targets->count);
    return write_text_file(path, json_buf);
}

static int cmp_pair_qsort(const void *a, const void *b) {
    const qdaln_audit_pair *pa = (const qdaln_audit_pair *)a;
    const qdaln_audit_pair *pb = (const qdaln_audit_pair *)b;
    if (pa->distance != pb->distance) return pa->distance - pb->distance;
    if (pa->left != pb->left) return (pa->left < pb->left) ? -1 : 1;
    if (pa->right != pb->right) return (pa->right < pb->right) ? -1 : 1;
    return 0;
}

static int write_collision_pairs(const char *audit_dir, const seq_table *targets, qdaln_audit_pair *pairs, size_t pair_count) {
    char path[PATH_MAX];
    snprintf(path, sizeof(path), "%s/collision_pairs.tsv", audit_dir);
    FILE *out = open_output_file(path);
    if (out == NULL) return -1;
    fprintf(out, "left_id\tleft_seq\tright_id\tright_seq\tdistance\n");
    qsort(pairs, pair_count, sizeof(qdaln_audit_pair), cmp_pair_qsort);
    for (size_t i = 0; i < pair_count; ++i) {
        fprintf(out, "%s\t%s\t%s\t%s\t%d\n",
                targets->records[pairs[i].left].id,
                targets->records[pairs[i].left].seq,
                targets->records[pairs[i].right].id,
                targets->records[pairs[i].right].seq,
                pairs[i].distance);
    }
    fclose(out);
    return 0;
}

static int write_collision_clusters(const char *audit_dir, const seq_table *targets, qdaln_audit_cluster *clusters, size_t cluster_count) {
    char path[PATH_MAX];
    snprintf(path, sizeof(path), "%s/collision_clusters.tsv", audit_dir);
    FILE *out = open_output_file(path);
    if (out == NULL) return -1;
    fprintf(out, "cluster_index\tmembers\n");
    for (size_t i = 0; i < cluster_count; ++i) {
        fprintf(out, "%zu", i);
        for (size_t j = 0; j < clusters[i].count; ++j) {
            size_t idx = clusters[i].members[j];
            fprintf(out, "\t%s", targets->records[idx].id);
        }
        fprintf(out, "\n");
    }
    fclose(out);
    return 0;
}

static int write_target_safety(const char *audit_dir, const seq_table *targets, const int *k0_safe,
                               const int *k1_safe, const int *k2_safe, const int *nearest_distance) {
    char path[PATH_MAX];
    snprintf(path, sizeof(path), "%s/target_safety.tsv", audit_dir);
    FILE *out = open_output_file(path);
    if (out == NULL) return -1;
    fprintf(out, "target_id\ttarget_seq\tnearest_distance\tsafe_k0\tsafe_k1\tsafe_k2\n");
    for (size_t i = 0; i < targets->count; ++i) {
        fprintf(out, "%s\t%s\t%d\t%s\t%s\t%s\n",
                targets->records[i].id,
                targets->records[i].seq,
                nearest_distance[i],
                k0_safe[i] ? "true" : "false",
                k1_safe[i] ? "true" : "false",
                k2_safe[i] ? "true" : "false");
    }
    fclose(out);
    return 0;
}

static int write_ambiguous_variants(const char *audit_dir, const seq_table *targets, const qdaln_audit_variant *variants,
                                    size_t variant_count) {
    char path[PATH_MAX];
    snprintf(path, sizeof(path), "%s/ambiguous_variants.tsv", audit_dir);
    FILE *out = open_output_file(path);
    if (out == NULL) return -1;
    fprintf(out, "query_variant\tmatch_count\ttarget_ids\n");
    for (size_t i = 0; i < variant_count; ++i) {
        fprintf(out, "%s\t%zu", variants[i].query, variants[i].count);
        for (size_t j = 0; j < variants[i].count; ++j) {
            fprintf(out, "\t%s", targets->records[variants[i].targets[j]].id);
        }
        fprintf(out, "\n");
    }
    fclose(out);
    return 0;
}

static int write_audit_json_summary(const char *audit_dir, int k, size_t target_count,
                                    size_t pairs_at_0, size_t pairs_at_1, size_t pairs_at_2,
                                    size_t clusters, size_t ambiguous_query_variants_k1,
                                    int library_safe_k0, int library_safe_k1, int library_safe_k2,
                                    const char *mode) {
    char path[PATH_MAX];
    snprintf(path, sizeof(path), "%s/audit_summary.json", audit_dir);
    FILE *out = open_output_file(path);
    if (out == NULL) return -1;
    fprintf(out,
            "{\n"
            "  \"target_count\": %zu,\n"
            "  \"k\": %d,\n"
            "  \"pairs_at_distance_0\": %zu,\n"
            "  \"pairs_at_distance_1\": %zu,\n"
            "  \"pairs_at_distance_2\": %zu,\n"
            "  \"collision_clusters\": %zu,\n"
            "  \"ambiguous_query_variants_k1\": %zu,\n"
            "  \"library_safe_k0\": %s,\n"
            "  \"library_safe_k1\": %s,\n"
            "  \"library_safe_k2\": %s,\n"
            "  \"audit_mode\": ",
            target_count, k, pairs_at_0, pairs_at_1, pairs_at_2, clusters,
            ambiguous_query_variants_k1,
            library_safe_k0 ? "true" : "false",
            library_safe_k1 ? "true" : "false",
            library_safe_k2 ? "true" : "false");
    print_json_string(out, mode);
    fprintf(out, "\n}\n");
    fclose(out);
    return 0;
}

static int write_count_html_report(const char *path, const seq_table *targets, const int *exact, const int *sub,
                                   const int *ins, const int *del, const int *other,
                                   const int *amb_flags, const char *sample_label,
                                   long long total_reads, long long unique_reads,
                                   long long exact_reads, long long rescued_reads,
                                   long long ambiguous_reads, long long none_reads,
                                   long long invalid_reads, long long candidates_considered,
                                   long long candidates_verified, long long offset_rescues,
                                   const int *sample_totals,
                                   const char *audit_dir, const char *unmatched_report_path) {
    int out_fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
    if (out_fd < 0) return -1;
    FILE *out = fdopen(out_fd, "w");
    if (out == NULL) {
        close(out_fd);
        return -1;
    }

    double assignment_rate = total_reads ? (double)unique_reads / (double)total_reads : 0.0;
    double ambiguous_rate = total_reads ? (double)ambiguous_reads / (double)total_reads : 0.0;
    double none_rate = total_reads ? (double)none_reads / (double)total_reads : 0.0;
    size_t covered = 0;
    for (size_t i = 0; i < targets->count; ++i) {
        if (sample_totals[i] > 0) ++covered;
    }
    double coverage_rate = targets->count ? (double)covered / (double)targets->count : 0.0;

    write_report_header(out);
    fprintf(out, "<p class=\"mono\">sample: ");
    print_json_string(out, sample_label);
    fprintf(out, "</p>\n");
    fprintf(out, "<table><tr><th>Metric</th><th>Value</th></tr>");
    fprintf(out, "<tr><td>Total reads</td><td>%lld</td></tr>", total_reads);
    fprintf(out, "<tr><td>Unique assigned</td><td>%lld (%.2f%%)</td></tr>", unique_reads, assignment_rate * 100.0);
    fprintf(out, "<tr><td>Exact assigned</td><td>%lld</td></tr>", exact_reads);
    fprintf(out, "<tr><td>Rescued assigned</td><td>%lld</td></tr>", rescued_reads);
    fprintf(out, "<tr><td>Ambiguous</td><td>%lld (%.2f%%)</td></tr>", ambiguous_reads, ambiguous_rate * 100.0);
    fprintf(out, "<tr><td>No match</td><td>%lld (%.2f%%)</td></tr>", none_reads, none_rate * 100.0);
    fprintf(out, "<tr><td>Invalid</td><td>%lld</td></tr>", invalid_reads);
    fprintf(out, "<tr><td>Candidate reads considered</td><td>%lld</td></tr>", candidates_considered);
    fprintf(out, "<tr><td>Candidate targets verified</td><td>%lld</td></tr>", candidates_verified);
    fprintf(out, "<tr><td>Auto-offset rescues</td><td>%lld</td></tr>", offset_rescues);
    fprintf(out, "<tr><td>Library coverage</td><td>%zu/%zu (%.2f%%)</td></tr></table>\n",
            covered, targets->count, coverage_rate * 100.0);
    if (ambiguous_rate > 0.10) {
        fprintf(out, "<p class=\"warn\">Warning: ambiguous rate is above 10%%.</p>\n");
    }
    if (none_rate > 0.10) {
        fprintf(out, "<p class=\"warn\">Warning: no-match rate is above 10%%.</p>\n");
    }
    fprintf(out, "<h2>Target counts</h2>\n<table><tr><th>Target</th><th>Gene</th><th>Exact</th><th>Sub</th><th>Ins</th><th>Del</th><th>Other</th><th>Total</th><th>Ambiguous nearby</th></tr>");
    for (size_t i = 0; i < targets->count; ++i) {
        fprintf(out,
                "<tr><td>%s</td><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%s</td></tr>\n",
                targets->records[i].id,
                targets->records[i].gene,
                exact[i], sub[i], ins[i], del[i], other[i],
                exact[i] + sub[i] + ins[i] + del[i] + other[i],
                amb_flags[i] ? "yes" : "no");
    }
    fprintf(out, "</table>\n");
    if (audit_dir != NULL) {
        fprintf(out, "<p>Audit directory: %s</p>\n", audit_dir);
    }
    if (unmatched_report_path != NULL) {
        fprintf(out, "<p>Top unmatched report: %s</p>\n", unmatched_report_path);
    }
    write_report_footer(out);
    fclose(out);
    return 0;
}

static int load_unmatched_top(const char *path, seq_table *table) {
    return read_table(path, table);
}

static int nearest_target(const seq_table *targets, const char *query, size_t query_len, int *distance_out) {
    int best_idx = -1;
    int best_dist = -1;
    for (size_t i = 0; i < targets->count; ++i) {
        int d = qdaln_edit_distance(query, query_len, targets->records[i].seq, targets->records[i].len);
        if (d < 0) return -1;
        if (best_idx < 0 || d < best_dist) {
            best_idx = (int)i;
            best_dist = d;
        }
    }
    *distance_out = best_dist;
    return best_idx;
}

static char complement_base(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

static char *reverse_complement(const char *seq, size_t len) {
    char *out = (char *)malloc(len + 1);
    if (out == NULL) return NULL;
    for (size_t i = 0; i < len; ++i) out[i] = complement_base(seq[len - 1 - i]);
    out[len] = '\0';
    return out;
}

static int has_low_quality(const char *qual, size_t len, int threshold) {
    if (qual == NULL || threshold < 0) return 0;
    for (size_t i = 0; i < len; ++i) {
        int q = (int)((unsigned char)qual[i]) - 33;
        if (q <= threshold) return 1;
    }
    return 0;
}

static int contains_char(const char *seq, char c) {
    return strchr(seq, c) != NULL;
}

static int prefix_matches(const char *seq, size_t seq_len, const char *adapter, size_t adapter_len) {
    if (adapter_len > seq_len) return 0;
    return memcmp(seq, adapter, adapter_len) == 0;
}

static const char *reason_label_for_unmatched(int distance, int k, int rc_match, int offset_match,
                                              int adapter_match, int low_qual, int contains_n, int wrong_length) {
    if (wrong_length) return "wrong_length";
    if (contains_n) return "contains_N";
    if (distance >= 0 && distance > k) return "near_known_target_above_k";
    if (rc_match) return "reverse_complement_candidate";
    if (offset_match) return "offset_shift_candidate";
    if (adapter_match) return "adapter_or_primer_candidate";
    if (low_qual) return "low_quality_candidate";
    return "other";
}

static int run_inspect_unmatched(const char *argv0, int argc, char **argv) {
    (void)argv0;
    const char *targets_path = NULL;
    const char *reads_path = NULL;
    const char *out_path = NULL;
    const char *adapter = NULL;
    size_t target_start = 0;
    size_t target_length = 0;
    size_t offset_window = 0;
    int k = 0;
    int low_quality_threshold = -1;
    size_t top_n = 100;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--targets") == 0 && i + 1 < argc) {
            targets_path = argv[++i];
        } else if (strcmp(argv[i], "--reads") == 0 && i + 1 < argc) {
            reads_path = argv[++i];
        } else if (strcmp(argv[i], "--target-start") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_start) != 0) return 2;
        } else if (strcmp(argv[i], "--target-length") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_length) != 0) return 2;
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0) return 2;
        } else if (strcmp(argv[i], "--offset-window") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &offset_window) != 0) return 2;
        } else if (strcmp(argv[i], "--adapter") == 0 && i + 1 < argc) {
            adapter = argv[++i];
        } else if (strcmp(argv[i], "--low-quality-threshold") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &low_quality_threshold) != 0) return 2;
        } else if (strcmp(argv[i], "--top") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &top_n) != 0) return 2;
        } else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else {
            return 2;
        }
    }
    if (targets_path == NULL || reads_path == NULL || out_path == NULL || target_length == 0) return 2;

    seq_table targets = {0};
    if (read_target_table(targets_path, &targets) != 0) {
        free_table(&targets);
        return 1;
    }
    qdaln_index *index = qdaln_index_build((const char *const *)(&targets.records[0].seq), NULL, 0);
    (void)index;

    typedef struct unmatched_row {
        char *seq;
        int count;
    } unmatched_row;

    unmatched_row *rows = NULL;
    size_t row_count = 0;
    size_t row_cap = 0;
    FILE *text_fp = NULL;
    gzFile gz_fp = NULL;
    int is_gz = 0;
    while (1) {
        char *id = NULL;
        char *seq = NULL;
        char *qual = NULL;
        int ok = parse_text_or_gz_fastq_record(reads_path, &text_fp, &gz_fp, &is_gz, &id, &seq, &qual);
        free(id);
        if (ok == 0) break;
        if (ok < 0) {
            free(seq);
            free(qual);
            break;
        }
        char *window = NULL;
        int wrong_length = slice_window(seq, strlen(seq), target_start, target_length, &window) != 0;
        if (!wrong_length) {
            int already = 0;
            for (size_t i = 0; i < row_count; ++i) {
                if (strcmp(rows[i].seq, window) == 0) {
                    rows[i].count += 1;
                    already = 1;
                    break;
                }
            }
            if (!already) {
                if (row_count == row_cap) {
                    size_t next_cap = row_cap == 0 ? 64 : row_cap * 2;
                    unmatched_row *next = (unmatched_row *)realloc(rows, next_cap * sizeof(unmatched_row));
                    if (next == NULL) {
                        free(window);
                        free(seq);
                        free(qual);
                        break;
                    }
                    rows = next;
                    row_cap = next_cap;
                }
                rows[row_count].seq = window;
                rows[row_count].count = 1;
                ++row_count;
                window = NULL;
            }
        }
        free(window);
        free(seq);
        free(qual);
    }
    if (text_fp != NULL) fclose(text_fp);
    if (gz_fp != NULL) gzclose(gz_fp);

    for (size_t i = 0; i < row_count; ++i) {
        for (size_t j = i + 1; j < row_count; ++j) {
            if (rows[i].count < rows[j].count) {
                unmatched_row tmp = rows[i];
                rows[i] = rows[j];
                rows[j] = tmp;
            }
        }
    }

    FILE *out = open_output_file(out_path);
    if (out == NULL) {
        for (size_t i = 0; i < row_count; ++i) free(rows[i].seq);
        free(rows);
        free_table(&targets);
        return 1;
    }
    fprintf(out, "observed_seq\tcount\tnearest_target_id\tnearest_target_seq\tnearest_distance\tedit_class\treverse_complement_target_id\toffset_shift_hint\tadapter_hint\tlow_quality_hint\treason\n");
    for (size_t i = 0; i < row_count && i < top_n; ++i) {
        int nearest_d = -1;
        int nearest = nearest_target(&targets, rows[i].seq, strlen(rows[i].seq), &nearest_d);
        const char *nearest_id = nearest >= 0 ? targets.records[nearest].id : "";
        const char *nearest_seq = nearest >= 0 ? targets.records[nearest].seq : "";
        const char *edit_class = nearest_d == 0 ? "exact" : (nearest_d == 1 ? "one_edit" : "other");
        char *rc = reverse_complement(rows[i].seq, strlen(rows[i].seq));
        int rc_d = -1;
        int rc_nearest = rc == NULL ? -1 : nearest_target(&targets, rc, strlen(rc), &rc_d);
        const char *rc_id = (rc_nearest >= 0 && rc_d <= k) ? targets.records[rc_nearest].id : "";
        int offset_match = 0;
        if (offset_window > 0) {
            for (size_t off = 1; off <= offset_window; ++off) {
                if (target_start + off + strlen(rows[i].seq) <= strlen(rows[i].seq) + target_start) {
                    offset_match = 1;
                    break;
                }
            }
        }
        int adapter_match = adapter != NULL && prefix_matches(rows[i].seq, strlen(rows[i].seq), adapter, strlen(adapter));
        int low_qual = 0;
        int contains_n = contains_char(rows[i].seq, 'N');
        const char *reason = reason_label_for_unmatched(nearest_d, k, rc_nearest >= 0 && rc_d <= k,
                                                        offset_match, adapter_match, low_qual,
                                                        contains_n, 0);
        fprintf(out, "%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n",
                rows[i].seq,
                rows[i].count,
                nearest_id,
                nearest_seq,
                nearest_d,
                edit_class,
                rc_id,
                offset_match ? "possible" : "",
                adapter_match ? "possible" : "",
                low_qual ? "possible" : "",
                reason);
        free(rc);
    }
    fclose(out);
    for (size_t i = 0; i < row_count; ++i) free(rows[i].seq);
    free(rows);
    free_table(&targets);
    return 0;
}

static int run_audit(const char *argv0, int argc, char **argv) {
    (void)argv0;
    const char *targets_path = NULL;
    const char *out_dir = NULL;
    const char *audit_mode = "auto";
    int k = 1;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--targets") == 0 && i + 1 < argc) {
            targets_path = argv[++i];
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0) return 2;
        } else if (strcmp(argv[i], "--out-dir") == 0 && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (strcmp(argv[i], "--audit-mode") == 0 && i + 1 < argc) {
            audit_mode = argv[++i];
        } else {
            return 2;
        }
    }
    if (targets_path == NULL || out_dir == NULL || k < 0) return 2;

    seq_table targets = {0};
    if (read_target_table(targets_path, &targets) != 0) {
        free_table(&targets);
        return 1;
    }
    if (ensure_dir_exists(out_dir) != 0) {
        free_table(&targets);
        return 1;
    }

    qdaln_audit_result audit = {0};
    int use_fast = strcmp(audit_mode, "fast") == 0 || (strcmp(audit_mode, "auto") == 0 && targets.count > 2000);
    if (use_fast) {
        if (qdaln_audit_targets_fast((const char *const *)(&targets.records[0].seq), NULL, targets.count, &audit) != 0) {
            free_table(&targets);
            return 1;
        }
    } else {
        if (qdaln_audit_targets((const char *const *)(&targets.records[0].seq), NULL, targets.count, &audit) != 0) {
            free_table(&targets);
            return 1;
        }
    }

    int *k0_safe = (int *)calloc(targets.count, sizeof(int));
    int *k1_safe = (int *)calloc(targets.count, sizeof(int));
    int *k2_safe = (int *)calloc(targets.count, sizeof(int));
    int *nearest = (int *)calloc(targets.count, sizeof(int));
    if (k0_safe == NULL || k1_safe == NULL || k2_safe == NULL || nearest == NULL) {
        free(k0_safe); free(k1_safe); free(k2_safe); free(nearest);
        qdaln_audit_result_free(&audit);
        free_table(&targets);
        return 1;
    }
    for (size_t i = 0; i < targets.count; ++i) {
        k0_safe[i] = 1;
        k1_safe[i] = 1;
        k2_safe[i] = 1;
        nearest[i] = -1;
    }
    for (size_t i = 0; i < audit.pair_count; ++i) {
        qdaln_audit_pair pair = audit.pairs[i];
        if (pair.distance <= 0) {
            k0_safe[pair.left] = 0;
            k0_safe[pair.right] = 0;
        }
        if (pair.distance <= 1) {
            k1_safe[pair.left] = 0;
            k1_safe[pair.right] = 0;
        }
        if (pair.distance <= 2) {
            k2_safe[pair.left] = 0;
            k2_safe[pair.right] = 0;
        }
        if (nearest[pair.left] < 0 || pair.distance < nearest[pair.left]) nearest[pair.left] = pair.distance;
        if (nearest[pair.right] < 0 || pair.distance < nearest[pair.right]) nearest[pair.right] = pair.distance;
    }

    size_t pairs_at_0 = 0, pairs_at_1 = 0, pairs_at_2 = 0;
    for (size_t i = 0; i < audit.pair_count; ++i) {
        if (audit.pairs[i].distance == 0) ++pairs_at_0;
        if (audit.pairs[i].distance == 1) ++pairs_at_1;
        if (audit.pairs[i].distance == 2) ++pairs_at_2;
    }
    int library_safe_k0 = pairs_at_0 == 0;
    int library_safe_k1 = pairs_at_0 == 0 && pairs_at_1 == 0;
    int library_safe_k2 = pairs_at_0 == 0 && pairs_at_1 == 0 && pairs_at_2 == 0;

    int rc = 0;
    rc |= write_collision_pairs(out_dir, &targets, audit.pairs, audit.pair_count);
    rc |= write_collision_clusters(out_dir, &targets, audit.clusters, audit.cluster_count);
    rc |= write_target_safety(out_dir, &targets, k0_safe, k1_safe, k2_safe, nearest);
    rc |= write_ambiguous_variants(out_dir, &targets, audit.variants, audit.variant_count);
    rc |= write_audit_json_summary(out_dir, k, targets.count, pairs_at_0, pairs_at_1, pairs_at_2,
                                   audit.cluster_count, audit.variant_count,
                                   library_safe_k0, library_safe_k1, library_safe_k2,
                                   use_fast ? "fast" : "exact");

    free(k0_safe); free(k1_safe); free(k2_safe); free(nearest);
    qdaln_audit_result_free(&audit);
    free_table(&targets);
    return rc ? 1 : 0;
}

static qdaln_ambiguity_policy parse_ambiguity_policy(const char *value, int *ok) {
    if (strcmp(value, "best") == 0) {
        *ok = 1;
        return QDALN_POLICY_BEST;
    }
    if (strcmp(value, "radius") == 0) {
        *ok = 1;
        return QDALN_POLICY_RADIUS;
    }
    *ok = 0;
    return QDALN_POLICY_BEST;
}

static qdaln_metric parse_metric_value(const char *value, int *ok) {
    if (strcmp(value, "hamming") == 0) {
        *ok = 1;
        return QDALN_METRIC_HAMMING;
    }
    if (strcmp(value, "levenshtein") == 0) {
        *ok = 1;
        return QDALN_METRIC_LEVENSHTEIN;
    }
    *ok = 0;
    return QDALN_METRIC_HAMMING;
}

static int cmp_long_desc(const void *a, const void *b) {
    const long long ia = *(const long long *)a;
    const long long ib = *(const long long *)b;
    if (ia < ib) return 1;
    if (ia > ib) return -1;
    return 0;
}

static int cmp_pair_count_desc(const void *a, const void *b) {
    const qdaln_pair_count *pa = (const qdaln_pair_count *)a;
    const qdaln_pair_count *pb = (const qdaln_pair_count *)b;
    if (pa->count < pb->count) return 1;
    if (pa->count > pb->count) return -1;
    if (pa->left < pb->left) return -1;
    if (pa->left > pb->left) return 1;
    if (pa->right < pb->right) return -1;
    if (pa->right > pb->right) return 1;
    return 0;
}

static int find_csv_column(char **fields, size_t n, const char *a, const char *b, const char *c, const char *d) {
    for (size_t i = 0; i < n; ++i) {
        if (field_eq(fields[i], a) || (b && field_eq(fields[i], b)) || (c && field_eq(fields[i], c)) ||
            (d && field_eq(fields[i], d))) {
            return (int)i;
        }
    }
    return -1;
}

static int read_bcl_sample_sheet(const char *path, seq_table *samples) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return -1;
    char buf[8192];
    int in_data = 0;
    int sample_col = -1;
    int idx1_col = -1;
    int idx2_col = -1;
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        trim_line(buf);
        if (buf[0] == '\0') continue;
        if (buf[0] == '[') {
            in_data = field_eq(buf, "[data]") || field_eq(buf, "[Data]");
            continue;
        }
        if (!in_data) continue;
        char *fields[32];
        size_t nf = split_fields(buf, ',', fields, 32);
        if (sample_col < 0) {
            sample_col = find_csv_column(fields, nf, "Sample_ID", "SampleID", "sample_id", NULL);
            if (sample_col < 0) sample_col = find_csv_column(fields, nf, "Sample_Name", "SampleName", "sample_name", NULL);
            idx1_col = find_csv_column(fields, nf, "index", "Index", "index1", "i7_index_id");
            idx2_col = find_csv_column(fields, nf, "index2", "Index2", "index_2", "i5_index_id");
            continue;
        }
        if ((size_t)sample_col >= nf || (size_t)idx1_col >= nf) continue;
        char combined[512];
        const char *idx1 = fields[idx1_col];
        const char *idx2 = (idx2_col >= 0 && (size_t)idx2_col < nf) ? fields[idx2_col] : "";
        if (idx2[0] != '\0') {
            snprintf(combined, sizeof(combined), "%s+%s", idx1, idx2);
        } else {
            snprintf(combined, sizeof(combined), "%s", idx1);
        }
        if (push_record(samples, fields[sample_col], strlen(fields[sample_col]), combined, strlen(combined)) != 0) {
            fclose(fp);
            return -1;
        }
    }
    fclose(fp);
    return samples->count == 0 ? -1 : 0;
}

static int read_run_cycles(const char *run_folder, size_t *read1_cycles, size_t *index1_cycles, size_t *index2_cycles) {
    char path[PATH_MAX];
    snprintf(path, sizeof(path), "%s/RunInfo.xml", run_folder);
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return -1;
    char buf[8192];
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        char *read = strstr(buf, "<Read");
        if (read == NULL) continue;
        char *num = strstr(buf, "Number=\"");
        char *cycles = strstr(buf, "NumCycles=\"");
        char *indexed = strstr(buf, "IsIndexedRead=\"");
        if (num == NULL || cycles == NULL || indexed == NULL) continue;
        int number = atoi(num + 8);
        int num_cycles = atoi(cycles + 11);
        int is_indexed = *(indexed + 15) == 'Y';
        if (number == 1 && !is_indexed) *read1_cycles = (size_t)num_cycles;
        if (number == 2 && is_indexed) *index1_cycles = (size_t)num_cycles;
        if (number == 3 && is_indexed) *index2_cycles = (size_t)num_cycles;
    }
    fclose(fp);
    return *read1_cycles == 0 ? -1 : 0;
}

static int parse_lane_dir_name(const char *name) {
    if (strlen(name) != 4 || name[0] != 'L') return -1;
    return atoi(name + 1);
}

static int find_bcl_cycle_path(const char *run_folder, int lane, size_t cycle, char *out_path, size_t out_size) {
    const char *patterns[] = {
        "%s/Data/Intensities/BaseCalls/L%03d/C%zu.1/s_%d_1101.bcl.gz",
        "%s/Data/Intensities/BaseCalls/L%03d/C%zu.1/s_%d_1101.bcl",
        "%s/Data/Intensities/BaseCalls/L%03d/C%zu.1/s_%d.bcl.gz",
        "%s/Data/Intensities/BaseCalls/L%03d/C%zu.1/s_%d.bcl",
    };
    for (size_t i = 0; i < sizeof(patterns) / sizeof(patterns[0]); ++i) {
        snprintf(out_path, out_size, patterns[i], run_folder, lane, cycle, lane);
        struct stat st;
        if (stat(out_path, &st) == 0) return 0;
    }
    return -1;
}

static int find_bcl_filter_path(const char *run_folder, int lane, char *out_path, size_t out_size) {
    const char *patterns[] = {
        "%s/Data/Intensities/BaseCalls/L%03d/s_%d.filter",
        "%s/Data/Intensities/BaseCalls/L%03d/s_%d_1101.filter",
    };
    for (size_t i = 0; i < sizeof(patterns) / sizeof(patterns[0]); ++i) {
        snprintf(out_path, out_size, patterns[i], run_folder, lane, lane);
        struct stat st;
        if (stat(out_path, &st) == 0) return 0;
    }
    return -1;
}

static int decode_bcl_base(uint8_t byte, uint8_t *base, uint8_t *qual) {
    *base = byte & 0x03;
    *qual = byte >> 2;
    return 0;
}

static char bcl_base_char(uint8_t code) {
    switch (code) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

static int read_filter_flags(const char *path, unsigned char **flags_out, uint32_t *count_out) {
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    uint32_t zero = 0, count = 0, version = 0;
    if (fread(&zero, sizeof(uint32_t), 1, fp) != 1 ||
        fread(&version, sizeof(uint32_t), 1, fp) != 1 ||
        fread(&count, sizeof(uint32_t), 1, fp) != 1) {
        fclose(fp);
        return -1;
    }
    unsigned char *flags = (unsigned char *)malloc(count);
    if (flags == NULL) {
        fclose(fp);
        return -1;
    }
    if (fread(flags, 1, count, fp) != count) {
        free(flags);
        fclose(fp);
        return -1;
    }
    fclose(fp);
    *flags_out = flags;
    *count_out = count;
    return 0;
}

static int read_bcl_file_bytes(const char *path, unsigned char **data_out, uint32_t *count_out) {
    size_t len = strlen(path);
    if (len >= 3 && strcmp(path + len - 3, ".gz") == 0) {
        gzFile fp = gzopen(path, "rb");
        if (fp == NULL) return -1;
        uint32_t count = 0;
        if (gzread(fp, &count, sizeof(uint32_t)) != (int)sizeof(uint32_t)) {
            gzclose(fp);
            return -1;
        }
        unsigned char *data = (unsigned char *)malloc(count);
        if (data == NULL) {
            gzclose(fp);
            return -1;
        }
        if (gzread(fp, data, count) != (int)count) {
            free(data);
            gzclose(fp);
            return -1;
        }
        gzclose(fp);
        *data_out = data;
        *count_out = count;
        return 0;
    }
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    uint32_t count = 0;
    if (fread(&count, sizeof(uint32_t), 1, fp) != 1) {
        fclose(fp);
        return -1;
    }
    unsigned char *data = (unsigned char *)malloc(count);
    if (data == NULL) {
        fclose(fp);
        return -1;
    }
    if (fread(data, 1, count, fp) != count) {
        free(data);
        fclose(fp);
        return -1;
    }
    fclose(fp);
    *data_out = data;
    *count_out = count;
    return 0;
}

static int qdaln_hamming_mismatches(const char *a, size_t n, const char *b, size_t m) {
    if (n != m) return -1;
    int mismatches = 0;
    for (size_t i = 0; i < n; ++i) if (a[i] != b[i]) ++mismatches;
    return mismatches;
}

static int match_sample_index(const seq_table *samples, const char *index_seq, int mismatch_limit, int *ambiguous) {
    int best = -1;
    int best_dist = mismatch_limit + 1;
    int best_count = 0;
    *ambiguous = 0;
    for (size_t i = 0; i < samples->count; ++i) {
        int d = qdaln_hamming_mismatches(index_seq, strlen(index_seq), samples->records[i].seq, samples->records[i].len);
        if (d < 0 || d > mismatch_limit) continue;
        if (d < best_dist) {
            best = (int)i;
            best_dist = d;
            best_count = 1;
            *ambiguous = 0;
        } else if (d == best_dist) {
            ++best_count;
            *ambiguous = 1;
        }
    }
    if (best_count != 1) return -1;
    return best;
}

static int write_fastq_record(FILE *out, const char *id, const char *seq, const char *qual) {
    fprintf(out, "@%s\n%s\n+\n%s\n", id, seq, qual);
    return 0;
}

static int write_gz_fastq_record(gzFile out, const char *id, const char *seq, const char *qual) {
    gzprintf(out, "@%s\n%s\n+\n%s\n", id, seq, qual);
    return 0;
}

static int ensure_sample_out(gzFile *handles, const char *out_dir, const char *sample_id) {
    (void)handles;
    (void)out_dir;
    (void)sample_id;
    return 0;
}

static int run_bcl_demux(const char *argv0, int argc, char **argv) {
    (void)argv0;
    const char *run_folder = NULL;
    const char *sample_sheet = NULL;
    const char *out_dir = NULL;
    const char *summary_path = NULL;
    int barcode_mismatches = 1;
    int threads = 1;
    int gzip_level = 1;
    int emit_index_fastqs = 0;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--run-folder") == 0 && i + 1 < argc) {
            run_folder = argv[++i];
        } else if (strcmp(argv[i], "--sample-sheet") == 0 && i + 1 < argc) {
            sample_sheet = argv[++i];
        } else if (strcmp(argv[i], "--out-dir") == 0 && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (strcmp(argv[i], "--barcode-mismatches") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "1,1") == 0) {
                barcode_mismatches = 1;
            } else if (parse_int_value(value, &barcode_mismatches) != 0) {
                return 2;
            }
        } else if (strcmp(argv[i], "--summary") == 0 && i + 1 < argc) {
            summary_path = argv[++i];
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &threads) != 0) return 2;
        } else if (strcmp(argv[i], "--gzip-level") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &gzip_level) != 0) return 2;
        } else if (strcmp(argv[i], "--emit-index-fastqs") == 0) {
            emit_index_fastqs = 1;
        } else {
            return 2;
        }
    }
    if (run_folder == NULL || sample_sheet == NULL || out_dir == NULL) return 2;
    if (threads <= 0 || gzip_level < 0 || gzip_level > 9) return 2;

    if (ensure_dir_exists(out_dir) != 0) return 1;

    seq_table samples = {0};
    if (read_bcl_sample_sheet(sample_sheet, &samples) != 0) {
        free_table(&samples);
        return 1;
    }

    size_t read1_cycles = 0, index1_cycles = 0, index2_cycles = 0;
    if (read_run_cycles(run_folder, &read1_cycles, &index1_cycles, &index2_cycles) != 0) {
        free_table(&samples);
        return 1;
    }

    DIR *basecalls = opendir((snprintf((char[PATH_MAX]){0}, PATH_MAX, "%s/Data/Intensities/BaseCalls", run_folder), (char[PATH_MAX]){0}));
    (void)basecalls;

    char basecalls_path[PATH_MAX];
    snprintf(basecalls_path, sizeof(basecalls_path), "%s/Data/Intensities/BaseCalls", run_folder);
    DIR *dir = opendir(basecalls_path);
    if (dir == NULL) {
        free_table(&samples);
        return 1;
    }

    long long total_reads = 0;
    long long assigned_reads = 0;
    long long undetermined_reads = 0;

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        int lane = parse_lane_dir_name(entry->d_name);
        if (lane < 0) continue;

        unsigned char *filter_flags = NULL;
        uint32_t filter_count = 0;
        char filter_path[PATH_MAX];
        if (find_bcl_filter_path(run_folder, lane, filter_path, sizeof(filter_path)) != 0 ||
            read_filter_flags(filter_path, &filter_flags, &filter_count) != 0) {
            free(filter_flags);
            closedir(dir);
            free_table(&samples);
            return 1;
        }

        unsigned char **read1 = (unsigned char **)calloc(read1_cycles, sizeof(unsigned char *));
        unsigned char **index1 = (unsigned char **)calloc(index1_cycles, sizeof(unsigned char *));
        unsigned char **index2 = (unsigned char **)calloc(index2_cycles, sizeof(unsigned char *));
        uint32_t count = 0;
        if ((read1_cycles && read1 == NULL) || (index1_cycles && index1 == NULL) || (index2_cycles && index2 == NULL)) {
            free(filter_flags);
            free(read1);
            free(index1);
            free(index2);
            closedir(dir);
            free_table(&samples);
            return 1;
        }

        int lane_ok = 1;
        for (size_t cycle = 0; cycle < read1_cycles; ++cycle) {
            char path[PATH_MAX];
            if (find_bcl_cycle_path(run_folder, lane, cycle + 1, path, sizeof(path)) != 0 ||
                read_bcl_file_bytes(path, &read1[cycle], &count) != 0 || count != filter_count) {
                lane_ok = 0;
                break;
            }
        }
        for (size_t cycle = 0; lane_ok && cycle < index1_cycles; ++cycle) {
            char path[PATH_MAX];
            if (find_bcl_cycle_path(run_folder, lane, read1_cycles + cycle + 1, path, sizeof(path)) != 0 ||
                read_bcl_file_bytes(path, &index1[cycle], &count) != 0 || count != filter_count) {
                lane_ok = 0;
                break;
            }
        }
        for (size_t cycle = 0; lane_ok && cycle < index2_cycles; ++cycle) {
            char path[PATH_MAX];
            if (find_bcl_cycle_path(run_folder, lane, read1_cycles + index1_cycles + cycle + 1, path, sizeof(path)) != 0 ||
                read_bcl_file_bytes(path, &index2[cycle], &count) != 0 || count != filter_count) {
                lane_ok = 0;
                break;
            }
        }
        if (!lane_ok) {
            for (size_t i = 0; i < read1_cycles; ++i) free(read1[i]);
            for (size_t i = 0; i < index1_cycles; ++i) free(index1[i]);
            for (size_t i = 0; i < index2_cycles; ++i) free(index2[i]);
            free(read1); free(index1); free(index2); free(filter_flags);
            closedir(dir);
            free_table(&samples);
            return 1;
        }

        gzFile *sample_handles = (gzFile *)calloc(samples.count, sizeof(gzFile));
        char lane_stats_path[PATH_MAX];
        snprintf(lane_stats_path, sizeof(lane_stats_path), "%s/L%03d.fastq.gz", out_dir, lane);
        gzFile undetermined = NULL;
        char undetermined_path[PATH_MAX];
        snprintf(undetermined_path, sizeof(undetermined_path), "%s/Undetermined_L%03d.fastq.gz", out_dir, lane);
        undetermined = gzopen(undetermined_path, "wb");
        if (undetermined == NULL) {
            for (size_t i = 0; i < read1_cycles; ++i) free(read1[i]);
            for (size_t i = 0; i < index1_cycles; ++i) free(index1[i]);
            for (size_t i = 0; i < index2_cycles; ++i) free(index2[i]);
            free(read1); free(index1); free(index2); free(filter_flags); free(sample_handles);
            closedir(dir);
            free_table(&samples);
            return 1;
        }
        if (gzip_level != 6) gzsetparams(undetermined, gzip_level, Z_DEFAULT_STRATEGY);

        gzFile *index_fastqs = NULL;
        if (emit_index_fastqs) {
            index_fastqs = (gzFile *)calloc(samples.count + 1, sizeof(gzFile));
        }

        for (uint32_t cluster = 0; cluster < filter_count; ++cluster) {
            if (filter_flags[cluster] == 0) continue;
            char read_seq[1024];
            char read_qual[1024];
            char index_seq[1024];
            size_t read_len = read1_cycles < sizeof(read_seq) - 1 ? read1_cycles : sizeof(read_seq) - 1;
            size_t index_len = index1_cycles + index2_cycles;
            if (index_len >= sizeof(index_seq)) index_len = sizeof(index_seq) - 1;
            for (size_t cycle = 0; cycle < read_len; ++cycle) {
                uint8_t base = 0, qual = 0;
                decode_bcl_base(read1[cycle][cluster], &base, &qual);
                read_seq[cycle] = bcl_base_char(base);
                read_qual[cycle] = (char)(33 + (qual > 40 ? 40 : qual));
            }
            read_seq[read_len] = '\0';
            read_qual[read_len] = '\0';
            size_t pos = 0;
            for (size_t cycle = 0; cycle < index1_cycles && pos < index_len; ++cycle) {
                uint8_t base = 0, qual = 0;
                decode_bcl_base(index1[cycle][cluster], &base, &qual);
                index_seq[pos++] = bcl_base_char(base);
            }
            if (index2_cycles > 0 && pos < index_len) index_seq[pos++] = '+';
            for (size_t cycle = 0; cycle < index2_cycles && pos < index_len; ++cycle) {
                uint8_t base = 0, qual = 0;
                decode_bcl_base(index2[cycle][cluster], &base, &qual);
                index_seq[pos++] = bcl_base_char(base);
            }
            index_seq[pos] = '\0';

            int ambiguous = 0;
            int sample_idx = match_sample_index(&samples, index_seq, barcode_mismatches, &ambiguous);
            char id[128];
            snprintf(id, sizeof(id), "%d:%u", lane, cluster + 1);
            ++total_reads;
            if (sample_idx >= 0 && !ambiguous) {
                if (sample_handles[sample_idx] == NULL) {
                    char path[PATH_MAX];
                    snprintf(path, sizeof(path), "%s/%s_L%03d.fastq.gz", out_dir, samples.records[sample_idx].id, lane);
                    sample_handles[sample_idx] = gzopen(path, "ab");
                    if (sample_handles[sample_idx] == NULL) {
                        gzclose(undetermined);
                        for (size_t i = 0; i < samples.count; ++i) if (sample_handles[i] != NULL) gzclose(sample_handles[i]);
                        if (index_fastqs != NULL) {
                            for (size_t i = 0; i < samples.count + 1; ++i) if (index_fastqs[i] != NULL) gzclose(index_fastqs[i]);
                            free(index_fastqs);
                        }
                        for (size_t i = 0; i < read1_cycles; ++i) free(read1[i]);
                        for (size_t i = 0; i < index1_cycles; ++i) free(index1[i]);
                        for (size_t i = 0; i < index2_cycles; ++i) free(index2[i]);
                        free(read1); free(index1); free(index2); free(filter_flags); free(sample_handles);
                        closedir(dir);
                        free_table(&samples);
                        return 1;
                    }
                    if (gzip_level != 6) gzsetparams(sample_handles[sample_idx], gzip_level, Z_DEFAULT_STRATEGY);
                }
                write_gz_fastq_record(sample_handles[sample_idx], id, read_seq, read_qual);
                if (emit_index_fastqs) {
                    if (index_fastqs[sample_idx] == NULL) {
                        char path[PATH_MAX];
                        snprintf(path, sizeof(path), "%s/%s_index_L%03d.fastq.gz", out_dir, samples.records[sample_idx].id, lane);
                        index_fastqs[sample_idx] = gzopen(path, "ab");
                        if (index_fastqs[sample_idx] != NULL && gzip_level != 6) {
                            gzsetparams(index_fastqs[sample_idx], gzip_level, Z_DEFAULT_STRATEGY);
                        }
                    }
                    if (index_fastqs[sample_idx] != NULL) {
                        write_gz_fastq_record(index_fastqs[sample_idx], id, index_seq, read_qual);
                    }
                }
                ++assigned_reads;
            } else {
                write_gz_fastq_record(undetermined, id, read_seq, read_qual);
                if (emit_index_fastqs) {
                    size_t slot = samples.count;
                    if (index_fastqs[slot] == NULL) {
                        char path[PATH_MAX];
                        snprintf(path, sizeof(path), "%s/Undetermined_index_L%03d.fastq.gz", out_dir, lane);
                        index_fastqs[slot] = gzopen(path, "ab");
                        if (index_fastqs[slot] != NULL && gzip_level != 6) {
                            gzsetparams(index_fastqs[slot], gzip_level, Z_DEFAULT_STRATEGY);
                        }
                    }
                    if (index_fastqs[slot] != NULL) {
                        write_gz_fastq_record(index_fastqs[slot], id, index_seq, read_qual);
                    }
                }
                ++undetermined_reads;
            }
        }

        gzclose(undetermined);
        for (size_t i = 0; i < samples.count; ++i) if (sample_handles[i] != NULL) gzclose(sample_handles[i]);
        if (index_fastqs != NULL) {
            for (size_t i = 0; i < samples.count + 1; ++i) if (index_fastqs[i] != NULL) gzclose(index_fastqs[i]);
            free(index_fastqs);
        }
        free(sample_handles);
        for (size_t i = 0; i < read1_cycles; ++i) free(read1[i]);
        for (size_t i = 0; i < index1_cycles; ++i) free(index1[i]);
        for (size_t i = 0; i < index2_cycles; ++i) free(index2[i]);
        free(read1); free(index1); free(index2); free(filter_flags);
    }
    closedir(dir);

    char stats_path[PATH_MAX];
    snprintf(stats_path, sizeof(stats_path), "%s/Demultiplex_Stats.csv", out_dir);
    FILE *stats = open_output_file(stats_path);
    if (stats == NULL) {
        free_table(&samples);
        return 1;
    }
    fprintf(stats, "SampleID,Reads\n");
    for (size_t i = 0; i < samples.count; ++i) {
        fprintf(stats, "%s,%d\n", samples.records[i].id, 0);
    }
    fprintf(stats, "Undetermined,%lld\n", undetermined_reads);
    fclose(stats);

    char normalized_sheet[PATH_MAX];
    snprintf(normalized_sheet, sizeof(normalized_sheet), "%s/SampleSheet.normalized.csv", out_dir);
    FILE *norm = open_output_file(normalized_sheet);
    if (norm != NULL) {
        fprintf(norm, "Sample_ID,index\n");
        for (size_t i = 0; i < samples.count; ++i) fprintf(norm, "%s,%s\n", samples.records[i].id, samples.records[i].seq);
        fclose(norm);
    }

    if (summary_path != NULL) {
        FILE *summary = open_output_file(summary_path);
        if (summary == NULL) {
            free_table(&samples);
            return 1;
        }
        fprintf(summary,
                "{\n  \"total_reads\": %lld,\n  \"assigned_reads\": %lld,\n  \"undetermined_reads\": %lld,\n"
                "  \"threads\": %d,\n  \"gzip_level\": %d,\n  \"emit_index_fastqs\": %s\n}\n",
                total_reads, assigned_reads, undetermined_reads,
                threads, gzip_level, emit_index_fastqs ? "true" : "false");
        fclose(summary);
    }

    free_table(&samples);
    return 0;
}

static int compare_files_exact(const char *left, const char *right, int *same) {
    FILE *a = fopen(left, "rb");
    FILE *b = fopen(right, "rb");
    if (a == NULL || b == NULL) {
        if (a) fclose(a);
        if (b) fclose(b);
        return -1;
    }
    int rc = 0;
    while (1) {
        int ca = fgetc(a);
        int cb = fgetc(b);
        if (ca != cb) {
            *same = 0;
            break;
        }
        if (ca == EOF) {
            *same = 1;
            break;
        }
    }
    fclose(a);
    fclose(b);
    return rc;
}

static int run_bcl_validate(const char *argv0, int argc, char **argv) {
    (void)argv0;
    const char *dotmatch_out = NULL;
    const char *truth_out = NULL;
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--dotmatch-out") == 0 && i + 1 < argc) {
            dotmatch_out = argv[++i];
        } else if (strcmp(argv[i], "--truth-out") == 0 && i + 1 < argc) {
            truth_out = argv[++i];
        } else {
            return 2;
        }
    }
    if (dotmatch_out == NULL || truth_out == NULL) return 2;

    char left[PATH_MAX];
    char right[PATH_MAX];
    snprintf(left, sizeof(left), "%s/Demultiplex_Stats.csv", dotmatch_out);
    snprintf(right, sizeof(right), "%s/Demultiplex_Stats.csv", truth_out);
    int same = 0;
    if (compare_files_exact(left, right, &same) != 0 || !same) return 1;
    return 0;
}

static int run_count(const char *argv0, int argc, char **argv) {
    (void)argv0;
    const char *targets_path = NULL;
    const char *out_path = NULL;
    const char *target_counts_long_path = NULL;
    const char *sample_qc_path = NULL;
    const char *assignments_path = NULL;
    const char *summary_path = NULL;
    const char *report_path = NULL;
    const char *report_audit_dir = NULL;
    const char *report_unmatched = NULL;
    const char *format = "dotmatch";
    const char *hamming_index_mode = "auto";
    const char **reads = NULL;
    size_t reads_count = 0;
    char **sample_labels = NULL;
    size_t target_start = 0;
    size_t target_length = 0;
    int k = 0;
    qdaln_metric metric = QDALN_METRIC_HAMMING;
    qdaln_ambiguity_policy policy = QDALN_POLICY_BEST;
    qdaln_offset_mode offset_mode = QDALN_OFFSET_BEST;
    size_t offset_window = 0;
    int max_correction_qual = -1;
    int indel_window = 0;
    int threads = 1;
    int ambiguity_report = 0;
    int hamming_index_override = -1;

    if (parse_reads_option(argc, argv, &reads, &reads_count) != 0) return 2;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--targets") == 0 && i + 1 < argc) {
            targets_path = argv[++i];
        } else if (strcmp(argv[i], "--reads") == 0 && i + 1 < argc) {
            ++i;
        } else if (strcmp(argv[i], "--target-start") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_start) != 0) { free((void *)reads); return 2; }
        } else if (strcmp(argv[i], "--target-length") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_length) != 0) { free((void *)reads); return 2; }
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0) { free((void *)reads); return 2; }
        } else if (strcmp(argv[i], "--metric") == 0 && i + 1 < argc) {
            int ok = 0;
            metric = parse_metric_value(argv[++i], &ok);
            if (!ok) { free((void *)reads); return 2; }
        } else if (strcmp(argv[i], "--ambiguity-policy") == 0 && i + 1 < argc) {
            int ok = 0;
            policy = parse_ambiguity_policy(argv[++i], &ok);
            if (!ok) { free((void *)reads); return 2; }
        } else if (strcmp(argv[i], "--offset-mode") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "best") == 0) {
                offset_mode = QDALN_OFFSET_BEST;
            } else if (strcmp(value, "multi") == 0) {
                offset_mode = QDALN_OFFSET_MULTI;
            } else {
                free((void *)reads);
                return 2;
            }
        } else if (strcmp(argv[i], "--auto-offset") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &offset_window) != 0 || offset_window > MAX_AUTO_OFFSET) {
                free((void *)reads);
                return 2;
            }
        } else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else if (strcmp(argv[i], "--target-counts-long") == 0 && i + 1 < argc) {
            target_counts_long_path = argv[++i];
        } else if (strcmp(argv[i], "--sample-qc") == 0 && i + 1 < argc) {
            sample_qc_path = argv[++i];
        } else if (strcmp(argv[i], "--assignments") == 0 && i + 1 < argc) {
            assignments_path = argv[++i];
        } else if (strcmp(argv[i], "--summary") == 0 && i + 1 < argc) {
            summary_path = argv[++i];
        } else if (strcmp(argv[i], "--report") == 0 && i + 1 < argc) {
            report_path = argv[++i];
        } else if (strcmp(argv[i], "--report-audit-dir") == 0 && i + 1 < argc) {
            report_audit_dir = argv[++i];
        } else if (strcmp(argv[i], "--report-unmatched") == 0 && i + 1 < argc) {
            report_unmatched = argv[++i];
        } else if (strcmp(argv[i], "--format") == 0 && i + 1 < argc) {
            format = argv[++i];
        } else if (strcmp(argv[i], "--hamming-index") == 0 && i + 1 < argc) {
            hamming_index_mode = argv[++i];
        } else if (strcmp(argv[i], "--max-correction-qual") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &max_correction_qual) != 0) { free((void *)reads); return 2; }
        } else if (strcmp(argv[i], "--indel-window") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &indel_window) != 0) { free((void *)reads); return 2; }
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &threads) != 0) { free((void *)reads); return 2; }
        } else if (strcmp(argv[i], "--ambiguous") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "report") == 0) {
                ambiguity_report = 1;
            } else if (strcmp(value, "discard") == 0) {
                ambiguity_report = 0;
            } else {
                free((void *)reads);
                return 2;
            }
        }
    }

    if (targets_path == NULL || out_path == NULL || target_length == 0 || reads_count == 0 || k < 0) {
        free((void *)reads);
        return 2;
    }
    if (threads <= 0 || max_correction_qual > 93 || (metric == QDALN_METRIC_HAMMING && k > 1) || indel_window < 0 || indel_window > 1) {
        free((void *)reads);
        return 2;
    }
    if (metric == QDALN_METRIC_HAMMING && indel_window != 0) {
        free((void *)reads);
        return 2;
    }
    if (k > 1 && indel_window != 0) {
        free((void *)reads);
        return 2;
    }
    if (strcmp(hamming_index_mode, "auto") == 0) {
        hamming_index_override = -1;
    } else if (strcmp(hamming_index_mode, "query") == 0) {
        hamming_index_override = 0;
    } else if (strcmp(hamming_index_mode, "precompute") == 0) {
        hamming_index_override = 1;
    } else {
        free((void *)reads);
        return 2;
    }

    seq_table targets = {0};
    if (read_target_table(targets_path, &targets) != 0) {
        free((void *)reads);
        free_table(&targets);
        return 1;
    }

    if (parse_sample_labels(argc, argv, reads_count, &sample_labels) != 0) {
        free((void *)reads);
        free_table(&targets);
        return 1;
    }

    const char **target_ptrs = (const char **)malloc(targets.count * sizeof(char *));
    size_t *target_lens = (size_t *)malloc(targets.count * sizeof(size_t));
    if ((targets.count != 0) && (target_ptrs == NULL || target_lens == NULL)) {
        free((void *)reads);
        free(target_ptrs);
        free(target_lens);
        free_sample_labels(sample_labels, reads_count);
        free_table(&targets);
        return 1;
    }
    for (size_t i = 0; i < targets.count; ++i) {
        target_ptrs[i] = targets.records[i].seq;
        target_lens[i] = targets.records[i].len;
    }

    qdaln_index *index = qdaln_index_build_metric(target_ptrs, target_lens, targets.count, metric, hamming_index_override);
    if (index == NULL) {
        free((void *)reads);
        free(target_ptrs);
        free(target_lens);
        free_sample_labels(sample_labels, reads_count);
        free_table(&targets);
        return 1;
    }

    int *amb_flags = (int *)calloc(targets.count, sizeof(int));
    if (amb_flags == NULL || target_ambiguity_flags(&targets, k, amb_flags) != 0) {
        free((void *)reads);
        free(target_ptrs);
        free(target_lens);
        free(amb_flags);
        qdaln_index_free(index);
        free_sample_labels(sample_labels, reads_count);
        free_table(&targets);
        return 1;
    }

    FILE *assignments = NULL;
    if (assignments_path != NULL) {
        assignments = open_output_file(assignments_path);
        if (assignments == NULL) {
            free((void *)reads);
            free(target_ptrs);
            free(target_lens);
            free(amb_flags);
            qdaln_index_free(index);
            free_sample_labels(sample_labels, reads_count);
            free_table(&targets);
            return 1;
        }
        fprintf(assignments, "read_id\tobserved_seq\ttarget_id\ttarget_seq\tdistance\tstatus\tmatch_count\tsecond_best_distance\tcorrection\n");
    }

    FILE *target_counts_long = NULL;
    if (target_counts_long_path != NULL) {
        target_counts_long = open_output_file(target_counts_long_path);
        if (target_counts_long == NULL) {
            if (assignments) fclose(assignments);
            free((void *)reads);
            free(target_ptrs);
            free(target_lens);
            free(amb_flags);
            qdaln_index_free(index);
            free_sample_labels(sample_labels, reads_count);
            free_table(&targets);
            return 1;
        }
        fprintf(target_counts_long, "sample_label\ttarget_id\ttarget_seq\tgene\tcount_exact\tcount_corrected_substitution\tcount_corrected_insertion\tcount_corrected_deletion\tcount_corrected_other\tcount_total\tambiguous_nearby\n");
    }

    FILE *sample_qc = NULL;
    if (sample_qc_path != NULL) {
        sample_qc = open_output_file(sample_qc_path);
        if (sample_qc == NULL) {
            if (assignments) fclose(assignments);
            if (target_counts_long) fclose(target_counts_long);
            free((void *)reads);
            free(target_ptrs);
            free(target_lens);
            free(amb_flags);
            qdaln_index_free(index);
            free_sample_labels(sample_labels, reads_count);
            free_table(&targets);
            return 1;
        }
        fprintf(sample_qc, "sample_label\ttotal_reads\tassigned_unique\tassignment_rate\texact_rate\trescued_rate\tambiguous_rate\tnone_rate\ttarget_coverage\tzero_count_targets\tgini_index\ttop_1_percent_dominance\tcandidates_considered\tcandidates_verified\tauto_offset_rescues\n");
    }

    FILE *summary = NULL;
    if (summary_path != NULL) {
        summary = open_output_file(summary_path);
        if (summary == NULL) {
            if (assignments) fclose(assignments);
            if (target_counts_long) fclose(target_counts_long);
            if (sample_qc) fclose(sample_qc);
            free((void *)reads);
            free(target_ptrs);
            free(target_lens);
            free(amb_flags);
            qdaln_index_free(index);
            free_sample_labels(sample_labels, reads_count);
            free_table(&targets);
            return 1;
        }
        fprintf(summary, "[\n");
    }

    int *global_counts = (int *)calloc(reads_count * targets.count, sizeof(int));
    if (global_counts == NULL) {
        if (assignments) fclose(assignments);
        if (target_counts_long) fclose(target_counts_long);
        if (sample_qc) fclose(sample_qc);
        if (summary) fclose(summary);
        free((void *)reads);
        free(target_ptrs);
        free(target_lens);
        free(amb_flags);
        qdaln_index_free(index);
        free_sample_labels(sample_labels, reads_count);
        free_table(&targets);
        return 1;
    }

    int rc = 0;
    for (size_t s = 0; s < reads_count && rc == 0; ++s) {
        int *exact = (int *)calloc(targets.count, sizeof(int));
        int *sub = (int *)calloc(targets.count, sizeof(int));
        int *ins = (int *)calloc(targets.count, sizeof(int));
        int *del = (int *)calloc(targets.count, sizeof(int));
        int *other = (int *)calloc(targets.count, sizeof(int));
        int *totals = (int *)calloc(targets.count, sizeof(int));
        if (exact == NULL || sub == NULL || ins == NULL || del == NULL || other == NULL || totals == NULL) {
            free(exact); free(sub); free(ins); free(del); free(other); free(totals);
            rc = 1;
            break;
        }
        FILE *text_fp = NULL;
        gzFile gz_fp = NULL;
        int is_gz = 0;
        long long total_reads = 0, unique_reads = 0, exact_reads = 0, rescued_reads = 0, ambiguous_reads = 0,
                  none_reads = 0, invalid_reads = 0, candidates_considered = 0, candidates_verified = 0,
                  offset_rescues = 0;
        while (rc == 0) {
            char *id = NULL;
            char *seq = NULL;
            char *qual = NULL;
            int ok = parse_text_or_gz_fastq_record(reads[s], &text_fp, &gz_fp, &is_gz, &id, &seq, &qual);
            if (ok == 0) break;
            if (ok < 0) {
                free(id); free(seq); free(qual);
                rc = 1;
                break;
            }
            ++total_reads;
            char *window = NULL;
            int best_offset = (int)target_start;
            int auto_enabled = offset_window > 0;
            if (parse_target_window(seq, strlen(seq), target_start, target_length,
                                    auto_enabled, offset_window, index, metric,
                                    max_correction_qual, qual, &best_offset, &window) != 0) {
                ++invalid_reads;
                if (assignments != NULL) {
                    qdaln_match_result invalid = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
                    write_assign_result(assignments, "count", id, seq, &targets, invalid);
                }
                free(id); free(seq); free(qual);
                continue;
            }
            qdaln_match_result r = {0};
            qdaln_index_stats stats = {0};
            if (qdaln_index_match_read_stats(index, window, target_length, k, policy, metric, &r, &stats) != 0) {
                free(window); free(id); free(seq); free(qual);
                rc = 1;
                break;
            }
            candidates_considered += stats.candidates_considered;
            candidates_verified += stats.candidates_verified;
            if (qual != NULL && max_correction_qual >= 0 && r.status == QDALN_MATCH_UNIQUE &&
                r.target_index >= 0 && r.best_distance == 1) {
                int qok = qdaln_read_correction_quality_ok(
                    window, target_length,
                    targets.records[r.target_index].seq,
                    targets.records[r.target_index].len,
                    qual + best_offset, target_length,
                    max_correction_qual,
                    metric);
                if (qok == 0) {
                    r.target_index = -1;
                    r.best_distance = -1;
                    r.second_best_distance = -1;
                    r.match_count = 0;
                    r.status = QDALN_MATCH_NONE;
                } else if (qok < 0) {
                    free(window); free(id); free(seq); free(qual);
                    rc = 1;
                    break;
                }
            }
            if (r.status == QDALN_MATCH_UNIQUE && best_offset != (int)target_start && r.best_distance == 0) {
                ++offset_rescues;
            }
            if (r.status == QDALN_MATCH_UNIQUE && r.target_index >= 0) {
                ++unique_reads;
                count_provenance_fields(&r, window, target_length, &targets, exact, sub, ins, del, other);
                totals[r.target_index] += 1;
                global_counts[s * targets.count + r.target_index] += 1;
                if (r.best_distance == 0) ++exact_reads; else ++rescued_reads;
            } else if (r.status == QDALN_MATCH_AMBIGUOUS) {
                ++ambiguous_reads;
            } else if (r.status == QDALN_MATCH_NONE) {
                ++none_reads;
            } else {
                ++invalid_reads;
            }
            if (assignments != NULL && (r.status != QDALN_MATCH_AMBIGUOUS || ambiguity_report)) {
                const char *target_id = (r.target_index >= 0) ? targets.records[r.target_index].id : "";
                const char *target_seq = (r.target_index >= 0) ? targets.records[r.target_index].seq : "";
                const char *correction = r.status == QDALN_MATCH_UNIQUE && r.target_index >= 0
                                         ? (r.best_distance == 0 ? "exact"
                                            : (compute_edit_kind(window, target_length, target_seq, strlen(target_seq)) == 1 ? "substitution"
                                               : (compute_edit_kind(window, target_length, target_seq, strlen(target_seq)) == 2 ? "insertion"
                                                  : (compute_edit_kind(window, target_length, target_seq, strlen(target_seq)) == 3 ? "deletion" : "other"))))
                                         : (r.status == QDALN_MATCH_AMBIGUOUS ? "ambiguous" : (r.status == QDALN_MATCH_NONE ? "none" : "invalid"));
                fprintf(assignments, "%s\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\n",
                        id,
                        window,
                        target_id,
                        target_seq,
                        r.best_distance,
                        status_name(r.status),
                        r.match_count,
                        r.second_best_distance,
                        correction);
            }
            free(window); free(id); free(seq); free(qual);
        }
        if (text_fp != NULL) fclose(text_fp);
        if (gz_fp != NULL) gzclose(gz_fp);

        if (rc == 0) {
            if (strcmp(format, "mageck") == 0) {
                /* postpone writing until all samples processed */
            } else {
                FILE *out = (s == 0) ? open_output_file(out_path) : fopen(out_path, "a");
                if (out == NULL) {
                    rc = 1;
                } else {
                    if (s == 0) {
                        write_count_table(out, &targets, exact, sub, ins, del, other, amb_flags);
                    }
                    fclose(out);
                }
            }
            if (target_counts_long != NULL) {
                write_target_counts_long(target_counts_long, sample_labels[s], &targets, exact, sub, ins, del, other, amb_flags);
            }
            if (sample_qc != NULL) {
                write_sample_qc(sample_qc, sample_labels[s], total_reads, unique_reads, exact_reads, rescued_reads,
                                ambiguous_reads, none_reads, targets.count, totals,
                                candidates_considered, candidates_verified, offset_rescues);
            }
            if (summary != NULL) {
                if (s > 0) fprintf(summary, ",\n");
                write_summary_json(summary, total_reads, exact_reads, rescued_reads, ambiguous_reads,
                                   none_reads, invalid_reads, candidates_considered, candidates_verified,
                                   offset_rescues, sample_labels[s]);
            }
            if (report_path != NULL && s == 0) {
                if (report_audit_dir != NULL) write_audit_summary(report_audit_dir, &targets, k);
                write_count_html_report(report_path, &targets, exact, sub, ins, del, other, amb_flags,
                                        sample_labels[s], total_reads, unique_reads, exact_reads, rescued_reads,
                                        ambiguous_reads, none_reads, invalid_reads, candidates_considered,
                                        candidates_verified, offset_rescues, totals,
                                        report_audit_dir, report_unmatched);
            }
        }

        free(exact); free(sub); free(ins); free(del); free(other); free(totals);
    }

    if (rc == 0 && strcmp(format, "mageck") == 0) {
        FILE *out = open_output_file(out_path);
        if (out == NULL) {
            rc = 1;
        } else {
            write_mageck_table(out, &targets, sample_labels, reads_count, global_counts);
            fclose(out);
        }
    }

    if (summary != NULL) {
        fprintf(summary, "\n]\n");
        fclose(summary);
    }
    if (assignments != NULL) fclose(assignments);
    if (target_counts_long != NULL) fclose(target_counts_long);
    if (sample_qc != NULL) fclose(sample_qc);
    free(global_counts);
    free((void *)reads);
    free(target_ptrs);
    free(target_lens);
    free(amb_flags);
    qdaln_index_free(index);
    free_sample_labels(sample_labels, reads_count);
    free_table(&targets);
    return rc;
}

static int run_crispr_count(const char *argv0, int argc, char **argv) {
    const char *library = NULL;
    const char *samples_path = NULL;
    const char *out = NULL;
    const char *summary = NULL;
    size_t guide_start = 0;
    size_t guide_length = 0;
    int k = 0;
    qdaln_metric metric = QDALN_METRIC_HAMMING;
    int indel_window = 0;
    int max_correction_qual = -1;
    int ambiguity_report = 0;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--library") == 0 && i + 1 < argc) {
            library = argv[++i];
        } else if (strcmp(argv[i], "--samples") == 0 && i + 1 < argc) {
            samples_path = argv[++i];
        } else if (strcmp(argv[i], "--guide-start") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &guide_start) != 0) return 2;
        } else if (strcmp(argv[i], "--guide-length") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &guide_length) != 0) return 2;
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0) return 2;
        } else if (strcmp(argv[i], "--metric") == 0 && i + 1 < argc) {
            int ok = 0;
            metric = parse_metric_value(argv[++i], &ok);
            if (!ok) return 2;
        } else if (strcmp(argv[i], "--indel-window") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &indel_window) != 0) return 2;
        } else if (strcmp(argv[i], "--max-correction-qual") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &max_correction_qual) != 0) return 2;
        } else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out = argv[++i];
        } else if (strcmp(argv[i], "--summary") == 0 && i + 1 < argc) {
            summary = argv[++i];
        } else if (strcmp(argv[i], "--ambiguous") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "report") == 0) {
                ambiguity_report = 1;
            } else if (strcmp(value, "discard") == 0) {
                ambiguity_report = 0;
            } else {
                return 2;
            }
        } else {
            return 2;
        }
    }
    if (library == NULL || samples_path == NULL || out == NULL || guide_length == 0) return 2;

    seq_table samples = {0};
    if (load_sample_manifest(samples_path, &samples) != 0) {
        free_table(&samples);
        return 1;
    }

    int argc2 = 0;
    char **argv2 = (char **)calloc(32 + samples.count * 2, sizeof(char *));
    if (argv2 == NULL) {
        free_table(&samples);
        return 1;
    }
    argv2[argc2++] = argv[0];
    argv2[argc2++] = "count";
    argv2[argc2++] = "--targets";
    argv2[argc2++] = (char *)library;
    for (size_t i = 0; i < samples.count; ++i) {
        argv2[argc2++] = "--reads";
        argv2[argc2++] = samples.records[i].seq;
    }
    argv2[argc2++] = "--target-start";
    char guide_start_buf[32];
    snprintf(guide_start_buf, sizeof(guide_start_buf), "%zu", guide_start);
    argv2[argc2++] = guide_start_buf;
    argv2[argc2++] = "--target-length";
    char guide_length_buf[32];
    snprintf(guide_length_buf, sizeof(guide_length_buf), "%zu", guide_length);
    argv2[argc2++] = guide_length_buf;
    argv2[argc2++] = "--k";
    char k_buf[16];
    snprintf(k_buf, sizeof(k_buf), "%d", k);
    argv2[argc2++] = k_buf;
    argv2[argc2++] = "--metric";
    argv2[argc2++] = metric == QDALN_METRIC_HAMMING ? "hamming" : "levenshtein";
    argv2[argc2++] = "--indel-window";
    char indel_buf[16];
    snprintf(indel_buf, sizeof(indel_buf), "%d", indel_window);
    argv2[argc2++] = indel_buf;
    if (max_correction_qual >= 0) {
        argv2[argc2++] = "--max-correction-qual";
        char qual_buf[16];
        snprintf(qual_buf, sizeof(qual_buf), "%d", max_correction_qual);
        argv2[argc2++] = strdup(qual_buf);
    }
    argv2[argc2++] = "--format";
    argv2[argc2++] = "mageck";
    argv2[argc2++] = "--out";
    argv2[argc2++] = (char *)out;
    if (summary != NULL) {
        argv2[argc2++] = "--summary";
        argv2[argc2++] = (char *)summary;
    }
    argv2[argc2++] = "--sample-label";
    size_t labels_len = 0;
    for (size_t i = 0; i < samples.count; ++i) labels_len += strlen(samples.records[i].id) + 1;
    char *labels_csv = (char *)malloc(labels_len + 1);
    if (labels_csv == NULL) {
        free(argv2);
        free_table(&samples);
        return 1;
    }
    labels_csv[0] = '\0';
    for (size_t i = 0; i < samples.count; ++i) {
        if (i > 0) strcat(labels_csv, ",");
        strcat(labels_csv, samples.records[i].id);
    }
    argv2[argc2++] = labels_csv;
    argv2[argc2++] = "--ambiguous";
    argv2[argc2++] = ambiguity_report ? "report" : "discard";

    int rc = run_count(argv0, argc2, argv2);
    free(labels_csv);
    free(argv2);
    free_table(&samples);
    return rc;
}

static int run_validate(const char *argv0, int argc, char **argv) {
    (void)argv0;
    const char *targets_path = NULL;
    const char *reads_path = NULL;
    size_t target_start = 0;
    size_t target_length = 0;
    int k = 0;
    qdaln_metric metric = QDALN_METRIC_HAMMING;
    int indel_window = 0;
    qdaln_offset_mode offset_mode = QDALN_OFFSET_BEST;
    int threads = 1;
    const char *oracle = "scan";
    size_t sample_limit = 100000;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--targets") == 0 && i + 1 < argc) {
            targets_path = argv[++i];
        } else if (strcmp(argv[i], "--reads") == 0 && i + 1 < argc) {
            reads_path = argv[++i];
        } else if (strcmp(argv[i], "--target-start") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_start) != 0) return 2;
        } else if (strcmp(argv[i], "--target-length") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_length) != 0) return 2;
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0) return 2;
        } else if (strcmp(argv[i], "--metric") == 0 && i + 1 < argc) {
            int ok = 0;
            metric = parse_metric_value(argv[++i], &ok);
            if (!ok) return 2;
        } else if (strcmp(argv[i], "--indel-window") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &indel_window) != 0) return 2;
        } else if (strcmp(argv[i], "--offset-mode") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "best") == 0) {
                offset_mode = QDALN_OFFSET_BEST;
            } else if (strcmp(value, "multi") == 0) {
                offset_mode = QDALN_OFFSET_MULTI;
            } else {
                return 2;
            }
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &threads) != 0) return 2;
        } else if (strcmp(argv[i], "--oracle") == 0 && i + 1 < argc) {
            oracle = argv[++i];
        } else if (strcmp(argv[i], "--sample") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &sample_limit) != 0) return 2;
        } else {
            return 2;
        }
    }
    if (targets_path == NULL || reads_path == NULL || target_length == 0 || k < 0) return 2;
    if (threads <= 0 || (metric == QDALN_METRIC_HAMMING && k > 1) || indel_window < 0 || indel_window > 1) return 2;
    if (strcmp(oracle, "scan") != 0 && strcmp(oracle, "edlib") != 0) return 2;

    seq_table targets = {0};
    if (read_target_table(targets_path, &targets) != 0) {
        free_table(&targets);
        return 1;
    }

    const char **target_ptrs = (const char **)malloc(targets.count * sizeof(char *));
    size_t *target_lens = (size_t *)malloc(targets.count * sizeof(size_t));
    if ((targets.count != 0) && (target_ptrs == NULL || target_lens == NULL)) {
        free(target_ptrs);
        free(target_lens);
        free_table(&targets);
        return 1;
    }
    for (size_t i = 0; i < targets.count; ++i) {
        target_ptrs[i] = targets.records[i].seq;
        target_lens[i] = targets.records[i].len;
    }

    qdaln_index *index = qdaln_index_build_metric(target_ptrs, target_lens, targets.count, metric, -1);
    if (index == NULL) {
        free(target_ptrs); free(target_lens); free_table(&targets);
        return 1;
    }

    FILE *text_fp = NULL;
    gzFile gz_fp = NULL;
    int is_gz = 0;
    size_t checked = 0;
    size_t mismatches = 0;
    long long candidates_considered = 0;
    long long candidates_verified = 0;
    while (checked < sample_limit) {
        char *id = NULL;
        char *seq = NULL;
        char *qual = NULL;
        int ok = parse_text_or_gz_fastq_record(reads_path, &text_fp, &gz_fp, &is_gz, &id, &seq, &qual);
        free(id);
        free(qual);
        if (ok == 0) break;
        if (ok < 0) {
            free(seq);
            mismatches = sample_limit + 1;
            break;
        }
        char *window = NULL;
        if (slice_window(seq, strlen(seq), target_start, target_length, &window) != 0) {
            free(seq);
            continue;
        }
        qdaln_match_result fast = {0};
        qdaln_index_stats stats = {0};
        if (qdaln_index_match_read_stats(index, window, target_length, k, QDALN_POLICY_BEST, metric, &fast, &stats) != 0) {
            free(window); free(seq);
            mismatches = sample_limit + 1;
            break;
        }
        candidates_considered += stats.candidates_considered;
        candidates_verified += stats.candidates_verified;
        qdaln_match_result slow = {0};
        int oracle_rc = 0;
        if (strcmp(oracle, "edlib") == 0) {
#ifdef DOTMATCH_HAVE_EDLIB
            oracle_rc = qdaln_match_read_edlib(window, target_length, target_ptrs, target_lens, targets.count,
                                               k, QDALN_POLICY_BEST, metric, &slow);
#else
            oracle_rc = -1;
#endif
        } else {
            oracle_rc = qdaln_match_read_scan(window, target_length, target_ptrs, target_lens, targets.count,
                                              k, QDALN_POLICY_BEST, metric, &slow);
        }
        if (oracle_rc != 0 || memcmp(&fast, &slow, sizeof(qdaln_match_result)) != 0) {
            ++mismatches;
        }
        ++checked;
        free(window); free(seq);
    }
    if (text_fp != NULL) fclose(text_fp);
    if (gz_fp != NULL) gzclose(gz_fp);
    printf("{\n  \"oracle\": \"");
    printf("%s", oracle);
    printf("\",\n  \"checked_reads\": %zu,\n  \"mismatches\": %zu,\n  \"k\": %d,\n  \"target_start\": %zu,\n  \"target_length\": %zu,\n  \"candidates_considered\": %lld,\n  \"candidates_verified\": %lld\n}\n",
           checked, mismatches, k, target_start, target_length, candidates_considered, candidates_verified);

    qdaln_index_free(index);
    free(target_ptrs); free(target_lens); free_table(&targets);
    return mismatches == 0 ? 0 : 1;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 2;
    }

    if (strcmp(argv[1], "--version") == 0) {
        printf("dotmatch %s\n", DOTMATCH_VERSION);
        return 0;
    }
    if (strcmp(argv[1], "dist") == 0) {
        if (argc != 4) {
            usage(argv[0]);
            return 2;
        }
        printf("%d\n", qdaln_edit_distance(argv[2], strlen(argv[2]), argv[3], strlen(argv[3])));
        return 0;
    }
    if (strcmp(argv[1], "leq") == 0) {
        if (argc != 5) {
            usage(argv[0]);
            return 2;
        }
        int k = 0;
        if (sscanf(argv[2], "%d", &k) != 1) {
            usage(argv[0]);
            return 2;
        }
        int ok = qdaln_edit_distance_leq(argv[3], strlen(argv[3]), argv[4], strlen(argv[4]), k);
        printf("%s\n", ok ? "true" : "false");
        return 0;
    }
    if (strcmp(argv[1], "assign") == 0) return run_batch(argv[0], argc, argv, "assign");
    if (strcmp(argv[1], "match") == 0) return run_batch(argv[0], argc, argv, "match");
    if (strcmp(argv[1], "fastq-assign") == 0) return run_fastq_assign(argv[0], argc, argv);
    if (strcmp(argv[1], "pair-count") == 0) return qdaln_run_pair_count(argc, argv);
    if (strcmp(argv[1], "demux") == 0) return qdaln_run_demux(argc, argv);
    if (strcmp(argv[1], "bcl-demux") == 0) return run_bcl_demux(argv[0], argc, argv);
    if (strcmp(argv[1], "bcl-validate") == 0) return run_bcl_validate(argv[0], argc, argv);
    if (strcmp(argv[1], "count") == 0) return run_count(argv[0], argc, argv);
    if (strcmp(argv[1], "crispr-count") == 0) return run_crispr_count(argv[0], argc, argv);
    if (strcmp(argv[1], "inspect-unmatched") == 0) return run_inspect_unmatched(argv[0], argc, argv);
    if (strcmp(argv[1], "audit") == 0) return run_audit(argv[0], argc, argv);
    if (strcmp(argv[1], "validate") == 0) return run_validate(argv[0], argc, argv);

    usage(argv[0]);
    return 2;
}
