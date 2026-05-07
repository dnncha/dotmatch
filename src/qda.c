#include "qdalign.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>
#include <stdint.h>
#include <limits.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <unistd.h>
#include <zlib.h>

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
    fprintf(stderr, "  %s dist SEQ1 SEQ2\n", argv0);
    fprintf(stderr, "  %s leq K SEQ1 SEQ2\n", argv0);
    fprintf(stderr, "  %s assign K barcodes.txt reads.txt\n", argv0);
    fprintf(stderr, "  %s match K targets.txt reads.txt\n", argv0);
    fprintf(stderr, "  %s fastq-assign --barcodes barcodes.tsv --reads reads.fastq[.gz] --barcode-start N --barcode-length L --k 0|1 --out assignments.tsv\n", argv0);
    fprintf(stderr, "  %s demux --barcodes barcodes.tsv|barcodes.csv --reads reads.fastq[.gz] --barcode-start N --barcode-length L --k 0|1 --metric hamming|levenshtein --out-dir demux_dir [--summary qc.json]\n", argv0);
    fprintf(stderr, "  %s bcl-demux --run-folder RUN --sample-sheet SampleSheet.csv --out-dir demux_dir --barcode-mismatches 0|1|1,1 [--threads N] [--gzip-level 0..9] [--emit-index-fastqs] [--summary summary.json]\n", argv0);
    fprintf(stderr, "  %s bcl-validate --dotmatch-out DIR --truth-out DIR\n", argv0);
    fprintf(stderr, "  %s count --targets targets.tsv|targets.csv --reads reads.fastq[.gz] [--reads more.fastq.gz] --sample-label labels --target-start N --target-length L --k 0|1 --metric hamming|levenshtein [--hamming-index auto|query|precompute] --ambiguity-policy best|radius --offset-mode best|multi --out counts.tsv [--format dotmatch|mageck]\n", argv0);
    fprintf(stderr, "  %s crispr-count --library guides.tsv|guides.csv --samples samples.tsv --guide-start N --guide-length L --k 0|1 --out counts.tsv [--summary qc.json]\n", argv0);
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

static int read_samples_file(const char *path, string_list *labels, string_list *reads) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return -1;
    char buf[8192];
    size_t row = 0;
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        trim_line(buf);
        if (buf[0] == '\0' || buf[0] == '#') continue;
        char *tab = strchr(buf, '\t');
        if (tab == NULL) {
            fclose(fp);
            return -1;
        }
        *tab = '\0';
        char *sample = buf;
        char *path_field = tab + 1;
        char *extra = strchr(path_field, '\t');
        if (extra != NULL) *extra = '\0';
        if (row == 0 &&
            (strcmp(sample, "sample") == 0 || strcmp(sample, "sample_id") == 0 || strcmp(sample, "sample_id ") == 0)) {
            ++row;
            continue;
        }
        if (sample[0] == '\0' || path_field[0] == '\0') {
            fclose(fp);
            return -1;
        }
        if (push_string(labels, sample) != 0 || push_string(reads, path_field) != 0) {
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

static size_t code_hash_local(uint64_t code, size_t len, size_t cap) {
    uint64_t x = code ^ ((uint64_t)len * 0x9e3779b97f4a7c15ULL);
    x *= 0x9e3779b97f4a7c15ULL;
    x ^= x >> 32;
    return (size_t)x & (cap - 1);
}

static size_t seed_hash_local(uint64_t code, size_t len, unsigned char seed_id, size_t cap) {
    return code_hash_local(code ^ ((uint64_t)seed_id * 0x517cc1b727220a95ULL), len + seed_id * 37U, cap);
}

static uint64_t code_low_mask_local(size_t len) {
    return len == 0 ? 0 : ((1ULL << (2 * len)) - 1ULL);
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

static int build_hamming_lookup(const seq_table *targets, size_t target_len, hamming_lookup *lookup) {
    memset(lookup, 0, sizeof(*lookup));
    if (target_len == 0 || target_len > 32) return 0;
    for (size_t i = 0; i < targets->count; ++i) {
        if (targets->records[i].len != target_len) return 0;
        uint64_t code = 0;
        if (!dna2_code_local(targets->records[i].seq, target_len, &code)) return 0;
    }

    size_t exact_need = targets->count * 2 + 16;
    size_t mismatch_need = targets->count * target_len * 4 + 16;
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
    return (double)(((long double)n + 1.0L - 2.0L * weighted / (long double)sum) / (long double)n);
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
    FILE *out = fopen(path, "w");
    if (out == NULL) return -1;

    fprintf(out,
            "<!doctype html>\n<html><head><meta charset=\"utf-8\"><title>DotMatch Report</title>"
            "<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:32px;line-height:1.4;color:#17202a}"
            "table{border-collapse:collapse;width:100%%;margin:16px 0}th,td{border:1px solid #d8dee4;padding:6px 8px;text-align:right}"
            "th:first-child,td:first-child{text-align:left}th{background:#f6f8fa}.metric{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}"
            ".metric div{border:1px solid #d8dee4;padding:12px;border-radius:6px}.warn{color:#9a6700}.ok{color:#1a7f37}</style></head><body>\n");
    fprintf(out, "<h1>DotMatch Report</h1>\n");
    fprintf(out, "<p>Known-target assignment report for %zu target%s and %zu sample%s.</p>\n",
            targets->count, targets->count == 1 ? "" : "s", reads->count, reads->count == 1 ? "" : "s");
    fprintf(out, "<div class=\"metric\"><div><strong>k</strong><br>%d</div><div><strong>Metric</strong><br>%s</div>"
                 "<div><strong>Ambiguity policy</strong><br>%s</div><div><strong>Target length</strong><br>%zu</div></div>\n",
            k, metric_name(metric), ambiguity_policy_name(policy), target_len);

    fprintf(out, "<h2>Sample QC</h2><table><tr><th>Sample</th><th>Total reads</th><th>Assignment rate</th>"
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
        double denom = s->total == 0 ? 1.0 : (double)s->total;
        fprintf(out, "<tr><td>");
        html_escape(out, labels->items[sample]);
        fprintf(out, "</td><td>%llu</td><td>%.2f%%</td><td>%.2f%%</td><td>%.2f%%</td><td>%.2f%%</td>"
                     "<td>%.2f%%</td><td>%.2f%%</td><td>%llu</td></tr>\n",
                s->total, 100.0 * (double)s->unique / denom, 100.0 * (double)s->exact / denom,
                100.0 * (double)s->corrected / denom, 100.0 * (double)s->ambiguous / denom,
                100.0 * (double)s->unmatched / denom,
                targets->count == 0 ? 0.0 : 100.0 * (double)covered / (double)targets->count,
                s->candidates_verified);
    }
    fprintf(out, "</table>\n");

    fprintf(out, "<h2>Warnings</h2><ul>\n");
    for (size_t sample = 0; sample < reads->count; ++sample) {
        const count_stats *s = &stats_by_sample[sample];
        double denom = s->total == 0 ? 1.0 : (double)s->total;
        if ((double)s->ambiguous / denom > 0.01) {
            fprintf(out, "<li class=\"warn\">Sample ");
            html_escape(out, labels->items[sample]);
            fprintf(out, " has ambiguous assignments above 1%%.</li>\n");
        }
        if ((double)s->unmatched / denom > 0.10) {
            fprintf(out, "<li class=\"warn\">Sample ");
            html_escape(out, labels->items[sample]);
            fprintf(out, " has no-match reads above 10%%.</li>\n");
        }
    }
    fprintf(out, "<li class=\"ok\">Ambiguous reads are not silently counted.</li></ul>\n");

    fprintf(out, "<h2>Inputs</h2><table><tr><th>Sample</th><th>FASTQ</th><th>Selected start(s)</th></tr>\n");
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

    fprintf(out, "</body></html>\n");
    fclose(out);
    return 0;
}

typedef struct count_sample_job {
    const qdaln_index *index;
    const hamming_lookup *hlookup;
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
    int fused_offset_detection;
    size_t target_start;
    size_t auto_offset;
    size_t auto_offset_sample;
    offset_mode offsets_mode;
    double offset_min_fraction;
    size_t read_threads;
    int rc;
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
            int exact_rc = qdaln_index_assign_hamming_stats(index, &read_ptr, &read_len, 1, 0, &r, &s);
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
            memcpy(exact_observed, seq + offset, target_len);
            exact_observed[target_len] = '\0';
            uppercase_ascii(exact_observed);
            const char *read_ptr = exact_observed;
            size_t read_len = target_len;
            qdaln_match_result exact_one;
            qdaln_index_stats exact_stats = {0, 0};
            if (qdaln_index_assign_stats(index, &read_ptr, &read_len, 1, 0, &exact_one, &exact_stats) != 0) {
                merge_state_free(&exact_merge);
                return -1;
            }
            exact_stats_total.candidates_considered += exact_stats.candidates_considered;
            exact_stats_total.candidates_verified += exact_stats.candidates_verified;
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
            continue;
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
    if (buffer->count == 0 && len <= 8191) {
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

static void direct_hamming_visit_seed(const count_sample_job *job, unsigned char seed_id, uint64_t seed_code,
                                      uint64_t read_code, int *best_target, int *ambiguous) {
    const hamming_lookup *lookup = job->hlookup;
    size_t seed_len = seed_id == 0 ? lookup->seed0_len : lookup->target_len - lookup->seed0_len;
    size_t slot = seed_hash_local(seed_code, seed_len, seed_id, lookup->seed_hash_cap);
    for (int e = lookup->seed_heads[slot]; e >= 0; e = lookup->seeds[e].next) {
        const hamming_seed_entry *entry = &lookup->seeds[e];
        if (entry->seed_id != seed_id || entry->code != seed_code) continue;
        int target_index = entry->target_index;
        if (target_index < 0) continue;
        job->stats->candidates_considered += 1;
        job->stats->candidates_verified += 1;
        if (hamming_code_distance_local(read_code, lookup->target_codes[target_index], lookup->target_len) > 1) {
            continue;
        }
        direct_hamming_record_hit(target_index, 1, best_target, ambiguous);
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

    size_t n_offsets = job->selected_offsets == NULL || job->selected_offsets->count == 0
            ? 1 : job->selected_offsets->count;
    if (n_offsets == 1) {
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
                direct_hamming_record_hit(entry->target_index, entry->match_count, &exact_target,
                                          &exact_ambiguous);
            }
        }

        if (exact_target >= 0) {
            if (exact_ambiguous) {
                ++job->stats->ambiguous;
            } else {
                ++job->counts[((job->sample_index * job->targets->count + (size_t)exact_target) * 5) + 0];
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
                ++job->counts[((job->sample_index * job->targets->count + (size_t)mismatch_target) * 5) + 1];
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
            ++job->counts[((job->sample_index * job->targets->count + (size_t)exact_target) * 5) + 0];
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
            ++job->counts[((job->sample_index * job->targets->count + (size_t)mismatch_target) * 5) + 1];
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
    unsigned long long *local_counts;
    count_stats local_stats;
    int rc;
} direct_hamming_batch_job;

static void *direct_hamming_batch_worker(void *arg) {
    direct_hamming_batch_job *batch = (direct_hamming_batch_job *)arg;
    batch->job.counts = batch->local_counts;
    batch->job.sample_index = 0;
    batch->job.stats = &batch->local_stats;
    batch->rc = 0;
    for (size_t i = batch->start; i < batch->end; ++i) {
        if (direct_hamming_count_seq(&batch->job, batch->items[i], batch->lens[i]) != 0) {
            batch->rc = 1;
            break;
        }
    }
    return NULL;
}

static int process_direct_hamming_buffer(count_sample_job *job, const seq_buffer *buffer) {
    if (buffer->count == 0) return 0;
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
    unsigned long long *local_counts = (unsigned long long *)calloc(read_threads * (target_slots == 0 ? 1 : target_slots),
                                                                    sizeof(unsigned long long));
    if (thread_ids == NULL || jobs == NULL || local_counts == NULL) {
        free(thread_ids);
        free(jobs);
        free(local_counts);
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
        jobs[t].local_counts = local_counts + t * target_slots;
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
            for (size_t slot = 0; slot < target_slots; ++slot) {
                job->counts[dst_offset + slot] += jobs[t].local_counts[slot];
            }
        }
    }

    free(thread_ids);
    free(jobs);
    free(local_counts);
    return rc;
}

static void score_offsets_for_seq(const hamming_lookup *lookup, const char *seq, size_t seq_len,
                                  size_t target_start, size_t target_len, size_t range,
                                  unsigned long long *scores) {
    if (lookup == NULL || !lookup->ready || lookup->target_len != target_len) return;
    size_t n_offsets = range * 2 + 1;
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
        size_t n_offsets = job->auto_offset * 2 + 1;
        unsigned long long *scores = (unsigned long long *)calloc(n_offsets == 0 ? 1 : n_offsets,
                                                                  sizeof(unsigned long long));
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
        const size_t batch_reads = 262144;
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
                free_seq_buffer(&batch);
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

static int count_sample_sequence(count_sample_job *job, const char *seq, size_t seq_len, const char *read_id) {
    char observed[8192];
    qdaln_match_result result = {-1, -1, -1, 0, QDALN_MATCH_INVALID};
    qdaln_index_stats istats = {0, 0};
    observed[0] = '\0';
    ++job->stats->total;
    int best_exact_shortcut = job->k == 1 &&
            job->assignment_policy == AMBIGUITY_POLICY_BEST && job->assignments == NULL &&
            job->ambiguous_out == NULL && job->unmatched_out == NULL;
    int handled = 0;
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
    if (result.status == QDALN_MATCH_UNIQUE && result.target_index >= 0) {
        seq_record *target = &job->targets->records[result.target_index];
        int kind = correction_kind(observed, strlen(observed), target->seq, target->len, result.best_distance);
        size_t slot = ((job->sample_index * job->targets->count + (size_t)result.target_index) * 5) + (size_t)kind;
        ++job->counts[slot];
        ++job->stats->unique;
        if (result.best_distance == 0) ++job->stats->exact;
        else ++job->stats->corrected;
        correction = correction_name(kind);
    } else if (result.status == QDALN_MATCH_AMBIGUOUS) {
        ++job->stats->ambiguous;
        correction = "ambiguous";
    } else if (result.status == QDALN_MATCH_NONE) {
        ++job->stats->unmatched;
        correction = "none";
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
    unsigned long long *local_counts;
    count_stats local_stats;
    int rc;
} count_batch_job;

static void *count_batch_worker(void *arg) {
    count_batch_job *batch = (count_batch_job *)arg;
    batch->job.counts = batch->local_counts;
    batch->job.sample_index = 0;
    batch->job.stats = &batch->local_stats;
    batch->job.assignments = NULL;
    batch->job.ambiguous_out = NULL;
    batch->job.unmatched_out = NULL;
    batch->rc = 0;
    for (size_t i = batch->start; i < batch->end; ++i) {
        if (count_sample_sequence(&batch->job, batch->items[i], batch->lens[i], NULL) != 0) {
            batch->rc = 1;
            break;
        }
    }
    return NULL;
}

static int process_count_buffer(count_sample_job *job, const seq_buffer *buffer) {
    if (buffer->count == 0) return 0;
    size_t read_threads = job->read_threads;
    if (read_threads <= 1 || buffer->count < 1024) {
        for (size_t i = 0; i < buffer->count; ++i) {
            if (count_sample_sequence(job, buffer->items[i], buffer->lens[i], NULL) != 0) return 1;
        }
        return 0;
    }
    if (read_threads > buffer->count) read_threads = buffer->count;

    size_t target_slots = job->targets->count * 5;
    pthread_t *thread_ids = (pthread_t *)calloc(read_threads, sizeof(pthread_t));
    count_batch_job *jobs = (count_batch_job *)calloc(read_threads, sizeof(count_batch_job));
    unsigned long long *local_counts = (unsigned long long *)calloc(read_threads * (target_slots == 0 ? 1 : target_slots),
                                                                    sizeof(unsigned long long));
    if (thread_ids == NULL || jobs == NULL || local_counts == NULL) {
        free(thread_ids);
        free(jobs);
        free(local_counts);
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
        jobs[t].local_counts = local_counts + t * target_slots;
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
            for (size_t slot = 0; slot < target_slots; ++slot) {
                job->counts[dst_offset + slot] += jobs[t].local_counts[slot];
            }
        }
    }

    free(thread_ids);
    free(jobs);
    free(local_counts);
    return rc;
}

static void *count_sample_worker(void *arg) {
    count_sample_job *job = (count_sample_job *)arg;
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
    size_t seq_len = 0;
    if (job->read_threads > 1 && !need_read_id) {
        const size_t batch_reads = 262144;
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
                free_seq_buffer(&batch);
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
        while ((got = need_read_id
                ? fastq_read_record_len(&reader, header, seq, plus, qual, sizeof(header), &seq_len)
                : fastq_read_sequence_record_len(&reader, seq, sizeof(seq), &seq_len)) == 1) {
            read_id[0] = '\0';
            if (need_read_id) fastq_read_id(header, read_id, sizeof(read_id));
            if (count_sample_sequence(job, seq, seq_len, read_id) != 0) {
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

    size_t n_offsets = range * 2 + 1;

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

    size_t n_offsets = range * 2 + 1;
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
                char observed[8192];
                if (target_len >= sizeof(observed)) continue;
                memcpy(observed, seq + offset, target_len);
                observed[target_len] = '\0';
                uppercase_ascii(observed);
                const char *read_ptr = observed;
                size_t read_len = target_len;
                qdaln_match_result r;
                if (qdaln_index_assign(index, &read_ptr, &read_len, 1, 0, &r) != 0) {
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
    ambiguity_policy assignment_policy = AMBIGUITY_POLICY_BEST;
    count_metric metric = COUNT_METRIC_LEVENSHTEIN;
    hamming_index_strategy hamming_strategy = HAMMING_INDEX_AUTO;
    size_t target_start = 0;
    size_t target_len = 0;
    size_t indel_window = 0;
    size_t auto_offset = 0;
    size_t auto_offset_sample = 1000;
    offset_mode offsets_mode = OFFSET_MODE_BEST;
    double offset_min_fraction = 0.005;
    size_t threads = 1;
    int k = -1;
    string_list reads = {0};
    string_list labels = {0};

    for (int i = 2; i < argc; ++i) {
        if ((strcmp(argv[i], "--targets") == 0 || strcmp(argv[i], "--library") == 0) && i + 1 < argc) {
            targets_path = argv[++i];
        } else if (strcmp(argv[i], "--samples") == 0 && i + 1 < argc) {
            samples_path = argv[++i];
        } else if (strcmp(argv[i], "--reads") == 0 && i + 1 < argc) {
            if (push_string(&reads, argv[++i]) != 0) {
                fprintf(stderr, "out of memory\n");
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--sample-label") == 0 && i + 1 < argc) {
            if (split_string_list(&labels, argv[++i], ',') != 0) {
                fprintf(stderr, "out of memory\n");
                goto fail_args;
            }
        } else if ((strcmp(argv[i], "--target-start") == 0 || strcmp(argv[i], "--guide-start") == 0) && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_start) != 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if ((strcmp(argv[i], "--target-length") == 0 || strcmp(argv[i], "--guide-length") == 0) && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_len) != 0 || target_len == 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0 || (k != 0 && k != 1)) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--metric") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "hamming") == 0) {
                metric = COUNT_METRIC_HAMMING;
            } else if (strcmp(value, "levenshtein") == 0) {
                metric = COUNT_METRIC_LEVENSHTEIN;
            } else {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--hamming-index") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
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
        } else if (strcmp(argv[i], "--indel-window") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &indel_window) != 0 || indel_window > 1) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--auto-offset") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &auto_offset) != 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--auto-offset-sample") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &auto_offset_sample) != 0 || auto_offset_sample == 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--offset-mode") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "best") == 0) {
                offsets_mode = OFFSET_MODE_BEST;
            } else if (strcmp(value, "multi") == 0) {
                offsets_mode = OFFSET_MODE_MULTI;
            } else {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--offset-min-fraction") == 0 && i + 1 < argc) {
            if (parse_double_value(argv[++i], &offset_min_fraction) != 0 ||
                offset_min_fraction < 0.0 || offset_min_fraction > 1.0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &threads) != 0 || threads == 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else if (strcmp(argv[i], "--assignments") == 0 && i + 1 < argc) {
            assignments_path = argv[++i];
        } else if (strcmp(argv[i], "--summary") == 0 && i + 1 < argc) {
            summary_path = argv[++i];
        } else if (strcmp(argv[i], "--report") == 0 && i + 1 < argc) {
            report_path = argv[++i];
        } else if (strcmp(argv[i], "--report-audit-dir") == 0 && i + 1 < argc) {
            report_audit_dir = argv[++i];
        } else if (strcmp(argv[i], "--report-unmatched") == 0 && i + 1 < argc) {
            report_unmatched_path = argv[++i];
        } else if (strcmp(argv[i], "--qc") == 0 && i + 1 < argc) {
            sample_qc_path = argv[++i];
        } else if (strcmp(argv[i], "--sample-qc") == 0 && i + 1 < argc) {
            sample_qc_path = argv[++i];
        } else if (strcmp(argv[i], "--target-counts-long") == 0 && i + 1 < argc) {
            target_counts_long_path = argv[++i];
        } else if (strcmp(argv[i], "--ambiguous-out") == 0 && i + 1 < argc) {
            ambiguous_path = argv[++i];
        } else if (strcmp(argv[i], "--unmatched-out") == 0 && i + 1 < argc) {
            unmatched_path = argv[++i];
        } else if (strcmp(argv[i], "--ambiguous") == 0 && i + 1 < argc) {
            ambiguous_policy = argv[++i];
            if (strcmp(ambiguous_policy, "discard") != 0 && strcmp(ambiguous_policy, "report") != 0) {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--ambiguity-policy") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "best") == 0) {
                assignment_policy = AMBIGUITY_POLICY_BEST;
            } else if (strcmp(value, "radius") == 0) {
                assignment_policy = AMBIGUITY_POLICY_RADIUS;
            } else {
                usage(argv0);
                goto fail_args;
            }
        } else if (strcmp(argv[i], "--format") == 0 && i + 1 < argc) {
            format = argv[++i];
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
    if (metric == COUNT_METRIC_HAMMING && indel_window != 0) {
        fprintf(stderr, "--indel-window is only valid with --metric levenshtein\n");
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
    if (threads > 1 && (assignments_path != NULL || ambiguous_path != NULL || unmatched_path != NULL)) {
        fprintf(stderr, "--threads > 1 is not supported with row-level diagnostic outputs\n");
        goto fail_args;
    }
    int count_only = assignments_path == NULL && ambiguous_path == NULL && unmatched_path == NULL;

    seq_table targets = {0};
    qdaln_index *index = NULL;
    hamming_lookup hlookup = {0};
    hamming_lookup offset_lookup = {0};
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
    size_t effective_read_threads = 1;

    double phase_start_seconds = seconds_now();
    if (read_target_table(targets_path, &targets) != 0) {
        fprintf(stderr, "failed to read targets\n");
        goto done;
    }
    if (metric == COUNT_METRIC_HAMMING && !all_targets_have_length(&targets, target_len)) {
        fprintf(stderr, "--metric hamming requires every target to have --target-length bases\n");
        goto done;
    }
    target_index_seconds = seconds_now() - phase_start_seconds;

    int direct_hamming_counts = count_only && metric == COUNT_METRIC_HAMMING && indel_window == 0 &&
            assignment_policy == AMBIGUITY_POLICY_BEST && (k == 0 || k == 1) && target_len <= 32 &&
            (hamming_strategy == HAMMING_INDEX_PRECOMPUTE || hamming_strategy == HAMMING_INDEX_AUTO);
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

    size_t total_slots = reads.count * targets.count * 5;
    counts = (unsigned long long *)calloc(total_slots == 0 ? 1 : total_slots, sizeof(unsigned long long));
    stats_by_sample = (count_stats *)calloc(reads.count == 0 ? 1 : reads.count, sizeof(count_stats));
    ambiguous_nearby = (unsigned char *)calloc(targets.count == 0 ? 1 : targets.count, sizeof(unsigned char));
    selected_offsets = (offset_list *)calloc(reads.count == 0 ? 1 : reads.count, sizeof(offset_list));
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
        assignments = fopen(assignments_path, "w");
        if (assignments == NULL) {
            fprintf(stderr, "failed to open assignments output\n");
            goto done;
        }
        fprintf(assignments, "sample\tread_id\tobserved_seq\ttarget_index\ttarget_id\ttarget_seq\tbest_distance\tsecond_best_distance\tmatch_count\tstatus\tcorrection\n");
    }
    if (ambiguous_path != NULL) {
        ambiguous_out = fopen(ambiguous_path, "w");
        if (ambiguous_out == NULL) {
            fprintf(stderr, "failed to open ambiguous output\n");
            goto done;
        }
        fprintf(ambiguous_out, "sample\tread_id\tobserved_seq\ttarget_index\ttarget_id\ttarget_seq\tbest_distance\tsecond_best_distance\tmatch_count\tstatus\tcorrection\n");
    }
    if (unmatched_path != NULL) {
        unmatched_out = fopen(unmatched_path, "w");
        if (unmatched_out == NULL) {
            fprintf(stderr, "failed to open unmatched output\n");
            goto done;
        }
        fprintf(unmatched_out, "sample\tread_id\tobserved_seq\ttarget_index\ttarget_id\ttarget_seq\tbest_distance\tsecond_best_distance\tmatch_count\tstatus\tcorrection\n");
    }

    int fused_offset_detection = direct_hamming_counts && auto_offset != 0;
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
    if (direct_hamming_counts && !fused_offset_detection && max_selected_offsets <= 1) {
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

    phase_start_seconds = seconds_now();
    size_t sample_threads = threads;
    if (direct_hamming_counts && reads.count == 1 && threads > 1) {
        effective_read_threads = threads;
        sample_threads = 1;
    } else if (count_only && reads.count == 1 && threads > 1) {
        effective_read_threads = threads;
        sample_threads = 1;
    } else if (sample_threads > reads.count) {
        sample_threads = reads.count;
    }
    if (sample_threads <= 1 || reads.count <= 1) {
        for (size_t sample = 0; sample < reads.count; ++sample) {
            count_sample_job job = {
                index, &hlookup, &targets, target_ptrs, target_lens,
                reads.items[sample], labels.items[sample], sample, &selected_offsets[sample],
                target_len, k, metric, indel_window, counts, &stats_by_sample[sample],
                assignments, ambiguous_out, unmatched_out, ambiguous_policy, assignment_policy,
                direct_hamming_counts, fused_offset_detection, target_start, auto_offset, auto_offset_sample,
                offsets_mode, offset_min_fraction, effective_read_threads, 1
            };
            count_sample_worker(&job);
            if (job.rc != 0) goto done;
        }
    } else {
        pthread_t *thread_ids = (pthread_t *)calloc(sample_threads, sizeof(pthread_t));
        count_sample_job *jobs = (count_sample_job *)calloc(reads.count, sizeof(count_sample_job));
        if (thread_ids == NULL || jobs == NULL) {
            free(thread_ids);
            free(jobs);
            fprintf(stderr, "out of memory\n");
            goto done;
        }
        size_t next_sample = 0;
        while (next_sample < reads.count) {
            size_t batch = reads.count - next_sample;
            if (batch > sample_threads) batch = sample_threads;
            for (size_t i = 0; i < batch; ++i) {
                size_t sample = next_sample + i;
                jobs[sample] = (count_sample_job){
                    index, &hlookup, &targets, target_ptrs, target_lens,
                    reads.items[sample], labels.items[sample], sample, &selected_offsets[sample],
                    target_len, k, metric, indel_window, counts, &stats_by_sample[sample],
                    NULL, NULL, NULL, ambiguous_policy, assignment_policy,
                    direct_hamming_counts, fused_offset_detection, target_start, auto_offset, auto_offset_sample,
                    offsets_mode, offset_min_fraction, 1, 1
                };
                if (pthread_create(&thread_ids[i], NULL, count_sample_worker, &jobs[sample]) != 0) {
                    fprintf(stderr, "failed to create worker thread\n");
                    batch = i;
                    for (size_t j = 0; j < batch; ++j) pthread_join(thread_ids[j], NULL);
                    free(thread_ids);
                    free(jobs);
                    goto done;
                }
            }
            for (size_t i = 0; i < batch; ++i) {
                pthread_join(thread_ids[i], NULL);
                if (jobs[next_sample + i].rc != 0) {
                    free(thread_ids);
                    free(jobs);
                    goto done;
                }
            }
            next_sample += batch;
        }
        free(thread_ids);
        free(jobs);
    }
    counting_seconds = seconds_now() - phase_start_seconds;

    out = fopen(out_path, "w");
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
        FILE *long_out = fopen(target_counts_long_path, "w");
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
        FILE *qc = fopen(sample_qc_path, "w");
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
            double denom = s->total == 0 ? 1.0 : (double)s->total;
            fprintf(qc, "%s\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%llu\t%llu\t%.8f\t%.8f\t%llu\n",
                    labels.items[sample], reads.items[sample], s->total, valid, s->unique, s->exact, s->corrected,
                    sub, ins, del, s->ambiguous, s->unmatched, s->invalid,
                    (double)s->unique / denom, (double)s->exact / denom, (double)s->corrected / denom,
                    (double)s->ambiguous / denom, (double)s->unmatched / denom,
                    observed_targets, (unsigned long long)(targets.count - observed_targets),
                    gini_from_counts(target_totals, targets.count),
                    top_fraction_from_counts(target_totals, targets.count, 0.01),
                    s->candidates_verified);
            free(target_totals);
        }
        fclose(qc);
    }

    if (summary_path != NULL) {
        FILE *summary = fopen(summary_path, "w");
        if (summary == NULL) {
            fprintf(stderr, "failed to open summary output\n");
            goto done;
        }
        double total_before_summary_seconds = seconds_now() - run_start_seconds;
        fprintf(summary,
                "{\n  \"k\": %d,\n  \"metric\": \"%s\",\n  \"ambiguity_policy\": \"%s\",\n  \"indel_window\": %zu,\n  \"target_start\": %zu,\n  \"auto_offset\": %zu,\n  \"offset_mode\": \"%s\",\n  \"offset_min_fraction\": %.8f,\n  \"offset_detection_strategy\": \"%s\",\n  \"count_engine\": \"%s\",\n  \"hamming_index\": \"%s\",\n  \"target_length\": %zu,\n  \"n_targets\": %zu,\n  \"read_threads\": %zu,\n  \"phase_seconds\": {\"target_index\": %.6f, \"offset_detection\": %.6f, \"hamming_precompute\": %.6f, \"counting\": %.6f, \"total_before_summary\": %.6f},\n  \"samples\": [\n",
                k, metric_name(metric), ambiguity_policy_name(assignment_policy), indel_window, target_start,
                auto_offset, offset_mode_name(offsets_mode), offset_min_fraction,
                offset_detection_strategy, count_engine, hamming_lookup_kind(&hlookup), target_len, targets.count,
                effective_read_threads,
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
    rc = 0;

done:
    if (out != NULL) fclose(out);
    if (assignments != NULL) fclose(assignments);
    if (ambiguous_out != NULL) fclose(ambiguous_out);
    if (unmatched_out != NULL) fclose(unmatched_out);
    qdaln_index_free(index);
    free_hamming_lookup(&hlookup);
    free_hamming_lookup(&offset_lookup);
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

static int run_fastq_assign(const char *argv0, int argc, char **argv) {
    const char *barcodes_path = NULL;
    const char *reads_path = NULL;
    const char *out_path = NULL;
    size_t barcode_start = 0;
    size_t barcode_len = 0;
    int k = -1;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--barcodes") == 0 && i + 1 < argc) {
            barcodes_path = argv[++i];
        } else if (strcmp(argv[i], "--reads") == 0 && i + 1 < argc) {
            reads_path = argv[++i];
        } else if (strcmp(argv[i], "--barcode-start") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &barcode_start) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--barcode-length") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &barcode_len) != 0 || barcode_len == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0 || (k != 0 && k != 1)) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out_path = argv[++i];
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
    out = fopen(out_path, "w");
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
                                    const char *min_edit_distance_json,
                                    int safe_at_k0, int safe_at_k1, const char *safe_at_k2_json,
                                    unsigned long long pairs_d0, unsigned long long pairs_d1,
                                    unsigned long long pairs_d2, unsigned long long pairs_within_k,
                                    unsigned long long risk_pairs_k1, const char *risk_pairs_k2_json,
                                    unsigned long long ambiguous_query_variants_k1, int recommended_k) {
    char path[4096];
    if (path_join(path, sizeof(path), out_dir, "audit_summary.json") != 0) return -1;
    FILE *out = fopen(path, "w");
    if (out == NULL) return -1;
    fprintf(out,
            "{\n"
            "  \"audit_mode\": \"%s\",\n"
            "  \"k\": %d,\n"
            "  \"targets\": %zu,\n"
            "  \"unique_sequences\": %zu,\n"
            "  \"duplicate_sequences\": %zu,\n"
            "  \"min_edit_distance\": %s,\n"
            "  \"safe_at_k0\": %s,\n"
            "  \"safe_at_k1\": %s,\n"
            "  \"safe_at_k2\": %s,\n"
            "  \"pairs_distance_0\": %llu,\n"
            "  \"pairs_distance_1\": %llu,\n"
            "  \"pairs_distance_2\": %llu,\n"
            "  \"pairs_within_requested_k\": %llu,\n"
            "  \"risk_pairs_for_k1\": %llu,\n"
            "  \"risk_pairs_for_k2\": %s,\n"
            "  \"ambiguous_query_variants_k1\": %llu,\n"
            "  \"recommended_k\": %d\n"
            "}\n",
            audit_mode, k, n_targets, unique_sequences, n_targets - unique_sequences, min_edit_distance_json,
            safe_at_k0 ? "true" : "false", safe_at_k1 ? "true" : "false", safe_at_k2_json,
            pairs_d0, pairs_d1, pairs_d2, pairs_within_k, risk_pairs_k1, risk_pairs_k2_json,
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
    variants_out = fopen(path, "w");
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
    pairs = fopen(path, "w");
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
    safety = fopen(path, "w");
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
    clusters = fopen(path, "w");
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
    summary = fopen(path, "w");
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
    fprintf(summary, "pairs_distance_0\t%llu\n", pairs_d0);
    fprintf(summary, "pairs_distance_1\t%llu\n", pairs_d1);
    fprintf(summary, "pairs_distance_2\t%llu\n", pairs_d2);
    fprintf(summary, "pairs_within_requested_k\t%llu\n", pairs_within_k);
    fprintf(summary, "risk_pairs_for_k1\t%llu\n", risk_pairs_k1);
    fprintf(summary, "risk_pairs_for_k2\tnot_computed\n");
    fprintf(summary, "ambiguous_query_variants_k1\t%llu\n", ambiguous_query_variants_k1);
    fprintf(summary, "recommended_k\t%d\n", pairs_d0 == 0 && (k == 0 || risk_pairs_k1 == 0) ? k : 0);
    fclose(summary);
    summary = NULL;
    if (write_audit_summary_json(out_dir, "fast", k, targets->count, unique_sequences,
                                 min_dist < 0 ? "\">=3\"" : (min_dist == 0 ? "0" : (min_dist == 1 ? "1" : "2")),
                                 pairs_d0 == 0, risk_pairs_k1 == 0, "null",
                                 pairs_d0, pairs_d1, pairs_d2, pairs_within_k, risk_pairs_k1, "null",
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

    for (int i = 2; i < argc; ++i) {
        if ((strcmp(argv[i], "--targets") == 0 || strcmp(argv[i], "--library") == 0) && i + 1 < argc) {
            targets_path = argv[++i];
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0 || k < 0 || k > 2) {
                usage(argv0);
                return 2;
            }
        } else if ((strcmp(argv[i], "--out-dir") == 0 || strcmp(argv[i], "--out") == 0) && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (strcmp(argv[i], "--audit-mode") == 0 && i + 1 < argc) {
            audit_mode = argv[++i];
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
    pairs = fopen(path, "w");
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
    safety = fopen(path, "w");
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
    clusters = fopen(path, "w");
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
    variants_out = fopen(path, "w");
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
    summary = fopen(path, "w");
    if (summary == NULL) goto done;
    fprintf(summary, "metric\tvalue\n");
    fprintf(summary, "audit_mode\texact\n");
    fprintf(summary, "targets\t%zu\n", targets.count);
    fprintf(summary, "unique_sequences\t%zu\n", unique_sequences);
    fprintf(summary, "duplicate_sequences\t%zu\n", targets.count - unique_sequences);
    fprintf(summary, "min_edit_distance\t%d\n", min_dist);
    fprintf(summary, "safe_at_k0\t%s\n", pairs_d0 == 0 ? "yes" : "no");
    fprintf(summary, "safe_at_k1\t%s\n", risk_pairs_k1 == 0 ? "yes" : "no");
    fprintf(summary, "safe_at_k2\t%s\n", risk_pairs_k2 == 0 ? "yes" : "no");
    fprintf(summary, "pairs_distance_0\t%llu\n", pairs_d0);
    fprintf(summary, "pairs_distance_1\t%llu\n", pairs_d1);
    fprintf(summary, "pairs_distance_2\t%llu\n", pairs_d2);
    fprintf(summary, "pairs_within_requested_k\t%llu\n", pairs_within_k);
    fprintf(summary, "risk_pairs_for_k1\t%llu\n", risk_pairs_k1);
    fprintf(summary, "risk_pairs_for_k2\t%llu\n", risk_pairs_k2);
    fprintf(summary, "ambiguous_query_variants_k1\t%llu\n", ambiguous_query_variants_k1);
    fprintf(summary, "recommended_k\t%d\n", pairs_d0 == 0 && (k == 0 || risk_pairs_k1 == 0) ? k : 0);
    fclose(summary);
    summary = NULL;
    char min_dist_json[32];
    char risk_pairs_k2_json[32];
    snprintf(min_dist_json, sizeof(min_dist_json), "%d", min_dist);
    snprintf(risk_pairs_k2_json, sizeof(risk_pairs_k2_json), "%llu", risk_pairs_k2);
    if (write_audit_summary_json(out_dir, "exact", k, targets.count, unique_sequences,
                                 min_dist_json, pairs_d0 == 0, risk_pairs_k1 == 0,
                                 risk_pairs_k2 == 0 ? "true" : "false",
                                 pairs_d0, pairs_d1, pairs_d2, pairs_within_k, risk_pairs_k1,
                                 risk_pairs_k2_json, ambiguous_query_variants_k1,
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
            memcpy(observed, seq + offset, target_len);
            observed[target_len] = '\0';
            uppercase_ascii(observed);
            const char *read_ptr = observed;
            size_t read_len = target_len;
            qdaln_match_result r;
            qdaln_index_stats stats;
            if (qdaln_index_assign_stats(index, &read_ptr, &read_len, 1, k, &r, &stats) != 0) return 0;
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

    for (int i = 2; i < argc; ++i) {
        if ((strcmp(argv[i], "--targets") == 0 || strcmp(argv[i], "--library") == 0) && i + 1 < argc) {
            targets_path = argv[++i];
        } else if (strcmp(argv[i], "--reads") == 0 && i + 1 < argc) {
            reads_path = argv[++i];
        } else if ((strcmp(argv[i], "--target-start") == 0 || strcmp(argv[i], "--guide-start") == 0) && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_start) != 0) {
                usage(argv0);
                return 2;
            }
        } else if ((strcmp(argv[i], "--target-length") == 0 || strcmp(argv[i], "--guide-length") == 0) && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_len) != 0 || target_len == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0 || (k != 0 && k != 1)) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--top") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &top_n) != 0 || top_n == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--offset-window") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &offset_window) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--adapter") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            size_t n = strlen(value);
            if (n == 0 || n >= sizeof(adapter)) {
                usage(argv0);
                return 2;
            }
            memcpy(adapter, value, n + 1);
            uppercase_ascii(adapter);
        } else if (strcmp(argv[i], "--low-quality-threshold") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &low_quality_threshold) != 0 ||
                low_quality_threshold < 0 || low_quality_threshold > 93) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out_path = argv[++i];
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
            const char *read_ptr = observed;
            size_t read_len = target_len;
            qdaln_index_stats stats;
            if (qdaln_index_assign_stats(index, &read_ptr, &read_len, 1, k, &r, &stats) != 0) {
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
    out = fopen(out_path, "w");
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

static FILE *open_demux_target_file(FILE **files, const seq_table *targets, size_t target_index, const char *out_dir) {
    if (files[target_index] != NULL) return files[target_index];
    char safe_id[512];
    sanitize_filename(targets->records[target_index].id, safe_id, sizeof(safe_id));
    char path[4096];
    int n = snprintf(path, sizeof(path), "%s/%s.fastq", out_dir, safe_id);
    if (n < 0 || (size_t)n >= sizeof(path)) return NULL;
    files[target_index] = fopen(path, "w");
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
    size_t indel_window = 0;
    int k = -1;

    for (int i = 2; i < argc; ++i) {
        if ((strcmp(argv[i], "--barcodes") == 0 || strcmp(argv[i], "--targets") == 0) && i + 1 < argc) {
            barcodes_path = argv[++i];
        } else if (strcmp(argv[i], "--reads") == 0 && i + 1 < argc) {
            reads_path = argv[++i];
        } else if ((strcmp(argv[i], "--barcode-start") == 0 || strcmp(argv[i], "--target-start") == 0) && i + 1 < argc) {
            if (parse_size_value(argv[++i], &barcode_start) != 0) {
                usage(argv0);
                return 2;
            }
        } else if ((strcmp(argv[i], "--barcode-length") == 0 || strcmp(argv[i], "--target-length") == 0) && i + 1 < argc) {
            if (parse_size_value(argv[++i], &barcode_len) != 0 || barcode_len == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0 || (k != 0 && k != 1)) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--metric") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "hamming") == 0) {
                metric = COUNT_METRIC_HAMMING;
            } else if (strcmp(value, "levenshtein") == 0) {
                metric = COUNT_METRIC_LEVENSHTEIN;
            } else {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--indel-window") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &indel_window) != 0 || indel_window > 1) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--out-dir") == 0 && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (strcmp(argv[i], "--summary") == 0 && i + 1 < argc) {
            summary_path = argv[++i];
        } else if (strcmp(argv[i], "--qc") == 0 && i + 1 < argc) {
            summary_path = argv[++i];
        } else if (strcmp(argv[i], "--assignments") == 0 && i + 1 < argc) {
            assignments_path = argv[++i];
        } else if (strcmp(argv[i], "--ambiguous-out") == 0 && i + 1 < argc) {
            ambiguous_path = argv[++i];
        } else if (strcmp(argv[i], "--unmatched-out") == 0 && i + 1 < argc) {
            unmatched_path = argv[++i];
        } else {
            usage(argv0);
            return 2;
        }
    }

    if (barcodes_path == NULL || reads_path == NULL || out_dir == NULL || barcode_len == 0 || k < 0) {
        usage(argv0);
        return 2;
    }
    if (metric == COUNT_METRIC_HAMMING && indel_window != 0) {
        fprintf(stderr, "--indel-window is only valid with --metric levenshtein\n");
        return 2;
    }
    if (indel_window != 0 && k != 1) {
        fprintf(stderr, "--indel-window requires --k 1\n");
        return 2;
    }

    seq_table targets = {0};
    fastq_reader reader = {0};
    qdaln_index *index = NULL;
    const char **target_ptrs = NULL;
    size_t *target_lens = NULL;
    FILE **target_files = NULL;
    FILE *assignments = NULL;
    FILE *ambiguous_out = NULL;
    FILE *unmatched_out = NULL;
    count_stats stats = {0};
    unsigned long long *target_counts = NULL;
    int rc = 1;

    if (read_target_table(barcodes_path, &targets) != 0) {
        fprintf(stderr, "failed to read barcodes\n");
        goto done;
    }
    if (metric == COUNT_METRIC_HAMMING && !all_targets_have_length(&targets, barcode_len)) {
        fprintf(stderr, "--metric hamming requires every barcode to have --barcode-length bases\n");
        goto done;
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
        assignments = fopen(assignments_path, "w");
        if (assignments == NULL) {
            fprintf(stderr, "failed to open assignments output\n");
            goto done;
        }
        fprintf(assignments, "read_id\tobserved_barcode\ttarget_index\ttarget_id\ttarget_seq\tbest_distance\tsecond_best_distance\tmatch_count\tstatus\n");
    }
    if (ambiguous_path != NULL) {
        ambiguous_out = fopen(ambiguous_path, "w");
        if (ambiguous_out == NULL) {
            fprintf(stderr, "failed to open ambiguous FASTQ output\n");
            goto done;
        }
    }
    if (unmatched_path != NULL) {
        unmatched_out = fopen(unmatched_path, "w");
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

        if (assign_count_window(index, seq, seq_len, barcode_start, barcode_len, k, metric, indel_window,
                                &result, &istats, observed, sizeof(observed), 0) != 0) {
            fprintf(stderr, "FASTQ assignment failed\n");
            goto done;
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
        FILE *summary = fopen(summary_path, "w");
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
                "{\n  \"workflow\": \"demux\",\n  \"k\": %d,\n  \"metric\": \"%s\",\n  \"indel_window\": %zu,\n  \"barcode_start\": %zu,\n  \"barcode_length\": %zu,\n  \"n_barcodes\": %zu,\n  \"total_reads\": %llu,\n  \"assigned_unique\": %llu,\n  \"assigned_exact\": %llu,\n  \"assigned_corrected\": %llu,\n  \"ambiguous\": %llu,\n  \"unmatched\": %llu,\n  \"invalid\": %llu,\n  \"nonempty_outputs\": %zu,\n  \"top_barcode_id\": \"%s\",\n  \"top_barcode_count\": %llu,\n  \"candidates_considered\": %llu,\n  \"candidates_verified\": %llu\n}\n",
                k, metric_name(metric), indel_window, barcode_start, barcode_len, targets.count, stats.total,
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
    free(target_ptrs);
    free(target_lens);
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
    snprintf(path, sizeof(path), "%s/RunInfo.xml", run_folder);
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
            bcl_read_info *r = &info->reads[info->read_count++];
            r->number = atoi(number);
            r->cycles = (size_t)strtoul(cycles, NULL, 10);
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
    while (fgets(buf, sizeof(buf), fp) != NULL) {
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
            continue;
        }
        const char *name = (name_col >= 0 && (size_t)name_col < nf) ? fields[name_col] : fields[id_col];
        const char *index2 = (index2_col >= 0 && (size_t)index2_col < nf) ? fields[index2_col] : "";
        int lane = 0;
        if (lane_col >= 0 && (size_t)lane_col < nf && fields[lane_col][0] != '\0') lane = atoi(fields[lane_col]);
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
    if (n < 0) {
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
    if (count < 0) {
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
    char path[4096];
    int n = snprintf(path, sizeof(path), "%s/%s_S%zu_L%03d_%c%d_001.fastq.gz", out_dir, safe_id, sample_number,
                     lane, read_kind, read_number);
    if (n < 0 || (size_t)n >= sizeof(path)) return NULL;
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
        size_t n = sample_count * read_count;
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
    result->sample_buffers = (text_buffer *)calloc((sample_count == 0 ? 1 : sample_count) * (read_count == 0 ? 1 : read_count), sizeof(text_buffer));
    result->undetermined_buffers = (text_buffer *)calloc(read_count == 0 ? 1 : read_count, sizeof(text_buffer));
    result->sample_assigned = (unsigned long long *)calloc(sample_count == 0 ? 1 : sample_count, sizeof(unsigned long long));
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
        char index1[512] = "";
        char index2[512] = "";
        char seqs[16][8192];
        char quals[16][8192];
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
                if (indexed_seen == 0) snprintf(index1, sizeof(index1), "%s", seq_out);
                else if (indexed_seen == 1) snprintf(index2, sizeof(index2), "%s", seq_out);
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
                char header[1024];
                snprintf(header, sizeof(header), "@DOTMATCH:1:%s:%zu %d:N:0:%s%s%s",
                         job->tile, cluster + 1, bcl_output_number(run, r), index1, index2[0] ? "+" : "", index2);
                text_buffer *buf = &result->sample_buffers[out_i * run->read_count + r];
                if (append_fastq_record(buf, header, seqs[r], quals[r]) != 0) return -1;
            }
            ++result->sample_assigned[out_i];
        } else {
            char unknown_index[1024];
            snprintf(unknown_index, sizeof(unknown_index), "%s%s%s", index1, index2[0] ? "+" : "", index2);
            if (add_bcl_unknown(&result->unknowns, unknown_index) != 0) return -1;
            for (size_t r = 0; r < run->read_count; ++r) {
                const bcl_read_info *ri = &run->reads[r];
                if (!bcl_output_enabled(ri, job->emit_index_fastqs)) continue;
                char header[1024];
                snprintf(header, sizeof(header), "@DOTMATCH:1:%s:%zu %d:N:0:%s%s%s",
                         job->tile, cluster + 1, bcl_output_number(run, r), index1, index2[0] ? "+" : "", index2);
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
            if (gz != NULL && gzwrite(gz, buf->data, (unsigned int)buf->len) != (int)buf->len) return -1;
        }
    }
    for (size_t r = 0; r < read_count; ++r) {
        const text_buffer *buf = &result->undetermined_buffers[r];
        if (buf->len == 0) continue;
        if (undetermined_fastqs[r] != NULL && gzwrite(undetermined_fastqs[r], buf->data, (unsigned int)buf->len) != (int)buf->len) return -1;
    }
    return 0;
}

static int parse_mismatches(const char *s, int *k1, int *k2) {
    char *comma = strchr(s, ',');
    if (comma == NULL) {
        int k = atoi(s);
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
    int a = atoi(left);
    int b = atoi(comma + 1);
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
    int k1 = 1;
    int k2 = 1;
    int emit_index_fastqs = 0;
    size_t requested_threads = 1;
    int gzip_level = 1;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--run-folder") == 0 && i + 1 < argc) {
            run_folder = argv[++i];
        } else if (strcmp(argv[i], "--sample-sheet") == 0 && i + 1 < argc) {
            sample_sheet = argv[++i];
        } else if (strcmp(argv[i], "--out-dir") == 0 && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (strcmp(argv[i], "--summary") == 0 && i + 1 < argc) {
            summary_path = argv[++i];
        } else if (strcmp(argv[i], "--barcode-mismatches") == 0 && i + 1 < argc) {
            mismatches = argv[++i];
        } else if (strcmp(argv[i], "--emit-index-fastqs") == 0) {
            emit_index_fastqs = 1;
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &requested_threads) != 0 || requested_threads == 0) {
                fprintf(stderr, "invalid --threads value\n");
                return 2;
            }
        } else if (strcmp(argv[i], "--gzip-level") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &gzip_level) != 0 || gzip_level < 0 || gzip_level > 9) {
                fprintf(stderr, "invalid --gzip-level value\n");
                return 2;
            }
        } else if (strcmp(argv[i], "--lanes") == 0 && i + 1 < argc) {
            ++i;
        } else if (strcmp(argv[i], "--interop-dir") == 0 && i + 1 < argc) {
            ++i;
        } else {
            usage(argv0);
            return 2;
        }
    }
    if (run_folder == NULL || sample_sheet == NULL || out_dir == NULL || parse_mismatches(mismatches, &k1, &k2) != 0) {
        usage(argv0);
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
    if (ensure_dir(out_dir) != 0) {
        fprintf(stderr, "failed to create BCL output directory\n");
        goto done;
    }
    if (collect_tiles(basecalls, 1, &tiles, &tile_count) != 0) {
        fprintf(stderr, "failed to find classic BCL tiles; CBCL is not supported in this milestone\n");
        goto done;
    }

    size_t output_file_count = (samples.count == 0 ? 1 : samples.count) * (run.read_count == 0 ? 1 : run.read_count);
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
        if (threads == 0) threads = 1;
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
    snprintf(stats_path, sizeof(stats_path), "%s/Demultiplex_Stats.csv", out_dir);
    FILE *stats = fopen(stats_path, "w");
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
        snprintf(unknown_path, sizeof(unknown_path), "%s/Top_Unknown_Barcodes.csv", out_dir);
        FILE *unknown = fopen(unknown_path, "w");
        if (unknown == NULL) goto done;
        fprintf(unknown, "index,count\n");
        size_t n = unknowns.count < 100 ? unknowns.count : 100;
        for (size_t i = 0; i < n; ++i) fprintf(unknown, "%s,%llu\n", unknowns.items[i].index, unknowns.items[i].count);
        fclose(unknown);
    }

    char normalized_path[4096];
    snprintf(normalized_path, sizeof(normalized_path), "%s/SampleSheet.normalized.csv", out_dir);
    FILE *normalized = fopen(normalized_path, "w");
    if (normalized != NULL) {
        fprintf(normalized, "sample_id,sample_name,lane,index,index2\n");
        for (size_t i = 0; i < samples.count; ++i) {
            fprintf(normalized, "%s,%s,%d,%s,%s\n", samples.items[i].id, samples.items[i].name,
                    samples.items[i].lane, samples.items[i].index1, samples.items[i].index2);
        }
        fclose(normalized);
    }

    if (summary_path != NULL) {
        FILE *summary = fopen(summary_path, "w");
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
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--dotmatch-out") == 0 && i + 1 < argc) {
            dotmatch_out = argv[++i];
        } else if (strcmp(argv[i], "--truth-out") == 0 && i + 1 < argc) {
            truth_out = argv[++i];
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
        snprintf(truth_path, sizeof(truth_path), "%s/%s", truth_out, ent->d_name);
        snprintf(dotmatch_path, sizeof(dotmatch_path), "%s/%s", dotmatch_out, ent->d_name);
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
    size_t threads = 1;
    offset_mode offsets_mode = OFFSET_MODE_BEST;
    double offset_min_fraction = 0.005;
    count_metric metric = COUNT_METRIC_LEVENSHTEIN;
    int k = -1;

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--targets") == 0 && i + 1 < argc) {
            targets_path = argv[++i];
        } else if (strcmp(argv[i], "--reads") == 0 && i + 1 < argc) {
            reads_path = argv[++i];
        } else if (strcmp(argv[i], "--target-start") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_start) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--target-length") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &target_len) != 0 || target_len == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) {
            if (parse_int_value(argv[++i], &k) != 0 || (k != 0 && k != 1)) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--indel-window") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &indel_window) != 0 || indel_window > 1) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--metric") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "hamming") == 0) {
                metric = COUNT_METRIC_HAMMING;
            } else if (strcmp(value, "levenshtein") == 0) {
                metric = COUNT_METRIC_LEVENSHTEIN;
            } else {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--auto-offset") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &auto_offset) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--auto-offset-sample") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &auto_offset_sample) != 0 || auto_offset_sample == 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--offset-mode") == 0 && i + 1 < argc) {
            const char *value = argv[++i];
            if (strcmp(value, "best") == 0) {
                offsets_mode = OFFSET_MODE_BEST;
            } else if (strcmp(value, "multi") == 0) {
                offsets_mode = OFFSET_MODE_MULTI;
            } else {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--offset-min-fraction") == 0 && i + 1 < argc) {
            if (parse_double_value(argv[++i], &offset_min_fraction) != 0 ||
                offset_min_fraction < 0.0 || offset_min_fraction > 1.0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--oracle") == 0 && i + 1 < argc) {
            oracle = argv[++i];
        } else if (strcmp(argv[i], "--sample") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &sample_limit) != 0) {
                usage(argv0);
                return 2;
            }
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            if (parse_size_value(argv[++i], &threads) != 0 || threads == 0) {
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
    if (metric == COUNT_METRIC_HAMMING && indel_window != 0) {
        fprintf(stderr, "--indel-window is only valid with --metric levenshtein\n");
        return 2;
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
        if (sscanf(argv[2], "%d", &k) != 1) {
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
        return run_count(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "inspect-unmatched") == 0) {
        return run_inspect_unmatched(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "audit") == 0 || strcmp(argv[1], "audit-targets") == 0) {
        return run_audit(argv[0], argc, argv);
    }

    if (strcmp(argv[1], "validate") == 0) {
        return run_validate(argv[0], argc, argv);
    }

    usage(argv[0]);
    return 2;
}
