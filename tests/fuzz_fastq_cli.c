#define main dotmatch_cli_main
#include "../src/qda.c"
#undef main

#include <stdarg.h>

typedef struct fuzz_rng {
    uint64_t state;
} fuzz_rng;

static uint64_t fuzz_next(fuzz_rng *rng) {
    uint64_t x = rng->state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng->state = x;
    return x;
}

static size_t fuzz_below(fuzz_rng *rng, size_t n) {
    return n == 0 ? 0 : (size_t)(fuzz_next(rng) % n);
}

static int fuzz_chance(fuzz_rng *rng, unsigned percent) {
    return fuzz_below(rng, 100) < percent;
}

typedef struct fuzz_buf {
    char *data;
    size_t len;
    size_t cap;
} fuzz_buf;

static void fuzz_put(fuzz_buf *b, const char *s, size_t n) {
    if (n == 0 && b->data != NULL) return;
    if (b->len + n + 1 > b->cap) {
        size_t cap = b->cap == 0 ? 4096 : b->cap;
        while (cap < b->len + n + 1) cap *= 2;
        char *next = (char *)realloc(b->data, cap);
        if (next == NULL) abort();
        b->data = next;
        b->cap = cap;
    }
    if (n > 0) memcpy(b->data + b->len, s, n);
    b->len += n;
    b->data[b->len] = '\0';
}

static void fuzz_putc(fuzz_buf *b, char c) {
    fuzz_put(b, &c, 1);
}

static const char *const fuzz_targets[] = {
    "ACGTACGTAC", "ACGTACGTAA", "TTGCAAGGCT", "GGGGCCCCAA", "CATGCATGCA", "CATGCTTGCA", "CATGCNTGCA",
};
static const size_t fuzz_target_count = sizeof(fuzz_targets) / sizeof(fuzz_targets[0]);
static const size_t fuzz_target_len = 10;
static const size_t fuzz_target_start = 3;

static char fuzz_base(fuzz_rng *rng) {
    static const char alphabet[] = "ACGTACGTACGTNacgtRY";
    return alphabet[fuzz_below(rng, fuzz_chance(rng, 90) ? 4 : sizeof(alphabet) - 1)];
}

static void fuzz_sequence(fuzz_rng *rng, fuzz_buf *seq) {
    size_t len = fuzz_target_start + fuzz_target_len + 3;
    if (fuzz_chance(rng, 8)) len = fuzz_below(rng, len + 2);
    if (fuzz_chance(rng, 2)) len = 300 + fuzz_below(rng, 9000);
    for (size_t i = 0; i < len; ++i) fuzz_putc(seq, fuzz_base(rng));
    if (seq->len >= fuzz_target_start + fuzz_target_len && fuzz_chance(rng, 70)) {
        const char *t = fuzz_targets[fuzz_below(rng, fuzz_target_count)];
        memcpy(seq->data + fuzz_target_start, t, fuzz_target_len);
        size_t edits = fuzz_below(rng, 3);
        for (size_t e = 0; e < edits; ++e) {
            seq->data[fuzz_target_start + fuzz_below(rng, fuzz_target_len)] = "ACGT"[fuzz_below(rng, 4)];
        }
        if (fuzz_chance(rng, 10)) seq->data[fuzz_target_start + fuzz_below(rng, fuzz_target_len + 1)] = 'N';
        if (fuzz_chance(rng, 5)) {
            size_t at = fuzz_target_start + fuzz_below(rng, fuzz_target_len);
            if (seq->data[at] >= 'A' && seq->data[at] <= 'Z') seq->data[at] = (char)(seq->data[at] + ('a' - 'A'));
        }
        if (fuzz_chance(rng, 15)) {
            size_t at = fuzz_target_start + fuzz_below(rng, fuzz_target_len);
            memmove(seq->data + at, seq->data + at + 1, seq->len - at - 1);
            seq->data[seq->len - 1] = 'A';
        }
    }
}

enum fuzz_fault {
    FUZZ_CLEAN,
    FUZZ_STRAY_CR,
    FUZZ_NO_ID,
    FUZZ_LONG_HEADER,
    FUZZ_BAD_PLUS,
    FUZZ_PLUS_ID_MISMATCH,
    FUZZ_QUAL_LENGTH,
    FUZZ_QUAL_SPACE,
    FUZZ_TRUNCATE,
    FUZZ_BYTE,
    FUZZ_NUL_BOTH,
    FUZZ_BLANK_TAIL,
    FUZZ_NO_FINAL_NEWLINE,
    FUZZ_FAULT_COUNT
};

static void fuzz_newline(fuzz_buf *out, int crlf) {
    if (crlf) fuzz_putc(out, '\r');
    fuzz_putc(out, '\n');
}

static void fuzz_fastq(fuzz_rng *rng, fuzz_buf *out, size_t records) {
    int crlf = fuzz_chance(rng, 15);
    int fault = fuzz_chance(rng, 50) ? FUZZ_CLEAN : (int)(1 + fuzz_below(rng, FUZZ_FAULT_COUNT - 1));
    size_t fault_record = fuzz_below(rng, records);
    fuzz_buf seq = {0};
    char id[64];
    for (size_t r = 0; r < records; ++r) {
        int here = r == fault_record ? fault : FUZZ_CLEAN;
        seq.len = 0;
        fuzz_sequence(rng, &seq);
        snprintf(id, sizeof(id), "read%zu", r);
        fuzz_putc(out, '@');
        if (here != FUZZ_NO_ID) fuzz_put(out, id, strlen(id));
        if (fuzz_chance(rng, 20)) fuzz_put(out, " 1:N:0:1", 8);
        if (here == FUZZ_LONG_HEADER) {
            size_t extra = 200 + fuzz_below(rng, 9000);
            fuzz_putc(out, ' ');
            for (size_t i = 0; i < extra; ++i) fuzz_putc(out, 'x');
        }
        fuzz_newline(out, crlf || here == FUZZ_STRAY_CR);
        size_t seq_at = out->len;
        fuzz_put(out, seq.data, seq.len);
        fuzz_newline(out, crlf);
        fuzz_putc(out, here == FUZZ_BAD_PLUS ? '-' : '+');
        if (fuzz_chance(rng, 10) || here == FUZZ_PLUS_ID_MISMATCH) fuzz_put(out, id, strlen(id));
        if (here == FUZZ_PLUS_ID_MISMATCH) fuzz_put(out, "x", 1);
        fuzz_newline(out, crlf);
        size_t qlen = seq.len;
        if (here == FUZZ_QUAL_LENGTH) qlen = fuzz_below(rng, seq.len + 2);
        size_t qual_at = out->len;
        for (size_t i = 0; i < qlen; ++i) fuzz_putc(out, (char)('#' + fuzz_below(rng, 40)));
        if (here == FUZZ_QUAL_SPACE && qlen > 0) out->data[qual_at + fuzz_below(rng, qlen)] = ' ';
        if (here == FUZZ_NUL_BOTH && seq.len > 0 && qlen == seq.len) {
            size_t at = fuzz_below(rng, seq.len);
            out->data[seq_at + at] = '\0';
            out->data[qual_at + at] = '\0';
        }
        int last = r + 1 == records;
        if (!(last && fault == FUZZ_NO_FINAL_NEWLINE)) fuzz_newline(out, crlf);
    }
    if (fault == FUZZ_TRUNCATE && out->len > 0) {
        size_t cut = fuzz_below(rng, out->len);
        out->len = cut;
        out->data[cut] = '\0';
    }
    if (fault == FUZZ_BYTE && out->len > 0) out->data[fuzz_below(rng, out->len)] = (char)fuzz_below(rng, 256);
    if (fault == FUZZ_BLANK_TAIL) fuzz_put(out, "\n", 1);
    free(seq.data);
}

static char fuzz_dir[512];
static const char *fuzz_targets_file = "targets.tsv";
static int fuzz_report_fd = 2;

static void fuzz_path(char *out, size_t cap, const char *name) {
    snprintf(out, cap, "%s/%s", fuzz_dir, name);
}

static void fuzz_write_file(const char *name, const char *data, size_t len, int gz) {
    char path[600];
    fuzz_path(path, sizeof(path), name);
    if (gz) {
        gzFile f = gzopen(path, "wb1");
        if (f == NULL) abort();
        if (len > 0 && gzwrite(f, data, (unsigned)len) != (int)len) abort();
        gzclose(f);
        return;
    }
    FILE *f = fopen(path, "wb");
    if (f == NULL) abort();
    if (len > 0 && fwrite(data, 1, len, f) != len) abort();
    fclose(f);
}

static char *fuzz_slurp(const char *name, size_t *len_out) {
    char path[600];
    fuzz_path(path, sizeof(path), name);
    FILE *f = fopen(path, "rb");
    *len_out = 0;
    if (f == NULL) return NULL;
    fuzz_buf b = {0};
    char chunk[4096];
    size_t n;
    while ((n = fread(chunk, 1, sizeof(chunk), f)) > 0) fuzz_put(&b, chunk, n);
    fclose(f);
    *len_out = b.len;
    if (b.data == NULL) fuzz_put(&b, "", 0);
    return b.data;
}

static void fuzz_unlink(const char *name) {
    char path[600];
    fuzz_path(path, sizeof(path), name);
    unlink(path);
}

static int fuzz_run(int argc_hint, ...) {
    char *argv[64];
    int argc = 0;
    va_list ap;
    va_start(ap, argc_hint);
    const char *arg;
    while ((arg = va_arg(ap, const char *)) != NULL && argc < 63) argv[argc++] = (char *)arg;
    va_end(ap);
    argv[argc] = NULL;
    (void)argc_hint;
    fflush(stdout);
    return dotmatch_cli_main(argc, argv);
}

static void fuzz_fail(uint64_t seed, const char *what, const fuzz_buf *fastq) {
    dprintf(fuzz_report_fd, "fuzz_fastq_cli: seed %llu: %s\n", (unsigned long long)seed, what);
    char path[600];
    fuzz_path(path, sizeof(path), "failing.fastq");
    FILE *f = fopen(path, "wb");
    if (f != NULL) {
        fwrite(fastq->data, 1, fastq->len, f);
        fclose(f);
        dprintf(fuzz_report_fd, "fuzz_fastq_cli: input kept at %s\n", path);
    }
    exit(1);
}

static int fuzz_same_file(const char *a, const char *b) {
    size_t la = 0;
    size_t lb = 0;
    char *da = fuzz_slurp(a, &la);
    char *db = fuzz_slurp(b, &lb);
    int same = (da == NULL && db == NULL) || (da != NULL && db != NULL && la == lb && memcmp(da, db, la) == 0);
    free(da);
    free(db);
    return same;
}

static int fuzz_count(const char *reads, const char *out, const char *k, const char *metric, const char *index) {
    char targets_path[600];
    char reads_path[600];
    char out_path[600];
    fuzz_path(targets_path, sizeof(targets_path), fuzz_targets_file);
    fuzz_path(reads_path, sizeof(reads_path), reads);
    fuzz_path(out_path, sizeof(out_path), out);
    unlink(out_path);
    return fuzz_run(0, "dotmatch", "count", "--targets", targets_path, "--reads", reads_path, "--sample-label", "s",
                    "--target-start", "3", "--target-length", "10", "--k", k, "--metric", metric, "--hamming-index",
                    index, "--offset-mode", "best", "--out", out_path, (const char *)NULL);
}

static int fuzz_count_engine(const char *out, const char *k, const char *indel, const char *policy, int generic) {
    char targets_path[600];
    char reads_path[600];
    char out_path[600];
    char assignments_path[600];
    fuzz_path(targets_path, sizeof(targets_path), fuzz_targets_file);
    fuzz_path(reads_path, sizeof(reads_path), "reads.fastq");
    fuzz_path(out_path, sizeof(out_path), out);
    fuzz_path(assignments_path, sizeof(assignments_path), "engine_assignments.tsv");
    unlink(out_path);
    if (generic) {
        return fuzz_run(0, "dotmatch", "count", "--targets", targets_path, "--reads", reads_path, "--sample-label", "s",
                        "--target-start", "3", "--target-length", "10", "--k", k, "--metric", "levenshtein",
                        "--indel-window", indel, "--ambiguity-policy", policy, "--offset-mode", "best", "--out",
                        out_path, "--assignments", assignments_path, (const char *)NULL);
    }
    return fuzz_run(0, "dotmatch", "count", "--targets", targets_path, "--reads", reads_path, "--sample-label", "s",
                    "--target-start", "3", "--target-length", "10", "--k", k, "--metric", "levenshtein",
                    "--indel-window", indel, "--ambiguity-policy", policy, "--offset-mode", "best", "--out", out_path,
                    (const char *)NULL);
}

static int fuzz_assign(const char *reads, const char *out, const char *k) {
    char targets_path[600];
    char reads_path[600];
    char out_path[600];
    fuzz_path(targets_path, sizeof(targets_path), fuzz_targets_file);
    fuzz_path(reads_path, sizeof(reads_path), reads);
    fuzz_path(out_path, sizeof(out_path), out);
    unlink(out_path);
    return fuzz_run(0, "dotmatch", "fastq-assign", "--barcodes", targets_path, "--reads", reads_path,
                    "--barcode-start", "3", "--barcode-length", "10", "--k", k, "--out", out_path,
                    (const char *)NULL);
}

static int fuzz_validate(const char *reads, const char *k, const char *metric) {
    char targets_path[600];
    char reads_path[600];
    fuzz_path(targets_path, sizeof(targets_path), fuzz_targets_file);
    fuzz_path(reads_path, sizeof(reads_path), reads);
    return fuzz_run(0, "dotmatch", "validate", "--targets", targets_path, "--reads", reads_path, "--target-start", "3",
                    "--target-length", "10", "--k", k, "--metric", metric, "--threads", "1", "--oracle", "scan",
                    (const char *)NULL);
}

static void fuzz_one(uint64_t seed) {
    fuzz_rng rng = {seed * 0x9E3779B97F4A7C15ull + 1};
    fuzz_buf fastq = {0};
    fuzz_fastq(&rng, &fastq, 1 + fuzz_below(&rng, 40));
    fuzz_write_file("reads.fastq", fastq.data, fastq.len, 0);
    fuzz_write_file("reads.fastq.gz", fastq.data, fastq.len, 1);
    char what[256];

    static const char *const metrics[] = {"hamming", "levenshtein"};
    static const char *const ks[] = {"0", "1"};
    fuzz_targets_file = fuzz_chance(&rng, 50) ? "targets.tsv" : "targets_acgt.tsv";
    const char *metric = metrics[fuzz_below(&rng, 2)];
    const char *k = ks[fuzz_below(&rng, 2)];

    int plain_rc = fuzz_count("reads.fastq", "plain.tsv", k, metric, "auto");
    int gz_rc = fuzz_count("reads.fastq.gz", "gz.tsv", k, metric, "auto");
    if ((plain_rc == 0) != (gz_rc == 0)) {
        snprintf(what, sizeof(what), "count %s k=%s: plain rc %d but gzip rc %d", metric, k, plain_rc, gz_rc);
        fuzz_fail(seed, what, &fastq);
    }
    if (plain_rc == 0 && !fuzz_same_file("plain.tsv", "gz.tsv")) {
        snprintf(what, sizeof(what), "count %s k=%s: plain and gzip counts differ", metric, k);
        fuzz_fail(seed, what, &fastq);
    }
    if (plain_rc == 0 && strcmp(metric, "hamming") == 0) {
        static const char *const indexes[] = {"query", "precompute"};
        for (size_t i = 0; i < 2; ++i) {
            int rc = fuzz_count("reads.fastq", "index.tsv", k, metric, indexes[i]);
            if (rc != 0 || !fuzz_same_file("plain.tsv", "index.tsv")) {
                snprintf(what, sizeof(what), "count hamming k=%s: --hamming-index %s rc %d or counts differ from auto", k,
                         indexes[i], rc);
                fuzz_fail(seed, what, &fastq);
            }
        }
    }

    if (plain_rc == 0) {
        static const char *const indels[] = {"0", "1"};
        static const char *const policies[] = {"radius", "best"};
        const char *indel = strcmp(k, "0") == 0 ? "0" : indels[fuzz_below(&rng, 2)];
        const char *policy = policies[fuzz_below(&rng, 2)];
        int fast_rc = fuzz_count_engine("engine_fast.tsv", k, indel, policy, 0);
        int generic_rc = fuzz_count_engine("engine_generic.tsv", k, indel, policy, 1);
        if (fast_rc != 0 || generic_rc != 0 || !fuzz_same_file("engine_fast.tsv", "engine_generic.tsv")) {
            snprintf(what, sizeof(what),
                     "count levenshtein k=%s --indel-window %s --ambiguity-policy %s: lookup engine (rc %d) and "
                     "generic engine with --assignments (rc %d) give different counts",
                     k, indel, policy, fast_rc, generic_rc);
            fuzz_fail(seed, what, &fastq);
        }
    }

    int assign_plain = fuzz_assign("reads.fastq", "assign_plain.tsv", k);
    int assign_gz = fuzz_assign("reads.fastq.gz", "assign_gz.tsv", k);
    if ((assign_plain == 0) != (assign_gz == 0)) {
        snprintf(what, sizeof(what), "fastq-assign k=%s: plain rc %d but gzip rc %d", k, assign_plain, assign_gz);
        fuzz_fail(seed, what, &fastq);
    }
    if (assign_plain == 0 && !fuzz_same_file("assign_plain.tsv", "assign_gz.tsv")) {
        snprintf(what, sizeof(what), "fastq-assign k=%s: plain and gzip assignments differ", k);
        fuzz_fail(seed, what, &fastq);
    }
    if ((assign_plain == 0) != (plain_rc == 0)) {
        snprintf(what, sizeof(what), "count rc %d but fastq-assign rc %d on the same reads", plain_rc, assign_plain);
        fuzz_fail(seed, what, &fastq);
    }

    if (plain_rc == 0) {
        int rc = fuzz_validate("reads.fastq", k, metric);
        if (rc != 0) {
            snprintf(what, sizeof(what), "validate %s k=%s: lookup disagrees with scan oracle (rc %d)", metric, k, rc);
            fuzz_fail(seed, what, &fastq);
        }
    }
    free(fastq.data);
}

int main(int argc, char **argv) {
    uint64_t first = argc > 1 ? strtoull(argv[1], NULL, 10) : 1;
    uint64_t iterations = argc > 2 ? strtoull(argv[2], NULL, 10) : 300;
    const char *tmp = getenv("TMPDIR");
    snprintf(fuzz_dir, sizeof(fuzz_dir), "%s/dotmatch-fuzz-XXXXXX", tmp != NULL && *tmp ? tmp : "/tmp");
    if (mkdtemp(fuzz_dir) == NULL) {
        perror("mkdtemp");
        return 2;
    }
    for (int acgt_only = 0; acgt_only < 2; ++acgt_only) {
        fuzz_buf targets = {0};
        fuzz_put(&targets, "id\tseq\n", 7);
        for (size_t i = 0; i < fuzz_target_count; ++i) {
            if (acgt_only && strchr(fuzz_targets[i], 'N') != NULL) continue;
            char line[64];
            int n = snprintf(line, sizeof(line), "t%zu\t%s\n", i, fuzz_targets[i]);
            fuzz_put(&targets, line, (size_t)n);
        }
        fuzz_write_file(acgt_only ? "targets_acgt.tsv" : "targets.tsv", targets.data, targets.len, 0);
        free(targets.data);
    }

    fuzz_report_fd = dup(2);
    int devnull = open("/dev/null", O_WRONLY);
    if (devnull < 0 || fuzz_report_fd < 0) return 2;
    dup2(devnull, 1);
    if (getenv("FUZZ_SHOW_CLI_STDERR") == NULL) dup2(devnull, 2);
    close(devnull);

    for (uint64_t s = first; s < first + iterations; ++s) fuzz_one(s);

    static const char *const names[] = {"reads.fastq", "reads.fastq.gz", "plain.tsv", "gz.tsv", "index.tsv",
                                        "assign_plain.tsv", "assign_gz.tsv", "engine_fast.tsv",
                                        "engine_generic.tsv", "engine_assignments.tsv", "targets.tsv",
                                        "targets_acgt.tsv"};
    for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); ++i) fuzz_unlink(names[i]);
    rmdir(fuzz_dir);
    dprintf(fuzz_report_fd, "fuzz_fastq_cli: %llu seeds from %llu passed\n", (unsigned long long)iterations,
            (unsigned long long)first);
    return 0;
}
