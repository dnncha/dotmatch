#define _POSIX_C_SOURCE 199309L

#include "qdalign.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static uint64_t rng_state = 0x243f6a8885a308d3ULL;

static uint64_t xorshift64(void) {
    uint64_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng_state = x;
    return x;
}

static char rand_base(void) {
    static const char dna[] = "ACGT";
    return dna[xorshift64() & 3ULL];
}

static void rand_seq(char *s, size_t n) {
    for (size_t i = 0; i < n; ++i) s[i] = rand_base();
    s[n] = '\0';
}

static void mutate_seq(const char *src, char *dst, size_t n, unsigned per_thousand) {
    static const char dna[] = "ACGT";
    for (size_t i = 0; i < n; ++i) {
        char c = src[i];
        if ((xorshift64() % 1000ULL) < per_thousand) {
            char nc = c;
            while (nc == c) nc = dna[xorshift64() & 3ULL];
            c = nc;
        }
        dst[i] = c;
    }
    dst[n] = '\0';
}

static double seconds_now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static int parse_size_arg(const char *s, size_t *out) {
    char *end = NULL;
    unsigned long v = strtoul(s, &end, 10);
    if (end == s || *end != '\0') return -1;
    *out = (size_t)v;
    return 0;
}

static long checksum_results(const qdaln_match_result *results, size_t n) {
    long checksum = 0;
    for (size_t i = 0; i < n; ++i) {
        checksum += (long)(results[i].target_index + 1) * 17L;
        checksum += (long)(results[i].best_distance + 1) * 31L;
        checksum += (long)results[i].status * 43L;
        checksum += (long)results[i].match_count * 7L;
    }
    return checksum;
}

static void scan_assign(const char *const *reads, const size_t *read_lens, size_t n_reads,
                        const char *const *targets, const size_t *target_lens, size_t n_targets,
                        int k, qdaln_match_result *results) {
    if (qdaln_match_many(reads, read_lens, n_reads, targets, target_lens, n_targets, k, results) != 0) exit(1);
}

static void naive_assign(const char *const *reads, const size_t *read_lens, size_t n_reads,
                         const char *const *targets, const size_t *target_lens, size_t n_targets,
                         int k, qdaln_match_result *results) {
    for (size_t i = 0; i < n_reads; ++i) {
        qdaln_match_many(&reads[i], &read_lens[i], 1, targets, target_lens, n_targets, k, &results[i]);
    }
}

static void run_case(size_t n_reads, size_t n_targets, size_t len, int k, unsigned err_per_thousand) {
    char *target_buf = (char *)malloc(n_targets * (len + 1));
    char *read_buf = (char *)malloc(n_reads * (len + 1));
    const char **targets = (const char **)malloc(n_targets * sizeof(char *));
    const char **reads = (const char **)malloc(n_reads * sizeof(char *));
    size_t *target_lens = (size_t *)malloc(n_targets * sizeof(size_t));
    size_t *read_lens = (size_t *)malloc(n_reads * sizeof(size_t));
    qdaln_match_result *results = (qdaln_match_result *)malloc(n_reads * sizeof(qdaln_match_result));
    if (!target_buf || !read_buf || !targets || !reads || !target_lens || !read_lens || !results) exit(1);

    for (size_t i = 0; i < n_targets; ++i) {
        char *target = target_buf + i * (len + 1);
        rand_seq(target, len);
        targets[i] = target;
        target_lens[i] = len;
    }
    for (size_t i = 0; i < n_reads; ++i) {
        size_t idx = (size_t)(xorshift64() % n_targets);
        char *read = read_buf + i * (len + 1);
        mutate_seq(targets[idx], read, len, err_per_thousand);
        reads[i] = read;
        read_lens[i] = len;
    }

    qdaln_index *index = qdaln_index_build(targets, target_lens, n_targets);
    if (index == NULL) exit(1);

    double start = seconds_now();
    qdaln_index_stats stats;
    if (qdaln_index_assign_stats(index, reads, read_lens, n_reads, k, results, &stats) != 0) exit(1);
    double elapsed = seconds_now() - start;
    long checksum = checksum_results(results, n_reads);
    double reads_per_sec = (double)n_reads / elapsed;
    printf("dotmatch_indexed,synthetic_barcode,%zu,%zu,%zu,%d,%.3f,%.6f,%.1f,%.1f,%ld\n",
           n_reads, n_targets, len, k, (double)err_per_thousand / 1000.0, elapsed,
           reads_per_sec, reads_per_sec * (double)n_targets, checksum);

    start = seconds_now();
    scan_assign(reads, read_lens, n_reads, targets, target_lens, n_targets, k, results);
    elapsed = seconds_now() - start;
    checksum = checksum_results(results, n_reads);
    reads_per_sec = (double)n_reads / elapsed;
    printf("dotmatch_scan,synthetic_barcode,%zu,%zu,%zu,%d,%.3f,%.6f,%.1f,%.1f,%ld\n",
           n_reads, n_targets, len, k, (double)err_per_thousand / 1000.0, elapsed,
           reads_per_sec, reads_per_sec * (double)n_targets, checksum);

    size_t naive_reads = n_reads < 5000 ? n_reads : 5000;
    start = seconds_now();
    naive_assign(reads, read_lens, naive_reads, targets, target_lens, n_targets, k, results);
    elapsed = seconds_now() - start;
    checksum = checksum_results(results, naive_reads);
    reads_per_sec = (double)naive_reads / elapsed;
    printf("dotmatch_naive,synthetic_barcode,%zu,%zu,%zu,%d,%.3f,%.6f,%.1f,%.1f,%ld\n",
           naive_reads, n_targets, len, k, (double)err_per_thousand / 1000.0, elapsed,
           reads_per_sec, reads_per_sec * (double)n_targets, checksum);

    free(target_buf);
    free(read_buf);
    free(targets);
    free(reads);
    free(target_lens);
    free(read_lens);
    free(results);
    qdaln_index_free(index);
}

int main(int argc, char **argv) {
    size_t n_reads = 100000;
    if (argc == 2 && parse_size_arg(argv[1], &n_reads) != 0) {
        fprintf(stderr, "Usage: %s [n_reads]\n", argv[0]);
        return 2;
    }

    const size_t lens[] = {8, 12, 16, 24, 32};
    const size_t target_counts[] = {96, 384, 737, 4096};
    const int ks[] = {0, 1, 2, 3};
    const unsigned errs[] = {0, 5, 10, 30};

    printf("tool,workload,n_reads,n_targets,len,k,err,seconds,reads_per_sec,assignments_per_sec,checksum\n");
    for (size_t li = 0; li < sizeof(lens) / sizeof(lens[0]); ++li) {
        for (size_t ti = 0; ti < sizeof(target_counts) / sizeof(target_counts[0]); ++ti) {
            for (size_t ki = 0; ki < sizeof(ks) / sizeof(ks[0]); ++ki) {
                for (size_t ei = 0; ei < sizeof(errs) / sizeof(errs[0]); ++ei) {
                    run_case(n_reads, target_counts[ti], lens[li], ks[ki], errs[ei]);
                }
            }
        }
    }
    return 0;
}
