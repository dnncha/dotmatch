#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE
#endif
#define _POSIX_C_SOURCE 199309L

#include "qdalign.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

static uint64_t rng_state = UINT64_C(0x13198a2e03707344);

static uint64_t next_random(void) {
    uint64_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng_state = x;
    return x;
}

static void random_sequence(char *sequence, size_t length) {
    static const char bases[] = "ACGT";
    for (size_t i = 0; i < length; ++i) sequence[i] = bases[next_random() & 3U];
    sequence[length] = '\0';
}

static void one_substitution(const char *source, char *destination, size_t length) {
    static const char bases[] = "ACGT";
    memcpy(destination, source, length + 1U);
    if (length == 0) return;
    size_t position = (size_t)(next_random() % length);
    char replacement = source[position];
    while (replacement == source[position]) replacement = bases[next_random() & 3U];
    destination[position] = replacement;
}

static double monotonic_seconds(void) {
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return 0.0;
    return (double)now.tv_sec + (double)now.tv_nsec / 1e9;
}

static long peak_rss_kb(void) {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) return -1;
#ifdef __APPLE__
    return (long)(usage.ru_maxrss / 1024);
#else
    return (long)usage.ru_maxrss;
#endif
}

static int parse_size(const char *text, size_t *value) {
    char *end = NULL;
    unsigned long parsed = strtoul(text, &end, 10);
    if (end == text || *end != '\0' || parsed == 0) return -1;
    *value = (size_t)parsed;
    return 0;
}

static int parse_int(const char *text, int *value) {
    char *end = NULL;
    long parsed = strtol(text, &end, 10);
    if (end == text || *end != '\0' || parsed < 0 || parsed > 3) return -1;
    *value = (int)parsed;
    return 0;
}

static long checksum_results(const qdaln_match_result *results, size_t count) {
    long checksum = 0;
    for (size_t i = 0; i < count; ++i) {
        checksum += (long)(results[i].target_index + 1) * 17L;
        checksum += (long)(results[i].best_distance + 1) * 31L;
        checksum += (long)(results[i].second_best_distance + 1) * 13L;
        checksum += (long)results[i].match_count * 7L;
        checksum += (long)results[i].status * 43L;
    }
    return checksum;
}

static int same_result(qdaln_match_result lhs, qdaln_match_result rhs) {
    return lhs.target_index == rhs.target_index &&
           lhs.best_distance == rhs.best_distance &&
           lhs.second_best_distance == rhs.second_best_distance &&
           lhs.match_count == rhs.match_count &&
           lhs.status == rhs.status;
}

int main(int argc, char **argv) {
    size_t n_reads = 1000000;
    size_t n_targets = 4096;
    size_t length = 20;
    int k = 1;
    size_t repeats = 5;
    if (argc > 1 && parse_size(argv[1], &n_reads) != 0) goto usage;
    if (argc > 2 && parse_size(argv[2], &n_targets) != 0) goto usage;
    if (argc > 3 && parse_size(argv[3], &length) != 0) goto usage;
    if (argc > 4 && parse_int(argv[4], &k) != 0) goto usage;
    if (argc > 5 && parse_size(argv[5], &repeats) != 0) goto usage;
    if (argc > 6 || length == 0 || length > 32) goto usage;

    char *target_storage = (char *)calloc(n_targets, length + 1U);
    char *read_storage = (char *)calloc(n_reads, length + 1U);
    const char **targets = (const char **)calloc(n_targets, sizeof(*targets));
    const char **reads = (const char **)calloc(n_reads, sizeof(*reads));
    size_t *target_lens = (size_t *)calloc(n_targets, sizeof(*target_lens));
    size_t *read_lens = (size_t *)calloc(n_reads, sizeof(*read_lens));
    qdaln_match_result *results = (qdaln_match_result *)calloc(n_reads, sizeof(*results));
    if (target_storage == NULL || read_storage == NULL || targets == NULL || reads == NULL ||
        target_lens == NULL || read_lens == NULL || results == NULL) {
        fprintf(stderr, "allocation failed\n");
        return 1;
    }

    for (size_t i = 0; i < n_targets; ++i) {
        targets[i] = target_storage + i * (length + 1U);
        target_lens[i] = length;
        random_sequence((char *)targets[i], length);
    }
    for (size_t i = 0; i < n_reads; ++i) {
        size_t target_index = i % n_targets;
        char *read = read_storage + i * (length + 1U);
        if ((i & 3U) == 0U) {
            memcpy(read, targets[target_index], length + 1U);
        } else if ((i & 3U) == 1U) {
            one_substitution(targets[target_index], read, length);
        } else {
            random_sequence(read, length);
        }
        reads[i] = read;
        read_lens[i] = length;
    }

    qdaln_index *index = qdaln_index_build(targets, target_lens, n_targets);
    if (index == NULL) {
        fprintf(stderr, "index build failed\n");
        return 1;
    }

    size_t validation_reads = n_reads < 256U ? n_reads : 256U;
    qdaln_match_result *scan = (qdaln_match_result *)calloc(validation_reads, sizeof(*scan));
    qdaln_match_result *indexed = (qdaln_match_result *)calloc(validation_reads, sizeof(*indexed));
    if (scan == NULL || indexed == NULL ||
        qdaln_match_many(reads, read_lens, validation_reads, targets, target_lens, n_targets, k, scan) != 0 ||
        qdaln_index_assign_hamming_stats(index, reads, read_lens, validation_reads, k, indexed, NULL) != 0) {
        fprintf(stderr, "validation setup failed\n");
        qdaln_index_free(index);
        return 1;
    }
    for (size_t i = 0; i < validation_reads; ++i) {
        if (!same_result(scan[i], indexed[i])) {
            fprintf(stderr, "output mismatch at validation read %zu\n", i);
            qdaln_index_free(index);
            return 1;
        }
    }
    free(scan);
    free(indexed);

    printf("run,n_reads,n_targets,length,k,seconds,reads_per_sec,candidates_per_read,verified_per_read,peak_rss_kb,checksum\n");
    for (size_t repeat = 0; repeat < repeats; ++repeat) {
        qdaln_index_stats stats = {0, 0};
        double start = monotonic_seconds();
        if (qdaln_index_assign_hamming_stats(index, reads, read_lens, n_reads, k, results, &stats) != 0) {
            fprintf(stderr, "indexed assignment failed on repeat %zu\n", repeat);
            qdaln_index_free(index);
            return 1;
        }
        double elapsed = monotonic_seconds() - start;
        if (elapsed <= 0.0) {
            fprintf(stderr, "invalid timer result\n");
            qdaln_index_free(index);
            return 1;
        }
        printf("%zu,%zu,%zu,%zu,%d,%.9f,%.3f,%.3f,%.3f,%ld,%ld\n",
               repeat, n_reads, n_targets, length, k, elapsed, (double)n_reads / elapsed,
               (double)stats.candidates_considered / (double)n_reads,
               (double)stats.candidates_verified / (double)n_reads,
               peak_rss_kb(), checksum_results(results, n_reads));
    }

    qdaln_index_free(index);
    free(target_storage);
    free(read_storage);
    free(targets);
    free(reads);
    free(target_lens);
    free(read_lens);
    free(results);
    return 0;

usage:
    fprintf(stderr, "Usage: %s [n_reads] [n_targets] [length<=32] [k=0..3] [repeats]\n", argv[0]);
    return 2;
}
