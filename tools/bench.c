#define _POSIX_C_SOURCE 199309L

#include "qdalign.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static uint64_t rng_state = 0x123456789abcdef0ULL;

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

typedef int (*dist_fn)(const char *, size_t, const char *, size_t);
typedef int (*leq_fn)(const char *, size_t, const char *, size_t, int);

static double run_case(const char *name, dist_fn fn, size_t len,
                       unsigned err_per_thousand, size_t calls) {
    const size_t pair_count = 1024;
    char *a = (char *)malloc(pair_count * (len + 1));
    char *b = (char *)malloc(pair_count * (len + 1));
    if (!a || !b) exit(1);

    for (size_t i = 0; i < pair_count; ++i) {
        char *pa = a + i * (len + 1);
        char *pb = b + i * (len + 1);
        rand_seq(pa, len);
        mutate_seq(pa, pb, len, err_per_thousand);
    }

    volatile int checksum = 0;
    double start = seconds_now();
    for (size_t i = 0; i < calls; ++i) {
        size_t idx = i & (pair_count - 1);
        checksum += fn(a + idx * (len + 1), len, b + idx * (len + 1), len);
    }
    double elapsed = seconds_now() - start;
    double calls_per_sec = (double)calls / elapsed;
    printf("%-10s len=%-5zu err=%-4.1f%% calls=%-9zu calls/s=%12.1f ns/call=%8.1f checksum=%d\n",
           name, len, (double)err_per_thousand / 10.0, calls,
           calls_per_sec, 1e9 / calls_per_sec, checksum);

    free(a);
    free(b);
    return elapsed;
}

static double run_leq_case(const char *name, leq_fn fn, size_t len, int k,
                           unsigned err_per_thousand, size_t calls) {
    const size_t pair_count = 1024;
    char *a = (char *)malloc(pair_count * (len + 1));
    char *b = (char *)malloc(pair_count * (len + 1));
    if (!a || !b) exit(1);

    for (size_t i = 0; i < pair_count; ++i) {
        char *pa = a + i * (len + 1);
        char *pb = b + i * (len + 1);
        rand_seq(pa, len);
        mutate_seq(pa, pb, len, err_per_thousand);
    }

    volatile int checksum = 0;
    double start = seconds_now();
    for (size_t i = 0; i < calls; ++i) {
        size_t idx = i & (pair_count - 1);
        checksum += fn(a + idx * (len + 1), len, b + idx * (len + 1), len, k);
    }
    double elapsed = seconds_now() - start;
    double calls_per_sec = (double)calls / elapsed;
    printf("%-10s len=%-5zu k=%-2d err=%-4.1f%% calls=%-9zu calls/s=%12.1f ns/call=%8.1f checksum=%d\n",
           name, len, k, (double)err_per_thousand / 10.0, calls,
           calls_per_sec, 1e9 / calls_per_sec, checksum);

    free(a);
    free(b);
    return elapsed;
}

static int dp_wrapper(const char *a, size_t a_len, const char *b, size_t b_len) {
    return qdaln_edit_distance_dp(a, a_len, b, b_len);
}

static int fast_wrapper(const char *a, size_t a_len, const char *b, size_t b_len) {
    return qdaln_edit_distance(a, a_len, b, b_len);
}

static int leq_wrapper(const char *a, size_t a_len, const char *b, size_t b_len, int k) {
    return qdaln_edit_distance_leq(a, a_len, b, b_len, k);
}

int main(void) {
    const unsigned errs[] = {0, 10, 30, 100};
    const size_t lens[] = {16, 32, 64, 128, 256};
    const size_t leq_lens[] = {16, 24, 32, 64};
    const int ks[] = {0, 1, 2, 3};

    printf("DotMatch microbench\n");
    printf("Note: DP is the slow correctness oracle; fast uses Myers64 when one input <=64.\n");
    printf("Pairs are pre-generated outside the timed loop.\n\n");

    for (size_t li = 0; li < sizeof(lens) / sizeof(lens[0]); ++li) {
        for (size_t ei = 0; ei < sizeof(errs) / sizeof(errs[0]); ++ei) {
            size_t len = lens[li];
            size_t fast_calls = len <= 64 ? 500000 : 10000;
            size_t dp_calls = len <= 64 ? 10000 : 1000;
            run_case("fast", fast_wrapper, len, errs[ei], fast_calls);
            run_case("dp", dp_wrapper, len, errs[ei], dp_calls);
            printf("\n");
        }
    }

    printf("threshold queries\n");
    printf("leq returns whether edit distance <= k.\n\n");

    for (size_t li = 0; li < sizeof(leq_lens) / sizeof(leq_lens[0]); ++li) {
        for (size_t ki = 0; ki < sizeof(ks) / sizeof(ks[0]); ++ki) {
            for (size_t ei = 0; ei < sizeof(errs) / sizeof(errs[0]); ++ei) {
                run_leq_case("leq", leq_wrapper, leq_lens[li], ks[ki], errs[ei], 500000);
            }
            printf("\n");
        }
    }
    return 0;
}
