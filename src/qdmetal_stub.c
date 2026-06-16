#include "qdmetal.h"

int qdmetal_available(void) {
    return 0;
}

const char *qdmetal_device_name(void) {
    return NULL;
}

int qdmetal_hamming_assign(const uint64_t *read_codes, size_t n_reads, const uint64_t *target_codes,
                           size_t n_targets, size_t len, int k, qdmetal_match_result *results,
                           qdmetal_assign_stats *stats) {
    (void)read_codes;
    (void)n_reads;
    (void)target_codes;
    (void)n_targets;
    (void)len;
    (void)k;
    (void)results;
    (void)stats;
    return -1;
}