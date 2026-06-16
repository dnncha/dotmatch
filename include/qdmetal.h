#ifndef QDMETAL_H
#define QDMETAL_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define QDMETAL_STATUS_UNAVAILABLE 0
#define QDMETAL_STATUS_OK 1

#define QDMETAL_PATH_BRUTE_FORCE "brute_force"
#define QDMETAL_PATH_SEED_INDEX "seed_index"

typedef struct qdmetal_match_result {
    int target_index;
    int best_distance;
    int second_best_distance;
    int match_count;
    int status;
} qdmetal_match_result;

typedef struct qdmetal_assign_stats {
    size_t candidates_considered;
    size_t candidates_verified;
    const char *path;
    const char *device_name;
} qdmetal_assign_stats;

int qdmetal_available(void);
const char *qdmetal_device_name(void);

/*
 * Assign packed A/C/G/T read windows to packed targets under Hamming distance <= k.
 * Uses a brute-force Metal kernel for modest target counts and a CPU-seeded
 * candidate reduction for larger libraries. Semantics match best-distance Hamming
 * assignment (k=0 exact, k=1 one mismatch).
 *
 * Returns 0 on success, -1 on invalid input or Metal runtime failure.
 */
int qdmetal_hamming_assign(const uint64_t *read_codes, size_t n_reads, const uint64_t *target_codes,
                           size_t n_targets, size_t len, int k, qdmetal_match_result *results,
                           qdmetal_assign_stats *stats);

#ifdef __cplusplus
}
#endif

#endif /* QDMETAL_H */