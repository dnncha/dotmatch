#ifndef QDALIGN_H
#define QDALIGN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define QDALN_VERSION "0.2.1"
#define QDALN_ALPHABET_POLICY "literal-byte; A/C/G/T/N/IUPAC symbols are ordinary byte symbols; no wildcard expansion"

enum qdaln_match_status {
    QDALN_MATCH_INVALID = -1,
    QDALN_MATCH_NONE = 0,
    QDALN_MATCH_UNIQUE = 1,
    QDALN_MATCH_AMBIGUOUS = 2
};

typedef struct qdaln_match_result {
    int target_index;
    int best_distance;
    int second_best_distance;
    int match_count;
    int status;
} qdaln_match_result;

enum qdaln_ambiguity_policy {
    QDALN_POLICY_BEST = 0,
    QDALN_POLICY_RADIUS = 1
};

enum qdaln_edit_class {
    QDALN_EDIT_INVALID = -1,
    QDALN_EDIT_NONE = 0,
    QDALN_EDIT_EXACT = 1,
    QDALN_EDIT_K1_SUB = 2,
    QDALN_EDIT_K1_INS = 3,
    QDALN_EDIT_K1_DEL = 4,
    QDALN_EDIT_K2 = 5,
    QDALN_EDIT_OTHER = 6
};

typedef struct qdaln_assignment_result {
    int target_index;
    int distance;
    int second_best_distance;
    int num_best_targets;
    int num_targets_within_radius;
    int status;
    int edit_class;
} qdaln_assignment_result;

typedef struct qdaln_index qdaln_index;

typedef struct qdaln_index_stats {
    size_t candidates_considered;
    size_t candidates_verified;
} qdaln_index_stats;

/*
 * Public alphabet contract for reproducible known-target assignment.
 * DotMatch compares bytes exactly: A/C/G/T, N, and IUPAC ambiguity symbols
 * are ordinary byte symbols. There is no wildcard expansion.
 */
const char *qdaln_alphabet_policy(void);

/*
 * Exact Levenshtein edit distance between two byte strings.
 * Costs: substitution=1, insertion=1, deletion=1.
 * For DNA this follows qdaln_alphabet_policy().
 *
 * Returns -1 on invalid input.
 */
int qdaln_edit_distance(const char *a, size_t a_len, const char *b, size_t b_len);

/*
 * Myers bit-vector kernel (generalized from 64-bit to multi-word for pattern_len > 64,
 * using arrays of uint64_t with portable carry handling for + and <<1 in the
 * bit-parallel steps; no __int128). Computes exact Levenshtein distance.
 * Used by qdaln_edit_distance and qdaln_edit_distance_leq for fast exact/leq
 * even on longer patterns (primers etc). Falls back to DP only for extreme lengths.
 */
int qdaln_edit_distance_myers64(const char *pattern, size_t pattern_len,
                                const char *text, size_t text_len);

/*
 * Slow but simple two-row dynamic-programming oracle.
 * Intended for tests and correctness comparisons, not production speed.
 */
int qdaln_edit_distance_dp(const char *a, size_t a_len, const char *b, size_t b_len);

/*
 * Returns 1 if edit distance <= k, 0 if > k, -1 on invalid input.
 * Uses a thresholded dynamic-programming kernel with early rejection.
 */
int qdaln_edit_distance_leq(const char *a, size_t a_len, const char *b, size_t b_len, int k);

/*
 * Assign each read to the best target within edit-distance threshold k.
 * Inputs are flat arrays of string pointers and lengths.
 *
 * For each read:
 * - UNIQUE: exactly one target has the best distance <= k.
 * - AMBIGUOUS: two or more targets share the best distance <= k.
 * - NONE: no target is within k.
 * - INVALID: invalid input for that read or call.
 *
 * Returns 0 on successful processing, -1 on invalid call-level input.
 */
int qdaln_match_many(const char *const *reads, const size_t *read_lens, size_t n_reads,
                     const char *const *targets, const size_t *target_lens, size_t n_targets,
                     int k, qdaln_match_result *results);

/*
 * Stable assignment contract for trusted known-target workflows.
 *
 * Policy:
 * - QDALN_POLICY_BEST: unique if exactly one target has the minimum distance.
 * - QDALN_POLICY_RADIUS: unique only if exactly one target is within radius k.
 *
 * edit_class describes the best assignment when status is UNIQUE or AMBIGUOUS.
 * Insertions/deletions are named from the observed read relative to the target:
 * K1_INS means the read has one extra base; K1_DEL means it is missing one base.
 */
int qdaln_assign_many(const char *const *reads, const size_t *read_lens, size_t n_reads,
                      const char *const *targets, const size_t *target_lens, size_t n_targets,
                      int k, int policy, qdaln_assignment_result *results);

/*
 * Build/free a reusable target index for repeated short-DNA assignment.
 * The index owns a private copy of the target strings.
 */
qdaln_index *qdaln_index_build(const char *const *targets, const size_t *target_lens, size_t n_targets);
void qdaln_index_free(qdaln_index *index);

/*
 * Assign reads against a reusable target index.
 * Results are identical to qdaln_match_many for the same targets.
 */
int qdaln_index_assign(const qdaln_index *index, const char *const *reads, const size_t *read_lens,
                       size_t n_reads, int k, qdaln_match_result *results);

int qdaln_index_assign_stats(const qdaln_index *index, const char *const *reads, const size_t *read_lens,
                             size_t n_reads, int k, qdaln_match_result *results,
                             qdaln_index_stats *stats);

/*
 * Count/status-oriented indexed assignment.
 *
 * Preserves target_index for UNIQUE, best_distance, and UNIQUE/AMBIGUOUS/NONE/
 * INVALID status. For AMBIGUOUS results, match_count and second_best_distance
 * may be lower-bound values because the search can stop once ambiguity at the
 * best distance is proven. Use qdaln_index_assign_stats when exact match_count
 * reporting is required.
 */
int qdaln_index_assign_status_stats(const qdaln_index *index, const char *const *reads,
                                    const size_t *read_lens, size_t n_reads, int k,
                                    qdaln_match_result *results, qdaln_index_stats *stats);

/*
 * Hamming-distance indexed assignment for equal-length known-target workflows.
 * Supports k=0..3. Unlike Levenshtein assignment, this does not consider
 * insertion/deletion neighbors, which is the fair mode for one-mismatch guide
 * counters.
 */
int qdaln_index_assign_hamming_stats(const qdaln_index *index, const char *const *reads,
                                     const size_t *read_lens, size_t n_reads, int k,
                                     qdaln_match_result *results, qdaln_index_stats *stats);

/*
 * Direct exact lookup for a single same-length known-target window.
 * This is the low-overhead path for exact counting/demultiplexing loops.
 * Results match qdaln_index_assign_hamming_stats(..., k=0, n_reads=1).
 */
int qdaln_index_lookup_exact_stats(const qdaln_index *index, const char *read, size_t read_len,
                                   qdaln_match_result *result, qdaln_index_stats *stats);

/*
 * Batched direct exact lookup. Results match
 * qdaln_index_assign_hamming_stats(..., k=0) while avoiding the general
 * assignment dispatcher and repeated public API setup.
 */
int qdaln_index_lookup_exact_many_stats(const qdaln_index *index, const char *const *reads,
                                        const size_t *read_lens, size_t n_reads,
                                        qdaln_match_result *results, qdaln_index_stats *stats);

/*
 * Batched ASCII-folding exact lookup for buffered FASTQ windows.
 * This is the exact-table path used by count-mode readers that keep raw ASCII
 * windows and want to avoid per-read output plumbing.
 */
int qdaln_index_lookup_exact_ascii_many_stats(const qdaln_index *index, const char *const *reads,
                                              const size_t *read_lens, size_t n_reads,
                                              qdaln_match_result *results, qdaln_index_stats *stats);

/*
 * Direct exact lookup that folds ASCII a/c/g/t to A/C/G/T while preserving
 * literal-byte behavior for other symbols through the non-encodable fallback.
 * This matches DotMatch FASTQ window normalization without requiring callers
 * to copy and uppercase the read window before lookup.
 */
int qdaln_index_lookup_exact_ascii_stats(const qdaln_index *index, const char *read, size_t read_len,
                                         qdaln_match_result *result, qdaln_index_stats *stats);

#ifdef __cplusplus
}
#endif

#endif /* QDALIGN_H */
