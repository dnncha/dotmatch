#include "qdalign.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint64_t rng_state = 0x9e3779b97f4a7c15ULL;

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

static void check_pair(const char *a, const char *b) {
    int dp = qdaln_edit_distance_dp(a, strlen(a), b, strlen(b));
    int fast = qdaln_edit_distance(a, strlen(a), b, strlen(b));
    int fast_direct = qdaln_edit_distance_myers64(a, strlen(a), b, strlen(b));
    assert(dp == fast);
    if (strlen(a) <= 64) assert(dp == fast_direct);

    for (int k = 0; k <= 8; ++k) {
        int leq = qdaln_edit_distance_leq(a, strlen(a), b, strlen(b), k);
        assert(leq == (dp <= k ? 1 : 0));
    }
}

static void fixed_tests(void) {
    check_pair("", "");
    check_pair("A", "");
    check_pair("", "ACGT");
    check_pair("A", "A");
    check_pair("A", "C");
    check_pair("ACGT", "ACGT");
    check_pair("ACGT", "AGGT");
    check_pair("ACGT", "ACGTT");
    check_pair("GATTACA", "GCATGCU");
    check_pair("AAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAT");
    check_pair("ACGTACGTACGTACGT", "TGCATGCATGCATGCA");

    assert(qdaln_edit_distance_leq(NULL, 1, "A", 1, -1) == -1);
    assert(qdaln_edit_distance_leq("A", 1, NULL, 1, -1) == -1);
    assert(qdaln_edit_distance_leq("A", 1, "ACGT", 4, 2) == 0);
    assert(qdaln_edit_distance_leq("A", 1, "ACGT", 4, 3) == 1);
}

static void alphabet_policy_tests(void) {
    const char *policy = qdaln_alphabet_policy();

    assert(policy != NULL);
    assert(strcmp(policy, QDALN_ALPHABET_POLICY) == 0);
    assert(strstr(policy, "literal-byte") != NULL);
    assert(strstr(policy, "no wildcard expansion") != NULL);
    assert(qdaln_edit_distance("N", 1, "A", 1) == 1);
    assert(qdaln_edit_distance("R", 1, "A", 1) == 1);
    assert(qdaln_edit_distance("R", 1, "R", 1) == 0);
}

static qdaln_match_result oracle_one(const char *read, size_t read_len,
                                     const char *const *targets, const size_t *target_lens,
                                     size_t n_targets, int k) {
    qdaln_match_result r = {-1, -1, -1, 0, QDALN_MATCH_NONE};
    if ((read == NULL && read_len != 0) || k < 0) {
        r.status = QDALN_MATCH_INVALID;
        return r;
    }

    int tie_best = 0;
    for (size_t i = 0; i < n_targets; ++i) {
        int d = qdaln_edit_distance_dp(read, read_len, targets[i], target_lens[i]);
        assert(d >= 0);
        if (d <= k) {
            ++r.match_count;
            if (r.best_distance < 0 || d < r.best_distance) {
                r.second_best_distance = r.best_distance;
                r.best_distance = d;
                r.target_index = (int)i;
                tie_best = 1;
            } else if (d == r.best_distance) {
                ++tie_best;
            } else if (r.second_best_distance < 0 || d < r.second_best_distance) {
                r.second_best_distance = d;
            }
        }
    }

    if (r.match_count == 0) {
        r.status = QDALN_MATCH_NONE;
    } else if (tie_best > 1) {
        r.status = QDALN_MATCH_AMBIGUOUS;
    } else {
        r.status = QDALN_MATCH_UNIQUE;
    }
    return r;
}

static void assert_match_result(qdaln_match_result got, qdaln_match_result want) {
    assert(got.target_index == want.target_index);
    assert(got.best_distance == want.best_distance);
    assert(got.second_best_distance == want.second_best_distance);
    assert(got.match_count == want.match_count);
    assert(got.status == want.status);
}

static void batch_fixed_tests(void) {
    const char *targets[] = {"ACGT", "AGGT", "TTTT", "ACGA"};
    size_t target_lens[] = {4, 4, 4, 4};
    const char *reads[] = {"ACGT", "ACGC", "CCCC", "", "ACGTT"};
    size_t read_lens[] = {4, 4, 4, 0, 5};
    qdaln_match_result results[5];

    assert(qdaln_match_many(reads, read_lens, 5, targets, target_lens, 4, 0, results) == 0);
    assert(results[0].status == QDALN_MATCH_UNIQUE);
    assert(results[0].target_index == 0);
    assert(results[0].best_distance == 0);
    assert(results[1].status == QDALN_MATCH_NONE);
    assert(results[2].status == QDALN_MATCH_NONE);
    assert(results[3].status == QDALN_MATCH_NONE);
    assert(results[4].status == QDALN_MATCH_NONE);

    assert(qdaln_match_many(reads, read_lens, 5, targets, target_lens, 4, 1, results) == 0);
    assert(results[0].status == QDALN_MATCH_UNIQUE);
    assert(results[0].match_count == 3);
    assert(results[0].second_best_distance == 1);
    assert(results[1].status == QDALN_MATCH_AMBIGUOUS);
    assert(results[1].best_distance == 1);
    assert(results[1].match_count == 2);
    assert(results[4].status == QDALN_MATCH_UNIQUE);
    assert(results[4].target_index == 0);

    assert(qdaln_match_many(reads, read_lens, 0, targets, target_lens, 4, 1, results) == 0);
    assert(qdaln_match_many(reads, read_lens, 5, targets, target_lens, 0, 1, results) == 0);
    for (size_t i = 0; i < 5; ++i) assert(results[i].status == QDALN_MATCH_NONE);

    assert(qdaln_match_many(NULL, read_lens, 5, targets, target_lens, 4, 1, results) == -1);
    assert(qdaln_match_many(reads, read_lens, 5, NULL, target_lens, 4, 1, results) == -1);
    assert(qdaln_match_many(reads, read_lens, 5, targets, target_lens, 4, 1, NULL) == -1);
    assert(qdaln_match_many(reads, read_lens, 5, targets, target_lens, 4, -1, results) == -1);
}

static void assignment_contract_tests(void) {
    const char *targets[] = {"ACGT", "ACGA", "ACGTT"};
    size_t target_lens[] = {4, 4, 5};
    const char *reads[] = {"ACGT", "ACGC", "ACGTT", "ACG", "TTTT"};
    size_t read_lens[] = {4, 4, 5, 3, 4};
    qdaln_assignment_result best[5];
    qdaln_assignment_result radius[5];

    assert(qdaln_assign_many(reads, read_lens, 5, targets, target_lens, 3, 1, QDALN_POLICY_BEST, best) == 0);
    assert(qdaln_assign_many(reads, read_lens, 5, targets, target_lens, 3, 1, QDALN_POLICY_RADIUS, radius) == 0);

    assert(best[0].status == QDALN_MATCH_UNIQUE);
    assert(best[0].target_index == 0);
    assert(best[0].distance == 0);
    assert(best[0].num_best_targets == 1);
    assert(best[0].num_targets_within_radius == 3);
    assert(best[0].edit_class == QDALN_EDIT_EXACT);
    assert(radius[0].status == QDALN_MATCH_AMBIGUOUS);

    assert(best[1].status == QDALN_MATCH_AMBIGUOUS);
    assert(best[1].distance == 1);
    assert(best[1].num_best_targets == 2);
    assert(best[1].num_targets_within_radius == 2);
    assert(best[1].edit_class == QDALN_EDIT_K1_SUB);

    assert(best[2].status == QDALN_MATCH_UNIQUE);
    assert(best[2].target_index == 2);
    assert(best[2].edit_class == QDALN_EDIT_EXACT);
    assert(radius[2].status == QDALN_MATCH_AMBIGUOUS);

    assert(best[3].status == QDALN_MATCH_AMBIGUOUS);
    assert(best[3].edit_class == QDALN_EDIT_K1_DEL);

    assert(best[4].status == QDALN_MATCH_NONE);
    assert(best[4].edit_class == QDALN_EDIT_NONE);

    assert(qdaln_assign_many(reads, read_lens, 5, targets, target_lens, 3, 1, 99, best) == -1);
}

static void index_fixed_tests(void) {
    const char *targets[] = {"ACGT", "AGGT", "ACGA", "ACGTT", "NNNN"};
    size_t target_lens[] = {4, 4, 4, 5, 4};
    const char *reads[] = {"ACGT", "ACGC", "ACGTT", "ACG", "NNNN", "TTTT"};
    size_t read_lens[] = {4, 4, 5, 3, 4, 4};
    qdaln_match_result scan[6];
    qdaln_match_result indexed[6];

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 5);
    assert(idx != NULL);

    for (int k = 0; k <= 3; ++k) {
        assert(qdaln_match_many(reads, read_lens, 6, targets, target_lens, 5, k, scan) == 0);
        assert(qdaln_index_assign(idx, reads, read_lens, 6, k, indexed) == 0);
        for (size_t i = 0; i < 6; ++i) assert_match_result(indexed[i], scan[i]);
    }

    assert(qdaln_index_assign(NULL, reads, read_lens, 6, 1, indexed) == -1);
    assert(qdaln_index_assign(idx, reads, read_lens, 6, -1, indexed) == -1);
    assert(qdaln_index_assign(idx, reads, read_lens, 6, 1, NULL) == -1);
    qdaln_index_free(idx);
    qdaln_index_free(NULL);
    assert(qdaln_index_build(NULL, target_lens, 5) == NULL);
}

static void index_duplicate_exact_tests(void) {
    const char *targets[] = {"ACGT", "ACGT", "AGGT"};
    size_t target_lens[] = {4, 4, 4};
    const char *reads[] = {"ACGT", "AGGT", "TTTT"};
    size_t read_lens[] = {4, 4, 4};
    qdaln_match_result indexed[3];
    qdaln_index_stats stats;

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 3);
    assert(idx != NULL);
    assert(qdaln_index_assign_stats(idx, reads, read_lens, 3, 0, indexed, &stats) == 0);

    assert(indexed[0].status == QDALN_MATCH_AMBIGUOUS);
    assert(indexed[0].target_index == 0);
    assert(indexed[0].best_distance == 0);
    assert(indexed[0].match_count == 2);

    assert(indexed[1].status == QDALN_MATCH_UNIQUE);
    assert(indexed[1].target_index == 2);
    assert(indexed[1].best_distance == 0);
    assert(indexed[1].match_count == 1);

    assert(indexed[2].status == QDALN_MATCH_NONE);
    assert(stats.candidates_verified == 3);
    qdaln_index_free(idx);
}

static void index_stats_pruning_tests(void) {
    const char *targets[] = {"AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT"};
    size_t target_lens[] = {8, 8, 8, 8};
    const char *reads[] = {"AAAAAAAT", "CCCCCCCA", "GGGGGGGA", "TTTTTTTA"};
    size_t read_lens[] = {8, 8, 8, 8};
    qdaln_match_result scan[4];
    qdaln_match_result indexed[4];
    qdaln_index_stats stats;

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 4);
    assert(idx != NULL);
    assert(qdaln_match_many(reads, read_lens, 4, targets, target_lens, 4, 1, scan) == 0);
    assert(qdaln_index_assign_stats(idx, reads, read_lens, 4, 1, indexed, &stats) == 0);
    for (size_t i = 0; i < 4; ++i) assert_match_result(indexed[i], scan[i]);
    assert(stats.candidates_verified < 4 * 4);
    assert(stats.candidates_considered == stats.candidates_verified);
    qdaln_index_free(idx);
}

static void levenshtein_k1_avoids_false_deletion_seed_candidates_tests(void) {
    const char *targets[] = {
        "ACGT",   /* exact */
        "ACGA",   /* one substitution */
        "TACGT",  /* target has one inserted base */
        "CGT",    /* read has one inserted base */
        "AGCT",   /* shares a deletion seed with ACGT but edit distance is 2 */
        "ATGC",   /* shares a deletion seed with ACGT but edit distance is 2 */
    };
    size_t target_lens[] = {4, 4, 5, 3, 4, 4};
    const char *reads[] = {"ACGT"};
    size_t read_lens[] = {4};
    qdaln_match_result scan[1];
    qdaln_match_result indexed[1];
    qdaln_index_stats stats;

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 6);
    assert(idx != NULL);
    assert(qdaln_match_many(reads, read_lens, 1, targets, target_lens, 6, 1, scan) == 0);
    assert(qdaln_index_assign_stats(idx, reads, read_lens, 1, 1, indexed, &stats) == 0);
    assert_match_result(indexed[0], scan[0]);
    assert(indexed[0].match_count == 4);
    assert(stats.candidates_verified == 4);
    assert(stats.candidates_considered == stats.candidates_verified);
    qdaln_index_free(idx);
}

static void hamming_single_unknown_uses_index_tests(void) {
    const char *targets[] = {"ACGT", "ACCT", "TTTT", "CCCC"};
    size_t target_lens[] = {4, 4, 4, 4};
    const char *reads[] = {"ACNT", "NNNN"};
    size_t read_lens[] = {4, 4};
    qdaln_match_result scan[2];
    qdaln_match_result indexed[2];
    qdaln_index_stats stats;

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 4);
    assert(idx != NULL);
    assert(qdaln_match_many(reads, read_lens, 2, targets, target_lens, 4, 1, scan) == 0);
    assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 2, 1, indexed, &stats) == 0);
    for (size_t i = 0; i < 2; ++i) assert_match_result(indexed[i], scan[i]);

    assert(indexed[0].status == QDALN_MATCH_AMBIGUOUS);
    assert(indexed[0].best_distance == 1);
    assert(indexed[0].match_count == 2);
    assert(indexed[1].status == QDALN_MATCH_NONE);
    assert(stats.candidates_verified < 4 * 2);
    assert(stats.candidates_considered == stats.candidates_verified);
    qdaln_index_free(idx);
}

static void levenshtein_non_acgt_indel_uses_index_tests(void) {
    const char *targets[] = {"ACGT", "TGCA"};
    size_t target_lens[] = {4, 4};
    const char *reads[] = {"ACNGT", "ANNGT"};
    size_t read_lens[] = {5, 5};
    qdaln_match_result scan[2];
    qdaln_match_result indexed[2];
    qdaln_match_result status_only[2];
    qdaln_index_stats stats;
    qdaln_index_stats status_stats;

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 2);
    assert(idx != NULL);
    assert(qdaln_match_many(reads, read_lens, 2, targets, target_lens, 2, 1, scan) == 0);
    assert(qdaln_index_assign_stats(idx, reads, read_lens, 2, 1, indexed, &stats) == 0);
    assert(qdaln_index_assign_status_stats(idx, reads, read_lens, 2, 1, status_only, &status_stats) == 0);

    for (size_t i = 0; i < 2; ++i) {
        assert_match_result(indexed[i], scan[i]);
        assert(status_only[i].status == scan[i].status);
        assert(status_only[i].best_distance == scan[i].best_distance);
        if (scan[i].status == QDALN_MATCH_UNIQUE) assert(status_only[i].target_index == scan[i].target_index);
    }
    assert(indexed[0].status == QDALN_MATCH_UNIQUE);
    assert(indexed[0].target_index == 0);
    assert(indexed[0].best_distance == 1);
    assert(indexed[1].status == QDALN_MATCH_NONE);
    assert(status_stats.candidates_verified < 2 * 2);
    qdaln_index_free(idx);
}

static void index_status_shortcut_stops_after_ambiguity_tests(void) {
    char targets_buf[20][21];
    const char *targets[20];
    size_t target_lens[20];
    const char *reads[] = {"AAAAAAAAAAAAAAAAAAA"};
    size_t read_lens[] = {19};
    qdaln_match_result exhaustive[1];
    qdaln_match_result status_only[1];
    qdaln_index_stats exhaustive_stats;
    qdaln_index_stats status_stats;

    for (size_t i = 0; i < 20; ++i) {
        memset(targets_buf[i], 'A', 20);
        targets_buf[i][i] = 'C';
        targets_buf[i][20] = '\0';
        targets[i] = targets_buf[i];
        target_lens[i] = 20;
    }

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 20);
    assert(idx != NULL);
    assert(qdaln_index_assign_stats(idx, reads, read_lens, 1, 1, exhaustive, &exhaustive_stats) == 0);
    assert(qdaln_index_assign_status_stats(idx, reads, read_lens, 1, 1, status_only, &status_stats) == 0);

    assert(exhaustive[0].status == QDALN_MATCH_AMBIGUOUS);
    assert(exhaustive[0].best_distance == 1);
    assert(exhaustive[0].match_count == 20);
    assert(status_only[0].status == QDALN_MATCH_AMBIGUOUS);
    assert(status_only[0].best_distance == exhaustive[0].best_distance);
    assert(status_stats.candidates_verified < exhaustive_stats.candidates_verified);
    assert(status_stats.candidates_verified <= 2);
    qdaln_index_free(idx);
}

typedef struct large_panel {
    char (*storage)[34];
    const char **targets;
    size_t *lens;
    size_t count;
} large_panel;

static void free_large_panel(large_panel *panel) {
    free(panel->storage);
    free(panel->targets);
    free(panel->lens);
    panel->storage = NULL;
    panel->targets = NULL;
    panel->lens = NULL;
    panel->count = 0;
}

static void set_panel_target(large_panel *panel, size_t i, const char *seq) {
    size_t len = strlen(seq);
    assert(len < sizeof(panel->storage[i]));
    memcpy(panel->storage[i], seq, len + 1);
    panel->targets[i] = panel->storage[i];
    panel->lens[i] = len;
}

static void encode_background_target(char *out, uint32_t value) {
    for (size_t bit = 0; bit < 16; ++bit) {
        if ((value >> bit) & 1U) {
            out[2 * bit] = 'T';
            out[2 * bit + 1] = 'G';
        } else {
            out[2 * bit] = 'G';
            out[2 * bit + 1] = 'T';
        }
    }
    out[32] = '\0';
}

static void build_large_panel(large_panel *panel, size_t count) {
    static const char *special[] = {
        "ACACACACACACACACACAC",
        "CCCCAAAACCCCAAAACCCC",
        "AAAACCCCGGGGTTTTAAAA",
        "TTTTGGGGCCCCAAAATTTT",
        "AACCAACCAACCAACCAACC",
        "AACCAACCAACCAACCAACA",
        "CCCAAACCCAAACCCAAANC",
        "ACGTACGT",
        "AACCAACCAACCAACCAACCAACCAACCAACC",
    };
    const size_t n_special = sizeof(special) / sizeof(special[0]);
    assert(count > n_special);
    panel->storage = (char (*)[34])calloc(count, sizeof(*panel->storage));
    panel->targets = (const char **)calloc(count, sizeof(*panel->targets));
    panel->lens = (size_t *)calloc(count, sizeof(*panel->lens));
    assert(panel->storage != NULL && panel->targets != NULL && panel->lens != NULL);
    panel->count = count;

    for (size_t i = 0; i < n_special; ++i) set_panel_target(panel, i, special[i]);
    for (size_t i = n_special; i < count; ++i) {
        encode_background_target(panel->storage[i], (uint32_t)(i - n_special));
        panel->targets[i] = panel->storage[i];
        panel->lens[i] = 32;
    }
}

static void assert_large_panel_case(large_panel *panel, const char *const *reads,
                                    const size_t *read_lens, size_t n_reads, int k) {
    qdaln_match_result *scan = (qdaln_match_result *)calloc(n_reads, sizeof(qdaln_match_result));
    qdaln_match_result *indexed = (qdaln_match_result *)calloc(n_reads, sizeof(qdaln_match_result));
    assert(scan != NULL && indexed != NULL);
    qdaln_index_stats stats;

    qdaln_index *idx = qdaln_index_build(panel->targets, panel->lens, panel->count);
    assert(idx != NULL);
    assert(qdaln_match_many(reads, read_lens, n_reads, panel->targets, panel->lens, panel->count, k, scan) == 0);
    assert(qdaln_index_assign_stats(idx, reads, read_lens, n_reads, k, indexed, &stats) == 0);
    for (size_t i = 0; i < n_reads; ++i) assert_match_result(indexed[i], scan[i]);
    if (k == 1) assert(stats.candidates_verified < panel->count);
    qdaln_index_free(idx);
    free(scan);
    free(indexed);
}

static void large_panel_oracle_tests(void) {
    const size_t panel_sizes[] = {1024, 16384, 65536};
    char substitution[34];
    char insertion[34];
    char deletion[34];
    char ambiguous[34];
    char edge32_deletion[34];

    strcpy(substitution, "CCCCAAAACCCCAAAACCCA");
    strcpy(insertion, "AAAACCCCGGGGGTTTTAAAA");
    strcpy(deletion, "TTTGGGGCCCCAAAATTTT");
    strcpy(ambiguous, "AACCAACCAACCAACCAACG");
    strcpy(edge32_deletion, "AACCAACCAACCAACCAACCAACCAACCAAC");

    const char *reads[] = {
        "ACACACACACACACACACAC",
        substitution,
        insertion,
        deletion,
        ambiguous,
        "CCCAAACCCAAACCCAAANC",
        "ACGTACGT",
        edge32_deletion,
        "GGGGGGGG",
        NULL,
    };
    size_t read_lens[] = {
        20,
        strlen(substitution),
        strlen(insertion),
        strlen(deletion),
        strlen(ambiguous),
        20,
        8,
        strlen(edge32_deletion),
        8,
        4,
    };
    const size_t n_reads = sizeof(reads) / sizeof(reads[0]);

    for (size_t p = 0; p < sizeof(panel_sizes) / sizeof(panel_sizes[0]); ++p) {
        large_panel panel = {0};
        build_large_panel(&panel, panel_sizes[p]);
        assert_large_panel_case(&panel, reads, read_lens, n_reads, 0);
        assert_large_panel_case(&panel, reads, read_lens, n_reads, 1);

        qdaln_match_result scan[10];
        assert(qdaln_match_many(reads, read_lens, n_reads, panel.targets, panel.lens, panel.count, 1, scan) == 0);
        assert(scan[0].status == QDALN_MATCH_UNIQUE && scan[0].best_distance == 0);
        assert(scan[1].status == QDALN_MATCH_UNIQUE && scan[1].best_distance == 1);
        assert(scan[2].status == QDALN_MATCH_UNIQUE && scan[2].best_distance == 1);
        assert(scan[3].status == QDALN_MATCH_UNIQUE && scan[3].best_distance == 1);
        assert(scan[4].status == QDALN_MATCH_AMBIGUOUS && scan[4].best_distance == 1);
        assert(scan[5].status == QDALN_MATCH_UNIQUE && scan[5].best_distance == 0);
        assert(scan[6].status == QDALN_MATCH_UNIQUE && scan[6].best_distance == 0);
        assert(scan[7].status == QDALN_MATCH_UNIQUE && scan[7].best_distance == 1);
        assert(scan[8].status == QDALN_MATCH_NONE);
        assert(scan[9].status == QDALN_MATCH_INVALID);
        free_large_panel(&panel);
    }
}

static void fuzz_tests(void) {
    char a[129];
    char b[129];

    for (size_t t = 0; t < 50000; ++t) {
        size_t a_len = (size_t)(xorshift64() % 65ULL);
        size_t b_len = (size_t)(xorshift64() % 129ULL);
        rand_seq(a, a_len);
        rand_seq(b, b_len);
        check_pair(a, b);
    }
}

static void batch_fuzz_tests(void) {
    char reads_buf[8][33];
    char targets_buf[8][33];
    const char *reads[8];
    const char *targets[8];
    size_t read_lens[8];
    size_t target_lens[8];
    qdaln_match_result got[8];

    for (size_t t = 0; t < 5000; ++t) {
        size_t n_reads = 1 + (size_t)(xorshift64() % 8ULL);
        size_t n_targets = 1 + (size_t)(xorshift64() % 8ULL);
        int k = (int)(xorshift64() % 4ULL);
        for (size_t i = 0; i < n_reads; ++i) {
            read_lens[i] = (size_t)(xorshift64() % 33ULL);
            rand_seq(reads_buf[i], read_lens[i]);
            reads[i] = reads_buf[i];
        }
        for (size_t i = 0; i < n_targets; ++i) {
            target_lens[i] = (size_t)(xorshift64() % 33ULL);
            rand_seq(targets_buf[i], target_lens[i]);
            targets[i] = targets_buf[i];
        }

        assert(qdaln_match_many(reads, read_lens, n_reads, targets, target_lens, n_targets, k, got) == 0);
        for (size_t i = 0; i < n_reads; ++i) {
            qdaln_match_result want = oracle_one(reads[i], read_lens[i], targets, target_lens, n_targets, k);
            assert_match_result(got[i], want);
        }
    }
}

static void index_fuzz_tests(void) {
    char reads_buf[8][33];
    char targets_buf[8][33];
    const char *reads[8];
    const char *targets[8];
    size_t read_lens[8];
    size_t target_lens[8];
    qdaln_match_result scan[8];
    qdaln_match_result indexed[8];

    for (size_t t = 0; t < 5000; ++t) {
        size_t n_reads = 1 + (size_t)(xorshift64() % 8ULL);
        size_t n_targets = 1 + (size_t)(xorshift64() % 8ULL);
        int k = (int)(xorshift64() % 4ULL);
        for (size_t i = 0; i < n_reads; ++i) {
            read_lens[i] = (size_t)(xorshift64() % 33ULL);
            rand_seq(reads_buf[i], read_lens[i]);
            reads[i] = reads_buf[i];
        }
        for (size_t i = 0; i < n_targets; ++i) {
            target_lens[i] = (size_t)(xorshift64() % 33ULL);
            rand_seq(targets_buf[i], target_lens[i]);
            targets[i] = targets_buf[i];
        }

        qdaln_index *idx = qdaln_index_build(targets, target_lens, n_targets);
        assert(idx != NULL);
        assert(qdaln_match_many(reads, read_lens, n_reads, targets, target_lens, n_targets, k, scan) == 0);
        assert(qdaln_index_assign(idx, reads, read_lens, n_reads, k, indexed) == 0);
        for (size_t i = 0; i < n_reads; ++i) assert_match_result(indexed[i], scan[i]);
        qdaln_index_free(idx);
    }
}

int main(void) {
    fixed_tests();
    alphabet_policy_tests();
    batch_fixed_tests();
    assignment_contract_tests();
    index_fixed_tests();
    index_duplicate_exact_tests();
    index_stats_pruning_tests();
    levenshtein_k1_avoids_false_deletion_seed_candidates_tests();
    hamming_single_unknown_uses_index_tests();
    levenshtein_non_acgt_indel_uses_index_tests();
    index_status_shortcut_stops_after_ambiguity_tests();
    large_panel_oracle_tests();
    fuzz_tests();
    batch_fuzz_tests();
    index_fuzz_tests();
    puts("qdalign tests passed");
    return 0;
}
