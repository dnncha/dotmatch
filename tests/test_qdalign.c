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
    /* Now myers64 generalized to >64 via multi-word; always compare for oracle coverage */
    assert(dp == fast_direct);

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
    check_pair("ACG", "CGA");
    check_pair("ABCDEFGHIJKLMNOPQ", "QBCDEFGHIJKLMNOPA");

    /* Explicit longer patterns (>64) for multi-word Myers coverage + oracle */
    {
        char longp[80];
        char longt[80];
        memset(longp, 'A', 79); longp[79] = '\0';
        memset(longt, 'A', 78); longt[78] = 'T'; longt[79] = '\0';
        check_pair(longp, longt); /* dist 2: one sub + one del? wait len diff, actually 1sub +1del? but dp will tell */
        check_pair(longp, longp);
        /* 70bp with 1 edit */
        char lp[71]; char lt[71];
        memset(lp, 'C', 70); lp[70]='\0';
        memcpy(lt, lp, 71);
        lt[35] = 'G';
        check_pair(lp, lt);
    }

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

static int hamming_oracle_distance(const char *a, size_t a_len, const char *b, size_t b_len, int k) {
    if (a_len != b_len) return -1;
    int d = 0;
    for (size_t i = 0; i < a_len; ++i) {
        if (a[i] != b[i] && ++d > k) return -1;
    }
    return d;
}

static qdaln_match_result hamming_oracle_one(const char *read, size_t read_len,
                                             const char *const *targets, const size_t *target_lens,
                                             size_t n_targets, int k) {
    qdaln_match_result r = {-1, -1, -1, 0, QDALN_MATCH_NONE};
    if ((read == NULL && read_len != 0) || k < 0) {
        r.status = QDALN_MATCH_INVALID;
        return r;
    }

    int tie_best = 0;
    for (size_t i = 0; i < n_targets; ++i) {
        int d = hamming_oracle_distance(read, read_len, targets[i], target_lens[i], k);
        if (d < 0) continue;
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

static void empty_index_tests(void) {
    const char *reads[] = {"ACGT", ""};
    size_t read_lens[] = {4, 0};
    qdaln_match_result results[2];
    qdaln_index_stats stats;

    qdaln_index *idx = qdaln_index_build(NULL, NULL, 0);
    assert(idx != NULL);

    for (int k = 0; k <= 2; ++k) {
        assert(qdaln_index_assign_stats(idx, reads, read_lens, 2, k, results, &stats) == 0);
        assert(results[0].status == QDALN_MATCH_NONE);
        assert(results[1].status == QDALN_MATCH_NONE);
        assert(stats.candidates_considered == 0);
        assert(stats.candidates_verified == 0);
    }

    assert(qdaln_index_assign_status_stats(idx, reads, read_lens, 2, 1, results, &stats) == 0);
    assert(results[0].status == QDALN_MATCH_NONE);
    assert(results[1].status == QDALN_MATCH_NONE);
    assert(stats.candidates_considered == 0);
    assert(stats.candidates_verified == 0);

    assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 2, 1, results, &stats) == 0);
    assert(results[0].status == QDALN_MATCH_NONE);
    assert(results[1].status == QDALN_MATCH_NONE);
    assert(stats.candidates_considered == 0);
    assert(stats.candidates_verified == 0);
    qdaln_index_free(idx);
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

static void direct_exact_lookup_tests(void) {
    const char *targets[] = {"ACGT", "ACGT", "AGGT", "NNNN"};
    size_t target_lens[] = {4, 4, 4, 4};
    const char *reads[] = {"ACGT", "AGGT", "TTTT", "NNNN", "acgt", "nnnn"};
    size_t read_lens[] = {4, 4, 4, 4};
    qdaln_match_result direct;
    qdaln_match_result batch[1];
    qdaln_index_stats direct_stats;
    qdaln_index_stats batch_stats;

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 4);
    assert(idx != NULL);
    for (size_t i = 0; i < 4; ++i) {
        assert(qdaln_index_lookup_exact_stats(idx, reads[i], read_lens[i], &direct, &direct_stats) == 0);
        assert(qdaln_index_assign_hamming_stats(idx, &reads[i], &read_lens[i], 1, 0, batch, &batch_stats) == 0);
        assert_match_result(direct, batch[0]);
        assert(direct_stats.candidates_considered == batch_stats.candidates_considered);
        assert(direct_stats.candidates_verified == batch_stats.candidates_verified);
    }
    qdaln_match_result direct_many[4];
    qdaln_match_result batch_many[4];
    assert(qdaln_index_lookup_exact_many_stats(idx, reads, read_lens, 4, direct_many, &direct_stats) == 0);
    assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 4, 0, batch_many, &batch_stats) == 0);
    for (size_t i = 0; i < 4; ++i) assert_match_result(direct_many[i], batch_many[i]);
    assert(direct_stats.candidates_considered == batch_stats.candidates_considered);
    assert(direct_stats.candidates_verified == batch_stats.candidates_verified);
    assert(qdaln_index_lookup_exact_ascii_many_stats(idx, reads, read_lens, 4, direct_many, &direct_stats) == 0);
    assert(qdaln_index_lookup_exact_ascii_many_stats(idx, reads, read_lens, 4, batch_many, &batch_stats) == 0);
    for (size_t i = 0; i < 4; ++i) assert_match_result(direct_many[i], batch_many[i]);
    assert(qdaln_index_lookup_exact_ascii_stats(idx, reads[4], 4, &direct, &direct_stats) == 0);
    assert(direct.status == QDALN_MATCH_AMBIGUOUS);
    assert(direct.target_index == 0);
    assert(direct.match_count == 2);
    assert(qdaln_index_lookup_exact_ascii_stats(idx, reads[5], 4, &direct, &direct_stats) == 0);
    assert(direct.status == QDALN_MATCH_UNIQUE);
    assert(direct.target_index == 3);

    assert(qdaln_index_lookup_exact_stats(NULL, reads[0], read_lens[0], &direct, &direct_stats) == -1);
    assert(qdaln_index_lookup_exact_stats(idx, reads[0], read_lens[0], NULL, &direct_stats) == -1);
    assert(qdaln_index_lookup_exact_many_stats(NULL, reads, read_lens, 4, direct_many, &direct_stats) == -1);
    assert(qdaln_index_lookup_exact_many_stats(idx, NULL, read_lens, 4, direct_many, &direct_stats) == -1);
    assert(qdaln_index_lookup_exact_many_stats(idx, reads, NULL, 4, direct_many, &direct_stats) == -1);
    assert(qdaln_index_lookup_exact_many_stats(idx, reads, read_lens, 4, NULL, &direct_stats) == -1);
    assert(qdaln_index_lookup_exact_ascii_many_stats(NULL, reads, read_lens, 4, direct_many, &direct_stats) == -1);
    assert(qdaln_index_lookup_exact_ascii_stats(NULL, reads[0], read_lens[0], &direct, &direct_stats) == -1);
    assert(qdaln_index_lookup_exact_ascii_stats(idx, reads[0], read_lens[0], NULL, &direct_stats) == -1);
    assert(qdaln_index_lookup_exact_stats(idx, NULL, 1, &direct, &direct_stats) == 0);
    assert(direct.status == QDALN_MATCH_INVALID);
    assert(direct_stats.candidates_considered == 0);
    assert(direct_stats.candidates_verified == 0);
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

static void levenshtein_k2_uses_index_without_full_scan_tests(void) {
    const char *targets[] = {
        "ACGTACGT",    /* exact */
        "ACGTTCGA",    /* two substitutions */
        "ACGTTACGT",   /* target has one inserted base */
        "ACGTTTACGT",  /* target has two inserted bases */
        "ACGACGT",     /* read has one inserted base */
        "ACACGT",      /* read has two inserted bases */
        "TTTTTTTT",
        "CCCCCCCC",
        "GGGGGGGG",
        "TGTGTGTG",
        "CACACACA",
        "GATCGATC",
        "TTAACCGG",
        "CCGGTTAA",
        "GGCCAATT",
        "TCCGTCCG",
    };
    size_t target_lens[] = {8, 8, 9, 10, 7, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
    const char *reads[] = {"ACGTACGT"};
    size_t read_lens[] = {8};
    qdaln_match_result scan[1];
    qdaln_match_result indexed[1];
    qdaln_index_stats stats;

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 16);
    assert(idx != NULL);
    assert(qdaln_match_many(reads, read_lens, 1, targets, target_lens, 16, 2, scan) == 0);
    assert(qdaln_index_assign_stats(idx, reads, read_lens, 1, 2, indexed, &stats) == 0);

    assert(scan[0].match_count == 6);
    assert_match_result(indexed[0], scan[0]);
    assert(stats.candidates_verified < 16);
    assert(stats.candidates_considered == stats.candidates_verified);
    qdaln_index_free(idx);
}

static void levenshtein_k2_len32_insertion_uses_index_tests(void) {
    char target0[33];
    char target1[33];
    char target2[33];
    char target3[33];
    char read0[34];
    memset(target0, 'A', 32);
    target0[32] = '\0';
    strcpy(target1, target0);
    target1[0] = 'C';
    strcpy(target2, target0);
    target2[16] = 'G';
    strcpy(target3, target0);
    target3[31] = 'T';
    memcpy(read0, target0, 16);
    read0[16] = 'C';
    memcpy(read0 + 17, target0 + 16, 16);
    read0[33] = '\0';

    const char *targets[] = {
        target0,
        target1,
        target2,
        target3,
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",
    };
    size_t target_lens[] = {32, 32, 32, 32, 32, 32};
    const char *reads[] = {read0};
    size_t read_lens[] = {33};
    qdaln_match_result scan[1];
    qdaln_match_result indexed[1];
    qdaln_index_stats stats;

    qdaln_index *idx = qdaln_index_build(targets, target_lens, sizeof(targets) / sizeof(targets[0]));
    assert(idx != NULL);
    assert(qdaln_match_many(reads, read_lens, 1, targets, target_lens, sizeof(targets) / sizeof(targets[0]), 2, scan) == 0);
    assert(qdaln_index_assign_stats(idx, reads, read_lens, 1, 2, indexed, &stats) == 0);
    assert_match_result(indexed[0], scan[0]);
    assert(indexed[0].status == QDALN_MATCH_UNIQUE);
    assert(indexed[0].target_index == 0);
    assert(stats.candidates_verified < sizeof(targets) / sizeof(targets[0]));
    assert(stats.candidates_considered == stats.candidates_verified);
    qdaln_index_free(idx);
}

static void levenshtein_k2_single_unknown_uses_bounded_candidates_tests(void) {
    enum { N_TARGETS = 1024, TARGET_LEN = 12, MAX_TARGET_LEN = 13 };
    char target_buf[N_TARGETS][MAX_TARGET_LEN + 1];
    const char *targets[N_TARGETS];
    size_t target_lens[N_TARGETS];
    const char *reads[] = {"ACNTACGTACGT", "ACNGTACGTACGT", "NNNNNNNNNNNN"};
    size_t read_lens[] = {TARGET_LEN, TARGET_LEN, TARGET_LEN};
    qdaln_match_result scan[3];
    qdaln_match_result indexed[3];
    qdaln_index_stats stats;

    strcpy(target_buf[0], "ACGTACGTACGT");
    strcpy(target_buf[1], "ACTACGTACGT");
    strcpy(target_buf[2], "ACGTTACGTACGT");
    strcpy(target_buf[3], "ACCTACGTACGT");
    for (size_t i = 4; i < N_TARGETS; ++i) {
        size_t v = i * 11400714819323198485ULL;
        for (size_t pos = 0; pos < TARGET_LEN; ++pos) {
            target_buf[i][pos] = "TGCA"[v & 3U];
            v >>= 2;
        }
        target_buf[i][TARGET_LEN] = '\0';
    }
    for (size_t i = 0; i < N_TARGETS; ++i) {
        targets[i] = target_buf[i];
        target_lens[i] = strlen(target_buf[i]);
    }

    qdaln_index *idx = qdaln_index_build(targets, target_lens, N_TARGETS);
    assert(idx != NULL);
    assert(qdaln_match_many(reads, read_lens, 3, targets, target_lens, N_TARGETS, 2, scan) == 0);
    assert(qdaln_index_assign_stats(idx, reads, read_lens, 3, 2, indexed, &stats) == 0);
    for (size_t i = 0; i < 3; ++i) assert_match_result(indexed[i], scan[i]);
    assert(indexed[0].status != QDALN_MATCH_NONE);
    assert(indexed[1].status != QDALN_MATCH_NONE);
    assert(indexed[2].status == QDALN_MATCH_NONE);
    assert(stats.candidates_considered == stats.candidates_verified);
    assert(stats.candidates_verified < N_TARGETS * 3);

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

static void hamming_single_unknown_k2_k3_uses_bounded_candidates_tests(void) {
    enum { N_TARGETS = 512, TARGET_LEN = 12 };
    char target_buf[N_TARGETS][TARGET_LEN + 1];
    const char *targets[N_TARGETS];
    size_t target_lens[N_TARGETS];
    const char *reads[] = {"ACNTACGTACGT", "ACNTACGTACGA", "NNNNNNNNNNNN"};
    size_t read_lens[] = {TARGET_LEN, TARGET_LEN, TARGET_LEN};
    qdaln_match_result indexed[3];
    qdaln_index_stats stats;

    strcpy(target_buf[0], "ACGTACGTACGT");
    strcpy(target_buf[1], "ACGTACGTACGA");
    strcpy(target_buf[2], "ACCTACGTACGA");
    strcpy(target_buf[3], "TCCTACGTACGA");
    for (size_t i = 4; i < N_TARGETS; ++i) {
        size_t v = i * 2654435761U;
        for (size_t pos = 0; pos < TARGET_LEN; ++pos) {
            target_buf[i][pos] = "ACGT"[v & 3U];
            v >>= 2;
        }
        target_buf[i][TARGET_LEN] = '\0';
    }
    for (size_t i = 0; i < N_TARGETS; ++i) {
        targets[i] = target_buf[i];
        target_lens[i] = TARGET_LEN;
    }

    qdaln_index *idx = qdaln_index_build(targets, target_lens, N_TARGETS);
    assert(idx != NULL);

    assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 3, 2, indexed, &stats) == 0);
    for (size_t i = 0; i < 3; ++i) {
        qdaln_match_result want = hamming_oracle_one(reads[i], read_lens[i], targets, target_lens, N_TARGETS, 2);
        assert_match_result(indexed[i], want);
    }
    assert(stats.candidates_considered == stats.candidates_verified);
    assert(stats.candidates_verified < N_TARGETS * 3);

    assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 3, 3, indexed, &stats) == 0);
    for (size_t i = 0; i < 3; ++i) {
        qdaln_match_result want = hamming_oracle_one(reads[i], read_lens[i], targets, target_lens, N_TARGETS, 3);
        assert_match_result(indexed[i], want);
    }
    assert(stats.candidates_considered == stats.candidates_verified);
    assert(stats.candidates_verified < N_TARGETS * 3);

    qdaln_index_free(idx);
}

static void hamming_multi_unknown_k2_k3_uses_bounded_candidates_tests(void) {
    enum { N_TARGETS = 1024, TARGET_LEN = 12 };
    char target_buf[N_TARGETS][TARGET_LEN + 1];
    const char *targets[N_TARGETS];
    size_t target_lens[N_TARGETS];
    const char *reads[] = {"ACNTANGTACGT", "NCNTANGTACGT", "NNNNNNNNNNNN"};
    size_t read_lens[] = {TARGET_LEN, TARGET_LEN, TARGET_LEN};
    qdaln_match_result indexed[3];
    qdaln_index_stats stats;

    strcpy(target_buf[0], "ACGTACGTACGT");
    strcpy(target_buf[1], "ACCTACGTACGT");
    strcpy(target_buf[2], "TCCTACGTACGT");
    strcpy(target_buf[3], "TCCTAGGTACGT");
    for (size_t i = 4; i < N_TARGETS; ++i) {
        size_t v = i * 780291637U;
        for (size_t pos = 0; pos < TARGET_LEN; ++pos) {
            target_buf[i][pos] = "CATG"[v & 3U];
            v >>= 2;
        }
        target_buf[i][TARGET_LEN] = '\0';
    }
    for (size_t i = 0; i < N_TARGETS; ++i) {
        targets[i] = target_buf[i];
        target_lens[i] = TARGET_LEN;
    }

    qdaln_index *idx = qdaln_index_build(targets, target_lens, N_TARGETS);
    assert(idx != NULL);

    assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 3, 2, indexed, &stats) == 0);
    for (size_t i = 0; i < 3; ++i) {
        qdaln_match_result want = hamming_oracle_one(reads[i], read_lens[i], targets, target_lens, N_TARGETS, 2);
        assert_match_result(indexed[i], want);
    }
    assert(indexed[0].status != QDALN_MATCH_NONE);
    assert(indexed[1].status == QDALN_MATCH_NONE);
    assert(indexed[2].status == QDALN_MATCH_NONE);
    assert(stats.candidates_considered == stats.candidates_verified);
    assert(stats.candidates_verified < N_TARGETS * 3);

    assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 3, 3, indexed, &stats) == 0);
    for (size_t i = 0; i < 3; ++i) {
        qdaln_match_result want = hamming_oracle_one(reads[i], read_lens[i], targets, target_lens, N_TARGETS, 3);
        assert_match_result(indexed[i], want);
    }
    assert(indexed[1].status != QDALN_MATCH_NONE);
    assert(indexed[2].status == QDALN_MATCH_NONE);
    assert(stats.candidates_considered == stats.candidates_verified);
    assert(stats.candidates_verified < N_TARGETS * 3);

    qdaln_index_free(idx);
}

static void hamming_seed_index_semantics_tests(void) {
    const char *targets[] = {
        "ACGTACGT",
        "ACGTACGT",
        "ACGTACGA",
        "TTTTTTTT",
        "ACGTTCGT",
        "CCCCCCCN",
    };
    size_t target_lens[] = {8, 8, 8, 8, 8, 8};
    const char *reads[] = {"ACGTACGT", "ACGTACGC", "TTTTTTTA", "CCCCCCCC"};
    size_t read_lens[] = {8, 8, 8, 8};
    qdaln_match_result scan[4];
    qdaln_match_result indexed[4];
    qdaln_index_stats stats;

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 6);
    assert(idx != NULL);
    assert(qdaln_match_many(reads, read_lens, 4, targets, target_lens, 6, 1, scan) == 0);
    assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 4, 1, indexed, &stats) == 0);
    for (size_t i = 0; i < 4; ++i) assert_match_result(indexed[i], scan[i]);
    assert(indexed[0].status == QDALN_MATCH_AMBIGUOUS);
    assert(indexed[0].best_distance == 0);
    assert(indexed[0].second_best_distance == 1);
    assert(indexed[0].match_count == 4);
    assert(indexed[1].status == QDALN_MATCH_AMBIGUOUS);
    assert(indexed[1].best_distance == 1);
    assert(indexed[2].status == QDALN_MATCH_UNIQUE);
    assert(indexed[3].status == QDALN_MATCH_UNIQUE);
    assert(indexed[3].best_distance == 1);
    assert(indexed[3].target_index == 5);
    assert(stats.candidates_considered == stats.candidates_verified);
    qdaln_index_free(idx);
}

static void hamming_k2_k3_seed_index_semantics_tests(void) {
    char len32_read[33];
    char len32_exact[33];
    char len32_d2[33];
    char len32_d3[33];
    char len32_far[33];
    memset(len32_read, 'A', 32);
    len32_read[32] = '\0';
    strcpy(len32_exact, len32_read);
    strcpy(len32_d2, len32_read);
    len32_d2[0] = 'C';
    len32_d2[31] = 'G';
    strcpy(len32_d3, len32_read);
    len32_d3[0] = 'C';
    len32_d3[15] = 'G';
    len32_d3[31] = 'T';
    strcpy(len32_far, len32_read);
    len32_far[0] = 'C';
    len32_far[8] = 'G';
    len32_far[16] = 'T';
    len32_far[24] = 'C';

    const char *targets[] = {
        "ACGT",                 /* length 4 exact */
        "ACGA",                 /* length 4 distance 1 */
        "TCGA",                 /* length 4 distance 2 */
        "TTGA",                 /* length 4 distance 3 */
        "TGCA",                 /* length 4 distance 4 */
        "ACGTACGTACGTACGTACGT", /* length 20 exact */
        "ACGTACGTACGTACGTACGT",
        "TCGTACGTACGTACGTACGA", /* length 20 distance 2 */
        "TCGTACGTACGTACGTTCGA", /* length 20 distance 3 */
        "ACGTACGTACGTACGTACGN", /* non-ACGT target */
        len32_exact,
        len32_d2,
        len32_d3,
        len32_far,
    };
    size_t target_lens[] = {
        4, 4, 4, 4, 4,
        20, 20, 20, 20, 20,
        32, 32, 32, 32,
    };
    const char *reads[] = {
        "ACGT",                 /* exact plus duplicate/near ambiguity accounting */
        "TCGT",                 /* k=2 ambiguous at distance 1 */
        "NNNN",                 /* no match at k=2/k=3 */
        "ACGTACGTACGTACGTACGT", /* length 20 duplicate exact */
        "ACGTACGTACGTACGTACGA", /* length 20 k=3 candidates */
        "ACGTACGTACGTACGTACNN", /* non-ACGT read fallback */
        len32_read,
        "NNNN",                 /* no match at k=3 */
    };
    size_t read_lens[] = {4, 4, 4, 20, 20, 20, 32, 4};
    qdaln_match_result indexed[8];
    qdaln_index_stats stats;

    qdaln_index *idx = qdaln_index_build(targets, target_lens, sizeof(targets) / sizeof(targets[0]));
    assert(idx != NULL);

    assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 8, 2, indexed, &stats) == 0);
    for (size_t i = 0; i < 8; ++i) {
        qdaln_match_result want = hamming_oracle_one(reads[i], read_lens[i], targets, target_lens,
                                                     sizeof(targets) / sizeof(targets[0]), 2);
        assert_match_result(indexed[i], want);
    }
    assert(indexed[0].status == QDALN_MATCH_UNIQUE);
    assert(indexed[0].target_index == 0);
    assert(indexed[0].match_count == 3);
    assert(indexed[1].status == QDALN_MATCH_AMBIGUOUS);
    assert(indexed[1].best_distance == 1);
    assert(indexed[2].status == QDALN_MATCH_NONE);
    assert(indexed[3].status == QDALN_MATCH_AMBIGUOUS);
    assert(indexed[3].best_distance == 0);
    assert(indexed[5].status == QDALN_MATCH_UNIQUE);
    assert(indexed[5].target_index == 9);
    assert(indexed[6].status == QDALN_MATCH_UNIQUE);
    assert(indexed[6].match_count == 2);
    assert(stats.candidates_considered == stats.candidates_verified);

    assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 8, 3, indexed, &stats) == 0);
    for (size_t i = 0; i < 8; ++i) {
        qdaln_match_result want = hamming_oracle_one(reads[i], read_lens[i], targets, target_lens,
                                                     sizeof(targets) / sizeof(targets[0]), 3);
        assert_match_result(indexed[i], want);
    }
    assert(indexed[0].match_count == 4);
    assert(indexed[2].status == QDALN_MATCH_NONE);
    assert(indexed[6].match_count == 3);
    assert(indexed[7].status == QDALN_MATCH_NONE);
    assert(stats.candidates_considered == stats.candidates_verified);

    qdaln_index_free(idx);
}

static void hamming_seed_seen_hash_semantics_tests(void) {
    enum { N_TARGETS = 192, TARGET_LEN = 8 };
    char target_buf[N_TARGETS][TARGET_LEN + 1];
    const char *targets[N_TARGETS];
    size_t target_lens[N_TARGETS];
    const char *reads[] = {"AACCCCCC", "AAGGGGGG"};
    size_t read_lens[] = {TARGET_LEN, TARGET_LEN};
    qdaln_match_result indexed[2];
    qdaln_index_stats stats;

    for (size_t i = 0; i < N_TARGETS; ++i) {
        target_buf[i][0] = 'A';
        target_buf[i][1] = 'A';
        size_t v = i;
        for (size_t pos = 2; pos < TARGET_LEN; ++pos) {
            target_buf[i][pos] = "ACGT"[v & 3U];
            v >>= 2;
        }
        target_buf[i][TARGET_LEN] = '\0';
        targets[i] = target_buf[i];
        target_lens[i] = TARGET_LEN;
    }

    qdaln_index *idx = qdaln_index_build(targets, target_lens, N_TARGETS);
    assert(idx != NULL);

    assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 2, 3, indexed, &stats) == 0);
    for (size_t i = 0; i < 2; ++i) {
        qdaln_match_result want = hamming_oracle_one(reads[i], read_lens[i], targets, target_lens, N_TARGETS, 3);
        assert_match_result(indexed[i], want);
    }
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

static void index_status_shortcut_stops_unknown_after_ambiguity_tests(void) {
    char targets_buf[20][21];
    const char *targets[20];
    size_t target_lens[20];
    const char *reads[] = {"AAAAAAAAAANAAAAAAAAA"};
    size_t read_lens[] = {20};
    qdaln_match_result exhaustive[1];
    qdaln_match_result status_only[1];
    qdaln_index_stats exhaustive_stats;
    qdaln_index_stats status_stats;

    for (size_t i = 0; i < 20; ++i) {
        memset(targets_buf[i], 'A', 20);
        targets_buf[i][10] = "ACGT"[i & 3U];
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

static void index_status_shortcut_k2_stops_after_ambiguity_tests(void) {
    char targets_buf[20][21];
    const char *targets[20];
    size_t target_lens[20];
    const char *reads[] = {"ACGTACGTACGTACGTACGT"};
    size_t read_lens[] = {20};
    qdaln_match_result exhaustive[1];
    qdaln_match_result status_only[1];
    qdaln_index_stats exhaustive_stats;
    qdaln_index_stats status_stats;

    for (size_t i = 0; i < 20; ++i) {
        strcpy(targets_buf[i], reads[0]);
        targets[i] = targets_buf[i];
        target_lens[i] = read_lens[0];
    }

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 20);
    assert(idx != NULL);
    assert(qdaln_index_assign_stats(idx, reads, read_lens, 1, 2, exhaustive, &exhaustive_stats) == 0);
    assert(qdaln_index_assign_status_stats(idx, reads, read_lens, 1, 2, status_only, &status_stats) == 0);

    assert(exhaustive[0].status == QDALN_MATCH_AMBIGUOUS);
    assert(exhaustive[0].best_distance == 0);
    assert(exhaustive[0].match_count == 20);
    assert(status_only[0].status == QDALN_MATCH_AMBIGUOUS);
    assert(status_only[0].best_distance == exhaustive[0].best_distance);
    assert(status_stats.candidates_verified < exhaustive_stats.candidates_verified);
    assert(status_stats.candidates_verified <= 2);
    qdaln_index_free(idx);
}

static void index_status_shortcut_k2_unknown_stops_after_ambiguity_tests(void) {
    char targets_buf[20][21];
    const char *targets[20];
    size_t target_lens[20];
    const char *reads[] = {"ACGTACGTACNTACGTACGT"};
    size_t read_lens[] = {20};
    qdaln_match_result exhaustive[1];
    qdaln_match_result status_only[1];
    qdaln_index_stats exhaustive_stats;
    qdaln_index_stats status_stats;

    for (size_t i = 0; i < 20; ++i) {
        strcpy(targets_buf[i], "ACGTACGTACATACGTACGT");
        target_lens[i] = read_lens[0];
        targets_buf[i][10] = "ACGT"[i & 3U];
        targets[i] = targets_buf[i];
    }

    qdaln_index *idx = qdaln_index_build(targets, target_lens, 20);
    assert(idx != NULL);
    assert(qdaln_index_assign_stats(idx, reads, read_lens, 1, 2, exhaustive, &exhaustive_stats) == 0);
    assert(qdaln_index_assign_status_stats(idx, reads, read_lens, 1, 2, status_only, &status_stats) == 0);

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
        /* Extend a_len range to >64 to exercise multi-word Myers as pattern (oracle via dp) */
        size_t a_len = (size_t)(xorshift64() % 129ULL);
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

static void random_literal_sequence(char *sequence, size_t length) {
    static const char alphabet[] = "ACGTNRYacgtX-";
    for (size_t i = 0; i < length; ++i) {
        sequence[i] = alphabet[xorshift64() % (sizeof(alphabet) - 1U)];
    }
    sequence[length] = '\0';
}

static void literal_alphabet_index_fuzz_tests(void) {
    char reads_buf[8][41];
    char targets_buf[8][41];
    const char *reads[8];
    const char *targets[8];
    size_t read_lens[8];
    size_t target_lens[8];
    qdaln_match_result scan[8];
    qdaln_match_result indexed[8];

    for (size_t trial = 0; trial < 10000; ++trial) {
        size_t n_reads = 1 + (size_t)(xorshift64() % 8ULL);
        size_t n_targets = 1 + (size_t)(xorshift64() % 8ULL);
        int k = (int)(xorshift64() % 4ULL);
        for (size_t i = 0; i < n_reads; ++i) {
            read_lens[i] = (size_t)(xorshift64() % 41ULL);
            random_literal_sequence(reads_buf[i], read_lens[i]);
            reads[i] = reads_buf[i];
        }
        for (size_t i = 0; i < n_targets; ++i) {
            target_lens[i] = (size_t)(xorshift64() % 41ULL);
            random_literal_sequence(targets_buf[i], target_lens[i]);
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
    empty_index_tests();
    index_duplicate_exact_tests();
    direct_exact_lookup_tests();
    index_stats_pruning_tests();
    levenshtein_k1_avoids_false_deletion_seed_candidates_tests();
    levenshtein_k2_uses_index_without_full_scan_tests();
    levenshtein_k2_len32_insertion_uses_index_tests();
    levenshtein_k2_single_unknown_uses_bounded_candidates_tests();
    hamming_single_unknown_uses_index_tests();
    hamming_single_unknown_k2_k3_uses_bounded_candidates_tests();
    hamming_multi_unknown_k2_k3_uses_bounded_candidates_tests();
    hamming_seed_index_semantics_tests();
    hamming_k2_k3_seed_index_semantics_tests();
    hamming_seed_seen_hash_semantics_tests();
    levenshtein_non_acgt_indel_uses_index_tests();
    index_status_shortcut_stops_after_ambiguity_tests();
    index_status_shortcut_stops_unknown_after_ambiguity_tests();
    index_status_shortcut_k2_stops_after_ambiguity_tests();
    index_status_shortcut_k2_unknown_stops_after_ambiguity_tests();
    large_panel_oracle_tests();
    fuzz_tests();
    batch_fuzz_tests();
    index_fuzz_tests();
    literal_alphabet_index_fuzz_tests();
    puts("qdalign tests passed");
    return 0;
}
