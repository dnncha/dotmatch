#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

static size_t test_malloc_calls = 0;
static size_t test_calloc_calls = 0;
static size_t test_free_calls = 0;

static void *test_malloc(size_t n) {
    ++test_malloc_calls;
    return malloc(n);
}

static void *test_calloc(size_t n, size_t size) {
    ++test_calloc_calls;
    return calloc(n, size);
}

static void test_free(void *p) {
    if (p != NULL) ++test_free_calls;
    free(p);
}

#define malloc test_malloc
#define calloc test_calloc
#define free test_free
#include "../src/qdalign.c"
#undef malloc
#undef calloc
#undef free

static void reset_alloc_counts(void) {
    test_malloc_calls = 0;
    test_calloc_calls = 0;
    test_free_calls = 0;
}

static void assert_leq_no_heap(const char *a, size_t a_len, const char *b, size_t b_len, int k, int expected) {
    reset_alloc_counts();
    assert(qdaln_edit_distance_leq(a, a_len, b, b_len, k) == expected);
    assert(test_malloc_calls == 0);
    assert(test_calloc_calls == 0);
    assert(test_free_calls == 0);
}

static void assert_k1_leq_no_heap(const char *a, size_t a_len, const char *b, size_t b_len, int expected) {
    assert_leq_no_heap(a, a_len, b, b_len, 1, expected);
}

static void packed_hamming_tests(void) {
    uint64_t a = 0;
    uint64_t b = 0;

    assert(dna2_code("ACGTACGT", 8, &a) == 1);
    assert(dna2_code("ACGTTCGA", 8, &b) == 1);
    assert(code_hamming_distance_qd(a, b, 8) == 2);
    assert(same_length_hamming_distance_within_k("ACGTACGTACGTACGT",
                                                 "ACGTACGTACGTACGA", 16, 1) == 1);
    assert(same_length_hamming_distance_within_k("ACGTACGTACGTACGT",
                                                 "TCGTACGTACGTACGA", 16, 1) == -1);
}

static void hamming_seed_index_build_tests(void) {
    const char *targets[] = {"ACGTACGT", "ACGTTCGT", "A", "NNNN"};
    size_t target_lens[] = {8, 8, 1, 4};
    qdaln_index *idx = qdaln_index_build(targets, target_lens, 4);

    assert(idx != NULL);
    assert(idx->hamming_seed_ready == 1);
    assert(idx->n_hamming_seeds == 27);
    assert(idx->hamming_seed_hash_cap >= 8);
    assert(idx->hamming_seed_heads[0] == 0);
    assert(idx->hamming_seed_heads[idx->hamming_seed_hash_cap] == idx->n_hamming_seeds);
    for (size_t slot = 0; slot < idx->hamming_seed_hash_cap; ++slot) {
        assert(idx->hamming_seed_heads[slot] <= idx->hamming_seed_heads[slot + 1]);
        for (size_t e = idx->hamming_seed_heads[slot]; e < idx->hamming_seed_heads[slot + 1]; ++e) {
            const qdaln_hamming_seed_entry *seed = &idx->hamming_seeds[e];
            assert(hamming_seed_hash(seed->code, seed->seed_len, seed->seed_id,
                                     idx->hamming_seed_hash_cap) == slot);
            assert(seed->target_code == idx->codes[seed->target_index]);
            assert(seed->target_len == idx->target_lens[seed->target_index]);
        }
    }
    qdaln_index_free(idx);
}

static void hamming_dense_query_no_heap_tests(void) {
    enum { N_TARGETS = 512, LEN = 8 };
    char words[N_TARGETS][LEN + 1];
    const char *targets[N_TARGETS];
    size_t lens[N_TARGETS];
    const char *reads[] = {"AACCCCCC", "AAGGGGGG", "AAAAAAAA"};
    size_t read_lens[] = {LEN, LEN, LEN};
    qdaln_match_result results[3];
    for (size_t i = 0; i < N_TARGETS; ++i) {
        words[i][0] = words[i][1] = 'A';
        size_t v = i;
        for (size_t p = 2; p < LEN; ++p) {
            words[i][p] = "ACGT"[v & 3U];
            v >>= 2;
        }
        words[i][LEN] = '\0';
        targets[i] = words[i];
        lens[i] = LEN;
    }
    words[N_TARGETS - 1][LEN - 1] = 'N';
    qdaln_index *idx = qdaln_index_build(targets, lens, N_TARGETS);
    assert(idx != NULL);
    for (int k = 1; k <= 3; ++k) {
        reset_alloc_counts();
        qdaln_index_stats stats;
        assert(qdaln_index_assign_hamming_stats(idx, reads, read_lens, 3, k, results, &stats) == 0);
        assert(test_malloc_calls == 0);
        assert(test_calloc_calls == 0);
        assert(test_free_calls == 0);
        assert(stats.candidates_considered == stats.candidates_verified);
        assert(results[2].best_distance == 0);
    }
    qdaln_index_free(idx);
}

int main(void) {
    packed_hamming_tests();
    hamming_seed_index_build_tests();
    hamming_dense_query_no_heap_tests();
    assert_k1_leq_no_heap("ACGT", 4, "ACGT", 4, 1);
    assert_k1_leq_no_heap("ACGT", 4, "ACGA", 4, 1);
    assert_k1_leq_no_heap("ACGT", 4, "ACGTT", 5, 1);
    assert_k1_leq_no_heap("ACGT", 4, "ACG", 3, 1);
    assert_k1_leq_no_heap("ACGT", 4, "TGCA", 4, 0);
    assert_k1_leq_no_heap("ACGT", 4, "ACGTTT", 6, 0);
    assert_k1_leq_no_heap("N", 1, "A", 1, 1);
    assert_k1_leq_no_heap("N", 1, "N", 1, 1);
    assert_leq_no_heap("ACGTACGT", 8, "ACGTTCGA", 8, 2, 1);
    assert_leq_no_heap("ACGTACGT", 8, "TCGTTCGA", 8, 3, 1);
    assert_leq_no_heap("ACG", 3, "CGA", 3, 2, 1);

    /* Extend coverage to longer patterns (>64) with oracle comparison.
     * For k>=2, leq now routes to generalized multi-word myers (no heap),
     * and we assert exact match to dp oracle. */
    {
        const char *longa = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; /* 72bp */
        const char *longb = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA"; /* 1 mismatch */
        size_t lena = 72, lenb = 72;
        assert_leq_no_heap(longa, lena, longb, lenb, 2, 1);
        assert(qdaln_edit_distance_myers64(longa, lena, longb, lenb) == qdaln_edit_distance_dp(longa, lena, longb, lenb));
        assert(qdaln_edit_distance(longa, lena, longb, lenb) == qdaln_edit_distance_dp(longa, lena, longb, lenb));
        assert(qdaln_edit_distance_leq(longa, lena, longb, lenb, 0) == 0);
        assert(qdaln_edit_distance_leq(longa, lena, longb, lenb, 1) == 1);
    }

    puts("qdalign threshold allocation tests passed");
    return 0;
}
