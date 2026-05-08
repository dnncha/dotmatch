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

static void assert_k1_leq_no_heap(const char *a, size_t a_len, const char *b, size_t b_len, int expected) {
    reset_alloc_counts();
    assert(qdaln_edit_distance_leq(a, a_len, b, b_len, 1) == expected);
    assert(test_malloc_calls == 0);
    assert(test_calloc_calls == 0);
    assert(test_free_calls == 0);
}

int main(void) {
    assert_k1_leq_no_heap("ACGT", 4, "ACGT", 4, 1);
    assert_k1_leq_no_heap("ACGT", 4, "ACGA", 4, 1);
    assert_k1_leq_no_heap("ACGT", 4, "ACGTT", 5, 1);
    assert_k1_leq_no_heap("ACGT", 4, "ACG", 3, 1);
    assert_k1_leq_no_heap("ACGT", 4, "TGCA", 4, 0);
    assert_k1_leq_no_heap("ACGT", 4, "ACGTTT", 6, 0);
    assert_k1_leq_no_heap("N", 1, "A", 1, 1);
    assert_k1_leq_no_heap("N", 1, "N", 1, 1);

    puts("qdalign threshold allocation tests passed");
    return 0;
}
