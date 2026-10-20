/* Exercise ownership on the allocation-error paths that GCC 13 misdiagnoses. */
#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE 1
#define _GNU_SOURCE 1
#include <assert.h>
#include <stdlib.h>

static int fail_malloc = 0;
static int fail_realloc = 0;
static size_t live_allocations = 0;

static void *checked_malloc(size_t size) {
    if (fail_malloc) return NULL;
    void *result = malloc(size);
    if (result != NULL) ++live_allocations;
    return result;
}

static void *checked_realloc(void *pointer, size_t size) {
    if (fail_realloc) return NULL;
    int new_allocation = pointer == NULL;
    void *result = realloc(pointer, size);
    if (result != NULL && new_allocation) ++live_allocations;
    return result;
}

static void checked_free(void *pointer) {
    if (pointer != NULL) {
        assert(live_allocations > 0);
        --live_allocations;
    }
    free(pointer);
}

#define malloc checked_malloc
#define realloc checked_realloc
#define free checked_free
#define main dotmatch_cli_main
#include "../src/qda.c"
#undef main
#undef free
#undef realloc
#undef malloc

int main(void) {
    string_list list = {0};
    fail_realloc = 1;
    assert(push_string(&list, "first") == -1);
    assert(list.count == 0 && list.items == NULL);
    fail_realloc = 0;

    assert(push_string(&list, "first") == 0);
    fail_malloc = 1;
    assert(push_string(&list, "not stored") == -1);
    assert(list.count == 1 && strcmp(list.items[0], "first") == 0);
    fail_malloc = 0;

    /* Fill the initial capacity, then fail growth without losing its strings. */
    while (list.count < list.cap) assert(push_string(&list, "retained") == 0);
    size_t count = list.count;
    fail_realloc = 1;
    assert(push_string(&list, "not stored") == -1);
    assert(list.count == count && strcmp(list.items[0], "first") == 0);
    fail_realloc = 0;
    assert(push_string(&list, "after growth") == 0);

    free_string_list(&list);
    assert(list.count == 0 && list.cap == 0 && list.items == NULL);
    assert(live_allocations == 0);
    free_string_list(&list);
    puts("native string-list allocation failure tests passed");
    return 0;
}
