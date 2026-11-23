#include "qdalign.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint64_t rng_state;
static uint64_t rnd(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}
static size_t rn(size_t n) { return n ? (size_t)(rnd() % n) : 0; }

static int ref_lev(const char *a, size_t n, const char *b, size_t m) {
    int *d = malloc((n + 1) * (m + 1) * sizeof(int));
    for (size_t i = 0; i <= n; ++i) d[i * (m + 1)] = (int)i;
    for (size_t j = 0; j <= m; ++j) d[j] = (int)j;
    for (size_t i = 1; i <= n; ++i)
        for (size_t j = 1; j <= m; ++j) {
            int s = d[(i - 1) * (m + 1) + j - 1] + ((unsigned char)a[i - 1] != (unsigned char)b[j - 1]);
            int x = d[(i - 1) * (m + 1) + j] + 1, y = d[i * (m + 1) + j - 1] + 1;
            int r = s < x ? s : x;
            d[i * (m + 1) + j] = r < y ? r : y;
        }
    int r = d[n * (m + 1) + m];
    free(d);
    return r;
}

static int ref_ham(const char *a, size_t n, const char *b, size_t m) {
    if (n != m) return 1 << 20;
    int d = 0;
    for (size_t i = 0; i < n; ++i) d += a[i] != b[i];
    return d;
}

static char fold_ascii(char c) {
    return c >= 'a' && c <= 'z' ? (char)(c - ('a' - 'A')) : c;
}

static int ref_exact_ascii(const char *a, size_t n, const char *b, size_t m) {
    if (n != m) return 1 << 20;
    for (size_t i = 0; i < n; ++i)
        if (fold_ascii(a[i]) != fold_ascii(b[i])) return 1 << 20;
    return 0;
}

typedef int (*metric_fn)(const char *, size_t, const char *, size_t);

static qdaln_match_result ref_match(const char *r, size_t rl, char **t, size_t *tl, size_t nt, int k, metric_fn f) {
    qdaln_match_result res = {-1, -1, -1, 0, QDALN_MATCH_NONE};
    int ties = 0;
    for (size_t j = 0; j < nt; ++j) {
        int d = f(r, rl, t[j], tl[j]);
        if (d > k) continue;
        res.match_count++;
        if (res.best_distance < 0 || d < res.best_distance) {
            res.second_best_distance = res.best_distance;
            res.best_distance = d;
            res.target_index = (int)j;
            ties = 1;
        } else if (d == res.best_distance) {
            ties++;
        } else if (res.second_best_distance < 0 || d < res.second_best_distance) {
            res.second_best_distance = d;
        }
    }
    if (res.match_count) res.status = ties > 1 ? QDALN_MATCH_AMBIGUOUS : QDALN_MATCH_UNIQUE;
    return res;
}

static const char ALPHA[] = "ACGTACGTACGTACGTACGTNacgtn\x00\xff-";
static char rand_base(void) {
    size_t x = rn(100);
    if (x < 88) return "ACGT"[rn(4)];
    return ALPHA[rn(sizeof(ALPHA) - 1)];
}

static long long_every;
static size_t rand_len(void) {
    size_t x = rn(1000);
    if (x < 5) return 0;
    if (x < 600) return 1 + rn(12);
    if (x < 950) return 13 + rn(24);
    if (x < 995) return 37 + rn(100);
    return 500 + rn(600);
}

static void mutate(const char *src, size_t sl, char *dst, size_t *dl, int edits) {
    memcpy(dst, src, sl);
    size_t l = sl;
    for (int e = 0; e < edits; ++e) {
        int kind = (int)rn(3);
        if (kind == 0 && l > 0) dst[rn(l)] = rand_base();
        else if (kind == 1) {
            size_t p = rn(l + 1);
            memmove(dst + p + 1, dst + p, l - p);
            dst[p] = rand_base();
            ++l;
        } else if (l > 0) {
            size_t p = rn(l);
            memmove(dst + p, dst + p + 1, l - p - 1);
            --l;
        }
    }
    *dl = l;
}

static long failures;
static const char *only;
static const char *cat_names[64];
static long cat_counts[64];
static int n_cats;
static int note(const char *cat) {
    if (only && strcmp(only, cat) != 0) return 0;
    ++failures;
    int i = 0;
    while (i < n_cats && strcmp(cat_names[i], cat) != 0) ++i;
    if (i == n_cats) cat_names[n_cats++] = cat;
    return ++cat_counts[i] <= 3;
}

static void dump(const char *s, size_t l) {
    putchar('"');
    for (size_t i = 0; i < l; ++i) {
        unsigned char c = (unsigned char)s[i];
        if (c >= 32 && c < 127) putchar(c);
        else printf("\\x%02x", c);
    }
    putchar('"');
}

static void report_case(const char *what, const char *r, size_t rl, char **t, size_t *tl, size_t nt, int k,
                        qdaln_match_result got, qdaln_match_result exp) {
    if (!note(what)) return;
    printf("MISMATCH %s k=%d read=", what, k);
    dump(r, rl);
    printf(" targets=[");
    for (size_t j = 0; j < nt; ++j) {
        if (j) printf(",");
        dump(t[j], tl[j]);
    }
    printf("]\n  got  idx=%d best=%d second=%d count=%d status=%d\n", got.target_index, got.best_distance,
           got.second_best_distance, got.match_count, got.status);
    printf("  want idx=%d best=%d second=%d count=%d status=%d\n", exp.target_index, exp.best_distance,
           exp.second_best_distance, exp.match_count, exp.status);
}

static int same_full(qdaln_match_result a, qdaln_match_result b) {
    return a.target_index == b.target_index && a.best_distance == b.best_distance &&
           a.second_best_distance == b.second_best_distance && a.match_count == b.match_count && a.status == b.status;
}

static int same_status(qdaln_match_result a, qdaln_match_result b) {
    if (a.status != b.status || a.best_distance != b.best_distance) return 0;
    if (a.status == QDALN_MATCH_UNIQUE && a.target_index != b.target_index) return 0;
    return 1;
}

static int ref_edit_class(size_t rl, size_t tl, int d) {
    if (d == 0) return QDALN_EDIT_EXACT;
    if (d == 1) {
        if (rl == tl) return QDALN_EDIT_K1_SUB;
        if (rl == tl + 1) return QDALN_EDIT_K1_INS;
        if (tl == rl + 1) return QDALN_EDIT_K1_DEL;
    }
    if (d == 2) return QDALN_EDIT_K2;
    return QDALN_EDIT_OTHER;
}

int main(int argc, char **argv) {
    long iters = argc > 1 ? strtol(argv[1], NULL, 10) : 100000;
    rng_state = argc > 2 ? strtoull(argv[2], 0, 10) : 0x1234567ULL;
    only = argc > 3 ? argv[3] : NULL;
    (void)long_every;
    if (rng_state == 0) rng_state = 1;
    static char buf[16][4096];
    char *t[16];
    size_t tl[16];
    for (int i = 0; i < 16; ++i) t[i] = buf[i];
    static char rb[4096];

    for (long it = 0; it < iters; ++it) {
        size_t base_len = rand_len();
        char base[2048];
        for (size_t i = 0; i < base_len; ++i) base[i] = rand_base();
        size_t nt = rn(7);
        for (size_t j = 0; j < nt; ++j) {
            if (j > 0 && rn(8) == 0) {
                size_t src = rn(j);
                memcpy(t[j], t[src], tl[src]);
                tl[j] = tl[src];
            } else {
                mutate(base, base_len, t[j], &tl[j], (int)rn(4));
            }
        }
        size_t rl;
        mutate(base, base_len, rb, &rl, (int)rn(4));
        const char *read = rb;
        int k = (int)rn(5);
        if (rn(50) == 0) k = (int)(rl + rn(3));

        for (size_t j = 0; j < nt; ++j) {
            int ref = ref_lev(read, rl, t[j], tl[j]);
            int a = qdaln_edit_distance(read, rl, t[j], tl[j]);
            int b = qdaln_edit_distance_dp(read, rl, t[j], tl[j]);
            int c = qdaln_edit_distance_myers64(read, rl, t[j], tl[j]);
            int c2 = qdaln_edit_distance_myers64(t[j], tl[j], read, rl);
            if ((a != ref || b != ref || c != ref || c2 != ref) && note("dist")) {
                printf("MISMATCH dist ref=%d ed=%d dp=%d my=%d my2=%d a=", ref, a, b, c, c2);
                dump(read, rl);
                printf(" b=");
                dump(t[j], tl[j]);
                printf("\n");
            }
            for (int kk = -1; kk <= 4; ++kk) {
                int l = qdaln_edit_distance_leq(read, rl, t[j], tl[j], kk);
                int want = ref <= kk ? 1 : 0;
                if (l != want && note("leq")) {
                    printf("MISMATCH leq k=%d got=%d want=%d a=", kk, l, want);
                    dump(read, rl);
                    printf(" b=");
                    dump(t[j], tl[j]);
                    printf("\n");
                }
            }
        }

        const char *rp[1] = {read};
        size_t rlp[1] = {rl};
        qdaln_match_result exp = ref_match(read, rl, t, tl, nt, k, ref_lev);
        qdaln_match_result got;
        qdaln_match_many(rp, rlp, 1, (const char *const *)t, tl, nt, k, &got);
        if (!same_full(got, exp)) report_case("match_many", read, rl, t, tl, nt, k, got, exp);

        for (int policy = 0; policy < 2; ++policy) {
            qdaln_assignment_result ar;
            qdaln_assign_many(rp, rlp, 1, (const char *const *)t, tl, nt, k, policy, &ar);
            int want_status = exp.status;
            if (policy == 1 && exp.match_count > 1) want_status = QDALN_MATCH_AMBIGUOUS;
            int ties = 0;
            for (size_t j = 0; j < nt; ++j) ties += exp.match_count && ref_lev(read, rl, t[j], tl[j]) == exp.best_distance;
            int want_class = exp.match_count ? ref_edit_class(rl, tl[exp.target_index], exp.best_distance) : QDALN_EDIT_NONE;
            if ((ar.status != want_status || ar.distance != exp.best_distance || ar.target_index != exp.target_index ||
                 ar.num_targets_within_radius != exp.match_count || ar.num_best_targets != ties ||
                 ar.second_best_distance != exp.second_best_distance || ar.edit_class != want_class) &&
                note("assign_many")) {
                printf("MISMATCH assign_many policy=%d k=%d status %d/%d dist %d/%d idx %d/%d n %d/%d ties %d/%d class %d/%d read=",
                       policy, k, ar.status, want_status, ar.distance, exp.best_distance, ar.target_index,
                       exp.target_index, ar.num_targets_within_radius, exp.match_count, ar.num_best_targets, ties,
                       ar.edit_class, want_class);
                dump(read, rl);
                printf("\n");
            }
        }

        qdaln_index *idx = qdaln_index_build((const char *const *)t, tl, nt);
        if (!idx) {
            printf("index build failed\n");
            return 1;
        }
        qdaln_index_stats st;
        qdaln_index_assign_stats(idx, rp, rlp, 1, k, &got, &st);
        if (!same_full(got, exp)) report_case("index_assign", read, rl, t, tl, nt, k, got, exp);
        qdaln_index_assign_status_stats(idx, rp, rlp, 1, k, &got, &st);
        if (!same_status(got, exp)) report_case("index_status", read, rl, t, tl, nt, k, got, exp);

        for (int hk = 0; hk <= 3; ++hk) {
            qdaln_match_result hexp = ref_match(read, rl, t, tl, nt, hk, ref_ham);
            qdaln_index_assign_hamming_stats(idx, rp, rlp, 1, hk, &got, &st);
            if (!same_full(got, hexp)) report_case("hamming", read, rl, t, tl, nt, hk, got, hexp);
        }
        qdaln_match_result eexp = ref_match(read, rl, t, tl, nt, 0, ref_ham);
        qdaln_index_lookup_exact_stats(idx, read, rl, &got, &st);
        if (!same_full(got, eexp)) report_case("lookup_exact", read, rl, t, tl, nt, 0, got, eexp);
        qdaln_index_lookup_exact_many_stats(idx, rp, rlp, 1, &got, &st);
        if (!same_full(got, eexp)) report_case("lookup_exact_many", read, rl, t, tl, nt, 0, got, eexp);
        qdaln_match_result aexp = ref_match(read, rl, t, tl, nt, 0, ref_exact_ascii);
        qdaln_index_lookup_exact_ascii_stats(idx, read, rl, &got, &st);
        if (!same_full(got, aexp)) report_case("lookup_ascii", read, rl, t, tl, nt, 0, got, aexp);
        qdaln_index_lookup_exact_ascii_many_stats(idx, rp, rlp, 1, &got, &st);
        if (!same_full(got, aexp)) report_case("lookup_ascii_many", read, rl, t, tl, nt, 0, got, aexp);
        qdaln_index_free(idx);
    }
    for (int i = 0; i < n_cats; ++i) printf("CATEGORY %s %ld\n", cat_names[i], cat_counts[i]);
    printf("iters=%ld failures=%ld\n", iters, failures);
    return failures != 0;
}
