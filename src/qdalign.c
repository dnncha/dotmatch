#include "qdalign.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(__AVX2__)
#include <immintrin.h>
#elif defined(__ARM_NEON) || defined(__aarch64__)
#include <arm_neon.h>
#endif

/* Max words for stack-allocated multi-word Myers bitvectors.
 * Supports patterns up to 512 bp with the fast bit-parallel path.
 * (Roadmap target was 128; 512 is safe and sufficient for primers/amplicons.)
 */
#define QDALN_MYERS_MAX_WORDS 8

typedef struct qdaln_hamming_seed_entry {
    uint64_t code;
    size_t target_index;
    size_t seed_len;
    unsigned char seed_id;
    size_t next;
} qdaln_hamming_seed_entry;

struct qdaln_index {
    char **targets;
    size_t *target_lens;
    uint64_t *codes;
    unsigned char *encodable;
    size_t *hash_heads;
    size_t *hash_next;
    uint64_t *exact_codes;
    size_t *exact_lens;
    size_t *exact_first_targets;
    size_t *exact_counts;
    size_t exact_hash_cap;
    uint64_t *deletion_codes;
    size_t *deletion_lens;
    size_t *deletion_targets;
    size_t *deletion_heads;
    size_t *deletion_next;
    size_t deletion_hash_cap;
    size_t n_deletions;
    qdaln_hamming_seed_entry *hamming_seeds;
    size_t *hamming_seed_heads;
    size_t hamming_seed_hash_cap;
    size_t n_hamming_seeds;
    int hamming_seed_ready;
    size_t *nonencodable_targets;
    size_t hash_cap;
    size_t n_targets;
    size_t n_nonencodable;
};

#define QDALN_HASH_EMPTY SIZE_MAX

static inline int min3_int(int x, int y, int z) {
    int m = x < y ? x : y;
    return m < z ? m : z;
}

static int dna2_code(const char *s, size_t len, uint64_t *code_out) {
    if ((s == NULL && len != 0) || len > 32) return 0;
    uint64_t code = 0;
    for (size_t i = 0; i < len; ++i) {
        uint64_t v;
        switch (s[i]) {
            case 'A':
                v = 0;
                break;
            case 'C':
                v = 1;
                break;
            case 'G':
                v = 2;
                break;
            case 'T':
                v = 3;
                break;
            default:
                return 0;
        }
        code |= v << (2 * i);
    }
    *code_out = code;
    return 1;
}

static int dna2_code_ascii_upper(const char *s, size_t len, uint64_t *code_out) {
    if ((s == NULL && len != 0) || len > 32) return 0;
    uint64_t code = 0;
    for (size_t i = 0; i < len; ++i) {
        uint64_t v;
        switch (s[i]) {
            case 'A':
            case 'a':
                v = 0;
                break;
            case 'C':
            case 'c':
                v = 1;
                break;
            case 'G':
            case 'g':
                v = 2;
                break;
            case 'T':
            case 't':
                v = 3;
                break;
            default:
                return 0;
        }
        code |= v << (2 * i);
    }
    *code_out = code;
    return 1;
}

static size_t next_pow2_size(size_t n) {
    size_t x = 1;
    while (x < n) x <<= 1;
    return x;
}

static inline size_t code_hash(uint64_t code, size_t len, size_t cap) {
    uint64_t x = code ^ ((uint64_t)len * 0x9e3779b97f4a7c15ULL);
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    return (size_t)x & (cap - 1);
}

static inline size_t exact_code_hash(uint64_t code, size_t len, size_t cap) {
    uint64_t x = code * UINT64_C(0x9e3779b97f4a7c15);
    x ^= (uint64_t)len * UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 32;
    return (size_t)x & (cap - 1);
}

static inline size_t hamming_seed_hash(uint64_t code, size_t seed_len, unsigned char seed_id, size_t cap) {
    return code_hash(code ^ ((uint64_t)seed_id * 0x517cc1b727220a95ULL),
                     seed_len + (size_t)seed_id * 37U, cap);
}

const char *qdaln_alphabet_policy(void) {
    return QDALN_ALPHABET_POLICY;
}

static uint64_t code_delete_base(uint64_t code, size_t len, size_t drop_pos) {
    if (drop_pos + 1 >= len) {
        uint64_t keep_bits = (uint64_t)2 * drop_pos;
        return keep_bits == 0 ? 0 : (code & ((1ULL << keep_bits) - 1ULL));
    }
    uint64_t lower_bits = (uint64_t)2 * drop_pos;
    uint64_t lower = lower_bits == 0 ? 0 : (code & ((1ULL << lower_bits) - 1ULL));
    uint64_t upper = code >> ((uint64_t)2 * (drop_pos + 1));
    return lower | (upper << lower_bits);
}

static uint64_t code_insert_base(uint64_t code, size_t len, size_t insert_pos, uint64_t base) {
    uint64_t lower_bits = (uint64_t)2 * insert_pos;
    uint64_t lower = lower_bits == 0 ? 0 : (code & ((1ULL << lower_bits) - 1ULL));
    uint64_t upper = code >> lower_bits;
    uint64_t shifted_upper = lower_bits + 2 >= 64 ? 0 : (upper << (lower_bits + 2));
    (void)len;
    return lower | ((base & 3ULL) << lower_bits) | shifted_upper;
}

static uint64_t code_low_mask_qd(size_t len) {
    if (len == 0) return 0;
    if (len >= 32) return UINT64_MAX;
    return (1ULL << (2 * len)) - 1ULL;
}

static uint64_t code_segment_qd(uint64_t code, size_t start, size_t len) {
    return (code >> (2 * start)) & code_low_mask_qd(len);
}

static void hamming_seed_partition(size_t total_len, size_t n_seeds, size_t seed_id,
                                   size_t *start_out, size_t *len_out) {
    size_t base = total_len / n_seeds;
    size_t rem = total_len % n_seeds;
    size_t seed_len = base + (seed_id < rem ? 1U : 0U);
    size_t start = seed_id * base + (seed_id < rem ? seed_id : rem);
    *start_out = start;
    *len_out = seed_len;
}

static unsigned char hamming_seed_key(size_t n_seeds, size_t seed_id) {
    return (unsigned char)((n_seeds - 2U) * 4U + seed_id);
}

static int popcount64_qd(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    int count = 0;
    while (x != 0) {
        x &= x - 1;
        ++count;
    }
    return count;
#endif
}

static int code_hamming_distance_qd(uint64_t a, uint64_t b, size_t len) {
    uint64_t diff = a ^ b;
    diff |= (diff >> 1);
    diff &= (code_low_mask_qd(len) & 0x5555555555555555ULL);
    return popcount64_qd(diff);
}

/* Optimized nonzero byte count for the 8-byte fast path (used by same_length hamming).
   Keeps exact same semantics. */
static int nonzero_bytes_u64(uint64_t x) {
    x |= x >> 4;
    x |= x >> 2;
    x |= x >> 1;
    x &= 0x0101010101010101ULL;
    return popcount64_qd(x);
}

/* SIMD-accelerated same-length Hamming (bytewise diff count) within k.
   Uses 256-bit AVX2 on x86 (when __AVX2__ defined, e.g. -mavx2) or 128-bit NEON on AArch64.
   Falls back to original 64-bit chunk + scalar for other/ tails.
   Semantics identical to original: returns d if d<=k else -1, exact distances for tests/DP crosschecks. */
static int same_length_hamming_distance_within_k(const char *a, const char *b, size_t len, int k) {
#if defined(__AVX2__)
    int d = 0;
    size_t i = 0;
    for (; i + 32 <= len; i += 32) {
        __m256i va = _mm256_loadu_si256((const __m256i *)(a + i));
        __m256i vb = _mm256_loadu_si256((const __m256i *)(b + i));
        __m256i diff = _mm256_xor_si256(va, vb);
        __m256i zero = _mm256_setzero_si256();
        __m256i eq = _mm256_cmpeq_epi8(diff, zero); /* 0xff where bytes equal */
        int eqmask = _mm256_movemask_epi8(eq);
        unsigned eqm = (unsigned)eqmask & 0xFFFFFFFFU;
        d += 32 - __builtin_popcount(eqm);
        if (d > k) return -1;
    }
    /* 8-byte chunks for remainder before byte tail */
    while (i + 8 <= len) {
        uint64_t wa = 0;
        uint64_t wb = 0;
        memcpy(&wa, a + i, sizeof(wa));
        memcpy(&wb, b + i, sizeof(wb));
        uint64_t diff = wa ^ wb;
        if (diff != 0) {
            d += nonzero_bytes_u64(diff);
            if (d > k) return -1;
        }
        i += 8;
    }
    for (; i < len; ++i) {
        if (a[i] != b[i] && ++d > k) return -1;
    }
    return d;
#elif defined(__aarch64__)
    int d = 0;
    size_t i = 0;
    for (; i + 16 <= len; i += 16) {
        uint8x16_t va = vld1q_u8((const uint8_t *)(a + i));
        uint8x16_t vb = vld1q_u8((const uint8_t *)(b + i));
        uint8x16_t diff = veorq_u8(va, vb);
        uint8x16_t nz = vcgtq_u8(diff, vdupq_n_u8(0));
        uint8x16_t one = vshrq_n_u8(nz, 7); /* 0 or 1 per byte lane */
        d += (int)vaddvq_u8(one);
        if (d > k) return -1;
    }
    /* 8-byte chunks for remainder */
    while (i + 8 <= len) {
        uint64_t wa = 0;
        uint64_t wb = 0;
        memcpy(&wa, a + i, sizeof(wa));
        memcpy(&wb, b + i, sizeof(wb));
        uint64_t diff = wa ^ wb;
        if (diff != 0) {
            d += nonzero_bytes_u64(diff);
            if (d > k) return -1;
        }
        i += 8;
    }
    for (; i < len; ++i) {
        if (a[i] != b[i] && ++d > k) return -1;
    }
    return d;
#else
    /* Scalar fallback (original impl) */
    int d = 0;
    size_t i = 0;
#if defined(__AVX2__)
    for (; i + 32 <= len; i += 32) {
        __m256i va = _mm256_loadu_si256((const __m256i *)(a + i));
        __m256i vb = _mm256_loadu_si256((const __m256i *)(b + i));
        __m256i diff = _mm256_xor_si256(va, vb);
        d += count_nonzero_bytes_avx2(diff);
        if (d > k) return -1;
    }
#elif defined(__ARM_NEON) || defined(__aarch64__)
    for (; i + 16 <= len; i += 16) {
        uint8x16_t va = vld1q_u8((const uint8_t *)(a + i));
        uint8x16_t vb = vld1q_u8((const uint8_t *)(b + i));
        uint8x16_t diff = veorq_u8(va, vb);
        uint8x16_t nz = vcgtq_u8(diff, vdupq_n_u8(0));
        uint8x16_t one = vshrq_n_u8(nz, 7);
        d += (int)vaddvq_u8(one);
        if (d > k) return -1;
    }
    for (; i + 8 <= len; i += 8) {
        uint64_t wa = 0;
        uint64_t wb = 0;
        memcpy(&wa, a + i, 8);
        memcpy(&wb, b + i, 8);
        uint64_t diff = wa ^ wb;
        if (diff != 0) {
            d += nonzero_bytes_u64(diff);
            if (d > k) return -1;
        }
    }
#else
    while (i + sizeof(uint64_t) <= len) {
        uint64_t wa = 0;
        uint64_t wb = 0;
        memcpy(&wa, a + i, sizeof(wa));
        memcpy(&wb, b + i, sizeof(wb));
        uint64_t diff = wa ^ wb;
        if (diff != 0) {
            d += nonzero_bytes_u64(diff);
            if (d > k) return -1;
        }
        i += sizeof(uint64_t);
    }
#endif
    for (; i < len; ++i) {
        if (a[i] != b[i] && ++d > k) return -1;
    }
    return d;
#endif
}

int qdaln_edit_distance_dp(const char *a, size_t a_len, const char *b, size_t b_len) {
    if ((a == NULL && a_len != 0) || (b == NULL && b_len != 0)) return -1;
    if (a_len > (size_t)INT32_MAX || b_len > (size_t)INT32_MAX) return -1;

    if (a_len == 0) return (int)b_len;
    if (b_len == 0) return (int)a_len;

    int *prev = (int *)malloc((b_len + 1) * sizeof(int));
    int *curr = (int *)malloc((b_len + 1) * sizeof(int));
    if (prev == NULL || curr == NULL) {
        free(prev);
        free(curr);
        return -1;
    }

    for (size_t j = 0; j <= b_len; ++j) prev[j] = (int)j;

    for (size_t i = 1; i <= a_len; ++i) {
        curr[0] = (int)i;
        unsigned char ca = (unsigned char)a[i - 1];
        for (size_t j = 1; j <= b_len; ++j) {
            int cost = ca == (unsigned char)b[j - 1] ? 0 : 1;
            int del = prev[j] + 1;
            int ins = curr[j - 1] + 1;
            int sub = prev[j - 1] + cost;
            curr[j] = min3_int(del, ins, sub);
        }
        int *tmp = prev;
        prev = curr;
        curr = tmp;
    }

    int result = prev[b_len];
    free(prev);
    free(curr);
    return result;
}

int qdaln_edit_distance_myers64(const char *pattern, size_t pattern_len,
                                const char *text, size_t text_len) {
    if ((pattern == NULL && pattern_len != 0) || (text == NULL && text_len != 0)) return -1;
    if (pattern_len == 0) return (int)text_len;
    if (text_len == 0) return (int)pattern_len;

    int nwords = (int)((pattern_len + 63) / 64);
    if (nwords > QDALN_MYERS_MAX_WORDS) {
        /* Extreme length: fall back to DP (still correct, though slower). */
        return qdaln_edit_distance_dp(pattern, pattern_len, text, text_len);
    }
    uint64_t peq[256][QDALN_MYERS_MAX_WORDS];
    memset(peq, 0, sizeof(peq));

    for (size_t i = 0; i < pattern_len; ++i) {
        unsigned char c = (unsigned char)pattern[i];
        int w = (int)(i / 64);
        int b = (int)(i % 64);
        peq[c][w] |= (uint64_t)1 << b;
    }

    uint64_t valid_masks[QDALN_MYERS_MAX_WORDS];
    for (int w = 0; w < nwords; ++w) {
        if ((w + 1) * 64 <= (int)pattern_len) {
            valid_masks[w] = UINT64_MAX;
        } else {
            int rem = (int)(pattern_len % 64);
            valid_masks[w] = (rem == 0 ? UINT64_MAX : (((uint64_t)1 << rem) - 1));
        }
    }

    int last_w = (int)((pattern_len - 1) / 64);
    int last_b = (int)((pattern_len - 1) % 64);
    uint64_t last_bit = (uint64_t)1 << last_b;

    uint64_t pv[QDALN_MYERS_MAX_WORDS];
    uint64_t mv[QDALN_MYERS_MAX_WORDS];
    for (int w = 0; w < nwords; ++w) {
        pv[w] = valid_masks[w];
        mv[w] = 0;
    }
    int score = (int)pattern_len;

    /* temps to avoid repeated stack in loop */
    uint64_t eqq[QDALN_MYERS_MAX_WORDS];
    uint64_t xvv[QDALN_MYERS_MAX_WORDS];
    uint64_t xhh[QDALN_MYERS_MAX_WORDS];
    uint64_t phh[QDALN_MYERS_MAX_WORDS];
    uint64_t mhh[QDALN_MYERS_MAX_WORDS];
    uint64_t suml[QDALN_MYERS_MAX_WORDS];

    for (size_t j = 0; j < text_len; ++j) {
        unsigned char c = (unsigned char)text[j];
        for (int w = 0; w < nwords; ++w) {
            eqq[w] = peq[c][w];
        }

        for (int w = 0; w < nwords; ++w) {
            xvv[w] = (eqq[w] | mv[w]) & valid_masks[w];
        }

        /* (eq & pv) + pv  with carry propagation across words (carry can be 0-2).
           Portable carry calc avoids __int128 for broad compiler compat (MSVC etc). */
        uint64_t carry = 0;
        for (int w = 0; w < nwords; ++w) {
            uint64_t eqpv = eqq[w] & pv[w];
            uint64_t lo = eqpv + pv[w];
            uint64_t c1 = (lo < eqpv) ? 1ULL : 0ULL;
            lo += carry;
            uint64_t c2 = (lo < carry) ? 1ULL : 0ULL;
            suml[w] = lo;
            carry = c1 + c2;
        }
        for (int w = 0; w < nwords; ++w) {
            xhh[w] = (((suml[w] ^ pv[w]) | eqq[w]) & valid_masks[w]);
        }

        for (int w = 0; w < nwords; ++w) {
            uint64_t not_term = (~(xhh[w] | pv[w])) & valid_masks[w];
            phh[w] = (mv[w] | not_term) & valid_masks[w];
        }
        for (int w = 0; w < nwords; ++w) {
            mhh[w] = (pv[w] & xhh[w]) & valid_masks[w];
        }

        /* score delta from the MSB of the *pattern* (last bit) before the shift */
        if (phh[last_w] & last_bit) {
            ++score;
        } else if (mhh[last_w] & last_bit) {
            --score;
        }

        /* ph = ((ph << 1) | 1) & mask ; mh = (mh << 1) & mask ; cross-word carry for << */
        uint64_t new_ph[QDALN_MYERS_MAX_WORDS];
        uint64_t new_mh[QDALN_MYERS_MAX_WORDS];
        uint64_t carry_ph = 1; /* lsb set for ph */
        uint64_t carry_mh = 0;
        for (int w = 0; w < nwords; ++w) {
            uint64_t carry_out_ph = (phh[w] >> 63) & 1ULL;
            uint64_t tmp_ph = (phh[w] << 1) | carry_ph;
            new_ph[w] = tmp_ph & valid_masks[w];
            carry_ph = carry_out_ph;

            uint64_t carry_out_mh = (mhh[w] >> 63) & 1ULL;
            uint64_t tmp_mh = (mhh[w] << 1) | carry_mh;
            new_mh[w] = tmp_mh & valid_masks[w];
            carry_mh = carry_out_mh;
        }
        for (int w = 0; w < nwords; ++w) {
            phh[w] = new_ph[w];
            mhh[w] = new_mh[w];
        }

        for (int w = 0; w < nwords; ++w) {
            uint64_t not_term = (~(xvv[w] | phh[w])) & valid_masks[w];
            pv[w] = (mhh[w] | not_term) & valid_masks[w];
        }
        for (int w = 0; w < nwords; ++w) {
            mv[w] = (phh[w] & xvv[w]) & valid_masks[w];
        }
    }

    return score;
}

int qdaln_edit_distance(const char *a, size_t a_len, const char *b, size_t b_len) {
    if ((a == NULL && a_len != 0) || (b == NULL && b_len != 0)) return -1;

    if (a_len == b_len) {
        int d = same_length_hamming_distance_within_k(a, b, a_len, 2);
        if (d >= 0) return d;
    }

    /* Levenshtein distance is symmetric. Myers now handles arbitrary lengths
       (multi-word), so use it for all cases for speed vs DP. Pick shorter as
       "pattern" to bound the word count in the inner loop (though asymptotic
       cost is identical either way). */
    if (a_len <= b_len) return qdaln_edit_distance_myers64(a, a_len, b, b_len);
    return qdaln_edit_distance_myers64(b, b_len, a, a_len);
}

static int one_delete_matches_qd(const char *longer, size_t longer_len,
                                 const char *shorter, size_t shorter_len) {
    if (longer_len != shorter_len + 1) return 0;
    size_t i = 0;
    size_t j = 0;
    int edits = 0;
    while (i < longer_len && j < shorter_len) {
        if (longer[i] == shorter[j]) {
            ++i;
            ++j;
        } else {
            ++edits;
            if (edits > 1) return 0;
            ++i;
        }
    }
    return 1;
}

int qdaln_edit_distance_leq(const char *a, size_t a_len, const char *b, size_t b_len, int k) {
    if ((a == NULL && a_len != 0) || (b == NULL && b_len != 0)) return -1;
    if (a_len > (size_t)INT32_MAX || b_len > (size_t)INT32_MAX) return -1;
    if (k < 0) return 0;
    if (a_len > b_len) {
        if (a_len - b_len > (size_t)k) return 0;
    } else if (b_len - a_len > (size_t)k) {
        return 0;
    }
    if (k == 0) return a_len == b_len && memcmp(a, b, a_len) == 0 ? 1 : 0;
    if (a_len == 0) return b_len <= (size_t)k ? 1 : 0;
    if (b_len == 0) return a_len <= (size_t)k ? 1 : 0;
    if ((size_t)k >= (a_len > b_len ? a_len : b_len)) return 1;

    if (k == 1) {
        if (a_len == b_len) {
            return same_length_hamming_distance_within_k(a, b, a_len, 1) >= 0 ? 1 : 0;
        }
        if (a_len == b_len + 1) return one_delete_matches_qd(a, a_len, b, b_len);
        if (b_len == a_len + 1) return one_delete_matches_qd(b, b_len, a, a_len);
        return 0;
    }

    if (a_len == b_len && same_length_hamming_distance_within_k(a, b, a_len, k) >= 0) return 1;

    if (k >= 2) {
        /* Myers now generalized to any length via multi-word; use for leq too
           (avoids the DP alloc path and is faster). */
        if (a_len <= b_len) return qdaln_edit_distance_myers64(a, a_len, b, b_len) <= k ? 1 : 0;
        return qdaln_edit_distance_myers64(b, b_len, a, a_len) <= k ? 1 : 0;
    }

    int *prev = (int *)malloc((b_len + 1) * sizeof(int));
    int *curr = (int *)malloc((b_len + 1) * sizeof(int));
    if (prev == NULL || curr == NULL) {
        free(prev);
        free(curr);
        return -1;
    }

    const int inf = k + 1;
    for (size_t j = 0; j <= b_len; ++j) prev[j] = (j <= (size_t)k) ? (int)j : inf;

    for (size_t i = 1; i <= a_len; ++i) {
        size_t start = i > (size_t)k ? i - (size_t)k : 1;
        size_t end = i + (size_t)k < b_len ? i + (size_t)k : b_len;
        int row_min = inf;

        curr[0] = i <= (size_t)k ? (int)i : inf;
        if (start > 1) curr[start - 1] = inf;

        unsigned char ca = (unsigned char)a[i - 1];
        for (size_t j = start; j <= end; ++j) {
            int cost = ca == (unsigned char)b[j - 1] ? 0 : 1;
            int del = prev[j] + 1;
            int ins = curr[j - 1] + 1;
            int sub = prev[j - 1] + cost;
            curr[j] = min3_int(del, ins, sub);
            if (curr[j] < row_min) row_min = curr[j];
        }

        if (end < b_len) curr[end + 1] = inf;
        if (row_min > k) {
            free(prev);
            free(curr);
            return 0;
        }

        int *tmp = prev;
        prev = curr;
        curr = tmp;
    }

    int result = prev[b_len] <= k ? 1 : 0;
    free(prev);
    free(curr);
    return result;
}

static qdaln_match_result empty_match_result(int status) {
    qdaln_match_result r;
    r.target_index = -1;
    r.best_distance = -1;
    r.second_best_distance = -1;
    r.match_count = 0;
    r.status = status;
    return r;
}

static qdaln_assignment_result empty_assignment_result(int status) {
    qdaln_assignment_result r;
    r.target_index = -1;
    r.distance = -1;
    r.second_best_distance = -1;
    r.num_best_targets = 0;
    r.num_targets_within_radius = 0;
    r.status = status;
    r.edit_class = status == QDALN_MATCH_NONE ? QDALN_EDIT_NONE : QDALN_EDIT_INVALID;
    return r;
}

static int assignment_edit_class(const char *read, size_t read_len,
                                 const char *target, size_t target_len, int distance) {
    if (distance < 0) return QDALN_EDIT_INVALID;
    if (distance == 0) return QDALN_EDIT_EXACT;
    if (distance == 1) {
        if (read_len == target_len) return QDALN_EDIT_K1_SUB;
        if (one_delete_matches_qd(read, read_len, target, target_len)) return QDALN_EDIT_K1_INS;
        if (one_delete_matches_qd(target, target_len, read, read_len)) return QDALN_EDIT_K1_DEL;
        return QDALN_EDIT_OTHER;
    }
    if (distance == 2) return QDALN_EDIT_K2;
    return QDALN_EDIT_OTHER;
}

static int candidate_distance_within_k(const char *read, size_t read_len,
                                       const char *target, size_t target_len, int k) {
    if (k == 0) {
        return read_len == target_len && memcmp(read, target, read_len) == 0 ? 0 : -1;
    }

    if (k == 1) {
        if ((read == NULL && read_len != 0) || (target == NULL && target_len != 0)) return -2;
        if (read_len == target_len) {
            return same_length_hamming_distance_within_k(read, target, read_len, 1);
        }
        if (read_len == target_len + 1) {
            return one_delete_matches_qd(read, read_len, target, target_len) ? 1 : -1;
        }
        if (target_len == read_len + 1) {
            return one_delete_matches_qd(target, target_len, read, read_len) ? 1 : -1;
        }
        return -1;
    }

    if (read_len == target_len) {
        int d = same_length_hamming_distance_within_k(read, target, read_len, k);
        if (d >= 0 && d <= 2) return d;
    }

    if (k >= 2) {
        /* Use generalized myers for any length when k>=2 (no more 64 limit).
           This is exact and fast; falls to leq/edit only for k<2 (already handled). */
        int d = read_len <= target_len
                    ? qdaln_edit_distance_myers64(read, read_len, target, target_len)
                    : qdaln_edit_distance_myers64(target, target_len, read, read_len);
        return d <= k ? d : -1;
    }

    int within = qdaln_edit_distance_leq(read, read_len, target, target_len, k);
    if (within <= 0) return within < 0 ? -2 : -1;

    int d = qdaln_edit_distance(read, read_len, target, target_len);
    if (d < 0) return -2;
    return d <= k ? d : -1;
}

int qdaln_match_many(const char *const *reads, const size_t *read_lens, size_t n_reads,
                     const char *const *targets, const size_t *target_lens, size_t n_targets,
                     int k, qdaln_match_result *results) {
    if (results == NULL) return -1;
    if (k < 0) return -1;
    if (n_reads != 0 && (reads == NULL || read_lens == NULL)) return -1;
    if (n_targets != 0 && (targets == NULL || target_lens == NULL)) return -1;

    for (size_t i = 0; i < n_reads; ++i) {
        results[i] = empty_match_result(QDALN_MATCH_NONE);

        const char *read = reads[i];
        size_t read_len = read_lens[i];
        if (read == NULL && read_len != 0) {
            results[i] = empty_match_result(QDALN_MATCH_INVALID);
            continue;
        }

        int best_tie_count = 0;
        for (size_t j = 0; j < n_targets; ++j) {
            const char *target = targets[j];
            size_t target_len = target_lens[j];
            if (target == NULL && target_len != 0) {
                results[i] = empty_match_result(QDALN_MATCH_INVALID);
                best_tie_count = 0;
                break;
            }

            int d = candidate_distance_within_k(read, read_len, target, target_len, k);
            if (d == -2) {
                results[i] = empty_match_result(QDALN_MATCH_INVALID);
                best_tie_count = 0;
                break;
            }
            if (d < 0) continue;

            ++results[i].match_count;
            if (results[i].best_distance < 0 || d < results[i].best_distance) {
                results[i].second_best_distance = results[i].best_distance;
                results[i].best_distance = d;
                results[i].target_index = (int)j;
                best_tie_count = 1;
            } else if (d == results[i].best_distance) {
                ++best_tie_count;
            } else if (results[i].second_best_distance < 0 || d < results[i].second_best_distance) {
                results[i].second_best_distance = d;
            }
        }

        if (results[i].status == QDALN_MATCH_INVALID) continue;
        if (results[i].match_count == 0) {
            results[i].status = QDALN_MATCH_NONE;
        } else if (best_tie_count > 1) {
            results[i].status = QDALN_MATCH_AMBIGUOUS;
        } else {
            results[i].status = QDALN_MATCH_UNIQUE;
        }
    }

    return 0;
}

int qdaln_assign_many(const char *const *reads, const size_t *read_lens, size_t n_reads,
                      const char *const *targets, const size_t *target_lens, size_t n_targets,
                      int k, int policy, qdaln_assignment_result *results) {
    if (results == NULL) return -1;
    if (k < 0) return -1;
    if (policy != QDALN_POLICY_BEST && policy != QDALN_POLICY_RADIUS) return -1;
    if (n_reads != 0 && (reads == NULL || read_lens == NULL)) return -1;
    if (n_targets != 0 && (targets == NULL || target_lens == NULL)) return -1;

    for (size_t i = 0; i < n_reads; ++i) {
        results[i] = empty_assignment_result(QDALN_MATCH_NONE);
        const char *read = reads[i];
        size_t read_len = read_lens[i];
        if (read == NULL && read_len != 0) {
            results[i] = empty_assignment_result(QDALN_MATCH_INVALID);
            continue;
        }

        for (size_t j = 0; j < n_targets; ++j) {
            const char *target = targets[j];
            size_t target_len = target_lens[j];
            if (target == NULL && target_len != 0) {
                results[i] = empty_assignment_result(QDALN_MATCH_INVALID);
                break;
            }

            int d = candidate_distance_within_k(read, read_len, target, target_len, k);
            if (d == -2) {
                results[i] = empty_assignment_result(QDALN_MATCH_INVALID);
                break;
            }
            if (d < 0) continue;

            ++results[i].num_targets_within_radius;
            if (results[i].distance < 0 || d < results[i].distance) {
                results[i].second_best_distance = results[i].distance;
                results[i].distance = d;
                results[i].target_index = (int)j;
                results[i].num_best_targets = 1;
                results[i].edit_class = assignment_edit_class(read, read_len, target, target_len, d);
            } else if (d == results[i].distance) {
                ++results[i].num_best_targets;
                if (results[i].target_index < 0 || (int)j < results[i].target_index) {
                    results[i].target_index = (int)j;
                    results[i].edit_class = assignment_edit_class(read, read_len, target, target_len, d);
                }
            } else if (results[i].second_best_distance < 0 || d < results[i].second_best_distance) {
                results[i].second_best_distance = d;
            }
        }

        if (results[i].status == QDALN_MATCH_INVALID) continue;
        if (results[i].num_targets_within_radius == 0) {
            results[i].status = QDALN_MATCH_NONE;
            results[i].edit_class = QDALN_EDIT_NONE;
        } else if (policy == QDALN_POLICY_RADIUS && results[i].num_targets_within_radius > 1) {
            results[i].status = QDALN_MATCH_AMBIGUOUS;
        } else if (results[i].num_best_targets > 1) {
            results[i].status = QDALN_MATCH_AMBIGUOUS;
        } else {
            results[i].status = QDALN_MATCH_UNIQUE;
        }
    }
    return 0;
}

static int index_hamming_seed_insert(qdaln_index *index, unsigned char seed_id, uint64_t code,
                                     size_t seed_len, size_t target_index) {
    size_t slot = hamming_seed_hash(code, seed_len, seed_id, index->hamming_seed_hash_cap);
    size_t e = index->n_hamming_seeds++;
    index->hamming_seeds[e].code = code;
    index->hamming_seeds[e].target_index = target_index;
    index->hamming_seeds[e].seed_len = seed_len;
    index->hamming_seeds[e].seed_id = seed_id;
    index->hamming_seeds[e].next = index->hamming_seed_heads[slot];
    index->hamming_seed_heads[slot] = e;
    return 0;
}

static int index_exact_table_insert(qdaln_index *index, uint64_t code, size_t len, size_t target_index) {
    if (index->exact_hash_cap == 0) return 0;
    size_t slot = exact_code_hash(code, len, index->exact_hash_cap);
    for (;;) {
        size_t first = index->exact_first_targets[slot];
        if (first == QDALN_HASH_EMPTY) {
            index->exact_codes[slot] = code;
            index->exact_lens[slot] = len;
            index->exact_first_targets[slot] = target_index;
            index->exact_counts[slot] = 1;
            return 0;
        }
        if (index->exact_codes[slot] == code && index->exact_lens[slot] == len) {
            if (target_index < index->exact_first_targets[slot]) index->exact_first_targets[slot] = target_index;
            ++index->exact_counts[slot];
            return 0;
        }
        slot = (slot + 1U) & (index->exact_hash_cap - 1U);
    }
}

static int index_exact_table_find(const qdaln_index *index, uint64_t code, size_t len,
                                  size_t *first_target, size_t *count) {
    if (index->exact_hash_cap == 0) return 0;
    size_t slot = exact_code_hash(code, len, index->exact_hash_cap);
    for (;;) {
        size_t first = index->exact_first_targets[slot];
        if (first == QDALN_HASH_EMPTY) return 0;
        if (index->exact_codes[slot] == code && index->exact_lens[slot] == len) {
            *first_target = first;
            *count = index->exact_counts[slot];
            return 1;
        }
        slot = (slot + 1U) & (index->exact_hash_cap - 1U);
    }
}

qdaln_index *qdaln_index_build(const char *const *targets, const size_t *target_lens, size_t n_targets) {
    if (n_targets != 0 && (targets == NULL || target_lens == NULL)) return NULL;

    qdaln_index *idx = (qdaln_index *)calloc(1, sizeof(qdaln_index));
    if (idx == NULL) return NULL;
    idx->n_targets = n_targets;

    if (n_targets == 0) return idx;

    idx->targets = (char **)calloc(n_targets, sizeof(char *));
    idx->target_lens = (size_t *)calloc(n_targets, sizeof(size_t));
    idx->codes = (uint64_t *)calloc(n_targets, sizeof(uint64_t));
    idx->encodable = (unsigned char *)calloc(n_targets, sizeof(unsigned char));
    idx->hash_cap = next_pow2_size(n_targets * 2 + 1);
    idx->hash_heads = (size_t *)malloc(idx->hash_cap * sizeof(size_t));
    idx->hash_next = (size_t *)malloc(n_targets * sizeof(size_t));
    idx->nonencodable_targets = (size_t *)malloc(n_targets * sizeof(size_t));
    if (idx->targets == NULL || idx->target_lens == NULL || idx->codes == NULL ||
        idx->encodable == NULL || idx->hash_heads == NULL || idx->hash_next == NULL ||
        idx->nonencodable_targets == NULL) {
        qdaln_index_free(idx);
        return NULL;
    }
    for (size_t i = 0; i < idx->hash_cap; ++i) {
        idx->hash_heads[i] = QDALN_HASH_EMPTY;
    }
    for (size_t i = 0; i < n_targets; ++i) idx->hash_next[i] = QDALN_HASH_EMPTY;

    idx->exact_hash_cap = next_pow2_size(n_targets * 2 + 1);
    idx->exact_codes = (uint64_t *)calloc(idx->exact_hash_cap, sizeof(uint64_t));
    idx->exact_lens = (size_t *)calloc(idx->exact_hash_cap, sizeof(size_t));
    idx->exact_first_targets = (size_t *)malloc(idx->exact_hash_cap * sizeof(size_t));
    idx->exact_counts = (size_t *)calloc(idx->exact_hash_cap, sizeof(size_t));
    if (idx->exact_codes == NULL || idx->exact_lens == NULL || idx->exact_first_targets == NULL ||
        idx->exact_counts == NULL) {
        qdaln_index_free(idx);
        return NULL;
    }
    for (size_t i = 0; i < idx->exact_hash_cap; ++i) idx->exact_first_targets[i] = QDALN_HASH_EMPTY;

    for (size_t i = 0; i < n_targets; ++i) {
        if (targets[i] == NULL && target_lens[i] != 0) {
            qdaln_index_free(idx);
            return NULL;
        }
        idx->target_lens[i] = target_lens[i];
        idx->targets[i] = (char *)malloc(target_lens[i] + 1);
        if (idx->targets[i] == NULL) {
            qdaln_index_free(idx);
            return NULL;
        }
        if (target_lens[i] != 0) memcpy(idx->targets[i], targets[i], target_lens[i]);
        idx->targets[i][target_lens[i]] = '\0';
        idx->encodable[i] = (unsigned char)dna2_code(idx->targets[i], target_lens[i], &idx->codes[i]);
        if (idx->encodable[i]) {
            size_t slot = code_hash(idx->codes[i], target_lens[i], idx->hash_cap);
            idx->hash_next[i] = idx->hash_heads[slot];
            idx->hash_heads[slot] = i;
            if (index_exact_table_insert(idx, idx->codes[i], target_lens[i], i) != 0) {
                qdaln_index_free(idx);
                return NULL;
            }
        } else {
            idx->nonencodable_targets[idx->n_nonencodable] = i;
            ++idx->n_nonencodable;
        }
    }

    size_t deletion_need = 0;
    for (size_t i = 0; i < n_targets; ++i) {
        if (idx->encodable[i] && idx->target_lens[i] > 0 && idx->target_lens[i] <= 32) {
            deletion_need += idx->target_lens[i];
        }
    }
    if (deletion_need != 0) {
        idx->deletion_hash_cap = next_pow2_size(deletion_need * 2 + 1);
        idx->deletion_codes = (uint64_t *)malloc(deletion_need * sizeof(uint64_t));
        idx->deletion_lens = (size_t *)malloc(deletion_need * sizeof(size_t));
        idx->deletion_targets = (size_t *)malloc(deletion_need * sizeof(size_t));
        idx->deletion_next = (size_t *)malloc(deletion_need * sizeof(size_t));
        idx->deletion_heads = (size_t *)malloc(idx->deletion_hash_cap * sizeof(size_t));
        if (idx->deletion_codes == NULL || idx->deletion_lens == NULL || idx->deletion_targets == NULL ||
            idx->deletion_next == NULL || idx->deletion_heads == NULL) {
            qdaln_index_free(idx);
            return NULL;
        }
        for (size_t i = 0; i < idx->deletion_hash_cap; ++i) idx->deletion_heads[i] = SIZE_MAX;

        for (size_t i = 0; i < n_targets; ++i) {
            if (!idx->encodable[i] || idx->target_lens[i] == 0 || idx->target_lens[i] > 32) continue;
            size_t del_len = idx->target_lens[i] - 1;
            for (size_t pos = 0; pos < idx->target_lens[i]; ++pos) {
                size_t e = idx->n_deletions++;
                uint64_t del_code = code_delete_base(idx->codes[i], idx->target_lens[i], pos);
                idx->deletion_codes[e] = del_code;
                idx->deletion_lens[e] = del_len;
                idx->deletion_targets[e] = i;
                size_t slot = code_hash(del_code, del_len, idx->deletion_hash_cap);
                idx->deletion_next[e] = idx->deletion_heads[slot];
                idx->deletion_heads[slot] = e;
            }
        }
    }

    size_t hamming_seed_need = 0;
    for (size_t i = 0; i < n_targets; ++i) {
        if (idx->encodable[i] && idx->target_lens[i] <= 32) {
            hamming_seed_need += 2 + 3 + 4;
        }
    }
    if (hamming_seed_need != 0) {
        idx->hamming_seed_hash_cap = next_pow2_size(hamming_seed_need * 2 + 1);
        idx->hamming_seeds = (qdaln_hamming_seed_entry *)malloc(hamming_seed_need * sizeof(qdaln_hamming_seed_entry));
        idx->hamming_seed_heads = (size_t *)malloc(idx->hamming_seed_hash_cap * sizeof(size_t));
        if (idx->hamming_seeds == NULL || idx->hamming_seed_heads == NULL) {
            qdaln_index_free(idx);
            return NULL;
        }
        for (size_t i = 0; i < idx->hamming_seed_hash_cap; ++i) idx->hamming_seed_heads[i] = SIZE_MAX;

        for (size_t i = 0; i < n_targets; ++i) {
            if (!idx->encodable[i] || idx->target_lens[i] > 32) continue;
            for (size_t n_seeds = 2; n_seeds <= 4; ++n_seeds) {
                for (size_t seed_id = 0; seed_id < n_seeds; ++seed_id) {
                    size_t start = 0;
                    size_t seed_len = 0;
                    hamming_seed_partition(idx->target_lens[i], n_seeds, seed_id, &start, &seed_len);
                    uint64_t seed = code_segment_qd(idx->codes[i], start, seed_len);
                    if (index_hamming_seed_insert(idx, hamming_seed_key(n_seeds, seed_id), seed, seed_len, i) != 0) {
                        qdaln_index_free(idx);
                        return NULL;
                    }
                }
            }
        }
        idx->hamming_seed_ready = 1;
    }

    return idx;
}

void qdaln_index_free(qdaln_index *index) {
    if (index == NULL) return;
    if (index->targets != NULL) {
        for (size_t i = 0; i < index->n_targets; ++i) free(index->targets[i]);
    }
    free(index->targets);
    free(index->target_lens);
    free(index->codes);
    free(index->encodable);
    free(index->hash_heads);
    free(index->hash_next);
    free(index->exact_codes);
    free(index->exact_lens);
    free(index->exact_first_targets);
    free(index->exact_counts);
    free(index->deletion_codes);
    free(index->deletion_lens);
    free(index->deletion_targets);
    free(index->deletion_heads);
    free(index->deletion_next);
    free(index->hamming_seeds);
    free(index->hamming_seed_heads);
    free(index->nonencodable_targets);
    free(index);
}

static void index_consider_candidate(qdaln_match_result *r, int target_index, int distance) {
    ++r->match_count;
    if (r->best_distance < 0 || distance < r->best_distance) {
        r->second_best_distance = r->best_distance;
        r->best_distance = distance;
        r->target_index = target_index;
    } else if (distance > r->best_distance &&
               (r->second_best_distance < 0 || distance < r->second_best_distance)) {
        r->second_best_distance = distance;
    }
}

static void index_finalize_result(qdaln_match_result *r, int best_tie_count) {
    if (r->match_count == 0) {
        r->status = QDALN_MATCH_NONE;
    } else if (best_tie_count > 1) {
        r->status = QDALN_MATCH_AMBIGUOUS;
    } else {
        r->status = QDALN_MATCH_UNIQUE;
    }
}

static void index_update_verified(qdaln_match_result *result, int target_index, int distance, int *best_tie_count) {
    if (result->best_distance < 0 || distance < result->best_distance) {
        index_consider_candidate(result, target_index, distance);
        *best_tie_count = 1;
    } else if (distance == result->best_distance) {
        index_consider_candidate(result, target_index, distance);
        if (result->target_index < 0 || target_index < result->target_index) result->target_index = target_index;
        ++*best_tie_count;
    } else {
        index_consider_candidate(result, target_index, distance);
    }
}

static int index_visit_hamming_code_candidates(const qdaln_index *index, uint64_t code, size_t len,
                                               unsigned char *seen, int distance,
                                               qdaln_match_result *result, int *best_tie_count,
                                               qdaln_index_stats *stats) {
    size_t slot = code_hash(code, len, index->hash_cap);
    for (size_t j = index->hash_heads[slot]; j != SIZE_MAX; ) {
        size_t nextj = index->hash_next[j];
        if (!index->encodable[j] || index->target_lens[j] != len || index->codes[j] != code) {
            j = nextj;
            continue;
        }
        if (seen != NULL) {
            if (seen[j]) {
                j = nextj;
                continue;
            }
            seen[j] = 1;
        }
        if (stats != NULL) {
            ++stats->candidates_considered;
            ++stats->candidates_verified;
        }
        index_update_verified(result, (int)j, distance, best_tie_count);
        j = nextj;
    }
    return 0;
}

typedef struct candidate_seen candidate_seen;
static int candidate_seen_add(candidate_seen *seen, size_t target_index);

static int index_verify_seed_candidate(const qdaln_index *index, size_t target_index, candidate_seen *seen,
                                       const char *read, size_t read_len, int k,
                                       qdaln_match_result *result, int *best_tie_count,
                                       qdaln_index_stats *stats) {
    int seen_rc = candidate_seen_add(seen, target_index);
    if (seen_rc <= 0) return seen_rc;
    if (stats != NULL) ++stats->candidates_considered;
    int d = candidate_distance_within_k(read, read_len, index->targets[target_index], index->target_lens[target_index], k);
    if (stats != NULL) ++stats->candidates_verified;
    if (d < 0) return 0;
    index_update_verified(result, (int)target_index, d, best_tie_count);
    return 0;
}

static int index_visit_exact_seed_candidates(const qdaln_index *index, uint64_t code, size_t len,
                                             candidate_seen *seen, const char *read, size_t read_len, int k,
                                             qdaln_match_result *result, int *best_tie_count,
                                             qdaln_index_stats *stats, int stop_on_ambiguous) {
    size_t slot = code_hash(code, len, index->hash_cap);
    for (size_t j = index->hash_heads[slot]; j != SIZE_MAX; j = index->hash_next[j]) {
        if (!index->encodable[j] || index->target_lens[j] != len || index->codes[j] != code) continue;
        if (index_verify_seed_candidate(index, j, seen, read, read_len, k, result, best_tie_count, stats) != 0) {
            return -1;
        }
        if (stop_on_ambiguous && *best_tie_count > 1) return 0;
    }
    return 0;
}

static int index_visit_deletion_seed_candidates(const qdaln_index *index, uint64_t code, size_t len,
                                                candidate_seen *seen, const char *read, size_t read_len, int k,
                                                qdaln_match_result *result, int *best_tie_count,
                                                qdaln_index_stats *stats, int stop_on_ambiguous) {
    if (index->deletion_hash_cap == 0) return 0;
    size_t slot = code_hash(code, len, index->deletion_hash_cap);
    for (size_t e = index->deletion_heads[slot]; e != SIZE_MAX; e = index->deletion_next[e]) {
        if (index->deletion_lens[e] != len || index->deletion_codes[e] != code) continue;
        size_t j = index->deletion_targets[e];
        if (index_verify_seed_candidate(index, j, seen, read, read_len, k, result, best_tie_count, stats) != 0) {
            return -1;
        }
        if (stop_on_ambiguous && *best_tie_count > 1) return 0;
    }
    return 0;
}

static int code_already_seen(const uint64_t *codes, size_t n_codes, uint64_t code) {
    for (size_t i = 0; i < n_codes; ++i) {
        if (codes[i] == code) return 1;
    }
    return 0;
}

struct candidate_seen {
    size_t inline_ids[128];
    size_t *ids;
    size_t *hash;
    size_t count;
    size_t cap;
    size_t hash_cap;
};

static void candidate_seen_init(candidate_seen *seen) {
    seen->ids = seen->inline_ids;
    seen->hash = NULL;
    seen->count = 0;
    seen->cap = sizeof(seen->inline_ids) / sizeof(seen->inline_ids[0]);
    seen->hash_cap = 0;
}

static void candidate_seen_free(candidate_seen *seen) {
    if (seen->ids != seen->inline_ids) free(seen->ids);
    free(seen->hash);
    candidate_seen_init(seen);
}

static size_t candidate_seen_hash_value(size_t value) {
    uint64_t x = (uint64_t)value;
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return (size_t)x;
}

static int candidate_seen_rehash(candidate_seen *seen, size_t min_cap) {
    size_t cap = 128;
    while (cap < min_cap) cap <<= 1;
    size_t *hash = (size_t *)calloc(cap, sizeof(size_t));
    if (hash == NULL) return -1;

    for (size_t i = 0; i < seen->count; ++i) {
        size_t stored = seen->ids[i] + 1U;
        size_t pos = candidate_seen_hash_value(seen->ids[i]) & (cap - 1U);
        while (hash[pos] != 0) pos = (pos + 1U) & (cap - 1U);
        hash[pos] = stored;
    }

    free(seen->hash);
    seen->hash = hash;
    seen->hash_cap = cap;
    return 0;
}

static int candidate_seen_hash_add(candidate_seen *seen, size_t target_index) {
    if (seen->hash_cap == 0 || ((seen->count + 1U) * 2U) > seen->hash_cap) {
        if (candidate_seen_rehash(seen, (seen->count + 1U) * 4U) != 0) return -1;
    }

    size_t stored = target_index + 1U;
    size_t pos = candidate_seen_hash_value(target_index) & (seen->hash_cap - 1U);
    while (seen->hash[pos] != 0) {
        if (seen->hash[pos] == stored) return 0;
        pos = (pos + 1U) & (seen->hash_cap - 1U);
    }
    seen->hash[pos] = stored;
    return 1;
}

static int candidate_seen_add(candidate_seen *seen, size_t target_index) {
    if (seen->hash != NULL) {
        int hash_rc = candidate_seen_hash_add(seen, target_index);
        if (hash_rc <= 0) return hash_rc;
    } else if (seen->count == seen->cap) {
        if (candidate_seen_rehash(seen, seen->cap * 4U) != 0) return -1;
        int hash_rc = candidate_seen_hash_add(seen, target_index);
        if (hash_rc <= 0) return hash_rc;
    } else {
        for (size_t i = 0; i < seen->count; ++i) {
            if (seen->ids[i] == target_index) return 0;
        }
    }
    if (seen->count == seen->cap) {
        size_t next_cap = seen->cap * 2;
        size_t *next = (size_t *)malloc(next_cap * sizeof(size_t));
        if (next == NULL) return -1;
        memcpy(next, seen->ids, seen->count * sizeof(size_t));
        if (seen->ids != seen->inline_ids) free(seen->ids);
        seen->ids = next;
        seen->cap = next_cap;
    }
    seen->ids[seen->count++] = target_index;
    return 1;
}

static size_t count_non_acgt(const char *read, size_t read_len, size_t *last_bad_pos) {
    size_t n_bad = 0;
    for (size_t i = 0; i < read_len; ++i) {
        switch (read[i]) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                break;
            default:
                if (last_bad_pos != NULL) *last_bad_pos = i;
                ++n_bad;
                break;
        }
    }
    return n_bad;
}

static size_t collect_non_acgt_positions(const char *read, size_t read_len, size_t *positions, size_t max_positions) {
    size_t n_bad = 0;
    for (size_t i = 0; i < read_len; ++i) {
        switch (read[i]) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                break;
            default:
                if (n_bad < max_positions) positions[n_bad] = i;
                ++n_bad;
                break;
        }
    }
    return n_bad;
}

static int find_single_non_acgt(const char *read, size_t read_len, size_t *bad_pos) {
    if ((read == NULL && read_len != 0) || read_len > 32) return 0;
    return count_non_acgt(read, read_len, bad_pos) == 1;
}

static int code_without_base(const char *read, size_t read_len, size_t drop_pos, uint64_t *code_out) {
    if ((read == NULL && read_len != 0) || read_len == 0 || drop_pos >= read_len || read_len - 1 > 32) return 0;
    uint64_t code = 0;
    size_t out = 0;
    for (size_t i = 0; i < read_len; ++i) {
        if (i == drop_pos) continue;
        uint64_t v;
        switch (read[i]) {
            case 'A':
                v = 0;
                break;
            case 'C':
                v = 1;
                break;
            case 'G':
                v = 2;
                break;
            case 'T':
                v = 3;
                break;
            default:
                return 0;
        }
        code |= v << (2 * out);
        ++out;
    }
    *code_out = code;
    return 1;
}

static int code_without_two_bases(const char *read, size_t read_len, size_t drop_a, size_t drop_b,
                                  uint64_t *code_out) {
    if ((read == NULL && read_len != 0) || read_len < 2 || drop_a >= read_len || drop_b >= read_len ||
        drop_a == drop_b || read_len - 2 > 32) {
        return 0;
    }
    if (drop_a > drop_b) {
        size_t tmp = drop_a;
        drop_a = drop_b;
        drop_b = tmp;
    }
    uint64_t code = 0;
    size_t out = 0;
    for (size_t i = 0; i < read_len; ++i) {
        if (i == drop_a || i == drop_b) continue;
        uint64_t v;
        switch (read[i]) {
            case 'A':
                v = 0;
                break;
            case 'C':
                v = 1;
                break;
            case 'G':
                v = 2;
                break;
            case 'T':
                v = 3;
                break;
            default:
                return 0;
        }
        code |= v << (2 * out);
        ++out;
    }
    *code_out = code;
    return 1;
}

static int index_assign_single_unknown_one(const qdaln_index *index, const char *read, size_t read_len,
                                           qdaln_match_result *result, qdaln_index_stats *stats,
                                           int stop_on_ambiguous) {
    size_t bad_pos = 0;
    if ((read == NULL && read_len != 0) || read_len > 32) return 0;
    size_t n_bad = count_non_acgt(read, read_len, &bad_pos);
    if (n_bad == 0) return 0;

    static const char bases[] = {'A', 'C', 'G', 'T'};
    char tmp[32];
    int best_tie_count = 0;
    *result = empty_match_result(QDALN_MATCH_NONE);
    candidate_seen seen;
    candidate_seen_init(&seen);

    if (n_bad == 1) {
        memcpy(tmp, read, read_len);
        for (size_t b = 0; b < 4; ++b) {
            uint64_t code = 0;
            tmp[bad_pos] = bases[b];
            if (dna2_code(tmp, read_len, &code)) {
                if (index_visit_exact_seed_candidates(index, code, read_len, &seen, read, read_len, 1,
                                                      result, &best_tie_count, stats, stop_on_ambiguous) != 0) {
                    candidate_seen_free(&seen);
                    return 0;
                }
                if (stop_on_ambiguous && best_tie_count > 1) goto done;
            }
        }
        uint64_t code = 0;
        if (code_without_base(read, read_len, bad_pos, &code)) {
            if (index_visit_exact_seed_candidates(index, code, read_len - 1, &seen, read, read_len, 1,
                                                  result, &best_tie_count, stats, stop_on_ambiguous) != 0) {
                candidate_seen_free(&seen);
                return 0;
            }
            if (stop_on_ambiguous && best_tie_count > 1) goto done;
        }
    }

    if (index->n_nonencodable != 0) {
        for (size_t ni = 0; ni < index->n_nonencodable; ++ni) {
            size_t j = index->nonencodable_targets[ni];
            size_t target_len = index->target_lens[j];
            if (target_len > read_len + 1 || read_len > target_len + 1) continue;
            int seen_rc = candidate_seen_add(&seen, j);
            if (seen_rc < 0) {
                candidate_seen_free(&seen);
                return 0;
            }
            if (seen_rc == 0) continue;
            if (stats != NULL) ++stats->candidates_considered;
            int d = candidate_distance_within_k(read, read_len, index->targets[j], target_len, 1);
            if (stats != NULL) ++stats->candidates_verified;
            if (d < 0) continue;
            index_update_verified(result, (int)j, d, &best_tie_count);
            if (stop_on_ambiguous && best_tie_count > 1) goto done;
        }
    }

done:
    index_finalize_result(result, best_tie_count);
    candidate_seen_free(&seen);
    return 1;
}

static int index_assign_neighbor_one_impl(const qdaln_index *index, const char *read, size_t read_len,
                                          uint64_t read_code, qdaln_match_result *result,
                                          qdaln_index_stats *stats, int stop_on_ambiguous) {
    if (index->n_targets == 0) {
        *result = empty_match_result(QDALN_MATCH_NONE);
        return 1;
    }

    candidate_seen seen;
    candidate_seen_init(&seen);
    int best_tie_count = 0;
    *result = empty_match_result(QDALN_MATCH_NONE);

    if (index_visit_exact_seed_candidates(index, read_code, read_len, &seen, read, read_len, 1, result,
                                          &best_tie_count, stats, stop_on_ambiguous) != 0) {
        candidate_seen_free(&seen);
        return 0;
    }
    if (stop_on_ambiguous && result->best_distance == 0) goto done;

    uint64_t seen_sub_codes[128];
    size_t n_seen_sub_codes = 0;
    seen_sub_codes[n_seen_sub_codes++] = read_code;
    for (size_t pos = 0; pos < read_len; ++pos) {
        uint64_t shift = (uint64_t)2 * pos;
        uint64_t old_base = (read_code >> shift) & 3ULL;
        uint64_t mask = 3ULL << shift;
        for (uint64_t b = 0; b < 4; ++b) {
            if (b == old_base) continue;
            uint64_t code = (read_code & ~mask) | (b << shift);
            if (code_already_seen(seen_sub_codes, n_seen_sub_codes, code)) continue;
            if (n_seen_sub_codes < sizeof(seen_sub_codes) / sizeof(seen_sub_codes[0])) {
                seen_sub_codes[n_seen_sub_codes++] = code;
            }
            if (index_visit_exact_seed_candidates(index, code, read_len, &seen, read, read_len, 1, result,
                                                  &best_tie_count, stats, stop_on_ambiguous) != 0) {
                candidate_seen_free(&seen);
                return 0;
            }
            if (stop_on_ambiguous && best_tie_count > 1) goto done;
        }
    }

    if (index_visit_deletion_seed_candidates(index, read_code, read_len, &seen, read, read_len, 1, result,
                                             &best_tie_count, stats, stop_on_ambiguous) != 0) {
        candidate_seen_free(&seen);
        return 0;
    }
    if (stop_on_ambiguous && best_tie_count > 1) goto done;

    uint64_t seen_del_codes[64];
    size_t n_seen_del_codes = 0;
    for (size_t pos = 0; pos < read_len; ++pos) {
        uint64_t code = code_delete_base(read_code, read_len, pos);
        if (code_already_seen(seen_del_codes, n_seen_del_codes, code)) continue;
        if (n_seen_del_codes < sizeof(seen_del_codes) / sizeof(seen_del_codes[0])) {
            seen_del_codes[n_seen_del_codes++] = code;
        }
        if (index_visit_exact_seed_candidates(index, code, read_len - 1, &seen, read, read_len, 1, result,
                                              &best_tie_count, stats, stop_on_ambiguous) != 0) {
            candidate_seen_free(&seen);
            return 0;
        }
        if (stop_on_ambiguous && best_tie_count > 1) goto done;
    }

    for (size_t ni = 0; ni < index->n_nonencodable; ++ni) {
        size_t j = index->nonencodable_targets[ni];
        int seen_rc = candidate_seen_add(&seen, j);
        if (seen_rc < 0) {
            candidate_seen_free(&seen);
            return 0;
        }
        if (seen_rc == 0) continue;
        size_t target_len = index->target_lens[j];
        if (target_len > read_len + 1 || read_len > target_len + 1) continue;
        if (stats != NULL) ++stats->candidates_considered;
        int d = candidate_distance_within_k(read, read_len, index->targets[j], target_len, 1);
        if (stats != NULL) ++stats->candidates_verified;
        if (d < 0) continue;
        index_update_verified(result, (int)j, d, &best_tie_count);
        if (stop_on_ambiguous && best_tie_count > 1) goto done;
    }

done:
    index_finalize_result(result, best_tie_count);
    candidate_seen_free(&seen);
    return 1;
}

static int index_assign_neighbor_one(const qdaln_index *index, const char *read, size_t read_len,
                                     uint64_t read_code, qdaln_match_result *result,
                                     qdaln_index_stats *stats) {
    return index_assign_neighbor_one_impl(index, read, read_len, read_code, result, stats, 0);
}

typedef int (*edit_code_visitor)(uint64_t code, size_t len, void *ctx);

static int index_visit_one_edit_codes(uint64_t code, size_t len, edit_code_visitor visit, void *ctx) {
    for (size_t pos = 0; pos < len; ++pos) {
        uint64_t shift = (uint64_t)2 * pos;
        uint64_t old_base = (code >> shift) & 3ULL;
        uint64_t mask = 3ULL << shift;
        for (uint64_t b = 0; b < 4; ++b) {
            if (b == old_base) continue;
            int rc = visit((code & ~mask) | (b << shift), len, ctx);
            if (rc != 0) return rc;
        }
    }

    for (size_t pos = 0; pos < len; ++pos) {
        int rc = visit(code_delete_base(code, len, pos), len - 1, ctx);
        if (rc != 0) return rc;
    }

    if (len < 32) {
        for (size_t pos = 0; pos <= len; ++pos) {
            for (uint64_t b = 0; b < 4; ++b) {
                int rc = visit(code_insert_base(code, len, pos, b), len + 1, ctx);
                if (rc != 0) return rc;
            }
        }
    }

    return 0;
}

typedef struct k2_exact_visit_ctx {
    const qdaln_index *index;
    candidate_seen *seen;
    const char *read;
    size_t read_len;
    qdaln_match_result *result;
    int *best_tie_count;
    qdaln_index_stats *stats;
    int stop_on_ambiguous;
} k2_exact_visit_ctx;

static int visit_k2_exact_code(uint64_t code, size_t len, void *ctx) {
    k2_exact_visit_ctx *v = (k2_exact_visit_ctx *)ctx;
    if (index_visit_exact_seed_candidates(v->index, code, len, v->seen, v->read, v->read_len, 2,
                                          v->result, v->best_tie_count, v->stats,
                                          v->stop_on_ambiguous) != 0) {
        return -1;
    }
    return v->stop_on_ambiguous && *v->best_tie_count > 1 ? 1 : 0;
}

static int visit_k2_first_edit_code(uint64_t code, size_t len, void *ctx) {
    int rc = visit_k2_exact_code(code, len, ctx);
    if (rc != 0) return rc;
    return index_visit_one_edit_codes(code, len, visit_k2_exact_code, ctx);
}

static int lengths_within_radius(size_t a_len, size_t b_len, int k) {
    return a_len > b_len ? a_len - b_len <= (size_t)k : b_len - a_len <= (size_t)k;
}

static int index_scan_nonencodable_seed_candidates(const qdaln_index *index, candidate_seen *seen,
                                                   const char *read, size_t read_len, int k,
                                                   qdaln_match_result *result, int *best_tie_count,
                                                   qdaln_index_stats *stats, int stop_on_ambiguous) {
    if (index->n_nonencodable == 0) return 0;
    for (size_t ni = 0; ni < index->n_nonencodable; ++ni) {
        size_t j = index->nonencodable_targets[ni];
        if (!lengths_within_radius(read_len, index->target_lens[j], k)) continue;
        if (index_verify_seed_candidate(index, j, seen, read, read_len, k, result, best_tie_count, stats) != 0) {
            return -1;
        }
        if (stop_on_ambiguous && *best_tie_count > 1) return 1;
    }
    return 0;
}

static int index_assign_neighbor_two_impl(const qdaln_index *index, const char *read, size_t read_len,
                                          uint64_t read_code, qdaln_match_result *result,
                                          qdaln_index_stats *stats, int stop_on_ambiguous) {
    if (index->n_targets == 0) {
        *result = empty_match_result(QDALN_MATCH_NONE);
        return 1;
    }

    candidate_seen seen;
    candidate_seen_init(&seen);
    int best_tie_count = 0;
    *result = empty_match_result(QDALN_MATCH_NONE);

    k2_exact_visit_ctx ctx = {index, &seen, read, read_len, result, &best_tie_count, stats, stop_on_ambiguous};
    int rc = visit_k2_exact_code(read_code, read_len, &ctx);
    if (rc < 0) {
        candidate_seen_free(&seen);
        return 0;
    }
    if (rc > 0) goto done;
    rc = index_visit_one_edit_codes(read_code, read_len, visit_k2_first_edit_code, &ctx);
    if (rc < 0) {
        candidate_seen_free(&seen);
        return 0;
    }
    if (rc > 0) goto done;
    rc = index_scan_nonencodable_seed_candidates(index, &seen, read, read_len, 2, result, &best_tie_count,
                                                stats, stop_on_ambiguous);
    if (rc < 0) {
        candidate_seen_free(&seen);
        return 0;
    }
    if (rc > 0) goto done;

done:
    index_finalize_result(result, best_tie_count);
    candidate_seen_free(&seen);
    return 1;
}

static int index_assign_neighbor_two(const qdaln_index *index, const char *read, size_t read_len,
                                     uint64_t read_code, qdaln_match_result *result,
                                     qdaln_index_stats *stats) {
    return index_assign_neighbor_two_impl(index, read, read_len, read_code, result, stats, 0);
}

static int index_assign_single_unknown_two(const qdaln_index *index, const char *read, size_t read_len,
                                           qdaln_match_result *result, qdaln_index_stats *stats,
                                           int stop_on_ambiguous) {
    size_t bad_pos = 0;
    if ((read == NULL && read_len != 0) || !find_single_non_acgt(read, read_len, &bad_pos)) return 0;
    if (read_len > 32) return 0;
    if (index->n_targets == 0) {
        *result = empty_match_result(QDALN_MATCH_NONE);
        return 1;
    }

    candidate_seen seen;
    candidate_seen_init(&seen);
    int best_tie_count = 0;
    *result = empty_match_result(QDALN_MATCH_NONE);

    static const char bases[] = {'A', 'C', 'G', 'T'};
    char tmp[32];
    memcpy(tmp, read, read_len);
    k2_exact_visit_ctx ctx = {index, &seen, read, read_len, result, &best_tie_count, stats, stop_on_ambiguous};
    for (size_t b = 0; b < 4; ++b) {
        uint64_t code = 0;
        tmp[bad_pos] = bases[b];
        if (!dna2_code(tmp, read_len, &code)) {
            candidate_seen_free(&seen);
            return 0;
        }
        int rc = visit_k2_exact_code(code, read_len, &ctx);
        if (rc < 0) {
            candidate_seen_free(&seen);
            return 0;
        }
        if (rc > 0) goto done;
        rc = index_visit_one_edit_codes(code, read_len, visit_k2_exact_code, &ctx);
        if (rc < 0) {
            candidate_seen_free(&seen);
            return 0;
        }
        if (rc > 0) goto done;
    }
    if (read_len > 0) {
        uint64_t deleted_code = 0;
        if (code_without_base(read, read_len, bad_pos, &deleted_code)) {
            int rc = visit_k2_exact_code(deleted_code, read_len - 1, &ctx);
            if (rc < 0) {
                candidate_seen_free(&seen);
                return 0;
            }
            if (rc > 0) goto done;
            rc = index_visit_one_edit_codes(deleted_code, read_len - 1, visit_k2_exact_code, &ctx);
            if (rc < 0) {
                candidate_seen_free(&seen);
                return 0;
            }
            if (rc > 0) goto done;
        }
    }

    int rc = index_scan_nonencodable_seed_candidates(index, &seen, read, read_len, 2, result, &best_tie_count,
                                                     stats, stop_on_ambiguous);
    if (rc < 0) {
        candidate_seen_free(&seen);
        return 0;
    }
    if (rc > 0) goto done;

done:
    index_finalize_result(result, best_tie_count);
    candidate_seen_free(&seen);
    return 1;
}

static int index_assign_neighbor_two_long_read(const qdaln_index *index, const char *read, size_t read_len,
                                               qdaln_match_result *result, qdaln_index_stats *stats) {
    if (read_len < 33 || read_len > 34) return 0;
    if (count_non_acgt(read, read_len, NULL) != 0) return 0;
    if (index->n_targets == 0) {
        *result = empty_match_result(QDALN_MATCH_NONE);
        return 1;
    }

    candidate_seen seen;
    candidate_seen_init(&seen);
    int best_tie_count = 0;
    *result = empty_match_result(QDALN_MATCH_NONE);

    k2_exact_visit_ctx ctx = {index, &seen, read, read_len, result, &best_tie_count, stats, 0};
    if (read_len == 33) {
        for (size_t drop = 0; drop < read_len; ++drop) {
            uint64_t code = 0;
            if (!code_without_base(read, read_len, drop, &code)) {
                candidate_seen_free(&seen);
                return 0;
            }
            if (visit_k2_exact_code(code, read_len - 1, &ctx) != 0 ||
                index_visit_one_edit_codes(code, read_len - 1, visit_k2_exact_code, &ctx) != 0) {
                candidate_seen_free(&seen);
                return 0;
            }
        }
    } else {
        for (size_t a = 0; a < read_len; ++a) {
            for (size_t b = a + 1; b < read_len; ++b) {
                uint64_t code = 0;
                if (!code_without_two_bases(read, read_len, a, b, &code)) {
                    candidate_seen_free(&seen);
                    return 0;
                }
                if (visit_k2_exact_code(code, read_len - 2, &ctx) != 0) {
                    candidate_seen_free(&seen);
                    return 0;
                }
            }
        }
    }

    if (index_scan_nonencodable_seed_candidates(index, &seen, read, read_len, 2, result, &best_tie_count, stats, 0) != 0) {
        candidate_seen_free(&seen);
        return 0;
    }

    index_finalize_result(result, best_tie_count);
    candidate_seen_free(&seen);
    return 1;
}

typedef int (*qdaln_encoder)(const char *s, size_t len, uint64_t *code_out);

static int bytes_equal_ascii_upper(const char *a, const char *b, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        unsigned char ca = (unsigned char)a[i];
        unsigned char cb = (unsigned char)b[i];
        if (ca >= 'a' && ca <= 'z') ca = (unsigned char)(ca - ('a' - 'A'));
        if (cb >= 'a' && cb <= 'z') cb = (unsigned char)(cb - ('a' - 'A'));
        if (ca != cb) return 0;
    }
    return 1;
}

static int index_assign_exact_one_impl(const qdaln_index *index, const char *read, size_t read_len,
                                       qdaln_match_result *result, qdaln_index_stats *stats,
                                       qdaln_encoder encode, int fold_ascii_fallback) {
    uint64_t read_code = 0;
    *result = empty_match_result(QDALN_MATCH_NONE);
    if (index->n_targets == 0) return 1;
    if (!encode(read, read_len, &read_code)) {
        if (index->n_nonencodable == 0) return 1;
        for (size_t ni = 0; ni < index->n_nonencodable; ++ni) {
            size_t j = index->nonencodable_targets[ni];
            if (index->target_lens[j] != read_len) continue;
            if (stats != NULL) ++stats->candidates_considered;
            if (stats != NULL) ++stats->candidates_verified;
            if (fold_ascii_fallback) {
                if (!bytes_equal_ascii_upper(index->targets[j], read, read_len)) continue;
            } else if (memcmp(index->targets[j], read, read_len) != 0) {
                continue;
            }
            ++result->match_count;
            result->best_distance = 0;
            if (result->target_index < 0 || (int)j < result->target_index) result->target_index = (int)j;
        }
        if (result->match_count == 0) {
            result->status = QDALN_MATCH_NONE;
        } else if (result->match_count > 1) {
            result->status = QDALN_MATCH_AMBIGUOUS;
        } else {
            result->status = QDALN_MATCH_UNIQUE;
        }
        return 1;
    }

    size_t exact_first = QDALN_HASH_EMPTY;
    size_t exact_count = 0;
    if (index_exact_table_find(index, read_code, read_len, &exact_first, &exact_count)) {
        if (stats != NULL) {
            stats->candidates_considered += exact_count;
            stats->candidates_verified += exact_count;
        }
        result->target_index = (int)exact_first;
        result->best_distance = 0;
        result->match_count = (int)exact_count;
        result->status = exact_count > 1 ? QDALN_MATCH_AMBIGUOUS : QDALN_MATCH_UNIQUE;
        return 1;
    }

    if (result->match_count == 0) {
        result->status = QDALN_MATCH_NONE;
    } else if (result->match_count > 1) {
        result->status = QDALN_MATCH_AMBIGUOUS;
    } else {
        result->status = QDALN_MATCH_UNIQUE;
    }
    return 1;
}

static int index_assign_exact_one(const qdaln_index *index, const char *read, size_t read_len,
                                  qdaln_match_result *result, qdaln_index_stats *stats) {
    return index_assign_exact_one_impl(index, read, read_len, result, stats, dna2_code, 0);
}

int qdaln_index_lookup_exact_stats(const qdaln_index *index, const char *read, size_t read_len,
                                   qdaln_match_result *result, qdaln_index_stats *stats) {
    if (index == NULL || result == NULL) return -1;
    if (read == NULL && read_len != 0) {
        *result = empty_match_result(QDALN_MATCH_INVALID);
        if (stats != NULL) {
            stats->candidates_considered = 0;
            stats->candidates_verified = 0;
        }
        return 0;
    }
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }
    return index_assign_exact_one(index, read, read_len, result, stats) ? 0 : -1;
}

int qdaln_index_lookup_exact_many_stats(const qdaln_index *index, const char *const *reads,
                                        const size_t *read_lens, size_t n_reads,
                                        qdaln_match_result *results, qdaln_index_stats *stats) {
    if (index == NULL || (n_reads != 0 && (reads == NULL || read_lens == NULL || results == NULL))) return -1;
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }
    for (size_t i = 0; i < n_reads; ++i) {
        if (reads[i] == NULL && read_lens[i] != 0) {
            results[i] = empty_match_result(QDALN_MATCH_INVALID);
            continue;
        }
        if (!index_assign_exact_one(index, reads[i], read_lens[i], &results[i], stats)) return -1;
    }
    return 0;
}

int qdaln_index_lookup_exact_ascii_many_stats(const qdaln_index *index, const char *const *reads,
                                              const size_t *read_lens, size_t n_reads,
                                              qdaln_match_result *results, qdaln_index_stats *stats) {
    if (index == NULL || (n_reads != 0 && (reads == NULL || read_lens == NULL || results == NULL))) return -1;
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }
    for (size_t i = 0; i < n_reads; ++i) {
        if (reads[i] == NULL && read_lens[i] != 0) {
            results[i] = empty_match_result(QDALN_MATCH_INVALID);
            continue;
        }
        if (!index_assign_exact_one_impl(index, reads[i], read_lens[i], &results[i], stats,
                                         dna2_code_ascii_upper, 1)) {
            return -1;
        }
    }
    return 0;
}

int qdaln_index_lookup_exact_ascii_stats(const qdaln_index *index, const char *read, size_t read_len,
                                         qdaln_match_result *result, qdaln_index_stats *stats) {
    if (index == NULL || result == NULL) return -1;
    if (read == NULL && read_len != 0) {
        *result = empty_match_result(QDALN_MATCH_INVALID);
        if (stats != NULL) {
            stats->candidates_considered = 0;
            stats->candidates_verified = 0;
        }
        return 0;
    }
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }
    return index_assign_exact_one_impl(index, read, read_len, result, stats, dna2_code_ascii_upper, 1) ? 0 : -1;
}

static int hamming_distance_within_k(const char *a, size_t a_len, const char *b, size_t b_len, int k) {
    if (a_len != b_len) return -1;
    if (a_len <= 32) {
        uint64_t a_code = 0;
        uint64_t b_code = 0;
        if (dna2_code(a, a_len, &a_code) && dna2_code(b, b_len, &b_code)) {
            int d = code_hamming_distance_qd(a_code, b_code, a_len);
            return d <= k ? d : -1;
        }
    }
    return same_length_hamming_distance_within_k(a, b, a_len, k);
}

static int index_visit_hamming_exact_seen_candidates(const qdaln_index *index, uint64_t code, size_t len,
                                                     candidate_seen *seen, int distance,
                                                     qdaln_match_result *result, int *best_tie_count,
                                                     qdaln_index_stats *stats) {
    size_t slot = code_hash(code, len, index->hash_cap);
    for (size_t j = index->hash_heads[slot]; j != SIZE_MAX; ) {
        size_t nextj = index->hash_next[j];
        if (!index->encodable[j] || index->target_lens[j] != len || index->codes[j] != code) {
            j = nextj;
            continue;
        }
        int seen_rc = candidate_seen_add(seen, j);
        if (seen_rc < 0) return -1;
        if (seen_rc == 0) {
            j = nextj;
            continue;
        }
        if (stats != NULL) {
            ++stats->candidates_considered;
            ++stats->candidates_verified;
        }
        index_update_verified(result, (int)j, distance, best_tie_count);
        j = nextj;
    }
    return 0;
}

static int index_visit_hamming_seed_candidates(const qdaln_index *index, unsigned char seed_id,
                                               uint64_t seed_code, size_t seed_len, candidate_seen *seen,
                                               uint64_t read_code, size_t read_len, int k,
                                               qdaln_match_result *result, int *best_tie_count,
                                               qdaln_index_stats *stats) {
    if (!index->hamming_seed_ready || index->hamming_seed_hash_cap == 0) return 0;
    size_t slot = hamming_seed_hash(seed_code, seed_len, seed_id, index->hamming_seed_hash_cap);
    for (size_t e = index->hamming_seed_heads[slot]; e != SIZE_MAX; ) {
        size_t nexte = index->hamming_seeds[e].next;
        const qdaln_hamming_seed_entry *entry = &index->hamming_seeds[e];
        size_t j = entry->target_index;
        if (entry->seed_id != seed_id || entry->seed_len != seed_len || entry->code != seed_code ||
            index->target_lens[j] != read_len) {
            e = nexte;
            continue;
        }
        int seen_rc = candidate_seen_add(seen, j);
        if (seen_rc < 0) return -1;
        if (seen_rc == 0) {
            e = nexte;
            continue;
        }
        if (stats != NULL) ++stats->candidates_considered;
        int d = code_hamming_distance_qd(read_code, index->codes[j], read_len);
        if (stats != NULL) ++stats->candidates_verified;
        if (d > k) {
            e = nexte;
            continue;
        }
        index_update_verified(result, (int)j, d, best_tie_count);
        e = nexte;
    }
    return 0;
}

static int index_assign_hamming_seed_one(const qdaln_index *index, const char *read, uint64_t read_code, size_t read_len,
                                         int k, qdaln_match_result *result, qdaln_index_stats *stats) {
    candidate_seen seen;
    candidate_seen_init(&seen);
    int best_tie_count = 0;
    *result = empty_match_result(QDALN_MATCH_NONE);

    if (index_visit_hamming_exact_seen_candidates(index, read_code, read_len, &seen, 0, result,
                                                  &best_tie_count, stats) != 0) {
        candidate_seen_free(&seen);
        return 0;
    }

    size_t n_seeds = (size_t)k + 1U;
    for (size_t seed_id = 0; seed_id < n_seeds; ++seed_id) {
        size_t start = 0;
        size_t seed_len = 0;
        hamming_seed_partition(read_len, n_seeds, seed_id, &start, &seed_len);
        uint64_t seed = code_segment_qd(read_code, start, seed_len);
        if (index_visit_hamming_seed_candidates(index, hamming_seed_key(n_seeds, seed_id), seed, seed_len,
                                                &seen, read_code, read_len, k, result, &best_tie_count,
                                                stats) != 0) {
            candidate_seen_free(&seen);
            return 0;
        }
    }

    if (index->n_nonencodable != 0) {
        for (size_t ni = 0; ni < index->n_nonencodable; ++ni) {
            size_t j = index->nonencodable_targets[ni];
            if (index->target_lens[j] != read_len) continue;
            int seen_rc = candidate_seen_add(&seen, j);
            if (seen_rc < 0) {
                candidate_seen_free(&seen);
                return 0;
            }
            if (seen_rc == 0) continue;
            if (stats != NULL) ++stats->candidates_considered;
            int d = hamming_distance_within_k(read, read_len, index->targets[j], index->target_lens[j], k);
            if (stats != NULL) ++stats->candidates_verified;
            if (d < 0) continue;
            index_update_verified(result, (int)j, d, &best_tie_count);
        }
    }

    index_finalize_result(result, best_tie_count);
    candidate_seen_free(&seen);
    return 1;
}

static int index_assign_hamming_scan_one(const qdaln_index *index, const char *read, size_t read_len, int k,
                                         qdaln_match_result *result, qdaln_index_stats *stats) {
    int best_tie_count = 0;
    uint64_t read_code = 0;
    int read_encodable = dna2_code(read, read_len, &read_code);
    *result = empty_match_result(QDALN_MATCH_NONE);
    for (size_t j = 0; j < index->n_targets; ++j) {
        if (index->target_lens[j] != read_len) continue;
        if (stats != NULL) ++stats->candidates_considered;
        int d;
        if (read_encodable && index->encodable[j]) {
            d = code_hamming_distance_qd(read_code, index->codes[j], read_len);
            if (d > k) d = -1;
        } else {
            d = hamming_distance_within_k(read, read_len, index->targets[j], index->target_lens[j], k);
        }
        if (stats != NULL) ++stats->candidates_verified;
        if (d < 0) continue;
        index_update_verified(result, (int)j, d, &best_tie_count);
    }
    index_finalize_result(result, best_tie_count);
    return 1;
}

typedef int (*hamming_code_visitor)(uint64_t code, void *ctx);

static int visit_hamming_codes_within_radius(uint64_t code, size_t len, int radius,
                                             hamming_code_visitor visit, void *ctx) {
    if (visit(code, ctx) != 0) return -1;
    if (radius <= 0) return 0;

    for (size_t pos = 0; pos < len; ++pos) {
        uint64_t shift = (uint64_t)2 * pos;
        uint64_t old_base = (code >> shift) & 3ULL;
        uint64_t mask = 3ULL << shift;
        for (uint64_t b = 0; b < 4; ++b) {
            if (b == old_base) continue;
            uint64_t code1 = (code & ~mask) | (b << shift);
            if (visit(code1, ctx) != 0) return -1;
            if (radius <= 1) continue;
            for (size_t pos2 = pos + 1; pos2 < len; ++pos2) {
                uint64_t shift2 = (uint64_t)2 * pos2;
                uint64_t old_base2 = (code1 >> shift2) & 3ULL;
                uint64_t mask2 = 3ULL << shift2;
                for (uint64_t b2 = 0; b2 < 4; ++b2) {
                    if (b2 == old_base2) continue;
                    if (visit((code1 & ~mask2) | (b2 << shift2), ctx) != 0) return -1;
                }
            }
        }
    }
    return 0;
}

typedef struct hamming_unknown_visit_ctx {
    const qdaln_index *index;
    candidate_seen *seen;
    const char *read;
    size_t read_len;
    int k;
    qdaln_match_result *result;
    int *best_tie_count;
    qdaln_index_stats *stats;
} hamming_unknown_visit_ctx;

typedef struct hamming_unknown_replacement_ctx {
    hamming_unknown_visit_ctx *visit;
    char *tmp;
    const size_t *bad_positions;
    size_t n_bad;
    int remaining_radius;
} hamming_unknown_replacement_ctx;

static int visit_hamming_unknown_code(uint64_t code, void *ctx) {
    hamming_unknown_visit_ctx *v = (hamming_unknown_visit_ctx *)ctx;
    return index_visit_exact_seed_candidates(v->index, code, v->read_len, v->seen, v->read, v->read_len,
                                             v->k, v->result, v->best_tie_count, v->stats, 0);
}

static int visit_hamming_unknown_replacements(size_t depth, hamming_unknown_replacement_ctx *ctx) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    if (depth == ctx->n_bad) {
        uint64_t code = 0;
        if (!dna2_code(ctx->tmp, ctx->visit->read_len, &code)) return -1;
        return visit_hamming_codes_within_radius(code, ctx->visit->read_len, ctx->remaining_radius,
                                                 visit_hamming_unknown_code, ctx->visit);
    }
    size_t pos = ctx->bad_positions[depth];
    for (size_t b = 0; b < 4; ++b) {
        ctx->tmp[pos] = bases[b];
        if (visit_hamming_unknown_replacements(depth + 1U, ctx) != 0) return -1;
    }
    return 0;
}

static int index_assign_hamming_unknowns(const qdaln_index *index, const char *read, size_t read_len,
                                         int k, qdaln_match_result *result, qdaln_index_stats *stats) {
    *result = empty_match_result(QDALN_MATCH_NONE);
    if ((read == NULL && read_len != 0) || k < 1 || k > 3 || read_len > 32) return 0;
    size_t bad_positions[3];
    size_t n_bad = collect_non_acgt_positions(read, read_len, bad_positions,
                                              sizeof(bad_positions) / sizeof(bad_positions[0]));
    if (n_bad == 0) return 0;

    char tmp[32];
    memcpy(tmp, read, read_len);
    int best_tie_count = 0;
    candidate_seen seen;
    candidate_seen_init(&seen);
    hamming_unknown_visit_ctx ctx = {index, &seen, read, read_len, k, result, &best_tie_count, stats};

    if (n_bad <= (size_t)k) {
        hamming_unknown_replacement_ctx replacement_ctx = {
            &ctx,
            tmp,
            bad_positions,
            n_bad,
            k - (int)n_bad,
        };
        if (visit_hamming_unknown_replacements(0, &replacement_ctx) != 0) {
            candidate_seen_free(&seen);
            return 0;
        }
    }

    if (index->n_nonencodable != 0) {
        for (size_t ni = 0; ni < index->n_nonencodable; ++ni) {
            size_t j = index->nonencodable_targets[ni];
            if (index->target_lens[j] != read_len) continue;
            int seen_rc = candidate_seen_add(&seen, j);
            if (seen_rc < 0) {
                candidate_seen_free(&seen);
                return 0;
            }
            if (seen_rc == 0) continue;
            if (stats != NULL) ++stats->candidates_considered;
            int d = hamming_distance_within_k(read, read_len, index->targets[j], index->target_lens[j], k);
            if (stats != NULL) ++stats->candidates_verified;
            if (d < 0) continue;
            index_update_verified(result, (int)j, d, &best_tie_count);
        }
    }

    index_finalize_result(result, best_tie_count);
    candidate_seen_free(&seen);
    return 1;
}

static int index_assign_hamming_one(const qdaln_index *index, const char *read, size_t read_len, int k,
                                    qdaln_match_result *result, qdaln_index_stats *stats) {
    if (index->n_targets == 0) {
        *result = empty_match_result(QDALN_MATCH_NONE);
        return 1;
    }
    if (k == 0) return index_assign_exact_one(index, read, read_len, result, stats);
    if (k < 1 || k > 3) return 0;

    uint64_t read_code = 0;
    if (!dna2_code(read, read_len, &read_code)) {
        if (index_assign_hamming_unknowns(index, read, read_len, k, result, stats)) return 1;
        return index_assign_hamming_scan_one(index, read, read_len, k, result, stats);
    }

    if (read_len > 32) return index_assign_hamming_scan_one(index, read, read_len, k, result, stats);

    int best_tie_count = 0;
    uint64_t seen_codes[128];
    size_t n_seen_codes = 0;
    *result = empty_match_result(QDALN_MATCH_NONE);
    if (index->hamming_seed_ready && read_len <= 32) {
        return index_assign_hamming_seed_one(index, read, read_code, read_len, k, result, stats);
    }
    seen_codes[n_seen_codes++] = read_code;
    index_visit_hamming_code_candidates(index, read_code, read_len, NULL, 0, result, &best_tie_count, stats);

    for (size_t pos = 0; pos < read_len; ++pos) {
        uint64_t shift = (uint64_t)2 * pos;
        uint64_t old_base = (read_code >> shift) & 3ULL;
        uint64_t mask = 3ULL << shift;
        for (uint64_t b = 0; b < 4; ++b) {
            if (b == old_base) continue;
            uint64_t code = (read_code & ~mask) | (b << shift);
            if (code_already_seen(seen_codes, n_seen_codes, code)) continue;
            if (n_seen_codes < sizeof(seen_codes) / sizeof(seen_codes[0])) seen_codes[n_seen_codes++] = code;
            index_visit_hamming_code_candidates(index, code, read_len, NULL, 1, result, &best_tie_count, stats);
        }
    }

    if (index->n_nonencodable != 0) {
        for (size_t ni = 0; ni < index->n_nonencodable; ++ni) {
            size_t j = index->nonencodable_targets[ni];
            if (index->target_lens[j] != read_len) continue;
            if (stats != NULL) ++stats->candidates_considered;
            int d = hamming_distance_within_k(read, read_len, index->targets[j], index->target_lens[j], k);
            if (stats != NULL) ++stats->candidates_verified;
            if (d < 0) continue;
            index_update_verified(result, (int)j, d, &best_tie_count);
        }
    }

    index_finalize_result(result, best_tie_count);
    return 1;
}

static int index_assign_fast_one(const qdaln_index *index, const char *read, size_t read_len,
                                 int k, qdaln_match_result *result, qdaln_index_stats *stats) {
    if (index->n_targets == 0) {
        *result = empty_match_result(QDALN_MATCH_NONE);
        return 1;
    }
    if (k == 0) return index_assign_exact_one(index, read, read_len, result, stats);

    uint64_t read_code = 0;
    if (!dna2_code(read, read_len, &read_code)) {
        if (k == 2 && index_assign_single_unknown_two(index, read, read_len, result, stats, 0)) return 1;
        if (k == 2 && index_assign_neighbor_two_long_read(index, read, read_len, result, stats)) return 1;
        if (k == 1) return index_assign_single_unknown_one(index, read, read_len, result, stats, 0);
        return 0;
    }

    if (k == 1) return index_assign_neighbor_one(index, read, read_len, read_code, result, stats);
    if (k == 2) return index_assign_neighbor_two(index, read, read_len, read_code, result, stats);

    return 0;
}

int qdaln_index_assign(const qdaln_index *index, const char *const *reads, const size_t *read_lens,
                       size_t n_reads, int k, qdaln_match_result *results) {
    return qdaln_index_assign_stats(index, reads, read_lens, n_reads, k, results, NULL);
}

int qdaln_index_assign_stats(const qdaln_index *index, const char *const *reads, const size_t *read_lens,
                             size_t n_reads, int k, qdaln_match_result *results,
                             qdaln_index_stats *stats) {
    if (index == NULL || results == NULL) return -1;
    if (k < 0) return -1;
    if (n_reads != 0 && (reads == NULL || read_lens == NULL)) return -1;
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }

    for (size_t i = 0; i < n_reads; ++i) {
        if (reads[i] == NULL && read_lens[i] != 0) {
            results[i] = empty_match_result(QDALN_MATCH_INVALID);
            continue;
        }
        if (!index_assign_fast_one(index, reads[i], read_lens[i], k, &results[i], stats)) {
            if (stats != NULL) {
                stats->candidates_considered += index->n_targets;
                stats->candidates_verified += index->n_targets;
            }
            int rc = qdaln_match_many(&reads[i], &read_lens[i], 1, (const char *const *)index->targets,
                                      index->target_lens, index->n_targets, k, &results[i]);
            if (rc != 0) return rc;
        }
    }

    return 0;
}

int qdaln_index_assign_status_stats(const qdaln_index *index, const char *const *reads,
                                    const size_t *read_lens, size_t n_reads, int k,
                                    qdaln_match_result *results, qdaln_index_stats *stats) {
    if (index == NULL || results == NULL) return -1;
    if (k < 0) return -1;
    if (n_reads != 0 && (reads == NULL || read_lens == NULL)) return -1;
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }

    if (k != 1 && k != 2) {
        return qdaln_index_assign_stats(index, reads, read_lens, n_reads, k, results, stats);
    }

    for (size_t i = 0; i < n_reads; ++i) {
        if (reads[i] == NULL && read_lens[i] != 0) {
            results[i] = empty_match_result(QDALN_MATCH_INVALID);
            continue;
        }
        uint64_t read_code = 0;
        if (!dna2_code(reads[i], read_lens[i], &read_code)) {
            if (k == 2) {
                if (!index_assign_single_unknown_two(index, reads[i], read_lens[i], &results[i], stats, 1) &&
                    !index_assign_fast_one(index, reads[i], read_lens[i], k, &results[i], stats)) {
                    if (stats != NULL) {
                        stats->candidates_considered += index->n_targets;
                        stats->candidates_verified += index->n_targets;
                    }
                    int rc = qdaln_match_many(&reads[i], &read_lens[i], 1, (const char *const *)index->targets,
                                              index->target_lens, index->n_targets, k, &results[i]);
                    if (rc != 0) return rc;
                }
                continue;
            }
            if (!index_assign_single_unknown_one(index, reads[i], read_lens[i], &results[i], stats, 1)) {
                if (stats != NULL) {
                    stats->candidates_considered += index->n_targets;
                    stats->candidates_verified += index->n_targets;
                }
                int rc = qdaln_match_many(&reads[i], &read_lens[i], 1, (const char *const *)index->targets,
                                          index->target_lens, index->n_targets, k, &results[i]);
                if (rc != 0) return rc;
            }
            continue;
        }
        if (k == 2) {
            if (!index_assign_neighbor_two_impl(index, reads[i], read_lens[i], read_code, &results[i], stats, 1)) {
                if (stats != NULL) {
                    stats->candidates_considered += index->n_targets;
                    stats->candidates_verified += index->n_targets;
                }
                int rc = qdaln_match_many(&reads[i], &read_lens[i], 1, (const char *const *)index->targets,
                                          index->target_lens, index->n_targets, k, &results[i]);
                if (rc != 0) return rc;
            }
            continue;
        }
        if (!index_assign_neighbor_one_impl(index, reads[i], read_lens[i], read_code, &results[i], stats, 1)) {
            if (stats != NULL) {
                stats->candidates_considered += index->n_targets;
                stats->candidates_verified += index->n_targets;
            }
            int rc = qdaln_match_many(&reads[i], &read_lens[i], 1, (const char *const *)index->targets,
                                      index->target_lens, index->n_targets, k, &results[i]);
            if (rc != 0) return rc;
        }
    }

    return 0;
}

int qdaln_index_assign_hamming_stats(const qdaln_index *index, const char *const *reads,
                                     const size_t *read_lens, size_t n_reads, int k,
                                     qdaln_match_result *results, qdaln_index_stats *stats) {
    if (index == NULL || results == NULL) return -1;
    if (k < 0 || k > 3) return -1;
    if (n_reads != 0 && (reads == NULL || read_lens == NULL)) return -1;
    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
    }
    for (size_t i = 0; i < n_reads; ++i) {
        if (reads[i] == NULL && read_lens[i] != 0) {
            results[i] = empty_match_result(QDALN_MATCH_INVALID);
            continue;
        }
        if (!index_assign_hamming_one(index, reads[i], read_lens[i], k, &results[i], stats)) {
            return -1;
        }
    }
    return 0;
}
