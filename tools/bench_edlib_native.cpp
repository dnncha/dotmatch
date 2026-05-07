#include "edlib.h"
#include "qdalign.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

static uint64_t rng_state = 0x13198a2e03707344ULL;

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

static void mutate_seq(const char *src, char *dst, size_t n, unsigned per_thousand) {
    static const char dna[] = "ACGT";
    for (size_t i = 0; i < n; ++i) {
        char c = src[i];
        if ((xorshift64() % 1000ULL) < per_thousand) {
            char nc = c;
            while (nc == c) nc = dna[xorshift64() & 3ULL];
            c = nc;
        }
        dst[i] = c;
    }
    dst[n] = '\0';
}

static void make_one_substitution(const char *src, char *dst, size_t n) {
    static const char dna[] = "ACGT";
    memcpy(dst, src, n);
    if (n == 0) {
        dst[0] = '\0';
        return;
    }
    size_t pos = (size_t)(xorshift64() % n);
    char nc = dst[pos];
    while (nc == dst[pos]) nc = dna[xorshift64() & 3ULL];
    dst[pos] = nc;
    dst[n] = '\0';
}

static size_t make_one_insertion(const char *src, char *dst, size_t n) {
    size_t pos = (size_t)(xorshift64() % (n + 1));
    memcpy(dst, src, pos);
    dst[pos] = rand_base();
    memcpy(dst + pos + 1, src + pos, n - pos);
    dst[n + 1] = '\0';
    return n + 1;
}

static size_t make_one_deletion(const char *src, char *dst, size_t n) {
    if (n == 0) {
        dst[0] = '\0';
        return 0;
    }
    size_t pos = (size_t)(xorshift64() % n);
    memcpy(dst, src, pos);
    memcpy(dst + pos, src + pos + 1, n - pos - 1);
    dst[n - 1] = '\0';
    return n - 1;
}

static void make_no_match(char *dst, size_t n) {
    memset(dst, 'N', n);
    dst[n] = '\0';
}

static long peak_rss_kb(void) {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) return -1;
#ifdef __APPLE__
    return (long)(usage.ru_maxrss / 1024);
#else
    return (long)usage.ru_maxrss;
#endif
}

static double seconds_now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static int parse_size_arg(const char *s, size_t *out) {
    char *end = NULL;
    unsigned long v = strtoul(s, &end, 10);
    if (end == s || *end != '\0') return -1;
    *out = (size_t)v;
    return 0;
}

static qdaln_match_result empty_result(int status) {
    qdaln_match_result r;
    r.target_index = -1;
    r.best_distance = -1;
    r.second_best_distance = -1;
    r.match_count = 0;
    r.status = status;
    return r;
}

static void consider(qdaln_match_result *r, int target_index, int d, int *best_ties) {
    ++r->match_count;
    if (r->best_distance < 0 || d < r->best_distance) {
        r->second_best_distance = r->best_distance;
        r->best_distance = d;
        r->target_index = target_index;
        *best_ties = 1;
    } else if (d == r->best_distance) {
        ++(*best_ties);
    } else if (r->second_best_distance < 0 || d < r->second_best_distance) {
        r->second_best_distance = d;
    }
}

static void finalize(qdaln_match_result *r, int best_ties) {
    if (r->match_count == 0) r->status = QDALN_MATCH_NONE;
    else if (best_ties > 1) r->status = QDALN_MATCH_AMBIGUOUS;
    else r->status = QDALN_MATCH_UNIQUE;
}

struct ExactBucket {
    int first_index;
    int count;
};

static std::string key_for(const char *s, size_t len) {
    return std::string(s, len);
}

static void exact_hash_assign(const char *const *reads, const size_t *read_lens, size_t n_reads,
                              const char *const *targets, const size_t *target_lens, size_t n_targets,
                              qdaln_match_result *results, double *verified_per_read) {
    std::unordered_map<std::string, ExactBucket> table;
    table.reserve(n_targets * 2 + 1);
    for (size_t i = 0; i < n_targets; ++i) {
        std::string key = key_for(targets[i], target_lens[i]);
        auto it = table.find(key);
        if (it == table.end()) {
            table.emplace(key, ExactBucket{(int)i, 1});
        } else {
            ++it->second.count;
            if ((int)i < it->second.first_index) it->second.first_index = (int)i;
        }
    }

    size_t verified = 0;
    for (size_t i = 0; i < n_reads; ++i) {
        results[i] = empty_result(QDALN_MATCH_NONE);
        auto it = table.find(key_for(reads[i], read_lens[i]));
        if (it == table.end()) continue;
        verified += (size_t)it->second.count;
        results[i].target_index = it->second.first_index;
        results[i].best_distance = 0;
        results[i].match_count = it->second.count;
        results[i].status = it->second.count > 1 ? QDALN_MATCH_AMBIGUOUS : QDALN_MATCH_UNIQUE;
    }
    *verified_per_read = (double)verified / (double)n_reads;
}

struct BkNode {
    int target_index;
    std::unordered_map<int, int> children;
};

struct BkTree {
    std::vector<BkNode> nodes;
    const char *const *targets;
    const size_t *target_lens;
};

static void bk_insert(BkTree *tree, int target_index) {
    if (tree->nodes.empty()) {
        tree->nodes.push_back(BkNode{target_index, std::unordered_map<int, int>()});
        return;
    }
    int node_index = 0;
    for (;;) {
        BkNode &node = tree->nodes[(size_t)node_index];
        int d = qdaln_edit_distance(tree->targets[target_index], tree->target_lens[target_index],
                                    tree->targets[node.target_index], tree->target_lens[node.target_index]);
        auto it = node.children.find(d);
        if (it == node.children.end()) {
            int next = (int)tree->nodes.size();
            node.children.emplace(d, next);
            tree->nodes.push_back(BkNode{target_index, std::unordered_map<int, int>()});
            return;
        }
        node_index = it->second;
    }
}

static void bk_query(const BkTree *tree, int node_index, const char *read, size_t read_len, int k,
                     qdaln_match_result *result, int *best_ties, size_t *visited) {
    const BkNode &node = tree->nodes[(size_t)node_index];
    ++(*visited);
    int d = qdaln_edit_distance(read, read_len, tree->targets[node.target_index], tree->target_lens[node.target_index]);
    if (d <= k) consider(result, node.target_index, d, best_ties);

    int lo = d - k;
    int hi = d + k;
    for (const auto &child : node.children) {
        if (child.first >= lo && child.first <= hi) {
            bk_query(tree, child.second, read, read_len, k, result, best_ties, visited);
        }
    }
}

static void bk_tree_assign(const char *const *reads, const size_t *read_lens, size_t n_reads,
                           const char *const *targets, const size_t *target_lens, size_t n_targets,
                           int k, qdaln_match_result *results, double *visited_per_read) {
    BkTree tree;
    tree.targets = targets;
    tree.target_lens = target_lens;
    tree.nodes.reserve(n_targets);
    for (size_t i = 0; i < n_targets; ++i) bk_insert(&tree, (int)i);

    size_t visited = 0;
    for (size_t i = 0; i < n_reads; ++i) {
        results[i] = empty_result(QDALN_MATCH_NONE);
        int best_ties = 0;
        if (!tree.nodes.empty()) bk_query(&tree, 0, reads[i], read_lens[i], k, &results[i], &best_ties, &visited);
        finalize(&results[i], best_ties);
    }
    *visited_per_read = (double)visited / (double)n_reads;
}

static void neighbor_add_candidate(std::vector<int> &hits, const std::unordered_map<std::string, std::vector<int>> &table,
                                   const std::string &key) {
    auto it = table.find(key);
    if (it == table.end()) return;
    for (int idx : it->second) hits.push_back(idx);
}

static void neighbor_assign(const char *const *reads, const size_t *read_lens, size_t n_reads,
                            const char *const *targets, const size_t *target_lens, size_t n_targets,
                            int k, qdaln_match_result *results, double *visited_per_read) {
    std::unordered_map<std::string, std::vector<int>> table;
    table.reserve(n_targets * 2 + 1);
    for (size_t i = 0; i < n_targets; ++i) table[key_for(targets[i], target_lens[i])].push_back((int)i);

    static const char dna[] = "ACGT";
    size_t visited = 0;
    for (size_t i = 0; i < n_reads; ++i) {
        results[i] = empty_result(QDALN_MATCH_NONE);
        int best_ties = 0;
        std::string read = key_for(reads[i], read_lens[i]);
        std::vector<int> hits;
        neighbor_add_candidate(hits, table, read);
        if (k == 1) {
            for (size_t pos = 0; pos < read.size(); ++pos) {
                char old = read[pos];
                for (char base : dna) {
                    if (base == old) continue;
                    read[pos] = base;
                    neighbor_add_candidate(hits, table, read);
                }
                read[pos] = old;
            }
            for (size_t pos = 0; pos < read.size(); ++pos) {
                std::string deleted = read.substr(0, pos) + read.substr(pos + 1);
                neighbor_add_candidate(hits, table, deleted);
            }
            for (size_t pos = 0; pos <= read.size(); ++pos) {
                for (char base : dna) {
                    std::string inserted = read.substr(0, pos) + base + read.substr(pos);
                    neighbor_add_candidate(hits, table, inserted);
                }
            }
        }
        std::sort(hits.begin(), hits.end());
        hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
        for (int idx : hits) {
            ++visited;
            int d = qdaln_edit_distance(reads[i], read_lens[i], targets[idx], target_lens[idx]);
            if (d <= k) consider(&results[i], idx, d, &best_ties);
        }
        finalize(&results[i], best_ties);
    }
    *visited_per_read = (double)visited / (double)n_reads;
}

static void edlib_assign(const char *const *reads, const size_t *read_lens, size_t n_reads,
                         const char *const *targets, const size_t *target_lens, size_t n_targets,
                         int k, qdaln_match_result *results) {
    EdlibAlignConfig config = edlibNewAlignConfig(k, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
    for (size_t i = 0; i < n_reads; ++i) {
        results[i] = empty_result(QDALN_MATCH_NONE);
        int best_ties = 0;
        for (size_t j = 0; j < n_targets; ++j) {
            EdlibAlignResult r = edlibAlign(reads[i], (int)read_lens[i], targets[j], (int)target_lens[j], config);
            if (r.status != EDLIB_STATUS_OK) exit(1);
            int d = r.editDistance;
            edlibFreeAlignResult(r);
            if (d >= 0) consider(&results[i], (int)j, d, &best_ties);
        }
        finalize(&results[i], best_ties);
    }
}

static long checksum_results(const qdaln_match_result *results, size_t n) {
    long checksum = 0;
    for (size_t i = 0; i < n; ++i) {
        checksum += (long)(results[i].target_index + 1) * 17L;
        checksum += (long)(results[i].best_distance + 1) * 31L;
        checksum += (long)results[i].status * 43L;
        checksum += (long)results[i].match_count * 7L;
    }
    return checksum;
}

static int same_result(qdaln_match_result a, qdaln_match_result b) {
    return a.target_index == b.target_index &&
           a.best_distance == b.best_distance &&
           a.second_best_distance == b.second_best_distance &&
           a.match_count == b.match_count &&
           a.status == b.status;
}

static void assert_same(const qdaln_match_result *a, const qdaln_match_result *b, size_t n,
                        const char *lhs, const char *rhs) {
    for (size_t i = 0; i < n; ++i) {
        if (!same_result(a[i], b[i])) {
            fprintf(stderr, "assignment mismatch %s vs %s at read %zu\n", lhs, rhs, i);
            fprintf(stderr, "%s: idx=%d best=%d second=%d count=%d status=%d\n",
                    lhs, a[i].target_index, a[i].best_distance, a[i].second_best_distance,
                    a[i].match_count, a[i].status);
            fprintf(stderr, "%s: idx=%d best=%d second=%d count=%d status=%d\n",
                    rhs, b[i].target_index, b[i].best_distance, b[i].second_best_distance,
                    b[i].match_count, b[i].status);
            exit(1);
        }
    }
}

static void print_row(const char *tool, const char *error_mode, size_t n_reads, size_t n_targets, size_t len, int k,
                      unsigned err_per_thousand, double indel_rate, double elapsed, double candidates_per_read,
                      double verified_per_read, long checksum) {
    double reads_per_sec = (double)n_reads / elapsed;
    printf("%s,synthetic_barcode,%s,%zu,%zu,%zu,%d,%.3f,%.3f,%.6f,%.1f,%.2f,%.2f,%ld,%ld,0\n",
           tool, error_mode, n_reads, n_targets, len, k, (double)err_per_thousand / 1000.0,
           indel_rate, elapsed, reads_per_sec, candidates_per_read, verified_per_read, peak_rss_kb(), checksum);
}

static void run_case(size_t n_reads, size_t n_targets, size_t len, int k, unsigned err_per_thousand,
                     const char *error_mode) {
    char *target_buf = (char *)malloc(n_targets * (len + 1));
    char *read_buf = (char *)malloc(n_reads * (len + 2));
    const char **targets = (const char **)malloc(n_targets * sizeof(char *));
    const char **reads = (const char **)malloc(n_reads * sizeof(char *));
    size_t *target_lens = (size_t *)malloc(n_targets * sizeof(size_t));
    size_t *read_lens = (size_t *)malloc(n_reads * sizeof(size_t));
    qdaln_match_result *indexed = (qdaln_match_result *)malloc(n_reads * sizeof(qdaln_match_result));
    qdaln_match_result *scan = (qdaln_match_result *)malloc(n_reads * sizeof(qdaln_match_result));
    qdaln_match_result *edlib = (qdaln_match_result *)malloc(n_reads * sizeof(qdaln_match_result));
    qdaln_match_result *baseline = (qdaln_match_result *)malloc(n_reads * sizeof(qdaln_match_result));
    if (!target_buf || !read_buf || !targets || !reads || !target_lens || !read_lens || !indexed || !scan || !edlib || !baseline) exit(1);

    for (size_t i = 0; i < n_targets; ++i) {
        char *target = target_buf + i * (len + 1);
        rand_seq(target, len);
        targets[i] = target;
        target_lens[i] = len;
    }
    if (strcmp(error_mode, "ambiguous") == 0 && n_targets >= 2 && len > 0) {
        memset(target_buf, 'A', len);
        target_buf[len] = '\0';
        char *target1 = target_buf + (len + 1);
        memset(target1, 'A', len);
        target1[0] = 'C';
        target1[len] = '\0';
    }
    for (size_t i = 0; i < n_reads; ++i) {
        size_t idx = (size_t)(xorshift64() % n_targets);
        char *read = read_buf + i * (len + 2);
        if (strcmp(error_mode, "exact") == 0) {
            memcpy(read, targets[idx], len + 1);
            read_lens[i] = len;
        } else if (strcmp(error_mode, "one_substitution") == 0) {
            make_one_substitution(targets[idx], read, len);
            read_lens[i] = len;
        } else if (strcmp(error_mode, "one_insertion") == 0) {
            read_lens[i] = make_one_insertion(targets[idx], read, len);
        } else if (strcmp(error_mode, "one_deletion") == 0) {
            read_lens[i] = make_one_deletion(targets[idx], read, len);
        } else if (strcmp(error_mode, "no_match") == 0) {
            make_no_match(read, len);
            read_lens[i] = len;
        } else if (strcmp(error_mode, "ambiguous") == 0 && len > 0) {
            memset(read, 'A', len);
            read[0] = 'G';
            read[len] = '\0';
            read_lens[i] = len;
        } else {
            mutate_seq(targets[idx], read, len, err_per_thousand);
            read_lens[i] = len;
        }
        reads[i] = read;
    }
    double indel_rate = (strcmp(error_mode, "one_insertion") == 0 || strcmp(error_mode, "one_deletion") == 0) ? 1.0 / (double)len : 0.0;

    qdaln_index *index = qdaln_index_build(targets, target_lens, n_targets);
    if (index == NULL) exit(1);

    qdaln_index_stats stats;
    double start = seconds_now();
    if (qdaln_index_assign_stats(index, reads, read_lens, n_reads, k, indexed, &stats) != 0) exit(1);
    double indexed_elapsed = seconds_now() - start;

    start = seconds_now();
    if (qdaln_match_many(reads, read_lens, n_reads, targets, target_lens, n_targets, k, scan) != 0) exit(1);
    double scan_elapsed = seconds_now() - start;

    start = seconds_now();
    edlib_assign(reads, read_lens, n_reads, targets, target_lens, n_targets, k, edlib);
    double edlib_elapsed = seconds_now() - start;

    assert_same(indexed, edlib, n_reads, "dotmatch_indexed", "edlib_native_scan");
    assert_same(scan, edlib, n_reads, "dotmatch_scan", "edlib_native_scan");

    print_row("dotmatch_indexed", error_mode, n_reads, n_targets, len, k, err_per_thousand, indel_rate, indexed_elapsed,
              (double)stats.candidates_considered / (double)n_reads,
              (double)stats.candidates_verified / (double)n_reads,
              checksum_results(indexed, n_reads));
    print_row("dotmatch_scan", error_mode, n_reads, n_targets, len, k, err_per_thousand, indel_rate, scan_elapsed,
              (double)n_targets, (double)n_targets, checksum_results(scan, n_reads));
    print_row("edlib_native_scan", error_mode, n_reads, n_targets, len, k, err_per_thousand, indel_rate, edlib_elapsed,
              (double)n_targets, (double)n_targets, checksum_results(edlib, n_reads));

    if (k == 0) {
        double verified_per_read = 0.0;
        start = seconds_now();
        exact_hash_assign(reads, read_lens, n_reads, targets, target_lens, n_targets, baseline, &verified_per_read);
        double hash_elapsed = seconds_now() - start;
        assert_same(baseline, edlib, n_reads, "exact_hash_lookup", "edlib_native_scan");
        print_row("exact_hash_lookup", error_mode, n_reads, n_targets, len, k, err_per_thousand, indel_rate, hash_elapsed,
                  verified_per_read, verified_per_read, checksum_results(baseline, n_reads));
    } else if (k == 1) {
        double visited_per_read = 0.0;
        start = seconds_now();
        bk_tree_assign(reads, read_lens, n_reads, targets, target_lens, n_targets, k, baseline, &visited_per_read);
        double bk_elapsed = seconds_now() - start;
        assert_same(baseline, edlib, n_reads, "bk_tree", "edlib_native_scan");
        print_row("bk_tree", error_mode, n_reads, n_targets, len, k, err_per_thousand, indel_rate, bk_elapsed,
                  visited_per_read, visited_per_read, checksum_results(baseline, n_reads));

        start = seconds_now();
        neighbor_assign(reads, read_lens, n_reads, targets, target_lens, n_targets, k, baseline, &visited_per_read);
        double neighbor_elapsed = seconds_now() - start;
        assert_same(baseline, edlib, n_reads, "neighbor_lookup", "edlib_native_scan");
        print_row("neighbor_lookup", error_mode, n_reads, n_targets, len, k, err_per_thousand, indel_rate, neighbor_elapsed,
                  visited_per_read, visited_per_read, checksum_results(baseline, n_reads));
    }

    qdaln_index_free(index);
    free(target_buf);
    free(read_buf);
    free(targets);
    free(reads);
    free(target_lens);
    free(read_lens);
    free(indexed);
    free(scan);
    free(edlib);
    free(baseline);
}

int main(int argc, char **argv) {
    size_t n_reads = 10000;
    if (argc == 2 && parse_size_arg(argv[1], &n_reads) != 0) {
        fprintf(stderr, "Usage: %s [n_reads]\n", argv[0]);
        return 2;
    }

    const char *matrix = getenv("DOTMATCH_NATIVE_MATRIX");
    const size_t smoke_lens[] = {16, 24, 32};
    const size_t paper_lens[] = {12, 16, 20, 24, 32};
    const size_t smoke_target_counts[] = {96, 737, 4096};
    const size_t paper_target_counts[] = {96, 737, 4096, 16384, 65536};
    const size_t *lens = matrix && strcmp(matrix, "paper") == 0 ? paper_lens : smoke_lens;
    const size_t n_lens = matrix && strcmp(matrix, "paper") == 0 ? sizeof(paper_lens) / sizeof(paper_lens[0]) : sizeof(smoke_lens) / sizeof(smoke_lens[0]);
    const size_t *target_counts = matrix && strcmp(matrix, "paper") == 0 ? paper_target_counts : smoke_target_counts;
    const size_t n_target_counts = matrix && strcmp(matrix, "paper") == 0 ? sizeof(paper_target_counts) / sizeof(paper_target_counts[0]) : sizeof(smoke_target_counts) / sizeof(smoke_target_counts[0]);
    const int ks[] = {0, 1};
    const unsigned errs[] = {0, 5, 10, 30};
    const char *smoke_error_modes[] = {"exact", "one_substitution"};
    const char *paper_error_modes[] = {"exact", "one_substitution", "one_insertion", "one_deletion", "no_match", "ambiguous"};
    const char **error_modes = matrix && strcmp(matrix, "paper") == 0 ? paper_error_modes : smoke_error_modes;
    const size_t n_error_modes = matrix && strcmp(matrix, "paper") == 0 ? sizeof(paper_error_modes) / sizeof(paper_error_modes[0]) : sizeof(smoke_error_modes) / sizeof(smoke_error_modes[0]);

    printf("tool,workload,error_mode,n_reads,n_targets,len,k,err,indel_rate,seconds,reads_per_sec,candidates_per_read,verified_per_read,peak_rss_kb,checksum,mismatches\n");
    for (size_t li = 0; li < n_lens; ++li) {
        for (size_t ti = 0; ti < n_target_counts; ++ti) {
            for (size_t ki = 0; ki < sizeof(ks) / sizeof(ks[0]); ++ki) {
                for (size_t em = 0; em < n_error_modes; ++em) {
                    size_t n_errs = strcmp(error_modes[em], "one_substitution") == 0 ? sizeof(errs) / sizeof(errs[0]) : 1;
                    for (size_t ei = 0; ei < n_errs; ++ei) {
                        run_case(n_reads, target_counts[ti], lens[li], ks[ki], errs[ei], error_modes[em]);
                    }
                }
            }
        }
    }
    return 0;
}
