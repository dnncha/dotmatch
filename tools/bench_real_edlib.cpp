#include "edlib.h"
#include "qdalign.h"

#include <zlib.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

struct Target {
    std::string id;
    std::string seq;
};

static void uppercase(std::string &s) {
    for (char &c : s) {
        if (c >= 'a' && c <= 'z') c = (char)(c - 'a' + 'A');
    }
}

static void trim(std::string &s) {
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
}

static std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> out;
    size_t start = 0;
    for (;;) {
        size_t pos = s.find(delim, start);
        if (pos == std::string::npos) {
            out.push_back(s.substr(start));
            return out;
        }
        out.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
}

static bool ieq(std::string a, std::string b) {
    uppercase(a);
    uppercase(b);
    return a == b;
}

static int find_col(const std::vector<std::string> &cols, const char *a, const char *b, const char *c) {
    for (size_t i = 0; i < cols.size(); ++i) {
        if (ieq(cols[i], a) || (b && ieq(cols[i], b)) || (c && ieq(cols[i], c))) return (int)i;
    }
    return -1;
}

static std::vector<Target> read_targets(const char *path) {
    FILE *fp = std::fopen(path, "r");
    if (!fp) std::exit(1);
    std::vector<Target> targets;
    char buf[65536];
    bool first = true;
    int id_col = 0;
    int seq_col = 1;
    while (std::fgets(buf, sizeof(buf), fp)) {
        std::string line(buf);
        trim(line);
        if (line.empty() || line[0] == '#') continue;
        char delim = line.find(',') != std::string::npos && line.find('\t') == std::string::npos ? ',' : '\t';
        std::vector<std::string> cols = split(line, delim);
        if (first) {
            int maybe_id = find_col(cols, "id", "target_id", "sgRNA");
            int maybe_seq = find_col(cols, "gRNA.sequence", "target_seq", "sequence");
            if (maybe_seq < 0) maybe_seq = find_col(cols, "seq", "barcode_seq", "guide_seq");
            if (maybe_id >= 0 && maybe_seq >= 0) {
                id_col = maybe_id;
                seq_col = maybe_seq;
                first = false;
                continue;
            }
        }
        first = false;
        if (cols.size() == 1) {
            targets.push_back(Target{std::to_string(targets.size()), cols[0]});
        } else if ((size_t)id_col < cols.size() && (size_t)seq_col < cols.size()) {
            targets.push_back(Target{cols[(size_t)id_col], cols[(size_t)seq_col]});
        }
        uppercase(targets.back().seq);
    }
    std::fclose(fp);
    return targets;
}

static int read_fastq_record(gzFile gz, std::string &seq) {
    char h[8192], s[8192], p[8192], q[8192];
    char *got = gzgets(gz, h, sizeof(h));
    if (!got) return gzeof(gz) ? 0 : -1;
    if (!gzgets(gz, s, sizeof(s)) || !gzgets(gz, p, sizeof(p)) || !gzgets(gz, q, sizeof(q))) return -1;
    std::string header(h), plus(p), qual(q);
    seq = s;
    trim(header);
    trim(seq);
    trim(plus);
    trim(qual);
    if (header.empty() || plus.empty() || header[0] != '@' || plus[0] != '+' || seq.size() != qual.size()) return -1;
    uppercase(seq);
    return 1;
}

static std::vector<std::string> read_observed(const std::vector<std::string> &paths, size_t start, size_t len, size_t limit) {
    std::vector<std::string> reads;
    for (const std::string &path : paths) {
        gzFile gz = gzopen(path.c_str(), "rb");
        if (!gz) std::exit(1);
        std::string seq;
        int got = 0;
        while ((limit == 0 || reads.size() < limit) && (got = read_fastq_record(gz, seq)) == 1) {
            if (start <= seq.size() && len <= seq.size() - start) reads.push_back(seq.substr(start, len));
        }
        gzclose(gz);
        if (got < 0) std::exit(1);
        if (limit != 0 && reads.size() >= limit) break;
    }
    return reads;
}

static double seconds_now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static qdaln_match_result empty_result(void) {
    qdaln_match_result r;
    r.target_index = -1;
    r.best_distance = -1;
    r.second_best_distance = -1;
    r.match_count = 0;
    r.status = QDALN_MATCH_NONE;
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

static void edlib_assign(const std::vector<std::string> &reads, const std::vector<Target> &targets, int k,
                         std::vector<qdaln_match_result> &results) {
    EdlibAlignConfig config = edlibNewAlignConfig(k, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
    for (size_t i = 0; i < reads.size(); ++i) {
        results[i] = empty_result();
        int best_ties = 0;
        for (size_t j = 0; j < targets.size(); ++j) {
            EdlibAlignResult r = edlibAlign(reads[i].data(), (int)reads[i].size(), targets[j].seq.data(), (int)targets[j].seq.size(), config);
            if (r.status != EDLIB_STATUS_OK) std::exit(1);
            int d = r.editDistance;
            edlibFreeAlignResult(r);
            if (d >= 0) consider(&results[i], (int)j, d, &best_ties);
        }
        finalize(&results[i], best_ties);
    }
}

static long checksum(const std::vector<qdaln_match_result> &results) {
    long out = 0;
    for (const qdaln_match_result &r : results) {
        out += (long)(r.target_index + 1) * 17L;
        out += (long)(r.best_distance + 1) * 31L;
        out += (long)r.status * 43L;
        out += (long)r.match_count * 7L;
    }
    return out;
}

static bool same(qdaln_match_result a, qdaln_match_result b) {
    return a.target_index == b.target_index && a.best_distance == b.best_distance &&
           a.second_best_distance == b.second_best_distance && a.match_count == b.match_count &&
           a.status == b.status;
}

int main(int argc, char **argv) {
    const char *targets_path = NULL;
    std::vector<std::string> read_paths;
    size_t target_start = 23;
    size_t target_len = 19;
    size_t limit = 50;
    int k = 1;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--targets") == 0 && i + 1 < argc) targets_path = argv[++i];
        else if (strcmp(argv[i], "--reads") == 0 && i + 1 < argc) read_paths.push_back(argv[++i]);
        else if (strcmp(argv[i], "--target-start") == 0 && i + 1 < argc) target_start = (size_t)strtoull(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--target-length") == 0 && i + 1 < argc) target_len = (size_t)strtoull(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--limit") == 0 && i + 1 < argc) limit = (size_t)strtoull(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--k") == 0 && i + 1 < argc) k = atoi(argv[++i]);
        else return 2;
    }
    if (!targets_path || read_paths.empty() || k < 0) return 2;

    std::vector<Target> targets = read_targets(targets_path);
    std::vector<std::string> reads = read_observed(read_paths, target_start, target_len, limit);
    std::vector<const char *> target_ptrs, read_ptrs;
    std::vector<size_t> target_lens, read_lens;
    for (const Target &t : targets) {
        target_ptrs.push_back(t.seq.data());
        target_lens.push_back(t.seq.size());
    }
    for (const std::string &r : reads) {
        read_ptrs.push_back(r.data());
        read_lens.push_back(r.size());
    }

    std::vector<qdaln_match_result> indexed(reads.size());
    std::vector<qdaln_match_result> edlib(reads.size());
    qdaln_index *index = qdaln_index_build(target_ptrs.data(), target_lens.data(), targets.size());
    qdaln_index_stats stats;
    double start = seconds_now();
    if (qdaln_index_assign_stats(index, read_ptrs.data(), read_lens.data(), reads.size(), k, indexed.data(), &stats) != 0) return 1;
    double indexed_seconds = seconds_now() - start;

    start = seconds_now();
    edlib_assign(reads, targets, k, edlib);
    double edlib_seconds = seconds_now() - start;
    size_t mismatches = 0;
    for (size_t i = 0; i < reads.size(); ++i) {
        if (!same(indexed[i], edlib[i])) ++mismatches;
    }
    qdaln_index_free(index);

    printf("tool,workflow,n_reads,n_targets,len,k,seconds,reads_per_sec,candidates_per_read,verified_per_read,checksum,mismatches\n");
    printf("dotmatch_indexed,public_crispr_yusa,%zu,%zu,%zu,%d,%.6f,%.1f,%.2f,%.2f,%ld,%zu\n",
           reads.size(), targets.size(), target_len, k, indexed_seconds, (double)reads.size() / indexed_seconds,
           (double)stats.candidates_considered / (double)reads.size(), (double)stats.candidates_verified / (double)reads.size(),
           checksum(indexed), mismatches);
    printf("edlib_native_scan,public_crispr_yusa,%zu,%zu,%zu,%d,%.6f,%.1f,%.2f,%.2f,%ld,%zu\n",
           reads.size(), targets.size(), target_len, k, edlib_seconds, (double)reads.size() / edlib_seconds,
           (double)targets.size(), (double)targets.size(), checksum(edlib), mismatches);
    return mismatches == 0 ? 0 : 1;
}
