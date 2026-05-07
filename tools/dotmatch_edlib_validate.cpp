#include "edlib.h"
#include "qdalign.h"

#include <zlib.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

struct Target {
    std::string id;
    std::string seq;
    std::string gene;
};

struct ValidationStats {
    size_t edlib_alignments = 0;
    size_t bounded_windows = 0;
    size_t fallback_windows = 0;
};

struct EdlibCandidateOracle {
    std::unordered_map<std::string, std::vector<int>> exact;
    bool all_targets_bounded = true;
};

static bool ends_with(const std::string &s, const char *suffix) {
    size_t n = s.size();
    size_t m = std::strlen(suffix);
    return n >= m && s.compare(n - m, m, suffix) == 0;
}

static void trim(std::string &s) {
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
}

static void uppercase(std::string &s) {
    for (char &c : s) {
        if (c >= 'a' && c <= 'z') c = (char)(c - 'a' + 'A');
    }
}

static bool is_acgt(const std::string &s) {
    for (char c : s) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
    }
    return true;
}

static std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> out;
    size_t start = 0;
    for (;;) {
        size_t pos = s.find(delim, start);
        if (pos == std::string::npos) {
            out.push_back(s.substr(start));
            break;
        }
        out.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return out;
}

static bool ieq(std::string a, std::string b) {
    uppercase(a);
    uppercase(b);
    return a == b;
}

static bool parse_double_arg(const char *s, double &out) {
    char *end = nullptr;
    out = std::strtod(s, &end);
    return end != s && *end == '\0';
}

static int find_col(const std::vector<std::string> &cols, const char *a, const char *b, const char *c) {
    for (size_t i = 0; i < cols.size(); ++i) {
        if (ieq(cols[i], a) || (b && ieq(cols[i], b)) || (c && ieq(cols[i], c))) return (int)i;
    }
    return -1;
}

static std::vector<Target> read_targets(const char *path) {
    FILE *fp = std::fopen(path, "r");
    if (!fp) {
        std::perror("open targets");
        std::exit(1);
    }
    std::vector<Target> targets;
    char buf[16384];
    bool first = true;
    int id_col = 0;
    int seq_col = 1;
    int gene_col = 2;
    bool header = false;
    while (std::fgets(buf, sizeof(buf), fp)) {
        std::string line(buf);
        trim(line);
        if (line.empty() || line[0] == '#') continue;
        char delim = line.find(',') != std::string::npos && line.find('\t') == std::string::npos ? ',' : '\t';
        std::vector<std::string> cols = split(line, delim);
        if (first) {
            int maybe_id = find_col(cols, "id", "target_id", "sgRNA");
            if (maybe_id < 0) maybe_id = find_col(cols, "sgRNAID", "sgRNA_ID", "guide_id");
            int maybe_seq = find_col(cols, "gRNA.sequence", "target_seq", "sequence");
            if (maybe_seq < 0) maybe_seq = find_col(cols, "seq", "barcode_seq", "guide_seq");
            int maybe_gene = find_col(cols, "Gene", "gene", nullptr);
            if (maybe_id >= 0 && maybe_seq >= 0) {
                id_col = maybe_id;
                seq_col = maybe_seq;
                gene_col = maybe_gene;
                header = true;
                first = false;
                continue;
            }
        }
        first = false;
        Target t;
        if (cols.size() == 1) {
            t.id = std::to_string(targets.size());
            t.seq = cols[0];
        } else {
            if ((size_t)id_col >= cols.size() || (size_t)seq_col >= cols.size()) std::exit(1);
            t.id = cols[(size_t)id_col];
            t.seq = cols[(size_t)seq_col];
            if (header && gene_col >= 0 && (size_t)gene_col < cols.size()) t.gene = cols[(size_t)gene_col];
            if (!header && cols.size() > 2) t.gene = cols[2];
        }
        uppercase(t.seq);
        targets.push_back(t);
    }
    std::fclose(fp);
    if (targets.empty()) {
        std::fprintf(stderr, "no targets\n");
        std::exit(1);
    }
    return targets;
}

struct FastqReader {
    FILE *fp = nullptr;
    gzFile gz = nullptr;
    bool is_gz = false;
};

static bool open_fastq(FastqReader &r, const char *path) {
    r.is_gz = ends_with(path, ".gz");
    if (r.is_gz) {
        r.gz = gzopen(path, "rb");
        return r.gz != nullptr;
    }
    r.fp = std::fopen(path, "r");
    return r.fp != nullptr;
}

static void close_fastq(FastqReader &r) {
    if (r.gz) gzclose(r.gz);
    if (r.fp) std::fclose(r.fp);
}

static int getline_fastq(FastqReader &r, char *buf, size_t cap) {
    if (r.is_gz) {
        char *got = gzgets(r.gz, buf, (int)cap);
        if (!got) return gzeof(r.gz) ? 0 : -1;
        return 1;
    }
    if (!std::fgets(buf, (int)cap, r.fp)) return std::ferror(r.fp) ? -1 : 0;
    return 1;
}

static int read_fastq_record(FastqReader &r, std::string &seq) {
    char h[8192], s[8192], p[8192], q[8192];
    int got = getline_fastq(r, h, sizeof(h));
    if (got <= 0) return got;
    if (getline_fastq(r, s, sizeof(s)) != 1 || getline_fastq(r, p, sizeof(p)) != 1 || getline_fastq(r, q, sizeof(q)) != 1) return -1;
    std::string hs(h), ps(p), qs(q);
    seq = s;
    trim(hs);
    trim(seq);
    trim(ps);
    trim(qs);
    if (hs.empty() || ps.empty() || hs[0] != '@' || ps[0] != '+' || seq.size() != qs.size()) return -1;
    uppercase(seq);
    return 1;
}

static qdaln_match_result empty_result() {
    qdaln_match_result r;
    r.target_index = -1;
    r.best_distance = -1;
    r.second_best_distance = -1;
    r.match_count = 0;
    r.status = QDALN_MATCH_NONE;
    return r;
}

static qdaln_match_result invalid_result() {
    qdaln_match_result r;
    r.target_index = -1;
    r.best_distance = -1;
    r.second_best_distance = -1;
    r.match_count = 0;
    r.status = QDALN_MATCH_INVALID;
    return r;
}

static EdlibCandidateOracle build_edlib_candidate_oracle(const std::vector<Target> &targets) {
    EdlibCandidateOracle oracle;
    oracle.exact.reserve(targets.size() * 2 + 1);
    for (size_t i = 0; i < targets.size(); ++i) {
        if (targets[i].seq.size() > 32 || !is_acgt(targets[i].seq)) {
            oracle.all_targets_bounded = false;
        }
        oracle.exact[targets[i].seq].push_back((int)i);
    }
    return oracle;
}

static void add_exact_bucket(const EdlibCandidateOracle &oracle, const std::string &seq,
                             std::vector<int> &candidates) {
    auto it = oracle.exact.find(seq);
    if (it == oracle.exact.end()) return;
    for (int target_index : it->second) {
        bool seen = false;
        for (int existing : candidates) {
            if (existing == target_index) {
                seen = true;
                break;
            }
        }
        if (!seen) candidates.push_back(target_index);
    }
}

static std::vector<int> bounded_candidate_targets(const EdlibCandidateOracle &oracle,
                                                  const std::string &read, int k) {
    std::vector<int> candidates;
    add_exact_bucket(oracle, read, candidates);
    if (k == 0) {
        std::sort(candidates.begin(), candidates.end());
        return candidates;
    }

    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::string edited = read;
    for (size_t pos = 0; pos < read.size(); ++pos) {
        char old = edited[pos];
        for (char base : bases) {
            if (base == old) continue;
            edited[pos] = base;
            add_exact_bucket(oracle, edited, candidates);
        }
        edited[pos] = old;
    }

    if (!read.empty()) {
        for (size_t pos = 0; pos < read.size(); ++pos) {
            edited = read;
            edited.erase(pos, 1);
            add_exact_bucket(oracle, edited, candidates);
        }
    }

    if (read.size() < 32) {
        for (size_t pos = 0; pos <= read.size(); ++pos) {
            for (char base : bases) {
                edited = read;
                edited.insert(edited.begin() + (std::ptrdiff_t)pos, base);
                add_exact_bucket(oracle, edited, candidates);
            }
        }
    }

    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    return candidates;
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

static void merge_window_result(qdaln_match_result *best, qdaln_match_result r) {
    if (r.status == QDALN_MATCH_INVALID) return;
    if (r.match_count == 0) {
        if (best->status == QDALN_MATCH_INVALID) *best = r;
        return;
    }
    if (best->match_count == 0 || best->best_distance < 0 || r.best_distance < best->best_distance) {
        *best = r;
        return;
    }
    if (r.best_distance == best->best_distance) {
        if (r.target_index != best->target_index || r.status == QDALN_MATCH_AMBIGUOUS ||
            best->status == QDALN_MATCH_AMBIGUOUS) {
            if (r.target_index >= 0 && (best->target_index < 0 || r.target_index < best->target_index)) {
                best->target_index = r.target_index;
            }
            best->status = QDALN_MATCH_AMBIGUOUS;
            best->match_count += r.match_count;
        }
    } else if (best->second_best_distance < 0 || r.best_distance < best->second_best_distance) {
        best->second_best_distance = r.best_distance;
    }
}

static qdaln_match_result edlib_assign_targets(const std::string &read, const std::vector<Target> &targets,
                                               const std::vector<int> &candidate_indices, int k,
                                               ValidationStats *stats) {
    qdaln_match_result result = empty_result();
    int best_ties = 0;
    EdlibAlignConfig config = edlibNewAlignConfig(k, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0);
    for (int target_index : candidate_indices) {
        const Target &target = targets[(size_t)target_index];
        EdlibAlignResult r = edlibAlign(read.data(), (int)read.size(), target.seq.data(), (int)target.seq.size(), config);
        if (r.status != EDLIB_STATUS_OK) std::exit(1);
        if (stats != nullptr) ++stats->edlib_alignments;
        int d = r.editDistance;
        edlibFreeAlignResult(r);
        if (d >= 0) consider(&result, target_index, d, &best_ties);
    }
    finalize(&result, best_ties);
    return result;
}

static qdaln_match_result edlib_assign_full_scan(const std::string &read, const std::vector<Target> &targets,
                                                 int k, ValidationStats *stats) {
    std::vector<int> all_targets;
    all_targets.reserve(targets.size());
    for (size_t i = 0; i < targets.size(); ++i) all_targets.push_back((int)i);
    return edlib_assign_targets(read, targets, all_targets, k, stats);
}

static qdaln_match_result edlib_assign_one(const std::string &read, const std::vector<Target> &targets,
                                           const EdlibCandidateOracle &oracle, int k,
                                           ValidationStats *stats) {
    if (k <= 1 && oracle.all_targets_bounded && read.size() <= 32 && is_acgt(read)) {
        if (stats != nullptr) ++stats->bounded_windows;
        std::vector<int> candidates = bounded_candidate_targets(oracle, read, k);
        return edlib_assign_targets(read, targets, candidates, k, stats);
    }
    if (stats != nullptr) ++stats->fallback_windows;
    return edlib_assign_full_scan(read, targets, k, stats);
}

static bool same_result(qdaln_match_result a, qdaln_match_result b) {
    return a.target_index == b.target_index && a.best_distance == b.best_distance &&
           a.second_best_distance == b.second_best_distance && a.match_count == b.match_count &&
           a.status == b.status;
}

static size_t offset_distance(size_t a, size_t b) {
    return a > b ? a - b : b - a;
}

static std::vector<size_t> detect_offsets(qdaln_index *index, const char *reads_path, size_t target_start,
                                          size_t target_len, size_t range, size_t sample_limit,
                                          const std::string &mode, double min_fraction) {
    if (range == 0) return {target_start};
    std::vector<unsigned long long> scores(range * 2 + 1, 0);
    FastqReader reader;
    if (!open_fastq(reader, reads_path)) std::exit(1);
    size_t checked = 0;
    std::string seq;
    int got = 0;
    while (checked < sample_limit && (got = read_fastq_record(reader, seq)) == 1) {
        for (size_t oi = 0; oi < scores.size(); ++oi) {
            long delta = (long)oi - (long)range;
            if (delta < 0 && target_start < (size_t)(-delta)) continue;
            size_t offset = delta < 0 ? target_start - (size_t)(-delta) : target_start + (size_t)delta;
            if (!(offset <= seq.size() && target_len <= seq.size() - offset)) continue;
            std::string observed = seq.substr(offset, target_len);
            const char *read_ptr = observed.data();
            size_t read_len = observed.size();
            qdaln_match_result r;
            if (qdaln_index_assign(index, &read_ptr, &read_len, 1, 0, &r) != 0) std::exit(1);
            if (r.status == QDALN_MATCH_UNIQUE) ++scores[oi];
        }
        ++checked;
    }
    close_fastq(reader);
    if (got < 0) std::exit(1);

    size_t best_i = range;
    for (size_t oi = 0; oi < scores.size(); ++oi) {
        if (scores[oi] > scores[best_i] ||
            (scores[oi] == scores[best_i] && offset_distance(oi, range) < offset_distance(best_i, range))) {
            best_i = oi;
        }
    }

    std::vector<size_t> offsets;
    if (mode == "multi" && checked != 0) {
        for (size_t oi = 0; oi < scores.size(); ++oi) {
            double fraction = (double)scores[oi] / (double)checked;
            if (fraction + 1e-12 < min_fraction) continue;
            long delta = (long)oi - (long)range;
            if (delta < 0 && target_start < (size_t)(-delta)) continue;
            size_t offset = delta < 0 ? target_start - (size_t)(-delta) : target_start + (size_t)delta;
            bool seen = false;
            for (size_t existing : offsets) {
                if (existing == offset) {
                    seen = true;
                    break;
                }
            }
            if (!seen) offsets.push_back(offset);
        }
    }
    if (offsets.empty() && scores[best_i] != 0) {
        long delta = (long)best_i - (long)range;
        offsets.push_back(delta < 0 ? target_start - (size_t)(-delta) : target_start + (size_t)delta);
    }
    if (offsets.empty()) offsets.push_back(target_start);
    return offsets;
}

static bool validate_sequence(const std::string &seq, const std::vector<size_t> &offsets,
                              qdaln_index *index, const std::vector<Target> &targets,
                              const EdlibCandidateOracle &oracle, size_t target_len,
                              size_t indel_window, int k, ValidationStats *stats) {
    size_t min_len = target_len;
    size_t max_len = target_len;
    if (indel_window != 0 && k == 1) {
        min_len = target_len > indel_window ? target_len - indel_window : 0;
        max_len = target_len + indel_window;
    }
    qdaln_match_result indexed = invalid_result();
    qdaln_match_result edlib = invalid_result();
    for (size_t offset : offsets) {
        for (size_t len = min_len; len <= max_len; ++len) {
            if (!(offset <= seq.size() && len <= seq.size() - offset)) continue;
            std::string observed = seq.substr(offset, len);
            const char *read_ptr = observed.data();
            size_t read_len = observed.size();
            qdaln_match_result indexed_one;
            qdaln_index_stats indexed_stats;
            if (qdaln_index_assign_stats(index, &read_ptr, &read_len, 1, k, &indexed_one, &indexed_stats) != 0) std::exit(1);
            qdaln_match_result edlib_one = edlib_assign_one(observed, targets, oracle, k, stats);
            merge_window_result(&indexed, indexed_one);
            merge_window_result(&edlib, edlib_one);
        }
    }
    return same_result(indexed, edlib);
}

int main(int argc, char **argv) {
    const char *targets_path = nullptr;
    const char *reads_path = nullptr;
    size_t target_start = 0;
    size_t target_len = 0;
    size_t indel_window = 0;
    size_t sample_limit = 100000;
    size_t auto_offset = 0;
    size_t auto_offset_sample = 1000;
    size_t threads = 1;
    std::string offset_mode = "best";
    double offset_min_fraction = 0.005;
    int k = -1;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--targets") == 0 && i + 1 < argc) targets_path = argv[++i];
        else if (std::strcmp(argv[i], "--reads") == 0 && i + 1 < argc) reads_path = argv[++i];
        else if (std::strcmp(argv[i], "--target-start") == 0 && i + 1 < argc) target_start = (size_t)std::strtoull(argv[++i], nullptr, 10);
        else if (std::strcmp(argv[i], "--target-length") == 0 && i + 1 < argc) target_len = (size_t)std::strtoull(argv[++i], nullptr, 10);
        else if (std::strcmp(argv[i], "--indel-window") == 0 && i + 1 < argc) indel_window = (size_t)std::strtoull(argv[++i], nullptr, 10);
        else if (std::strcmp(argv[i], "--auto-offset") == 0 && i + 1 < argc) auto_offset = (size_t)std::strtoull(argv[++i], nullptr, 10);
        else if (std::strcmp(argv[i], "--auto-offset-sample") == 0 && i + 1 < argc) auto_offset_sample = (size_t)std::strtoull(argv[++i], nullptr, 10);
        else if (std::strcmp(argv[i], "--offset-mode") == 0 && i + 1 < argc) offset_mode = argv[++i];
        else if (std::strcmp(argv[i], "--offset-min-fraction") == 0 && i + 1 < argc) {
            if (!parse_double_arg(argv[++i], offset_min_fraction)) return 2;
        }
        else if (std::strcmp(argv[i], "--k") == 0 && i + 1 < argc) k = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--sample") == 0 && i + 1 < argc) sample_limit = (size_t)std::strtoull(argv[++i], nullptr, 10);
        else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) threads = (size_t)std::strtoull(argv[++i], nullptr, 10);
        else {
            std::fprintf(stderr, "Usage: %s --targets targets.tsv --reads reads.fastq[.gz] --target-start N --target-length L --k 0|1 [--indel-window N] [--auto-offset N] [--offset-mode best|multi] [--sample N] [--threads N]\n", argv[0]);
            return 2;
        }
    }
    if (!targets_path || !reads_path || target_len == 0 || (k != 0 && k != 1)) return 2;
    if (offset_mode != "best" && offset_mode != "multi") return 2;
    if (auto_offset_sample == 0 || offset_min_fraction < 0.0 || offset_min_fraction > 1.0) return 2;
    if (threads == 0) return 2;

    std::vector<Target> targets = read_targets(targets_path);
    std::vector<const char *> target_ptrs;
    std::vector<size_t> target_lens;
    for (const Target &t : targets) {
        target_ptrs.push_back(t.seq.data());
        target_lens.push_back(t.seq.size());
    }
    qdaln_index *index = qdaln_index_build(target_ptrs.data(), target_lens.data(), targets.size());
    if (!index) return 1;
    EdlibCandidateOracle oracle = build_edlib_candidate_oracle(targets);
    std::vector<size_t> offsets = detect_offsets(index, reads_path, target_start, target_len, auto_offset,
                                                 auto_offset_sample, offset_mode, offset_min_fraction);

    FastqReader reader;
    if (!open_fastq(reader, reads_path)) return 1;
    std::vector<std::string> seqs;
    if (sample_limit != 0) seqs.reserve(sample_limit);
    std::string seq;
    int got = 0;
    while ((sample_limit == 0 || seqs.size() < sample_limit) && (got = read_fastq_record(reader, seq)) == 1) {
        seqs.push_back(seq);
    }
    close_fastq(reader);
    if (got < 0) return 1;

    size_t checked = seqs.size();
    size_t worker_count = threads;
    if (worker_count > checked && checked != 0) worker_count = checked;
    if (worker_count == 0) worker_count = 1;
    std::vector<size_t> worker_mismatches(worker_count, 0);
    std::vector<ValidationStats> worker_stats(worker_count);
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (size_t tid = 0; tid < worker_count; ++tid) {
        workers.emplace_back([&, tid]() {
            size_t local_mismatches = 0;
            for (size_t i = tid; i < seqs.size(); i += worker_count) {
                if (!validate_sequence(seqs[i], offsets, index, targets, oracle, target_len,
                                       indel_window, k, &worker_stats[tid])) {
                    ++local_mismatches;
                }
            }
            worker_mismatches[tid] = local_mismatches;
        });
    }
    for (std::thread &worker : workers) worker.join();
    size_t mismatches = 0;
    for (size_t value : worker_mismatches) mismatches += value;
    ValidationStats stats;
    for (const ValidationStats &worker : worker_stats) {
        stats.edlib_alignments += worker.edlib_alignments;
        stats.bounded_windows += worker.bounded_windows;
        stats.fallback_windows += worker.fallback_windows;
    }
    qdaln_index_free(index);
    std::printf("{\n  \"oracle\": \"edlib_native\",\n  \"checked_reads\": %zu,\n  \"mismatches\": %zu,\n  \"k\": %d,\n  \"target_start\": %zu,\n  \"target_length\": %zu,\n  \"indel_window\": %zu,\n  \"offset_mode\": \"%s\",\n  \"selected_target_starts\": [",
                checked, mismatches, k, target_start, target_len, indel_window, offset_mode.c_str());
    for (size_t i = 0; i < offsets.size(); ++i) {
        if (i != 0) std::printf(", ");
        std::printf("%zu", offsets[i]);
    }
    std::printf("],\n  \"oracle_strategy\": \"bounded_edlib_candidates\",\n  \"edlib_alignments\": %zu,\n  \"bounded_windows\": %zu,\n  \"fallback_windows\": %zu\n}\n",
                stats.edlib_alignments, stats.bounded_windows, stats.fallback_windows);
    return mismatches == 0 ? 0 : 1;
}
