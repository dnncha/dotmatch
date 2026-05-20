#include "qdalign.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <zlib.h>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

struct GpuResult {
    int target_index;
    int best_distance;
    int second_best_distance;
    int match_count;
    int status;
};

struct GpuParams {
    uint32_t n_targets;
    uint32_t len;
    uint32_t k;
    uint32_t reserved;
};

struct Counts {
    size_t assigned_unique = 0;
    size_t assigned_exact = 0;
    size_t assigned_corrected = 0;
    size_t ambiguous = 0;
    size_t unmatched = 0;
};

struct Inputs {
    std::vector<std::string> targets;
    std::vector<uint64_t> target_codes;
    std::vector<std::string> reads;
    std::vector<uint64_t> read_codes;
    size_t total_reads = 0;
    size_t invalid_windows = 0;
    size_t non_acgt_windows = 0;
    size_t skipped_targets = 0;
};

struct Args {
    std::string targets = "examples/crispr_guides/data/yusa_library.csv";
    std::vector<std::string> reads;
    size_t target_start = 23;
    size_t target_length = 19;
    int k = 1;
    size_t max_reads_per_sample = 10000;
};

static double seconds_now(void) {
    using clock = std::chrono::steady_clock;
    static const auto start = clock::now();
    auto now = clock::now();
    return std::chrono::duration<double>(now - start).count();
}

static std::string clean_csv(const char *value) {
    std::string s = value == NULL ? "" : value;
    for (char &c : s) {
        if (c == ',' || c == '\n' || c == '\r') c = ' ';
    }
    return s;
}

static std::string trim_line(const char *line) {
    std::string s(line == NULL ? "" : line);
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return s;
}

static std::vector<std::string> split_csv_simple(const std::string &line) {
    std::vector<std::string> out;
    std::string cell;
    for (char c : line) {
        if (c == ',') {
            out.push_back(cell);
            cell.clear();
        } else {
            cell.push_back(c);
        }
    }
    out.push_back(cell);
    return out;
}

static bool pack_dna2(const std::string &s, uint64_t *out) {
    if (s.size() > 32) return false;
    uint64_t code = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        uint64_t bits = 0;
        switch (s[i]) {
            case 'A': bits = 0; break;
            case 'C': bits = 1; break;
            case 'G': bits = 2; break;
            case 'T': bits = 3; break;
            default: return false;
        }
        code |= bits << (2 * i);
    }
    *out = code;
    return true;
}

static int parse_args(int argc, char **argv, Args *args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--targets" && i + 1 < argc) {
            args->targets = argv[++i];
        } else if (arg == "--reads" && i + 1 < argc) {
            args->reads.push_back(argv[++i]);
        } else if (arg == "--target-start" && i + 1 < argc) {
            args->target_start = (size_t)strtoull(argv[++i], NULL, 10);
        } else if (arg == "--target-length" && i + 1 < argc) {
            args->target_length = (size_t)strtoull(argv[++i], NULL, 10);
        } else if (arg == "--k" && i + 1 < argc) {
            args->k = atoi(argv[++i]);
        } else if (arg == "--max-reads-per-sample" && i + 1 < argc) {
            args->max_reads_per_sample = (size_t)strtoull(argv[++i], NULL, 10);
        } else {
            return -1;
        }
    }
    if (args->reads.empty()) {
        args->reads.push_back("examples/crispr_guides/data/ERR376998.fastq.gz");
        args->reads.push_back("examples/crispr_guides/data/ERR376999.fastq.gz");
    }
    if (args->target_length == 0 || args->target_length > 32 || args->k < 0 || args->k > 1) return -1;
    return 0;
}

static int load_targets(const Args &args, Inputs *inputs) {
    FILE *fh = fopen(args.targets.c_str(), "r");
    if (fh == NULL) return -1;
    char buf[8192];
    if (fgets(buf, sizeof(buf), fh) == NULL) {
        fclose(fh);
        return -1;
    }
    std::vector<std::string> header = split_csv_simple(trim_line(buf));
    size_t seq_col = SIZE_MAX;
    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == "gRNA.sequence" || header[i] == "target_seq" || header[i] == "sequence") {
            seq_col = i;
            break;
        }
    }
    if (seq_col == SIZE_MAX) {
        fclose(fh);
        return -1;
    }
    while (fgets(buf, sizeof(buf), fh) != NULL) {
        std::vector<std::string> row = split_csv_simple(trim_line(buf));
        if (seq_col >= row.size()) continue;
        std::string seq = row[seq_col];
        if (seq.size() != args.target_length) {
            ++inputs->skipped_targets;
            continue;
        }
        uint64_t code = 0;
        if (!pack_dna2(seq, &code)) {
            ++inputs->skipped_targets;
            continue;
        }
        inputs->targets.push_back(seq);
        inputs->target_codes.push_back(code);
    }
    fclose(fh);
    return inputs->targets.empty() ? -1 : 0;
}

static int load_reads(const Args &args, Inputs *inputs) {
    char header[8192];
    char seq_buf[8192];
    char plus[8192];
    char qual[8192];
    for (const std::string &path : args.reads) {
        gzFile fh = gzopen(path.c_str(), "rb");
        if (fh == NULL) return -1;
        size_t sample_reads = 0;
        while (gzgets(fh, header, sizeof(header)) != NULL) {
            if (gzgets(fh, seq_buf, sizeof(seq_buf)) == NULL ||
                gzgets(fh, plus, sizeof(plus)) == NULL ||
                gzgets(fh, qual, sizeof(qual)) == NULL) {
                gzclose(fh);
                return -1;
            }
            if (args.max_reads_per_sample != 0 && sample_reads >= args.max_reads_per_sample) break;
            ++sample_reads;
            ++inputs->total_reads;
            std::string seq = trim_line(seq_buf);
            if (seq.size() < args.target_start + args.target_length) {
                ++inputs->invalid_windows;
                continue;
            }
            std::string window = seq.substr(args.target_start, args.target_length);
            uint64_t code = 0;
            if (!pack_dna2(window, &code)) {
                ++inputs->non_acgt_windows;
                continue;
            }
            inputs->reads.push_back(window);
            inputs->read_codes.push_back(code);
        }
        gzclose(fh);
    }
    return 0;
}

static long checksum_cpu(const qdaln_match_result *results, size_t n) {
    long checksum = 0;
    for (size_t i = 0; i < n; ++i) {
        checksum += (long)(results[i].target_index + 1) * 17L;
        checksum += (long)(results[i].best_distance + 1) * 31L;
        checksum += (long)results[i].status * 43L;
        checksum += (long)results[i].match_count * 7L;
    }
    return checksum;
}

static long checksum_gpu(const GpuResult *results, size_t n) {
    long checksum = 0;
    for (size_t i = 0; i < n; ++i) {
        checksum += (long)(results[i].target_index + 1) * 17L;
        checksum += (long)(results[i].best_distance + 1) * 31L;
        checksum += (long)results[i].status * 43L;
        checksum += (long)results[i].match_count * 7L;
    }
    return checksum;
}

static Counts count_cpu(const qdaln_match_result *results, size_t n) {
    Counts c;
    for (size_t i = 0; i < n; ++i) {
        if (results[i].status == QDALN_MATCH_UNIQUE) {
            ++c.assigned_unique;
            if (results[i].best_distance == 0) ++c.assigned_exact;
            if (results[i].best_distance == 1) ++c.assigned_corrected;
        } else if (results[i].status == QDALN_MATCH_AMBIGUOUS) {
            ++c.ambiguous;
        } else {
            ++c.unmatched;
        }
    }
    return c;
}

static Counts count_gpu(const GpuResult *results, size_t n) {
    Counts c;
    for (size_t i = 0; i < n; ++i) {
        if (results[i].status == QDALN_MATCH_UNIQUE) {
            ++c.assigned_unique;
            if (results[i].best_distance == 0) ++c.assigned_exact;
            if (results[i].best_distance == 1) ++c.assigned_corrected;
        } else if (results[i].status == QDALN_MATCH_AMBIGUOUS) {
            ++c.ambiguous;
        } else {
            ++c.unmatched;
        }
    }
    return c;
}

static size_t result_mismatches(const qdaln_match_result *cpu, const GpuResult *gpu, size_t n) {
    size_t mismatches = 0;
    for (size_t i = 0; i < n; ++i) {
        if (cpu[i].target_index != gpu[i].target_index ||
            cpu[i].best_distance != gpu[i].best_distance ||
            cpu[i].second_best_distance != gpu[i].second_best_distance ||
            cpu[i].match_count != gpu[i].match_count ||
            cpu[i].status != gpu[i].status) {
            ++mismatches;
        }
    }
    return mismatches;
}

static long count_delta(const Counts &a, const Counts &b) {
    long delta = 0;
    delta += labs((long)a.assigned_unique - (long)b.assigned_unique);
    delta += labs((long)a.assigned_exact - (long)b.assigned_exact);
    delta += labs((long)a.assigned_corrected - (long)b.assigned_corrected);
    delta += labs((long)a.ambiguous - (long)b.ambiguous);
    delta += labs((long)a.unmatched - (long)b.unmatched);
    return delta;
}

static NSString *kernel_source(void) {
    return @"#include <metal_stdlib>\n"
            "using namespace metal;\n"
            "struct GpuResult { int target_index; int best_distance; int second_best_distance; int match_count; int status; };\n"
            "struct GpuParams { uint n_targets; uint len; uint k; uint reserved; };\n"
            "kernel void hamming_assign(device const ulong *reads [[buffer(0)]],\n"
            "                           device const ulong *targets [[buffer(1)]],\n"
            "                           constant GpuParams &params [[buffer(2)]],\n"
            "                           device GpuResult *results [[buffer(3)]],\n"
            "                           uint gid [[thread_position_in_grid]]) {\n"
            "    ulong read = reads[gid];\n"
            "    ulong mask = params.len == 32 ? 0xffffffffffffffffUL : ((1UL << (params.len * 2)) - 1UL);\n"
            "    int best = -1;\n"
            "    int second = -1;\n"
            "    int best_idx = -1;\n"
            "    int match_count = 0;\n"
            "    int best_ties = 0;\n"
            "    for (uint j = 0; j < params.n_targets; ++j) {\n"
            "        ulong diff = (read ^ targets[j]) & mask;\n"
            "        ulong pair_bits = (diff | (diff >> 1)) & 0x5555555555555555UL;\n"
            "        int d = (int)popcount(pair_bits);\n"
            "        if ((uint)d > params.k) continue;\n"
            "        ++match_count;\n"
            "        if (best < 0 || d < best) {\n"
            "            second = best;\n"
            "            best = d;\n"
            "            best_idx = (int)j;\n"
            "            best_ties = 1;\n"
            "        } else if (d == best) {\n"
            "            if (best_idx < 0 || (int)j < best_idx) best_idx = (int)j;\n"
            "            ++best_ties;\n"
            "        } else if (second < 0 || d < second) {\n"
            "            second = d;\n"
            "        }\n"
            "    }\n"
            "    GpuResult out;\n"
            "    out.target_index = best_idx;\n"
            "    out.best_distance = best;\n"
            "    out.second_best_distance = second;\n"
            "    out.match_count = match_count;\n"
            "    out.status = match_count == 0 ? 0 : (best_ties > 1 ? 2 : 1);\n"
            "    results[gid] = out;\n"
            "}\n";
}

static void print_row(const char *tool, const char *backend, const char *status, const Args &args,
                      const Inputs &inputs, double input_seconds, double prep_seconds, double seconds,
                      const Counts &counts, long checksum, size_t mismatches, long delta,
                      const char *device, const char *notes) {
    double total_seconds = input_seconds + prep_seconds + seconds;
    double reads_per_sec = seconds > 0.0 ? (double)inputs.total_reads / seconds : 0.0;
    double total_reads_per_sec = total_seconds > 0.0 ? (double)inputs.total_reads / total_seconds : 0.0;
    std::cout << tool << ','
              << backend << ','
              << status << ",public_crispr_yusa_hamming,"
              << inputs.total_reads << ','
              << inputs.reads.size() << ','
              << inputs.targets.size() << ','
              << args.target_start << ','
              << args.target_length << ','
              << args.k << ','
              << input_seconds << ','
              << prep_seconds << ','
              << seconds << ','
              << total_seconds << ','
              << reads_per_sec << ','
              << total_reads_per_sec << ','
              << counts.assigned_unique << ','
              << counts.assigned_exact << ','
              << counts.assigned_corrected << ','
              << counts.ambiguous << ','
              << counts.unmatched << ','
              << inputs.invalid_windows << ','
              << inputs.non_acgt_windows << ','
              << inputs.skipped_targets << ','
              << checksum << ','
              << mismatches << ','
              << delta << ','
              << clean_csv(device) << ','
              << clean_csv(notes) << '\n';
}

static void print_unavailable(const char *device, const char *notes) {
    Args args;
    Inputs inputs;
    Counts counts;
    print_row("dotmatch_gpu_metal", "metal", "unavailable", args, inputs, 0.0, 0.0, 0.0,
              counts, 0, 0, 0, device, notes);
}

int main(int argc, char **argv) {
    @autoreleasepool {
        Args args;
        if (parse_args(argc, argv, &args) != 0) {
            std::cerr << "usage: " << argv[0] << " [--targets PATH] [--reads PATH ...] "
                      << "[--target-start N] [--target-length N] [--k 0|1] [--max-reads-per-sample N]\n";
            return 2;
        }

        std::cout << "tool,backend,status,workload,total_reads,packable_reads,n_targets,target_start,target_length,k,input_seconds,prep_seconds,seconds,total_seconds,reads_per_sec,total_reads_per_sec,assigned_unique,assigned_exact,assigned_corrected,ambiguous,unmatched,invalid_windows,non_acgt_windows,skipped_targets,checksum,mismatches,count_delta,device,notes\n";

        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (device == nil) {
            print_unavailable("no_metal_device", "MTLCreateSystemDefaultDevice returned nil");
            return 0;
        }

        double input_start = seconds_now();
        Inputs inputs;
        if (load_targets(args, &inputs) != 0 || load_reads(args, &inputs) != 0) {
            print_unavailable([[device name] UTF8String], "failed to load public CRISPR targets or reads");
            return 0;
        }
        double input_seconds = seconds_now() - input_start;
        if (inputs.reads.empty()) {
            print_unavailable([[device name] UTF8String], "no packable extracted CRISPR guide windows");
            return 0;
        }

        std::vector<const char *> target_ptrs(inputs.targets.size());
        std::vector<size_t> target_lens(inputs.targets.size(), args.target_length);
        for (size_t i = 0; i < inputs.targets.size(); ++i) target_ptrs[i] = inputs.targets[i].c_str();
        std::vector<const char *> read_ptrs(inputs.reads.size());
        std::vector<size_t> read_lens(inputs.reads.size(), args.target_length);
        for (size_t i = 0; i < inputs.reads.size(); ++i) read_ptrs[i] = inputs.reads[i].c_str();

        double cpu_prep_start = seconds_now();
        qdaln_index *index = qdaln_index_build(target_ptrs.data(), target_lens.data(), inputs.targets.size());
        double cpu_prep_seconds = seconds_now() - cpu_prep_start;
        if (index == NULL) return 2;
        std::vector<qdaln_match_result> cpu_results(inputs.reads.size());
        qdaln_index_stats stats;
        double cpu_start = seconds_now();
        int cpu_rc = qdaln_index_assign_hamming_stats(index, read_ptrs.data(), read_lens.data(), inputs.reads.size(),
                                                      args.k, cpu_results.data(), &stats);
        Counts cpu_counts = count_cpu(cpu_results.data(), inputs.reads.size());
        long cpu_sum = checksum_cpu(cpu_results.data(), inputs.reads.size());
        double cpu_seconds = seconds_now() - cpu_start;
        if (cpu_rc != 0) return 2;
        print_row("dotmatch_cpu_index", "cpu", "ok", args, inputs, input_seconds, cpu_prep_seconds,
                  cpu_seconds, cpu_counts, cpu_sum, 0, 0, [[device name] UTF8String],
                  "CPU indexed Hamming k=1 public CRISPR baseline");

        NSError *error = nil;
        id<MTLLibrary> library = [device newLibraryWithSource:kernel_source() options:nil error:&error];
        if (library == nil) {
            print_unavailable([[device name] UTF8String],
                              error == nil ? "Metal library compile failed" : [[error localizedDescription] UTF8String]);
            qdaln_index_free(index);
            return 0;
        }
        id<MTLFunction> function = [library newFunctionWithName:@"hamming_assign"];
        id<MTLComputePipelineState> pipeline = [device newComputePipelineStateWithFunction:function error:&error];
        id<MTLCommandQueue> queue = [device newCommandQueue];
        if (pipeline == nil || queue == nil) {
            print_unavailable([[device name] UTF8String],
                              error == nil ? "Metal pipeline unavailable" : [[error localizedDescription] UTF8String]);
            qdaln_index_free(index);
            return 0;
        }

        double gpu_prep_start = seconds_now();
        id<MTLBuffer> read_buffer = [device newBufferWithLength:inputs.read_codes.size() * sizeof(uint64_t)
                                                        options:MTLResourceStorageModeShared];
        id<MTLBuffer> target_buffer = [device newBufferWithLength:inputs.target_codes.size() * sizeof(uint64_t)
                                                          options:MTLResourceStorageModeShared];
        id<MTLBuffer> result_buffer = [device newBufferWithLength:inputs.reads.size() * sizeof(GpuResult)
                                                          options:MTLResourceStorageModeShared];
        id<MTLBuffer> params_buffer = [device newBufferWithLength:sizeof(GpuParams)
                                                          options:MTLResourceStorageModeShared];
        if (read_buffer == nil || target_buffer == nil || result_buffer == nil || params_buffer == nil) return 2;
        memcpy([read_buffer contents], inputs.read_codes.data(), inputs.read_codes.size() * sizeof(uint64_t));
        memcpy([target_buffer contents], inputs.target_codes.data(), inputs.target_codes.size() * sizeof(uint64_t));
        GpuParams params{(uint32_t)inputs.targets.size(), (uint32_t)args.target_length, (uint32_t)args.k, 0};
        memcpy([params_buffer contents], &params, sizeof(params));
        double gpu_prep_seconds = seconds_now() - gpu_prep_start;

        double gpu_start = seconds_now();
        id<MTLCommandBuffer> command_buffer = [queue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoder];
        [encoder setComputePipelineState:pipeline];
        [encoder setBuffer:read_buffer offset:0 atIndex:0];
        [encoder setBuffer:target_buffer offset:0 atIndex:1];
        [encoder setBuffer:params_buffer offset:0 atIndex:2];
        [encoder setBuffer:result_buffer offset:0 atIndex:3];
        NSUInteger threads_per_group = pipeline.maxTotalThreadsPerThreadgroup;
        if (threads_per_group > 256) threads_per_group = 256;
        [encoder dispatchThreads:MTLSizeMake(inputs.reads.size(), 1, 1)
            threadsPerThreadgroup:MTLSizeMake(threads_per_group, 1, 1)];
        [encoder endEncoding];
        [command_buffer commit];
        [command_buffer waitUntilCompleted];
        GpuResult *gpu_results = (GpuResult *)[result_buffer contents];
        Counts gpu_counts = count_gpu(gpu_results, inputs.reads.size());
        long gpu_sum = checksum_gpu(gpu_results, inputs.reads.size());
        size_t mismatches = result_mismatches(cpu_results.data(), gpu_results, inputs.reads.size());
        long delta = count_delta(cpu_counts, gpu_counts);
        double gpu_seconds = seconds_now() - gpu_start;

        print_row("dotmatch_gpu_metal", "metal", "ok", args, inputs, input_seconds, gpu_prep_seconds,
                  gpu_seconds, gpu_counts, gpu_sum, mismatches, delta, [[device name] UTF8String],
                  "Metal public CRISPR FASTQ extract-pack-dispatch-readback-count lane");

        qdaln_index_free(index);
    }
    return 0;
}
