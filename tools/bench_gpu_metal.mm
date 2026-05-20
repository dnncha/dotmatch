#include "qdalign.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

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

struct CaseSpec {
    size_t n_reads;
    size_t n_targets;
    size_t len;
    int k;
    unsigned err_per_thousand;
};

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

static void rand_seq(std::string &s, size_t n) {
    s.resize(n);
    for (size_t i = 0; i < n; ++i) s[i] = rand_base();
}

static void mutate_seq(const std::string &src, std::string &dst, unsigned per_thousand) {
    static const char dna[] = "ACGT";
    dst.resize(src.size());
    for (size_t i = 0; i < src.size(); ++i) {
        char c = src[i];
        if ((xorshift64() % 1000ULL) < per_thousand) {
            char nc = c;
            while (nc == c) nc = dna[xorshift64() & 3ULL];
            c = nc;
        }
        dst[i] = c;
    }
}

static uint64_t pack_dna2(const std::string &s) {
    uint64_t code = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        uint64_t bits = 0;
        switch (s[i]) {
            case 'A': bits = 0; break;
            case 'C': bits = 1; break;
            case 'G': bits = 2; break;
            case 'T': bits = 3; break;
            default: bits = 0; break;
        }
        code |= bits << (2 * i);
    }
    return code;
}

static double seconds_now(void) {
    using clock = std::chrono::steady_clock;
    static const auto start = clock::now();
    auto now = clock::now();
    return std::chrono::duration<double>(now - start).count();
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

static size_t mismatch_count(const qdaln_match_result *cpu, const GpuResult *gpu, size_t n) {
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

static std::string csv_clean(const char *value) {
    std::string s = value == NULL ? "" : value;
    for (char &c : s) {
        if (c == ',' || c == '\n' || c == '\r') c = ' ';
    }
    return s;
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

static void print_row(const char *tool, const char *backend, const char *status, const CaseSpec &spec,
                      double prep_seconds, double seconds, long checksum, size_t mismatches,
                      const char *device, const char *notes) {
    double total_seconds = prep_seconds + seconds;
    double reads_per_sec = seconds > 0.0 ? (double)spec.n_reads / seconds : 0.0;
    double total_reads_per_sec = total_seconds > 0.0 ? (double)spec.n_reads / total_seconds : 0.0;
    double pairs_per_sec = seconds > 0.0 ? ((double)spec.n_reads * (double)spec.n_targets) / seconds : 0.0;
    std::cout << tool << ','
              << backend << ','
              << status << ",synthetic_hamming,"
              << spec.n_reads << ','
              << spec.n_targets << ','
              << spec.len << ','
              << spec.k << ','
              << ((double)spec.err_per_thousand / 1000.0) << ','
              << prep_seconds << ','
              << seconds << ','
              << total_seconds << ','
              << reads_per_sec << ','
              << total_reads_per_sec << ','
              << pairs_per_sec << ','
              << checksum << ','
              << mismatches << ','
              << csv_clean(device) << ','
              << csv_clean(notes) << '\n';
}

int main(void) {
    @autoreleasepool {
        std::cout << "tool,backend,status,workload,n_reads,n_targets,len,k,error_rate,prep_seconds,seconds,total_seconds,reads_per_sec,total_reads_per_sec,pairs_per_sec,checksum,mismatches,device,notes\n";

        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (device == nil) {
            CaseSpec spec{0, 0, 0, 1, 0};
            print_row("dotmatch_gpu_metal", "metal", "unavailable", spec, 0.0, 0.0, 0, 0,
                      "no_metal_device", "MTLCreateSystemDefaultDevice returned nil");
            return 0;
        }

        NSError *error = nil;
        id<MTLLibrary> library = [device newLibraryWithSource:kernel_source() options:nil error:&error];
        if (library == nil) {
            CaseSpec spec{0, 0, 0, 1, 0};
            print_row("dotmatch_gpu_metal", "metal", "unavailable", spec, 0.0, 0.0, 0, 0,
                      [[device name] UTF8String],
                      error == nil ? "Metal library compile failed" : [[error localizedDescription] UTF8String]);
            return 0;
        }
        id<MTLFunction> function = [library newFunctionWithName:@"hamming_assign"];
        id<MTLComputePipelineState> pipeline = [device newComputePipelineStateWithFunction:function error:&error];
        if (pipeline == nil) {
            CaseSpec spec{0, 0, 0, 1, 0};
            print_row("dotmatch_gpu_metal", "metal", "unavailable", spec, 0.0, 0.0, 0, 0,
                      [[device name] UTF8String],
                      error == nil ? "Metal pipeline compile failed" : [[error localizedDescription] UTF8String]);
            return 0;
        }
        id<MTLCommandQueue> queue = [device newCommandQueue];
        if (queue == nil) {
            CaseSpec spec{0, 0, 0, 1, 0};
            print_row("dotmatch_gpu_metal", "metal", "unavailable", spec, 0.0, 0.0, 0, 0,
                      [[device name] UTF8String], "Metal command queue unavailable");
            return 0;
        }

        const CaseSpec cases[] = {
            {20000, 737, 20, 1, 10},
            {50000, 4096, 20, 1, 10},
        };

        for (const CaseSpec &spec : cases) {
            std::vector<std::string> target_strings(spec.n_targets);
            std::vector<std::string> read_strings(spec.n_reads);
            std::vector<uint64_t> target_codes(spec.n_targets);
            std::vector<uint64_t> read_codes(spec.n_reads);
            std::vector<const char *> target_ptrs(spec.n_targets);
            std::vector<const char *> read_ptrs(spec.n_reads);
            std::vector<size_t> target_lens(spec.n_targets, spec.len);
            std::vector<size_t> read_lens(spec.n_reads, spec.len);

            for (size_t i = 0; i < spec.n_targets; ++i) {
                rand_seq(target_strings[i], spec.len);
                target_codes[i] = pack_dna2(target_strings[i]);
                target_ptrs[i] = target_strings[i].c_str();
            }
            for (size_t i = 0; i < spec.n_reads; ++i) {
                size_t target_idx = (size_t)(xorshift64() % spec.n_targets);
                mutate_seq(target_strings[target_idx], read_strings[i], spec.err_per_thousand);
                read_codes[i] = pack_dna2(read_strings[i]);
                read_ptrs[i] = read_strings[i].c_str();
            }

            double cpu_prep_start = seconds_now();
            qdaln_index *index = qdaln_index_build(target_ptrs.data(), target_lens.data(), spec.n_targets);
            double cpu_prep_seconds = seconds_now() - cpu_prep_start;
            if (index == NULL) return 2;

            std::vector<qdaln_match_result> cpu_results(spec.n_reads);
            qdaln_index_stats stats;
            double cpu_start = seconds_now();
            int cpu_rc = qdaln_index_assign_hamming_stats(index, read_ptrs.data(), read_lens.data(), spec.n_reads,
                                                          spec.k, cpu_results.data(), &stats);
            double cpu_seconds = seconds_now() - cpu_start;
            if (cpu_rc != 0) return 2;
            long cpu_sum = checksum_cpu(cpu_results.data(), spec.n_reads);
            print_row("dotmatch_cpu_index", "cpu", "ok", spec, cpu_prep_seconds, cpu_seconds, cpu_sum, 0,
                      [[device name] UTF8String], "CPU indexed Hamming k=1 baseline; prep is index build");

            double gpu_prep_start = seconds_now();
            id<MTLBuffer> read_buffer = [device newBufferWithLength:read_codes.size() * sizeof(uint64_t)
                                                            options:MTLResourceStorageModeShared];
            id<MTLBuffer> target_buffer = [device newBufferWithLength:target_codes.size() * sizeof(uint64_t)
                                                              options:MTLResourceStorageModeShared];
            id<MTLBuffer> result_buffer = [device newBufferWithLength:spec.n_reads * sizeof(GpuResult)
                                                              options:MTLResourceStorageModeShared];
            id<MTLBuffer> params_buffer = [device newBufferWithLength:sizeof(GpuParams)
                                                              options:MTLResourceStorageModeShared];
            if (read_buffer == nil || target_buffer == nil || result_buffer == nil || params_buffer == nil) {
                qdaln_index_free(index);
                return 2;
            }
            memcpy([read_buffer contents], read_codes.data(), read_codes.size() * sizeof(uint64_t));
            memcpy([target_buffer contents], target_codes.data(), target_codes.size() * sizeof(uint64_t));
            GpuParams params{(uint32_t)spec.n_targets, (uint32_t)spec.len, (uint32_t)spec.k, 0};
            memcpy([params_buffer contents], &params, sizeof(params));
            double gpu_prep_seconds = seconds_now() - gpu_prep_start;

            id<MTLCommandBuffer> command_buffer = [queue commandBuffer];
            id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoder];
            [encoder setComputePipelineState:pipeline];
            [encoder setBuffer:read_buffer offset:0 atIndex:0];
            [encoder setBuffer:target_buffer offset:0 atIndex:1];
            [encoder setBuffer:params_buffer offset:0 atIndex:2];
            [encoder setBuffer:result_buffer offset:0 atIndex:3];
            NSUInteger threads_per_group = pipeline.maxTotalThreadsPerThreadgroup;
            if (threads_per_group > 256) threads_per_group = 256;
            MTLSize threads = MTLSizeMake(spec.n_reads, 1, 1);
            MTLSize group = MTLSizeMake(threads_per_group, 1, 1);
            double gpu_start = seconds_now();
            [encoder dispatchThreads:threads threadsPerThreadgroup:group];
            [encoder endEncoding];
            [command_buffer commit];
            [command_buffer waitUntilCompleted];
            double gpu_seconds = seconds_now() - gpu_start;

            GpuResult *gpu_results = (GpuResult *)[result_buffer contents];
            size_t mismatches = mismatch_count(cpu_results.data(), gpu_results, spec.n_reads);
            long gpu_sum = checksum_gpu(gpu_results, spec.n_reads);
            print_row("dotmatch_gpu_metal", "metal", "ok", spec, gpu_prep_seconds, gpu_seconds, gpu_sum, mismatches,
                      [[device name] UTF8String],
                      "Metal brute-force packed Hamming k=1; prep is shared-buffer allocation and copy");

            qdaln_index_free(index);
        }
    }
    return 0;
}
