#include "qdmetal.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <algorithm>
#include <cstring>
#include <unordered_map>
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

struct MetalContext {
    id<MTLDevice> device;
    id<MTLComputePipelineState> brute_pipeline;
    id<MTLComputePipelineState> seed_pipeline;
    id<MTLCommandQueue> queue;
    char device_name[256];
};

static MetalContext g_ctx;
static int g_ctx_ready = 0;

static uint64_t seed_code_local(uint64_t packed, size_t offset, size_t length) {
    if (length == 0) return 0;
    uint64_t mask = length == 32 ? UINT64_MAX : ((1ULL << (2 * length)) - 1ULL);
    return (packed >> (2 * offset)) & mask;
}

static void build_seed_candidates(const std::vector<uint64_t> &reads, const std::vector<uint64_t> &targets, size_t len,
                                  int k, std::vector<uint32_t> *offsets, std::vector<uint32_t> *indices,
                                  size_t *candidate_count) {
    offsets->clear();
    indices->clear();
    *candidate_count = 0;
    offsets->reserve(reads.size() + 1);
    offsets->push_back(0);

    if (k > 1 || len == 0 || len > 32 || targets.empty()) {
        offsets->resize(reads.size() + 1, 0);
        return;
    }

    const size_t first_len = len / 2;
    const size_t second_offset = first_len;
    const size_t second_len = len - first_len;
    std::unordered_map<uint64_t, std::vector<uint32_t> > first_seed;
    std::unordered_map<uint64_t, std::vector<uint32_t> > second_seed;
    first_seed.reserve(targets.size() * 2);
    second_seed.reserve(targets.size() * 2);
    for (size_t i = 0; i < targets.size(); ++i) {
        first_seed[seed_code_local(targets[i], 0, first_len)].push_back((uint32_t)i);
        second_seed[seed_code_local(targets[i], second_offset, second_len)].push_back((uint32_t)i);
    }

    std::vector<uint32_t> seen(targets.size(), 0);
    uint32_t stamp = 1;
    for (uint64_t read : reads) {
        uint32_t before = (uint32_t)indices->size();
        const uint64_t s0 = seed_code_local(read, 0, first_len);
        const uint64_t s1 = seed_code_local(read, second_offset, second_len);
        std::unordered_map<uint64_t, std::vector<uint32_t> >::const_iterator it0 = first_seed.find(s0);
        if (it0 != first_seed.end()) {
            for (uint32_t target_idx : it0->second) {
                if (seen[target_idx] == stamp) continue;
                seen[target_idx] = stamp;
                indices->push_back(target_idx);
            }
        }
        std::unordered_map<uint64_t, std::vector<uint32_t> >::const_iterator it1 = second_seed.find(s1);
        if (it1 != second_seed.end()) {
            for (uint32_t target_idx : it1->second) {
                if (seen[target_idx] == stamp) continue;
                seen[target_idx] = stamp;
                indices->push_back(target_idx);
            }
        }
        *candidate_count += (size_t)((uint32_t)indices->size() - before);
        offsets->push_back((uint32_t)indices->size());
        ++stamp;
        if (stamp == 0) {
            std::fill(seen.begin(), seen.end(), 0);
            stamp = 1;
        }
    }
}

static NSString *kernel_source(void) {
    return @"#include <metal_stdlib>\n"
            "using namespace metal;\n"
            "struct GpuResult { int target_index; int best_distance; int second_best_distance; int match_count; int status; };\n"
            "struct GpuParams { uint n_targets; uint len; uint k; uint reserved; };\n"
            "static inline GpuResult assign_from_candidates(ulong read, device const ulong *targets, constant GpuParams &params, uint begin, uint end, device const uint *candidate_indices, bool use_candidates) {\n"
            "    ulong mask = params.len == 32 ? 0xffffffffffffffffUL : ((1UL << (params.len * 2)) - 1UL);\n"
            "    int best = -1;\n"
            "    int second = -1;\n"
            "    int best_idx = -1;\n"
            "    int match_count = 0;\n"
            "    int best_ties = 0;\n"
            "    for (uint p = begin; p < end; ++p) {\n"
            "        uint j = use_candidates ? candidate_indices[p] : p;\n"
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
            "    return out;\n"
            "}\n"
            "kernel void hamming_assign(device const ulong *reads [[buffer(0)]],\n"
            "                           device const ulong *targets [[buffer(1)]],\n"
            "                           constant GpuParams &params [[buffer(2)]],\n"
            "                           device GpuResult *results [[buffer(3)]],\n"
            "                           uint gid [[thread_position_in_grid]]) {\n"
            "    results[gid] = assign_from_candidates(reads[gid], targets, params, 0, params.n_targets, nullptr, false);\n"
            "}\n"
            "kernel void hamming_seed_assign(device const ulong *reads [[buffer(0)]],\n"
            "                                device const ulong *targets [[buffer(1)]],\n"
            "                                constant GpuParams &params [[buffer(2)]],\n"
            "                                device GpuResult *results [[buffer(3)]],\n"
            "                                device const uint *candidate_offsets [[buffer(4)]],\n"
            "                                device const uint *candidate_indices [[buffer(5)]],\n"
            "                                uint gid [[thread_position_in_grid]]) {\n"
            "    uint begin = candidate_offsets[gid];\n"
            "    uint end = candidate_offsets[gid + 1];\n"
            "    results[gid] = assign_from_candidates(reads[gid], targets, params, begin, end, candidate_indices, true);\n"
            "}\n";
}

static int ensure_metal_context(void) {
    if (g_ctx_ready) return 1;
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (device == nil) return 0;
        NSError *error = nil;
        id<MTLLibrary> library = [device newLibraryWithSource:kernel_source() options:nil error:&error];
        if (library == nil) return 0;
        id<MTLFunction> brute_fn = [library newFunctionWithName:@"hamming_assign"];
        id<MTLFunction> seed_fn = [library newFunctionWithName:@"hamming_seed_assign"];
        id<MTLComputePipelineState> brute_pipeline = [device newComputePipelineStateWithFunction:brute_fn error:&error];
        id<MTLComputePipelineState> seed_pipeline = [device newComputePipelineStateWithFunction:seed_fn error:&error];
        id<MTLCommandQueue> queue = [device newCommandQueue];
        if (brute_pipeline == nil || seed_pipeline == nil || queue == nil) return 0;
        g_ctx.device = device;
        g_ctx.brute_pipeline = brute_pipeline;
        g_ctx.seed_pipeline = seed_pipeline;
        g_ctx.queue = queue;
        strncpy(g_ctx.device_name, [[device name] UTF8String], sizeof(g_ctx.device_name) - 1);
        g_ctx.device_name[sizeof(g_ctx.device_name) - 1] = '\0';
        g_ctx_ready = 1;
        return 1;
    }
}

static int dispatch_brute(id<MTLBuffer> read_buffer, id<MTLBuffer> target_buffer, id<MTLBuffer> params_buffer,
                          id<MTLBuffer> result_buffer, size_t n_reads) {
    id<MTLCommandBuffer> command_buffer = [g_ctx.queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoder];
    [encoder setComputePipelineState:g_ctx.brute_pipeline];
    [encoder setBuffer:read_buffer offset:0 atIndex:0];
    [encoder setBuffer:target_buffer offset:0 atIndex:1];
    [encoder setBuffer:params_buffer offset:0 atIndex:2];
    [encoder setBuffer:result_buffer offset:0 atIndex:3];
    NSUInteger threads_per_group = g_ctx.brute_pipeline.maxTotalThreadsPerThreadgroup;
    if (threads_per_group > 256) threads_per_group = 256;
    [encoder dispatchThreads:MTLSizeMake(n_reads, 1, 1) threadsPerThreadgroup:MTLSizeMake(threads_per_group, 1, 1)];
    [encoder endEncoding];
    [command_buffer commit];
    [command_buffer waitUntilCompleted];
    return command_buffer.status == MTLCommandBufferStatusCompleted ? 0 : -1;
}

static int dispatch_seed(id<MTLBuffer> read_buffer, id<MTLBuffer> target_buffer, id<MTLBuffer> params_buffer,
                         id<MTLBuffer> result_buffer, id<MTLBuffer> candidate_offsets_buffer,
                         id<MTLBuffer> candidate_indices_buffer, size_t n_reads) {
    id<MTLCommandBuffer> command_buffer = [g_ctx.queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoder];
    [encoder setComputePipelineState:g_ctx.seed_pipeline];
    [encoder setBuffer:read_buffer offset:0 atIndex:0];
    [encoder setBuffer:target_buffer offset:0 atIndex:1];
    [encoder setBuffer:params_buffer offset:0 atIndex:2];
    [encoder setBuffer:result_buffer offset:0 atIndex:3];
    [encoder setBuffer:candidate_offsets_buffer offset:0 atIndex:4];
    [encoder setBuffer:candidate_indices_buffer offset:0 atIndex:5];
    NSUInteger threads_per_group = g_ctx.seed_pipeline.maxTotalThreadsPerThreadgroup;
    if (threads_per_group > 256) threads_per_group = 256;
    [encoder dispatchThreads:MTLSizeMake(n_reads, 1, 1) threadsPerThreadgroup:MTLSizeMake(threads_per_group, 1, 1)];
    [encoder endEncoding];
    [command_buffer commit];
    [command_buffer waitUntilCompleted];
    return command_buffer.status == MTLCommandBufferStatusCompleted ? 0 : -1;
}

int qdmetal_available(void) {
    return ensure_metal_context();
}

const char *qdmetal_device_name(void) {
    if (!ensure_metal_context()) return NULL;
    return g_ctx.device_name;
}

int qdmetal_hamming_assign(const uint64_t *read_codes, size_t n_reads, const uint64_t *target_codes, size_t n_targets,
                           size_t len, int k, qdmetal_match_result *results, qdmetal_assign_stats *stats) {
    if (read_codes == NULL || target_codes == NULL || results == NULL || n_reads == 0 || n_targets == 0 || len == 0 ||
        len > 32 || (k != 0 && k != 1)) {
        return -1;
    }
    if (!ensure_metal_context()) return -1;

    if (stats != NULL) {
        stats->candidates_considered = 0;
        stats->candidates_verified = 0;
        stats->path = QDMETAL_PATH_BRUTE_FORCE;
        stats->device_name = g_ctx.device_name;
    }

    @autoreleasepool {
        id<MTLDevice> device = g_ctx.device;
        id<MTLBuffer> read_buffer =
                [device newBufferWithBytes:read_codes length:n_reads * sizeof(uint64_t) options:MTLResourceStorageModeShared];
        id<MTLBuffer> target_buffer = [device newBufferWithBytes:target_codes length:n_targets * sizeof(uint64_t)
                                                       options:MTLResourceStorageModeShared];
        id<MTLBuffer> result_buffer =
                [device newBufferWithLength:n_reads * sizeof(GpuResult) options:MTLResourceStorageModeShared];
        id<MTLBuffer> params_buffer = [device newBufferWithLength:sizeof(GpuParams) options:MTLResourceStorageModeShared];
        if (read_buffer == nil || target_buffer == nil || result_buffer == nil || params_buffer == nil) return -1;

        GpuParams params{(uint32_t)n_targets, (uint32_t)len, (uint32_t)k, 0};
        memcpy([params_buffer contents], &params, sizeof(params));

        const size_t seed_threshold = 1024;
        int rc = 0;
        if (n_targets >= seed_threshold && k <= 1) {
            std::vector<uint64_t> reads(read_codes, read_codes + n_reads);
            std::vector<uint64_t> targets(target_codes, target_codes + n_targets);
            std::vector<uint32_t> candidate_offsets;
            std::vector<uint32_t> candidate_indices;
            size_t candidate_count = 0;
            build_seed_candidates(reads, targets, len, k, &candidate_offsets, &candidate_indices, &candidate_count);
            id<MTLBuffer> candidate_offsets_buffer =
                    [device newBufferWithBytes:candidate_offsets.data()
                                        length:candidate_offsets.size() * sizeof(uint32_t)
                                       options:MTLResourceStorageModeShared];
            size_t candidate_bytes = std::max<size_t>(candidate_indices.size() * sizeof(uint32_t), sizeof(uint32_t));
            id<MTLBuffer> candidate_indices_buffer =
                    [device newBufferWithLength:candidate_bytes options:MTLResourceStorageModeShared];
            if (candidate_offsets_buffer == nil || candidate_indices_buffer == nil) return -1;
            if (!candidate_indices.empty()) {
                memcpy([candidate_indices_buffer contents], candidate_indices.data(),
                       candidate_indices.size() * sizeof(uint32_t));
            }
            rc = dispatch_seed(read_buffer, target_buffer, params_buffer, result_buffer, candidate_offsets_buffer,
                               candidate_indices_buffer, n_reads);
            if (stats != NULL) {
                stats->path = QDMETAL_PATH_SEED_INDEX;
                stats->candidates_considered = candidate_count;
                stats->candidates_verified = candidate_count;
            }
        } else {
            rc = dispatch_brute(read_buffer, target_buffer, params_buffer, result_buffer, n_reads);
            if (stats != NULL) {
                stats->candidates_considered = n_reads * n_targets;
                stats->candidates_verified = n_reads * n_targets;
            }
        }
        if (rc != 0) return -1;

        const GpuResult *gpu_results = (const GpuResult *)[result_buffer contents];
        for (size_t i = 0; i < n_reads; ++i) {
            results[i].target_index = gpu_results[i].target_index;
            results[i].best_distance = gpu_results[i].best_distance;
            results[i].second_best_distance = gpu_results[i].second_best_distance;
            results[i].match_count = gpu_results[i].match_count;
            results[i].status = gpu_results[i].status;
        }
    }
    return 0;
}