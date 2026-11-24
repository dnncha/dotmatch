#define main dotmatch_cli_main
#include "../src/qda.c"
#undef main

#include <unistd.h>

static char libfuzz_path[64];

static int libfuzz_write(const uint8_t *data, size_t size) {
    if (libfuzz_path[0] == '\0') {
        snprintf(libfuzz_path, sizeof(libfuzz_path), "/tmp/dotmatch-libfuzz-%ld.fastq", (long)getpid());
    }
    FILE *f = fopen(libfuzz_path, "wb");
    if (f == NULL) return -1;
    size_t wrote = size == 0 ? 0 : fwrite(data, 1, size, f);
    return fclose(f) == 0 && wrote == size ? 0 : -1;
}

typedef struct libfuzz_trace {
    int rc;
    size_t records;
    uint64_t digest;
} libfuzz_trace;

static uint64_t libfuzz_mix(uint64_t h, const char *s, size_t n) {
    for (size_t i = 0; i < n; ++i) h = (h ^ (unsigned char)s[i]) * 0x100000001b3ULL;
    return (h ^ n) * 0x100000001b3ULL;
}

static libfuzz_trace libfuzz_read(size_t buffer_bytes, int view) {
    libfuzz_trace trace = {0, 0, 0xcbf29ce484222325ULL};
    fastq_reader reader;
    fastq_reader_buffer_bytes = buffer_bytes;
    int opened = fastq_reader_open(&reader, libfuzz_path);
    fastq_reader_buffer_bytes = (size_t)1 << 20;
    if (opened != 0) abort();
    char seq[128];
    for (;;) {
        const char *seq_out = seq;
        size_t seq_len = 0;
        int got = view ? fastq_read_sequence_view(&reader, seq, sizeof(seq), &seq_out, &seq_len)
                       : fastq_read_sequence_record_len(&reader, seq, sizeof(seq), &seq_len);
        if (got != 1) {
            trace.rc = got;
            break;
        }
        ++trace.records;
        trace.digest = libfuzz_mix(trace.digest, seq_out, seq_len);
    }
    fastq_reader_close(&reader);
    return trace;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 1) return 0;
    size_t small = 1 + data[0] % 64;
    if (libfuzz_write(data + 1, size - 1) != 0) abort();
    libfuzz_trace refill = libfuzz_read(small, 0);
    libfuzz_trace whole_view = libfuzz_read((size_t)1 << 20, 1);
    libfuzz_trace small_view = libfuzz_read(small, 1);
    if (refill.rc != whole_view.rc || refill.records != whole_view.records || refill.digest != whole_view.digest ||
        refill.rc != small_view.rc || refill.records != small_view.records || refill.digest != small_view.digest) {
        abort();
    }
    return 0;
}
