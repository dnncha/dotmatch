CC ?= cc
CXX ?= c++
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -Wpedantic -Iinclude
CXXFLAGS ?= -O3 -std=c++11 -Wall -Wextra -Iinclude
LDFLAGS ?=
ZLIB_LIBS ?= -lz
PTHREAD_LIBS ?= -pthread
COVERAGE_CC ?= clang
COVERAGE_MIN ?= 75
LLVM_PROFDATA ?= $(shell command -v llvm-profdata 2>/dev/null || xcrun --find llvm-profdata 2>/dev/null)
LLVM_COV ?= $(shell command -v llvm-cov 2>/dev/null || xcrun --find llvm-cov 2>/dev/null)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
SHARED_EXT := dylib
DOTMATCH_SHARED_FLAGS := -dynamiclib -install_name @rpath/libdotmatch.dylib
QDALIGN_SHARED_FLAGS := -dynamiclib -install_name @rpath/libqdalign.dylib
else
SHARED_EXT := so
DOTMATCH_SHARED_FLAGS := -shared
QDALIGN_SHARED_FLAGS := -shared
endif

.PHONY: all clean test cli-test coverage bench bench-batch bench-small bench-paper figures bench-real-report bench-barcode-demux bench-barcode-demux-competitors bench-barcode-sota barcode-sota-report barcode-sota-gate barcode-demux-report barcode-competitor-env fetch-barcode-demo fetch-barcode-demo-claim fetch-sanson-crispr fetch-10x-bcl-demo bench-bcl-small bench-bcl-10x bench-bcl-real bench-bcl-real-repeated bcl-figures bcl-competitor-env bcl-linux-env bcl-sota-gate bench-public-crispr-small bench-public-crispr bench-public-crispr-competitors bench-public-crispr-repeated bench-public-crispr-scaling bench-crispr-sota crispr-sota-report crispr-sota-gate atlas-crispr-sota atlas-barcode-sota count-agreement count-agreement-sota validate-public-crispr-edlib validate-crispr-sota-edlib public-crispr-report public-crispr-claim-gate public-crispr-smoke-gate competitor-env edlib edlib-tools bench-edlib-native benchmark-report benchmark-report-native asan shared python-test python-package-test publication-ready

all: dotmatch libdotmatch.a qda libqdalign.a

build:
	mkdir -p build

build/qdalign.o: src/qdalign.c include/qdalign.h | build
	$(CC) $(CFLAGS) -c src/qdalign.c -o $@

libdotmatch.a: build/qdalign.o
	ar rcs $@ $^

libqdalign.a: build/qdalign.o
	ar rcs $@ $^

build/qdalign.pic.o: src/qdalign.c include/qdalign.h | build
	$(CC) $(CFLAGS) -fPIC -c src/qdalign.c -o $@

libdotmatch.$(SHARED_EXT): build/qdalign.pic.o
	$(CC) $(DOTMATCH_SHARED_FLAGS) $^ -o $@ $(LDFLAGS)

libqdalign.$(SHARED_EXT): build/qdalign.pic.o
	$(CC) $(QDALIGN_SHARED_FLAGS) $^ -o $@ $(LDFLAGS)

shared: libdotmatch.$(SHARED_EXT) libqdalign.$(SHARED_EXT)

dotmatch: src/qda.c build/qdalign.o include/qdalign.h
	$(CC) $(CFLAGS) src/qda.c build/qdalign.o -o $@ $(LDFLAGS) $(ZLIB_LIBS) $(PTHREAD_LIBS)

qda: src/qda.c build/qdalign.o include/qdalign.h
	$(CC) $(CFLAGS) src/qda.c build/qdalign.o -o $@ $(LDFLAGS) $(ZLIB_LIBS) $(PTHREAD_LIBS)

build/test_qdalign: tests/test_qdalign.c build/qdalign.o include/qdalign.h | build
	$(CC) $(CFLAGS) tests/test_qdalign.c build/qdalign.o -o $@ $(LDFLAGS)

test: build/test_qdalign
	./build/test_qdalign

cli-test: dotmatch
	sh tests/test_cli_fastq.sh

coverage:
	test -n "$(LLVM_PROFDATA)"
	test -n "$(LLVM_COV)"
	rm -rf build/coverage
	mkdir -p build/coverage
	$(COVERAGE_CC) -O0 -g -std=c11 -Wall -Wextra -Wpedantic -Iinclude -fprofile-instr-generate -fcoverage-mapping -c src/qdalign.c -o build/coverage/qdalign.o
	$(COVERAGE_CC) -O0 -g -std=c11 -Wall -Wextra -Wpedantic -Iinclude -fprofile-instr-generate -fcoverage-mapping tests/test_qdalign.c build/coverage/qdalign.o -o build/coverage/test_qdalign
	LLVM_PROFILE_FILE=build/coverage/test-%p.profraw ./build/coverage/test_qdalign
	$(COVERAGE_CC) -O0 -g -std=c11 -Wall -Wextra -Wpedantic -Iinclude -fprofile-instr-generate -fcoverage-mapping src/qda.c src/qdalign.c -o build/coverage/dotmatch $(ZLIB_LIBS) $(PTHREAD_LIBS)
	LLVM_PROFILE_FILE=build/coverage/cli-%p.profraw DOTMATCH_BIN="$(CURDIR)/build/coverage/dotmatch" sh tests/test_cli_fastq.sh
	$(LLVM_PROFDATA) merge -sparse build/coverage/*.profraw -o build/coverage/coverage.profdata
	$(LLVM_COV) report build/coverage/test_qdalign -object build/coverage/dotmatch -instr-profile=build/coverage/coverage.profdata --sources src/qdalign.c --sources src/qda.c --show-branch-summary | tee build/coverage/coverage.txt
	$(LLVM_COV) export build/coverage/test_qdalign -object build/coverage/dotmatch -instr-profile=build/coverage/coverage.profdata --sources src/qdalign.c --sources src/qda.c > build/coverage/coverage.json
	$(LLVM_COV) show build/coverage/test_qdalign -object build/coverage/dotmatch -instr-profile=build/coverage/coverage.profdata --sources src/qdalign.c --sources src/qda.c -format=html -output-dir=build/coverage/html
	python3 scripts/check_coverage_threshold.py build/coverage/coverage.json --min-lines "$(COVERAGE_MIN)"

build/bench: tools/bench.c build/qdalign.o include/qdalign.h | build
	$(CC) $(CFLAGS) tools/bench.c build/qdalign.o -o $@ $(LDFLAGS)

bench: build/bench
	./build/bench

build/bench_batch: tools/bench_batch.c build/qdalign.o include/qdalign.h | build
	$(CC) $(CFLAGS) tools/bench_batch.c build/qdalign.o -o $@ $(LDFLAGS)

bench-batch: build/bench_batch
	./build/bench_batch

edlib:
	python3 scripts/fetch_edlib.py

build/bench_edlib_native: tools/bench_edlib_native.cpp build/qdalign.o include/qdalign.h | build edlib
	$(CXX) $(CXXFLAGS) -Ibuild/edlib/edlib/include tools/bench_edlib_native.cpp build/qdalign.o build/edlib/edlib/src/edlib.cpp -o $@ $(LDFLAGS)

build/dotmatch_edlib_validate: tools/dotmatch_edlib_validate.cpp build/qdalign.o include/qdalign.h | build edlib
	$(CXX) $(CXXFLAGS) -Ibuild/edlib/edlib/include tools/dotmatch_edlib_validate.cpp build/qdalign.o build/edlib/edlib/src/edlib.cpp -o $@ $(LDFLAGS) $(ZLIB_LIBS) $(PTHREAD_LIBS)

build/bench_real_edlib: tools/bench_real_edlib.cpp build/qdalign.o include/qdalign.h | build edlib
	$(CXX) $(CXXFLAGS) -Ibuild/edlib/edlib/include tools/bench_real_edlib.cpp build/qdalign.o build/edlib/edlib/src/edlib.cpp -o $@ $(LDFLAGS) $(ZLIB_LIBS)

edlib-tools: build/dotmatch_edlib_validate

bench-edlib-native: build/bench_edlib_native
	./build/bench_edlib_native

bench-small: build/bench_edlib_native
	./build/bench_edlib_native 100

bench-paper: build/bench_edlib_native
	DOTMATCH_NATIVE_MATRIX=paper DOTMATCH_NATIVE_REPORT_READS=$${DOTMATCH_NATIVE_REPORT_READS:-1000} python3 scripts/generate_native_benchmark_report.py

figures: benchmark-report-native

bench-real-report: dotmatch build/bench_real_edlib
	python3 scripts/generate_real_benchmark_report.py

bench-barcode-demux: dotmatch
	python3 scripts/bench_barcode_demux.py --run-hash-splitter
	python3 scripts/generate_barcode_demux_report.py

barcode-competitor-env:
	sh scripts/install_barcode_competitors.sh

bench-barcode-demux-competitors: dotmatch barcode-competitor-env
	PATH="$(CURDIR)/build/barcode-competitors/bin:$$PATH" python3 scripts/bench_barcode_demux.py --run-cutadapt --run-hash-splitter
	python3 scripts/generate_barcode_demux_report.py

fetch-barcode-demo:
	python3 scripts/fetch_srp009896_barcode_demo.py --metadata-only

fetch-barcode-demo-claim:
	test -n "$$DOTMATCH_BARCODE_SOTA_BARCODES" || test -n "$$DOTMATCH_BARCODE_SOTA_BARCODES_URL"
	python3 scripts/fetch_srp009896_barcode_demo.py --require-barcodes --subsample "$${DOTMATCH_BARCODE_SOTA_SUBSAMPLE:-100000}" $${DOTMATCH_BARCODE_SOTA_BARCODES:+--barcodes-file "$$DOTMATCH_BARCODE_SOTA_BARCODES"} $${DOTMATCH_BARCODE_SOTA_BARCODES_URL:+--barcodes-url "$$DOTMATCH_BARCODE_SOTA_BARCODES_URL"} --barcode-start "$${DOTMATCH_BARCODE_START:-0}" $${DOTMATCH_BARCODE_LENGTH:+--barcode-length "$$DOTMATCH_BARCODE_LENGTH"}

bench-barcode-sota: dotmatch barcode-competitor-env fetch-barcode-demo-claim
	PATH="$(CURDIR)/build/barcode-competitors/bin:$$PATH" python3 scripts/bench_barcode_demux.py --reads "$$(python3 -c 'import json;print(json.load(open("examples/barcode_demux/data/metadata.json"))["runs"][0]["local_fastq"])')" --barcodes "$$(python3 -c 'import json;print(json.load(open("examples/barcode_demux/data/metadata.json"))["barcodes"])')" --barcode-start "$${DOTMATCH_BARCODE_START:-0}" --barcode-length "$$(python3 -c 'import json;print(json.load(open("examples/barcode_demux/data/metadata.json")).get("barcode_length") or 8)')" --k "$${DOTMATCH_BARCODE_K:-1}" --metric "$${DOTMATCH_BARCODE_METRIC:-hamming}" --workflow-name real_srp009896_inline_barcode --run-cutadapt --run-hash-splitter --repeats "$${DOTMATCH_BARCODE_REPEATS:-5}"
	python3 scripts/generate_barcode_demux_report.py

barcode-sota-report:
	python3 scripts/generate_barcode_demux_report.py

barcode-sota-gate:
	python3 scripts/check_barcode_sota_gate.py

fetch-sanson-crispr:
	python3 scripts/fetch_sanson_brunello_demo.py --subsample "$${DOTMATCH_SANSON_SUBSAMPLE:-100000}"

fetch-10x-bcl-demo:
	python3 scripts/fetch_10x_tiny_bcl.py --extract

barcode-demux-report:
	python3 scripts/generate_barcode_demux_report.py

bench-bcl-small: dotmatch
	python3 scripts/bench_bcl_demux.py
	python3 scripts/generate_bcl_demux_report.py

bench-bcl-10x: dotmatch fetch-10x-bcl-demo
	python3 scripts/bench_bcl_demux.py --run-folder examples/bcl_demux/data/cellranger-tiny-bcl-1.2.0 --sample-sheet examples/bcl_demux/data/cellranger-tiny-bcl-samplesheet.normalized.csv --workflow-name public_10x_tiny_bcl --detect-competitors --run-installed-competitors --threads "$${DOTMATCH_BCL_THREADS:-1}" --gzip-level "$${DOTMATCH_BCL_GZIP_LEVEL:-1}"
	python3 scripts/generate_bcl_demux_report.py

bench-bcl-real: dotmatch
	test -n "$$DOTMATCH_BCL_RUN_FOLDER"
	test -n "$$DOTMATCH_BCL_SAMPLE_SHEET"
	python3 scripts/bench_bcl_demux.py --run-folder "$$DOTMATCH_BCL_RUN_FOLDER" --sample-sheet "$$DOTMATCH_BCL_SAMPLE_SHEET" --workflow-name real_bcl_user_supplied --detect-competitors --run-installed-competitors --threads "$${DOTMATCH_BCL_THREADS:-1}" --gzip-level "$${DOTMATCH_BCL_GZIP_LEVEL:-1}"
	python3 scripts/generate_bcl_demux_report.py

bench-bcl-real-repeated: dotmatch
	test -n "$$DOTMATCH_BCL_RUN_FOLDER"
	test -n "$$DOTMATCH_BCL_SAMPLE_SHEET"
	python3 scripts/run_bcl_repeated.py --run-folder "$$DOTMATCH_BCL_RUN_FOLDER" --sample-sheet "$$DOTMATCH_BCL_SAMPLE_SHEET" --workflow-name real_bcl_user_supplied --detect-competitors --run-installed-competitors --threads "$${DOTMATCH_BCL_THREADS:-1}" --gzip-level "$${DOTMATCH_BCL_GZIP_LEVEL:-1}" --repeats "$${DOTMATCH_BCL_REPEATS:-5}"

bcl-figures:
	python3 scripts/generate_bcl_demux_report.py

bcl-competitor-env:
	sh scripts/check_bcl_competitors.sh

bcl-linux-env:
	sh scripts/check_bcl_linux_env.sh

bcl-sota-gate:
	python3 scripts/check_bcl_sota_gate.py

bench-public-crispr-small: dotmatch
	python3 scripts/run_public_crispr_benchmark.py --small

bench-public-crispr: dotmatch
	python3 scripts/run_public_crispr_benchmark.py

competitor-env:
	sh scripts/install_competitors.sh

bench-public-crispr-competitors: dotmatch competitor-env
	PATH="$(CURDIR)/build/guide-counter/bin:$(CURDIR)/build/competitor-env/bin:$$PATH" python3 scripts/run_public_crispr_benchmark.py --small --run-mageck --run-cutadapt --run-bowtie2 --run-guide-counter

bench-public-crispr-repeated: dotmatch competitor-env
	PATH="$(CURDIR)/build/guide-counter/bin:$(CURDIR)/build/competitor-env/bin:$$PATH" python3 scripts/run_public_crispr_repeated.py --run-mageck --run-guide-counter

bench-public-crispr-scaling: dotmatch competitor-env
	PATH="$(CURDIR)/build/guide-counter/bin:$(CURDIR)/build/competitor-env/bin:$$PATH" python3 scripts/bench_public_crispr_sample_scaling.py --run-guide-counter

count-agreement:
	python3 scripts/compare_count_tables.py

validate-public-crispr-edlib: dotmatch edlib-tools
	python3 scripts/validate_public_crispr_edlib.py

validate-crispr-sota-edlib: dotmatch edlib-tools
	python3 scripts/validate_crispr_sota_edlib.py

bench-crispr-sota: dotmatch competitor-env
	PATH="$(CURDIR)/build/guide-counter/bin:$(CURDIR)/build/competitor-env/bin:$$PATH" python3 scripts/run_crispr_sota_repeated.py --run-mageck --run-guide-counter $${DOTMATCH_SOTA_FULL:+--full}

count-agreement-sota:
	python3 scripts/compare_crispr_sota_counts.py

crispr-sota-report:
	python3 scripts/generate_crispr_sota_report.py

crispr-sota-gate:
	python3 scripts/check_crispr_sota_gate.py

atlas-crispr-sota:
	sh scripts/atlas_run_crispr_sota.sh

atlas-barcode-sota:
	sh scripts/atlas_run_barcode_sota.sh

public-crispr-report:
	python3 scripts/generate_public_crispr_report.py

public-crispr-claim-gate:
	python3 scripts/check_public_crispr_claim_gate.py

public-crispr-smoke-gate:
	python3 scripts/check_public_crispr_claim_gate.py --smoke

benchmark-report: shared build/bench_batch
	python3 scripts/generate_benchmark_report.py

benchmark-report-native: build/bench_edlib_native
	python3 scripts/generate_native_benchmark_report.py

python-test: shared
	DOTMATCH_LIB=$(CURDIR)/libdotmatch.$(SHARED_EXT) PYTHONPATH=$(CURDIR)/python python3 -m pytest python/tests

python-package-test:
	python3 scripts/check_python_wheel.py

publication-ready:
	python3 scripts/check_publication_ready.py

asan:
	$(MAKE) clean
	$(MAKE) CFLAGS='-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -std=c11 -Wall -Wextra -Wpedantic -Iinclude' LDFLAGS='-fsanitize=address,undefined' test
	$(MAKE) clean

clean:
	rm -rf build dotmatch qda libdotmatch.a libdotmatch.so libdotmatch.dylib libqdalign.a libqdalign.so libqdalign.dylib
