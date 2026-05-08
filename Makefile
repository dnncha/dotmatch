CC ?= cc
CXX ?= c++
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -Wpedantic -Iinclude
CXXFLAGS ?= -O3 -std=c++11 -Wall -Wextra -Iinclude
LDFLAGS ?=
ZLIB_LIBS ?= -lz
PTHREAD_LIBS ?= -pthread
DOTMATCH_VERSION ?= $(shell python3 -c 'import pathlib; print(next(line.split("\"")[1] for line in pathlib.Path("pyproject.toml").read_text(encoding="utf-8").splitlines() if line.startswith("version = ")))')
DOTMATCH_VERSION_CFLAGS := -DDOTMATCH_VERSION=\"$(DOTMATCH_VERSION)\"
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

.PHONY: all clean test cli-test coverage bench bench-batch bench-small bench-native-matrix figures bench-real-report bench-barcode-demux bench-barcode-demux-competitors bench-barcode-comparison barcode-comparison-report barcode-comparison-gate barcode-demux-report barcode-competitor-env fetch-barcode-demo fetch-barcode-demo-claim fetch-sanson-crispr fetch-10x-bcl-demo bench-bcl-small bench-bcl-10x bench-bcl-real bench-bcl-real-repeated bcl-figures bcl-competitor-env bcl-linux-env bcl-tiny-public-gate bcl-comparison-gate fetch-oligo-adapter-demo bench-oligo-adapter bench-oligo-adapter-public oligo-adapter-smoke-gate oligo-adapter-public-gate fetch-amplicon-panel-demo bench-amplicon-panel bench-amplicon-panel-public amplicon-panel-smoke-gate amplicon-panel-public-gate fetch-feature-barcode-demo bench-feature-barcode bench-feature-barcode-public feature-barcode-smoke-gate feature-barcode-public-gate fetch-perturb-seq-demo bench-perturb-seq bench-perturb-seq-public perturb-seq-smoke-gate perturb-seq-public-gate bench-public-crispr-small bench-public-crispr bench-public-crispr-competitors bench-public-crispr-repeated bench-public-crispr-scaling bench-crispr-comparison crispr-comparison-report crispr-comparison-gate count-agreement count-agreement-comparison validate-public-crispr-edlib validate-crispr-comparison-edlib public-crispr-report public-crispr-evidence-gate public-crispr-smoke-gate competitor-env edlib edlib-tools bench-edlib-native benchmark-report benchmark-report-native asan shared python-test python-package-test repository-ready release-ready pretag-ready assay-evidence-ready alphabet-policy-ready citation-metadata-ready native-comparator-scope-ready workflow-examples-ready workflow-adoption-status distribution-record-ready distribution-submission-ready bioconda-recipe-ready distribution-channels

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

dotmatch: src/qda.c build/qdalign.o include/qdalign.h Makefile
	$(CC) $(CFLAGS) $(DOTMATCH_VERSION_CFLAGS) src/qda.c build/qdalign.o -o $@ $(LDFLAGS) $(ZLIB_LIBS) $(PTHREAD_LIBS)

qda: src/qda.c build/qdalign.o include/qdalign.h Makefile
	$(CC) $(CFLAGS) $(DOTMATCH_VERSION_CFLAGS) src/qda.c build/qdalign.o -o $@ $(LDFLAGS) $(ZLIB_LIBS) $(PTHREAD_LIBS)

build/test_qdalign: tests/test_qdalign.c build/qdalign.o include/qdalign.h | build
	$(CC) $(CFLAGS) tests/test_qdalign.c build/qdalign.o -o $@ $(LDFLAGS)

build/test_qdalign_threshold_alloc: tests/test_qdalign_threshold_alloc.c src/qdalign.c include/qdalign.h | build
	$(CC) $(CFLAGS) tests/test_qdalign_threshold_alloc.c -o $@ $(LDFLAGS)

test: build/test_qdalign build/test_qdalign_threshold_alloc
	./build/test_qdalign
	./build/test_qdalign_threshold_alloc

cli-test: dotmatch
	sh tests/test_cli_fastq.sh
	sh tests/test_crispr_example_expected.sh

coverage:
	test -n "$(LLVM_PROFDATA)"
	test -n "$(LLVM_COV)"
	rm -rf build/coverage
	mkdir -p build/coverage
	$(COVERAGE_CC) -O0 -g -std=c11 -Wall -Wextra -Wpedantic -Iinclude -fprofile-instr-generate -fcoverage-mapping -c src/qdalign.c -o build/coverage/qdalign.o
	$(COVERAGE_CC) -O0 -g -std=c11 -Wall -Wextra -Wpedantic -Iinclude -fprofile-instr-generate -fcoverage-mapping tests/test_qdalign.c build/coverage/qdalign.o -o build/coverage/test_qdalign
	LLVM_PROFILE_FILE=build/coverage/test-%p.profraw ./build/coverage/test_qdalign
	$(COVERAGE_CC) -O0 -g -std=c11 -Wall -Wextra -Wpedantic -Iinclude $(DOTMATCH_VERSION_CFLAGS) -fprofile-instr-generate -fcoverage-mapping src/qda.c src/qdalign.c -o build/coverage/dotmatch $(ZLIB_LIBS) $(PTHREAD_LIBS)
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

bench-native-matrix: build/bench_edlib_native
	DOTMATCH_NATIVE_MATRIX=full DOTMATCH_NATIVE_REPORT_READS=$${DOTMATCH_NATIVE_REPORT_READS:-1000} python3 scripts/generate_native_benchmark_report.py

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
	python3 scripts/fetch_srp009896_barcode_demo.py --metadata-only --use-public-example-barcodes

fetch-barcode-demo-claim:
	python3 scripts/fetch_srp009896_barcode_demo.py --require-barcodes --subsample "$${DOTMATCH_BARCODE_COMPARISON_SUBSAMPLE:-100000}" $${DOTMATCH_BARCODE_COMPARISON_BARCODES:+--barcodes-file "$$DOTMATCH_BARCODE_COMPARISON_BARCODES"} $${DOTMATCH_BARCODE_COMPARISON_BARCODES_URL:+--barcodes-url "$$DOTMATCH_BARCODE_COMPARISON_BARCODES_URL"} $${DOTMATCH_BARCODE_COMPARISON_USE_PUBLIC_EXAMPLE:+--use-public-example-barcodes} --barcode-start "$${DOTMATCH_BARCODE_START:-1}" $${DOTMATCH_BARCODE_LENGTH:+--barcode-length "$$DOTMATCH_BARCODE_LENGTH"}

bench-barcode-comparison: dotmatch barcode-competitor-env fetch-barcode-demo-claim
	PATH="$(CURDIR)/build/barcode-competitors/bin:$$PATH" python3 scripts/bench_barcode_demux.py --reads "$$(python3 -c 'import json;print(json.load(open("examples/barcode_demux/data/metadata.json"))["runs"][0]["local_fastq"])')" --barcodes "$$(python3 -c 'import json;print(json.load(open("examples/barcode_demux/data/metadata.json"))["barcodes"])')" --barcode-start "$${DOTMATCH_BARCODE_START:-1}" --barcode-length "$$(python3 -c 'import json;m=json.load(open("examples/barcode_demux/data/metadata.json"));print(m.get("barcode_length") or ("auto" if m.get("barcode_length_mode") == "auto" else 8))')" --k "$${DOTMATCH_BARCODE_K:-0}" --metric "$${DOTMATCH_BARCODE_METRIC:-hamming}" --workflow-name real_srp009896_inline_barcode --run-cutadapt --run-hash-splitter --repeats "$${DOTMATCH_BARCODE_REPEATS:-5}"
	python3 scripts/generate_barcode_demux_report.py

barcode-comparison-report:
	python3 scripts/generate_barcode_demux_report.py

barcode-comparison-gate:
	python3 scripts/check_barcode_comparison_gate.py

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

bcl-tiny-public-gate:
	python3 scripts/check_bcl_tiny_public_gate.py

bcl-comparison-gate:
	python3 scripts/check_bcl_comparison_gate.py

bench-oligo-adapter: dotmatch
	python3 scripts/bench_oligo_adapter.py

fetch-oligo-adapter-demo:
	python3 scripts/fetch_public_oligo_adapter_demo.py --records "$${DOTMATCH_OLIGO_ADAPTER_SUBSAMPLE:-10000}"

bench-oligo-adapter-public: dotmatch fetch-oligo-adapter-demo
	python3 scripts/bench_oligo_adapter.py --include-public --metadata examples/oligo_adapter/data/metadata.json

oligo-adapter-smoke-gate:
	python3 scripts/check_oligo_adapter_gate.py

oligo-adapter-public-gate:
	python3 scripts/check_oligo_adapter_gate.py --public

bench-amplicon-panel: dotmatch
	python3 scripts/bench_amplicon_panel.py

fetch-amplicon-panel-demo:
	python3 scripts/fetch_nfcore_amplicon_panel_demo.py --records "$${DOTMATCH_AMPLICON_PANEL_SUBSAMPLE:-20000}"

bench-amplicon-panel-public: dotmatch fetch-amplicon-panel-demo
	python3 scripts/bench_amplicon_panel.py --include-public --metadata examples/amplicon_panel/data/metadata.json

amplicon-panel-smoke-gate:
	python3 scripts/check_amplicon_panel_gate.py

amplicon-panel-public-gate:
	python3 scripts/check_amplicon_panel_gate.py --public

bench-feature-barcode: dotmatch
	python3 scripts/bench_feature_barcode.py

fetch-feature-barcode-demo:
	python3 scripts/fetch_10x_feature_barcode_demo.py --records "$${DOTMATCH_FEATURE_BARCODE_SUBSAMPLE:-20000}"

bench-feature-barcode-public: dotmatch fetch-feature-barcode-demo
	python3 scripts/bench_feature_barcode.py --include-public --metadata examples/feature_barcode/data/metadata.json

feature-barcode-smoke-gate:
	python3 scripts/check_feature_barcode_gate.py

feature-barcode-public-gate:
	python3 scripts/check_feature_barcode_gate.py

bench-perturb-seq: dotmatch
	python3 scripts/bench_perturb_seq.py

fetch-perturb-seq-demo:
	python3 scripts/fetch_10x_crispr_guide_demo.py --records "$${DOTMATCH_PERTURB_SEQ_SUBSAMPLE:-20000}"

bench-perturb-seq-public: dotmatch fetch-perturb-seq-demo
	python3 scripts/bench_perturb_seq.py --include-public --metadata examples/perturb_seq/data/metadata.json

perturb-seq-smoke-gate:
	python3 scripts/check_perturb_seq_gate.py

perturb-seq-public-gate:
	python3 scripts/check_perturb_seq_gate.py --public

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

validate-crispr-comparison-edlib: dotmatch edlib-tools
	python3 scripts/validate_crispr_comparison_edlib.py

bench-crispr-comparison: dotmatch competitor-env
	PATH="$(CURDIR)/build/guide-counter/bin:$(CURDIR)/build/competitor-env/bin:$$PATH" python3 scripts/run_crispr_comparison_repeated.py --run-mageck --run-guide-counter $${DOTMATCH_COMPARISON_FULL:+--full}

count-agreement-comparison:
	python3 scripts/compare_crispr_comparison_counts.py

crispr-comparison-report:
	python3 scripts/generate_crispr_comparison_report.py

crispr-comparison-gate:
	python3 scripts/check_crispr_comparison_gate.py

public-crispr-report:
	python3 scripts/generate_public_crispr_report.py

public-crispr-evidence-gate:
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

repository-ready:
	python3 scripts/check_repository_ready.py

release-ready: assay-evidence-ready alphabet-policy-ready citation-metadata-ready native-comparator-scope-ready workflow-examples-ready distribution-record-ready distribution-submission-ready bioconda-recipe-ready public-crispr-evidence-gate crispr-comparison-gate barcode-comparison-gate feature-barcode-public-gate perturb-seq-public-gate amplicon-panel-public-gate bcl-tiny-public-gate oligo-adapter-public-gate
	python3 scripts/check_release_readiness.py

pretag-ready:
	$(MAKE) test
	$(MAKE) cli-test
	$(MAKE) python-test
	$(MAKE) python-package-test
	$(MAKE) repository-ready
	$(MAKE) release-ready
	$(MAKE) coverage
	npm run lint
	npm audit --audit-level=moderate
	npm run build
	NEXT_OUTPUT=export NEXT_PUBLIC_BASE_PATH=/dotmatch NEXT_PUBLIC_SITE_URL=https://dnncha.github.io/dotmatch npm run build

assay-evidence-ready:
	python3 scripts/check_assay_evidence.py

alphabet-policy-ready:
	python3 scripts/check_alphabet_policy.py

citation-metadata-ready:
	python3 scripts/check_citation_metadata.py

native-comparator-scope-ready:
	python3 scripts/check_native_comparator_scope.py

workflow-examples-ready:
	python3 scripts/check_workflow_examples.py

workflow-adoption-status:
	python3 scripts/check_workflow_adoption.py

distribution-record-ready:
	python3 scripts/check_distribution_record.py

distribution-submission-ready:
	python3 scripts/check_distribution_submission.py

bioconda-recipe-ready:
	python3 scripts/check_bioconda_recipe.py

distribution-channels:
	python3 scripts/check_distribution_channels.py

asan:
	$(MAKE) clean
	$(MAKE) CFLAGS='-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -std=c11 -Wall -Wextra -Wpedantic -Iinclude' LDFLAGS='-fsanitize=address,undefined' test
	$(MAKE) clean

clean:
	rm -rf build dotmatch qda libdotmatch.a libdotmatch.so libdotmatch.dylib libqdalign.a libqdalign.so libqdalign.dylib
