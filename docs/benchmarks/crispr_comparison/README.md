# CRISPR Comparison Evidence

This report is generated from raw CSV artifacts. It is intentionally stricter than the public smoke report: comparison rows require both MAGeCK/Yusa and Sanson/Brunello real-data rows, competitor rows, count agreement, and Edlib validation.

## Evidence Boundary

- Hamming `k=1` rows are the fair guide-counter lane: one mismatch, no indels.
- Levenshtein `k=1` rows are the DotMatch differentiator lane: substitutions plus single-base insertions/deletions, with Edlib validation.
- Full FASTQ rows are reported separately from repeated subsamples.
- guide-counter speed ratios are reported when present; they are not universal replacement gates.
- Broad comparisons require `make crispr-comparison-gate` to pass.

## GuideCounter Compatibility Note

The compatibility mode documented for this release is:

```bash
dotmatch guide-counter count \
  --input sample.fastq.gz \
  --library guides.tsv \
  --output guide_counts
```

It also accepts the `guide-counter-count` and `guide-count` aliases. By default
the wrapper uses Hamming `k=1`, no indels, best-distance assignment, automatic
multi-offset guide-window detection, `--offset-sample-size 100000`, and
`--offset-min-fraction 0.0025`; `--exact-match` switches the delegated count
engine to exact matching. For an output prefix such as `guide_counts`, it writes
`guide_counts.counts.txt`, `guide_counts.extended-counts.txt`, and
`guide_counts.stats.txt`.

This is a command-line and file-format compatibility bridge. The assignment and
counting step remains DotMatch's deterministic CPU path. The GPU/backend
optimizer row below is advisory evidence only and does not authorize unchecked
production GPU assignment.

## Throughput Figure

![CRISPR comparison throughput](../../../benchmarks/figures/crispr_comparison_throughput.svg)

## Repeated Subsample Rows

|dataset|tool|records_per_sample|repeats|mean_reads_per_sec|p50_reads_per_sec|p95_reads_per_sec|cv|max_peak_rss_mb|mean_verified_per_read|
|---|---|---|---|---|---|---|---|---|---|
|mageck_yusa|dotmatch_exact_k0|10000|5|160359.4|158280.7|165730.3|0.0304|111.8|0.895|
|mageck_yusa|dotmatch_exact_k0|100000|5|887206.3|960437.0|1011090.6|0.1482|113.4|0.894|
|mageck_yusa|dotmatch_hamming_k1|10000|5|213575.6|217085.7|228121.1|0.0706|34.6|0.997|
|mageck_yusa|dotmatch_hamming_k1|100000|5|754902.5|756958.1|785171.6|0.0464|49.3|0.996|
|mageck_yusa|dotmatch_levenshtein_k1|10000|5|6921.9|7016.6|7108.6|0.0419|111.9|2.828|
|mageck_yusa|dotmatch_levenshtein_k1|100000|5|6635.0|6592.2|7371.6|0.0806|113.5|2.822|
|mageck_yusa|guide_counter_one_mismatch|10000|5|22849.6|23567.2|27945.6|0.1782|528.7||
|mageck_yusa|guide_counter_one_mismatch|100000|5|184061.4|183607.3|202375.6|0.0625|528.7||
|mageck_yusa|mageck_count_exact|10000|5|21856.7|23041.7|25733.9|0.1830|133.6||
|mageck_yusa|mageck_count_exact|100000|5|127848.7|124634.5|147025.8|0.0901|152.9||
|sanson_brunello|dotmatch_exact_k0|10000|5|103809.6|100934.5|122058.7|0.1384|120.2|0.803|
|sanson_brunello|dotmatch_exact_k0|100000|5|167486.7|166245.7|193132.4|0.1198|122.2|0.805|
|sanson_brunello|dotmatch_hamming_k1|10000|5|136650.6|152481.8|157797.5|0.2752|170.2|0.875|
|sanson_brunello|dotmatch_hamming_k1|100000|5|721462.1|719108.1|741919.6|0.0169|196.8|0.874|
|sanson_brunello|dotmatch_levenshtein_k1|10000|5|878.1|863.9|937.8|0.0626|113.2|3.889|
|sanson_brunello|dotmatch_levenshtein_k1|100000|5|803.4|813.7|853.8|0.0561|112.7|3.963|
|sanson_brunello|guide_counter_one_mismatch|10000|5|43986.2|45605.9|49264.6|0.1114|527.7||
|sanson_brunello|guide_counter_one_mismatch|100000|5|302166.4|298310.9|327193.6|0.0759|527.7||
|sanson_brunello|mageck_count_exact|10000|5|37236.0|41500.2|41965.8|0.1666|114.8||
|sanson_brunello|mageck_count_exact|100000|5|219199.4|232963.0|247063.8|0.1956|115.1||


## Full FASTQ Rows

|dataset|tool|records_per_sample|repeats|mean_reads_per_sec|mean_seconds|max_peak_rss_mb|mean_verified_per_read|
|---|---|---|---|---|---|---|---|
|sanson_brunello|dotmatch_exact_k0|full|1|261396.8|944.7340|170.7|0.837|
|sanson_brunello|dotmatch_hamming_k1|full|1|634950.2|388.9288|198.0|0.905|
|sanson_brunello|guide_counter_one_mismatch|full|1|473609.5|521.4220|528.7||


## Guide-Counter-Style Public Paper-Data Lane

DotMatch `dotmatch_hamming_k1` versus `guide_counter_one_mismatch` on the public paper-data inputs. This lane uses best-distance Hamming assignment with guide-counter's offset threshold, is limited to one mismatch and no indels, and keeps Levenshtein rows as a separate DotMatch capability lane.

|dataset|records_per_sample|dotmatch_hamming_reads_per_sec|guide_counter_reads_per_sec|speedup|count_agreement_status|count_total_delta|semantics|
|---|---|---|---|---|---|---|---|
|mageck_yusa|10000|213575.6|22849.6|9.35|ok|-24533|one mismatch, no indels|
|mageck_yusa|100000|754902.5|184061.4|4.10|ok|-24533|one mismatch, no indels|
|sanson_brunello|10000|136650.6|43986.2|3.11|ok|-1190|one mismatch, no indels|
|sanson_brunello|100000|721462.1|302166.4|2.39|ok|-1190|one mismatch, no indels|
|sanson_brunello|full|634950.2|473609.5|1.34|ok|-1190|one mismatch, no indels|


## Full Hamming Guide-Counter Ratio

|dataset|dotmatch_hamming_reads_per_sec|guide_counter_reads_per_sec|speedup|status|
|---|---|---|---|---|
|sanson_brunello|634950.2|473609.5|1.34|reported|


## Backend Optimizer

The optimizer is advisory and CPU-authoritative: it records the fastest eligible candidate backend, but speed claims still require CPU checksum agreement.

|dataset|optimizer|authority|selected_backend|candidate_backend|recommendation|expected_speedup_band|estimated_total_speedup|
|---|---|---|---|---|---|---|---|
|sanson_brunello|local_benchmark_informed_scorer_v1|cpu|cpu|gpu-metal-experimental|gpu_candidate_requires_cpu_validation|1.5-3x|1.92|


## Edlib Oracle Validation

|dataset|sample|checked_reads|mismatches|oracle_strategy|edlib_alignments|bounded_windows|fallback_windows|selected_target_start|stratum_exact|stratum_corrected|stratum_ambiguous|stratum_unmatched|stratum_contains_n|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|mageck_yusa|plasmid|10000|0|bounded_edlib_candidates|6586091|29925|75|23|8486|755|518|241|25|
|mageck_yusa|ESC1|10000|0|bounded_edlib_candidates|5536726|29937|63|23|8450|861|478|211|21|
|sanson_brunello|plasmid|10000|0|bounded_edlib_candidates|43488|290000|0|21|8684|754|152|410|0|
|sanson_brunello|RepA|10000|0|bounded_edlib_candidates|7005889|229910|90|22|7722|1053|112|1113|5|
|sanson_brunello|RepB|10000|0|bounded_edlib_candidates|10489798|259865|135|22|7129|970|115|1786|5|
|sanson_brunello|RepC|10000|0|bounded_edlib_candidates|8559115|289890|110|21|8048|988|131|833|4|


## Count Agreement

|dataset|comparison|status|n_guides|total_delta|differing_guides|max_abs_delta|pearson|spearman|
|---|---|---|---|---|---|---|---|---|
|mageck_yusa|mageck_yusa:dotmatch_hamming_vs_guide_counter|ok|87437|-24533|13537|26|0.94176266|0.95124146|
|mageck_yusa|mageck_yusa:dotmatch_exact_vs_mageck_exact|ok|87437|0|0|0|1.00000000|1.00000000|
|sanson_brunello|sanson_brunello:dotmatch_hamming_vs_guide_counter|ok|77441|-1190|255|21|0.99639346|0.99659479|
|sanson_brunello|sanson_brunello:dotmatch_exact_vs_mageck_exact|non_comparable|77441|321536|67253|71|nan|nan|


## Raw Inputs

- `benchmarks/raw/crispr_comparison_repeated.csv`
- `benchmarks/raw/crispr_comparison_full_sanson_atlas_latest_dotmatch.csv`
- `benchmarks/raw/crispr_comparison_edlib_validation.csv`
- `benchmarks/raw/crispr_comparison_count_agreement_summary.csv`
- `benchmarks/raw/crispr_sanson_brunello_backend_optimization_atlas_latest_dotmatch.json`
