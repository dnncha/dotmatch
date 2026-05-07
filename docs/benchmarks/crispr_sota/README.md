# CRISPR SOTA Evidence

This report is generated from raw CSV artifacts. It is intentionally stricter than the public smoke report: SOTA claims require both MAGeCK/Yusa and Sanson/Brunello real-data rows, competitor rows, count agreement, and Edlib validation.

## Claim Discipline

- Hamming `k=1` rows are the fair guide-counter lane: one mismatch, no indels.
- Levenshtein `k=1` rows are the DotMatch differentiator lane: substitutions plus single-base insertions/deletions, with Edlib validation.
- Full FASTQ rows are reported separately from repeated subsamples.
- Broad SOTA wording is blocked unless `make crispr-sota-gate` passes.

## Throughput Figure

![CRISPR SOTA throughput](../../../benchmarks/figures/crispr_sota_throughput.svg)

## Repeated Subsample Rows

|dataset|tool|records_per_sample|repeats|mean_reads_per_sec|p50_reads_per_sec|p95_reads_per_sec|cv|max_peak_rss_mb|mean_verified_per_read|
|---|---|---|---|---|---|---|---|---|---|
|mageck_yusa|dotmatch_exact_k0|100000|5|1226604.5|1236776.2|1307857.7|0.0676|29.9|0.894|
|mageck_yusa|dotmatch_exact_k0|5|1|484667.9|484667.9|484667.9|0.0000|27.7|0.894|
|mageck_yusa|dotmatch_hamming_k1|100000|5|648130.4|646694.0|667286.6|0.0214|49.4|0.996|
|mageck_yusa|dotmatch_hamming_k1|5|1|280311.2|280311.2|280311.2|0.0000|44.5|0.996|
|mageck_yusa|dotmatch_levenshtein_k1|100000|5|107288.7|107211.0|108555.2|0.0079|121.8|0.993|
|mageck_yusa|dotmatch_levenshtein_k1|5|1|39516.0|39516.0|39516.0|0.0000|113.0|1.039|
|mageck_yusa|guide_counter_one_mismatch|100000|5|183944.9|183334.5|197431.9|0.0475|528.7||
|mageck_yusa|mageck_count_exact|100000|5|133813.8|139958.8|143631.3|0.0782|146.1||
|sanson_brunello|dotmatch_exact_k0|100000|5|877105.2|874365.8|897482.5|0.0160|68.9|0.805|
|sanson_brunello|dotmatch_exact_k0|5|1|134.7|134.7|134.7|0.0000|28.8|0.800|
|sanson_brunello|dotmatch_hamming_k1|100000|5|651211.1|647760.6|661037.9|0.0137|196.9|0.873|
|sanson_brunello|dotmatch_hamming_k1|5|1|57.6|57.6|57.6|0.0000|156.8|0.800|
|sanson_brunello|dotmatch_levenshtein_k1|100000|5|14135.6|14017.0|14633.9|0.0198|130.9|8.872|
|sanson_brunello|dotmatch_levenshtein_k1|5|1|83.9|83.9|83.9|0.0000|110.3|0.800|
|sanson_brunello|guide_counter_one_mismatch|100000|5|289641.0|297252.2|311829.6|0.0808|527.7||
|sanson_brunello|mageck_count_exact|100000|5|272922.1|261472.2|293751.4|0.0638|114.2||


## Full FASTQ Rows

|dataset|tool|records_per_sample|repeats|mean_reads_per_sec|mean_seconds|max_peak_rss_mb|mean_verified_per_read|
|---|---|---|---|---|---|---|---|
|mageck_yusa|dotmatch_exact_k0|full|1|5136342.8|3.9707|29.9|0.838|
|mageck_yusa|dotmatch_hamming_k1|full|1|4701710.0|4.3377|49.4|1.000|
|mageck_yusa|dotmatch_levenshtein_k1|full|1|77265.8|263.9547|117.8|21.718|
|mageck_yusa|guide_counter_one_mismatch|full|1|2974449.2|6.8566|516.2||
|mageck_yusa|mageck_count_exact|full|1|479694.0|42.5160|138.7||
|sanson_brunello|dotmatch_exact_k0|full|1|3838096.0|64.3419|78.4|0.835|
|sanson_brunello|dotmatch_hamming_k1|full|1|3568639.0|69.2002|206.4|0.903|
|sanson_brunello|dotmatch_levenshtein_k1|full|1|35750.2|6907.6665|146.3|0.955|
|sanson_brunello|guide_counter_one_mismatch|full|1|2320384.8|106.4265|527.7||
|sanson_brunello|mageck_count_exact|full|1|520740.1|18.8599|104.1||


## Full Hamming Speed Check

|dataset|dotmatch_hamming_reads_per_sec|guide_counter_reads_per_sec|speedup|status|
|---|---|---|---|---|
|mageck_yusa|4701710.0|2974449.2|1.58|pass|
|sanson_brunello|3568639.0|2320384.8|1.54|pass|


## Edlib Oracle Validation

|dataset|sample|checked_reads|mismatches|oracle_strategy|edlib_alignments|bounded_windows|fallback_windows|selected_target_start|stratum_exact|stratum_corrected|stratum_ambiguous|stratum_unmatched|stratum_contains_n|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|mageck_yusa|plasmid|10000|0|bounded_edlib_candidates|6586091|29925|75|23|8984|755|20|241|25|
|mageck_yusa|ESC1|10000|0|bounded_edlib_candidates|5536726|29937|63|23|8910|861|18|211|21|
|sanson_brunello|plasmid|10000|0|bounded_edlib_candidates|43488|290000|0|21|8826|754|10|410|0|
|sanson_brunello|RepA|10000|0|bounded_edlib_candidates|7005889|229910|90|22|7830|1053|4|1113|5|
|sanson_brunello|RepB|10000|0|bounded_edlib_candidates|10489798|259865|135|22|7237|970|7|1786|5|
|sanson_brunello|RepC|10000|0|bounded_edlib_candidates|8559115|289890|110|21|8173|988|6|833|4|


## Count Agreement

|dataset|comparison|status|n_guides|total_delta|differing_guides|max_abs_delta|pearson|spearman|
|---|---|---|---|---|---|---|---|---|
|mageck_yusa|mageck_yusa:dotmatch_hamming_vs_guide_counter|ok|87437|-24533|13537|26|0.94176266|0.95124146|
|mageck_yusa|mageck_yusa:dotmatch_exact_vs_mageck_exact|ok|87437|0|0|0|1.00000000|1.00000000|
|sanson_brunello|sanson_brunello:dotmatch_hamming_vs_guide_counter|ok|77441|-1709|764|21|0.99622960|0.99642069|
|sanson_brunello|sanson_brunello:dotmatch_exact_vs_mageck_exact|ok|77441|321536|67253|71|nan|nan|


## Raw Inputs

- `benchmarks/raw/crispr_sota_repeated.csv`
- `benchmarks/raw/crispr_sota_edlib_validation.csv`
- `benchmarks/raw/crispr_sota_count_agreement_summary.csv`
