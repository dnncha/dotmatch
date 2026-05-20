# CRISPR Comparison Evidence

This report is generated from raw CSV artifacts. It is intentionally stricter than the public smoke report: comparison rows require both MAGeCK/Yusa and Sanson/Brunello real-data rows, competitor rows, count agreement, and Edlib validation.

## Evidence Boundary

- Hamming `k=1` rows are the fair guide-counter lane: one mismatch, no indels.
- Levenshtein `k=1` rows are the DotMatch differentiator lane: substitutions plus single-base insertions/deletions, with Edlib validation.
- Full FASTQ rows are reported separately from repeated subsamples.
- guide-counter speed ratios are reported when present; they are not universal replacement gates.
- Broad comparisons require `make crispr-comparison-gate` to pass.

## Throughput Figure

![CRISPR comparison throughput](../../../benchmarks/figures/crispr_comparison_throughput.svg)

## Repeated Subsample Rows

|dataset|tool|records_per_sample|repeats|mean_reads_per_sec|p50_reads_per_sec|p95_reads_per_sec|cv|max_peak_rss_mb|mean_verified_per_read|
|---|---|---|---|---|---|---|---|---|---|
|mageck_yusa|dotmatch_exact_k0|10000|5|160359.4|158280.7|165730.3|0.0304|111.8|0.895|
|mageck_yusa|dotmatch_exact_k0|100000|5|887206.3|960437.0|1011090.6|0.1482|113.4|0.894|
|mageck_yusa|dotmatch_hamming_k1|10000|5|94076.1|96138.0|100448.8|0.0601|115.8|0.929|
|mageck_yusa|dotmatch_hamming_k1|100000|5|210421.0|216553.9|226452.1|0.1099|121.4|0.928|
|mageck_yusa|dotmatch_levenshtein_k1|10000|5|6921.9|7016.6|7108.6|0.0419|111.9|2.828|
|mageck_yusa|dotmatch_levenshtein_k1|100000|5|6635.0|6592.2|7371.6|0.0806|113.5|2.822|
|mageck_yusa|guide_counter_one_mismatch|10000|5|26684.9|26003.1|30334.1|0.0787|528.7||
|mageck_yusa|guide_counter_one_mismatch|100000|5|191862.5|210721.3|235253.3|0.2325|528.7||
|mageck_yusa|mageck_count_exact|10000|5|21856.7|23041.7|25733.9|0.1830|133.6||
|mageck_yusa|mageck_count_exact|100000|5|127848.7|124634.5|147025.8|0.0901|152.9||
|sanson_brunello|dotmatch_exact_k0|10000|5|103809.6|100934.5|122058.7|0.1384|120.2|0.803|
|sanson_brunello|dotmatch_exact_k0|100000|5|167486.7|166245.7|193132.4|0.1198|122.2|0.805|
|sanson_brunello|dotmatch_hamming_k1|10000|5|75635.5|76519.1|80993.8|0.0587|245.2|0.873|
|sanson_brunello|dotmatch_hamming_k1|100000|5|141126.5|147940.8|161900.9|0.1213|251.2|0.875|
|sanson_brunello|dotmatch_levenshtein_k1|10000|5|878.1|863.9|937.8|0.0626|113.2|3.889|
|sanson_brunello|dotmatch_levenshtein_k1|100000|5|803.4|813.7|853.8|0.0561|112.7|3.963|
|sanson_brunello|guide_counter_one_mismatch|10000|5|46552.4|48594.2|52361.9|0.1547|527.7||
|sanson_brunello|guide_counter_one_mismatch|100000|5|332365.0|360760.8|386608.9|0.2047|527.7||
|sanson_brunello|mageck_count_exact|10000|5|37236.0|41500.2|41965.8|0.1666|114.8||
|sanson_brunello|mageck_count_exact|100000|5|219199.4|232963.0|247063.8|0.1956|115.1||


## Full FASTQ Rows

_No rows available._


## Full Hamming Guide-Counter Ratio

_No rows available._


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
|mageck_yusa|mageck_yusa:dotmatch_hamming_vs_guide_counter|ok|87437|-26043|13852|27|0.93492501|0.94108416|
|mageck_yusa|mageck_yusa:dotmatch_exact_vs_mageck_exact|ok|87437|0|0|0|1.00000000|1.00000000|
|sanson_brunello|sanson_brunello:dotmatch_hamming_vs_guide_counter|ok|77441|-2352|817|24|0.99330729|0.99268340|
|sanson_brunello|sanson_brunello:dotmatch_exact_vs_mageck_exact|ok|77441|321536|67253|71|nan|nan|


## Raw Inputs

- `benchmarks/raw/crispr_comparison_repeated.csv`
- `benchmarks/raw/crispr_comparison_edlib_validation.csv`
- `benchmarks/raw/crispr_comparison_count_agreement_summary.csv`
