# Correlated-guide null calibration

This known-truth simulation isolates calibration of frozen scores. It does not simulate read counts, fit MAGeCK/JACKS models or measure biological accuracy.
Each of 200 seeds supplies 600 test genes, 4,000 independently sampled control genes and four guides per gene.
Guide residual SD is 0.35. The correlated scenario adds one shared N(0,1) artifact to all guides of each gene.
Mixed scenarios shift 20% of test genes by -2.5. All callers use the same controls and test family at directional q<=0.05.
The pooled null resamples individual control guides into 12,000 four-guide averages; the gene null retains the 4,000 control gene averages.
The reported mean false discovery proportion estimates FDR in this specified simulation. Intervals are normal Monte Carlo intervals for the replicate mean, not a guarantee across unknown biological scenarios.

| Shared artifact | True effects | Method | Mean FDP [MC 95% interval] | Mean power | P(any false call) |
| --- | --- | --- | --- | --- | --- |
| False | False | pooled_guide_bh | 0.070 [0.035, 0.105] | — | 0.070 |
| False | False | whole_gene_bh | 0.000 [0.000, 0.000] | — | 0.000 |
| False | False | whole_gene_by | 0.000 [0.000, 0.000] | — | 0.000 |
| False | True | pooled_guide_bh | 0.040 [0.038, 0.043] | 1.000 | 0.985 |
| False | True | whole_gene_bh | 0.040 [0.038, 0.043] | 1.000 | 0.985 |
| False | True | whole_gene_by | 0.005 [0.004, 0.006] | 1.000 | 0.400 |
| True | False | pooled_guide_bh | 1.000 [1.000, 1.000] | — | 1.000 |
| True | False | whole_gene_bh | 0.000 [0.000, 0.000] | — | 0.000 |
| True | False | whole_gene_by | 0.000 [0.000, 0.000] | — | 0.000 |
| True | True | pooled_guide_bh | 0.356 [0.351, 0.361] | 0.909 | 1.000 |
| True | True | whole_gene_bh | 0.038 [0.034, 0.041] | 0.438 | 0.855 |
| True | True | whole_gene_by | 0.003 [0.001, 0.004] | 0.068 | 0.075 |

Under shared artifacts with true effects, pooled guide averaging understates null variance: mean FDP rises above the nominal 0.05 level.
The whole-gene null preserves that correlation and reduces false discoveries, with a substantial power cost. BY is more conservative still.
Under null-only scenarios, finite p-value resolution can prevent all calls; zero observed false calls does not prove zero population risk.
A shared finite control set induces dependence between empirical test p-values. BH validity requires appropriate independence or positive dependence; this simulation alone does not establish that for all assays.
BY accommodates arbitrary dependence of valid marginal p-values. Both methods still require exchangeable, representative control genes within guide-count strata.
