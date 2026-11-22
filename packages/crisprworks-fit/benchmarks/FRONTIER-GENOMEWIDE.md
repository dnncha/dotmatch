# Genome-wide frontier: DepMap 24Q4

24 cell lines screened in both Avana and KY, selected as the first sorted eligible ModelIDs without looking at effects (176 were eligible).
Every UsedByChronos single-alignment gene with the modal guide count is fitted: Avana 14843 genes / 59372 guides / 61 samples, KY 13445 genes / 67225 guides / 66 samples.
Evaluation genes are Achilles common-essential and nonessential controls, held out from negative-control training in five folds.

## Average precision across cell lines

| Library | Method | Mean AP | Median AP | Min AP | Essential call rate | Nonessential call rate |
| --- | --- | --- | --- | --- | --- | --- |
| Avana | fit | 0.9629 | 0.9697 | 0.9052 | 0.714 | 0.028 |
| Avana | fit_gene_bh | 0.9629 | 0.9697 | 0.9052 | 0.226 | 0.002 |
| Avana | chronos | 0.9789 | 0.9809 | 0.9424 | 0.617 | 0.004 |
| Avana | jacks | 0.9787 | 0.9842 | 0.9384 | — | — |
| Avana | joint | 0.9787 | 0.9842 | 0.9384 | — | — |
| Avana | joint_gene_bh | 0.9787 | 0.9842 | 0.9384 | 0.616 | 0.003 |
| KY | fit | 0.9679 | 0.9705 | 0.9322 | 0.732 | 0.028 |
| KY | fit_gene_bh | 0.9679 | 0.9705 | 0.9322 | 0.124 | 0.002 |
| KY | chronos | 0.9767 | 0.9769 | 0.9617 | 0.422 | 0.003 |
| KY | jacks | 0.9800 | 0.9795 | 0.9641 | — | — |
| KY | joint | 0.9800 | 0.9795 | 0.9641 | — | — |
| KY | joint_gene_bh | 0.9800 | 0.9795 | 0.9641 | 0.647 | 0.002 |

## Paired differences across independent cell lines

Each cell line is one paired observation; intervals are 10,000 cell-level bootstrap draws (seed 1729) of the mean AP difference.

| Library | Candidate | Comparator | Mean AP difference | 95% interval | Cells where candidate is higher |
| --- | --- | --- | --- | --- | --- |
| Avana | joint | chronos | -0.0002 | [-0.0026, +0.0019] | 17/24 |
| Avana | fit | chronos | -0.0159 | [-0.0236, -0.0093] | 2/24 |
| Avana | joint | fit | +0.0157 | [+0.0107, +0.0215] | 24/24 |
| Avana | jacks | chronos | -0.0002 | [-0.0026, +0.0019] | 17/24 |
| KY | joint | chronos | +0.0033 | [+0.0019, +0.0046] | 21/24 |
| KY | fit | chronos | -0.0088 | [-0.0124, -0.0058] | 1/24 |
| KY | joint | fit | +0.0121 | [+0.0093, +0.0152] | 24/24 |
| KY | jacks | chronos | +0.0033 | [+0.0019, +0.0046] | 21/24 |

## Methods and limits

Chronos (commit d85e55f, TensorFlow 2.20) trains 1001 epochs with default regularization; its calls use the official empirical p-value and two-stage BH.
Chronos did not use the GPU meaningfully at this scale (0% utilisation, one CPU thread, about 29 minutes per fold); two folds ran with an L4 visible and eight on CPU only.
JACKS (commit dd5c4be) and Fit joint learn guide efficacy jointly across all panel cell lines in each library; Fit MLE fits each cell separately.
Fit gene BH calibrates frozen Fit effects with whole control genes, excluding training controls from the testing family.
Call rates use negative effects and directional q<=0.05; nonessential call rates are class-specific false-positive frequencies, not genome-wide FDR.
No copy-number correction is applied to any method. The Achilles reference sets derive partly from DepMap screens, so absolute AP is optimistic for every method; paired differences are the quantity of interest.
No hyperparameters were selected on evaluation annotations.
Fit joint hit its 50-iteration limit for 8 gene-fold fits (OR5D14 (219436), PRRC1 (133619), TMPRSS13 (84000)); those genes are left out of the joint gene-BH testing family and are never called, but keep their effects for ranking.
