# Gene-level frontier comparison

A frozen panel compares Fit with Chronos 2.3.15 and JACKS 0.2 on three cell lines and two libraries.
Avana (Broad) and KY (Sanger) are curated Chronos vignette subsets enriched for reference controls; these are not full genome-wide screens.
Cell IDs are the first three sorted shared IDs with >=2 replicates: T47D, L363 and NCI-H1915.
Every evaluation gene is excluded from negative-control training in its external fold.

| Library | Cell | Method | AP | Essential calls / class size | Nonessential calls / class size |
| --- | --- | --- | --- | --- | --- |
| Avana | ACH-000147 | fit | 0.9857 | 288/303 | 26/322 |
| Avana | ACH-000147 | fit_gene_bh | 0.9857 | 276/303 | 9/322 |
| Avana | ACH-000147 | fit_gene_by | 0.9857 | 0/303 | 0/322 |
| Avana | ACH-000147 | chronos | 0.9874 | 270/303 | 3/322 |
| Avana | ACH-000147 | jacks | 0.9860 | — | — |
| Avana | ACH-000147 | joint | 0.9860 | — | — |
| Avana | ACH-000147 | joint_gene_bh | 0.9860 | 275/303 | 7/322 |
| Avana | ACH-000147 | joint_gene_by | 0.9860 | 0/303 | 0/322 |
| Avana | ACH-000183 | fit | 0.9900 | 288/303 | 12/322 |
| Avana | ACH-000183 | fit_gene_bh | 0.9900 | 285/303 | 9/322 |
| Avana | ACH-000183 | fit_gene_by | 0.9900 | 0/303 | 0/322 |
| Avana | ACH-000183 | chronos | 0.9882 | 264/303 | 1/322 |
| Avana | ACH-000183 | jacks | 0.9921 | — | — |
| Avana | ACH-000183 | joint | 0.9921 | — | — |
| Avana | ACH-000183 | joint_gene_bh | 0.9921 | 290/303 | 7/322 |
| Avana | ACH-000183 | joint_gene_by | 0.9921 | 0/303 | 0/322 |
| Avana | ACH-000434 | fit | 0.9557 | 251/303 | 18/322 |
| Avana | ACH-000434 | fit_gene_bh | 0.9557 | 213/303 | 6/322 |
| Avana | ACH-000434 | fit_gene_by | 0.9557 | 0/303 | 0/322 |
| Avana | ACH-000434 | chronos | 0.9675 | 254/303 | 8/322 |
| Avana | ACH-000434 | jacks | 0.9675 | — | — |
| Avana | ACH-000434 | joint | 0.9675 | — | — |
| Avana | ACH-000434 | joint_gene_bh | 0.9675 | 250/303 | 7/322 |
| Avana | ACH-000434 | joint_gene_by | 0.9675 | 0/303 | 0/322 |
| KY | ACH-000147 | fit | 0.9616 | 227/277 | 13/261 |
| KY | ACH-000147 | fit_gene_bh | 0.9616 | 206/277 | 5/261 |
| KY | ACH-000147 | fit_gene_by | 0.9616 | 0/277 | 0/261 |
| KY | ACH-000147 | chronos | 0.9668 | 213/277 | 1/261 |
| KY | ACH-000147 | jacks | 0.9750 | — | — |
| KY | ACH-000147 | joint | 0.9750 | — | — |
| KY | ACH-000147 | joint_gene_bh | 0.9750 | 235/277 | 5/261 |
| KY | ACH-000147 | joint_gene_by | 0.9750 | 0/277 | 0/261 |
| KY | ACH-000183 | fit | 0.9815 | 239/277 | 5/261 |
| KY | ACH-000183 | fit_gene_bh | 0.9815 | 243/277 | 7/261 |
| KY | ACH-000183 | fit_gene_by | 0.9815 | 0/277 | 0/261 |
| KY | ACH-000183 | chronos | 0.9769 | 228/277 | 2/261 |
| KY | ACH-000183 | jacks | 0.9861 | — | — |
| KY | ACH-000183 | joint | 0.9861 | — | — |
| KY | ACH-000183 | joint_gene_bh | 0.9861 | 244/277 | 6/261 |
| KY | ACH-000183 | joint_gene_by | 0.9861 | 0/277 | 0/261 |
| KY | ACH-000434 | fit | 0.9056 | 176/277 | 10/261 |
| KY | ACH-000434 | fit_gene_bh | 0.9056 | 76/277 | 3/261 |
| KY | ACH-000434 | fit_gene_by | 0.9056 | 0/277 | 0/261 |
| KY | ACH-000434 | chronos | 0.9231 | 134/277 | 3/261 |
| KY | ACH-000434 | jacks | 0.9233 | — | — |
| KY | ACH-000434 | joint | 0.9233 | — | — |
| KY | ACH-000434 | joint_gene_bh | 0.9233 | 156/277 | 4/261 |
| KY | ACH-000434 | joint_gene_by | 0.9233 | 0/277 | 0/261 |

## Paired ranking differences

Positive AP differences favor production Fit joint. Percentile intervals use 2,000 paired, class-stratified gene bootstrap draws (seed 1729), conditional on frozen fits and folds.
They do not measure uncertainty across independent screens, labels or refitted training sets; the six comparisons are correlated.

| Library | Cell | Comparator | AP difference | Conditional 95% interval |
| --- | --- | --- | --- | --- |
| Avana | ACH-000147 | fit | +0.0003 | [-0.0041, +0.0052] |
| Avana | ACH-000147 | chronos | -0.0014 | [-0.0045, +0.0019] |
| Avana | ACH-000183 | fit | +0.0021 | [-0.0013, +0.0058] |
| Avana | ACH-000183 | chronos | +0.0038 | [+0.0005, +0.0072] |
| Avana | ACH-000434 | fit | +0.0118 | [+0.0059, +0.0185] |
| Avana | ACH-000434 | chronos | -0.0000 | [-0.0044, +0.0044] |
| KY | ACH-000147 | fit | +0.0134 | [+0.0051, +0.0229] |
| KY | ACH-000147 | chronos | +0.0083 | [+0.0007, +0.0165] |
| KY | ACH-000183 | fit | +0.0046 | [-0.0007, +0.0102] |
| KY | ACH-000183 | chronos | +0.0092 | [+0.0033, +0.0157] |
| KY | ACH-000434 | fit | +0.0177 | [+0.0079, +0.0303] |
| KY | ACH-000434 | chronos | +0.0002 | [-0.0059, +0.0059] |

## Methods and limits

Fit fits each cell separately; Chronos and JACKS learn guide efficacy jointly across the three cells within each library.
All methods see the same finite-read complete-guide cohort (four guides Avana; five guides KY) and the same training controls.
Chronos trains for a predeclared 1001 epochs with default model regularization and its official empirical p-value / two-stage BH caller.
JACKS retains unmodified upstream inference and control-guide normalization with pseudocount 32; only effect ranking is evaluated here.
Fit joint adapts that tested JACKS model. Every fitted effect is checked against independent unmodified JACKS at 1e-10 absolute/relative tolerance in the same NumPy 1.26.4 environment.
Reported Fit joint runs use the production NumPy 2.3.5 dependency set. NumPy version changes alter variance-window ordering at near-tied means in some guides; cross-version differences are recorded rather than described as exact parity.
Fit uses updated efficiency, control normalization, two pooled-guide permutation rounds, and median/variance conventions inherited from MAGeCK2.
The gene-level Fit callers calibrate those frozen effects using whole control genes within guide-count strata, excluding training controls from the testing family.
Fit MLE families cover one cell per fit; joint families cover all three conditions together. These call counts therefore compare complete analysis choices, not only effect estimators.
BH assumes independent or suitable positively dependent valid p-values; BY tolerates arbitrary dependence of valid p-values, at a substantial power cost.
Calls use negative effects and directional q<=0.05. Nonessential call rates are class-specific false-positive frequencies, not measured genome-wide FDR.
Fit gene callers have identical AP to Fit because calibration does not refit effects. No hyperparameters were selected on the evaluation annotations.

This panel does not establish overall SOTA: it is enriched for controls, has only three cell lines, has no copy-number correction, and does not independently validate biological truth.
The annotation collections and source screen populations may overlap; gene holdout cannot remove that collection-level bias.
