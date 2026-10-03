# CRISPRWorks Fit: held-out gene essentiality

This evaluates ranking against CEGv2/NEGv1 annotations, not biological ground truth or calibrated FDR.
Both screens are HAP1 from the Hart Lab BAGEL repository; they are not independent cell lines or labs.
Five external stratified gene folds exclude each evaluation gene from BAGEL2's training labels.
All methods receive complete four-guide labels, selected before normalization; incomplete labels and large control bins are excluded.
Fit uses a T0 baseline and three T18 replicates. TKOv3 T3 and starvation samples are excluded.

| Screen | Method | AP (95% conditional bootstrap CI) | ROC AUC | Recall at 95% precision | Essential / nonessential call rate at nominal 5% FDR |
| --- | --- | --- | --- | --- | --- |
| hap1-tko | fit_beta | 0.9984 (0.9973–0.9994) | 0.9986 | 0.9935 | — |
| hap1-tko | fit_wald_z | 0.9987 (0.9975–0.9997) | 0.9988 | 0.9968 | — |
| hap1-tko | fit_finite_depletion_p | 0.9983 (0.9969–0.9994) | 0.9985 | 0.9935 | 0.0162 / 0.0000 |
| hap1-tko | fit_control_beta | 0.9983 (0.9969–0.9993) | 0.9984 | 0.9935 | — |
| hap1-tko | fit_control_depletion_p | 0.9983 (0.9970–0.9993) | 0.9984 | 0.9935 | 0.9451 / 0.0067 |
| hap1-tko | bagel2 | 0.9980 (0.9964–0.9993) | 0.9980 | 0.9919 | — |
| hap1-tko | mean_log2fc | 0.9980 (0.9965–0.9993) | 0.9980 | 0.9919 | — |
| hap1-tkov3 | fit_beta | 0.9974 (0.9947–0.9991) | 0.9967 | 0.9918 | — |
| hap1-tkov3 | fit_wald_z | 0.9947 (0.9915–0.9971) | 0.9946 | 0.9803 | — |
| hap1-tkov3 | fit_finite_depletion_p | 0.9974 (0.9947–0.9991) | 0.9967 | 0.9918 | 0.5721 / 0.0000 |
| hap1-tkov3 | fit_control_beta | 0.9984 (0.9970–0.9994) | 0.9985 | 0.9934 | — |
| hap1-tkov3 | fit_control_depletion_p | 0.9984 (0.9970–0.9994) | 0.9985 | 0.9934 | 0.9869 / 0.0123 |
| hap1-tkov3 | bagel2 | 0.9985 (0.9971–0.9994) | 0.9985 | 0.9918 | — |
| hap1-tkov3 | mean_log2fc | 0.9985 (0.9971–0.9995) | 0.9985 | 0.9918 | — |

AP is threshold-block average precision; all ties enter together. Larger scores mean greater depletion.
Confidence intervals use 500 class-stratified gene bootstrap samples with seed 1729, conditional on the fitted scores and training folds.
They do not capture variability across screens or refitting the training sets. Recall at 95% precision is a retrospective ranking diagnostic, not a deployable FDR threshold.
Fit scores are negative beta, negative Wald z, or negative finite depletion p. BAGEL2 is a supervised Bayes-factor classifier; mean log2FC is a simple baseline.
BAGEL2 uses its own sum-read normalization and pseudocount 5, no network boost or multi-target correction; Fit uses median normalization and updated guide efficiency.
The Fit control lane instead normalizes on nonessential training genes and draws its permutation background exclusively from their guides. Each scored gene is excluded from those controls.
Call rates use negative fitted effects with directional permutation FDR <=0.05. They measure sensitivity and false-positive frequency in the annotated classes; the nonessential call rate is not the false discovery rate among all genome-wide hits.
Finite tails repair tie/zero-p behavior. They do not change the fitted beta or Wald z and do not establish exchangeability of pooled pseudo-genes.

Chronos 2.3 hit calling and JACKS, copy-number effects, other cell lines, multi-condition contrasts, and independent validation remain untested here.
These measurements cannot support an overall state-of-the-art claim.
