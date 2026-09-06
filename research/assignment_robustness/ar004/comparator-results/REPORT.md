# Actual ReCo native-workflow comparison

The unchanged ReCo source at `e2daf48b610f8db29ad014bff5be8bb983aaa76f`, Cutadapt 2.8 and Bowtie 2.3.0 processed all seven complete original archives. This comparison changes extraction/context rules and assignment; it does not isolate one algorithmic choice or establish true biological accuracy.

The same baseline-only guide population and MAGeCK 0.5.9.5 settings used in the matched-position study were applied to its counts. Technical files were combined within the three biological replicates.

| Comparison | Common genes | Effect changes >=0.5 | Only left at FDR0.05 | Only ReCo at FDR0.05 | Rank correlation |
|---|---:|---:|---:|---:|---:|
| native_workflow: joint_best vs actual_ReCo | 19112 | 371 | 42 | 72 | 0.987690 |
| native_workflow: event_best vs actual_ReCo | 19112 | 384 | 46 | 71 | 0.987158 |
| native_workflow: joint_exact vs actual_ReCo | 19112 | 19 | 9 | 15 | 0.997539 |
| identical_input_repeat: actual_ReCo vs actual_ReCo | 19112 | 0 | 0 | 0 | 1.000000 |

All gene rows, FDR0.01/0.05/0.10 comparisons, individual archive count budgets, original commands and repeated-input checks are retained. Discordant calls are not labelled true/false discoveries. No current-version replacement or handwritten ReCo imitation was used. Original release acquisition/dependency failures are retained in their earlier workflow logs and do not enter these results.
