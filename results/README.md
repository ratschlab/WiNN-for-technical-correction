# Results

`final/` contains the compact publication results, including aggregate clinical statistics. Large corrected matrices, candidate fits, private clinical inputs and scheduler files remain outside the public repository. Current WiNN outputs use development version 0.1.7.9999 with permutation maxT drift calibration; Raw and comparator outputs retain their preceding published values.

## Reading the comparison

The primary WiNN configuration uses supplied acquisition segments and fixed settings without QC identities. QC-tuned and inferred-segment modes are separate supplementary configurations; a mode is not chosen independently for each reported endpoint. Held-out QCs, biological labels and replicate identities are excluded from configuration selection.

Technical cleanup and preservation measure different properties. Held-out QC CV, residual acquisition-order structure and batch association do not establish biological truth in empirical data. Read them together with coverage, cross-batch matching and each dataset's biological endpoints. Uncertainty is conditional on the fitted correction; the bootstrap does not refit the full correction or selection pipeline.

## Calibrated results and limitations

Across ten simulation seeds, fixed WiNN and QC-tuned WiNN retain lower mean truth RMSE and higher mean ICC and sample-profile Pearson correlation than the included comparator methods. Fixed WiNN's mean RMSE is 0.07444, ICC 0.91863 and Pearson correlation 0.91993. Inferred segmentation is less reliable in this design with known acquisition batches.

The forced-correction comparison changes after calibration. Selective/selective correction has worse truth RMSE, ICC and Pearson correlation than both forced-drift policies in all ten seeds, while beating selective-drift/forced-batch correction on those endpoints. Its held-out QC CV is lower than both policies that force batch correction in all ten seeds, and lower than forced-drift/selective-batch correction in seven. Thus valid drift-gate calibration does not imply that selective drift correction maximizes every preservation metric.

The partial-confounding experiment is a stylized sensitivity analysis. Under complete joint confounding, mean recovery relative to the clean effect is 0.740, sensitivity is 0.507 and the mean false-discovery proportion is 0.359. These results show a recovery/error tradeoff; they do not establish a general identifiability boundary or validate a confounded empirical design.

At fixed segments and a prespecified threshold, the permutation gate's error guarantee requires joint exchangeability of true-null profiles. Scheduled sample-type shifts, data-selected segmentation, QC-guided tuning and multiple segments jointly require separate interpretation. Calibration experiments and their Monte Carlo uncertainty are in `manuscript/generated/calibration_revision/`.

## Files

Compact primary and simulation tables are in `final/tables/`. Publication figure sources, detailed estimates and support tables are in `manuscript/generated/`, with vector figures in `manuscript/figures/`. Current supplementary and public derivative archives are `final/revision_inputs/supplement_20261114/` and `derivative_20261114/`; each has an uncompressed SHA-256 manifest. Historical archives remain version-labelled. The manuscript sources provide the full empirical interpretation.
