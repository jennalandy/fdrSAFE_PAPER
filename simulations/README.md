# Simulation Results

These tests compare `fdrSAFE`, partial implementations, and benchmarks on Symmetric, Asymmetric, and Curated Ovarian Data-Based simulation studies. Use the following steps to create Figure 1, Supplementary Figure A.1, and Supplementary Table B.1.

1. Run slurm jobs for model comparisons `test_symmetric`, `test_asymmetric`, `test_COD_based` and for hyperparameter choices `test_symmetric_hyperparameters` and `test_asymmetric_hyperparameters`.
2. Run through `1_simulations.qmd` to create Figure 1, Supplementary Figure A.1, and Supplementary Table B.1