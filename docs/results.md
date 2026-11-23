# Synthetic verification

60 artificial samples, one binary target, three outer folds and seed 42. Each entry gives the outer-fold mean and sample standard deviation.

| Features | Model | AUROC | Balanced accuracy |
|---|---|---:|---:|
| combined | tensorflow | 0.513 ± 0.119 | 0.550 ± 0.180 |
| combined | xgboost | 0.997 ± 0.006 | 0.950 ± 0.050 |
| protein | tensorflow | 0.790 ± 0.130 | 0.683 ± 0.058 |
| protein | xgboost | 0.997 ± 0.006 | 0.967 ± 0.029 |
| rna | tensorflow | 0.877 ± 0.150 | 0.750 ± 0.180 |
| rna | xgboost | 1.000 ± 0.000 | 1.000 ± 0.000 |

These scores are from the replacement archive’s verification run. They measure recovery of artificial signal. Real-cohort predictive performance and raw-search QC remain unverified.
