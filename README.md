# Mapping Heterogeneity in the Neuroanatomical Correlates of Depression

Analysis code for a study examining how the depth of depression phenotyping influences neuroanatomic profiles in UK Biobank neuroimaging data (n = 30,122). Methods include multivariable logistic regression, random forest classifiers with Bayesian hyperparameter tuning, super learner ensembles, and cross-phenotype feature transfer. This repository contains all analysis scripts for the regression and machine learning pipelines.

## Citation

Watts, D., et al. Mapping heterogeneity in the neuroanatomical correlates of depression. *Molecular Psychiatry*. [DOI forthcoming]

## Repository structure

```
watts-depression-neuroanatomical-correlates/
├── ml_scripts/
│   ├── deep_vs_controls/
│   │   └── deep_vs_control_classifier_nonoverlap_control.R
│   ├── feature_transfer/
│   │   ├── bootstrapped_auc_f1_mcnemars_feature_transfer_deep_vs_controls_task.R
│   │   ├── bootstrapped_auc_f1_mcnemars_feature_transfer_shallow_vs_controls_task.R
│   │   ├── deep_vs_control_task_shallow_vs_control_feature_transfer.R
│   │   └── shallow_vs_control_task_deep_vs_shallow_feature_transfer.R
│   ├── generalizability_test/
│   │   └── portability_test.R
│   ├── shallow_vs_controls/
│   │   ├── downsampled_cases_to_deep_task/
│   │   │   └── top30percent_shallow_vs_controls_downsampled_cases_to_match_deep_vs_control_ratio_model.R
│   │   └── shallow_vs_control_rf_classifier_top30percent_removing_overlap_duplicate_controls.R
│   └── super_learners/
│       └── deep_vs_controls_super_learner_residualized_covariates_all_features.R
└── regression_analyses/
    ├── S1_neuroimaging_feature_extraction.R
    ├── S2_calculating_odds_ratios_and_pairwise_associations.R
    ├── S3_sample_characterization.R
    ├── S4_sensitivity_analyses.R
    ├── S4b_sensitivity_analyses_within_cases.R
    └── S5_symptom_dimension_logistic_regression_for_cases.R
```

- `regression_analyses/` - Logistic regression pipelines: feature extraction, association testing, sample characterization, sensitivity analyses, within-cases symptom dimension mapping (S1-S5)
- `ml_scripts/deep_vs_controls/` - Random forest classification for deep depression phenotypes vs. controls
- `ml_scripts/shallow_vs_controls/` - Random forest classification for shallow phenotypes vs. controls, including downsampled variants
- `ml_scripts/feature_transfer/` - Cross-phenotype feature transfer analyses and bootstrapped AUC/F1 comparison
- `ml_scripts/generalizability_test/` - Portability testing across phenotype definitions
- `ml_scripts/super_learners/` - Super learner ensemble classification

## Data availability

Analyses use UK Biobank data (Application 32568). Data are available to approved researchers via the UK Biobank access process: https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access

## Dependencies

R packages used across scripts:

boot, Boruta, car, caret, data.table, DiceKriging, DMwR2, dplyr, furrr, future, future.apply, ggplot2, ggrepel, glmnet, gridExtra, here, kknn, lightgbm, limma, magrittr, MASS, mlr3, mlr3extralearners, mlr3learners, mlr3mbo, mlr3pipelines, mlr3torch, mlr3tuning, mlr3viz, openxlsx, parallelly, pbapply, PerfMeas, pheatmap, pROC, progressr, PRROC, ranger, RColorBrewer, readr, readxl, remotes, rgenoud, rsample, shapviz, stringr, themis, tidymodels, tidyr, torch, vctrs, xgboost, yardstick

## Contact

Devon Watts, dwatts1@mgh.harvard.edu
