#!/usr/bin/env Rscript
# S4b_sensitivity_within_cases.R
# Models 4 & 5: Within-Cases Analyses (Cases Only)
#
# Model 4: Episode Count (Negative Binomial Regression)
#   For each IDP:
#     episode_count ~ IDP + age + sex + assessment_center + eTIV + head_motion + PC1:PC20
#   Asks: among people with depression, does brain structure relate to recurrence burden?
#
# Model 5: Age at First Episode (Linear Regression)
#   For each IDP:
#     age_at_first_episode ~ IDP + age + sex + assessment_center + eTIV + head_motion + PC1:PC20
#   Asks: among people with depression, is brain structure at imaging associated with earlier versus later illness onset?
#
#   Age at imaging is retained as a covariate so the model tests whether,
#   holding current age constant, earlier onset predicts different brain
#   structure. Without this, onset age confounds with current age (older
#   participants have had more time to develop depression at a younger age).
#
# Samples (both models):
#   - LifetimeMDD cases: Full range. Primary within-cases analysis.
#   - MDDRecur cases: Left-truncated (>=2 episodes / earlier onset typical).
#
# DESIGN RATIONALE FOR REGRESSION CHOICES:
# Model 4 uses negative binomial regression because the outcome (lifetime episode count) is 
# a non-negative integer count variable that is typically right-skewed and overdispersed (variance > mean). 
# Negative binomial handles overdispersion via an additional dispersion parameter (theta).
# If variance/mean is close to 1, the model reduces to Poisson. 
# Linear regression on raw counts would violate normality assumptions; log-transforming counts compresses 
# the low end (log(1)=0) and loses the natural interpretability of count data. 
# Coefficients are exponentiated to incidence rate ratios (IRRs).
#
# Model 5 uses ordinary least squares (linear regression) because age at
# first episode is a continuous variable measured in years and is
# approximately normally distributed, and the assumptions 
# (linearity, normality of residuals, homoscedasticity) are reasonable for this outcome
#  Coefficients are unstandardized betas (years per 1-SD increase in IDP).
# A residual normality check is included in the output.

# MHQ timing exclusion: Fields 20442 and 20433 are from the MHQ.
# Sentinel handling: Negative values in field 20442 set to NA; cases with NA excluded.
#
# Author: Devon Watts

# --- User-defined paths (modify for your environment) ---
UKBB_DATA_PATH   = "/path/to/ukbiobank/ukb674571.csv"
ZSCORE_DATA_PATH = "/path/to/preprocessing_objects/z_score_neuroimaging_objects.RData"

library(data.table)
library(MASS)  # for glm.nb
library(openxlsx)
library(here)

# 0. CONFIGURATION

csv_file_path  = UKBB_DATA_PATH
z_score_path   = ZSCORE_DATA_PATH
output_dir     = here("output", "sensitivity_results")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1. LOAD DATA

cat("Loading z_score objects...\n")
load(z_score_path)

cat("Loading UKB CSV...\n")
csv_data = fread(csv_file_path, showProgress = TRUE, nThread = 20)

# 2. PHENOTYPE REGISTRY (cases-only targets)

phenotype_registry_m4 = list(
  LifetimeMDD = list(
    data  = z_score_completed_LifetimeMDD,
    col   = "LifetimeMDD",
    label = "Lifetime MDD (CIDI-SF)",
    note  = "Full episode range (1+). Primary within-cases analysis."
  ),
  MDDRecur = list(
    data  = z_score_completed_MDDRecurr,
    col   = "MDDRecur",
    label = "Recurrent MDD",
    note  = "Episode count left-truncated at >=2 (definitional). Tests variation among already-recurrent individuals. Narrower question than LifetimeMDD."
  )
)

# 3. EXTRACT EPISODE COUNT AND APPLY EXCLUSIONS

id_col = ifelse("subject_num" %in% names(csv_data), "subject_num", "eid")

# --- 3a. Episode count (field 20442) ---
num_ep_col = "20442-0.0"

if (!num_ep_col %in% names(csv_data)) {
  stop(sprintf("Field %s not found in csv_data. Check column naming.", num_ep_col))
}

episode_dt = csv_data[, c(id_col, num_ep_col), with = FALSE]
setnames(episode_dt, c("subject_num", "episode_count_raw"))
episode_dt[, subject_num := as.integer(subject_num)]

# Recode: negative sentinels to NA
# UKB coding: -1 = "Do not know", -3 = "Prefer not to answer"
episode_dt[, episode_count := {
  x = as.numeric(episode_count_raw)
  x[x < 0] = NA  # Remove ALL negative sentinels
  x
}]

cat("\n=== Episode count distribution (before phenotype restriction) ===\n")
cat("Raw values including sentinels:\n")
print(table(episode_dt$episode_count_raw, useNA = "always"))
cat("\nAfter sentinel removal:\n")
print(summary(episode_dt$episode_count[!is.na(episode_dt$episode_count)]))

episode_dt = episode_dt[, .(subject_num, episode_count)]

# --- 3b. Age at first episode (field 20433) ---
onset_col = "20433-0.0"

if (!onset_col %in% names(csv_data)) {
  stop(sprintf("Field %s not found in csv_data. Check column naming.", onset_col))
}

onset_dt = csv_data[, c(id_col, onset_col), with = FALSE]
setnames(onset_dt, c("subject_num", "onset_age_raw"))
onset_dt[, subject_num := as.integer(subject_num)]

# Recode: negative sentinels to NA
# UKB coding: -1 = "Do not know", -3 = "Prefer not to answer"
onset_dt[, age_at_first_episode := {
  x = as.numeric(onset_age_raw)
  x[x < 0] = NA
  x
}]

cat("\n=== Age at first episode distribution (before phenotype restriction) ===\n")
cat("After sentinel removal:\n")
print(summary(onset_dt$age_at_first_episode[!is.na(onset_dt$age_at_first_episode)]))

onset_dt = onset_dt[, .(subject_num, age_at_first_episode)]

# --- 3c. MHQ timing exclusion ---
mhq_date_col    = "20400-0.0"
imaging_date_col = "53-2.0"

mhq_timing = csv_data[, c(id_col, mhq_date_col, imaging_date_col), with = FALSE]
setnames(mhq_timing, c("subject_num", "mhq_date", "imaging_date"))
mhq_timing[, subject_num := as.integer(subject_num)]

if (is.character(mhq_timing$mhq_date))    mhq_timing[, mhq_date := as.Date(mhq_date)]
if (is.character(mhq_timing$imaging_date)) mhq_timing[, imaging_date := as.Date(imaging_date)]

mri_before_mhq_ids = mhq_timing[!is.na(mhq_date) & !is.na(imaging_date) & imaging_date < mhq_date, subject_num]
mhq_valid_ids      = mhq_timing[!is.na(mhq_date) & !is.na(imaging_date) & mhq_date <= imaging_date, subject_num]

cat(sprintf("\nMHQ timing: %d valid, %d excluded (MRI before MHQ)\n",
            length(mhq_valid_ids), length(mri_before_mhq_ids)))

# 4. COVARIATES (same as primary analysis)

base_covariates = c("bl_31",    # Sex
                    "bl_54",    # Age at imaging
                    "bl_21022", # Assessment center
                    "bl_26521", # eTIV
                    "bl_24419", # Head motion
                    paste0("PC", 1:20))

# 5. NEGATIVE BINOMIAL REGRESSION FUNCTIONS

#' Run negative binomial regression for a single IDP
#'
#' Model: episode_count ~ IDP + covariates
#' Returns exponentiated coefficient (IRR) and profile likelihood CI.
#'
#' @param data data.table with all required columns
#' @param feature Character: column name of the IDP predictor
#' @param outcome_col Character: column name of the episode count outcome
#' @param covariates Character vector of covariate column names
#' @return data.table row with feature, IRR, CI, p-value; or NULL on failure
run_nb_single_feature = function(data, feature, outcome_col, covariates) {
  
  # Check feature exists and has variance
  if (!feature %in% names(data)) return(NULL)
  if (sd(data[[feature]], na.rm = TRUE) == 0) {
    warning(sprintf("Feature %s has zero variance. Skipping.", feature))
    return(NULL)
  }
  
  formula_str = paste(outcome_col, "~", feature, "+", paste(covariates, collapse = " + "))
  
  result = tryCatch({
    model = glm.nb(as.formula(formula_str), data = data)
    
    # Extract coefficient for the IDP feature
    coefs   = summary(model)$coefficients
    if (!feature %in% rownames(coefs)) return(NULL)
    
    beta    = coefs[feature, "Estimate"]
    se      = coefs[feature, "Std. Error"]
    z_val   = coefs[feature, "z value"]
    p_val   = coefs[feature, "Pr(>|z|)"]
    
    # Profile likelihood CIs on log scale
    ci_log = tryCatch({
      ci = confint(model, parm = feature, level = 0.95)
      as.numeric(ci)
    }, error = function(e) {
      # Fallback to Wald CIs if profiling fails
      # Warning: Profile likelihood CI failed for this feature.
      # Falling back to Wald CI (beta +/- 1.96*SE). This may occur for
      # features near boundary conditions or with sparse data.
      cat(sprintf("    [WARNING] Profile CI failed for %s: %s. Using Wald CI.\n",
                  feature, e$message))
      c(beta - 1.96 * se, beta + 1.96 * se)
    })
    
    # Convergence check
    converged = model$converged
    theta     = model$theta  # dispersion parameter
    
    data.table(
      feature    = feature,
      IRR        = exp(beta),
      lower_CI   = exp(ci_log[1]),
      upper_CI   = exp(ci_log[2]),
      beta_log   = beta,
      se_log     = se,
      z_value    = z_val,
      p_value    = p_val,
      theta      = theta,
      converged  = converged
    )
    
  }, error = function(e) {
    cat(sprintf("    [ERROR] NB regression failed for %s: %s\n", feature, e$message))
    return(NULL)
  })
  
  return(result)
}

#' Run negative binomial regression for all IDP features
#'
#' Iterates over all IDP features, computes FDR-corrected p-values,
#' and returns a compiled results table.
#'
#' @param analysis_data data.table with outcome, features, and covariates
#' @param feature_cols Character vector of IDP column names
#' @param outcome_col Character: episode count column name
#' @param covariates Character vector of covariate column names
#' @return data.table with results for all features
run_nb_all_features = function(analysis_data, feature_cols, outcome_col, covariates) {
  
  n_features = length(feature_cols)
  cat(sprintf("  Running negative binomial regression for %d features...\n", n_features))
  
  results_list = vector("list", n_features)
  
  for (i in seq_along(feature_cols)) {
    feature = feature_cols[i]
    
    if (i %% 50 == 0 || i == n_features) {
      cat(sprintf("    Progress: %d / %d features\n", i, n_features))
    }
    
    results_list[[i]] = run_nb_single_feature(analysis_data, feature, outcome_col, covariates)
  }
  
  # Combine results, dropping NULLs
  results_dt = rbindlist(results_list[!sapply(results_list, is.null)])
  
  if (nrow(results_dt) == 0) {
    warning("All features failed. Returning empty data.table.")
    return(data.table())
  }
  
  # FDR correction (within-phenotype, matching primary analysis)
  results_dt[, fdr_p_value       := p.adjust(p_value, method = "fdr")]
  results_dt[, bonferroni_p_value := p.adjust(p_value, method = "bonferroni")]
  results_dt[, significant_fdr   := fdr_p_value < 0.05]
  
  # Check convergence
  n_failed = sum(!results_dt$converged)
  if (n_failed > 0) {
    cat(sprintf("  Warning: %d / %d models did not converge. Review these features.\n",
                n_failed, nrow(results_dt)))
  }
  
  return(results_dt)
}

# 6. BUILD ANALYSIS DATASETS AND RUN

model5_results = list()

for (pheno_name in names(phenotype_registry_m4)) {
  
  reg = phenotype_registry_m4[[pheno_name]]
  dt  = as.data.table(reg$data)
  outcome_col_pheno = reg$col
  
  cat(sprintf("\n{'='*60}\n"))
  cat(sprintf("=== Model 4: %s ===\n", reg$label))
  cat(sprintf("Note: %s\n", reg$note))
  cat(sprintf("{'='*60}\n"))
  
  # --- Step 1: Restrict to CASES only ---
  cases_dt = dt[dt[[outcome_col_pheno]] == 1]
  n_cases_start = nrow(cases_dt)
  cat(sprintf("  Starting cases: n = %d\n", n_cases_start))
  
  # --- Step 2: Merge episode count ---
  cases_dt = merge(cases_dt, episode_dt, by = "subject_num", all.x = TRUE)
  
  # --- Step 3: Apply MHQ timing exclusion ---
  n_before_mhq = nrow(cases_dt)
  cases_dt = cases_dt[!subject_num %in% mri_before_mhq_ids]
  n_after_mhq = nrow(cases_dt)
  cat(sprintf("  After MHQ timing exclusion: n = %d (removed %d)\n",
              n_after_mhq, n_before_mhq - n_after_mhq))
  
  # --- Step 4: Remove cases with missing/invalid episode count ---
  n_before_ep = nrow(cases_dt)
  n_ep_missing = sum(is.na(cases_dt$episode_count))
  cases_dt = cases_dt[!is.na(episode_count)]
  n_after_ep = nrow(cases_dt)
  cat(sprintf("  After removing missing episode count: n = %d (removed %d)\n",
              n_after_ep, n_before_ep - n_after_ep))
  
  pct_total_dropped = 100 * (n_cases_start - n_after_ep) / n_cases_start
  cat(sprintf("  Total attrition from starting cases: %.1f%%\n", pct_total_dropped))
  
  if (pct_total_dropped > 20) {
    cat(sprintf("  Warning: >20%% of %s cases were excluded. Selection bias is possible.\n", pheno_name))
  }
  
  # --- Step 5: Report episode count distribution ---
  cat(sprintf("\n  Episode count distribution among %s cases (n = %d):\n", pheno_name, n_after_ep))
  print(table(cases_dt$episode_count, useNA = "always"))
  cat(sprintf("  Mean = %.1f, Median = %.0f, SD = %.1f, Range = %d-%d\n",
              mean(cases_dt$episode_count), median(cases_dt$episode_count),
              sd(cases_dt$episode_count),
              min(cases_dt$episode_count), max(cases_dt$episode_count)))
  
  # Sanity check for MDDRecur: all episodes should be >= 2
  if (pheno_name == "MDDRecur") {
    n_below_2 = sum(cases_dt$episode_count < 2)
    if (n_below_2 > 0) {
      cat(sprintf("  Warning: %d MDDRecur cases have episode_count < 2. This contradicts the phenotype definition (>=2 episodes). Investigate.\n", n_below_2))
    } else {
      cat("  Sanity check passed: all MDDRecur cases have episode_count >= 2.\n")
    }
  }
  
  # --- Step 6: Check for overdispersion ---
  # Variance >> mean suggests overdispersion, supporting NB over Poisson
  ep_mean = mean(cases_dt$episode_count)
  ep_var  = var(cases_dt$episode_count)
  cat(sprintf("  Overdispersion check: mean = %.2f, variance = %.2f, ratio = %.2f\n",
              ep_mean, ep_var, ep_var / ep_mean))
  if (ep_var / ep_mean < 1.5) {
    cat("  Warning: Variance/mean ratio close to 1. Poisson may be adequate.\n")
    cat("  Proceeding with negative binomial (still valid; reduces to Poisson when theta -> Inf).\n")
  }
  
  # --- Step 7: Identify IDP features ---
  bl_cols      = grep("^bl_", names(cases_dt), value = TRUE)
  feature_cols = setdiff(bl_cols, base_covariates)
  cat(sprintf("  Number of IDP features: %d\n", length(feature_cols)))
  
  # --- Step 8: Complete-case on covariates ---
  analysis_cols = c("episode_count", feature_cols, base_covariates)
  n_before_cc = nrow(cases_dt)
  cases_dt = cases_dt[complete.cases(cases_dt[, ..base_covariates])]
  n_after_cc = nrow(cases_dt)
  if (n_before_cc != n_after_cc) {
    cat(sprintf("  Complete-case on covariates: removed %d, remaining n = %d\n",
                n_before_cc - n_after_cc, n_after_cc))
  }
  
  # --- Step 9: Run negative binomial regression ---
  results = run_nb_all_features(
    analysis_data = cases_dt,
    feature_cols  = feature_cols,
    outcome_col   = "episode_count",
    covariates    = base_covariates
  )
  
  if (nrow(results) > 0) {
    # Add metadata
    results[, phenotype := pheno_name]
    results[, model     := "Model 4: Within-cases NB (episode count)"]
    results[, n_cases   := n_after_cc]
    
    # Classify IDP category from feature name
    results[, idp_category := fcase(
      grepl("^bl_25(0[5-9]|10[0-3])", feature), "FA",
      grepl("^bl_25(1[0-4]|15[01])", feature),  "MD",
      grepl("^bl_2(67[2-5]|68[2-5])", feature), "SA",
      grepl("^bl_2(67[5-8]|68[5-8])", feature), "CT",
      grepl("^bl_250[0-9]$|^bl_2501[0]$", feature), "Global Volume",
      grepl("^bl_265[5-9]|^bl_2659[0-5]", feature), "Subcortical Volume",
      default = "Other"
    )]
    
    # Save
    saveRDS(results, file = file.path(output_dir,
                                      sprintf("model5_NB_%s.RDS", pheno_name)))
    fwrite(results, file = file.path(output_dir,
                                     sprintf("model5_NB_%s.csv", pheno_name)))
    
    # Summary
    n_sig = sum(results$significant_fdr)
    cat(sprintf("\n  Results: %d features tested, %d FDR-significant\n", nrow(results), n_sig))
    
    if (n_sig > 0) {
      cat("  FDR-significant features:\n")
      sig_results = results[significant_fdr == TRUE][order(fdr_p_value)]
      print(sig_results[, .(feature, IRR, lower_CI, upper_CI, p_value, fdr_p_value)])
    }
  }
  
  model5_results[[pheno_name]] = results
}

# 7. MODEL 5: AGE AT FIRST EPISODE (LINEAR REGRESSION)
# For each IDP: age_at_first_episode ~ IDP + covariates
# Linear regression (continuous outcome, approximately normal).
# Coefficient = change in onset age (years) per 1-SD increase in IDP.

cat("=== MODEL 5: Within-Cases Linear Regression (Age at First Episode) ===\n")

#' Run linear regression for a single IDP predicting onset age
run_lm_single_feature = function(data, feature, outcome_col, covariates) {
  if (!feature %in% names(data)) return(NULL)
  if (sd(data[[feature]], na.rm = TRUE) == 0) return(NULL)
  
  formula_str = paste(outcome_col, "~", feature, "+", paste(covariates, collapse = " + "))
  
  result = tryCatch({
    model = lm(as.formula(formula_str), data = data)
    coefs = summary(model)$coefficients
    if (!feature %in% rownames(coefs)) return(NULL)
    
    beta  = coefs[feature, "Estimate"]
    se    = coefs[feature, "Std. Error"]
    t_val = coefs[feature, "t value"]
    p_val = coefs[feature, "Pr(>|t|)"]
    
    ci = confint(model, parm = feature, level = 0.95)
    
    data.table(
      feature      = feature,
      beta         = beta,
      se           = se,
      lower_CI     = ci[1],
      upper_CI     = ci[2],
      t_value      = t_val,
      p_value      = p_val,
      r_squared    = summary(model)$r.squared,
      adj_r_squared = summary(model)$adj.r.squared,
      residual_df  = model$df.residual
    )
  }, error = function(e) {
    cat(sprintf("    [ERROR] LM failed for %s: %s\n", feature, e$message))
    return(NULL)
  })
  return(result)
}

#' Run linear regression for all IDP features
run_lm_all_features = function(analysis_data, feature_cols, outcome_col, covariates) {
  n_features = length(feature_cols)
  cat(sprintf("  Running linear regression for %d features...\n", n_features))
  
  results_list = vector("list", n_features)
  for (i in seq_along(feature_cols)) {
    feature = feature_cols[i]
    if (i %% 50 == 0 || i == n_features) {
      cat(sprintf("    Progress: %d / %d features\n", i, n_features))
    }
    results_list[[i]] = run_lm_single_feature(analysis_data, feature, outcome_col, covariates)
  }
  
  results_dt = rbindlist(results_list[!sapply(results_list, is.null)])
  if (nrow(results_dt) == 0) {
    warning("All features failed. Returning empty data.table.")
    return(data.table())
  }
  
  results_dt[, fdr_p_value       := p.adjust(p_value, method = "fdr")]
  results_dt[, bonferroni_p_value := p.adjust(p_value, method = "bonferroni")]
  results_dt[, significant_fdr   := fdr_p_value < 0.05]
  
  return(results_dt)
}

model6_results = list()

for (pheno_name in names(phenotype_registry_m4)) {
  
  reg = phenotype_registry_m4[[pheno_name]]
  dt  = as.data.table(reg$data)
  outcome_col_pheno = reg$col
  
  cat(sprintf("\n{'='*60}\n"))
  cat(sprintf("=== Model 5: %s ===\n", reg$label))
  cat(sprintf("{'='*60}\n"))
  
  # --- Step 1: Restrict to CASES only ---
  cases_dt = dt[dt[[outcome_col_pheno]] == 1]
  n_cases_start = nrow(cases_dt)
  cat(sprintf("  Starting cases: n = %d\n", n_cases_start))
  
  # --- Step 2: Merge onset age ---
  cases_dt = merge(cases_dt, onset_dt, by = "subject_num", all.x = TRUE)
  
  # --- Step 3: Apply MHQ timing exclusion ---
  n_before_mhq = nrow(cases_dt)
  cases_dt = cases_dt[!subject_num %in% mri_before_mhq_ids]
  n_after_mhq = nrow(cases_dt)
  cat(sprintf("  After MHQ timing exclusion: n = %d (removed %d)\n",
              n_after_mhq, n_before_mhq - n_after_mhq))
  
  # --- Step 4: Remove cases with missing onset age ---
  n_before = nrow(cases_dt)
  n_missing = sum(is.na(cases_dt$age_at_first_episode))
  cases_dt = cases_dt[!is.na(age_at_first_episode)]
  n_after = nrow(cases_dt)
  cat(sprintf("  After removing missing onset age: n = %d (removed %d)\n",
              n_after, n_before - n_after))
  
  pct_total_dropped = 100 * (n_cases_start - n_after) / n_cases_start
  cat(sprintf("  Total attrition from starting cases: %.1f%%\n", pct_total_dropped))
  
  if (pct_total_dropped > 20) {
    cat(sprintf("  Warning: >20%% of %s cases were excluded. Selection bias possible.\n", pheno_name))
  }
  
  # --- Step 5: Report onset age distribution ---
  cat(sprintf("\n  Age at first episode among %s cases (n = %d):\n", pheno_name, n_after))
  print(summary(cases_dt$age_at_first_episode))
  cat(sprintf("  Mean = %.1f, Median = %.0f, SD = %.1f, Range = %.0f-%.0f\n",
              mean(cases_dt$age_at_first_episode), median(cases_dt$age_at_first_episode),
              sd(cases_dt$age_at_first_episode),
              min(cases_dt$age_at_first_episode), max(cases_dt$age_at_first_episode)))
  
  # --- Step 6: Check correlation with age at imaging ---
  cor_val = cor(cases_dt$age_at_first_episode, cases_dt$bl_54, use = "complete.obs")
  cat(sprintf("  Correlation between onset age and age at imaging (z-scored): r = %.3f\n", cor_val))
  if (abs(cor_val) > 0.7) {
    cat("  Warning: High correlation between onset age and age at imaging.\n")
    cat("  Coefficients may be imprecise due to collinearity, but model is still identifiable.\n")
  }
  
  # --- Step 7: Identify IDP features ---
  bl_cols      = grep("^bl_", names(cases_dt), value = TRUE)
  feature_cols = setdiff(bl_cols, base_covariates)
  cat(sprintf("  Number of IDP features: %d\n", length(feature_cols)))
  
  # --- Step 8: Complete-case on covariates ---
  n_before_cc = nrow(cases_dt)
  cases_dt = cases_dt[complete.cases(cases_dt[, ..base_covariates])]
  n_after_cc = nrow(cases_dt)
  if (n_before_cc != n_after_cc) {
    cat(sprintf("  Complete-case on covariates: removed %d, remaining n = %d\n",
                n_before_cc - n_after_cc, n_after_cc))
  }
  
  # --- Step 9: Run linear regression ---
  results = run_lm_all_features(
    analysis_data = cases_dt,
    feature_cols  = feature_cols,
    outcome_col   = "age_at_first_episode",
    covariates    = base_covariates
  )
  
  if (nrow(results) > 0) {
    results[, phenotype := pheno_name]
    results[, model     := "Model 5: Within-cases LM (onset age)"]
    results[, n_cases   := n_after_cc]
    
    saveRDS(results, file = file.path(output_dir,
                                      sprintf("model6_LM_onset_%s.RDS", pheno_name)))
    fwrite(results, file = file.path(output_dir,
                                     sprintf("model6_LM_onset_%s.csv", pheno_name)))
    
    n_sig = sum(results$significant_fdr)
    cat(sprintf("\n  Results: %d features tested, %d FDR-significant\n", nrow(results), n_sig))
    
    if (n_sig > 0) {
      cat("  FDR-significant features:\n")
      sig_results = results[significant_fdr == TRUE][order(fdr_p_value)]
      print(sig_results[, .(feature, beta, lower_CI, upper_CI, p_value, fdr_p_value)])
    }
  }
  
  model6_results[[pheno_name]] = results
}

# 8. COMBINED SUMMARY (MODELS 5 AND 6)

cat("Model 6 summary (Episode Count, Negative Binomial)\n")

for (pheno_name in names(model5_results)) {
  results = model5_results[[pheno_name]]
  reg     = phenotype_registry_m4[[pheno_name]]
  
  cat(sprintf("--- %s (n = %d cases) ---\n", reg$label,
              ifelse(nrow(results) > 0, results$n_cases[1], 0)))
  
  if (nrow(results) == 0) {
    cat("  No results (all features failed).\n\n")
    next
  }
  
  n_sig = sum(results$significant_fdr)
  n_nom = sum(results$p_value < 0.05)
  cat(sprintf("  Features tested: %d\n", nrow(results)))
  cat(sprintf("  Nominally significant (p < 0.05): %d\n", n_nom))
  cat(sprintf("  FDR-significant: %d\n", n_sig))
  cat(sprintf("  Convergence failures: %d\n", sum(!results$converged)))
  
  if (n_sig > 0) {
    n_irr_above_1 = sum(results[significant_fdr == TRUE, IRR > 1])
    n_irr_below_1 = sum(results[significant_fdr == TRUE, IRR < 1])
    cat(sprintf("  Among FDR-significant: %d with IRR > 1 (more episodes), %d with IRR < 1 (fewer episodes)\n",
                n_irr_above_1, n_irr_below_1))
  }
  cat("\n")
}

cat("Model 6 summary (Age at First Episode, Linear Regression)\n")

for (pheno_name in names(model6_results)) {
  results = model6_results[[pheno_name]]
  reg     = phenotype_registry_m4[[pheno_name]]
  
  cat(sprintf("--- %s (n = %d cases) ---\n", reg$label,
              ifelse(nrow(results) > 0, results$n_cases[1], 0)))
  
  if (nrow(results) == 0) {
    cat("  No results (all features failed).\n\n")
    next
  }
  
  n_sig = sum(results$significant_fdr)
  n_nom = sum(results$p_value < 0.05)
  cat(sprintf("  Features tested: %d\n", nrow(results)))
  cat(sprintf("  Nominally significant (p < 0.05): %d\n", n_nom))
  cat(sprintf("  FDR-significant: %d\n", n_sig))
  
  if (n_sig > 0) {
    n_pos = sum(results[significant_fdr == TRUE, beta > 0])
    n_neg = sum(results[significant_fdr == TRUE, beta < 0])
    cat(sprintf("  Among FDR-significant: %d with beta > 0 (later onset), %d with beta < 0 (earlier onset)\n",
                n_pos, n_neg))
  }
  cat("\n")
}

# 9. INTERPRETATION GUIDE

cat("MODEL 5 (Episode Count)\n")
cat("- Negative binomial regression: episode_count ~ IDP + covariates\n")
cat("- IRR > 1: 1-SD increase in IDP associated with MORE episodes.\n")
cat("- IRR < 1: 1-SD increase in IDP associated with FEWER episodes.\n")
cat("- IDPs are z-scored, so coefficients reflect standardized effects.\n")
cat("- MDDRecur is left-truncated at 2 episodes. Range restriction reduces power.\n\n")

cat("MODEL 6 (Age at First Episode)\n")
cat("- Linear regression: age_at_first_episode ~ IDP + covariates\n")
cat("- Beta > 0: 1-SD increase in IDP associated with LATER onset (in years).\n")
cat("- Beta < 0: 1-SD increase in IDP associated with EARLIER onset (in years).\n")
cat("- Age at imaging is in the model, so this tests: holding current age constant,\n")
cat("  does someone with earlier onset have different brain structure?\n")
cat("- This separates onset-timing effects from age-related brain changes.\n\n")

# 10. SUPPLEMENTARY TABLE GENERATION (Models 5-6)
#
# Generates one .xlsx per model with tabs by modality and a Sample_Sizes tab.
# Within each modality tab, phenotypes are stacked (LifetimeMDD, MDDRecur).
#
# Output:
#   Supplementary_Table_Dd_Episode_Count.xlsx   (Model 5, NB)
#   Supplementary_Table_De_Onset_Age.xlsx       (Model 6, LM)
#
cat("\n\n=== Section 11: Generating Supplementary Tables Dd-De ===\n")

# --- Load feature-to-region and modality lookup from S4 ---

lookup_path = file.path(output_dir, "feature_region_lookup.RDS")
if (file.exists(lookup_path)) {
  feature_region_lookup = readRDS(lookup_path)
  cat(sprintf("Loaded feature_region_lookup: %d features\n", nrow(feature_region_lookup)))
} else {
  cat("Warning: feature_region_lookup.RDS not found. Region names will be blank.\n")
  cat("  Run S4_sensitivity_analyses.R Section 9c first to generate this file.\n")
  feature_region_lookup = data.table(feature = character(), region_name = character(),
                                     modality = character())
}

# --- Configuration ---

phenotype_order = c("LifetimeMDD", "MDDRecur")

modality_order = c("FA", "MD", "SA", "CT", "global_volume", "subcortical_volume")

modality_tab_names = c(
  FA                 = "FA",
  MD                 = "MD",
  SA                 = "SA",
  CT                 = "CT",
  global_volume      = "Global_Volume",
  subcortical_volume = "Subcortical_Volume"
)

# --- Styles ---

st_header = createStyle(
  fontName = "Arial", fontSize = 11, textDecoration = "bold",
  fgFill = "#4472C4", fontColour = "#FFFFFF",
  halign = "center", border = "TopBottom", borderColour = "#2F5496"
)

st_pheno_label = createStyle(
  fontName = "Arial", fontSize = 11, textDecoration = "bold",
  fgFill = "#D9E2F3", border = "Bottom", borderColour = "#4472C4"
)

st_sig = createStyle(
  fontName = "Arial", fontSize = 10, textDecoration = "bold"
)

# --- Assign modality to results using lookup ---
# This replaces the regex-based idp_category with the validated
# modality classification from the primary analysis.

assign_modality_from_lookup = function(dt, lookup) {
  dt_out = merge(dt, lookup[, .(feature, region_name, modality)],
                 by = "feature", all.x = TRUE)
  n_missing = sum(is.na(dt_out$modality))
  if (n_missing > 0) {
    cat(sprintf("  Warning: %d features without modality assignment\n", n_missing))
  }
  return(dt_out)
}


# 11a. HELPER: Write one modality tab for within-cases models

write_within_tab = function(wb, tab_name, results_dt, mod_name,
                            phenotypes, model_type) {
  
  addWorksheet(wb, tab_name)
  
  if (model_type == "nb") {
    col_headers = c("Feature", "Region", "IRR", "Lower CI", "Upper CI",
                    "p-value", "FDR p-value", "Theta", "Converged")
  } else {
    col_headers = c("Feature", "Region", "Beta (years)", "SE",
                    "Lower CI", "Upper CI", "p-value", "FDR p-value")
  }
  
  current_row = 1
  
  for (pheno in phenotypes) {
    sub = results_dt[modality == mod_name & phenotype == pheno]
    if (nrow(sub) == 0) next
    
    setorder(sub, feature)
    n_cases = sub$n_cases[1]
    
    # Phenotype label row
    writeData(wb, tab_name,
              x = sprintf("%s (n_cases=%s)", pheno, format(n_cases, big.mark = ",")),
              startRow = current_row, startCol = 1)
    mergeCells(wb, tab_name, cols = 1:length(col_headers), rows = current_row)
    addStyle(wb, tab_name, st_pheno_label, rows = current_row,
             cols = 1:length(col_headers), gridExpand = TRUE)
    current_row = current_row + 1
    
    # Column headers
    writeData(wb, tab_name, x = as.data.frame(t(col_headers)),
              startRow = current_row, startCol = 1, colNames = FALSE)
    addStyle(wb, tab_name, st_header, rows = current_row,
             cols = 1:length(col_headers), gridExpand = TRUE)
    current_row = current_row + 1
    
    # Data rows
    if (model_type == "nb") {
      out_df = data.frame(
        feature   = sub$feature,
        region    = fifelse(is.na(sub$region_name), "", sub$region_name),
        irr       = sub$IRR,
        lower_ci  = sub$lower_CI,
        upper_ci  = sub$upper_CI,
        p_val     = sub$p_value,
        fdr_p     = sub$fdr_p_value,
        theta     = sub$theta,
        converged = fifelse(sub$converged, "TRUE", "FALSE"),
        stringsAsFactors = FALSE
      )
    } else {
      out_df = data.frame(
        feature   = sub$feature,
        region    = fifelse(is.na(sub$region_name), "", sub$region_name),
        beta      = sub$beta,
        se        = sub$se,
        lower_ci  = sub$lower_CI,
        upper_ci  = sub$upper_CI,
        p_val     = sub$p_value,
        fdr_p     = sub$fdr_p_value,
        stringsAsFactors = FALSE
      )
    }
    
    writeData(wb, tab_name, x = out_df,
              startRow = current_row, startCol = 1, colNames = FALSE)
    
    # Bold FDR-significant rows
    sig_rows = which(sub$fdr_p_value < 0.05)
    if (length(sig_rows) > 0) {
      effect_cols = if (model_type == "nb") 3:7 else 3:8
      addStyle(wb, tab_name, st_sig,
               rows = current_row + sig_rows - 1,
               cols = effect_cols, gridExpand = TRUE, stack = TRUE)
    }
    
    current_row = current_row + nrow(sub) + 1  # +1 blank separator
  }
  
  # Column widths
  setColWidths(wb, tab_name, cols = 1, widths = 14)
  setColWidths(wb, tab_name, cols = 2, widths = 55)
  n_remaining = length(col_headers) - 2
  setColWidths(wb, tab_name, cols = 3:(2 + n_remaining), widths = 14)
}


# 11b. HELPER: Write Sample_Sizes tab for within-cases

write_within_ss_tab = function(wb, results_list, model_label, model_type) {
  addWorksheet(wb, "Sample_Sizes")
  
  writeData(wb, "Sample_Sizes",
            x = sprintf("Sample sizes: %s", model_label),
            startRow = 1, startCol = 1)
  addStyle(wb, "Sample_Sizes",
           createStyle(fontSize = 12, textDecoration = "bold"),
           rows = 1, cols = 1)
  
  ss_cols = c("Phenotype", "N Cases", "N Features Tested",
              "N FDR-Significant", "N Nominally Significant (p<0.05)")
  
  if (model_type == "nb") {
    ss_cols = c(ss_cols, "Convergence Failures")
  }
  
  writeData(wb, "Sample_Sizes", x = as.data.frame(t(ss_cols)),
            startRow = 3, startCol = 1, colNames = FALSE)
  addStyle(wb, "Sample_Sizes", st_header, rows = 3,
           cols = 1:length(ss_cols), gridExpand = TRUE)
  
  row_i = 4
  for (pheno in names(results_list)) {
    r = results_list[[pheno]]
    if (is.null(r) || nrow(r) == 0) next
    
    n_cases = r$n_cases[1]
    n_tested = nrow(r)
    n_fdr = sum(r$significant_fdr, na.rm = TRUE)
    n_nom = sum(r$p_value < 0.05, na.rm = TRUE)
    
    row_data = data.frame(
      pheno = pheno,
      n_cases = format(n_cases, big.mark = ","),
      n_tested = n_tested,
      n_fdr = n_fdr,
      n_nom = n_nom,
      stringsAsFactors = FALSE
    )
    
    if (model_type == "nb") {
      n_fail = sum(!r$converged, na.rm = TRUE)
      row_data$n_fail = n_fail
    }
    
    writeData(wb, "Sample_Sizes", x = row_data,
              startRow = row_i, startCol = 1, colNames = FALSE)
    row_i = row_i + 1
  }
  
  setColWidths(wb, "Sample_Sizes", cols = 1, widths = 18)
  setColWidths(wb, "Sample_Sizes", cols = 2:length(ss_cols), widths = 16)
}


# 11c. GENERATE TABLE 4e (Model 4: Episode Count)

cat("\n--- Generating Supplementary_Table_Dd_Episode_Count.xlsx ---\n")

# Combine Model 4 results and assign modality from lookup
m4_all = rbindlist(model5_results, fill = TRUE)
m4_all = assign_modality_from_lookup(m4_all, feature_region_lookup)

# Drop the regex-based idp_category to avoid confusion
if ("idp_category" %in% names(m4_all)) m4_all[, idp_category := NULL]

wb_dd = createWorkbook()

for (mod in modality_order) {
  tab_name = modality_tab_names[mod]
  if (nrow(m4_all[modality == mod]) == 0) next
  
  write_within_tab(
    wb         = wb_dd,
    tab_name   = tab_name,
    results_dt = m4_all,
    mod_name   = mod,
    phenotypes = phenotype_order,
    model_type = "nb"
  )
  cat(sprintf("  %s: written\n", tab_name))
}

write_within_ss_tab(wb_dd, model5_results,
                    "Model 4: Within-cases NB (Episode Count)", "nb")
cat("  Sample_Sizes: written\n")

out_dd = file.path(output_dir, "Supplementary_Table_Dd_Episode_Count.xlsx")
saveWorkbook(wb_dd, out_dd, overwrite = TRUE)
cat(sprintf("  Saved: %s\n", out_dd))


# 11d. GENERATE TABLE 4f (Model 5: Age at First Episode)

cat("\n--- Generating Supplementary_Table_De_Onset_Age.xlsx ---\n")

# Combine Model 5 results and assign modality from lookup
m5_all = rbindlist(model6_results, fill = TRUE)
m5_all = assign_modality_from_lookup(m5_all, feature_region_lookup)

wb_de = createWorkbook()

for (mod in modality_order) {
  tab_name = modality_tab_names[mod]
  if (nrow(m5_all[modality == mod]) == 0) next
  
  write_within_tab(
    wb         = wb_de,
    tab_name   = tab_name,
    results_dt = m5_all,
    mod_name   = mod,
    phenotypes = phenotype_order,
    model_type = "lm"
  )
  cat(sprintf("  %s: written\n", tab_name))
}

write_within_ss_tab(wb_de, model6_results,
                    "Model 5: Within-cases LM (Age at First Episode)", "lm")
cat("  Sample_Sizes: written\n")

out_de = file.path(output_dir, "Supplementary_Table_De_Onset_Age.xlsx")
saveWorkbook(wb_de, out_de, overwrite = TRUE)
cat(sprintf("  Saved: %s\n", out_de))

cat("Done. Proceed to S5\n")
