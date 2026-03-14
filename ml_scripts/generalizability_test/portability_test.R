# Direct model portability test
# Reviewer 2.6: cross-definition generalizability
#
# Tests whether the decision boundary learned from one phenotype definition
# generalizes to the other by applying each trained classifier directly to
# the opposite task's test set without retraining.
#
# Terminology:
#   deep_model    = classifier for deep depression vs controls (76 RFE features)
#   shallow_model = classifier for shallow depression vs controls (76 RFE features)
#   "ported"      = a trained model applied to the opposite task's test set

# --- User-defined paths (modify for your environment) ---
DEEP_MODEL_WORKSPACE_PATH = "/path/to/your/data/deep_vs_controls_residualized_removing_dup_control_deep_cases_top30.RData"  # Deep depression classifier workspace
SHALLOW_MODEL_WORKSPACE_PATH = "/path/to/your/data/shallow_vs_controls_rf_rfe_removing_dup_control_deep_cases.RData"  # Shallow depression classifier workspace

options(timeout = max(10000, getOption("timeout")))

package_install = function(package_list) {
  for (pkg in package_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      tryCatch({
        install.packages(pkg)
      }, error = function(e) {
        message(paste0("Could not install ", pkg, ": ", e))
      })
    }
    tryCatch({
      library(pkg, character.only = TRUE)
    }, error = function(e) {
      message(paste0("Could not load ", pkg, ": ", e))  
    })
  }
}

packages = c('dplyr', 'data.table', 'magrittr', 'rsample', 'readr', 'stringr', 
             'car', 'future', 'future.apply', 'parallel', 'furrr', 'ggplot2', 'pbapply', 
             'progressr', 'gridExtra', 'pbapply', 'vctrs', 'boot', 'themis', 'rgenoud',
             'tidymodels', 'mlr3', 'mlr3tuning', 'mlr3learners', 'mlr3pipelines', 'limma',
             'mlr3mbo', 'ranger', 'parallel', 'future', 'future.apply', 'glmnet', 
             'pROC', 'RColorBrewer', 'remotes', 'lightgbm', 'xgboost', 'shapviz', 'readxl',
             'PerfMeas', 'caret')

package_install(packages)
library(here)

output_dir = here("output", "generalizability_test")

# Step 1: Load workspaces and extract objects
cat("--- STEP 1: Loading workspaces ---\n\n")

deep_env = new.env()
shallow_env = new.env()

cat("Loading deep depression classifier workspace...\n")
load(DEEP_MODEL_WORKSPACE_PATH,
     envir = deep_env)

cat("Loading shallow depression classifier workspace...\n")
load(SHALLOW_MODEL_WORKSPACE_PATH,
     envir = shallow_env)

# Extract trained models
# deep_model: ranger random forest trained to classify deep depression vs controls
# shallow_model: ranger random forest trained to classify shallow depression vs controls
deep_model = deep_env$auc_model
shallow_model = shallow_env$auc_model

# Extract test data
# Full test sets (256 features each) — needed because each model's 76 features
# are a different subset, and we need all features available for cross-task prediction
deep_test_data = deep_env$test_data
shallow_test_data = shallow_env$test_data

# RFE-filtered test sets (76 features each) — used for same-task baseline predictions
deep_testing_rfe = deep_env$testing_rfe
shallow_testing_rfe = shallow_env$testing_rfe

deep_features = deep_model$forest$independent.variable.names
shallow_features = shallow_model$forest$independent.variable.names

cat("Deep depression classifier: ", length(deep_features), " features\n")
cat("Shallow depression classifier: ", length(shallow_features), " features\n")
cat("Shared features between classifiers: ", length(intersect(deep_features, shallow_features)), "\n")

# Step 2: Verify feature alignment
cat("\n--- STEP 2: Feature alignment verification ---\n\n")

# For cross-task prediction, each model's 76 features must be present
# in the opposite task's full test_data (256 features)
deep_in_shallow = all(deep_features %in% names(shallow_test_data))
cat("All deep classifier features present in shallow test_data: ", deep_in_shallow, "\n")

shallow_in_deep = all(shallow_features %in% names(deep_test_data))
cat("All shallow classifier features present in deep test_data: ", shallow_in_deep, "\n")

if (!deep_in_shallow | !shallow_in_deep) {
  stop("CRITICAL: Feature mismatch. Cannot proceed with portability test.")
}

# Verify outcome factor levels
cat("\nDeep test outcome levels: ", levels(deep_test_data$outcome), "\n")
cat("Deep test outcome distribution:\n")
print(table(deep_test_data$outcome))
cat("\nShallow test outcome levels: ", levels(shallow_test_data$outcome), "\n")
cat("Shallow test outcome distribution:\n")
print(table(shallow_test_data$outcome))

# Step 3: Apply models to opposite-task test sets
cat("\n--- STEP 3: Cross-task predictions ---\n\n")

# --- Deep depression classifier applied to shallow test set (ported) ---
cat("Applying deep depression classifier to shallow test set...\n")
# Extract the 76 features the deep model expects from the shallow full test data
shallow_test_with_deep_features = shallow_test_data[, ..deep_features]
cat("  Prepared data: ", nrow(shallow_test_with_deep_features), " subjects x ",
    ncol(shallow_test_with_deep_features), " features\n")

set.seed(42)
pred_deep_on_shallow = predict(deep_model, data = shallow_test_with_deep_features, type = "response")
# Column 2 = P(deep_depression), the case probability from the deep classifier
prob_deep_on_shallow = pred_deep_on_shallow$predictions[, 2]
cat("  Prediction columns: ", colnames(pred_deep_on_shallow$predictions), "\n")
cat("  Predicted probability range: ", range(prob_deep_on_shallow), "\n")
cat("  Predicted probability mean: ", mean(prob_deep_on_shallow), "\n")

# --- Shallow depression classifier applied to deep test set (ported) ---
cat("\nApplying shallow depression classifier to deep test set...\n")
# Extract the 76 features the shallow model expects from the deep full test data
deep_test_with_shallow_features = deep_test_data[, ..shallow_features]
cat("  Prepared data: ", nrow(deep_test_with_shallow_features), " subjects x ",
    ncol(deep_test_with_shallow_features), " features\n")

set.seed(42)
pred_shallow_on_deep = predict(shallow_model, data = deep_test_with_shallow_features, type = "response")
# Column 2 = P(shallow_depression), the case probability from the shallow classifier
prob_shallow_on_deep = pred_shallow_on_deep$predictions[, 2]
cat("  Prediction columns: ", colnames(pred_shallow_on_deep$predictions), "\n")
cat("  Predicted probability range: ", range(prob_shallow_on_deep), "\n")
cat("  Predicted probability mean: ", mean(prob_shallow_on_deep), "\n")

# --- Outcome-specific baseline predictions (same-task, for comparison) ---
cat("\nRe-deriving outcome-specific baseline predictions...\n")

# Shallow classifier on its own shallow test set (baseline)
set.seed(42)
pred_shallow_on_shallow = predict(shallow_model, data = shallow_testing_rfe[, ..shallow_features], type = "response")
prob_shallow_on_shallow = pred_shallow_on_shallow$predictions[, 2]
cat("  Shallow classifier on shallow test - prob range: ", range(prob_shallow_on_shallow), "\n")

# Deep classifier on its own deep test set (baseline)
# Note: deep_testing_rfe has an extra subject_num column; exclude it
deep_rfe_features = intersect(names(deep_testing_rfe), deep_features)
set.seed(42)
pred_deep_on_deep = predict(deep_model, data = deep_testing_rfe[, ..deep_rfe_features], type = "response")
prob_deep_on_deep = pred_deep_on_deep$predictions[, 2]
cat("  Deep classifier on deep test - prob range: ", range(prob_deep_on_deep), "\n")

# Verify against stored AUCs from the original workspaces
cat("\nVerification against stored ROC objects:\n")
stored_deep_auc = as.numeric(deep_env$roc_obj$auc)
stored_shallow_auc = as.numeric(shallow_env$roc_obj$auc)
cat("  Stored deep classifier AUC: ", stored_deep_auc, "\n")
cat("  Stored shallow classifier AUC: ", stored_shallow_auc, "\n")

# Step 4: Compute AUC with 95% CIs
cat("\n--- STEP 4: AUC computation ---\n\n")

# Convert outcome labels to numeric: 1 = case, 0 = control
deep_true_numeric = ifelse(deep_test_data$outcome == "deep_depression", 1, 0)
shallow_true_numeric = ifelse(shallow_test_data$outcome == "shallow_depression", 1, 0)

cat("Deep test set: ", sum(deep_true_numeric == 1), " cases (deep depression), ",
    sum(deep_true_numeric == 0), " controls\n")
cat("Shallow test set: ", sum(shallow_true_numeric == 1), " cases (shallow depression), ",
    sum(shallow_true_numeric == 0), " controls\n")

# --- Ported model AUCs ---

# Deep classifier on shallow test (ported)
roc_deep_on_shallow = roc(shallow_true_numeric, prob_deep_on_shallow,
                            levels = c(0, 1), direction = "<", quiet = TRUE)
auc_deep_on_shallow = as.numeric(auc(roc_deep_on_shallow))
ci_deep_on_shallow = ci.auc(roc_deep_on_shallow, conf.level = 0.95, method = "delong")

cat("Deep classifier on SHALLOW test (ported): AUC = ", sprintf("%.3f", auc_deep_on_shallow),
    " (95% CI: ", sprintf("%.3f", ci_deep_on_shallow[1]), "-",
    sprintf("%.3f", ci_deep_on_shallow[3]), ")\n")

# Shallow classifier on deep test (ported)
roc_shallow_on_deep = roc(deep_true_numeric, prob_shallow_on_deep,
                            levels = c(0, 1), direction = "<", quiet = TRUE)
auc_shallow_on_deep = as.numeric(auc(roc_shallow_on_deep))
ci_shallow_on_deep = ci.auc(roc_shallow_on_deep, conf.level = 0.95, method = "delong")

cat("Shallow classifier on DEEP test (ported): AUC = ", sprintf("%.3f", auc_shallow_on_deep),
    " (95% CI: ", sprintf("%.3f", ci_shallow_on_deep[1]), "-",
    sprintf("%.3f", ci_shallow_on_deep[3]), ")\n")

# --- Outcome-specific baseline AUCs (recomputed for verification) ---

# Deep classifier on deep test (baseline — same task)
roc_deep_on_deep = roc(deep_true_numeric, prob_deep_on_deep,
                         levels = c(0, 1), direction = "<", quiet = TRUE)
auc_deep_on_deep = as.numeric(auc(roc_deep_on_deep))
ci_deep_on_deep = ci.auc(roc_deep_on_deep, conf.level = 0.95, method = "delong")

# Shallow classifier on shallow test (baseline — same task)
roc_shallow_on_shallow = roc(shallow_true_numeric, prob_shallow_on_shallow,
                               levels = c(0, 1), direction = "<", quiet = TRUE)
auc_shallow_on_shallow = as.numeric(auc(roc_shallow_on_shallow))
ci_shallow_on_shallow = ci.auc(roc_shallow_on_shallow, conf.level = 0.95, method = "delong")

cat("\nDeep classifier on DEEP test (baseline): AUC = ", sprintf("%.3f", auc_deep_on_deep),
    " (95% CI: ", sprintf("%.3f", ci_deep_on_deep[1]), "-",
    sprintf("%.3f", ci_deep_on_deep[3]), ")\n")
cat("Shallow classifier on SHALLOW test (baseline): AUC = ", sprintf("%.3f", auc_shallow_on_shallow),
    " (95% CI: ", sprintf("%.3f", ci_shallow_on_shallow[1]), "-",
    sprintf("%.3f", ci_shallow_on_shallow[3]), ")\n")

cat("\nVerification: Recomputed deep baseline AUC = ", sprintf("%.4f", auc_deep_on_deep),
    " vs stored = ", sprintf("%.4f", stored_deep_auc), "\n")
cat("Verification: Recomputed shallow baseline AUC = ", sprintf("%.4f", auc_shallow_on_shallow),
    " vs stored = ", sprintf("%.4f", stored_shallow_auc), "\n")

# Step 5: DeLong's paired test
cat("\n--- STEP 5: DeLong's paired test for AUC comparison ---\n\n")

# Compare ported vs baseline on the SHALLOW test set
# (deep classifier [ported] vs shallow classifier [baseline], same subjects)
delong_shallow = roc.test(roc_deep_on_shallow, roc_shallow_on_shallow,
                            method = "delong", paired = TRUE)
cat("SHALLOW test set: deep classifier (ported) vs shallow classifier (baseline)\n")
cat("  Deep classifier (ported) AUC: ", sprintf("%.4f", auc_deep_on_shallow), "\n")
cat("  Shallow classifier (baseline) AUC: ", sprintf("%.4f", auc_shallow_on_shallow), "\n")
cat("  DeLong Z: ", sprintf("%.4f", delong_shallow$statistic), "\n")
cat("  DeLong p: ", sprintf("%.6f", delong_shallow$p.value), "\n")

# Compare ported vs baseline on the DEEP test set
# (shallow classifier [ported] vs deep classifier [baseline], same subjects)
delong_deep = roc.test(roc_shallow_on_deep, roc_deep_on_deep,
                         method = "delong", paired = TRUE)
cat("\nDEEP test set: shallow classifier (ported) vs deep classifier (baseline)\n")
cat("  Shallow classifier (ported) AUC: ", sprintf("%.4f", auc_shallow_on_deep), "\n")
cat("  Deep classifier (baseline) AUC: ", sprintf("%.4f", auc_deep_on_deep), "\n")
cat("  DeLong Z: ", sprintf("%.4f", delong_deep$statistic), "\n")
cat("  DeLong p: ", sprintf("%.6f", delong_deep$p.value), "\n")

# FDR correction across the two DeLong p-values
delong_p_values = c(delong_shallow$p.value, delong_deep$p.value)
delong_p_fdr = p.adjust(delong_p_values, method = "fdr")
cat("\nFDR-corrected DeLong p-values:\n")
cat("  Shallow test comparison: ", sprintf("%.6f", delong_p_fdr[1]), "\n")
cat("  Deep test comparison: ", sprintf("%.6f", delong_p_fdr[2]), "\n")

# Step 6: Threshold-dependent metrics at 0.5
cat("\n--- STEP 6: Threshold-dependent metrics (threshold = 0.5) ---\n\n")

compute_threshold_metrics = function(true_labels, predicted_probs, threshold = 0.5) {
  pred_class = ifelse(predicted_probs > threshold, 1, 0)

  TP = sum(true_labels == 1 & pred_class == 1)
  TN = sum(true_labels == 0 & pred_class == 0)
  FP = sum(true_labels == 0 & pred_class == 1)
  FN = sum(true_labels == 1 & pred_class == 0)

  accuracy = (TP + TN) / (TP + TN + FP + FN)
  sensitivity = TP / (TP + FN)
  specificity = TN / (TN + FP)
  ppv = ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  npv = ifelse((TN + FN) > 0, TN / (TN + FN), NA)
  f1 = ifelse((2 * TP + FP + FN) > 0, (2 * TP) / (2 * TP + FP + FN), NA)

  data.frame(
    TP = TP, TN = TN, FP = FP, FN = FN,
    Accuracy = round(accuracy, 4),
    Sensitivity = round(sensitivity, 4),
    Specificity = round(specificity, 4),
    PPV = round(ppv, 4),
    NPV = round(npv, 4),
    F1 = round(f1, 4)
  )
}

# Ported models
metrics_deep_on_shallow = compute_threshold_metrics(shallow_true_numeric, prob_deep_on_shallow)
cat("Deep classifier on SHALLOW test (ported):\n")
print(metrics_deep_on_shallow)

metrics_shallow_on_deep = compute_threshold_metrics(deep_true_numeric, prob_shallow_on_deep)
cat("\nShallow classifier on DEEP test (ported):\n")
print(metrics_shallow_on_deep)

# Baselines
metrics_shallow_on_shallow = compute_threshold_metrics(shallow_true_numeric, prob_shallow_on_shallow)
cat("\nShallow classifier on SHALLOW test (baseline):\n")
print(metrics_shallow_on_shallow)

metrics_deep_on_deep = compute_threshold_metrics(deep_true_numeric, prob_deep_on_deep)
cat("\nDeep classifier on DEEP test (baseline):\n")
print(metrics_deep_on_deep)

# Step 7: Summary table
cat("\n============================================================\n")
cat("SUMMARY TABLE\n")
cat("============================================================\n\n")

summary_table = data.frame(
  Model = c("Deep depression classifier", "Shallow depression classifier",
            "Shallow depression classifier", "Deep depression classifier"),
  Applied_To = c("Shallow test (ported)", "Shallow test (baseline)",
                 "Deep test (ported)", "Deep test (baseline)"),
  AUC = c(sprintf("%.3f", auc_deep_on_shallow),
          sprintf("%.3f", auc_shallow_on_shallow),
          sprintf("%.3f", auc_shallow_on_deep),
          sprintf("%.3f", auc_deep_on_deep)),
  CI_95 = c(sprintf("%.3f-%.3f", ci_deep_on_shallow[1], ci_deep_on_shallow[3]),
            sprintf("%.3f-%.3f", ci_shallow_on_shallow[1], ci_shallow_on_shallow[3]),
            sprintf("%.3f-%.3f", ci_shallow_on_deep[1], ci_shallow_on_deep[3]),
            sprintf("%.3f-%.3f", ci_deep_on_deep[1], ci_deep_on_deep[3])),
  DeLong_p = c("", sprintf("%.6f", delong_shallow$p.value),
               "", sprintf("%.6f", delong_deep$p.value)),
  F1 = c(sprintf("%.3f", metrics_deep_on_shallow$F1),
         sprintf("%.3f", metrics_shallow_on_shallow$F1),
         sprintf("%.3f", metrics_shallow_on_deep$F1),
         sprintf("%.3f", metrics_deep_on_deep$F1)),
  Sensitivity = c(sprintf("%.3f", metrics_deep_on_shallow$Sensitivity),
                  sprintf("%.3f", metrics_shallow_on_shallow$Sensitivity),
                  sprintf("%.3f", metrics_shallow_on_deep$Sensitivity),
                  sprintf("%.3f", metrics_deep_on_deep$Sensitivity)),
  Specificity = c(sprintf("%.3f", metrics_deep_on_shallow$Specificity),
                  sprintf("%.3f", metrics_shallow_on_shallow$Specificity),
                  sprintf("%.3f", metrics_shallow_on_deep$Specificity),
                  sprintf("%.3f", metrics_deep_on_deep$Specificity)),
  PPV = c(sprintf("%.3f", metrics_deep_on_shallow$PPV),
          sprintf("%.3f", metrics_shallow_on_shallow$PPV),
          sprintf("%.3f", metrics_shallow_on_deep$PPV),
          sprintf("%.3f", metrics_deep_on_deep$PPV)),
  NPV = c(sprintf("%.3f", metrics_deep_on_shallow$NPV),
          sprintf("%.3f", metrics_shallow_on_shallow$NPV),
          sprintf("%.3f", metrics_shallow_on_deep$NPV),
          sprintf("%.3f", metrics_deep_on_deep$NPV))
)

print(summary_table, row.names = FALSE)

cat("\n--- DeLong paired test results ---\n\n")
cat("Shallow test set: deep classifier (ported) vs shallow classifier (baseline)\n")
cat("  AUC diff: ", sprintf("%.4f", auc_deep_on_shallow - auc_shallow_on_shallow),
    "  DeLong p = ", sprintf("%.6f", delong_shallow$p.value),
    "  p(FDR) = ", sprintf("%.6f", delong_p_fdr[1]), "\n")

cat("\nDeep test set: shallow classifier (ported) vs deep classifier (baseline)\n")
cat("  AUC diff: ", sprintf("%.4f", auc_shallow_on_deep - auc_deep_on_deep),
    "  DeLong p = ", sprintf("%.6f", delong_deep$p.value),
    "  p(FDR) = ", sprintf("%.6f", delong_p_fdr[2]), "\n")

# Step 8: Save results
cat("\n--- Saving results ---\n")

portability_results = list(
  # AUC results — ported models
  auc_deep_on_shallow = auc_deep_on_shallow,
  ci_deep_on_shallow = ci_deep_on_shallow,
  auc_shallow_on_deep = auc_shallow_on_deep,
  ci_shallow_on_deep = ci_shallow_on_deep,

  # AUC results — baselines (same-task)
  auc_deep_on_deep = auc_deep_on_deep,
  ci_deep_on_deep = ci_deep_on_deep,
  auc_shallow_on_shallow = auc_shallow_on_shallow,
  ci_shallow_on_shallow = ci_shallow_on_shallow,

  # Predicted probabilities
  prob_deep_on_shallow = prob_deep_on_shallow,
  prob_shallow_on_deep = prob_shallow_on_deep,
  prob_deep_on_deep = prob_deep_on_deep,
  prob_shallow_on_shallow = prob_shallow_on_shallow,

  # True labels (numeric: 1 = case, 0 = control)
  shallow_true_numeric = shallow_true_numeric,
  deep_true_numeric = deep_true_numeric,

  # ROC objects
  roc_deep_on_shallow = roc_deep_on_shallow,
  roc_shallow_on_deep = roc_shallow_on_deep,
  roc_deep_on_deep = roc_deep_on_deep,
  roc_shallow_on_shallow = roc_shallow_on_shallow,

  # DeLong paired test results
  delong_shallow = delong_shallow,
  delong_deep = delong_deep,
  delong_p_fdr = delong_p_fdr,

  # Threshold-dependent metrics
  metrics_deep_on_shallow = metrics_deep_on_shallow,
  metrics_shallow_on_deep = metrics_shallow_on_deep,
  metrics_shallow_on_shallow = metrics_shallow_on_shallow,
  metrics_deep_on_deep = metrics_deep_on_deep,

  # Summary table
  summary_table = summary_table,

  # Feature info
  deep_features = deep_features,
  shallow_features = shallow_features,
  shared_features = intersect(deep_features, shallow_features),

  # Metadata
  date_run = Sys.time(),
  seed = 42
)

save(portability_results,
     file = file.path(output_dir, "portability_test_results.RData"))
cat("Saved to: ", file.path(output_dir, "portability_test_results.RData"), "\n")

write.csv(summary_table,
          file = file.path(output_dir, "portability_test_summary.csv"),
          row.names = FALSE)
cat("Summary CSV saved to: ", file.path(output_dir, "portability_test_summary.csv"), "\n")
