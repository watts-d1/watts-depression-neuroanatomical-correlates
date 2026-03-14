# --- User-defined paths (modify for your environment) ---
FEATURE_TRANSFER_RESULTS_PATH = "/path/to/shallow_task_deep_control_features_binary_classifier_significance.RData"  # Feature transfer classifier results
ORIGINAL_MODEL_RESULTS_PATH = "/path/to/shallow_controls_binary_classifier_rfe_significance_residualized_covariates.RData"  # Original shallow vs controls classifier results

library(pROC)
library(boot)
library(yardstick)
library(stats)
library(here)

# Feature Transfer Objects
load(FEATURE_TRANSFER_RESULTS_PATH)

# Rename the objects to avoid confusion
roc_obj_feature_transfer = roc_obj
auc_feature_transfer = auc
ci_feature_transfer = ci_auc
f1_score_feature_transfer = f1_score_sens
precision_feature_transfer = ppv_sens
recall_feature_transfer = npv_sens
true_labels_feature_transfer = true_labels_numeric

rm(roc_obj, auc, ci_auc, f1_score_sens, ppv_sens, npv_sens, true_labels_numeric)

# Outcome_specific objects
load(ORIGINAL_MODEL_RESULTS_PATH)

roc_obj_data_driven = roc_obj
auc_data_driven = auc
ci_data_driven = ci_auc
f1_score_data_driven = f1_score_sens
precision_data_driven = ppv_sens
recall_data_driven = npv_sens
true_labels_data_driven = true_labels_numeric

rm(roc_obj, auc, ci_auc, f1_score_sens, ppv_sens, npv_sens, true_labels_numeric)

# Bootstrap test for F1 scores
set.seed(123)

# Ensure true labels and predicted probabilities are aligned
scores_df_feature_transfer = data.frame(
  predictions = roc_obj_feature_transfer$predictor,
  true_labels = true_labels_feature_transfer
)

scores_df_data_driven = data.frame(
  predictions = roc_obj_data_driven$predictor,
  true_labels = true_labels_data_driven
)

f1_score_func = function(data, indices) {
  d = data[indices, ]
  pred_labels = ifelse(d$predictions > 0.5, 1, 0)
  
  tp = sum(d$true_labels == 1 & pred_labels == 1)
  fp = sum(d$true_labels == 0 & pred_labels == 1)
  fn = sum(d$true_labels == 1 & pred_labels == 0)
  
  precision = tp / (tp + fp)
  recall = tp / (tp + fn)
  
  f1 = 2 * (precision * recall) / (precision + recall)
  return(f1)
}

set.seed(163)

combined_scores = rbind(
  cbind(scores_df_feature_transfer, model = "Feature Transfer"),
  cbind(scores_df_data_driven, model = "Data Driven")
)

stratified_bootstrap_test = function(true_labels_1, predicted_probs_1, 
                                     true_labels_2, predicted_probs_2, 
                                     n_iterations = 10000) {
  calculate_f1 = function(y_true, y_pred) {
    confusion = table(y_true, y_pred)
    precision = confusion[2,2] / sum(confusion[,2])
    recall = confusion[2,2] / sum(confusion[2,])
    f1 = 2 * (precision * recall) / (precision + recall)
    return(f1)
  }
  
  # Original metrics
  original_roc_1 = roc(true_labels_1, predicted_probs_1)
  original_auc_1 = auc(original_roc_1)
  original_f1_1 = calculate_f1(true_labels_1, ifelse(predicted_probs_1 > 0.5, 1, 0))
  
  original_roc_2 = roc(true_labels_2, predicted_probs_2)
  original_auc_2 = auc(original_roc_2)
  original_f1_2 = calculate_f1(true_labels_2, ifelse(predicted_probs_2 > 0.5, 1, 0))
  
  # Observed differences
  observed_diff_auc = original_auc_1 - original_auc_2
  observed_diff_f1 = original_f1_1 - original_f1_2
  
  data_1 = data.frame(true_labels = true_labels_1, predicted_probs = predicted_probs_1)
  data_2 = data.frame(true_labels = true_labels_2, predicted_probs = predicted_probs_2)
  class_distribution_1 = table(true_labels_1) / length(true_labels_1)
  class_distribution_2 = table(true_labels_2) / length(true_labels_2)
  
  # Bootstrap
  bootstrap_results = replicate(n_iterations, {
    # Stratified sampling for both datasets
    boot_data_1 = do.call(rbind, lapply(names(class_distribution_1), function(label) {
      label_data = data_1[data_1$true_labels == label,]
      sample_size = round(nrow(data_1) * class_distribution_1[label])
      label_data[sample(nrow(label_data), size = sample_size, replace = TRUE),]
    }))
    
    boot_data_2 = do.call(rbind, lapply(names(class_distribution_2), function(label) {
      label_data = data_2[data_2$true_labels == label,]
      sample_size = round(nrow(data_2) * class_distribution_2[label])
      label_data[sample(nrow(label_data), size = sample_size, replace = TRUE),]
    }))
    
    boot_roc_1 = roc(boot_data_1$true_labels, boot_data_1$predicted_probs)
    boot_auc_1 = auc(boot_roc_1)
    boot_f1_1 = calculate_f1(boot_data_1$true_labels, ifelse(boot_data_1$predicted_probs > 0.5, 1, 0))
    
    boot_roc_2 = roc(boot_data_2$true_labels, boot_data_2$predicted_probs)
    boot_auc_2 = auc(boot_roc_2)
    boot_f1_2 = calculate_f1(boot_data_2$true_labels, ifelse(boot_data_2$predicted_probs > 0.5, 1, 0))
    
    c(auc_diff = boot_auc_1 - boot_auc_2, f1_diff = boot_f1_1 - boot_f1_2)
  })
  
  p_value_auc = mean(abs(bootstrap_results["auc_diff",]) >= abs(observed_diff_auc))
  p_value_f1 = mean(abs(bootstrap_results["f1_diff",]) >= abs(observed_diff_f1))
  
  permuted_diff_mean_auc = mean(bootstrap_results["auc_diff",])
  permuted_diff_mean_f1 = mean(bootstrap_results["f1_diff",])
  
  ci_auc = quantile(bootstrap_results["auc_diff",], c(0.025, 0.975))
  ci_f1 = quantile(bootstrap_results["f1_diff",], c(0.025, 0.975))
  
  # P-value correction
  p_values = c(p_value_auc, p_value_f1)
  p_values_fdr = p.adjust(p_values, method = "fdr")
  
  list(
    observed_diff_auc = observed_diff_auc,
    observed_diff_f1 = observed_diff_f1,
    permuted_diff_mean_auc = permuted_diff_mean_auc,
    permuted_diff_mean_f1 = permuted_diff_mean_f1,
    ci_auc = ci_auc,
    ci_f1 = ci_f1,
    p_value_auc = p_value_auc,
    p_value_f1 = p_value_f1,
    p_value_auc_fdr = p_values_fdr[1],
    p_value_f1_fdr = p_values_fdr[2],
    bootstrap_auc_diff = bootstrap_results["auc_diff",],
    bootstrap_f1_diff = bootstrap_results["f1_diff",]
  )
}

results = stratified_bootstrap_test(true_labels_feature_transfer, roc_obj_feature_transfer$predictor,
                                    true_labels_data_driven, roc_obj_data_driven$predictor)

print(paste("Observed AUC Difference:", results$observed_diff_auc))
print(paste("Permuted AUC Difference (Mean):", results$permuted_diff_mean_auc))
print(paste("Permuted AUC Difference (95% CI):",results$ci_auc[1], "-", results$ci_auc[2]))
print(paste("AUC P-value:", results$p_value_auc))

print(paste("Observed F1 Difference:", results$observed_diff_f1))
print(paste("Permuted F1 Difference (Mean):", results$permuted_diff_mean_f1))
print(paste("Permuted F1 Difference (95% CI):",   results$ci_f1[1], "-", results$ci_f1[2]))
print(paste("F1 P-value:", results$p_value_f1))

# McNemar's test

comparison_df = data.frame(
  pred_labels_feature_transfer = ifelse(scores_df_feature_transfer$predictions > 0.5, 1, 0),
  pred_labels_data_driven = ifelse(scores_df_data_driven$predictions > 0.5, 1, 0),
  true_labels = scores_df_feature_transfer$true_labels
)

contingency_df = data.frame(
  feature_transfer_correct = comparison_df$pred_labels_feature_transfer == comparison_df$true_labels,
  data_driven_correct = comparison_df$pred_labels_data_driven == comparison_df$true_labels
)

contingency_table = table(contingency_df$feature_transfer_correct, contingency_df$data_driven_correct)
print(contingency_table)

mcnemar_result = mcnemar.test(contingency_table)
cat("McNemar's test results:\n")
print(mcnemar_result)

# FALSE/FALSE (2408): Both models incorrectly classified these instances
# FALSE/TRUE (731): The feature transfer model misclassified, but the data-driven model correctly classified
# TRUE/FALSE (637): The feature transfer model correctly classified, but the data-driven model misclassified
# TRUE/TRUE (26127): Both models correctly classified these instances



