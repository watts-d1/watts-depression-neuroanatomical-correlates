# Feature transfer - Deep vs Controls Features applied to Shallow vs Controls

# --- User-defined paths (modify for your environment) ---
DEEP_MODEL_PATH = "/path/to/deep_vs_controls_rf_rfe_80percent_downsampled_residualized_covariates.RData"  # Deep vs controls RFE model workspace
SHALLOW_MODEL_PATH = "/path/to/shallow_vs_controls_rf_rfe_80percent_downsampled_residualized_covariates.RData"  # Shallow vs controls RFE model workspace
LOOKUP_TABLE_PATH = "/path/to/lookup_table_neuroimaging_ukbb.xlsx"  # Neuroimaging feature lookup table

package_install = function(package_list) {
  for (pkg in package_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      tryCatch({
        install.packages(pkg)
      }, error = function(e) {
        message(paste0("Could not install ", pkg, ": ", e))
      })
    }
    library(pkg, character.only = TRUE)
  }
}

packages = c('dplyr', 'data.table', 'magrittr', 'rsample', 'readr', 'stringr', 
             'car', 'future', 'future.apply', 'parallel', 'furrr', 'ggplot2', 'pbapply', 
             'progressr', 'gridExtra', 'pbapply', 'vctrs', 'boot', 'themis', 'rgenoud',
             'DMwR2', 'tidymodels', 'mlr3', 'mlr3tuning', 'mlr3learners', 'mlr3pipelines', 
             'mlr3mbo', 'ranger', 'parallel', 'future', 'future.apply', 'glmnet', 'Boruta', 
             'pROC', 'RColorBrewer', 'remotes', 'lightgbm', 'xgboost', 'shapviz', 'readxl',
             'PerfMeas')

package_install(packages)
library(here)

attach(DEEP_MODEL_PATH)

training_features_deep = training_rfe
testing_features_deep  = testing_rfe
training_data_deep     = train_data_downsampled
testing_data_deep      = test_data
var_importance_df_deep = var_importance_df

detach()

attach(SHALLOW_MODEL_PATH)

training_features_shallow = training_rfe
testing_features_shallow  = testing_rfe
training_data_shallow     = train_data_downsampled
testing_data_shallow      = test_data
var_importance_df_shallow = var_importance_df

detach()

training_data_deep_subset = training_data_deep %>%
  select(all_of(names(training_features_shallow)))

testing_data_deep_subset = testing_data_deep %>%
  select(all_of(names(testing_features_shallow)))

classes_training_downsampled = table(training_data_deep_subset$outcome)
classes_testing_downsampled = table(testing_data_deep_subset$outcome)
print(classes_training_downsampled)
print(classes_testing_downsampled)

# Hyperparameter tuning - Bayesian optimization

bayes_tune = mlr3tuning::tnr("mbo")
task_train = TaskClassif$new("Outcome Train", 
                             backend = training_data_deep_subset, 
                             target = "outcome", 
                             positive = "deep_depression")

search_space_prediction = ps(# mtry: Number of variables randomly sampled as candidates at each split
  mtry = p_int(lower = 1, upper = floor(sqrt(ncol(training_data_deep_subset)))),
  # min.node.size: Minimum size of terminal nodes
  min.node.size = p_int(lower = 5, upper = 20),
  # sample.fraction: Fraction of observations to sample without replacement
  sample.fraction = p_dbl(lower = 0.1, upper = 0.8),
  # splitrule: Splitting rule to use (gini, extratrees, or hellinger)
  # hellinger distance is a measure of convergence b/w probability distributions - useful with imbalanced classes
  splitrule = p_fct(levels = c("gini", "extratrees", "hellinger")),
  # num.trees: Number of trees to grow in the forest
  num.trees = p_int(lower = 2500, upper = 2500),
  # alpha: Significance level for the test statistic in conditional inference trees
  alpha = p_dbl(lower = 0, upper = 1),
  # minprop: Minimum proportion of observations per tree that must reach a terminal node
  minprop = p_dbl(lower = 0, upper = 0.5),
  # max.depth: Maximum depth of any node in the tree
  max.depth = p_int(lower = 1, upper = 20)
)

learner_ranger = set_threads(lrn("classif.ranger"), n = 12)
learner_ranger$predict_type = "prob"

set.seed(62)
tuning_ranger_prediction = mlr3tuning::tune(tuner = bayes_tune,
                                            search_space = search_space_prediction,
                                            task = task_train,
                                            learner = learner_ranger,
                                            resampling = rsmp("cv", folds = 5),
                                            measure = msr("classif.auc"),
                                            term_evals = 25)

best_params_prediction = tuning_ranger_prediction$result_learner_param_vals
print(best_params_prediction)

f_controls = table(training_data_deep_subset$outcome)[1]
f_cases    = table(training_data_deep_subset$outcome)[2]
w_cases    = f_controls / (f_controls + f_cases)
w_controls = f_cases / (f_controls + f_cases)

weights = ifelse(training_data_deep_subset$outcome == '1', w_cases, w_controls)

set.seed(150)
per_model = ranger::ranger(dependent.variable.name = "outcome",
                           data = training_data_deep_subset,
                           mtry = best_params_prediction$mtry,
                           min.node.size   = best_params_prediction$min.node.size,
                           sample.fraction = best_params_prediction$sample.fraction,
                           max.depth       = best_params_prediction$max.depth,
                           minprop         = best_params_prediction$minprop,
                           alpha           = best_params_prediction$alpha,
                           splitrule       = best_params_prediction$splitrule,
                           case.weights    = weights, 
                           num.trees       = best_params_prediction$num.trees,
                           probability     = FALSE)

predictions = predict(per_model, data = testing_data_deep_subset, type = "response")

# Performance metrics

confusion_matrix = table(testing_data_deep_subset$outcome, predictions$predictions)

set.seed(130)
auc_model = ranger::ranger(dependent.variable.name = "outcome",
                           data = training_data_deep_subset,
                           mtry = best_params_prediction$mtry,
                           min.node.size   = best_params_prediction$min.node.size,
                           sample.fraction = best_params_prediction$sample.fraction,
                           max.depth       = best_params_prediction$max.depth,
                           minprop         = best_params_prediction$minprop,
                           alpha           = best_params_prediction$alpha,
                           splitrule       = best_params_prediction$splitrule,
                           case.weights    = weights, 
                           num.trees       = best_params_prediction$num.trees,
                           probability     = TRUE
)

prediction_probabilities = predict(auc_model, data = testing_data_deep_subset, type = "response")

roc_curve = roc(testing_data_deep_subset$outcome, prediction_probabilities$predictions[,2])

TP = confusion_matrix[2, 2]
TN = confusion_matrix[1, 1]
FP = confusion_matrix[1, 2]
FN = confusion_matrix[2, 1]

f1_score = (2 * TP) / (2 * TP + FP + FN)
cat(paste0("F1 Score: ", f1_score, "\n"))

accuracy = (TP + TN) / (TP + TN + FP + FN)
sensitivity = TP / (TP + FN)
specificity = TN / (TN + FP)
ppv = TP / (TP + FP)
npv = TN / (TN + FN)

cat(paste0("Accuracy: ", accuracy, "\n"))
cat(paste0("Sensitivity: ", sensitivity, "\n"))
cat(paste0("Specificity: ", specificity, "\n"))
cat(paste0("Positive Predictive Value: ", ppv, "\n"))
cat(paste0("Negative Predictive Value: ", npv))

accuracy_ci = binom.test(x = TP + TN,
                         n = TP + TN + FP + FN,
                         conf.level = 0.95)$conf.int

sensitivity_ci = prop.test(x = confusion_matrix[2,2],
                           n = sum(confusion_matrix[2,]),
                           conf.level = 0.95)$conf.int

specificity_ci = prop.test(x = confusion_matrix[1,1],
                           n = sum(confusion_matrix[1,]),
                           conf.level = 0.95)$conf.int

ppv_ci = prop.test(x = TP,
                   n = TP + FP,
                   conf.level = 0.95)$conf.int

npv_ci = prop.test(x = TN,
                   n = TN + FN,
                   conf.level = 0.95)$conf.int

cat(paste0("Accuracy: ", round(accuracy * 100, 2), "% (", round(accuracy_ci[1] * 100, 2), "% - ", round(accuracy_ci[2] * 100, 2), "%)\n"))

cat(paste0("Sensitivity: ", round(sensitivity * 100, 2), "% (", round(sensitivity_ci[1] * 100, 2), "% - ", round(sensitivity_ci[2] * 100, 2), "%)\n"))
cat(paste0("Specificity: ", round(specificity * 100, 2), "% (", round(specificity_ci[1] * 100, 2), "% - ", round(specificity_ci[2] * 100, 2), "%)\n"))

cat(paste0("Positive Predictive Value: ", round(ppv * 100, 2), "% (", round(ppv_ci[1] * 100, 2), "% - ", round(ppv_ci[2] * 100, 2), "%)\n"))
cat(paste0("Negative Predictive Value: ", round(npv * 100, 2), "% (", round(npv_ci[1] * 100, 2), "% - ", round(npv_ci[2] * 100, 2), "%)\n"))

predictions_auc = predict(auc_model, data = testing_data_deep_subset, type = "response") 
predicted_probabilities_auc = predictions_auc$predictions[, 2]

roc_obj = roc(testing_data_deep_subset$outcome, predicted_probabilities_auc)
auc = pROC::auc(roc_obj)

cat(paste0("AUC: ", auc, "\n"))

ci_auc = ci(roc_obj, conf.level = 0.95)
print(ci_auc)

roc_df = data.frame(
  Sensitivity = rev(roc_curve$sensitivities),
  Specificity = rev(1 - roc_curve$specificities)
)

roc_plot = ggplot(roc_df, aes(x = Specificity, y = Sensitivity)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_smooth(color = "darkorange", linewidth = 1, linetype = 2, se = FALSE, method = "loess") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "ROC Curve",
       x = "1 - Specificity",
       y = "Sensitivity") +
  annotate("text", x = 0.6, y = 0.4, label = paste0("AUC: ", round(auc, 3)), size = 6, color = "black") +
  theme(plot.title   = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(roc_plot)

# PR-AUC calculation

binarized_predictions = ifelse(predicted_probabilities_auc >= 0.5, 1, 0)

true_labels = testing_data_deep_subset$outcome

true_labels_numeric = as.numeric(true_labels) - 1

res = precision.at.all.recall.levels(predicted_probabilities_auc, true_labels_numeric)

auprc = AUPRC(list(res), comp.precision = TRUE)

cat(paste0("AUPRC: ", round(auprc, 4), "\n"))

# Bootstrap 95% CI for AUPRC
n_bootstraps = 10000
bootstrapped_auprc = numeric(n_bootstraps)

for (i in 1:n_bootstraps) {
  bootstrap_indices = sample(1:length(true_labels_numeric), replace = TRUE)
  bootstrap_predicted_probs = predicted_probabilities_auc[bootstrap_indices]
  bootstrap_true_labels = true_labels_numeric[bootstrap_indices]
  
  bootstrap_res = precision.at.all.recall.levels(bootstrap_predicted_probs, bootstrap_true_labels)
  bootstrapped_auprc[i] = AUPRC(list(bootstrap_res), comp.precision = TRUE)
}

ci_auprc = quantile(bootstrapped_auprc, c(0.025, 0.975))

cat(paste0("95% CI for AUPRC: (", round(ci_auprc[1], 3), ", ", round(ci_auprc[2], 3), ")\n"))

# Performance metrics at 80% sensitivity threshold
desired_sensitivity = 0.80

sens_thresholds = data.frame(
  sensitivity = roc_curve$sensitivities,
  threshold = roc_curve$thresholds
)

optimal_threshold = sens_thresholds$threshold[which.min(abs(sens_thresholds$sensitivity - desired_sensitivity))]

adjusted_predictions_sens = as.numeric(prediction_probabilities$predictions[,2] > optimal_threshold)

adjusted_confusion_matrix_sens = table(testing_data_deep_subset$outcome, adjusted_predictions_sens)

TP_sens = adjusted_confusion_matrix_sens[2, 2]
TN_sens = adjusted_confusion_matrix_sens[1, 1]
FP_sens = adjusted_confusion_matrix_sens[1, 2]
FN_sens = adjusted_confusion_matrix_sens[2, 1]

f1_score_sens = (2 * TP_sens) / (2 * TP_sens + FP_sens + FN_sens)

cat(paste0("F1 Score (Sensitivity 80%): ", f1_score_sens, "\n"))

# Bootstrap 95% CI for F1 score

calculate_f1_from_confusion_sens = function(confusion_matrix) {
  TP_sens = confusion_matrix[2, 2]
  FP_sens = confusion_matrix[1, 2]
  FN_sens = confusion_matrix[2, 1]
  
  f1_score_sens = (2 * TP_sens) / (2 * TP_sens + FP_sens + FN_sens)
  return(f1_score_sens)
}

n_bootstraps = 10000
bootstrapped_f1_sens = numeric(n_bootstraps)

for (i in 1:n_bootstraps) {
  bootstrap_indices = sample(1:length(testing_data_deep_subset$outcome), replace = TRUE)
  bootstrap_true_labels = testing_data_deep_subset$outcome[bootstrap_indices]
  bootstrap_predicted_labels = adjusted_predictions_sens[bootstrap_indices]
  
  bootstrap_confusion = table(bootstrap_true_labels, bootstrap_predicted_labels)
  
  bootstrapped_f1_sens[i] = calculate_f1_from_confusion_sens(bootstrap_confusion)
}

ci_f1_sens = quantile(bootstrapped_f1_sens, c(0.025, 0.975))

cat(paste0("95% CI for F1 Score (Sensitivity 80%): (", round(ci_f1_sens[1], 4), ", ", round(ci_f1_sens[2], 4), ")\n"))

accuracy_sens = (TP_sens + TN_sens) / (TP_sens + TN_sens + FP_sens + FN_sens)

sensitivity_sens = TP_sens / (TP_sens + FN_sens)

specificity_sens = TN_sens / (TN_sens + FP_sens)

ppv_sens = TP_sens / (TP_sens + FP_sens)

npv_sens = TN_sens / (TN_sens + FN_sens)

cat(paste0("Confusion Matrix (Sensitivity 80%):\n"))
print(adjusted_confusion_matrix_sens)
cat(paste0("Accuracy (Sensitivity 80%): ", accuracy_sens, "\n"))
cat(paste0("Sensitivity (Sensitivity 80%): ", sensitivity_sens, "\n"))
cat(paste0("Specificity (Sensitivity 80%): ", specificity_sens, "\n"))
cat(paste0("Positive Predictive Value (Sensitivity 80%): ", ppv_sens, "\n"))
cat(paste0("Negative Predictive Value (Sensitivity 80%): ", npv_sens, "\n"))

accuracy_ci_sens = binom.test(x = TP_sens + TN_sens,
                              n = TP_sens + TN_sens + FP_sens + FN_sens,
                              conf.level = 0.95)$conf.int

sensitivity_ci_sens = prop.test(x = TP_sens,
                                n = TP_sens + FN_sens,
                                conf.level = 0.95)$conf.int

specificity_ci_sens = prop.test(x = TN_sens,
                                n = TN_sens + FP_sens,
                                conf.level = 0.95)$conf.int

ppv_ci_sens = prop.test(x = TP_sens,
                        n = TP_sens + FP_sens,
                        conf.level = 0.95)$conf.int

npv_ci_sens = prop.test(x = TN_sens,
                        n = TN_sens + FN_sens,
                        conf.level = 0.95)$conf.int

cat(paste0("Accuracy (Sensitivity 80%): ", round(accuracy_sens * 100, 2), "% (", round(accuracy_ci_sens[1] * 100, 2), "% - ", round(accuracy_ci_sens[2] * 100, 2), "%)\n"))

cat(paste0("Sensitivity (Sensitivity 80%): ", round(sensitivity_sens * 100, 2), "% (", round(sensitivity_ci_sens[1] * 100, 2), "% - ", round(sensitivity_ci_sens[2] * 100, 2), "%)\n"))
cat(paste0("Specificity (Sensitivity 80%): ", round(specificity_sens * 100, 2), "% (", round(specificity_ci_sens[1] * 100, 2), "% - ", round(specificity_ci_sens[2] * 100, 2), "%)\n"))

cat(paste0("Positive Predictive Value (Sensitivity 80%): ", round(ppv_sens * 100, 2), "% (", round(ppv_ci_sens[1] * 100, 2), "% - ", round(ppv_ci_sens[2] * 100, 2), "%)\n"))
cat(paste0("Negative Predictive Value (Sensitivity 80%): ", round(npv_sens * 100, 2), "% (", round(npv_ci_sens[1] * 100, 2), "% - ", round(npv_ci_sens[2] * 100, 2), "%)\n"))

# Variable importance - permutation importance

set.seed(140)
var_model = ranger::ranger(dependent.variable.name = "outcome",
                           data = training_data_deep_subset,
                           mtry = best_params_prediction$mtry,
                           min.node.size   = best_params_prediction$min.node.size,
                           sample.fraction = best_params_prediction$sample.fraction,
                           minprop         = best_params_prediction$minprop,
                           alpha           = best_params_prediction$alpha,
                           num.trees       = best_params_prediction$num.trees,
                           splitrule       = best_params_prediction$splitrule,
                           max.depth       = best_params_prediction$max.depth,
                           case.weights    = weights, 
                           importance      = "permutation")

var_importance = var_model$variable.importance

var_importance_df = data.frame(Variable = names(var_importance),
                               Importance = var_importance)

var_importance_df = var_importance_df[order(-var_importance_df$Importance),]
var_importance_df$Variable = factor(var_importance_df$Variable, levels = var_importance_df$Variable)

top_20_df = var_importance_df[1:20,]
top_20_df$Variable = factor(top_20_df$Variable, levels = top_20_df$Variable)

palette = brewer.pal(10, "Set3")

top_20_plot = ggplot(top_20_df, aes(x = Variable, y = Importance, fill = Importance)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_distiller(palette = "Set3") +
  theme_minimal() +
  labs(title = "Top 20 Variable Importance - Deep Depression vs Controls",
       x = "Variables",
       y = "Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(10, 10, 30, 10), "pt"))

top_20_plot

bottom_20_df = var_importance_df[(nrow(var_importance_df)-19):nrow(var_importance_df),]
bottom_20_df$Variable = factor(bottom_20_df$Variable, levels = bottom_20_df$Variable)

bottom_20_plot = ggplot(bottom_20_df, aes(x = Variable, y = Importance, fill = Importance)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_distiller(palette = "Set3") +
  theme_minimal() +
  labs(title = "Bottom 20 Variable Importance",
       x = "Variables",
       y = "Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(10, 10, 30, 10), "pt"))

bottom_20_plot

# Map feature IDs to descriptive names

# Map feature IDs to descriptive names via lookup table
lookup_table_path = LOOKUP_TABLE_PATH

read_lookup_table = function(filepath) {
  read_excel(filepath)
}

lookup_table_df = read_lookup_table(lookup_table_path)

lookup_table_df$corrected_id = gsub("baseline_", "bl_", lookup_table_df$corrected_id)

var_importance_lookup = var_importance_df %>%
  left_join(lookup_table_df, by = c("Variable" = "corrected_id")) %>%
  select(Variable, Importance, processed_name_abbreviated)

head(var_importance_lookup)

save(var_importance_df, var_importance_lookup, file = here("output", "feature_transfer", "var_importance_deep_task_shallow_control_features_feature_transfer.RData"))

save(binarized_predictions, true_labels_numeric, f1_score_sens, auc, roc_df, auc_model, roc_obj, ci_auc, ppv_sens, npv_sens, file = here("output", "feature_transfer", "deep_task_shallow_control_features_binary_classifier_significance.RData"))

