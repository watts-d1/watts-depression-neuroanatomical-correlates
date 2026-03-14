# --- User-defined paths (modify for your environment) ---
INPUT_DATA_PATH = "/path/to/your/data/non_overlapping_control_df.RData"  # Non-overlapping control dataframe

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
    library(pkg, character.only = TRUE)
  }
}

packages = c('dplyr', 'data.table', 'magrittr', 'rsample', 'readr', 'stringr', 
             'car', 'future', 'future.apply', 'parallel', 'furrr', 'ggplot2', 'pbapply', 
             'progressr', 'gridExtra', 'pbapply', 'vctrs', 'boot', 'themis', 'rgenoud',
             'DMwR2', 'tidymodels', 'mlr3', 'mlr3tuning', 'mlr3learners', 'mlr3pipelines', 
             'mlr3mbo', 'ranger', 'parallel', 'future', 'future.apply', 'glmnet', 'Boruta', 
             'pROC', 'RColorBrewer', 'remotes', 'lightgbm', 'xgboost', 'shapviz', 'readxl',
             'PerfMeas', 'caret', 'DiceKriging', 'kknn', 'mlr3extralearners', 'torch', 
             'mlr3torch', 'mlr3viz', 'PRROC')

package_install(packages)
library(here)

setDTthreads(threads = 10)

load(INPUT_DATA_PATH)

# Shallow depression phenotypes
outcomes_shallow = list(
  GPpsy              = non_overlap_controls_df$cleaned_datasets$GPpsy$GPpsy, 
  Psypsy             = non_overlap_controls_df$cleaned_datasets$Psypsy$Psypsy,
  SelfRepDep         = non_overlap_controls_df$cleaned_datasets$SelfRepDep$SelfRepDep, 
  DepAll             = non_overlap_controls_df$cleaned_datasets$DepAll$DepAll)

# Deep depression phenotypes
outcomes_deep = list(
  ICD10Dep           = non_overlap_controls_df$cleaned_datasets$ICD10Dep$ICD10Dep, 
  ICD10Dep.exclpsych = non_overlap_controls_df$cleaned_datasets$ICD10Dep_exclpsych$ICD10Dep.exclpsych, 
  LifetimeMDD        = non_overlap_controls_df$cleaned_datasets$LifetimeMDD$LifetimeMDD, 
  MDDRecurr          = non_overlap_controls_df$cleaned_datasets$MDDRecur$MDDRecur)

# Phenotype dataframes for shallow depression
shallow_phenotypes = list(GPpsy       = non_overlap_controls_df$cleaned_datasets$GPpsy, 
                          Psypsy      = non_overlap_controls_df$cleaned_datasets$Psypsy,
                          SelfRepDep  = non_overlap_controls_df$cleaned_datasets$SelfRepDep, 
                          DepAll      = non_overlap_controls_df$cleaned_datasets$DepAll)

# Phenotype dataframes for deep depression
deep_phenotypes = list(ICD10Dep           = non_overlap_controls_df$cleaned_datasets$ICD10Dep, 
                       ICD10Dep.exclpsych = non_overlap_controls_df$cleaned_datasets$ICD10Dep_exclpsych, 
                       LifetimeMDD        = non_overlap_controls_df$cleaned_datasets$LifetimeMDD,
                       MDDRecurr          = non_overlap_controls_df$cleaned_datasets$MDDRecur)

#' Modify Outcomes Based on Given Label
#'
#' This function iterates over a list of outcomes and modifies each element
#' based on a given condition: if the element is 0, it is replaced with "control",
#' otherwise, it is replaced with the specified label.
#'
#' @param outcomes A list of numerical outcomes (0 or 1).
#' @param label A string to replace non-zero outcomes.
#' @return A modified list with "control" or the specified label.
#' @examples
#' outcomes = list(example1 = c(0, 1, 0), example2 = c(1, 1, 0))
#' modify_outcomes(outcomes, "positive")
modify_outcomes = function(outcomes, label) {
  sapply(outcomes, function(x) ifelse(x == 0, "control", label), simplify = FALSE)
}

outcomes_shallow_modified = modify_outcomes(outcomes_shallow, "shallow_depression")
outcomes_deep_modified    = modify_outcomes(outcomes_deep, "deep_depression")

# Function to add outcomes to phenotypes
add_outcomes_to_phenotypes = function(phenotypes, modified_outcomes) {
  mapply(function(phenotype, outcome) {
    if (nrow(phenotype) != length(outcome)) {
      stop("Mismatch in number of rows in phenotype and length of outcome")
    }
    phenotype$outcome = outcome
    return(phenotype)
  }, phenotypes, modified_outcomes, SIMPLIFY = FALSE)
}

shallow_phenotypes_modified = add_outcomes_to_phenotypes(shallow_phenotypes, outcomes_shallow_modified)
deep_phenotypes_modified = add_outcomes_to_phenotypes(deep_phenotypes, outcomes_deep_modified)

# Function to remove specified columns and combine dataframes
process_and_combine_dataframes = function(df_list, columns_to_remove) {
  df_list %>% 
    lapply(., function(df) {
      # Ensure only existing columns are selected for removal
      cols_to_remove = intersect(names(df), columns_to_remove)
      select(df, !all_of(cols_to_remove))
    }) %>% 
    bind_rows()
}

original_outcomes_shallow = names(outcomes_shallow)
original_outcomes_deep    = names(outcomes_deep)

combined_shallow_phenotype = process_and_combine_dataframes(shallow_phenotypes_modified, original_outcomes_shallow)
combined_deep_phenotype    = process_and_combine_dataframes(deep_phenotypes_modified, original_outcomes_deep)

## Identifying and removing covariates
covariates = c("bl_31",    # Sex
               "bl_54",    # Assessment Center
               "bl_21022", # Age
               "bl_26521", # Estimated Total Intracranial Volume 
               "bl_24419", # Head motion 
               paste0("PC", 1:20))

covariate_data_deep = combined_deep_phenotype[, ..covariates]
covariate_data_shallow = combined_shallow_phenotype[, ..covariates]

# Function to residualize features
residualize_features = function(data, covariates) {
  residualized_data = data
  covariate_data = data[, ..covariates]
  
  # Get neuroimaging feature names (everything except covariates, outcome, subject_num)
  neuro_features = setdiff(names(data), c(covariates, "outcome", "subject_num"))
  
  for (feature in neuro_features) {
    model = lm(data[[feature]] ~ ., data = covariate_data)
    residualized_data[[feature]] = residuals(model)
  }
  
  return(residualized_data)
}

# Residualize the features in the shallow and deep phenotype dataframes
combined_phenotypes_shallow_residualized = residualize_features(combined_shallow_phenotype, covariates)
combined_phenotypes_deep_residualized = residualize_features(combined_deep_phenotype, covariates)

# Creating a copy of dataframes prior to residualization 
combined_shallow_phenotype_pre_resid = combined_shallow_phenotype
combined_deep_phenotype_pre_resid = combined_deep_phenotype

# Define a function to remove specified columns from a dataframe
remove_covariates = function(df, covariates) {
  # Ensure only existing columns are selected for removal
  cols_to_remove = intersect(names(df), covariates)
  select(df, !all_of(cols_to_remove))
}

combined_shallow_phenotype = remove_covariates(combined_phenotypes_shallow_residualized, covariates)
combined_deep_phenotype = remove_covariates(combined_phenotypes_deep_residualized, covariates)

combined_shallow_phenotype$outcome = as.factor(combined_shallow_phenotype$outcome)
combined_deep_phenotype$outcome = as.factor(combined_deep_phenotype$outcome)

combined_shallow_phenotype$outcome = as.factor(combined_shallow_phenotype$outcome)
combined_deep_phenotype$outcome    = as.factor(combined_deep_phenotype$outcome)

# Print statement to show counts before deduplication for shallow phenotype
print("Counts before deduplication (Shallow Phenotype):")
print(paste("Total rows:", nrow(combined_shallow_phenotype)))
print(table(combined_shallow_phenotype$outcome))

# Check for duplicate subjects in shallow phenotype
shallow_subject_counts = combined_shallow_phenotype %>%
  group_by(subject_num) %>%
  summarise(n = n(), 
            has_multiple_outcomes = n_distinct(outcome) > 1)

print(paste("Shallow - Total unique subjects:", nrow(shallow_subject_counts)))
print(paste("Shallow - Subjects appearing multiple times:", sum(shallow_subject_counts$n > 1)))
print(paste("Shallow - Subjects with conflicting outcomes:", sum(shallow_subject_counts$has_multiple_outcomes)))

# Print statement to show counts before deduplication for deep phenotype
print("Counts before deduplication (Deep Phenotype):")
print(paste("Total rows:", nrow(combined_deep_phenotype)))
print(table(combined_deep_phenotype$outcome))

# Check for duplicate subjects in deep phenotype
deep_subject_counts = combined_deep_phenotype %>%
  group_by(subject_num) %>%
  summarise(n = n(), 
            has_multiple_outcomes = n_distinct(outcome) > 1)

print(paste("Deep - Total unique subjects:", nrow(deep_subject_counts)))
print(paste("Deep - Subjects appearing multiple times:", sum(deep_subject_counts$n > 1)))
print(paste("Deep - Subjects with conflicting outcomes:", sum(deep_subject_counts$has_multiple_outcomes)))

# Deduplicate shallow phenotype
shallow_cases = combined_shallow_phenotype %>% 
  filter(outcome != "control")
shallow_controls = combined_shallow_phenotype %>% 
  filter(outcome == "control")

shallow_cases_unique = shallow_cases %>%
  distinct(subject_num, .keep_all = TRUE)
shallow_controls_unique = shallow_controls %>%
  distinct(subject_num, .keep_all = TRUE)

combined_shallow_phenotype_deduplicated = bind_rows(shallow_cases_unique, shallow_controls_unique) %>%
  arrange(subject_num)

# Deduplicate deep phenotype
deep_cases = combined_deep_phenotype %>% 
  filter(outcome != "control")
deep_controls = combined_deep_phenotype %>% 
  filter(outcome == "control")

deep_cases_unique = deep_cases %>%
  distinct(subject_num, .keep_all = TRUE)
deep_controls_unique = deep_controls %>%
  distinct(subject_num, .keep_all = TRUE)

combined_deep_phenotype_deduplicated = bind_rows(deep_cases_unique, deep_controls_unique) %>%
  arrange(subject_num)

# Print statement to show counts after deduplication
print("Counts after deduplication (Shallow Phenotype):")
print(paste("Total rows:", nrow(combined_shallow_phenotype_deduplicated)))
print(table(combined_shallow_phenotype_deduplicated$outcome))

print("Counts after deduplication (Deep Phenotype):")
print(paste("Total rows:", nrow(combined_deep_phenotype_deduplicated)))
print(table(combined_deep_phenotype_deduplicated$outcome))

# Replace original dataframes with deduplicated versions
combined_shallow_phenotype = combined_shallow_phenotype_deduplicated
combined_deep_phenotype = combined_deep_phenotype_deduplicated

# Now merge deduplicated dataframes
combined_phenotypes = bind_rows(combined_shallow_phenotype, combined_deep_phenotype)

# Check for subjects that appear in both shallow and deep datasets
subjects_in_both = inner_join(
  combined_shallow_phenotype %>% select(subject_num, outcome) %>% rename(shallow_outcome = outcome),
  combined_deep_phenotype %>% select(subject_num, outcome) %>% rename(deep_outcome = outcome),
  by = "subject_num"
)

print(paste("Subjects appearing in both datasets:", nrow(subjects_in_both)))

# Display examples of subjects in both datasets (if any)
if(nrow(subjects_in_both) > 0) {
  print("Examples of subjects in both datasets:")
  print(head(subjects_in_both, 10))
}

# Removing controls from the phenotype definition for a binary classifier (deep vs shallow)
combined_phenotypes_deep_shallow = combined_phenotypes %>% filter(outcome != "control")
combined_phenotypes_deep_shallow$outcome = droplevels(combined_phenotypes_deep_shallow$outcome)

# Final verification after filtering
print("Final binary classifier dataset:")
print(paste("Total rows:", nrow(combined_phenotypes_deep_shallow)))
print(table(combined_phenotypes_deep_shallow$outcome))

# Stratified sampling 
set.seed(404)  
data_split = initial_split(combined_phenotypes_deep_shallow, prop = 0.7, strata = "outcome")
train_data = training(data_split)
test_data  = testing(data_split)

train_data = train_data %>% select(-subject_num, -MDDRecur)
test_data  = test_data %>% select(-subject_num, -MDDRecur)

classes_training = table(train_data$outcome)
print(classes_training)

#' Downsample majority class by a specified percentage
#'
#' @param data A data frame containing the data to be downsampled (training set)
#' @param outcome_name The name of the outcome variable (e.g. LifetimeDepress)
#' @param downsample_percent The percentage of the majority class to retain after downsampling
#' 
#' @return A data frame where the majority class has been downsampled
downsample_majority_class = function(data, outcome_name, downsample_percent) {
  
  # Error handling: Check if outcome_name is a valid column
  if (!outcome_name %in% names(data)) {
    stop("The specified outcome_name is not a column in the provided data frame.")
  }
  
  # Identify unique classes
  classes = unique(data[[outcome_name]])
  
  # Determine the majority class
  class_counts = table(data[[outcome_name]])
  majority_class_label = names(which.max(class_counts))
  
  # Initialize an empty data frame to store the downsampled data
  result = data.frame()
  
  # Loop over each class
  for (class in classes) {
    # Extract data corresponding to the current class
    class_data = data[data[[outcome_name]] == class, ]
    
    # If the current class is the majority class, downsample it
    if (as.character(class) == majority_class_label) {
      
      class_sample_size = floor(nrow(class_data) * downsample_percent)
      class_data = class_data[sample(nrow(class_data), class_sample_size), ]
    }
    
    # Combine the downsampled class data back into the result data frame
    result = rbind(result, class_data)
  }
  
  return(result)
}

# Specify the percentage of majority class to retain (e.g. 0.30 retains 30% of the sample)
downsample_percent = 0.50

# Apply downsampling to the training data
set.seed(62)
train_data_downsampled = downsample_majority_class(train_data, "outcome", downsample_percent)

classes_training_downsampled = table(train_data_downsampled$outcome)
classes_training = table(train_data$outcome)
classes_testing  = table(test_data$outcome)

print(classes_training)
print(classes_training_downsampled)
print(classes_testing)

# Hyperparameter tuning - Bayesian optimization

bayes_tune = mlr3tuning::tnr("mbo")

task_train = TaskClassif$new("Deep Depression vs Shallow Train", 
                             backend = train_data_downsampled, 
                             target = "outcome",
                             positive = "deep_depression")

task_test  = TaskClassif$new("Deep Depression vs Shallow Test", 
                             backend = test_data, 
                             target = "outcome",
                             positive = "deep_depression")

# Base learners


# Elastic net
learner_glmnet = set_threads(lrn("classif.glmnet"), n = 14)

learner_glmnet$predict_type = "prob"

search_space_glmnet = ps(alpha = p_dbl(lower = 0, upper = 1),
                         lambda = p_dbl(lower = 0.01, upper = 1))

#learner_glmnet = custom_cv_glmnet
set.seed(120)
tuning_glmnet = mlr3tuning::tune(tuner = bayes_tune,
                                 search_space = search_space_glmnet,
                                 task = task_train,
                                 learner = learner_glmnet,
                                 resampling = rsmp("cv", folds = 5),
                                 measure = msr("classif.auc"),
                                 term_evals = 25
)

# K-nearest neighbors
learner_knn = set_threads(lrn("classif.kknn"), n = 14)
learner_knn$predict_type = "prob"

search_space_knn = ps(k = p_int(lower = 1, upper = 10))  # Number of nearest neighbors

set.seed(130)
tuning_knn = mlr3tuning::tune(tuner = bayes_tune,
                              search_space = search_space_knn,
                              task = task_train,
                              learner = learner_knn,
                              resampling = rsmp("cv", folds = 5),
                              measure = msr("classif.auc"),
                              term_evals = 10
)

# Naive Bayes

learner_nb = lrn("classif.naive_bayes", predict_type = "prob")

# Random forest (ranger)
learner_ranger = set_threads(lrn('classif.ranger', predict_type = 'prob'), n = 14)

search_space_ranger = ps(# mtry: Number of variables randomly sampled as candidates at each split
  mtry = p_int(lower = 1, upper = floor(sqrt(ncol(train_data_downsampled)))),
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


set.seed(180)
tuning_ranger = mlr3tuning::tune(tuner = bayes_tune,
                                 search_space = search_space_ranger,
                                 task = task_train,
                                 learner = learner_ranger,
                                 resampling = rsmp("cv", folds = 5),
                                 measure = msr("classif.auc"),
                                 term_evals = 25
)

# SVM with radial kernel
learner_svm = set_threads(lrn("classif.svm", type = "C-classification", kernel = c("radial"), predict_type = "prob"), n = 14)

search_space_svm = ps(
  # Cost: Controls the trade-off between achieving a low training error and a low testing error
  cost = p_dbl(lower = 0.1, upper = 50),
  
  # Kernel: Specifies the kernel function to be used in the algorithm
  kernel = p_fct(levels = c("radial")),
  
  # Gamma: Defines how far the influence of a single training example reaches
  # low values result in smoother, less complex decision boundary
  # high value means influence is limited to points close to the decision boundary
  gamma = p_dbl(lower = 0.0001, upper = 10),
  
  # Scale: Whether to scale the features
  scale = p_fct(levels = c(FALSE))
)

set.seed(150)
tuning_svm = mlr3tuning::tune(tuner = bayes_tune,
                              search_space = search_space_svm,
                              task = task_train,
                              learner = learner_svm,
                              resampling = rsmp("cv", folds = 5),
                              measure = msr("classif.auc"),
                              term_evals = 25)

learner_glmnet$param_set$values   = tuning_glmnet$result_learner_param_vals
#learner_nb$param_set$values       = tuning_nb$result_learner_param_vals
learner_knn$param_set$values      = tuning_knn$result_learner_param_vals
learner_ranger$param_set$values   = tuning_ranger$result_learner_param_vals
learner_svm$param_set$values      = tuning_svm$result_learner_param_vals

base_learners = list(learner_glmnet, 
                     learner_nb, 
                     learner_knn,
                     learner_ranger, 
                     learner_svm)

outer_resampling = rsmp("cv", folds = 5)
inner_resampling = rsmp("cv", folds = 5)

# Super learner - XGBoost

super_learner_xgb = set_threads(lrn("classif.xgboost", predict_type = "prob"), n = 12)

search_space_xgb = ps(# eta: Learning rate (also called shrinkage). Smaller values make the model more robust to overfitting.
  eta = p_dbl(lower = 0.01, upper = 0.3), 
  # max_depth: Maximum depth of a tree. Increasing this value makes the model more complex and more likely to overfit.
  max_depth = p_int(lower = 1, upper = 10),
  # gamma: Minimum loss reduction required to make a further partition on a leaf node of the tree.
  gamma = p_dbl(lower = 0, upper = 1),
  # colsample_bytree: The fraction of features to choose for each boosting round.
  colsample_bytree = p_dbl(lower = 0.5, upper = 1),
  # min_child_weight: Minimum sum of instance weight needed in a child.
  min_child_weight = p_dbl(lower = 1, upper = 10),
  # subsample: The fraction of observations to be randomly sampled for each boosting round.
  subsample = p_dbl(lower = 0.5, upper = 1),
  # nrounds: Number of boosting rounds.
  nrounds = p_int(lower = 1500, upper = 2000),
  # alpha: L1 regularization 
  alpha = p_dbl(lower = 0, upper = 1),
  # lambda: L2 regularization 
  lambda = p_dbl(lower = 0.001, upper = 1), 
  # control the balance of positive and negative weights, which can be useful for unbalanced classes
  scale_pos_weight = p_dbl(lower = 1, upper = 10)
)

set.seed(210)
tuning_xgb = mlr3tuning::tune(tuner = bayes_tune,
                              search_space = search_space_xgb,
                              task = task_train,
                              learner = super_learner_xgb,
                              resampling = inner_resampling,
                              measure = msr("classif.auc"),
                              term_evals = 25
)

super_learner_xgb$param_set$values = tuning_xgb$result_learner_param_vals

# mlr3 pipeline stacking

set.seed(250)
graph_stack = pipeline_stacking(base_learners, super_learner_xgb, method = "cv", folds = 5, use_features = FALSE)

ensemble = as_learner(graph_stack)

set.seed(260)
rr = resample(task_train, ensemble, outer_resampling, store_models = TRUE)

best_model = rr$learners[[which.min(rr$aggregate(msr("classif.auc")))]]

best_model$train(task_train)

ensemble$train(task_train)

best_model$predict_type = "prob"
test_preds_best_model = best_model$predict(task_test)
test_perform_best_model = test_preds_best_model$score(msr("classif.auc"))

ensemble$predict_type = "prob"
test_preds_ensemble_model = ensemble$predict(task_test)
test_perform_ensemble = test_preds_ensemble_model$score(msr("classif.auc"))

row_ids = test_preds_ensemble_model$row_ids
truth = test_preds_ensemble_model$truth
response = test_preds_ensemble_model$response
probabilities = test_preds_ensemble_model$prob

predictions_df = data.frame(
  row_ids = row_ids,
  truth = truth,
  response = response,
  prob_deep_depression = probabilities[, "deep_depression"],
  prob_control = probabilities[, "shallow_depression"]
)

head(predictions_df)

deep_depression_predictions = predictions_df[predictions_df$response == "deep_depression", ]

# Performance metrics

confusion_matrix = test_preds_ensemble_model$confusion

TP = confusion_matrix[1, 1]  
TN = confusion_matrix[2, 2]  
FP = confusion_matrix[1, 2]  
FN = confusion_matrix[2, 1]


f1_score = (2 * TP) / (2 * TP + FP + FN)

cat(paste0("F1 Score: ", f1_score, "\n"))

accuracy = (TP + TN) / (TP + TN + FP + FN)

sensitivity = TP / (TP + FN)

specificity = TN / (TN + FP)

# Balanced accuracy 
balanced_acc = (sensitivity + specificity) / 2

ppv = TP / (TP + FP)

npv = TN / (TN + FN)

cat(paste0("Accuracy: ", accuracy, "\n"))
cat(paste0("Sensitivity: ", sensitivity, "\n"))
cat(paste0("Specificity: ", specificity, "\n"))
cat(paste0("Balanced Accuracy: ", balanced_acc, "\n"))
cat(paste0("Positive Predictive Value: ", ppv, "\n"))
cat(paste0("Negative Predictive Value: ", npv))

set.seed(530)

n_bootstraps = 10000

bootstrapped_accuracy = numeric(n_bootstraps)
bootstrapped_sensitivity = numeric(n_bootstraps)
bootstrapped_specificity = numeric(n_bootstraps)
bootstrapped_ppv = numeric(n_bootstraps)
bootstrapped_npv = numeric(n_bootstraps)
bootstrapped_f1 = numeric(n_bootstraps)

for (i in 1:n_bootstraps) {
  bootstrap_indices = sample(1:nrow(predictions_df), replace = TRUE)
  bootstrap_predictions = predictions_df[bootstrap_indices, ]
  
  bootstrap_confusion = table(bootstrap_predictions$truth, bootstrap_predictions$response)
  
  bootstrap_TP = bootstrap_confusion["deep_depression", "deep_depression"]
  bootstrap_TN = bootstrap_confusion["shallow_depression", "shallow_depression"]
  bootstrap_FP = bootstrap_confusion["shallow_depression", "deep_depression"]
  bootstrap_FN = bootstrap_confusion["deep_depression", "shallow_depression"]
  
  bootstrapped_accuracy[i] = (bootstrap_TP + bootstrap_TN) / sum(bootstrap_confusion)
  bootstrapped_sensitivity[i] = bootstrap_TP / (bootstrap_TP + bootstrap_FN)
  bootstrapped_specificity[i] = bootstrap_TN / (bootstrap_TN + bootstrap_FP)
  bootstrapped_ppv[i] = bootstrap_TP / (bootstrap_TP + bootstrap_FP)
  bootstrapped_npv[i] = bootstrap_TN / (bootstrap_TN + bootstrap_FN)
  bootstrapped_f1[i] = (2 * bootstrap_TP) / (2 * bootstrap_TP + bootstrap_FP + bootstrap_FN)
}

ci_accuracy = quantile(bootstrapped_accuracy, c(0.025, 0.975))
ci_sensitivity = quantile(bootstrapped_sensitivity, c(0.025, 0.975))
ci_specificity = quantile(bootstrapped_specificity, c(0.025, 0.975))
ci_ppv = quantile(bootstrapped_ppv, c(0.025, 0.975))
ci_npv = quantile(bootstrapped_npv, c(0.025, 0.975))
ci_f1 = quantile(bootstrapped_f1, c(0.025, 0.975))

print("Super Learner Performance (with bootstrapped 95% CIs):")
cat(paste0("F1 Score: ", round(f1_score, 4), " (95% CI: ", round(ci_f1[1], 4), " - ", round(ci_f1[2], 4), ")\n"))
cat(paste0("Accuracy: ", round(accuracy, 4), " (95% CI: ", round(ci_accuracy[1], 4), " - ", round(ci_accuracy[2], 4), ")\n"))
cat(paste0("Sensitivity: ", round(sensitivity, 4), " (95% CI: ", round(ci_sensitivity[1], 4), " - ", round(ci_sensitivity[2], 4), ")\n"))
cat(paste0("Specificity: ", round(specificity, 4), " (95% CI: ", round(ci_specificity[1], 4), " - ", round(ci_specificity[2], 4), ")\n"))
cat(paste0("PPV: ", round(ppv, 4), " (95% CI: ", round(ci_ppv[1], 4), " - ", round(ci_ppv[2], 4), ")\n"))
cat(paste0("NPV: ", round(npv, 4), " (95% CI: ", round(ci_npv[1], 4), " - ", round(ci_npv[2], 4), ")\n"))

set.seed(550)

n_bootstraps = 10000

bootstrapped_auc = numeric(n_bootstraps)

pred_probs = test_preds_ensemble_model$prob[, "deep_depression"]
true_labels = test_preds_ensemble_model$truth

for (i in 1:n_bootstraps) {
  bootstrap_indices = sample(1:length(true_labels), replace = TRUE)
  bootstrap_pred_probs = pred_probs[bootstrap_indices]
  bootstrap_true_labels = true_labels[bootstrap_indices]
  
  bootstrap_roc = roc(bootstrap_true_labels, bootstrap_pred_probs)
  bootstrapped_auc[i] = auc(bootstrap_roc)
}

ci_auc = quantile(bootstrapped_auc, c(0.025, 0.975))

original_roc = roc(true_labels, pred_probs)
original_auc = auc(original_roc)

cat(paste0("AUC: ", round(original_auc, 3), 
           " (95% CI: ", round(ci_auc[1], 3), " - ", round(ci_auc[2], 3), ")\n"))

roc_df = data.frame(
  Sensitivity = original_roc$sensitivities,
  Specificity = 1 - original_roc$specificities
)

# PR-AUC calculation
pred_probs_deep_depression = test_preds_ensemble_model$prob[, "deep_depression"]
true_labels = test_preds_ensemble_model$truth

pr_obj = pr.curve(scores.class0 = pred_probs_deep_depression, 
                  weights.class0 = ifelse(true_labels == "deep_depression", 1, 0),
                  curve = TRUE)
auprc = pr_obj$auc.integral

cat(paste0("AUPRC: ", round(auprc, 3), "\n"))

# Calculate 95% confidence interval for AUPRC using bootstrapping
set.seed(520)  # for reproducibility
n_bootstraps = 10000
bootstrapped_auprc = numeric(n_bootstraps)

for (i in 1:n_bootstraps) {
  bootstrap_indices = sample(1:length(true_labels), replace = TRUE)
  bootstrap_pred_probs = pred_probs_deep_depression[bootstrap_indices]
  bootstrap_true_labels = true_labels[bootstrap_indices]
  
  bootstrap_pr = pr.curve(scores.class0 = bootstrap_pred_probs, 
                          weights.class0 = ifelse(bootstrap_true_labels == "deep_depression", 1, 0),
                          curve = TRUE)
  bootstrapped_auprc[i] = bootstrap_pr$auc.integral
}

ci_auprc = quantile(bootstrapped_auprc, c(0.025, 0.975))

cat(paste0("AUPRC: ", round(auprc, 3), 
           " (95% CI: ", round(ci_auprc[1], 3), " - ", round(ci_auprc[2], 3), ")\n"))

# Learning curves
# Function to generate learning curve data
generate_learning_curve = function(learner, task, measure, n_points = 7) {
  train_sizes = seq(0.1, 1, length.out = n_points)
  results = list()
  
  for (size in train_sizes) {
    
    # Create a subsampled task
    subtask = task$clone()
    subtask$filter(sample(task$nrow, size = floor(task$nrow * size)))
    
    # Create a resampling instance
    resampling = rsmp("cv", folds = 5)
    
    # Perform resampling
    rr = resample(subtask, learner, resampling)
    
    # Store results
    results[[as.character(size)]] = list(
      size = size,
      performance = rr$aggregate(measure)
    )
  }
  
  return(rbindlist(results))
}

# Generate learning curve data
set.seed(502)  # for reproducibility
learning_curve_data = generate_learning_curve(
  learner = ensemble,
  task = task_train,
  measure = msr("classif.auc"),
  n_points = 7
)

# Plot the learning curve
ggplot(learning_curve_data, aes(x = size, y = performance)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Learning Curve for Meta Learner",
       x = "Training Set Size (Proportion)",
       y = "Performance (AUC)")


# Identify best base learner

# Function to evaluate a single learner
evaluate_learner = function(learner, task, resampling) {
  set.seed(300)  # for reproducibility
  rr = resample(task, learner, resampling)
  return(list(
    learner = learner$id,
    auc = rr$aggregate(msr("classif.auc")),
    acc = rr$aggregate(msr("classif.acc")),
    bal_acc = rr$aggregate(msr("classif.bacc"))
  ))
}

base_learner_results = lapply(base_learners, function(lrn) {
  evaluate_learner(lrn, task_train, rsmp("cv", folds = 5))
})

base_learner_df = do.call(rbind, lapply(base_learner_results, function(x) {
  data.frame(
    learner = x$learner,
    algorithm = paste(x$algorithm, collapse = ", "),
    auc = x$auc,
    acc = x$acc,
    bal_acc = x$bal_acc
  )
}))

best_base_learner = base_learner_df[which.max(base_learner_df$auc), ]

print("Base Learner Performance:")
print(base_learner_df[, c("learner", "algorithm", "auc", "acc", "bal_acc")])

cat("\nBest Base Learner:", best_base_learner$learner, "\n")
cat("Algorithm:", best_base_learner$algorithm, "\n")
cat("AUC:", best_base_learner$auc, "\n")
cat("Accuracy:", best_base_learner$acc, "\n")
cat("Balanced Accuracy:", best_base_learner$bal_acc, "\n")

# Extended performance metrics for best base learner

best_base_learner_model = base_learners[[which(sapply(base_learners, function(x) x$id) == best_base_learner$learner)]]
best_base_learner_model$train(task_train)

best_base_preds = best_base_learner_model$predict(task_test)

confusion_matrix_base = best_base_preds$confusion
print("SVM Base Learner:")
print(confusion_matrix_base)

TP_base = confusion_matrix_base[1, 1]  
TN_base = confusion_matrix_base[2, 2]  
FP_base = confusion_matrix_base[1, 2]  
FN_base = confusion_matrix_base[2, 1]

accuracy_base = (TP_base + TN_base) / (TP_base + TN_base + FP_base + FN_base)
sensitivity_base = TP_base / (TP_base + FN_base)
specificity_base = TN_base / (TN_base + FP_base)
balanced_acc_base = (sensitivity_base + specificity_base) / 2
ppv_base = TP_base / (TP_base + FP_base)
npv_base = TN_base / (TN_base + FN_base)
f1_score_base = (2 * TP_base) / (2 * TP_base + FP_base + FN_base)

cat("\nBest Base Learner Performance Metrics:\n")
cat("\nBest Base Learner:", best_base_learner$learner, "\n")
cat(paste0("AUC: ", best_base_preds$score(msr("classif.auc")), "\n"))
cat(paste0("Accuracy: ", accuracy_base, "\n"))
cat(paste0("Sensitivity: ", sensitivity_base, "\n"))
cat(paste0("Specificity: ", specificity_base, "\n"))
cat(paste0("Balanced Accuracy: ", balanced_acc_base, "\n"))
cat(paste0("Positive Predictive Value: ", ppv_base, "\n"))
cat(paste0("Negative Predictive Value: ", npv_base, "\n"))
cat(paste0("F1 Score: ", f1_score_base, "\n"))

# Compare best base learner to ensemble

cat("\nPerformance Comparison:\n")
cat("Metric\t\t\tBest Base Learner\tEnsemble Model\n")
cat(paste0("F1 Score\t\t", round(f1_score_base, 4), "\t\t\t", round(f1_score, 4), "\n"))
cat(paste0("AUC\t\t\t", round(best_base_preds$score(msr("classif.auc")), 4), "\t\t\t", 
           round(test_perform_ensemble, 4), "\n"))
cat(paste0("Balanced Accuracy\t", round(balanced_acc_base, 4), "\t\t\t", round(balanced_acc, 4), "\n"))
cat(paste0("Sensitivity\t\t", round(sensitivity_base, 4), "\t\t\t", round(sensitivity, 4), "\n"))
cat(paste0("Specificity\t\t", round(specificity_base, 4), "\t\t\t", round(specificity, 4), "\n"))
cat(paste0("PPV\t\t\t", round(ppv_base, 4), "\t\t\t", round(ppv, 4), "\n"))
cat(paste0("NPV\t\t\t", round(npv_base, 4), "\t\t\t", round(npv, 4), "\n"))





