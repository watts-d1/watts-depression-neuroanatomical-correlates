# Neuroimaging Feature Extraction Pipeline - Heterogeneity of Depression Part 2
#
# Note: Information pertaining to neuroimaging preprocessing can be found here:
#nhttps://doi.org/10.1016/j.neuroimage.2017.10.034
# Author: Devon Watts

# --- User-defined paths (modify for your environment) ---
S1_RESULTS_DIR = "/path/to/S1/logistic_regression/results"  # Directory containing S1 per-phenotype .RData output files
S1_ICD_OBJECTS_PATH = "/path/to/all_completed_ICD_objects_z_score.RData"  # S1 output: ICD z-score objects
LOOKUP_TABLE_PATH = "/path/to/lookup_table_neuroimaging_ukbb.xlsx"  # Neuroimaging variable lookup table

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
             'car', 'future', 'future.apply', 'parallel', 'furrr', 'ggplot2', 
             'pbapply', 'progressr', 'gridExtra', 'pbapply', 'vctrs', 'readxl', 
             'tidyr', 'pheatmap', 'ggrepel', 'boot')

package_install(packages)
library(here)

logistic_regression_names = c("LR_GPNoDep_PL_95_w_MD.RDS", 
                              "LR_SelfRepDep_PL_95_w_MD.RDS", 
                              "LR_DepAll_PL_95_w_MD.RDS",
                              "LR_ICD10Dep_PL_95_w_MD.RDS",
                              "LR_ICD10Dep_exclpsych_PL_95_w_MD.RDS",
                              "LR_LifetimeMDD_PL_95_w_MD.RDS",
                              "LR_MDDRecur_PL_95_w_MD.RDS",
                              "LR_Psypsy_PL_95_w_MD.RDS",
                              "LR_GPpsy_PL_95_w_MD.RDS"
                              
)

for (file in logistic_regression_names) {
  file_path = file.path(S1_RESULTS_DIR, file)
  data = readRDS(file_path)
  obj_name = sub("\\.RDS$", "", file)
  assign(obj_name, data, envir = .GlobalEnv)
}

# Extracting feature names

lookup_table_path = LOOKUP_TABLE_PATH

read_lookup_table = function(filepath) {
  read_excel(filepath)
}

lookup_table_df = read_lookup_table(lookup_table_path)

# Harmonize prefix: "baseline_" -> "bl_" to match regression output naming
lookup_table_df$corrected_id = gsub("baseline_", "bl_", lookup_table_df$corrected_id)

#' Join results with lookup table to get descriptive feature names
#' @return List with 'matched' and 'no_match' data frames
merge_with_lookup_table = function(results_df, lookup_df) {
  results_df$feature = as.character(results_df$feature)
  lookup_df$corrected_id = as.character(lookup_df$corrected_id)

  lookup_df = lookup_df[, c("corrected_id", "field", "processed_name_abbreviated")]

  merged_df = merge(results_df, lookup_df, by.x = 'feature', by.y = 'corrected_id', all.x = TRUE)

  no_match_rows = is.na(merged_df$processed_name_abbreviated)
  merged_df$processed_name_abbreviated[no_match_rows] = "no_match"

  no_match_df = merged_df[no_match_rows, ]
  matched_df = merged_df[!no_match_rows, ]

  return(list(matched = matched_df, no_match = no_match_df))
}

# Function to Extract Significant Features

#' Extract Significant Features from Logistic Regression Results
#'
#' Processes logistic regression results to identify and describe statistically significant features.
#' This function filters features based on FDR-adjusted significance and combines this information 
#' with odds ratios, confidence intervals, and descriptive names from a lookup table.
#'
#' @param combined_df Data frame that combines logistic regression results with odds ratios 
#' and confidence intervals for each feature.
#' @param lookup_df Data frame from the lookup table containing feature descriptions.
#'
#' @return A list containing:
#'         - `significant_features`: Data frame with significant features, their descriptive names, 
#'           odds ratios, and confidence intervals, sorted by FDR p-value and odds ratio.
#'         - `number_significant`: Total number of significant features based on the FDR criterion.
#'         - `no_match_features`: Data frame with significant features that did not find a match in 
#'           the lookup table, requiring further investigation.
#'
#' @examples
#' # Assuming combined_df and lookup_df are predefined data frames
#' significant_results = extract_significant_features(combined_df, lookup_df)
#' @export
extract_significant_features = function(combined_df, lookup_df) {
  # Check if there are any significant results
  if (!any(combined_df$significant_fdr)) {
    return(list(
      significant_features = data.frame(),
      number_significant = 0,
      no_match_features = data.frame()
    ))
  }
  
  # Filter out only significant results
  significant_df = combined_df[combined_df$significant_fdr, ]
  
  # Merge with lookup table to get descriptive names
  merged_data = merge_with_lookup_table(significant_df, lookup_df)
  significant_df_with_names = merged_data$matched
  no_match_df = merged_data$no_match
  
  # Order by FDR p-value and odds ratio
  ordered_significant_df = significant_df_with_names[order(significant_df_with_names$fdr_p_value, -significant_df_with_names$odds_ratio), ]
  
  return(list(
    significant_features = ordered_significant_df,
    number_significant = nrow(ordered_significant_df),
    no_match_features = no_match_df
  ))
}

# Function to Process Results for All Datasets

##' Process Results for Multiple Datasets
#'
#' Applies the `extract_significant_features` function to a list of data frames containing logistic 
#' regression results for multiple datasets. This function is designed to process and extract 
#' significant features for each dataset in the list.
#'
#' @param results_list A list of data frames, each representing logistic regression results for 
#' different datasets. These data frames should include odds ratios, confidence intervals, and 
#' significance indicators.
#' @param lookup_df Data frame from the lookup table containing feature descriptions.
#'
#' @return A list where each element corresponds to the processed significant results for each dataset 
#' in `results_list`. Each element in the list is itself a list containing `significant_features`, 
#' `number_significant`, and `no_match_features`.
#'
#' @examples
#' # Assuming results_list and lookup_df are predefined
#' processed_results = process_results(results_list, lookup_df)
#' @export
process_results = function(results_list, lookup_df) {
  lapply(results_list, function(df) {
    extract_significant_features(df, lookup_df)
  })
}

# Create results_list - all cases and controls 
# 
# Including all features, volume, surface_area, thickness, and FA
results_list_all = list(
  GPNoDep            = LR_GPNoDep_PL_95_w_MD$All_Features$results_df,
  Psypsy             = LR_Psypsy_PL_95_w_MD$All_Features$results_df,
  SelfRepDep         = LR_SelfRepDep_PL_95_w_MD$All_Features$results_df,
  GPpsy              = LR_GPpsy_PL_95_w_MD$All_Features$results_df,
  DepAll             = LR_DepAll_PL_95_w_MD$All_Features$results_df,
  ICD10Dep           = LR_ICD10Dep_PL_95_w_MD$All_Features$results_df,
  ICD10Dep_exclpsych = LR_ICD10Dep_exclpsych_PL_95_w_MD$All_Features$results_df,
  LifetimeMDD        = LR_LifetimeMDD_PL_95_w_MD$All_Features$results_df,
  MDDRecur           = LR_MDDRecur_PL_95_w_MD$All_Features$results_df
)

results_list_global_volume = list(
  GPNoDep            = LR_GPNoDep_PL_95_w_MD$global_volume$results_df,
  Psypsy             = LR_Psypsy_PL_95_w_MD$global_volume$results_df,
  SelfRepDep         = LR_SelfRepDep_PL_95_w_MD$global_volume$results_df,
  GPpsy              = LR_GPpsy_PL_95_w_MD$global_volume$results_df,
  DepAll             = LR_DepAll_PL_95_w_MD$global_volume$results_df,
  ICD10Dep           = LR_ICD10Dep_PL_95_w_MD$global_volume$results_df,
  ICD10Dep_exclpsych = LR_ICD10Dep_exclpsych_PL_95_w_MD$global_volume$results_df,
  LifetimeMDD        = LR_LifetimeMDD_PL_95_w_MD$global_volume$results_df,
  MDDRecur           = LR_MDDRecur_PL_95_w_MD$global_volume$results_df
)


# Volume metrics 
results_list_subcortical_volume = list(
  GPNoDep            = LR_GPNoDep_PL_95_w_MD$subcortical_volume$results_df,
  Psypsy             = LR_Psypsy_PL_95_w_MD$subcortical_volume$results_df,
  SelfRepDep         = LR_SelfRepDep_PL_95_w_MD$subcortical_volume$results_df,
  GPpsy              = LR_GPpsy_PL_95_w_MD$subcortical_volume$results_df,
  DepAll             = LR_DepAll_PL_95_w_MD$subcortical_volume$results_df,
  ICD10Dep           = LR_ICD10Dep_PL_95_w_MD$subcortical_volume$results_df,
  ICD10Dep_exclpsych = LR_ICD10Dep_exclpsych_PL_95_w_MD$subcortical_volume$results_df,
  LifetimeMDD        = LR_LifetimeMDD_PL_95_w_MD$subcortical_volume$results_df,
  MDDRecur           = LR_MDDRecur_PL_95_w_MD$subcortical_volume$results_df
)

# Surface Area
results_list_surface_area = list(
  GPNoDep            = LR_GPNoDep_PL_95_w_MD$surface_area$results_df,
  Psypsy             = LR_Psypsy_PL_95_w_MD$surface_area$results_df,
  SelfRepDep         = LR_SelfRepDep_PL_95_w_MD$surface_area$results_df,
  GPpsy              = LR_GPpsy_PL_95_w_MD$surface_area$results_df,
  DepAll             = LR_DepAll_PL_95_w_MD$surface_area$results_df,
  ICD10Dep           = LR_ICD10Dep_PL_95_w_MD$surface_area$results_df,
  ICD10Dep_exclpsych = LR_ICD10Dep_exclpsych_PL_95_w_MD$surface_area$results_df,
  LifetimeMDD        = LR_LifetimeMDD_PL_95_w_MD$surface_area$results_df,
  MDDRecur           = LR_MDDRecur_PL_95_w_MD$surface_area$results_df
)

# Thickness 
results_list_thickness = list(
  GPNoDep            = LR_GPNoDep_PL_95_w_MD$thickness$results_df,
  Psypsy             = LR_Psypsy_PL_95_w_MD$thickness$results_df,
  SelfRepDep         = LR_SelfRepDep_PL_95_w_MD$thickness$results_df,
  GPpsy              = LR_GPpsy_PL_95_w_MD$thickness$results_df,
  DepAll             = LR_DepAll_PL_95_w_MD$thickness$results_df,
  ICD10Dep           = LR_ICD10Dep_PL_95_w_MD$thickness$results_df,
  ICD10Dep_exclpsych = LR_ICD10Dep_exclpsych_PL_95_w_MD$thickness$results_df,
  LifetimeMDD        = LR_LifetimeMDD_PL_95_w_MD$thickness$results_df,
  MDDRecur           = LR_MDDRecur_PL_95_w_MD$thickness$results_df
)


# Fractional anisotropy
results_list_FA = list(
  GPNoDep            = LR_GPNoDep_PL_95_w_MD$FA$results_df,
  Psypsy             = LR_Psypsy_PL_95_w_MD$FA$results_df,
  SelfRepDep         = LR_SelfRepDep_PL_95_w_MD$FA$results_df,
  GPpsy              = LR_GPpsy_PL_95_w_MD$FA$results_df,
  DepAll             = LR_DepAll_PL_95_w_MD$FA$results_df,
  ICD10Dep           = LR_ICD10Dep_PL_95_w_MD$FA$results_df,
  ICD10Dep_exclpsych = LR_ICD10Dep_exclpsych_PL_95_w_MD$FA$results_df,
  LifetimeMDD        = LR_LifetimeMDD_PL_95_w_MD$FA$results_df,
  MDDRecur           = LR_MDDRecur_PL_95_w_MD$FA$results_df
)

# Mean Diffusivity 
results_list_MD = list(
  GPNoDep            = LR_GPNoDep_PL_95_w_MD$MD$results_df,
  Psypsy             = LR_Psypsy_PL_95_w_MD$MD$results_df,
  SelfRepDep         = LR_SelfRepDep_PL_95_w_MD$MD$results_df,
  GPpsy              = LR_GPpsy_PL_95_w_MD$MD$results_df,
  DepAll             = LR_DepAll_PL_95_w_MD$MD$results_df,
  ICD10Dep           = LR_ICD10Dep_PL_95_w_MD$MD$results_df,
  ICD10Dep_exclpsych = LR_ICD10Dep_exclpsych_PL_95_w_MD$MD$results_df,
  LifetimeMDD        = LR_LifetimeMDD_PL_95_w_MD$MD$results_df,
  MDDRecur           = LR_MDDRecur_PL_95_w_MD$MD$results_df
)

# Saving surface area and thickness results lists for downstream analyses (e.g. ENIGMA)
save(results_list_surface_area, results_list_thickness, file = here("output", "results_lists_thickness_SA.RData"))
save(results_list_FA, results_list_MD, file = here("output", "results_lists_FA_MD.RData"))

# Accessing significant results across datasets with all cases and controls
# Note: To access the significant results for a specific dataset, use the $ operator
# e.g., significant_results$GPNoDep$significant_features

# Process the results
significant_results_all                = process_results(results_list_all, lookup_table_df)
significant_results_global_volume      = process_results(results_list_global_volume, lookup_table_df)
significant_results_subcortical_volume = process_results(results_list_subcortical_volume, lookup_table_df)
significant_results_thickness          = process_results(results_list_thickness, lookup_table_df)
significant_results_surface_area       = process_results(results_list_surface_area, lookup_table_df)
significant_results_fa                 = process_results(results_list_FA, lookup_table_df)
significant_results_md                 = process_results(results_list_MD, lookup_table_df)

post_process_results = function(significant_results_list) {
  lapply(significant_results_list, function(result) {
    if (!is.null(result$significant_features) && nrow(result$significant_features) > 0) {
      result$significant_features %>%
        select(feature, odds_ratio, lower_CI, upper_CI, fdr_p_value, processed_name_abbreviated) %>%
        mutate(
          feature = gsub("^bl_", "", feature),  # Remove "bl_" prefix from feature
          odds_ratio = round(odds_ratio, 3),
          lower_CI = round(lower_CI, 3),
          upper_CI = round(upper_CI, 3),
          fdr_p_value = fdr_p_value,
          processed_name_abbreviated = gsub("_", " ", processed_name_abbreviated)  # Replace underscores with spaces
        )
    } else {
      data.frame() # Return an empty data frame if there are no significant features
    }
  })
}

# Usage example:
significant_results_fa_processed = post_process_results(significant_results_fa)
significant_results_md_processed = post_process_results(significant_results_md)
significant_results_surface_area_processed = post_process_results(significant_results_surface_area)
significant_results_thickness_processed = post_process_results(significant_results_thickness)


# Saving significant result objects 
save(significant_results_all,
     significant_results_global_volume,
     significant_results_subcortical_volume,
     significant_results_thickness,
     significant_results_surface_area,
     significant_results_fa,
     significant_results_md,
     file = here("output", "significant_results.RData")
)

create_forest_plot = function(data, title) {
  data %>%
    mutate(
      log_odds_ratio = log(odds_ratio),
      log_lower_CI = log(lower_CI),
      log_upper_CI = log(upper_CI)
    ) %>%
    ggplot(aes(x = processed_name_abbreviated, y = log_odds_ratio, ymin = log_lower_CI, ymax = log_upper_CI)) +
    geom_pointrange(color = "blue", size = 0.8) +
    geom_errorbar(aes(ymin = log_lower_CI, ymax = log_lower_CI), color = "black", width = 0.2) +
    geom_errorbar(aes(ymin = log_upper_CI, ymax = log_upper_CI), color = "black", width = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    coord_flip() +
    labs(title = title, x = "Cortical Region", y = "Log Odds Ratio")
}

# Process significant_results_thickness
for (def in names(significant_results_thickness)) {
  sig_features = significant_results_thickness[[def]]$significant_features
  
  if (nrow(sig_features) > 0) {
    plot_title = paste("Cortical Thickness:", def, "- Significant Associations")
    forest_plot = create_forest_plot(sig_features, plot_title)
    print(forest_plot)
  }
}

# Process significant_results_surface_area
for (def in names(significant_results_surface_area)) {
  sig_features = significant_results_surface_area[[def]]$significant_features
  
  if (nrow(sig_features) > 0) {
    plot_title = paste("Surface Area:", def, "- Significant Associations")
    forest_plot = create_forest_plot(sig_features, plot_title)
    print(forest_plot)
  }
}

# Process significant_results_fa
for (def in names(significant_results_fa)) {
  sig_features = significant_results_fa[[def]]$significant_features
  
  if (nrow(sig_features) > 0) {
    plot_title = paste("Fractional Anisotropy (FA):", def, "- Significant Associations")
    forest_plot = create_forest_plot(sig_features, plot_title)
    print(forest_plot)
  }
}

# Process significant_results_md
for (def in names(significant_results_md)) {
  sig_features = significant_results_md[[def]]$significant_features
  
  if (nrow(sig_features) > 0) {
    plot_title = paste("Mean Diffusivity (MD):", def, "- Significant Associations")
    forest_plot = create_forest_plot(sig_features, plot_title)
    print(forest_plot)
  }
}

# Loading Depression Dataframes for Correlations across Depression Labels
load(S1_ICD_OBJECTS_PATH)

case_control_counts = list(
  DepAll             = z_score_completed_DepAll$DepAll,
  GPNoDep            = z_score_completed_GPNoDep$GPNoDep,
  Psypsy             = z_score_completed_Psypsy$Psypsy,
  SelfRepDep         = z_score_completed_SelfRepDep$SelfRepDep, 
  GPpsy              = z_score_completed_GPpsy$GPpsy, 
  ICD10Dep           = z_score_completed_ICD10Dep$ICD10Dep, 
  ICD10Dep_exclpsych = z_score_completed_ICD10Dep_exclpsych$ICD10Dep.exclpsych, 
  LifetimeMDD        = z_score_completed_LifetimeMDD$LifetimeMDD, 
  MDDRecurr          = z_score_completed_MDDRecurr$MDDRecur
)

# Function to Calculate Correlations Across Depression Labels

#' Calculate Correlations of Effect Sizes Across Depression Labels
#'
#' @param results_list List of odds ratio dataframes for different labels.
#'
#' @return List of correlation results between each pair of labels.
calculate_correlations_across_labels = function(results_list) {
  correlation_results = data.frame(
    Group1 = character(),
    Group2 = character(),
    SpearmanCorrelation = numeric(),
    Observations = integer(),
    PValue = numeric(),
    FDR_PValue = numeric(),
    stringsAsFactors = FALSE
  )
  
  p_values = numeric()
  scatter_plot_data = list()
  
  ors_list = results_list
  
  for (i in 1:(length(ors_list) - 1)) {
    for (j in (i + 1):length(ors_list)) {
      df1 = ors_list[[i]]
      df2 = ors_list[[j]]
      label1 = names(ors_list)[i]
      label2 = names(ors_list)[j]
      
      common_features = intersect(df1$feature, df2$feature)
      if (length(common_features) > 0) {
        odds_ratios1 = df1[df1$feature %in% common_features, "odds_ratio"]
        odds_ratios2 = df2[df2$feature %in% common_features, "odds_ratio"]
        
        spearman_correlation = cor(odds_ratios1, odds_ratios2, method = "spearman")
        p_value = cor.test(odds_ratios1, odds_ratios2, method = "spearman")$p.value
        
        correlation_results = rbind(correlation_results, data.frame(
          Group1 = label1, 
          Group2 = label2, 
          SpearmanCorrelation = spearman_correlation,
          Observations = length(common_features),
          PValue = p_value
        ))
        
        p_values = c(p_values, p_value)
        
        scatter_plot_data[[paste(label1, label2, sep = "_vs_")]] = list(
          Data = data.frame(Feature = common_features, OddsRatio1 = odds_ratios1, OddsRatio2 = odds_ratios2),
          Observations = length(common_features),
          PValue = p_value
        )
      }
    }
  }
  
  # Apply FDR correction to the gathered p-values
  fdr_corrected_p_values = p.adjust(p_values, method = "fdr")
  
  # Assign corrected p-values back to the correlation_results
  correlation_results$FDR_PValue = fdr_corrected_p_values
  
  return(list(Correlations = correlation_results, ScatterData = scatter_plot_data))
}

# Calculating correlations for each depression definition
# Separately considering all features, global volume, sub-cortical volume, surface area, thickness, FA, and MD
# Observations in this case refers to the number of available features 
correlation_results_all           = calculate_correlations_across_labels(results_list_all)
correlation_results_global_volume = calculate_correlations_across_labels(results_list_global_volume)
correlation_results_sub_volume    = calculate_correlations_across_labels(results_list_subcortical_volume)
correlation_results_surface_area  = calculate_correlations_across_labels(results_list_surface_area)
correlation_results_thickness     = calculate_correlations_across_labels(results_list_thickness)
correlation_results_fa            = calculate_correlations_across_labels(results_list_FA)
correlation_results_md            = calculate_correlations_across_labels(results_list_MD)

#' Plot a scatter plot with correlations and p-values
#'
#' This function creates a scatter plot of odds ratios between two groups. It also adds a linear model fit line 
#' and annotates the plot with the Spearman correlation coefficient, the Spearman p-value, and the FDR 
#' corrected p-value.
#'
#' @param data A data frame containing the variables `OddsRatio1` and `OddsRatio2` and optionally `Feature`.
#' @param group1_label A character string for the label of group 1.
#' @param group2_label A character string for the label of group 2.
#' @param spearman_corr The Spearman correlation coefficient between the two groups.
#' @param p_value The p-value corresponding to the Spearman correlation coefficient.
#' @param fdr_p_value The FDR corrected p-value.
#'
#' @return A ggplot object representing the scatter plot.
#' @export
#'
#' @examples
#' # Assuming `correlation_results_all` and `correlation_info` are predefined:
#' scatter_plot = plot_correlation_scatter(
#'   data = correlation_results_all$ScatterData[["GPNoDep_vs_Psypsy"]]$Data,
#'   group1_label = "GPNoDep",
#'   group2_label = "Psypsy",
#'   spearman_corr = correlation_info$SpearmanCorrelation,
#'   p_value = correlation_info$PValue,
#'   fdr_p_value = correlation_info$FDR_PValue
#' )
#' print(scatter_plot)
plot_correlation_scatter = function(data, group1_label, group2_label, spearman_corr, p_value, fdr_p_value) {
  # Construct title and subtitle text
  title_text = paste("Scatter Plot for", group1_label, "vs", group2_label)
  subtitle_text = paste("Spearman Correlation:", round(spearman_corr, 4),
                        "Spearman P-Value:", round(p_value, 4, scientific = TRUE))
  
  
  # Modify fdr_p_value based on condition
  fdr_p_value_text = ifelse(fdr_p_value < 0.001, "< 0.001", format(fdr_p_value, digits = 4))
  
  # Base plot with points and line of best fit
  p = ggplot(data, aes(x = OddsRatio1, y = OddsRatio2)) +
    geom_point() +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = paste("Odds Ratio -", group1_label),
      y = paste("Odds Ratio -", group2_label)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "none"
    )
  
  # Annotate the FDR corrected p-value on the plot within the plotting area
  p = p + annotate("text", x = max(data$OddsRatio1), y = min(data$OddsRatio2), 
                   label = paste("FDR P-Value:", round(fdr_p_value, 4)), 
                   hjust = 1, vjust = 0, size = 4, color = "red", fontface = "bold")
  
  return(p)
}

# Retrieve the Spearman correlation and the FDR p-value for GPNoDep vs Psypsy
correlation_info = correlation_results_all$Correlations[correlation_results_all$Correlations$Group1 == "ICD10Dep" & correlation_results_all$Correlations$Group2 == "LifetimeMDD", ]

plot_correlation_scatter = function(data, group1_label, group2_label, spearman_corr, p_value, fdr_p_value) {
  # Construct title and subtitle text
  title_text = paste("Scatter Plot for", group1_label, "vs", group2_label)
  subtitle_text = paste("Spearman Correlation:", round(spearman_corr, 6),
                        "Spearman P-Value:", format(p_value, scientific = TRUE))
  
  # Custom function to format the FDR P-value
  format_p_value = function(p_value) {
    if (p_value < 0.001) {
      return("< 0.001")
    } else {
      # Format to avoid scientific notation for small non-zero values
      return(formatC(p_value, format = "f", digits = 6))
    }
  }
  
  # Apply the custom formatting function to the FDR P-value
  fdr_p_value_text = format_p_value(fdr_p_value)
  
  # Base plot with points and line of best fit
  p = ggplot(data, aes(x = OddsRatio1, y = OddsRatio2)) +
    geom_point() +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = paste("Odds Ratio -", group1_label),
      y = paste("Odds Ratio -", group2_label)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "none"
    )
  
  # Annotate the FDR corrected p-value on the plot within the plotting area
  p = p + annotate("text", x = max(data$OddsRatio1), y = min(data$OddsRatio2), 
                   label = paste("FDR P-Value:", fdr_p_value_text), 
                   hjust = 1, vjust = 0, size = 4, color = "red", fontface = "bold")
  
  return(p)
}
# Generate the plot for GPNoDep vs Psypsy
scatter_plot = plot_correlation_scatter(
  data = correlation_results_all$ScatterData[["ICD10Dep_vs_LifetimeMDD"]]$Data,
  group1_label = "ICD10Dep",
  group2_label = "LifetimeMDD",
  spearman_corr = correlation_info$SpearmanCorrelation,
  p_value = correlation_info$PValue,
  fdr_p_value = correlation_info$FDR_PValue
)

# To view or save the plot
print(scatter_plot)

#' Reorder dataframes according to severity 
reorder_dataframes = function(sensitivity_list, order_sequence) {
  # Initialize an empty list to store the reordered dataframes
  reordered_list = list()
  
  # Iterate through the order_sequence and extract the corresponding dataframe
  for (group in order_sequence) {
    if (group %in% names(sensitivity_list)) {
      reordered_list[[group]] = sensitivity_list[[group]]
    }
  }
  
  return(reordered_list)
}

# Define the sequence
order_sequence = c("GPNoDep", "GPpsy", "Psypsy", "DepAll", "SelfRepDep", 
                   "MDD_noimpairment", "MDD2_w_impair", "ICD10Dep", 
                   "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur")

# Reorder sensitivity results list according to sequential shallow vs deep definitions
results_list_all                = reorder_dataframes(results_list_all, order_sequence)
results_list_global_volume      = reorder_dataframes(results_list_global_volume, order_sequence)
results_list_subcortical_volume = reorder_dataframes(results_list_subcortical_volume, order_sequence)
results_list_surface_area       = reorder_dataframes(results_list_surface_area, order_sequence)
results_list_thickness          = reorder_dataframes(results_list_thickness, order_sequence)
results_list_fa                 = reorder_dataframes(results_list_FA, order_sequence)
results_list_md                 = reorder_dataframes(results_list_MD, order_sequence)


#' Perform Dependency Analysis of Effect Sizes Across Depression Labels

#' Calculate Confidence Interval Differences Between Groups
calculate_ci_difference_pairwise = function(results_list, alpha = 0.05) {
  all_pairwise_results = data.frame()
  group_names = names(results_list)
  comparisons_performed = 0
  comparisons_skipped = 0
  
  for (i in 1:(length(group_names) - 1)) {
    for (j in (i + 1):length(group_names)) {
      group1 = group_names[i]
      group2 = group_names[j]
      
      df1 = results_list[[group1]]$significant_features
      df2 = results_list[[group2]]$significant_features
      
      if (nrow(df1) > 0 && nrow(df2) > 0) {
        merged_df = merge(df1, df2, by = "feature", suffixes = c(".x", ".y"))
        
        if (nrow(merged_df) > 0) {
          ci_diff_results = tryCatch({
            calculate_ci_difference(merged_df, group1, group2)
          }, error = function(e) {
            message(paste("Error in calculating CI difference for", group1, "and", group2, ":", e$message))
            return(NULL)
          })
          
          if (!is.null(ci_diff_results) && nrow(ci_diff_results) > 0) {
            all_pairwise_results = rbind(all_pairwise_results, ci_diff_results)
            comparisons_performed = comparisons_performed + 1
          } else {
            comparisons_skipped = comparisons_skipped + 1
            message(paste("No valid results for comparison between", group1, "and", group2))
          }
        } else {
          comparisons_skipped = comparisons_skipped + 1
          message(paste("No common features found between", group1, "and", group2))
        }
      } else {
        comparisons_skipped = comparisons_skipped + 1
        message(paste("Skipping comparison between", group1, "and", group2, "due to insufficient significant features"))
      }
    }
  }
  
  if (nrow(all_pairwise_results) > 0) {
    all_pairwise_results$adjusted_p_value = p.adjust(all_pairwise_results$p_value, method = "fdr")
    all_pairwise_results$is_significant_adjusted = all_pairwise_results$adjusted_p_value < alpha 
    all_pairwise_results$is_significant_unadjusted = all_pairwise_results$p_value < alpha 
    message(paste("Comparisons performed:", comparisons_performed))
    message(paste("Comparisons skipped:", comparisons_skipped))
    return(all_pairwise_results)
  } else {
    message("No pairwise comparisons could be performed.")
    message(paste("Comparisons skipped:", comparisons_skipped))
    return(NULL)
  }
}

# Calculate the CI differences pairwise across all groups
pairwise_results_all = calculate_ci_difference_pairwise(significant_results_all, alpha = 0.01)
pairwise_results_fa  = calculate_ci_difference_pairwise(significant_results_fa, alpha = 0.01)
pairwise_results_md  = calculate_ci_difference_pairwise(significant_results_md, alpha = 0.01)
pairwise_results_sa  = calculate_ci_difference_pairwise(significant_results_surface_area, alpha = 0.05, n_perm = 10000)
pairwise_results_ct  = calculate_ci_difference_pairwise(significant_results_thickness, alpha = 0.05, n_perm = 10000)

# Joining with lookup_table_df and rearranging columns while excluding 'name'
pairwise_results_all = pairwise_results_all %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

pairwise_results_fa = pairwise_results_fa %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

pairwise_results_md = pairwise_results_md %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

pairwise_results_sa = pairwise_results_sa %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

pairwise_results_ct = pairwise_results_ct %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

# Viewing a portion of the results
head(pairwise_results_all)
head(pairwise_results_fa)
head(pairwise_results_md)
head(pairwise_results_sa)
head(pairwise_results_ct)

# Gathering significant results for all features 
significant_ci_results_all = pairwise_results_all %>% filter(is_significant_adjusted)
significant_ci_results_fa  = pairwise_results_fa %>% filter(is_significant_adjusted)
significant_ci_results_md  = pairwise_results_md %>% filter(is_significant_adjusted)
significant_ci_results_sa  = pairwise_results_sa %>% filter(is_significant_adjusted)
significant_ci_results_ct  = pairwise_results_ct %>% filter(is_significant_adjusted)

# Viewing the updated pairwise results
head(significant_ci_results_all)
head(significant_ci_results_fa)
head(significant_ci_results_md)
head(significant_ci_results_sa)
head(significant_ci_results_ct)

#' Calculate Odds Ratio Differences Between Groups Using Permutation Testing
#'
#' This function compares the odds ratios between two depression labels for each matched feature.
#' It employs permutation testing to assess the statistical significance of the observed differences 
#' in odds ratios between the groups.
#'
#' @param merged_df A data frame with merged significant features from two groups.
#'                  It expects columns with odds ratios named 'odds_ratio_x' and 'odds_ratio_y' for the two groups.
#' @param group1_name The name of the first group.
#' @param group2_name The name of the second group.
#' @param num_permutations The number of permutations to use in the permutation test.
#' @param alpha The significance level to determine if the difference is statistically significant.
#' @return A data frame containing each feature, the difference in odds ratios, the odds ratios for each group,
#'         the calculated p-value, a boolean indicating significance, and the group names.
#' @export
calculate_or_difference = function(merged_df, group1_name, group2_name, num_permutations = 10000, alpha = 0.05) {
  
  # Check if merged_df is empty
  if (nrow(merged_df) == 0) {
    return(data.frame())
  }
  # Initialize an empty dataframe to store the results 
  or_diff_df = data.frame(feature = character(),
                          odds_ratio_difference = numeric(),
                          odds_ratio_group1 = numeric(),
                          odds_ratio_group2 = numeric(),
                          p_value = numeric(),
                          is_significant = logical(),
                          group1 = character(),
                          group2 = character(),
                          stringsAsFactors = FALSE)
  
  # Loop over each row in the merged dataframe
  for (i in 1:nrow(merged_df)) {
    row = merged_df[i, ]
    # Calculate the absolute difference in odds ratios between the two groups
    observed_difference = abs(row$odds_ratio_x - row$odds_ratio_y)
    
    # Initialize a vector to store differences from each permutation
    perm_diffs = numeric(num_permutations)
    for (j in 1:num_permutations) {
      
      # Randomly permute odds ratios and calculate the absolute difference
      permuted_odds_ratio_x = sample(merged_df$odds_ratio_x, 1)
      permuted_odds_ratio_y = sample(merged_df$odds_ratio_y, 1)
      perm_diffs[j] = abs(permuted_odds_ratio_x - permuted_odds_ratio_y)
    }
    
    # Calculate the p-value based on how often the permuted differences are greater than or equal to the observed difference
    p_value = mean(perm_diffs >= observed_difference)
    
    # Replace zero p-values with the minimum non-zero value
    p_values_corrected = ifelse(p_value == 0, 1/(num_permutations + 1), p_value)
    
    # Determine if the observed difference is statistically significant
    is_significant = p_value < alpha
    
    # Append the results for this feature to the results dataframe
    or_diff_df = rbind(or_diff_df, data.frame(feature = row$feature,
                                              odds_ratio_difference = observed_difference,
                                              odds_ratio_group1 = row$odds_ratio_x,
                                              odds_ratio_group2 = row$odds_ratio_y,
                                              p_value = p_value,
                                              is_significant = is_significant,
                                              group1 = group1_name,
                                              group2 = group2_name))
  }
  
  return(or_diff_df)
}

#' Pairwise Comparison of Odds Ratio Differences Across Groups
#'
#' This function performs pairwise comparisons of odds ratio differences across various groups.
#' It leverages the calculate_or_difference function for each pair of groups.
#'
#' @param results_list A list containing data frames of significant features for each group.
#' @param num_permutations The number of permutations to use in the permutation test.
#' @param alpha The significance level to determine if the difference is statistically significant.
#' @return A list where each element is a data frame containing all pairwise comparison results for a pair of groups.
calculate_or_difference_pairwise = function(results_list, num_permutations = 10000, alpha = 0.05) {
  
  # Initialize a list to store the results of pairwise comparisons
  pairwise_results = data.frame()
  
  # Extract the names of the groups from the results list
  group_names = names(results_list)
  
  # Iterate over each unique pair of groups
  for (i in 1:(length(group_names) - 1)) {
    for (j in (i + 1):length(group_names)) {
      # Extract the names of the two groups being compared
      group1 = group_names[i]
      group2 = group_names[j]
      
      # Retrieve the data frames of significant features for each group
      df1 = results_list[[group1]]$significant_features
      df2 = results_list[[group2]]$significant_features
      
      # Check if df1 and df2 are not NULL and have rows
      # Proceed only if both groups have at least one significant feature
      if (nrow(df1) > 0 && nrow(df2) > 0) {
        # Merge the data frames on the 'feature' column with appropriate suffixes for each group
        merged_df = merge(df1, df2, by = "feature", suffixes = c("_x", "_y"))
        
        # Calculate the odds ratio differences for the merged data frame
        or_diff_results = calculate_or_difference(merged_df, group1, group2, num_permutations, alpha)
        
        # Store the results in the list with a key representing the pair of groups
        pairwise_results = rbind(pairwise_results, or_diff_results)
      }
    }
  }
  
  # Apply FDR correction 
  p_values = pairwise_results$p_value
  adjusted_p_values = p.adjust(p_values, method = "fdr", n = length(p_values))
  pairwise_results$adjusted_p_values = adjusted_p_values
  pairwise_results$is_significant = adjusted_p_values <= alpha
  
  
  return(pairwise_results)
}

# Calculate the OR differences pairwise across all groups with a specific alpha level
pairwise_or_results_all = calculate_or_difference_pairwise(sensitivity_significant_results_all, 10000, 0.05)
pairwise_or_results_global = calculate_or_difference_pairwise(sensitivity_significant_results_global_volume, 10000, 0.05)
pairwise_or_results_subcortical = calculate_or_difference_pairwise(sensitivity_significant_results_subcortical_volume, 10000, 0.05)
pairwise_or_results_thickness = calculate_or_difference_pairwise(sensitivity_significant_results_thickness, 10000, 0.05)
pairwise_or_results_sa = calculate_or_difference_pairwise(sensitivity_significant_results_surface_area, 10000, 0.05)
pairwise_or_results_fa = calculate_or_difference_pairwise(sensitivity_significant_results_fa, 10000, 0.05)
pairwise_or_results_md = calculate_or_difference_pairwise(sensitivity_significant_results_md, 10000, 0.05)

# Gathering significant results for all features 
# Joining with lookup_table_df and rearranging columns while excluding 'name'
pairwise_or_results_all = pairwise_or_results_all %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

pairwise_or_results_global = pairwise_or_results_global %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

pairwise_or_results_subcortical = pairwise_or_results_subcortical %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

pairwise_or_results_thickness = pairwise_or_results_thickness %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

pairwise_or_results_sa = pairwise_or_results_sa %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

pairwise_or_results_fa = pairwise_or_results_fa %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

pairwise_or_results_md = pairwise_or_results_md %>%
  left_join(lookup_table_df %>% select(-name, -id), by = c("feature" = "corrected_id")) %>%
  select(feature, processed_name_abbreviated, everything())

significant_ors_results_all         = pairwise_or_results_all %>% filter(is_significant)
significant_ors_results_global      = pairwise_or_results_global %>% filter(is_significant)
significant_ors_results_subcortical = pairwise_or_results_subcortical  %>% filter(is_significant)
significant_ors_results_thickness   = pairwise_or_results_thickness %>% filter(is_significant)
significant_ors_results_sa          = pairwise_or_results_sa  %>% filter(is_significant)
significant_ors_results_fa          = pairwise_or_results_fa %>% filter(is_significant)
significant_ors_results_md          = pairwise_or_results_md %>% filter(is_significant)


# Viewing pairwise results
head(significant_ors_results_all)
head(significant_ors_results_sa)
head(significant_ors_results_global)
head(significant_ors_results_subcortical)
head(significant_ors_results_thickness)
head(significant_ors_results_fa)
head(significant_ors_results_md)

cat("Done. Proceed to S3\n")