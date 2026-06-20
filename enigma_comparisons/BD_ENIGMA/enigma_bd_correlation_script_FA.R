# --------------------------------------------------------
# ENIGMA - Calculating correlations with BD working group summary statistics (for FA) 
# --------------------------------------------------------

###### Package setup ######

# Function to install and load packages
package_install = function(package_list) {
  for (pkg in package_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      tryCatch({
        install.packages(pkg)
      }, error = function(e) {
        message(paste0("Could not install ", pkg, ": ", e))
      })
    }
    # Load the library
    library(pkg, character.only = TRUE)
  }
}

# List of packages
packages = c('dplyr', 'data.table', 'readr', 'stringr', 'ggplot2', 'readxl',
             'tidyr', 'pheatmap', 'ggrepel', 'gridExtra', 'here')

# Call the package_install function
package_install(packages)

###### ENIGMA BD effect sizes and lookup table ######

# Create the BD data frame with category information
BD_FA = data.frame(
  ROI = c("PTR", "ACR", "CR", "ALIC", "PCR", "SCR", "IC", "RLIC", "CST", "PLIC",
          "CGC", "UNC", "EC", "SLF", "SS", "IFO", "FXST", "SFO", "CGH",
          "CC", "BCC", "GCC", "SCC", "FX"),
  Full.tract.name = c("Posterior thalamic radiation", "Anterior corona radiata", 
                      "Corona radiata", "Anterior limb of internal capsule",
                      "Posterior corona radiata", "Superior corona radiata",
                      "Internal capsule", "Retrolenticular part of IC",
                      "Corticospinal tract", "Posterior limb of internal capsule",
                      "Cingulum (cingulate)", "Uncinate fasciculus",
                      "External capsule", "Superior longitudinal fasciculus",
                      "Sagittal stratum", "Inferior fronto-occipital fasciculus",
                      "Fornix/Stria terminalis", "Superior fronto-occipital fasciculus",
                      "Cingulum (hippocampal)", "Corpus callosum",
                      "Body of corpus callosum", "Genu of corpus callosum",
                      "Splenium of corpus callosum", "Fornix"),
  Cohens_D = c(-0.3028, -0.2431, -0.1969, -0.1532, -0.1524, -0.086, -0.0665, -0.0471,
               0.0002, 0.0363, -0.3892, -0.2479, -0.2298, -0.2266, -0.1954, -0.1885,
               -0.156, -0.1548, -0.0678, -0.4625, -0.4297, -0.3729, -0.3387, -0.2876),
  Category = c(rep("Projection fibers", 10), 
               rep("Association fibers", 9),
               rep("Commissural fibers", 5))
)

# Load results lists and matching tables for white matter analysis
# Cohen's d from the linear regressions
load(here("data", "results_lists_FA_MD_linear_regression.RData"))
matching_table_FA = read_excel(here("data", "enigma_sumstats", "mdd", "matching_table_UKBB_to_Schmaal_FA_MD.xlsx"),
                               sheet = 1)

# use per-feature Cohen's d from the results lists

# Function to average hemispheric effects
average_hemispheric_effects = function(results_list, matching_table) {
  averaged_results = data.frame()
  
  for (row in 1:nrow(matching_table)) {
    # map region to its left/right UKBB features
    region = matching_table$Region[row]
    left_feature = matching_table$left_UKBB_feature_name[row]
    right_feature = matching_table$right_UKBB_feature_name[row]
    full_tract_name = matching_table$`Full tract name`[row]

    if (!is.na(left_feature) && !is.na(right_feature)) {
      # pull the left and right hemisphere results
      left_results = results_list[results_list$feature == left_feature, ]
      right_results = results_list[results_list$feature == right_feature, ]

      if (nrow(left_results) > 0 && nrow(right_results) > 0) {
        # average effect sizes across hemispheres
        averaged_odds_ratio = mean(c(left_results$odds_ratio, right_results$odds_ratio))
        cohens_d = mean(c(left_results$cohen_d_corrected, right_results$cohen_d_corrected))
        # bilateral SE and 95% CI
        se_bilateral = sqrt(left_results$cohen_d_se^2 + right_results$cohen_d_se^2) / 2
        cohens_d_lower_CI = cohens_d - 1.96 * se_bilateral
        cohens_d_upper_CI = cohens_d + 1.96 * se_bilateral
        averaged_fdr_p_value = mean(c(left_results$fdr_p_value, right_results$fdr_p_value))

        # assemble the averaged row
        averaged_results = rbind(averaged_results, data.frame(
          Region = region,
          `Full tract name` = full_tract_name,
          Averaged_Odds_Ratio = averaged_odds_ratio,
          Cohens_D = cohens_d,
          Cohens_D_Lower_CI = cohens_d_lower_CI,
          Cohens_D_Upper_CI = cohens_d_upper_CI,
          FDR_P_Value = averaged_fdr_p_value
        ))
      }
    }
  }
  return(averaged_results)
}

# Function to get top 5 features by absolute Cohen's D
get_top_5_meta = function(df) {
  df %>%
    arrange(desc(abs(Cohens_D))) %>%
    slice_head(n = 5) %>%
    rename(d_enigma = Cohens_D)
}

# --------------------------------------------------------
# Process white matter tract data and create comparison dataframes
# --------------------------------------------------------

###### White matter tract processing ######

# Function to process white matter tract data
process_white_matter_data = function(data, measure_name) {
  # Extract unique tracts
  tracts = unique(data$ROI)
  
  # Create measurement list
  measurement_list = paste0(measure_name, "_", tolower(tracts))
  
  # Compare lengths
  cat("Number of", measure_name, "variables:", length(measurement_list), "\n")
  cat("Number of Tract names:", length(tracts), "\n")
  
  # Compare order
  tracts_lower = tolower(tracts)
  measurement_lower = tolower(gsub(paste0("^", tolower(measure_name), "_"), "", measurement_list))
  
  order_match = all(tracts_lower == measurement_lower)
  cat("Order matches:", order_match, "\n")
  
  # Print mismatches if any
  if (!order_match) {
    mismatches = which(tracts_lower != measurement_lower)
    cat("Mismatches at positions:", paste(mismatches, collapse = ", "), "\n")
    cat("Tract names at mismatches:", paste(tracts[mismatches], collapse = ", "), "\n")
    cat("Measurement names at mismatches:", paste(measurement_list[mismatches], collapse = ", "), "\n")
  }
  
  # Create a data frame with Tract and effect size
  result = data.frame(
    Tract = tracts,
    Effect_Size = data$Cohens_D
  )
  
  return(result)
}

# Process data for FA
cat("Fractional Anisotropy:\n")
FA_processed = process_white_matter_data(BD_FA, "fa")
print(head(FA_processed))

# Create comparison dataframe for FA
create_comparison_dataframe_fa = function(fa_data) {
  # Extract Tract names
  tract_names = fa_data$Tract
  
  # Create list of processed name abbreviated
  processed_name_abbreviated = sapply(tract_names, function(name) {
    # Convert to lowercase and append fa_
    paste0("fa_", tolower(name))
  })
  
  # Create dataframe
  comparison_df = data.frame(
    Tract = tract_names,
    processed_name_abbreviated = processed_name_abbreviated,
    Effect_Size = fa_data$Effect_Size,
    stringsAsFactors = FALSE
  )
  
  return(comparison_df)
}

# Create comparison dataframe
cat("\nCreating comparison dataframe:\n")
comparison_df_FA = create_comparison_dataframe_fa(FA_processed)

# Display first few rows of the comparison dataframe
cat("\nFirst few rows of the FA comparison dataframe:\n")
print(head(comparison_df_FA, 10))

# Check if the number of rows match
cat("\nNumber of rows in the FA comparison dataframe:", nrow(comparison_df_FA), "\n")
cat("Number of tracts in FA_processed:", nrow(FA_processed), "\n")

#------------
# FA OBJECTS
#------------
cat("\nTract names:\n")
cat(paste(comparison_df_FA$Tract, collapse = "\n"), "\n\n")

cat("FA Variable names:\n")
cat(paste(comparison_df_FA$processed_name_abbreviated, collapse = "\n"), "\n\n")

cat("Effect Size Variable names:\n")
cat(paste(comparison_df_FA$Effect_Size, collapse = "\n"), "\n\n")

cat("Effect Sizes:\n")
cat(paste(sprintf("%.3f", comparison_df_FA$Effect_Size), collapse = "\n"), "\n\n")

# --------------------------------------------------------
# Processing Lookup Tables and Matching Names
# --------------------------------------------------------

###### Lookup tables and name matching ######

# Process lookup table to match ENIGMA BD naming conventions
matching_table_FA$Region = toupper(matching_table_FA$Region)

# Create mapping between ENIGMA BD and UKBB features
BD_FA$processed_name_abbreviated = paste0("fa_", tolower(BD_FA$ROI))

# Merge BD FA data with matching table
BD_FA_matched = merge(BD_FA, 
                      matching_table_FA[, c("Region", "left_UKBB_feature_name", "right_UKBB_feature_name")],
                      by.x = "ROI",
                      by.y = "Region",
                      all.x = TRUE)

# use per-feature Cohen's d
results_list_FA = lapply(results_list_FA, function(df) {
  df$cohen_d       = df$cohen_d_corrected
  df$cohen_d_lower = df$cohen_d_lower
  df$cohen_d_upper = df$cohen_d_upper
  return(df)
})


# Process and merge results 
results_list_FA_merged = lapply(results_list_FA, function(df) {
  # Average the hemispheric effects
  averaged_results = average_hemispheric_effects(df, matching_table_FA)
  
  # Merge with BD data
  merged_results = merge(averaged_results,
                         BD_FA[, c("ROI", "Cohens_D", "Category")],
                         by.x = "Region",
                         by.y = "ROI",
                         all.x = TRUE)
  
  return(merged_results)
})

# rename columns
results_list_FA_final = lapply(results_list_FA_merged, function(df) {
  df %>%
    rename(
      cohen_d = Cohens_D.x,    # UKBB effect sizes
      d_enigma = Cohens_D.y    # ENIGMA effect sizes
    ) %>%
    filter(!is.na(Region)) %>%
    arrange(Region)
})


# Print summary of matched results
cat("\nSummary of matched results:\n")
for(pheno in names(results_list_FA_final)) {
  cat("\nPhenotype:", pheno, "\n")
  cat("Number of matched tracts:", nrow(results_list_FA_final[[pheno]]), "\n")
}

# Save the matched results if needed
# save(results_list_FA_final, file = "results_list_FA_final.RData")

# --------------------------------------------------------
# Calculate correlations between BD ENIGMA and UKBB effect sizes
# --------------------------------------------------------

###### Correlations: ENIGMA vs UKBB ######

calculate_correlations = function(df, col1, col2) {
  # pair up the two effect-size columns
  cor_data = data.frame(
    x = df[[col1]],
    y = df[[col2]]
  )

  # spearman rho, guard against degenerate input
  spearman_raw = tryCatch(
    cor(cor_data$x, cor_data$y, method = "spearman"),
    error = function(e) NA
  )

  # spearman p-value, non-exact for ties
  spearman_p = tryCatch(
    cor.test(cor_data$x, cor_data$y, method = "spearman", exact = FALSE)$p.value,
    error = function(e) NA
  )

  list(spearman_raw = spearman_raw, spearman_p = spearman_p)
}

# Calculating correlations for FA measures across phenotypes
correlations_FA_GPNoDep            = calculate_correlations(results_list_FA_final$GPNoDep, "cohen_d", "d_enigma")
correlations_FA_Psypsy             = calculate_correlations(results_list_FA_final$Psypsy, "cohen_d", "d_enigma")
correlations_FA_SelfRepDep         = calculate_correlations(results_list_FA_final$SelfRepDep, "cohen_d", "d_enigma")
correlations_FA_GPpsy              = calculate_correlations(results_list_FA_final$GPpsy, "cohen_d", "d_enigma")
correlations_FA_DepAll             = calculate_correlations(results_list_FA_final$DepAll, "cohen_d", "d_enigma")
correlations_FA_ICD10Dep           = calculate_correlations(results_list_FA_final$ICD10Dep, "cohen_d", "d_enigma")
correlations_FA_ICD10Dep_exclpsych = calculate_correlations(results_list_FA_final$ICD10Dep_exclpsych, "cohen_d", "d_enigma")
correlations_FA_LifetimeMDD        = calculate_correlations(results_list_FA_final$LifetimeMDD, "cohen_d", "d_enigma")
correlations_FA_MDDRecur           = calculate_correlations(results_list_FA_final$MDDRecur, "cohen_d", "d_enigma")

# Creating a list of dataframes for correlations FA
correlations_FA = list(
  GPNoDep = correlations_FA_GPNoDep,
  Psypsy = correlations_FA_Psypsy,
  SelfRepDep = correlations_FA_SelfRepDep,
  GPpsy = correlations_FA_GPpsy,
  DepAll = correlations_FA_DepAll,
  ICD10Dep = correlations_FA_ICD10Dep,
  ICD10Dep_exclpsych = correlations_FA_ICD10Dep_exclpsych,
  LifetimeMDD = correlations_FA_LifetimeMDD,
  MDDRecur = correlations_FA_MDDRecur
)

# Extract raw p-values for FA correlations
fa_pvals = unlist(lapply(correlations_FA, function(x) x$spearman_p))

# FDR correction for p-values
fa_fdr_pvals = p.adjust(fa_pvals, method = "fdr")

# Add FDR-corrected p-values to the correlation results
for (i in seq_along(correlations_FA)) {
  correlations_FA[[i]]$fdr_p = fa_fdr_pvals[i]
}

# Print correlation results
cat("\nCorrelation Results:\n")
for(pheno in names(correlations_FA)) {
  cat("\nPhenotype:", pheno)
  cat("\nSpearman correlation:", round(correlations_FA[[pheno]]$spearman_raw, 3))
  cat("\nRaw p-value:", format(correlations_FA[[pheno]]$spearman_p, digits = 3))
  cat("\nFDR-adjusted p-value:", format(correlations_FA[[pheno]]$fdr_p, digits = 3))
  cat("\n---")
}


# Extract top features based on absolute effect sizes
get_top_features = function(df, n = 5) {
  df %>%
    arrange(desc(abs(d_enigma))) %>%  # Sort by absolute ENIGMA effect size
    slice_head(n = n) %>%
    select(Region, Full.tract.name, d_enigma, cohen_d) %>%
    rename(
      ENIGMA_Effect = d_enigma,
      UKBB_Effect = cohen_d
    )
}


# Get top features for each phenotype
top_features_FA = lapply(results_list_FA_final, get_top_features)

# Print top features
cat("\nTop Features by Effect Size:\n")
for(pheno in names(top_features_FA)) {
  cat("\nPhenotype:", pheno, "\n")
  print(top_features_FA[[pheno]])
  cat("\n---")
}

# Get top features for each phenotype
top_features_FA = lapply(results_list_FA_final, get_top_features)

# get top features
get_top_features = function(df, n = 5) {
  df %>%
    arrange(desc(abs(d_enigma))) %>%
    slice_head(n = n) %>%
    select(Region, `Full.tract.name`, d_enigma, cohen_d) %>%
    rename(
      ENIGMA_Effect = d_enigma,
      UKBB_Effect = cohen_d
    )
}


# Now proceeding with visualization steps 9 and 10

###### Correlation heatmap ######

# Create correlation heatmap with fixed scale
plot_correlation_heatmap = function(correlations, title) {
  # Order definitions by phenotype granularity
  order_definitions = c("GPNoDep", "GPpsy", "Psypsy", "SelfRepDep", "DepAll", 
                        "ICD10Dep", "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur")
  
  # Convert correlations to a data frame if it's not already
  if (!is.data.frame(correlations)) {
    correlations_df = do.call(rbind, lapply(names(correlations), function(x) {
      data.frame(Definition = x, 
                 Spearman_Raw = correlations[[x]]$spearman_raw,
                 Spearman_P = correlations[[x]]$spearman_p,
                 FDR_P = correlations[[x]]$fdr_p,
                 stringsAsFactors = FALSE)
    }))
  } else {
    correlations_df = correlations
  }
  
  # Reorder the correlations based on the order of definitions
  # keep only known definitions, order by granularity
  correlations_df = correlations_df %>%
    filter(Definition %in% order_definitions) %>%
    arrange(factor(Definition, levels = order_definitions))

  # Create correlation matrix
  # single-row matrix of rho, phenotypes as columns
  corr_matrix = matrix(correlations_df$Spearman_Raw, nrow = 1)
  rownames(corr_matrix) = "ENIGMA BD"
  colnames(corr_matrix) = correlations_df$Definition
  
  # Create p-value matrix
  p_matrix = matrix(correlations_df$FDR_P, nrow = 1)
  rownames(p_matrix) = "ENIGMA BD"
  colnames(p_matrix) = correlations_df$Definition
  
  # Create labels
  labels = matrix(paste0(round(corr_matrix, 2), "\n(", format(p_matrix, digits = 3), ")"),
                  nrow = nrow(corr_matrix), ncol = ncol(corr_matrix))
  
  # Set fixed maximum correlation value for consistent scaling
  max_abs_corr = 0.8  # range -0.8 to 0.8

  # Create a custom color palette
  # blue-white-red diverging ramp
  custom_colors = colorRampPalette(c("#6666ff", "white", "#ff3333"))(100)

  # Add the scale limits to pheatmap
  # draw the heatmap with fixed breaks and annotated cells
  pheatmap(corr_matrix,
           main = title, 
           color = custom_colors,
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           display_numbers = labels, 
           number_format = "%.2f",
           fontsize_number = 8, 
           angle_col = 45, 
           cellwidth = 80, 
           cellheight = 35,
           fontcolor_number = "black", 
           border_color = "black",
           breaks = seq(-max_abs_corr, max_abs_corr, length.out = 101),
           legend_breaks = seq(-0.8, 0.8, by = 0.4),  # Add explicit legend breaks
           legend_labels = sprintf("%.1f", seq(-0.8, 0.8, by = 0.4)))  # Add explicit legend labels
}

# Generate heatmap for FA correlations
plot_correlation_heatmap(correlations_FA, "Spearman Correlations (Fractional Anisotropy)")

###### Scatter plots for top correlations ######

# Step 10: Scatter Plots for Top Correlations
get_top_correlations = function(correlations, n = 3) {
  corr_values = sapply(correlations, function(x) x$spearman_raw)
  top_indices = order(abs(corr_values), decreasing = TRUE)[1:n]
  names(correlations)[top_indices]
}

# Get top 3 strongest correlations
top_corr_FA = get_top_correlations(correlations_FA, n = 3)

# Calculate global ranges
calculate_global_range = function(data_list, top_corr) {
  all_x_values = unlist(lapply(top_corr, function(definition) data_list[[definition]]$cohen_d))
  all_y_values = unlist(lapply(top_corr, function(definition) data_list[[definition]]$d_enigma))
  x_range = range(all_x_values, na.rm = TRUE)
  y_range = range(all_y_values, na.rm = TRUE)
  list(x_range = x_range, y_range = y_range)
}

global_range_FA = calculate_global_range(results_list_FA_final, top_corr_FA)

# Create scatter plot function
create_scatter_plot = function(data, definition, corr, p_value, global_range, color) {
  # points plus linear fit, UKB vs ENIGMA d
  ggplot(data[[definition]], aes(x = cohen_d, y = d_enigma)) +
    geom_point(size = 2, alpha = 0.7, color = color) +
    geom_smooth(method = "lm", se = FALSE, color = color, linewidth = 1) +
    xlab("Cohen's d (UKB)") +
    ylab("Cohen's d (ENIGMA BD)") +
    # share axis ranges across plots for comparability
    scale_y_continuous(breaks = seq(min(global_range$y_range), max(global_range$y_range), length.out = 5),
                       limits = global_range$y_range,
                       expand = c(0.01, 0)) +
    scale_x_continuous(breaks = seq(min(global_range$x_range), max(global_range$x_range), length.out = 5),
                       limits = global_range$x_range,
                       labels = function(x) sprintf("%.3f", x),
                       expand = c(0.01, 0)) +
    theme_classic() +
    # match all text and axes to the plot color
    theme(
      axis.title = element_text(size = 12, color = color),
      axis.text = element_text(size = 10, color = color),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_line(color = color),
      axis.ticks = element_line(color = color),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 20, 10),
      plot.title = element_text(hjust = 0.5, size = 14, color = color, face = "bold")
    ) +
    # stamp r and p in the bottom-right corner
    annotate("text", x = Inf, y = -Inf,
             label = paste0("r = ", round(corr, 3), "\np = ", format.pval(p_value, digits = 3)),
             hjust = 1.1, vjust = -0.8, size = 3, color = color) +
    labs(title = definition)
}

# Generate scatter plots
plot_colors = c("#4a80ed", "#ed984a", "#45a049")

scatter_plots_FA = lapply(seq_along(top_corr_FA), function(i) {
  definition = top_corr_FA[i]
  corr = correlations_FA[[definition]]$spearman_raw
  p_value = correlations_FA[[definition]]$spearman_p
  print(paste("Plot", i, "corresponds to:", paste0("FA: ", definition)))
  create_scatter_plot(results_list_FA_final, definition, corr, p_value, 
                      global_range_FA, plot_colors[i])
})

# Arrange the FA plots in a grid layout
scatter_plot_grid_FA = do.call(grid.arrange, c(scatter_plots_FA, ncol = 3))

