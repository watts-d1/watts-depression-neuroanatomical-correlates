# --------------------------------------------------------
# ENIGMA - Calculating correlations with Schizophrenia working group summary statistics 
# --------------------------------------------------------

###### Package setup ######

# Function to install and load packages
package_install = function(package_list) {
  # install each package if missing, then load
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
             'tidyr', 'pheatmap', 'ggrepel', 'gridExtra', 'ggpubr', 'grid', 'here')

# Call the package_install function 
package_install(packages)

###### Load ENIGMA SCZ summary stats ######

# Loading csv files
CT   = fread(here("data", "enigma_sumstats", "schizophrenia", "schizophrenia_cortical_thickness.csv"))
SA   = fread(here("data", "enigma_sumstats", "schizophrenia", "schizophrenia_surface_area.csv"))

process_and_print_data = function(data, drop_temporalpole = FALSE) {
  # optionally drop temporal pole regions
  filtered_data = if (drop_temporalpole) {
    data[!endsWith(data$Structure, "_temporalpole"), ]
  } else {
    data
  }

  # build structure: effect-size strings
  combined_list = sprintf("%s: %.3f", filtered_data$Structure, filtered_data$d_icv)

  # echo structures and effect sizes to console
  cat("Combined list:\n", paste(combined_list, collapse = "\n"), "\n\n", sep = "")
  cat("Structure:\n",     paste(filtered_data$Structure, collapse = "\n"), "\n\n", sep = "")
  cat("d_icv:\n",         paste(sprintf("%.3f", filtered_data$d_icv), collapse = "\n"), "\n\n", sep = "")

  # return the pieces invisibly for downstream use
  invisible(list(combined_list = combined_list,
                 Structure     = filtered_data$Structure,
                 d_icv         = filtered_data$d_icv))
}

# Process and print data for SA and CT
cat("Surface Area:\n")
sa_results_ggseg = process_and_print_data(SA)

cat("\nCortical Thickness:\n")
ct_results_ggseg = process_and_print_data(CT)


print_results_sequentially = function(results) {
  cat("Combined list:\n")
  cat(paste(results$combined_list, collapse = "\n"), "\n\n")
  
  cat("Structure:\n")
  cat(paste(results$Structure, collapse = "\n"), "\n\n")
  
  cat("d_icv:\n")
  cat(paste(sprintf("%.3f", results$d_icv), collapse = "\n"), "\n\n")
}

# Print sa_results_ggseg sequentially
cat("Surface Area Results:\n")
print_results_sequentially(sa_results_ggseg)

# Print sa_results_ggseg sequentially
cat("Cortical Thickness Results:\n")
print_results_sequentially(ct_results_ggseg)

###### Load UKBB Cohen's d results lists ######

# Loading results lists (thickness & surface area)
# Cohen's d results lists
load(here("data", "results_lists_thickness_SA_linear_regression.RData"))

# Removing unused depression definitions 
results_list_surface_area = results_list_surface_area[!names(results_list_surface_area) %in% c("MDD_noimpairment", "MDD2_w_impair")]
results_list_thickness    = results_list_thickness[!names(results_list_thickness) %in% c("MDD_noimpairment", "MDD2_w_impair")]

###### Lookup table and effect-size helpers ######

# --------------------------------------------------------
# Extracting Feature Names
# --------------------------------------------------------

# Path to the neuroimaging lookup table
lookup_table_path = here("data", "lookup_table_neuroimaging_ukbb.xlsx")

# read the lookup table from excel
read_lookup_table = function(filepath) {
  read_excel(filepath)
}

# Lookup table mapping feature ids to names
lookup_table_df = read_lookup_table(lookup_table_path)

# Rename the prefixes "baseline_" to "bl_" for each feature value in the table
lookup_table_df$corrected_id = gsub("baseline_", "bl_", lookup_table_df$corrected_id)


# Matching table of UKBB field ids to sumstat Structure (region)
matching_table_path = here("data", "enigma_sumstats", "bipolar_disorder", "matching_table_bipolar_sumstats_Hibar2018.xlsx")
SA_matching_table = as.data.frame(read_excel(matching_table_path, sheet = 1))
CT_matching_table = as.data.frame(read_excel(matching_table_path, sheet = 2))

# Add the corrected_id to the SA dataframe based on the matching Structure column
SA$corrected_id = SA_matching_table$corrected_id[match(SA$Structure, SA_matching_table$Structure)]

# Add the corrected_id to the CT dataframe based on the matching Structure column
CT$corrected_id = CT_matching_table$corrected_id[match(CT$Structure, CT_matching_table$Structure)]

# Condense the "baseline_" prefix to "bl_" in the SA and CT dataframes
SA$corrected_id = gsub("baseline_", "bl_", SA$corrected_id)
CT$corrected_id = gsub("baseline_", "bl_", CT$corrected_id)

# Add the processed_name_abbreviated to the SA dataframe based on the matching Structure column
SA$processed_name_abbreviated = SA_matching_table$processed_name_abbreviated[match(SA$Structure, SA_matching_table$Structure)]

# Add the processed_name_abbreviated to the CT dataframe based on the matching Structure column
CT$processed_name_abbreviated = CT_matching_table$processed_name_abbreviated[match(CT$Structure, CT_matching_table$Structure)]
#stopifnot(all(grepl("^area_", SA$processed_name_abbreviated)),
#all(grepl("^mean_thickness_", CT$processed_name_abbreviated)))

# Function to add the corresponding d_icv column to each dataframe in results_list
add_d_icv = function(results_list, enigma_df) {
  lapply(results_list, function(df) {
    # join ENIGMA d_icv onto each result by feature id
    df$enigma_d = enigma_df$d_icv[match(df$feature, enigma_df$corrected_id)]
    df
  })
}

# Add cohen_d column to surface area results
results_list_surface_area_final = lapply(results_list_surface_area, function(df) {
  # Cohen's d from the linear regressions
  df$cohen_d = df$cohen_d_corrected
  df
})

# Add cohen_d column to thickness results
results_list_thickness_final = lapply(results_list_thickness, function(df) {
  # Cohen's d from the linear regressions
  df$cohen_d = df$cohen_d_corrected
  df
})

# Add the corresponding d_icv column to surface area results
results_list_surface_area_final = add_d_icv(results_list_surface_area_final, SA)

# Add the corresponding d_icv column to thickness results
results_list_thickness_final = add_d_icv(results_list_thickness_final, CT)

##---------------------------------
# Removing rows with missing values 
##---------------------------------
results_list_surface_area_final = lapply(results_list_surface_area_final, function(df) {
  df %>% 
    filter(!is.na(enigma_d))
})

results_list_thickness_final = lapply(results_list_thickness_final, function(df) {
  df %>% 
    filter(!is.na(enigma_d))
})

###### Spearman correlations and FDR ######

# Function to calculate Spearman's rank correlation
calculate_correlations = function(df) {
  # keep rows where both effect sizes are present
  ok = complete.cases(df$cohen_d, df$enigma_d)
  x = df$cohen_d[ok]; y = df$enigma_d[ok]
  # spearman rho and its p-value (UKB vs ENIGMA)
  list(spearman_raw = cor(x, y, method = "spearman"),
       spearman_p   = cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
}

# Calculating correlations for surface area measures across depression phenotypes 
correlations_surface_area = lapply(results_list_surface_area_final, calculate_correlations)

# Calculating correlations for thickness measures across depression phenotypes 
correlations_thickness = lapply(results_list_thickness_final, calculate_correlations)

# Extract raw p-values for surface area correlations
surface_area_pvals = sapply(correlations_surface_area, function(x) x$spearman_p)

# Extract raw p-values for cortical thickness correlations
thickness_pvals = sapply(correlations_thickness, function(x) x$spearman_p)

# FDR correction for surface area p-values
surface_area_fdr_pvals = p.adjust(surface_area_pvals, method = "fdr")

# FDR correction for cortical thickness p-values
thickness_fdr_pvals = p.adjust(thickness_pvals, method = "fdr")

# Add FDR-corrected p-values to the correlation results
for (i in seq_along(correlations_surface_area)) {
  correlations_surface_area[[i]]$fdr_p = surface_area_fdr_pvals[i]
}

for (i in seq_along(correlations_thickness)) {
  correlations_thickness[[i]]$fdr_p = thickness_fdr_pvals[i]
}

# --------------------------------------------------------
# Write cortical Spearman correlations to CSV
# per-definition cortical correlations
# --------------------------------------------------------
correlations_to_df = function(correlations, measure) {
  # one row per depression definition with r, p, and FDR q
  do.call(rbind, lapply(names(correlations), function(def) {
    data.frame(Definition   = def,
               Measure      = measure,
               Spearman_r   = correlations[[def]]$spearman_raw,
               p_value      = correlations[[def]]$spearman_p,
               FDR_p_value  = correlations[[def]]$fdr_p,
               stringsAsFactors = FALSE)
  }))
}

cortical_correlations_export = rbind(
  correlations_to_df(correlations_surface_area, "Surface Area"),
  correlations_to_df(correlations_thickness, "Thickness")
)

dir.create(here("output"), showWarnings = FALSE, recursive = TRUE)
write.csv(cortical_correlations_export,
          here("output", "scz_ctsa_cortical_enigma_correlations.csv"),
          row.names = FALSE)

###### Correlation heatmaps ######

# --------------------------------------------------------
# Step 9. Visualize Correlations
# --------------------------------------------------------

# Maximum absolute correlation across both datasets
max_abs_corr = max(abs(c(
  sapply(correlations_surface_area, function(x) x$spearman_raw),
  sapply(correlations_thickness, function(x) x$spearman_raw)
)))

# heatmap of Spearman correlations and p-values
plot_correlation_heatmap = function(correlations, title, max_abs_corr = NULL) {
  # Order definitions by phenotype granularity
  order_definitions = c("GPNoDep", "GPpsy", "Psypsy", "SelfRepDep", "DepAll", "ICD10Dep", "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur")
  
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
  # keep only ordered definitions and sort them
  correlations_df = correlations_df %>%
    filter(Definition %in% order_definitions) %>%
    arrange(factor(Definition, levels = order_definitions))

  # Create correlation matrix
  corr_matrix = matrix(correlations_df$Spearman_Raw, nrow = 1)
  rownames(corr_matrix) = "ENIGMA Schizophrenia"
  colnames(corr_matrix) = correlations_df$Definition
  
  # Create p-value matrix
  p_matrix = matrix(correlations_df$FDR_P, nrow = 1)
  rownames(p_matrix) = "ENIGMA Schizophrenia"
  colnames(p_matrix) = correlations_df$Definition
  
  # Create labels
  # annotate each cell with rho and its FDR q
  labels = matrix(paste0(round(corr_matrix, 2), "\n(", format(p_matrix, digits = 3), ")"),
                  nrow = nrow(corr_matrix), ncol = ncol(corr_matrix))

  # Determine color range
  if (is.null(max_abs_corr)) {
    max_abs_corr = max(abs(corr_matrix))
  }
  
  # Create a custom color palette
  # blue-white-red diverging scale centered on zero
  custom_colors = colorRampPalette(c("#6666ff", "white", "#ff3333"))(100)

  # Create the heatmap
  # draw single-row heatmap with symmetric color breaks
  pheatmap(corr_matrix, main = title, color = custom_colors,
           cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = labels, number_format = "%.2f",
           fontsize_number = 8, angle_col = 45, cellwidth = 80, cellheight = 35,
           fontcolor_number = "black", border_color = "black",
           breaks = seq(-max_abs_corr, max_abs_corr, length.out = 101))
}

# Create heatmaps
plot_correlation_heatmap(correlations_surface_area, "Spearman Correlations (Surface Area)", max_abs_corr)
plot_correlation_heatmap(correlations_thickness, "Spearman Correlations (Thickness)", max_abs_corr)

###### Spot-check thickness correlations ######

# --------------------------------------------------------
# check thickness correlations
# correlation between cohen_d and enigma_d
# --------------------------------------------------------
head(results_list_thickness_final$ICD10Dep$cohen_d)
head(results_list_thickness_final$ICD10Dep$enigma_d)

# Calculate Spearman correlation
cor_result = cor.test(results_list_thickness_final$ICD10Dep$cohen_d, 
                      results_list_thickness_final$ICD10Dep$enigma_d, 
                      method = "spearman", 
                      exact = FALSE)

# Print the result
print(cor_result)

head(results_list_thickness_final$ICD10Dep_exclpsych$cohen_d)
head(results_list_thickness_final$ICD10Dep_exclpsych$enigma_d)

# Calculate Spearman correlation
cor_result = cor.test(results_list_thickness_final$ICD10Dep_exclpsych$cohen_d, 
                      results_list_thickness_final$ICD10Dep_exclpsych$enigma_d, 
                      method = "spearman", 
                      exact = FALSE)

# Print the result
print(cor_result)

###### Scatter plots for top correlations ######

# --------------------------------------------------------
# Step 10. Scatter Plots for Top Correlations
# --------------------------------------------------------

# get the top n definitions by absolute Spearman correlation
get_top_correlations = function(correlations, n = 3) {
  corr_values = sapply(correlations, function(x) abs(x$spearman_raw))
  top_indices = order(corr_values, decreasing = TRUE)[1:n]
  names(correlations)[top_indices]
}

# Get the top 3 strongest correlations for surface area and thickness
top_corr_surface_area = get_top_correlations(correlations_surface_area, n = 3)
top_corr_thickness = get_top_correlations(correlations_thickness, n = 3)

# Define colors
plot_colors = c("#4a80ed", "#ed984a", "#45a049")

# Calculate global ranges for surface area and thickness
calculate_global_range = function(data_list, top_corr) {
  # pool x/y across the top definitions for shared axes
  all_x_values = unlist(lapply(top_corr, function(definition) data_list[[definition]]$cohen_d))
  all_y_values = unlist(lapply(top_corr, function(definition) data_list[[definition]]$enigma_d))
  # common axis ranges so panels are comparable
  x_range = range(all_x_values, na.rm = TRUE)
  y_range = range(all_y_values, na.rm = TRUE)
  list(x_range = x_range, y_range = y_range)
}

global_range_surface_area = calculate_global_range(results_list_surface_area_final, top_corr_surface_area)
global_range_thickness = calculate_global_range(results_list_thickness_final, top_corr_thickness)

# scatter plot of UKB cohen_d vs ENIGMA enigma_d for a definition
create_scatter_plot = function(data, definition, corr, p_value, measure, global_range, color) {
  # drop incomplete pairs before plotting
  data = data[complete.cases(data[, c("cohen_d", "enigma_d")]), ]
  pear = cor(data$cohen_d, data$enigma_d, method = "pearson")

  # points plus linear fit of UKB vs ENIGMA
  ggplot(data, aes(x = cohen_d, y = enigma_d)) +
    geom_point(size = 2, alpha = 0.7, color = color) +
    geom_smooth(method = "lm", se = FALSE, color = color, linewidth = 1) +
    xlab("UKB association (log-OR per SD, rescaled)") +
    ylab("Effect Size (ENIGMA)") +
    # lock axes to the shared global range
    scale_y_continuous(breaks = seq(min(global_range$y_range), max(global_range$y_range), length.out = 5),
                       limits = global_range$y_range,
                       expand = c(0.01, 0)) +
    scale_x_continuous(breaks = seq(min(global_range$x_range), max(global_range$x_range), length.out = 5),
                       limits = global_range$x_range,
                       labels = function(x) sprintf("%.3f", x),
                       expand = c(0.01, 0)) +
    theme_classic() +
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
    # annotate corner with rho and FDR q
    annotate("text", x = Inf, y = -Inf,
             label = paste0("r = ", round(corr, 3), "\nq = ", format.pval(p_value, digits = 3)),
             hjust = 1.1, vjust = -0.8, size = 3, color = color) +
    labs(title = definition)
}

# Scatter plots for surface area
scatter_plots_surface_area = lapply(seq_along(top_corr_surface_area), function(i) {
  definition = top_corr_surface_area[i]
  corr = correlations_surface_area[[definition]]$spearman_raw
  p_value = correlations_surface_area[[definition]]$fdr_p
  print(paste("Plot", i, "corresponds to:", paste0("Surface Area: ", definition)))
  create_scatter_plot(results_list_surface_area_final[[definition]], definition, corr, p_value, 
                      "Surface Area", global_range_surface_area, plot_colors[i])
})

# Arrange the surface area plots in a grid layout
# 1107 x 417 
# SCZ_Enigma_scatterplot_surface_area_top3
scatter_plot_grid_surface_area = do.call(grid.arrange, c(scatter_plots_surface_area, ncol = 3))

# Scatter plots for thickness
scatter_plots_thickness = lapply(seq_along(top_corr_thickness), function(i) {
  definition = top_corr_thickness[i]
  corr = correlations_thickness[[definition]]$spearman_raw
  p_value = correlations_thickness[[definition]]$fdr_p
  print(paste("Plot", i, "corresponds to:", paste0("Thickness: ", definition)))
  create_scatter_plot(results_list_thickness_final[[definition]], definition, corr, p_value, 
                      "Thickness", global_range_thickness, plot_colors[i])
})

# Arrange the thickness plots in a grid layout
# 1107 x 417
scatter_plot_grid_thickness = do.call(grid.arrange, c(scatter_plots_thickness, ncol = 3))
