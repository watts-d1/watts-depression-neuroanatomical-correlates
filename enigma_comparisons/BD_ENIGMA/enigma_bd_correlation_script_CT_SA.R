# --------------------------------------------------------
# ENIGMA - Calculating correlations with Bipolar working group summary statistics
# --------------------------------------------------------

###### Package setup and data loading ######

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
             'tidyr', 'pheatmap', 'ggrepel', 'gridExtra', 'ggpubr', 'grid', 'here')

# Call the package_install function 
package_install(packages)

# Loading csv files
CT   = fread(here("data", "enigma_sumstats", "bipolar_disorder", "CT.csv"))
SA   = fread(here("data", "enigma_sumstats", "bipolar_disorder", "SA.csv"))

###### Region/effect-size processing helpers ######

# Function to return Structure and corresponding effect size
process_brain_structure_data = function(data, measure_name) {
  # Extract unique structures
  structures = unique(data$Structure)
  
  # Create measurement list
  measurement_list = paste0(measure_name, "_", tolower(gsub("^[LR]_", "", structures)))
  
  # Compare lengths
  cat("Number of", measure_name, "variables:", length(measurement_list), "\n")
  cat("Number of Structure names:", length(structures), "\n")
  
  # Compare order
  structures_lower = tolower(structures)
  measurement_lower = tolower(gsub(paste0("^", tolower(measure_name), "_"), "", measurement_list))
  
  order_match = all(structures_lower == measurement_lower)
  cat("Order matches:", order_match, "\n")
  
  # Print mismatches if any
  if (!order_match) {
    mismatches = which(structures_lower != measurement_lower)
    cat("Mismatches at positions:", paste(mismatches, collapse = ", "), "\n")
    cat("Structure names at mismatches:", paste(structures[mismatches], collapse = ", "), "\n")
    cat("Measurement names at mismatches:", paste(measurement_list[mismatches], collapse = ", "), "\n")
  }
  
  # Create a data frame with Structure and effect size
  result = data.frame(
    Structure = structures,
    Effect_Size = data$d_icv
  )
  
  # Return the result
  return(result)
}

# Process and print data for SA and CT
cat("Surface Area:\n")
SA_processed_for_ggseg = process_brain_structure_data(SA, "surface_area")
print(head(SA_processed_for_ggseg))

cat("\nCortical Thickness:\n")
CT_processed_for_ggseg = process_brain_structure_data(CT, "cortical_thickness")
print(head(CT_processed_for_ggseg))

# Process and print data for SA and CT
cat("Surface Area:\n")
SA_processed_for_ggseg = process_brain_structure_data(SA, "surface_area")
print(head(SA_processed_for_ggseg))

cat("\nCortical Thickness:\n")
CT_processed_for_ggseg = process_brain_structure_data(CT, "cortical_thickness")
print(head(CT_processed_for_ggseg))

create_comparison_dataframe_ct = function(ct_data) {
  # Extract Structure names from ct_data
  structure_names = ct_data$Structure
  
  # Create list of processed name abbreviated
  processed_name_abbreviated = sapply(structure_names, function(name) {
    # Extract the hemisphere (L or R) and the structure name
    hemisphere = substr(name, 1, 1)
    structure = substr(name, 3, nchar(name))
    
    # Convert hemisphere to lf or rf
    hemisphere_suffix = ifelse(hemisphere == "L", "_lf", "_rf")
    
    # Construct the processed name
    paste0("mean_thickness_", tolower(structure), hemisphere_suffix)
  })
  
  # Create dataframe
  comparison_df = data.frame(
    Structure = structure_names,
    processed_name_abbreviated = processed_name_abbreviated,
    Effect_Size = ct_data$Effect_Size,
    stringsAsFactors = FALSE
  )
  
  return(comparison_df)
}

create_comparison_dataframe_sa = function(sa_data) {
  # Extract Structure names from sa_data
  structure_names = sa_data$Structure
  
  # Create list of processed name abbreviated
  processed_name_abbreviated = sapply(structure_names, function(name) {
    # Extract the hemisphere (L or R) and the structure name
    hemisphere = substr(name, 1, 1)
    structure = substr(name, 3, nchar(name))
    
    # Convert hemisphere to lf or rf
    hemisphere_suffix = ifelse(hemisphere == "L", "_lf", "_rf")
    
    # Construct the processed name
    paste0("area_", tolower(structure), hemisphere_suffix)
  })
  
  # Create dataframe
  comparison_df = data.frame(
    Structure = structure_names,
    processed_name_abbreviated = processed_name_abbreviated,
    Effect_Size = sa_data$Effect_Size,
    stringsAsFactors = FALSE
  )
  
  return(comparison_df)
}

# After processing SA and CT data
cat("\nCreating comparison dataframe:\n")
comparison_df_CT = create_comparison_dataframe_ct(CT_processed_for_ggseg)
comparison_df_SA = create_comparison_dataframe_sa(SA_processed_for_ggseg)

# Display first few rows of the comparison dataframe
cat("\nFirst few rows of the comparison dataframe:\n")
print(tail(comparison_df_CT, 10))
print(tail(comparison_df_SA, 10))

# Check if the number of rows match
cat("\nNumber of rows in comparison dataframe:", nrow(comparison_df_CT), "\n")
cat("Number of structures in CT_processed_for_ggseg:", length(CT_processed_for_ggseg), "\n")

#------------
# CT OBJECTS
#------------
cat("\nStructure names:\n")
cat(paste(comparison_df_CT$Structure, collapse = "\n"), "\n\n")

cat("Mean Thickness Variable names:\n")
cat(paste(comparison_df_CT$processed_name_abbreviated, collapse = "\n"), "\n\n")

cat("Effect Size Variable names:\n")
cat(paste(comparison_df_CT$Effect_Size, collapse = "\n"), "\n\n")

cat("Effect Sizes:\n")
cat(paste(sprintf("%.3f", comparison_df_SA$Effect_Size), collapse = "\n"), "\n\n")


#------------
# SA OBJECTS
#------------
cat("\nStructure names:\n")
cat(paste(comparison_df_SA$Structure, collapse = "\n"), "\n\n")

cat("Mean Thickness Variable names:\n")
cat(paste(comparison_df_SA$processed_name_abbreviated, collapse = "\n"), "\n\n")

cat("Effect Size Variable names:\n")
cat(paste(comparison_df_SA$Effect_Size, collapse = "\n"), "\n\n")

cat("Effect Sizes:\n")
cat(paste(sprintf("%.3f", comparison_df_SA$Effect_Size), collapse = "\n"), "\n\n")


###### UKB lookup table and feature matching ######

# Loading results lists (thickness & surface area)
load(here("data", "results_lists_thickness_SA_linear_regression.RData"))

# Removing unused depression definitions 
results_list_surface_area = results_list_surface_area[!names(results_list_surface_area) %in% c("MDD_noimpairment", "MDD2_w_impair")]
results_list_thickness    = results_list_thickness[!names(results_list_thickness) %in% c("MDD_noimpairment", "MDD2_w_impair")]

# --------------------------------------------------------
# Extracting Feature Names
# --------------------------------------------------------

# Path to the neuroimaging lookup table
lookup_table_path = here("data", "lookup_table_neuroimaging_ukbb.xlsx")

# Read the lookup table from excel
read_lookup_table = function(filepath) {
  read_excel(filepath)
}

# Read the lookup table into a dataframe
lookup_table_df = read_lookup_table(lookup_table_path)

# Rename the prefixes "baseline_" to "bl_" for each feature value in the table
lookup_table_df$corrected_id = gsub("baseline_", "bl_", lookup_table_df$corrected_id)


# Matching table of UKBB field IDs to sumstat Structure (region)
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


# Function to add the corresponding d_icv column to each dataframe in results_list
add_d_icv = function(results_list, enigma_df) {
  lapply(results_list, function(df) {
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

###### Spearman correlations and FDR correction ######

# Function to calculate Spearman's rank correlation
calculate_correlations = function(df) {
  # spearman rho between UKB and ENIGMA effect sizes
  spearman_raw = tryCatch(
    cor(df$cohen_d, df$enigma_d, method = "spearman"),
    error = function(e) NA
  )
  # two-sided p-value, non-exact to allow ties
  spearman_p = tryCatch(
    cor.test(df$cohen_d, df$enigma_d, method = "spearman", exact = FALSE)$p.value,
    error = function(e) NA
  )
  list(spearman_raw = spearman_raw, spearman_p = spearman_p)
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
# --------------------------------------------------------
correlations_to_df = function(correlations, measure_label) {
  # flatten the per-definition correlation list into one tidy frame
  do.call(rbind, lapply(names(correlations), function(def) {
    data.frame(
      Definition   = def,
      Measure      = measure_label,
      Spearman_r   = correlations[[def]]$spearman_raw,
      p_value      = correlations[[def]]$spearman_p,
      FDR_p_value  = correlations[[def]]$fdr_p,
      stringsAsFactors = FALSE
    )
  }))
}

cortical_correlations_corrected = rbind(
  correlations_to_df(correlations_surface_area, "Surface Area"),
  correlations_to_df(correlations_thickness, "Cortical Thickness")
)

write.csv(
  cortical_correlations_corrected,
  here("output", "bd_ctsa_cortical_enigma_correlations.csv"),
  row.names = FALSE
)

###### Correlation heatmaps ######

# --------------------------------------------------------
# Step 9. Visualize Correlations
# --------------------------------------------------------

# Calculating the maximum absolute correlation across both datasets
max_abs_corr = max(abs(c(
  sapply(correlations_surface_area, function(x) x$spearman_raw),
  sapply(correlations_thickness, function(x) x$spearman_raw)
)))

# heatmap of spearman correlations with fdr p-values
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
  correlations_df = correlations_df %>%
    filter(Definition %in% order_definitions) %>%
    arrange(factor(Definition, levels = order_definitions))
  
  # Create correlation matrix
  corr_matrix = matrix(correlations_df$Spearman_Raw, nrow = 1)
  rownames(corr_matrix) = "ENIGMA Bipolar"
  colnames(corr_matrix) = correlations_df$Definition
  
  # Create p-value matrix
  p_matrix = matrix(correlations_df$FDR_P, nrow = 1)
  rownames(p_matrix) = "ENIGMA Bipolar"
  colnames(p_matrix) = correlations_df$Definition
  
  # Create labels
  labels = matrix(paste0(round(corr_matrix, 2), "\n(", format(p_matrix, digits = 3), ")"),
                  nrow = nrow(corr_matrix), ncol = ncol(corr_matrix))
  
  # Determine color range
  if (is.null(max_abs_corr)) {
    max_abs_corr = max(abs(corr_matrix))
  }
  
  # Create a custom color palette
  custom_colors = colorRampPalette(c("#6666ff", "white", "#ff3333"))(100)
  
  # Create the heatmap
  pheatmap(corr_matrix, main = title, color = custom_colors,
           cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = labels, number_format = "%.2f",
           fontsize_number = 8, angle_col = 45, cellwidth = 80, cellheight = 35,
           fontcolor_number = "black", border_color = "black",
           breaks = seq(-max_abs_corr, max_abs_corr, length.out = 101))
}

# Create heatmaps
plot_correlation_heatmap(correlations_surface_area, "Spearman Correlations (Surface Area)", max_abs_corr)
plot_correlation_heatmap(correlations_thickness, "Spearman Correlations (Thickness)", max_abs_corr)

###### Scatter plots for top correlations ######

# --------------------------------------------------------
# Step 10. Scatter Plots for Top Correlations
# --------------------------------------------------------

# return the top n definitions by absolute spearman correlation
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
  # pool effect sizes across the top definitions
  all_x_values = unlist(lapply(top_corr, function(definition) data_list[[definition]]$cohen_d))
  all_y_values = unlist(lapply(top_corr, function(definition) data_list[[definition]]$enigma_d))
  # shared axis limits so panels are comparable
  x_range = range(all_x_values, na.rm = TRUE)
  y_range = range(all_y_values, na.rm = TRUE)
  list(x_range = x_range, y_range = y_range)
}

global_range_surface_area = calculate_global_range(results_list_surface_area_final, top_corr_surface_area)
global_range_thickness = calculate_global_range(results_list_thickness_final, top_corr_thickness)

# scatter plot of UKB cohen's d vs ENIGMA enigma_d for one definition
create_scatter_plot = function(data, definition, corr, p_value, measure, global_range, color) {
  # Remove rows with missing values
  data = data[complete.cases(data[, c("cohen_d", "enigma_d")]), ]

  # points plus a linear fit line
  ggplot(data, aes(x = cohen_d, y = enigma_d)) +
    geom_point(size = 2, alpha = 0.7, color = color) +
    geom_smooth(method = "lm", se = FALSE, color = color, linewidth = 1) +
    xlab("Effect Size (UKB)") +
    ylab("Effect Size (ENIGMA)") +
    # lock axes to the shared global range
    scale_y_continuous(breaks = seq(min(global_range$y_range), max(global_range$y_range), length.out = 5),
                       limits = global_range$y_range,
                       expand = c(0.01, 0)) +
    scale_x_continuous(breaks = seq(min(global_range$x_range), max(global_range$x_range), length.out = 5),
                       limits = global_range$x_range,
                       labels = function(x) sprintf("%.3f", x),
                       expand = c(0.01, 0)) +
    # match all plot elements to the panel color
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
    # annotate corner with r and fdr p-value
    annotate("text", x = Inf, y = -Inf,
             label = paste0("r = ", round(corr, 3), "\np = ", format.pval(p_value, digits = 3)),
             hjust = 1.1, vjust = -0.8, size = 3, color = color) +
    labs(title = definition)
}

# Scatter plots for surface area
# build one panel per top surface-area definition
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
# BD_Enigma_scatterplot_thickness_top3
scatter_plot_grid_surface_area = do.call(grid.arrange, c(scatter_plots_surface_area, ncol = 3))

# Scatter plots for thickness
# build one panel per top thickness definition
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