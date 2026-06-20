# --------------------------------------------------------
# ENIGMA - Calculating correlations with MDD working group summary statistics 
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
             'tidyr', 'pheatmap', 'ggrepel', 'gridExtra', 'here')

# Call the package_install function 
package_install(packages)

# Loading csv files
CT   = fread(here("data", "enigma_sumstats", "mdd", "CT.csv"))
SA   = fread(here("data", "enigma_sumstats", "mdd", "SA.csv"))

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
  # strip prefixes and check structures line up with measurement names
  structures_lower = tolower(structures)
  measurement_lower = tolower(gsub(paste0("^", tolower(measure_name), "_"), "", measurement_list))

  order_match = all(structures_lower == measurement_lower)
  cat("Order matches:", order_match, "\n")

  # flag any positions where the names diverge
  # Print mismatches if any
  if (!order_match) {
    mismatches = which(structures_lower != measurement_lower)
    cat("Mismatches at positions:", paste(mismatches, collapse = ", "), "\n")
    cat("Structure names at mismatches:", paste(structures[mismatches], collapse = ", "), "\n")
    cat("Measurement names at mismatches:", paste(measurement_list[mismatches], collapse = ", "), "\n")
  }
  
  # pair each structure with its d_icv effect size
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

###### Build UKBB-style abbreviated names (CT and SA) ######

create_comparison_dataframe_ct = function(ct_data) {
  # Extract Structure names from ct_data
  structure_names = ct_data$Structure

  # build a ukbb-style abbreviated name per structure
  # Create list of processed name abbreviated
  processed_name_abbreviated = sapply(structure_names, function(name) {
    # Extract the hemisphere (L or R) and the structure name
    hemisphere = substr(name, 1, 1)
    structure = substr(name, 3, nchar(name))

    # map hemisphere letter to lf/rf suffix
    # Convert hemisphere to lf or rf
    hemisphere_suffix = ifelse(hemisphere == "L", "_lf", "_rf")

    # Construct the processed name
    paste0("mean_thickness_", tolower(structure), hemisphere_suffix)
  })

  # bundle structure, abbreviated name, and effect size
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

  # build a ukbb-style abbreviated name per structure
  # Create list of processed name abbreviated
  processed_name_abbreviated = sapply(structure_names, function(name) {
    # Extract the hemisphere (L or R) and the structure name
    hemisphere = substr(name, 1, 1)
    structure = substr(name, 3, nchar(name))

    # map hemisphere letter to lf/rf suffix
    # Convert hemisphere to lf or rf
    hemisphere_suffix = ifelse(hemisphere == "L", "_lf", "_rf")

    # Construct the processed name
    paste0("area_", tolower(structure), hemisphere_suffix)
  })

  # bundle structure, abbreviated name, and effect size
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
cat("\nFirst few rows of the CT comparison dataframe:\n")
print(tail(comparison_df_CT, 10))

# Check if the number of rows match
cat("\nNumber of rows in the CTcomparison dataframe:", nrow(comparison_df_CT), "\n")
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


# Loading results lists (thickness & surface area)
load(here("data", "results_lists_thickness_SA_linear_regression.RData"))

###### Lookup table and feature-name extraction ######

# --------------------------------------------------------
# Extracting Feature Names
# --------------------------------------------------------

# Path to the neuroimaging lookup table
lookup_table_path = here("data", "lookup_table_neuroimaging_ukbb.xlsx")

# read the lookup table from excel
read_lookup_table = function(filepath) {
  read_excel(filepath)
}

# Read the lookup table
lookup_table_df = read_lookup_table(lookup_table_path)

# Rename the prefixes "baseline_" to "bl_" for each feature value in the table
lookup_table_df$corrected_id = gsub("baseline_", "bl_", lookup_table_df$corrected_id)

# merge results with lookup table, splitting matched from no_match
merge_with_lookup_table = function(results_df, lookup_df) {
  # coerce join keys to character so the merge lines up
  results_df$feature = as.character(results_df$feature)
  lookup_df$corrected_id = as.character(lookup_df$corrected_id)

  # Select only the necessary columns from the lookup table
  lookup_df = lookup_df[, c("corrected_id", "field", "processed_name_abbreviated")]

  # left-join results onto the lookup by feature id
  merged_df = merge(results_df, lookup_df, by.x = 'feature', by.y = 'corrected_id', all.x = TRUE)

  # Identifying 'no_match' rows
  no_match_rows = is.na(merged_df$processed_name_abbreviated)
  merged_df$processed_name_abbreviated[no_match_rows] = "no_match"

  # split matched from no_match for follow-up
  # Separating 'no_match' entries for further investigation
  no_match_df = merged_df[no_match_rows, ]
  matched_df = merged_df[!no_match_rows, ]
  
  return(list(matched = matched_df, no_match = no_match_df))
}

# --------------------------------------------------------
# Subsetting lookup tables - subcortical and thickness measures
# --------------------------------------------------------

# Subset for t1_subcortical_sMRI
subcortical_df = lookup_table_df %>%
  filter(field == "t1_thickness_area")

# Split into two tables based on prefix
lookup_table_area = subcortical_df %>%
  filter(grepl("^area_", processed_name_abbreviated))

lookup_table_thickness = subcortical_df %>%
  filter(grepl("^mean_thickness", processed_name_abbreviated))

###### Fuzzy-matching ENIGMA names to UKBB features ######

# --------------------------------------------------------
# Matching summary statistics names against lookup_table
# --------------------------------------------------------
# Step 1: Add 'enigma_name' column to lookup_table_area and lookup_table_thickness
lookup_table_area$enigma_name = NA
lookup_table_thickness$enigma_name = NA

# fuzzy-match ENIGMA region names to UKBB feature ids
match_enigma_names = function(enigma_name, lookup_table) {
  # Remove prefix and convert to lower case
  enigma_name_clean = tolower(gsub("^(L|R)_", "", enigma_name))
  
  # approximate-match against the abbreviated names
  # Find the best match in the lookup_table
  best_match = lookup_table$processed_name_abbreviated[agrep(enigma_name_clean, lookup_table$processed_name_abbreviated, max.distance = 0.2)]

  # keep only matches on the correct hemisphere
  # Filter out matches with incorrect hemisphere
  hemisphere = ifelse(grepl("^L_", enigma_name), "lf", "rf")
  best_match = best_match[grepl(hemisphere, best_match)]

  # special-case insula to avoid spurious fuzzy hits
  # If the enigma_name is "insula", filter the best_match to only include entries that contain "insula" in their name
  if (enigma_name_clean == "insula") {
    best_match = best_match[grepl("insula", best_match)]
  }

  # return first match (or NAs when nothing matched)
  if (length(best_match) > 0) {
    return(list(best_guess = enigma_name,
                processed_name_abbreviated = best_match[1],
                corrected_id = lookup_table$corrected_id[lookup_table$processed_name_abbreviated == best_match[1]]))
  } else {
    return(list(best_guess = NA, processed_name_abbreviated = NA, corrected_id = NA))  
  }
}

# --------------------------------------------------------
# Step 3: Apply the matching function to create separate objects with the best guess match
# --------------------------------------------------------

# Surface Area
best_guess_area = sapply(SA$Structure, function(x) match_enigma_names(x, lookup_table_area)$best_guess)
processed_name_abbreviated_area = sapply(SA$Structure, function(x) match_enigma_names(x, lookup_table_area)$processed_name_abbreviated)
corrected_id_area = sapply(SA$Structure, function(x) match_enigma_names(x, lookup_table_area)$corrected_id)

# Cortical Thickness
best_guess_thickness = sapply(CT$Structure, function(x) match_enigma_names(x, lookup_table_thickness)$best_guess)
processed_name_abbreviated_thickness = sapply(CT$Structure, function(x) match_enigma_names(x, lookup_table_thickness)$processed_name_abbreviated)
corrected_id_thickness = sapply(CT$Structure, function(x) match_enigma_names(x, lookup_table_thickness)$corrected_id)

# --------------------------------------------------------
# Step 4: Add relevant columns from lookup tables to SA and CT
# --------------------------------------------------------
SA$best_guess_area = best_guess_area
SA$processed_name_abbreviated = processed_name_abbreviated_area
SA$corrected_id = corrected_id_area

CT$matched_thickness = best_guess_thickness
CT$processed_name_abbreviated = processed_name_abbreviated_thickness
CT$corrected_id = corrected_id_thickness

# --------------------------------------------------------
# Step 5: Remove NA from SA and CT
# --------------------------------------------------------
SA = SA[complete.cases(SA$corrected_id), ]
CT = CT[complete.cases(CT$corrected_id), ]

# optionally save matched SA and CT name tables

###### Merge ENIGMA effects into results lists ######

# --------------------------------------------------------
# Step 6: Merge SA and CT with the results_list_surface_area and results_list_thickness
# --------------------------------------------------------
results_list_surface_area_merged = lapply(results_list_surface_area, function(df) {
  merged_df = merge(df, SA[, c("corrected_id", "d_icv", "se_icv", "processed_name_abbreviated")], by.x = "feature", by.y = "corrected_id", all.x = TRUE)
  merged_df
})

results_list_thickness_merged = lapply(results_list_thickness, function(df) {
  merged_df = merge(df, CT[, c("corrected_id", "d_icv", "se_icv", "processed_name_abbreviated")], by.x = "feature", by.y = "corrected_id", all.x = TRUE)
  merged_df
})

# --------------------------------------------------------
# Removing NA rows 
# --------------------------------------------------------
remove_na_rows = function(df, column_name) {
  # Check if the column exists in the dataframe
  if(!column_name %in% names(df)) {
    stop("Column not found in the dataframe")
  }
  # Remove rows where the column has NA
  df[!is.na(df[[column_name]]), ]
}

# Apply the function to each dataframe in the list
results_list_surface_area_merged = lapply(results_list_surface_area_merged, remove_na_rows, "processed_name_abbreviated")
results_list_thickness_merged    = lapply(results_list_thickness_merged, remove_na_rows, "processed_name_abbreviated")

# optionally save merged surface area and thickness results lists

# --------------------------------------------------------
# Step 7: Add cohen_d and se_d columns to the results_list objects
# --------------------------------------------------------
# Cohen's d from the linear regressions

# Add cohen_d column to surface area results
results_list_surface_area_final = lapply(results_list_surface_area_merged, function(df) {
  # use the Cohen's d
  df$cohen_d = df$cohen_d_corrected
  # Return the modified data frame
  df
})

# Add cohen_d column to thickness results
results_list_thickness_final = lapply(results_list_thickness_merged, function(df) {
  # use the Cohen's d
  df$cohen_d = df$cohen_d_corrected
  # Return the modified data frame
  df
})

# Removing MDD_noimpairment and MDD2_w_impair phenotypes 
results_list_surface_area_final$MDD_noimpairment = NULL
results_list_surface_area_final$MDD2_w_impair = NULL
results_list_thickness_final$MDD_noimpairment = NULL
results_list_thickness_final$MDD2_w_impair = NULL

# Extract significant cortical thickness regions and store as separate objects 
for (name in names(results_list_thickness_final)) {
  # Extract the dataframe
  df = results_list_thickness_final[[name]]
  
  # Filter rows where significant_fdr is TRUE
  significant_df = df[df$significant_fdr == TRUE, ]
  
  # Assign the filtered dataframe to a new object in the environment
  assign(paste0(name, "_cortical_thickness_significant"), significant_df)
}


# Extract significant surface area regions and store as separate objects 
for (name in names(results_list_surface_area_final)) {
  # Extract the dataframe
  df = results_list_surface_area_final[[name]]
  
  # Filter rows where significant_fdr is TRUE
  significant_df = df[df$significant_fdr == TRUE, ]
  
  # Assign the filtered dataframe to a new object in the environment
  assign(paste0(name, "_surface_area_significant"), significant_df)
}

###### Top-5 feature extraction ######

# --------------------------------------------------------
# Step 7: Extracting the top 5 features based on the absolute value of d_icv
# This is done for each depression definition
# --------------------------------------------------------

# Function to get top 5 features by absolute value of d_icv
get_top_5 = function(df) {
  df %>%
    mutate(abs_d_icv = abs(d_icv)) %>%
    top_n(5, abs_d_icv) %>%
    select(-abs_d_icv)
}

# Apply get_top_5 function to each depression definition for cortical thickness
top_5_list_thickness = lapply(results_list_thickness_final, get_top_5)

# Apply get_top_5 function to each depression definition for surface area
top_5_list_surface_area = lapply(results_list_surface_area_final, get_top_5)

# Combine top 5 lists into a single data frame for cortical thickness
top_5_combined_thickness = bind_rows(
  lapply(names(top_5_list_thickness), function(name) {
    top_5_list_thickness[[name]] %>%
      mutate(definition = name)
  })
)

# Cohen's d confidence interval
top_5_combined_thickness$d_lower_CI = top_5_combined_thickness$cohen_d_lower
top_5_combined_thickness$d_upper_CI = top_5_combined_thickness$cohen_d_upper
top_5_combined_thickness$cohen_d_se = top_5_combined_thickness$cohen_d_se

# Combine top 5 lists into a single data frame for surface area
top_5_combined_surface_area = bind_rows(
  lapply(names(top_5_list_surface_area), function(name) {
    top_5_list_surface_area[[name]] %>%
      mutate(definition = name)
  })
)

# Cohen's d confidence interval
top_5_combined_surface_area$d_lower_CI = top_5_combined_surface_area$cohen_d_lower
top_5_combined_surface_area$d_upper_CI = top_5_combined_surface_area$cohen_d_upper
top_5_combined_surface_area$cohen_d_se = top_5_combined_surface_area$cohen_d_se


# Extract top 5 effect sizes for surface area
top_5_surface_area = top_5_combined_surface_area %>%
  group_by(definition) %>%
  arrange(desc(abs(d_icv))) %>%
  slice_head(n = 5) %>%
  select(definition, processed_name_abbreviated, d_icv, cohen_d)

# Extract top 5 effect sizes for thickness
top_5_thickness = top_5_combined_thickness %>%
  group_by(definition) %>%
  arrange(desc(abs(d_icv))) %>%
  slice_head(n = 5) %>%
  select(definition, processed_name_abbreviated, d_icv, cohen_d)


###### Forest plots (ENIGMA-ordered and UKBB-ordered) ######

# Order depression definitions
definition_order = c("GPNoDep", "GPpsy", "Psypsy", "SelfRepDep", "DepAll", "ICD10Dep", "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur")

create_forest_plot = function(data, title) {
  # reshape to long: one row per cohort effect size
  data_long = data %>%
    select(processed_name_abbreviated, cohen_d, cohen_d_se, d_icv, se_icv, definition) %>%
    pivot_longer(cols = c(cohen_d, d_icv),
                 names_to = "Group",
                 values_to = "Effect_Size") %>%
    # label cohorts and compute 95% CIs from each SE
    mutate(Group = case_when(
      Group == "cohen_d" ~ "UK Biobank",
      Group == "d_icv"   ~ "Enigma MDD",
      TRUE               ~ Group
    ),
    lower_CI = case_when(
      Group == "UK Biobank" ~ Effect_Size - 1.96 * cohen_d_se,
      Group == "Enigma MDD" ~ Effect_Size - 1.96 * se_icv
    ),
    upper_CI = case_when(
      Group == "UK Biobank" ~ Effect_Size + 1.96 * cohen_d_se,
      Group == "Enigma MDD" ~ Effect_Size + 1.96 * se_icv
    ))

  # order regions by ENIGMA effect, definitions by predefined order
  data_long$processed_name_abbreviated = factor(data_long$processed_name_abbreviated, levels = unique(data$processed_name_abbreviated[order(data$d_icv)]))
  data_long$definition = factor(data_long$definition, levels = definition_order)

  # full-panel background rectangle per facet
  rect_data = data.frame(
    definition = unique(data_long$definition),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf
  )

  # draw point-ranges faceted by definition
  ggplot(data_long, aes(x = Effect_Size, y = processed_name_abbreviated, color = Group, shape = Group)) +
    geom_rect(
      data = rect_data,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey90",
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    geom_pointrange(aes(xmin = lower_CI, xmax = upper_CI), position = position_dodge(width = 0.5), linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_grid(cols = vars(definition), scales = "free_y", space = "free") +
    labs(title = title, x = "Effect Size (95% CI)", y = "Processed Name Abbreviated") +
    scale_color_manual(values = c("UK Biobank" = "blue", "Enigma MDD" = "red")) +
    scale_shape_manual(values = c("UK Biobank" = 16, "Enigma MDD" = 17)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey90", colour = NA),
      panel.spacing = unit(0.5, "lines")
    )
}

# Create forest plots
forest_plot_thickness = create_forest_plot(
  top_5_combined_thickness,
  "Top 5 Associations between Depression Phenotypes and Cortical Thickness"
)
forest_plot_surface_area = create_forest_plot(
  top_5_combined_surface_area,
  "Top 5 Associations between Depression Phenotypes and Surface Area"
)

print(forest_plot_thickness)
print(forest_plot_surface_area)

# --------------------------------------------------------
# Extracting the top 5 regions according to the UK Biobank Sample 
# --------------------------------------------------------
# Extract top 5 effect sizes for surface area based on cohen_d
top_5_surface_area_cohen_d = top_5_combined_surface_area %>%
  group_by(definition) %>%
  arrange(desc(abs(cohen_d))) %>%
  slice_head(n = 5) %>%
  select(definition, processed_name_abbreviated, d_icv, se_icv, cohen_d, cohen_d_se)

# Extract top 5 effect sizes for thickness based on cohen_d
top_5_thickness_cohen_d = top_5_combined_thickness %>%
  group_by(definition) %>%
  arrange(desc(abs(cohen_d))) %>%
  slice_head(n = 5) %>%
  select(definition, processed_name_abbreviated, d_icv, se_icv, cohen_d, cohen_d_se)

# Function to get top 5 features by absolute value of cohen_d across all definitions
get_top_5_overall = function(df) {
  df %>%
    mutate(abs_cohen_d = abs(cohen_d)) %>%
    group_by(processed_name_abbreviated) %>%
    summarise(max_abs_cohen_d = max(abs_cohen_d)) %>%
    top_n(5, max_abs_cohen_d) %>%
    select(processed_name_abbreviated)
}

# Apply get_top_5_overall function to surface area and thickness data
top_5_surface_area_overall = get_top_5_overall(top_5_combined_surface_area)
top_5_thickness_overall = get_top_5_overall(top_5_combined_thickness)

ukbb_forest_plot = function(data, title) {
  # Reshape the data to have separate rows for cohen_d and d_icv
  data_long = data %>%
    select(processed_name_abbreviated, cohen_d, cohen_d_se, d_icv, se_icv, definition) %>%
    pivot_longer(cols = c(cohen_d, d_icv),
                 names_to = "Group",
                 values_to = "Effect_Size") %>%
    # label cohorts and compute 95% CIs from each SE
    mutate(Group = case_when(
      Group == "cohen_d" ~ "UK Biobank",
      Group == "d_icv"   ~ "Enigma MDD",
      TRUE               ~ Group
    ),
    lower_CI = case_when(
      Group == "UK Biobank" ~ Effect_Size - 1.96 * cohen_d_se,
      Group == "Enigma MDD" ~ Effect_Size - 1.96 * se_icv
    ),
    upper_CI = case_when(
      Group == "UK Biobank" ~ Effect_Size + 1.96 * cohen_d_se,
      Group == "Enigma MDD" ~ Effect_Size + 1.96 * se_icv
    ))

  # here features are ordered by the UKB cohen_d instead
  # Order the features based on cohen_d effect size
  data_long$processed_name_abbreviated = factor(data_long$processed_name_abbreviated, levels = unique(data$processed_name_abbreviated[order(data$cohen_d)]))
  
  # Order the depression definitions
  data_long$definition = factor(data_long$definition, levels = definition_order)
  
  # Create the plot
  ggplot(data_long, aes(x = Effect_Size, y = processed_name_abbreviated, color = Group, shape = Group)) +
    geom_pointrange(aes(xmin = lower_CI, xmax = upper_CI), position = position_dodge(width = 0.5), linewidth = 0.5) +
    scale_color_manual(values = c("UK Biobank" = "blue", "Enigma MDD" = "red")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    facet_grid(cols = vars(definition), scales = "free_y", space = "free") +
    labs(title = title, x = "Effect Size (95% CI)", y = "Processed Name Abbreviated") +
    scale_shape_manual(values = c("UK Biobank" = 16, "Enigma MDD" = 17)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold"),
      panel.spacing = unit(0.5, "lines")
    )
}

# Create forest plots
ukbb_forest_plot_thickness = ukbb_forest_plot(
  top_5_thickness_cohen_d,
  "Top 5 Associations (UK Biobank) between Depression Phenotypes and Cortical Thickness"
)
ukbb_forest_plot_surface_area = ukbb_forest_plot(
  top_5_surface_area_cohen_d,
  "Top 5 Associations (UK Biobank) between Depression Phenotypes and Surface Area"
)

print(ukbb_forest_plot_thickness)
print(ukbb_forest_plot_surface_area)

###### Correlations between UKBB and ENIGMA effects ######

# --------------------------------------------------------
# Step 8: Calculate the correlations between effect sizes
# --------------------------------------------------------
# spearman correlation and p-value between two columns, NA on error
calculate_correlations = function(df, col1, col2) {
  # pull the two effect-size columns into x/y
  cor_data = data.frame(
    x = df[[col1]],
    y = df[[col2]]
  )

  # spearman rho, NA if it fails
  spearman_raw = tryCatch(
    cor(cor_data$x, cor_data$y, method = "spearman"),
    error = function(e) NA
  )

  # matching p-value, NA if it fails
  spearman_p = tryCatch(
    cor.test(cor_data$x, cor_data$y, method = "spearman", exact = FALSE)$p.value,
    error = function(e) NA
  )
  
  #fdr_p = p.adjust(spearman_p, method = "fdr")
  
  #list(spearman_raw = spearman_raw, spearman_p = spearman_p, fdr_p = fdr_p)
  list(spearman_raw = spearman_raw, spearman_p = spearman_p)
}

# Calculating correlations for surface area measures across depression phenotypes 
correlations_surface_area_GPNoDep            = calculate_correlations(results_list_surface_area_final$GPNoDep, "cohen_d", "d_icv")
correlations_surface_area_Psypsy             = calculate_correlations(results_list_surface_area_final$Psypsy, "cohen_d", "d_icv")
correlations_surface_area_SelfRepDep         = calculate_correlations(results_list_surface_area_final$SelfRepDep, "cohen_d", "d_icv")
correlations_surface_area_GPpsy              = calculate_correlations(results_list_surface_area_final$GPpsy, "cohen_d", "d_icv")
correlations_surface_area_DepAll             = calculate_correlations(results_list_surface_area_final$DepAll, "cohen_d", "d_icv")
correlations_surface_area_ICD10Dep           = calculate_correlations(results_list_surface_area_final$ICD10Dep, "cohen_d", "d_icv")
correlations_surface_area_ICD10Dep_exclpsych = calculate_correlations(results_list_surface_area_final$ICD10Dep_exclpsych, "cohen_d", "d_icv")
correlations_surface_area_LifetimeMDD        = calculate_correlations(results_list_surface_area_final$LifetimeMDD, "cohen_d", "d_icv")
#correlations_surface_area_MDD_noimpairment   = calculate_correlations(results_list_surface_area_final$MDD_noimpairment, "cohen_d", "d_icv")
#correlations_surface_area_MDD2_w_impair      = calculate_correlations(results_list_surface_area_final$MDD2_w_impair, "cohen_d", "d_icv")
correlations_surface_area_MDDRecur           = calculate_correlations(results_list_surface_area_final$MDDRecur, "cohen_d", "d_icv")


# Calculating correlations for surface area measures across depression phenotypes 
correlations_thickness_GPNoDep            = calculate_correlations(results_list_thickness_final$GPNoDep, "cohen_d", "d_icv")
correlations_thickness_Psypsy             = calculate_correlations(results_list_thickness_final$Psypsy, "cohen_d", "d_icv")
correlations_thickness_SelfRepDep         = calculate_correlations(results_list_thickness_final$SelfRepDep, "cohen_d", "d_icv")
correlations_thickness_GPpsy              = calculate_correlations(results_list_thickness_final$GPpsy, "cohen_d", "d_icv")
correlations_thickness_DepAll             = calculate_correlations(results_list_thickness_final$DepAll, "cohen_d", "d_icv")
correlations_thickness_ICD10Dep           = calculate_correlations(results_list_thickness_final$ICD10Dep, "cohen_d", "d_icv")
correlations_thickness_ICD10Dep_exclpsych = calculate_correlations(results_list_thickness_final$ICD10Dep_exclpsych, "cohen_d", "d_icv")
correlations_thickness_LifetimeMDD        = calculate_correlations(results_list_thickness_final$LifetimeMDD, "cohen_d", "d_icv")
#correlations_thickness_MDD_noimpairment   = calculate_correlations(results_list_thickness_final$MDD_noimpairment, "cohen_d", "d_icv")
#correlations_thickness_MDD2_w_impair      = calculate_correlations(results_list_thickness_final$MDD2_w_impair, "cohen_d", "d_icv")
correlations_thickness_MDDRecur           = calculate_correlations(results_list_thickness_final$MDDRecur, "cohen_d", "d_icv")

# --------------------------------------------------------
# Organize Correlation Results
# --------------------------------------------------------

# Creating a list of dataframes for correlations surface area
correlations_surface_area = list(
  GPNoDep = correlations_surface_area_GPNoDep,
  Psypsy = correlations_surface_area_Psypsy,
  SelfRepDep = correlations_surface_area_SelfRepDep,
  GPpsy = correlations_surface_area_GPpsy,
  DepAll = correlations_surface_area_DepAll,
  ICD10Dep = correlations_surface_area_ICD10Dep,
  ICD10Dep_exclpsych = correlations_surface_area_ICD10Dep_exclpsych,
  LifetimeMDD = correlations_surface_area_LifetimeMDD,
  #MDD_noimpairment = correlations_surface_area_MDD_noimpairment,
  #MDD2_w_impair = correlations_surface_area_MDD2_w_impair,
  MDDRecur = correlations_surface_area_MDDRecur
)

# Creating a list of dataframes for correlations thickness
correlations_thickness = list(
  GPNoDep = correlations_thickness_GPNoDep,
  Psypsy = correlations_thickness_Psypsy,
  SelfRepDep = correlations_thickness_SelfRepDep,
  GPpsy = correlations_thickness_GPpsy,
  DepAll = correlations_thickness_DepAll,
  ICD10Dep = correlations_thickness_ICD10Dep,
  ICD10Dep_exclpsych = correlations_thickness_ICD10Dep_exclpsych,
  LifetimeMDD = correlations_thickness_LifetimeMDD,
  #MDD_noimpairment = correlations_thickness_MDD_noimpairment,
  #MDD2_w_impair = correlations_thickness_MDD2_w_impair,
  MDDRecur = correlations_thickness_MDDRecur
)

# --------------------------------------------------------
# FDR Correction for Thickness and Surface Area effect sizes 
# --------------------------------------------------------

# Extract raw p-values for surface area correlations
surface_area_pvals = unlist(lapply(correlations_surface_area, function(x) x$spearman_p))

# Extract raw p-values for cortical thickness correlations
thickness_pvals = unlist(lapply(correlations_thickness, function(x) x$spearman_p))

# FDR correction for surface area p-values
surface_area_fdr_pvals = p.adjust(surface_area_pvals, method = "fdr")

# FDR correction for cortical thickness p-values
thickness_fdr_pvals = p.adjust(thickness_pvals, method = "fdr")

# Create new lists with FDR-corrected p-values for surface area
correlations_surface_area_fdr = lapply(names(correlations_surface_area), function(name) {
  list(
    spearman_raw = correlations_surface_area[[name]]$spearman_raw,
    spearman_p = correlations_surface_area[[name]]$spearman_p,
    fdr_p = surface_area_fdr_pvals[name]
  )
})

# Add FDR-corrected p-values to the respective data frames
for (i in seq_along(correlations_surface_area)) {
  correlations_surface_area[[i]]$fdr_p = surface_area_fdr_pvals[i]
}

for (i in seq_along(correlations_thickness)) {
  correlations_thickness[[i]]$fdr_p = thickness_fdr_pvals[i]
}

# --------------------------------------------------------
# Write the cortical Spearman correlations to CSV
# (per-definition, for surface area and cortical thickness)
# --------------------------------------------------------
cortical_correlations_corrected = bind_rows(
  bind_rows(lapply(names(correlations_surface_area), function(name) {
    data.frame(
      Definition   = name,
      Measure      = "Surface Area",
      Spearman_r   = correlations_surface_area[[name]]$spearman_raw,
      p_value      = correlations_surface_area[[name]]$spearman_p,
      FDR_p_value  = correlations_surface_area[[name]]$fdr_p,
      stringsAsFactors = FALSE
    )
  })),
  bind_rows(lapply(names(correlations_thickness), function(name) {
    data.frame(
      Definition   = name,
      Measure      = "Cortical Thickness",
      Spearman_r   = correlations_thickness[[name]]$spearman_raw,
      p_value      = correlations_thickness[[name]]$spearman_p,
      FDR_p_value  = correlations_thickness[[name]]$fdr_p,
      stringsAsFactors = FALSE
    )
  }))
)

write.csv(
  cortical_correlations_corrected,
  file = here("output", "mdd_sact_cortical_enigma_correlations.csv"),
  row.names = FALSE
)

# Calculating the maximum absolute correlation across both datasets
max_abs_corr = max(abs(c(
  sapply(correlations_surface_area, function(x) x$spearman_raw),
  sapply(correlations_thickness, function(x) x$spearman_raw)
)))

###### Correlation heatmap ######

# --------------------------------------------------------
# Step 9. Visualize Correlations
# --------------------------------------------------------

# heatmap of spearman correlations with p-value labels
plot_correlation_heatmap = function(correlations, title, max_abs_corr = NULL) {
  # Order definitions by phenotype granularity
  order_definitions = c("GPNoDep", "GPpsy", "Psypsy", "SelfRepDep", "DepAll", "ICD10Dep", "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur")

  # flatten the correlation list into a data frame
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
  
  # restrict to known definitions and apply the ordering
  # Reorder the correlations based on the order of definitions
  correlations_df = correlations_df %>%
    filter(Definition %in% order_definitions) %>%
    arrange(factor(Definition, levels = order_definitions))

  # single-row matrices of rho and FDR p
  # Create correlation matrix
  corr_matrix = matrix(correlations_df$Spearman_Raw, nrow = 1)
  rownames(corr_matrix) = "ENIGMA MDD"
  colnames(corr_matrix) = correlations_df$Definition

  # Create p-value matrix
  p_matrix = matrix(correlations_df$FDR_P, nrow = 1)
  rownames(p_matrix) = "ENIGMA MDD"
  colnames(p_matrix) = correlations_df$Definition

  # cell labels showing rho with p underneath
  # Create labels
  labels = matrix(paste0(round(corr_matrix, 2), "\n(", format(p_matrix, digits = 3), ")"),
                  nrow = nrow(corr_matrix), ncol = ncol(corr_matrix))

  # symmetric color range centered on zero
  # Determine color range
  if (is.null(max_abs_corr)) {
    max_abs_corr = max(abs(corr_matrix))
  }

  # Create a custom color palette
  custom_colors = colorRampPalette(c("#6666ff", "white", "#ff3333"))(100)

  # render the single-row heatmap
  # Create the heatmap
  pheatmap(corr_matrix, main = title, color = custom_colors,
           cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = labels, number_format = "%.2f",
           fontsize_number = 8, angle_col = 45, cellwidth = 80, cellheight = 35,
           fontcolor_number = "black", border_color = "black",
           breaks = seq(-max_abs_corr, max_abs_corr, length.out = 101))
}

plot_correlation_heatmap(correlations_surface_area, "Spearman Correlations (Surface Area)")
plot_correlation_heatmap(correlations_thickness, "Spearman Correlations (Thickness)")

###### Scatter plots for top correlations ######

# --------------------------------------------------------
# Step 10. Scatter Plots for Top Correlations
# --------------------------------------------------------

# get names of the top n strongest absolute correlations
get_top_correlations = function(correlations, n = 3) {
  corr_values = sapply(correlations, function(x) x$spearman_raw)
  top_indices = order(abs(corr_values), decreasing = TRUE)[1:n]
  names(correlations)[top_indices]
}

# Get the top 3 strongest correlations for surface area and thickness
top_corr_surface_area = get_top_correlations(correlations_surface_area, n = 3)
top_corr_thickness = get_top_correlations(correlations_thickness, n = 3)

# Define colors
plot_colors = c("#4a80ed", "#ed984a", "#45a049")

# Calculate global ranges for surface area and thickness
calculate_global_range = function(data_list, top_corr) {
  all_x_values = unlist(lapply(top_corr, function(definition) data_list[[definition]]$cohen_d))
  all_y_values = unlist(lapply(top_corr, function(definition) data_list[[definition]]$d_icv))
  x_range = range(all_x_values, na.rm = TRUE)
  y_range = range(all_y_values, na.rm = TRUE)
  list(x_range = x_range, y_range = y_range)
}

global_range_surface_area = calculate_global_range(results_list_surface_area_final, top_corr_surface_area)
global_range_thickness = calculate_global_range(results_list_thickness_final, top_corr_thickness)

# scatter plot of UKB cohen_d against ENIGMA d_icv for one definition
create_scatter_plot = function(data, definition, corr, p_value, measure, global_range, color) {
  # points plus a linear fit line
  ggplot(data[[definition]], aes(x = cohen_d, y = d_icv)) +
    geom_point(size = 2, alpha = 0.7, color = color) +
    geom_smooth(method = "lm", se = FALSE, color = color, linewidth = 1) +
    xlab("Cohen's d (UKB)") +
    ylab("Cohen's d (ENIGMA)") +
    # fix axes to the shared global range
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
    # annotate corner with r and p
    annotate("text", x = Inf, y = -Inf,
             label = paste0("r = ", round(corr, 3), "\np = ", format.pval(p_value, digits = 3)),
             hjust = 1.1, vjust = -0.8, size = 3, color = color) +
    labs(title = definition)
}

# Define colors
plot_colors = c("#4a80ed", "#ed984a", "#45a049")

# Scatter plots for surface area
scatter_plots_surface_area = lapply(seq_along(top_corr_surface_area), function(i) {
  definition = top_corr_surface_area[i]
  corr = correlations_surface_area[[definition]]$spearman_raw
  p_value = correlations_surface_area[[definition]]$spearman_p
  print(paste("Plot", i, "corresponds to:", paste0("Surface Area: ", definition)))
  create_scatter_plot(results_list_surface_area_final, definition, corr, p_value, 
                      "Surface Area", global_range_surface_area, plot_colors[i])
})

# Arrange the surface area plots in a grid layout
# 1107 x 417
# CT_EnigmaMDD_ukbb_scatterplot
scatter_plot_grid_surface_area = do.call(grid.arrange, c(scatter_plots_surface_area, ncol = 3))

# Scatter plots for thickness
scatter_plots_thickness = lapply(seq_along(top_corr_thickness), function(i) {
  definition = top_corr_thickness[i]
  corr = correlations_thickness[[definition]]$spearman_raw
  p_value = correlations_thickness[[definition]]$spearman_p
  print(paste("Plot", i, "corresponds to:", paste0("Thickness: ", definition)))
  create_scatter_plot(results_list_thickness_final, definition, corr, p_value, 
                      "Thickness", global_range_thickness, plot_colors[i])
})

# Arrange the thickness plots in a grid layout
scatter_plot_grid_thickness = do.call(grid.arrange, c(scatter_plots_thickness, ncol = 3))
