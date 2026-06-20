# --------------------------------------------------------
# ENIGMA - Calculating correlations with MDD working group white matter summary statistics
# ENIGMA MDD white matter (van Velzen et al. 2020)
# --------------------------------------------------------

###### Package setup ######

# Function to install and load packages
package_install = function(package_list) {
  for (pkg in package_list) {
    # install only if not already available
    if (!requireNamespace(pkg, quietly = TRUE)) {
      tryCatch({
        install.packages(pkg)
      }, error = function(e) {
        # warn but keep going on install failure
        message(paste0("Could not install ", pkg, ": ", e))
      })
    }
    # Load the library
    library(pkg, character.only = TRUE)
  }
}

# List of packages
packages = c('dplyr', 'data.table', 'readr', 'stringr', 'ggplot2', 'readxl',
             'tidyr', 'pheatmap', 'ggrepel', 'gridExtra', 'ggpubr', 'here')

# Call the package_install function 
package_install(packages)

###### Read ENIGMA summary statistics ######

# Specify the path to the xlsx file
xlsx_file = here("data", "enigma_sumstats", "mdd", "white_matter_sumstats.xlsx")

#-----------------
# Fractional Ansiotropy 
#-----------------
# Read in the first tab of the xlsx file as FA
FA = as.data.frame(read_excel(xlsx_file, sheet = 1))

# Renaming Cohens_D column
FA = FA %>% rename(Cohens_D = `Cohen’s d`)
FA = FA %>% rename(FDR_P_Value = `FDR P-value`)

# Removing FA regions not present in UKBB
FA = FA %>%
  filter(!Region %in% c("CC", "CR", "IC", "IFO", "AverageMD", "AverageFA"))

#-----------------
# Mean Diffusivity 
#-----------------
# Read in the second tab of the xlsx file as MD
MD = as.data.frame(read_excel(xlsx_file, sheet = 2))

MD = MD %>% rename(Cohens_D = `Cohen’s d`)

MD = MD %>% rename(FDR_P_Value = `FDR P-value`)

# Removing MD regions not present in UKBB
MD = MD %>%
  filter(!Region %in% c("CC", "CR", "IC", "IFO", "AverageMD", "AverageFA"))

###### Lookup table and effect-size helpers ######

# Creating an object to change column names and create the lower_CI and upper_CI columns for subsequent plotting
create_base_object = function(data, top_5_regions) {
  data %>%
    # standardise column names
    rename(Full.tract.name = `Full tract name`,
           Cohens_D_SE = SE) %>%
    # derive 95% CI bounds from SE
    mutate(
      Cohens_D_Lower_CI = Cohens_D - 1.96 * Cohens_D_SE,
      Cohens_D_Upper_CI = Cohens_D + 1.96 * Cohens_D_SE
    ) %>%
    # keep plotting columns and restrict to the top 5 regions
    select(Region, Full.tract.name, Cohens_D, Cohens_D_SE, Cohens_D_Lower_CI, Cohens_D_Upper_CI) %>%
    filter(Region %in% top_5_regions)
}

# build base objects after top_5_* are defined

# Loading results lists (FA & MD)
# Cohen's d from the linear regressions
load(here("data", "results_lists_FA_MD_linear_regression.RData"))

# --------------------------------------------------------
# Extracting Feature Names
# --------------------------------------------------------

# path to the neuroimaging lookup table
lookup_table_path = here("data", "lookup_table_neuroimaging_ukbb.xlsx")

# read the lookup table from xlsx
read_lookup_table = function(filepath) {
  read_excel(filepath)
}

# read the lookup table into a dataframe
lookup_table_df = read_lookup_table(lookup_table_path)

# Rename the prefixes "baseline_" to "bl_" for each feature value in the table
lookup_table_df$corrected_id = gsub("baseline_", "bl_", lookup_table_df$corrected_id)

# use Cohen's d from the results lists

# average hemispheric effects and extract FDR p-values
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

# Specify the path to the matching table Excel file
matching_table_path = here("data", "enigma_sumstats", "mdd", "matching_table_UKBB_to_Schmaal_FA_MD.xlsx")

# Read in the matching table for FA
matching_table_FA = read_excel(matching_table_path, sheet = 1)

# Read in the matching table for MD
matching_table_MD = read_excel(matching_table_path, sheet = 2)

# Apply the averaging and conversion function to each results list
averaged_results_lists_FA = lapply(results_list_FA, function(x) average_hemispheric_effects(x, matching_table_FA))
averaged_results_lists_MD = lapply(results_list_MD, function(x) average_hemispheric_effects(x, matching_table_MD))

# Combine the averaged results into a single data frame for FA and MD
averaged_results_df_FA = do.call(rbind, averaged_results_lists_FA)
averaged_results_df_MD = do.call(rbind, averaged_results_lists_MD)

# Adding rownames as a column - depression_definitions 
averaged_results_df_FA$Phenotypes = rownames(averaged_results_df_FA) 
averaged_results_df_MD$Phenotypes = rownames(averaged_results_df_MD) 

# Function to split the averaged_results_df_FA dataframe into separate dataframes for each depression definition
split_dataframe = function(averaged_results_df) {
  split(averaged_results_df, gsub("\\.[0-9]+$", "", averaged_results_df$Phenotypes))
}

# Split the averaged_results_df_FA dataframe
FA_split = split_dataframe(averaged_results_df_FA)
MD_split = split_dataframe(averaged_results_df_MD)

# Calculate standard error for MD data
FA_split = lapply(FA_split, function(df) {
  df$Cohens_D_SE = (df$Cohens_D_Upper_CI - df$Cohens_D_Lower_CI) / (2 * 1.96)
  return(df)
})

MD_split = lapply(MD_split, function(df) {
  df$Cohens_D_SE = (df$Cohens_D_Upper_CI - df$Cohens_D_Lower_CI) / (2 * 1.96)
  return(df)
})

# Function to get top 5 features by absolute value of Cohens_D from the meta-analytic results
get_top_5_meta = function(df) {
  df %>%
    arrange(desc(abs(Cohens_D))) %>%
    slice_head(n = 5) %>%
    rename(d_enigma = Cohens_D, se_enigma = SE)
}

# Apply get_top_5_meta function to the FA and MD data
top_5_FA_meta = get_top_5_meta(FA)
top_5_MD_meta = get_top_5_meta(MD)

# Function to rank top effect sizes across the nine phenotypes in the ukbb
rank_regions = function(df) {
  df %>%
    mutate(rank = rank(-abs(Cohens_D), ties.method = "first")) %>%
    mutate(points = max(rank) - rank + 1) %>%
    select(Region, Full.tract.name, Cohens_D, Cohens_D_SE, Cohens_D_Lower_CI, Cohens_D_Upper_CI, points)
}

get_top_5_ukbb = function(split_obj, exclude_regions = character(0)) {
  # drop excluded regions before ranking, default empty
  split_obj = lapply(split_obj, function(df) df[!df$Region %in% exclude_regions, ])

  # Apply ranking to each phenotype
  ranked_list = lapply(split_obj, rank_regions)
  
  # Combine all rankings
  all_rankings = do.call(rbind, ranked_list)
  
  # Sum points for each region
  top_regions = all_rankings %>%
    group_by(Region, Full.tract.name) %>%
    summarize(
      total_points = sum(points),
      mean_cohens_d = mean(Cohens_D),
      mean_se = mean(Cohens_D_SE),
      mean_lower_ci = mean(Cohens_D_Lower_CI),
      mean_upper_ci = mean(Cohens_D_Upper_CI)
    ) %>%
    ungroup() %>%
    arrange(desc(total_points)) %>%
    slice_head(n = 5)
  
  # Create a list similar to the meta objects
  top_5_list = list(
    Region = top_regions$Region,
    `Full tract name` = top_regions$Full.tract.name,
    d_ukbb = top_regions$mean_cohens_d,
    se_ukbb = top_regions$mean_se,
    lower_CI = top_regions$mean_lower_ci,
    upper_CI = top_regions$mean_upper_ci
  )
  
  return(top_5_list)
}

# apply to FA and MD split objects; exclude composite ROIs for FA
top_5_FA_ukbb = get_top_5_ukbb(FA_split, exclude_regions = c("CC", "CR", "IC", "IFO"))
top_5_MD_ukbb = get_top_5_ukbb(MD_split)

# create base objects for MD and FA
MD_base = create_base_object(MD, top_5_MD_ukbb$Region)
FA_base = create_base_object(FA, top_5_FA_ukbb$Region)

# Modify the extract_corresponding_regions function to include full tract names
extract_corresponding_regions = function(split_obj, top_5, is_ukbb = FALSE) {
  lapply(split_obj, function(df) {
    df_filtered = df[df$Region %in% top_5$Region, ]
    #df_filtered = add_full_tract_names(df_filtered)
    if (is_ukbb) {
      df_filtered %>%
        rename(d_ukbb = Cohens_D, se_ukbb = Cohens_D_SE)
    } else {
      df_filtered
    }
  })
}

# Extract the corresponding regions from FA_split and MD_split for both meta and UKBB top 5
FA_split_top_5_meta = extract_corresponding_regions(FA_split, top_5_FA_meta)
MD_split_top_5_meta = extract_corresponding_regions(MD_split, top_5_MD_meta)
FA_split_top_5_ukbb = extract_corresponding_regions(FA_split, top_5_FA_ukbb, is_ukbb = TRUE)
MD_split_top_5_ukbb = extract_corresponding_regions(MD_split, top_5_MD_ukbb, is_ukbb = TRUE)


# Calculate standard error for UKBB effect sizes
FA_split = lapply(FA_split, function(df) {
  df$Cohens_D_SE = (df$Cohens_D_Upper_CI - df$Cohens_D_Lower_CI) / (2 * 1.96)
  return(df)
})

MD_split = lapply(MD_split, function(df) {
  df$Cohens_D_SE = (df$Cohens_D_Upper_CI - df$Cohens_D_Lower_CI) / (2 * 1.96)
  return(df)
})

# Extract the corresponding regions from FA_split and MD_split
FA_split_top_5 = extract_corresponding_regions(FA_split, top_5_FA_meta)
MD_split_top_5 = extract_corresponding_regions(MD_split, top_5_MD_meta)

# For top_5_FA_meta
top_5_FA_meta$lower_CI = top_5_FA_meta$d_enigma - qnorm(0.975) * top_5_FA_meta$se_enigma
top_5_FA_meta$upper_CI = top_5_FA_meta$d_enigma + qnorm(0.975) * top_5_FA_meta$se_enigma

# For top_5_MD_meta
top_5_MD_meta$lower_CI = top_5_MD_meta$d_enigma - qnorm(0.975) * top_5_MD_meta$se_enigma
top_5_MD_meta$upper_CI = top_5_MD_meta$d_enigma + qnorm(0.975) * top_5_MD_meta$se_enigma

rename_columns = function(df) {
  if ("d_ukbb" %in% names(df)) {
    df = rename(df, Cohens_D = d_ukbb)
  }
  if ("se_ukbb" %in% names(df)) {
    df = rename(df, Cohens_D_SE = se_ukbb)
  }
  if ("lower_CI" %in% names(df)) {
    df = rename(df, Cohens_D_Lower_CI = lower_CI)
  }
  if ("upper_CI" %in% names(df)) {
    df = rename(df, Cohens_D_Upper_CI = upper_CI)
  }
  # Ensure "Full.tract.name" is consistent
  if ("Full tract name" %in% names(df)) {
    df = rename(df, Full.tract.name = `Full tract name`)
  }
  return(df)
}

# Rename columns in top_5_MD_ukbb
top_5_MD_ukbb = rename_columns(as.data.frame(top_5_MD_ukbb))

# Rename columns in MD_split_top_5_ukbb
MD_split_top_5_ukbb = lapply(MD_split_top_5_ukbb, rename_columns)

# Doing the same for FA
top_5_FA_ukbb = rename_columns(as.data.frame(top_5_FA_ukbb))
FA_split_top_5_ukbb = lapply(FA_split_top_5_ukbb, rename_columns)

###### Forest plots ######

# Function to create a forest plot
create_forest_plot = function(meta_df, split_obj, title) {
  # Combine the split objects into a single data frame
  combined_df = bind_rows(split_obj, .id = "Definition")
  
  # Merge with the meta-analytic results
  merged_df = left_join(combined_df, meta_df, by = "Region")
  
  # Reshape the data to have separate rows for d_ukbb and d_enigma
  data_long = merged_df %>%
    select(Region, Definition, d_ukbb, se_ukbb, d_enigma, se_enigma) %>%
    pivot_longer(cols = c(d_ukbb, se_ukbb, d_enigma, se_enigma),
                 names_to = c(".value", "Group"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(Group = case_when(
      Group == "ukbb"   ~ "UK Biobank",
      Group == "enigma" ~ "Enigma MDD",
      TRUE              ~ Group
    ),
    lower_CI = d - 1.96 * se,
    upper_CI = d + 1.96 * se)
  
  # Order the regions based on the meta-analytic effect size
  #data_long$Region = factor(data_long$Region, levels = meta_df$Region)
  data_long$Region = factor(meta_df$`Full tract name`[match(data_long$Region, meta_df$Region)], 
                            levels = meta_df$`Full tract name`)
  
  # Order the depression definitions
  data_long$Definition = factor(data_long$Definition, levels = definition_order)
  
  # Create a data frame for rectangle positions
  rect_data = data.frame(
    Definition = unique(data_long$Definition),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf
  )
  
  # Custom labeling function (vectorized)
  custom_labeller = function(values) {
    sapply(values, function(value) {
      if (value == "ICD10Dep_exclpsych") {
        return("ICD10Dep\n(exclpsych)")
      } else {
        return(as.character(value))
      }
    })
  }
  
  # Create the plot
  ggplot(data_long, aes(x = d, y = Region, color = Group, shape = Group)) +
    geom_pointrange(aes(xmin = lower_CI, xmax = upper_CI), position = position_dodge(width = 0.5), linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_rect(
      data = rect_data,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey90",
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    facet_grid(cols = vars(Definition), scales = "free_y", space = "free", 
               labeller = labeller(Definition = custom_labeller)) +
    labs(title = title, x = "Cohen's d (95% CI)", y = NULL) +
    scale_color_manual(values = c("UK Biobank" = "blue", "Enigma MDD" = "red")) +
    scale_shape_manual(values = c("UK Biobank" = 16, "Enigma MDD" = 17)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, color = color, face = "bold"),
      #plot.title = element_blank(),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 11),
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

# Function to create a forest plot for meta-analytic results
create_forest_plot_meta = function(meta_df, split_obj, title) {
  # Combine the split objects into a single data frame
  combined_df = bind_rows(split_obj, .id = "Definition")
  
  # Merge with the meta-analytic results
  merged_df = left_join(combined_df, meta_df, by = "Region")
  
  # Reshape the data to have separate rows for d_ukbb and d_enigma
  data_long = merged_df %>%
    select(Region, Definition, Cohens_D, Cohens_D_SE, d_enigma, se_enigma) %>%
    rename(d_ukbb = Cohens_D, se_ukbb = Cohens_D_SE) %>%
    pivot_longer(cols = c(d_ukbb, se_ukbb, d_enigma, se_enigma),
                 names_to = c(".value", "Group"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(Group = case_when(
      Group == "ukbb"   ~ "UK Biobank",
      Group == "enigma" ~ "Enigma MDD",
      TRUE              ~ Group
    ),
    lower_CI = d - 1.96 * se,
    upper_CI = d + 1.96 * se)
  
  # Order the regions based on the meta-analytic effect size
  data_long$Region = factor(meta_df$`Full tract name`[match(data_long$Region, meta_df$Region)], 
                            levels = meta_df$`Full tract name`)
  
  # Order the depression definitions
  data_long$Definition = factor(data_long$Definition, levels = definition_order)
  
  # Create a data frame for rectangle positions
  rect_data = data.frame(
    Definition = unique(data_long$Definition),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf
  )
  
  # Custom labeling function (vectorized)
  custom_labeller = function(values) {
    sapply(values, function(value) {
      if (value == "ICD10Dep_exclpsych") {
        return("ICD10Dep\n(exclpsych)")
      } else {
        return(as.character(value))
      }
    })
  }
  
  # Create the plot
  ggplot(data_long, aes(x = d, y = Region, color = Group, shape = Group)) +
    geom_pointrange(aes(xmin = lower_CI, xmax = upper_CI), position = position_dodge(width = 0.5), linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_rect(
      data = rect_data,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey90",
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    facet_grid(cols = vars(Definition), scales = "free_y", space = "free", 
               labeller = labeller(Definition = custom_labeller)) +
    labs(title = title, x = "Cohen's d (95% CI)", y = NULL) +
    scale_color_manual(values = c("UK Biobank" = "blue", "Enigma MDD" = "red")) +
    scale_shape_manual(values = c("UK Biobank" = 16, "Enigma MDD" = 17)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 11),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey90", colour = NA),
      panel.spacing = unit(0.5, "lines")
    )
}

create_forest_plot_ukbb = function(base_df, split_obj, title) {
  # Combine the split objects into a single data frame
  combined_df = bind_rows(split_obj, .id = "Definition")
  
  # Merge with the base results
  merged_df = left_join(combined_df, base_df, by = c("Region", "Full.tract.name"))
  
  # Reshape the data to have separate rows for UKBB and ENIGMA
  data_long = merged_df %>%
    select(Region, Full.tract.name, Definition, 
           d_ukbb = Cohens_D.x, se_ukbb = Cohens_D_SE.x,
           d_enigma = Cohens_D.y, se_enigma = Cohens_D_SE.y) %>%
    pivot_longer(cols = c(d_ukbb, se_ukbb, d_enigma, se_enigma),
                 names_to = c(".value", "Group"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(Group = case_when(
      Group == "ukbb"   ~ "UK Biobank",
      Group == "enigma" ~ "ENIGMA MDD",
      TRUE              ~ Group
    ),
    lower_CI = d - 1.96 * se,
    upper_CI = d + 1.96 * se)

  # Order the regions based on the effect size
  data_long$Full.tract.name = factor(data_long$Full.tract.name,
                                     levels = unique(base_df$Full.tract.name))
  
  # Order the depression definitions
  data_long$Definition = factor(data_long$Definition, levels = definition_order)
  
  # Custom labeling function (vectorized)
  custom_labeller = function(values) {
    sapply(values, function(value) {
      if (value == "ICD10Dep_exclpsych") {
        return("ICD10Dep\n(exclpsych)")
      } else {
        return(as.character(value))
      }
    })
  }
  
  # Create the plot
  ggplot(data_long, aes(x = d, y = Full.tract.name, color = Group, shape = Group)) +
    geom_pointrange(aes(xmin = lower_CI, xmax = upper_CI), position = position_dodge(width = 0.5), linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_grid(cols = vars(Definition), scales = "free_y", space = "free", 
               labeller = labeller(Definition = custom_labeller)) +
    labs(title = title, x = "Cohen's d (95% CI)", y = NULL) +
    scale_color_manual(values = c("UK Biobank" = "blue", "ENIGMA MDD" = "red")) +
    scale_shape_manual(values = c("UK Biobank" = 16, "ENIGMA MDD" = 17)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 11),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey90", colour = NA),
      panel.spacing = unit(0.5, "lines")
    )
}

###### Build and print forest plots ######

# Order depression definitions
definition_order = c("GPNoDep", "GPpsy", "Psypsy", "SelfRepDep", "DepAll", "ICD10Dep", "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur")


# Create forest plots for meta-analytic results
forest_plot_FA_meta = create_forest_plot_meta(
  top_5_FA_meta,
  FA_split_top_5_meta,
  "Top 5 Meta-Analytic Associations between Depression Phenotypes and Fractional Anisotropy (FA)"
)

forest_plot_MD_meta = create_forest_plot_meta(
  top_5_MD_meta,
  MD_split_top_5_meta,
  "Top 5 Meta-Analytic Associations between Depression Phenotypes and Mean Diffusivity (MD)"
)

# Create forest plots for UKBB results
# For MD
forest_plot_MD_ukbb = create_forest_plot_ukbb(
  MD_base,
  MD_split_top_5_ukbb,
  "Top 5 UK Biobank Associations between Depression Phenotypes and Mean Diffusivity (MD)"
)

# For FA
forest_plot_FA_ukbb = create_forest_plot_ukbb(
  FA_base,
  FA_split_top_5_ukbb,
  "Top 5 UK Biobank Associations between Depression Phenotypes and Fractional Anisotropy (FA)"
)

# Print the plots - for each plot, current version uses 1959 x 477
print(forest_plot_FA_meta)
print(forest_plot_MD_meta)
print(forest_plot_FA_ukbb)
print(forest_plot_MD_ukbb)


###### FA correlations ######

# Function to calculate correlations for each depression definition
calculate_correlations = function(FA_split, enigma_df) {
  correlations = lapply(FA_split, function(df) {
    # align UKB and ENIGMA by Region
    common = intersect(df$Region, enigma_df$Region)
    enigma_cohens_d = enigma_df$Cohens_D[match(common, enigma_df$Region)]
    df_cohens_d = df$Cohens_D[match(common, df$Region)]

    # spearman correlation and its p-value
    spearman_raw = cor(enigma_cohens_d, df_cohens_d, method = "spearman")
    spearman_p = cor.test(enigma_cohens_d, df_cohens_d, method = "spearman", exact = FALSE)$p.value

    list(spearman_raw = spearman_raw, spearman_p = spearman_p)
  })

  return(correlations)
}

# Calculate the correlations for FA
correlations_FA = calculate_correlations(FA_split, FA)

# Extract raw p-values for FA correlations
FA_pvals = sapply(correlations_FA, function(x) x$spearman_p)

# FDR correction for FA p-values
FA_fdr_pvals = p.adjust(FA_pvals, method = "fdr")

# Add FDR-corrected p-values to the correlation results
for (i in seq_along(correlations_FA)) {
  correlations_FA[[i]]$fdr_p = FA_fdr_pvals[i]
}

# Print the correlations
print(correlations_FA)

# Convert correlations_FA to a data frame
correlations_FA_df = do.call(rbind, lapply(names(correlations_FA), function(x) {
  data.frame(Definition = x, 
             Spearman_Raw = correlations_FA[[x]]$spearman_raw,
             Spearman_P = correlations_FA[[x]]$spearman_p,
             FDR_P = correlations_FA[[x]]$fdr_p,
             stringsAsFactors = FALSE)
}))

###### Correlation heatmap ######

# Create a heatmap for FA correlations
plot_correlation_heatmap = function(correlations_df, title, max_abs_corr) {
  # Reorder the rows of correlations_df based on order_definitions
  correlations_df = correlations_df[match(order_definitions, correlations_df$Definition), ]

  # build single-row correlation and p-value matrices
  corr_matrix = matrix(correlations_df$Spearman_Raw, nrow = 1)
  p_matrix = matrix(correlations_df$Spearman_P, nrow = 1)
  rownames(corr_matrix) = rownames(p_matrix) = "ENIGMA MDD"
  colnames(corr_matrix) = colnames(p_matrix) = correlations_df$Definition

  #labels = matrix(paste0(round(corr_matrix, 2), "\n(", format(p_matrix, digits = 3), ")"),
  #nrow = nrow(corr_matrix), ncol = ncol(corr_matrix))

  # cell labels show rounded correlations
  labels = matrix(round(corr_matrix, 2), nrow = nrow(corr_matrix), ncol = ncol(corr_matrix))

  # Create color breaks and custom colors
  # diverging blue-white-red scale centred on zero
  color_breaks = seq(-max_abs_corr, max_abs_corr, length.out = 101)  # Use odd number to ensure 0 is included
  neg_colors = colorRampPalette(c("#4a80ed", "white"))(51)  # 51 for negative values including 0
  pos_colors = colorRampPalette(c("white", "red"))(51)[-1]  # 50 for positive values, exclude first to avoid duplicate white
  custom_colors = c(neg_colors, pos_colors)

  # draw the heatmap, no clustering
  pheatmap(corr_matrix,
           #main = title, 
           color = custom_colors, 
           breaks = color_breaks,
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           display_numbers = labels, 
           number_format = "%.2f",
           fontsize_number = 10, 
           angle_col = 45, cellwidth = 80, cellheight = 35,
           fontcolor_number = "black", border_color = "black")
}

# Order definitions by phenotype granularity 
order_definitions = c("GPNoDep", "GPpsy", "Psypsy", "SelfRepDep", "DepAll", "ICD10Dep", "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur")

# Get the top 3 strongest correlations for FA
top_corr_FA = correlations_FA_df[order(-abs(correlations_FA_df$Spearman_Raw)), "Definition"][1:3]

# Calculate global x-range
calculate_global_x_range = function(top_corr_FA, FA_split) {
  all_x_values = unlist(lapply(top_corr_FA, function(definition) {
    FA_split[[definition]]$Cohens_D
  }))
  x_min = min(all_x_values, na.rm = TRUE)
  x_max = max(all_x_values, na.rm = TRUE)
  return(c(x_min, x_max))
}

# Calculate global y-range
calculate_global_y_range = function(top_corr_FA, FA, FA_split) {
  all_y_values = unlist(lapply(top_corr_FA, function(definition) {
    FA$Cohens_D[match(FA_split[[definition]]$Region, FA$Region)]
  }))
  y_min = min(all_y_values, na.rm = TRUE)
  y_max = max(all_y_values, na.rm = TRUE)
  return(c(y_min, y_max))
}

# Calculate global ranges
#global_x_range = c(-0.075, 0.025)  # Adding a buffer for x-range
global_x_range = calculate_global_x_range(top_corr_FA, FA_split)
global_y_range = calculate_global_y_range(top_corr_FA, FA, FA_split)

###### FA scatter plots ######

create_scatter_plot = function(definition, corr, p_value, measure, global_x_range, global_y_range, color, plot_title) {
  # pick UKB x and matched ENIGMA y for this measure
  data = if (measure == "FA") FA_split[[definition]] else MD_split[[definition]]
  y_data = if (measure == "FA") FA$Cohens_D[match(data$Region, FA$Region)] else MD$Cohens_D[match(data$Region, MD$Region)]

  # scatter with linear fit and fixed axis ranges
  ggplot(data, aes(x = Cohens_D, y = y_data)) +
    geom_point(size = 2, alpha = 0.7, color = color) +
    geom_smooth(method = "lm", se = FALSE, color = color, linewidth = 1) +
    xlab("Cohen's d (UKB)") +
    ylab("Cohen's d (ENIGMA)") +
    scale_y_continuous(breaks = seq(-0.25, -0.05, 0.05),
                       limits = c(-0.25, -0.05),
                       expand = c(0.01, 0)) +
    scale_x_continuous(breaks = seq(-0.075, 0.025, 0.025),
                       limits = global_x_range,
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
    # annotate r and p in the corner
    annotate("text", x = Inf, y = -Inf,
             label = paste0("r = ", round(corr, 3), "\np = ", format.pval(p_value, digits = 3)),
             hjust = 1.1, vjust = -0.8, size = 3, color = color) +
    labs(title = plot_title)  # Use labs() instead of ggtitle()
}

# define colors
plot_colors = c("#4a80ed", "#ed984a", "#45a049")

# Create scatter plots
scatter_plots_FA = lapply(seq_along(top_corr_FA), function(i) {
  definition = top_corr_FA[i]
  corr = correlations_FA_df$Spearman_Raw[correlations_FA_df$Definition == definition]
  p_value = correlations_FA_df$Spearman_P[correlations_FA_df$Definition == definition]
  print(paste("Plot", i, "corresponds to:", paste0("FA: ", definition)))
  create_scatter_plot(definition, corr, p_value, "FA", global_x_range, global_y_range, plot_colors[i], 
                      plot_title = paste0("(", letters[i], ") FA: ", definition))
})

# Arrange the FA scatter plots in a grid layout
scatter_plot_grid_FA = grid.arrange(
  grobs = scatter_plots_FA,
  ncol = 3
)

###### Mean diffusivity (MD) correlations and plots ######

##------------------------
# Calculating Heatmap and Scatterplots for MD Measures
##------------------------
# Function to split the averaged_results_df_MD dataframe into separate dataframes for each depression definition
split_dataframe_MD = function(averaged_results_df) {
  split(averaged_results_df, gsub("\\.[0-9]+$", "", averaged_results_df$Phenotypes))
}

# Split the averaged_results_df_MD dataframe
MD_split = split_dataframe_MD(averaged_results_df_MD)

# Function to calculate correlations for each depression definition
calculate_correlations_MD = function(MD_split, enigma_df) {
  correlations = lapply(MD_split, function(df) {
    # align UKB and ENIGMA by Region
    common = intersect(df$Region, enigma_df$Region)
    enigma_cohens_d = enigma_df$Cohens_D[match(common, enigma_df$Region)]
    df_cohens_d = df$Cohens_D[match(common, df$Region)]

    # spearman correlation and its p-value
    spearman_raw = cor(enigma_cohens_d, df_cohens_d, method = "spearman")
    spearman_p = cor.test(enigma_cohens_d, df_cohens_d, method = "spearman", exact = FALSE)$p.value

    list(spearman_raw = spearman_raw, spearman_p = spearman_p)
  })
  
  return(correlations)
}

# Calculate the correlations for MD
correlations_MD = calculate_correlations_MD(MD_split, MD)

# Extract raw p-values for MD correlations
MD_pvals = sapply(correlations_MD, function(x) x$spearman_p)

# FDR correction for MD p-values
MD_fdr_pvals = p.adjust(MD_pvals, method = "fdr")

# Add FDR-corrected p-values to the correlation results
for (i in seq_along(correlations_MD)) {
  correlations_MD[[i]]$fdr_p = MD_fdr_pvals[i]
}

# Print the correlations
print(correlations_MD)

# global MD ranges computed below, after top_corr_MD is defined


# Convert correlations_MD to a data frame
correlations_MD_df = do.call(rbind, lapply(names(correlations_MD), function(x) {
  data.frame(Definition = x, 
             Spearman_Raw = correlations_MD[[x]]$spearman_raw,
             Spearman_P = correlations_MD[[x]]$spearman_p,
             FDR_P = correlations_MD[[x]]$fdr_p,
             stringsAsFactors = FALSE)
}))


# Calculate the maximum absolute correlation across both FA and MD
max_abs_corr = max(abs(c(correlations_FA_df$Spearman_Raw, correlations_MD_df$Spearman_Raw)))

# Plot the heatmap for FA correlations
plot_correlation_heatmap(correlations_FA_df, "Spearman Correlations (FA)", max_abs_corr)

# Plot the heatmap for MD correlations
plot_correlation_heatmap(correlations_MD_df, "Spearman Correlations (MD)", max_abs_corr)

# Get the top 3 strongest correlations for MD
top_corr_MD = correlations_MD_df[order(-abs(correlations_MD_df$Spearman_Raw)), "Definition"][1:3]

create_scatter_plot = function(definition, corr, p_value, measure, color, plot_title) {
  # pick UKB x and matched ENIGMA y for this measure
  data = if (measure == "FA") FA_split[[definition]] else MD_split[[definition]]
  y_data = if (measure == "FA") FA$Cohens_D[match(data$Region, FA$Region)] else MD$Cohens_D[match(data$Region, MD$Region)]

  # Calculate the range of the data
  x_range = range(data$Cohens_D, na.rm = TRUE)
  y_range = range(y_data, na.rm = TRUE)

  # Add some padding to the ranges
  x_padding = diff(x_range) * 0.05
  y_padding = diff(y_range) * 0.05

  # scatter with linear fit, axes scaled to padded data range
  ggplot(data, aes(x = Cohens_D, y = y_data)) +
    geom_point(size = 2, alpha = 0.7, color = color) +
    geom_smooth(method = "lm", se = FALSE, color = color, linewidth = 1) +
    xlab("Cohen's d (UKB)") +
    ylab("Cohen's d (ENIGMA)") +
    scale_y_continuous(limits = c(y_range[1] - y_padding, y_range[2] + y_padding),
                       expand = c(0, 0)) +
    scale_x_continuous(limits = c(x_range[1] - x_padding, x_range[2] + x_padding),
                       labels = function(x) sprintf("%.3f", x),
                       expand = c(0, 0)) +
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
    # annotate r and p in the corner
    annotate("text", x = Inf, y = -Inf,
             label = paste0("r = ", round(corr, 3), "\np = ", format.pval(p_value, digits = 3)),
             hjust = 1.1, vjust = -0.8, size = 3, color = color) +
    labs(title = plot_title)
}

# Create scatter plots for MD
scatter_plots_MD = lapply(seq_along(top_corr_MD), function(i) {
  definition = top_corr_MD[i]
  corr = correlations_MD_df$Spearman_Raw[correlations_MD_df$Definition == definition]
  p_value = correlations_MD_df$Spearman_P[correlations_MD_df$Definition == definition]
  print(paste("Plot", i, "corresponds to:", paste0("MD: ", definition)))
  create_scatter_plot(definition, corr, p_value, "MD", plot_colors[i], 
                      plot_title = paste0("(", letters[i], ") MD: ", definition))
})

# Arrange the MD scatter plots in a grid layout
scatter_plot_grid_MD = grid.arrange(
  grobs = scatter_plots_MD,
  ncol = 3
)

# Combine FA and MD scatter plot grids vertically
combined_scatter_plot_grid = grid.arrange(
  scatter_plot_grid_FA,
  scatter_plot_grid_MD,
  ncol = 1
)


# Get the top 3 strongest correlations for FA and MD
top_corr_FA = correlations_FA_df[order(-abs(correlations_FA_df$Spearman_Raw)), "Definition"][1:3]
top_corr_MD = correlations_MD_df[order(-abs(correlations_MD_df$Spearman_Raw)), "Definition"][1:3]

# Calculate global x-range
calculate_global_x_range = function(top_corr, split_data) {
  all_x_values = unlist(lapply(top_corr, function(definition) {
    split_data[[definition]]$Cohens_D
  }))
  x_min = min(all_x_values, na.rm = TRUE)
  x_max = max(all_x_values, na.rm = TRUE)
  return(c(x_min, x_max))
}

# Calculate global y-range
calculate_global_y_range = function(top_corr, data, split_data) {
  all_y_values = unlist(lapply(top_corr, function(definition) {
    data$Cohens_D[match(split_data[[definition]]$Region, data$Region)]
  }))
  y_min = min(all_y_values, na.rm = TRUE)
  y_max = max(all_y_values, na.rm = TRUE)
  return(c(y_min, y_max))
}

# Calculate global ranges for FA and MD
global_x_range_FA = calculate_global_x_range(top_corr_FA, FA_split)
global_y_range_FA = calculate_global_y_range(top_corr_FA, FA, FA_split)
global_x_range_MD = calculate_global_x_range(top_corr_MD, MD_split)
global_y_range_MD = calculate_global_y_range(top_corr_MD, MD, MD_split)

create_scatter_plot = function(definition, corr, p_value, measure, global_x_range, global_y_range, color, plot_title) {
  # pick UKB x and matched ENIGMA y for this measure
  data = if (measure == "FA") FA_split[[definition]] else MD_split[[definition]]
  y_data = if (measure == "FA") FA$Cohens_D[match(data$Region, FA$Region)] else MD$Cohens_D[match(data$Region, MD$Region)]

  # scatter with linear fit on shared global axis ranges
  ggplot(data, aes(x = Cohens_D, y = y_data)) +
    geom_point(size = 2, alpha = 0.7, color = color) +
    geom_smooth(method = "lm", se = FALSE, color = color, linewidth = 1) +
    xlab("Cohen's d (UKB)") +
    ylab("Cohen's d (ENIGMA)") +
    scale_y_continuous(breaks = seq(min(global_y_range), max(global_y_range), length.out = 5),
                       limits = global_y_range,
                       expand = c(0.01, 0)) +
    scale_x_continuous(breaks = seq(min(global_x_range), max(global_x_range), length.out = 5),
                       limits = global_x_range,
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
    # annotate r and p in the corner
    annotate("text", x = Inf, y = -Inf,
             label = paste0("r = ", round(corr, 3), "\np = ", format.pval(p_value, digits = 3)),
             hjust = 1.1, vjust = -0.8, size = 3, color = color) +
    labs(title = plot_title)
}

# Define colors
plot_colors = c("#4a80ed", "#ed984a", "#45a049")

# Create scatter plots for FA
scatter_plots_FA = lapply(seq_along(top_corr_FA), function(i) {
  definition = top_corr_FA[i]
  corr = correlations_FA_df$Spearman_Raw[correlations_FA_df$Definition == definition]
  p_value = correlations_FA_df$Spearman_P[correlations_FA_df$Definition == definition]
  print(paste("Plot", i, "corresponds to:", paste0("FA: ", definition)))
  create_scatter_plot(definition, corr, p_value, "FA", global_x_range_FA, global_y_range_FA, plot_colors[i], 
                      plot_title = paste0("(", letters[i], ") FA: ", definition))
})

# Create scatter plots for MD
scatter_plots_MD = lapply(seq_along(top_corr_MD), function(i) {
  definition = top_corr_MD[i]
  corr = correlations_MD_df$Spearman_Raw[correlations_MD_df$Definition == definition]
  p_value = correlations_MD_df$Spearman_P[correlations_MD_df$Definition == definition]
  print(paste("Plot", i, "corresponds to:", paste0("MD: ", definition)))
  create_scatter_plot(definition, corr, p_value, "MD", global_x_range_MD, global_y_range_MD, plot_colors[i], 
                      plot_title = paste0("(", letters[i], ") MD: ", definition))
})

# Arrange the FA scatter plots in a grid layout
scatter_plot_grid_FA = grid.arrange(
  grobs = scatter_plots_FA,
  ncol = 3
)

# Arrange the MD scatter plots in a grid layout
scatter_plot_grid_MD = grid.arrange(
  grobs = scatter_plots_MD,
  ncol = 3
)
