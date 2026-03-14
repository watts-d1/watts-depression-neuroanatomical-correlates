# Neuroimaging Feature Extraction Pipeline - Heterogeneity of Depression
# Part 1 - Parallel Processing version
#
# Neuroimaging preprocessing: https://doi.org/10.1016/j.neuroimage.2017.10.034
# Author: Devon Watts

# --- User-defined paths (modify for your environment) ---
UKBB_DATA_PATH = "/path/to/ukbiobank/ukb674571.csv"  # UK Biobank phenotype CSV

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
             'car', 'future', 'future.apply', 'parallelly', 'parallel', 'furrr', 
             'ggplot2', 'pbapply', 'progressr', 'gridExtra', 'pbapply', 'vctrs')

package_install(packages)
library(here)

setDTthreads(threads = 20)

# Read a file and add a subject_num column from whichever ID column is present
read_and_examine = function(file_path, default_id_column = "fussIID", alt_id_column1 = "IIDs", alt_id_column2 = "IID", alt_id_column3 = "eid") {
  data = fread(file_path, header = TRUE, sep = "auto")
  
  if (default_id_column %in% names(data)) {
    id_column = default_id_column
  } else if (alt_id_column1 %in% names(data)) {
    id_column = alt_id_column1
  } else if (alt_id_column2 %in% names(data)) {
    id_column = alt_id_column2
  } else if (alt_id_column3 %in% names(data)) {
    id_column = alt_id_column3
  } else {
    id_column = NA
    warning(paste0("ID column not found in ", file_path))
  }
  
  if (!is.na(id_column)) {
    data[, subject_num := get(id_column)]
  } else {
    warning(paste0("ID column not found in ", file_path))
  }
  
  print(str(data))
  
  return(data)
}

file_paths = c('fuss.mappedDepressionPhenos.txt',
               'fuss.mappedImputedDepression.txt', 
               'rawExtracted_pcs_covs.txt', 
               'UKB.depressionphenotypes.txt', 
               'ukb_md_imputed.txt')

data_tables = lapply(file_paths, read_and_examine)

names(data_tables) = gsub("\\.txt$", "", file_paths)

fuss_mapped_depression_phenos  = data_tables$fuss.mappedDepressionPhenos
fuss_mapped_imputed_depression = data_tables$fuss.mappedImputedDepression
ukb_depression_phenotypes      = data_tables$UKB.depressionphenotypes
ukb_md_imputed                 = data_tables$ukb_md_imputed
ukbiobank_PCs                  = data_tables$rawExtracted_pcs_covs

ukbiobank_PCs = ukbiobank_PCs %>%
  rename_with(~ gsub("f.22009.0.", "PC", .), starts_with("f.22009.0."))

subset_depression_data = function(data, diagnosis_column) {
  if (!("subject_num" %in% names(data))) {
    stop("subject_num column not found in the dataset")
  }
  
  subset_data = na.omit(data[, c("subject_num", diagnosis_column), with = FALSE])
  
  cat("Original data (rows, columns):", dim(data)[1], ",", dim(data)[2], "\n")
  cat("Subset data for", diagnosis_column, "(rows, columns):", dim(subset_data)[1], ",", dim(subset_data)[2], "\n")
  
  return(subset_data)
}

# Depression diagnosis columns
depression_cols = c("LifetimeMDD", "MDDRecur", "GPpsy", "Psypsy", "DepAll", 
                    "SelfRepDep", "GPNoDep", "ICD10Dep", "ICD10Dep.exclpsych")

depression_subsets = lapply(depression_cols, function(col) subset_depression_data(fuss_mapped_depression_phenos, col))
names(depression_subsets) = depression_cols

for (col in depression_cols) {
  assign(paste0("df_", col), depression_subsets[[col]])
}

# List of all fields: https://biobank.ctsu.ox.ac.uk/showcase/label.cgi?id=100003
csv_file_path = UKBB_DATA_PATH
csv_data = fread(csv_file_path, showProgress = TRUE, nThread = 20)

csv_data$subject_num = csv_data$eid

# Subset UK Biobank data to columns matching specified field IDs across timepoints
subset_ukbb_data = function(data, field_ids, always_include = "subject_num") {
  pattern = paste0("^(", paste(field_ids, collapse = "|"), ")(-[0-9]+\\.[0-9]+)?$")
  subset_cols = grep(pattern, names(data), value = TRUE)
  if (!always_include %in% subset_cols) {
    subset_cols = c(always_include, subset_cols)
  }
  
  subset_data = data[, subset_cols, with = FALSE]
  
  return(subset_data)
}

# Define field IDs for T1 global structural MRI measurements
# https://biobank.ctsu.ox.ac.uk/showcase/label.cgi?id=110
# All volume metrics 
# Instance 2 is considered the baseline imaging visit 
t1_global_sMRI = c("25010", "25009", "25008", "25007", "25006", "25005", "25002", "25001", "25004", "25003")

t1_global_sMRI_data = subset_ukbb_data(csv_data, t1_global_sMRI)

# Define field IDs for T1 subcortical structural MRI measurements
# These fields cover volumes of key subcortical structures including:
# The accumbens, amygdala, caudate, hippocampus, pallidum, putamen, and thalamus.
# https://biobank.ctsu.ox.ac.uk/showcase/label.cgi?id=190
# Instance 2 is considered the baseline imaging visit 
# All volume metrics 
t1_subcortical_sMRI = c("26564", "26595", "26563", "26594", "26559", "26590", "26562", "26593", "26561", "26592", "26560", "26591", "26558", "26589")

t1_subcortical_sMRI_data = subset_ukbb_data(csv_data, t1_subcortical_sMRI)

# Define field IDs for T1 MRI thickness and area measurements
# https://biobank.ctsu.ox.ac.uk/showcase/label.cgi?id=192
# These fields include area and mean thickness measurements of various brain regions in both hemispheres.
# Instance 2 is considered the baseline imaging visit 
t1_thickness = c(
  
  # Mean Thickness 
  "26755", "26856", "26756", "26857", "26757", "26858", "26758", "26859", "26759", "26860", 
  "26760", "26861", "26786", "26887", "26761", "26862", "26762", "26863", "26763", "26864", 
  "26788", "26889", "26764", "26865", "26765", "26866", "26766", "26867", "26767", "26868", 
  "26768", "26869", "26769", "26870", "26771", "26872", "26770", "26871", "26772", "26873", 
  "26773", "26874", "26774", "26875", "26775", "26876", "26776", "26877", "26777", "26878", 
  "26778", "26879", "26779", "26880", "26780", "26881", "26781", "26882", "26782", "26883", 
  "26783", "26884", "26784", "26885", "26785", "26886", "26787", "26888")

# Surface Area
t1_surface_area = c("26721", "26822", "26722", "26823", "26723", "26824", "26724", "26825", "26725", "26826",
                    "26726", "26827", "26752", "26853", "26727", "26828", "26728", "26829", "26729", "26830",
                    "26754", "26855", "26730", "26831", "26731", "26832", "26732", "26833", "26733", "26834",
                    "26734", "26835", "26735", "26836", "26737", "26838", "26736", "26837", "26738", "26839",
                    "26739", "26840", "26740", "26841", "26741", "26842", "26742", "26843", "26743", "26844",
                    "26744", "26845", "26745", "26846", "26746", "26847", "26747", "26848", "26748", "26849",
                    "26749", "26850", "26750", "26851", "26751", "26852", "26753", "26854")

t1_thickness_sMRI_data = subset_ukbb_data(csv_data, t1_thickness)
t1_surface_area_sMRI_data = subset_ukbb_data(csv_data, t1_surface_area)

# Define field IDs for DTI (Diffusion Tensor Imaging) measurements
# https://biobank.ctsu.ox.ac.uk/showcase/label.cgi?id=134
# These fields represent Mean Fractional Anisotropy (FA) values across various brain white matter tracts.
# Instance 2 is considered the baseline visit 
dti_measures = c(
  # FA
  "25079", "25078", "25073", "25072", "25059", "25071", "25070", "25091", "25090", "25093",
  "25092", "25063", "25062", "25089", "25088", "25095", "25094", "25061", "25058", "25067",
  "25066", "25065", "25064", "25056", "25057", "25083", "25082", "25075", "25074", "25085",
  "25084", "25077", "25076", "25087", "25086", "25060", "25069", "25068", "25081", "25080",
  "25099", "25098", "25097", "25096", "25103", "25102", "25101", "25100",
  # MD measures 
  "25127", "25126", "25121", "25120", "25107", "25119", "25118", "25139", "25138", "25141", "25140", 
  "25111", "25110", "25137", "25136", "25143", "25142", "25109", "25106", "25115", "25114", "25113", 
  "25112", "25104", "25105", "25131", "25130", "25123", "25122", "25133", "25132", "25125", "25124", 
  "25135", "25134", "25108", "25117", "25116", "25129", "25128", "25147", "25146", "25145", "25144", 
  "25151", "25150", "25149", "25148"
)

dti_measures_data = subset_ukbb_data(csv_data, dti_measures)

# Merge a list of dataframes by a common column
merge_data_by_column = function(data_list, by_column = "subject_num") {
  Reduce(function(x, y) merge(x, y, by = by_column, all = TRUE), data_list)
}

combined_neuro_data = merge_data_by_column(list(t1_global_sMRI_data, t1_subcortical_sMRI_data, t1_thickness_sMRI_data, t1_surface_area_sMRI_data, dti_measures_data))
cat("The initial neuroimaging features across baseline and follow-up has", length(combined_neuro_data), "columns in total.\n")

# For these fields, 2.0 is baseline and 3.0 is follow-up
baseline_combined_neuro = combined_neuro_data %>% select(ends_with("2.0"))
cat("The initial neuroimaging features at baseline has", length(baseline_combined_neuro), "columns in total.\n")

baseline_combined_neuro$subject_num = combined_neuro_data$subject_num

# Covariates for baseline encounter (2.0)
t1_covariates = c("21022",  # Age
                  "31",     # Sex
                  "54",     # Assessment Center
                  "26521",  # Estimated Total Intracranial Volume 
                  "24419")  # Head Motion 

t1_covariates_data = csv_data %>%
  select(subject_num, starts_with(t1_covariates))

t1_covariates_data = t1_covariates_data %>% select(c("54-0.0", "31-0.0", "26521-2.0", "21022-0.0", "24419-2.0", "subject_num"))

PC_data = ukbiobank_PCs %>%
  select(subject_num, matches("^PC[1-9]$|^PC1[0-9]$|^PC20$"))

PC_data$subject_num = as.integer(PC_data$subject_num)
t1_covariates_data = left_join(t1_covariates_data, PC_data, by = "subject_num")
baseline_neuroimaging_and_covariates = merge(baseline_combined_neuro, t1_covariates_data, by = "subject_num")
baseline_neuroimaging_and_covariates$subject_num = as.integer(baseline_neuroimaging_and_covariates$subject_num)

# Standardize column names: prefix numeric IDs with "bl_" and strip timepoint suffixes
colnames(baseline_neuroimaging_and_covariates) = gsub("^([0-9]+)", "bl_\\1", colnames(baseline_neuroimaging_and_covariates))
colnames(baseline_neuroimaging_and_covariates) = gsub("-[0-9]+\\.[0-9]+$", "", colnames(baseline_neuroimaging_and_covariates))

colnames(t1_covariates_data) = gsub("^([0-9]+)", "bl_\\1", colnames(t1_covariates_data))
colnames(t1_covariates_data) = gsub("-[0-9]+\\.[0-9]+$", "", colnames(t1_covariates_data))

colnames(baseline_combined_neuro) = gsub("^([0-9]+)", "bl_\\1", colnames(baseline_combined_neuro))
colnames(baseline_combined_neuro) = gsub("-[0-9]+\\.[0-9]+$", "", colnames(baseline_combined_neuro))

cat("The combined neuroimaging and covariates data has", nrow(baseline_neuroimaging_and_covariates), "rows and", ncol(baseline_neuroimaging_and_covariates) -1, "columns.\n")

### Approximately 90% of rows do not have Estimated Total Intracranial Volume 
### Note: This is consistent with overall percentage of complete cases (~9.4%)
t1_covariates_data %>%
  summarise(
    missing_rows = sum(is.na(bl_26521)),
    complete_rows = sum(!is.na(bl_26521)),
    percent_missing = (sum(is.na(bl_26521)) / nrow(t1_covariates_data)) * 100
  )

# Subset neuroimaging data for each depression definition
subset_and_combine = function(base_data, smaller_data, by_column = "subject_num") {
  # Filter the base data for subject numbers in the smaller data
  filtered_base = base_data %>%
    filter((!!sym(by_column)) %in% smaller_data[[by_column]])
  
  # Combine the filtered base data with the depression dataframes 
  combined_data = filtered_base %>%
    left_join(smaller_data, by = by_column)
  
  return(combined_data)
}

neuro_GPNoDep            = subset_and_combine(baseline_neuroimaging_and_covariates, df_GPNoDep)
neuro_GPpsy              = subset_and_combine(baseline_neuroimaging_and_covariates, df_GPpsy)
neuro_Psypsy             = subset_and_combine(baseline_neuroimaging_and_covariates, df_Psypsy)
neuro_self_rep_dep       = subset_and_combine(baseline_neuroimaging_and_covariates, df_SelfRepDep)
neuro_DepAll             = subset_and_combine(baseline_neuroimaging_and_covariates, df_DepAll)
neuro_ICD10Dep           = subset_and_combine(baseline_neuroimaging_and_covariates, df_ICD10Dep)
neuro_ICD10Dep_exclpsych = subset_and_combine(baseline_neuroimaging_and_covariates, df_ICD10Dep.exclpsych)
neuro_lifetime_mdd       = subset_and_combine(baseline_neuroimaging_and_covariates, df_LifetimeMDD)
neuro_mdd_recurr         = subset_and_combine(baseline_neuroimaging_and_covariates, df_MDDRecur)

# Calculate and print per-column missingness percentages
calculate_missingness = function(data) {
  
  # Exclude PC columns
  #data_no_PC = data %>% select(-starts_with("PC"))
  
  # Total number of rows
  total_rows = nrow(data)
  
  # Number of complete cases
  complete_cases = sum(complete.cases(data))
  
  # Calculate the percentage of missing data for each column
  missing_percentage = colSums(is.na(data)) / total_rows * 100
  
  # Print results
  cat("Total number of cases:", total_rows, "\n")
  cat("Number of complete cases:", complete_cases, "\n")
  cat("Percentage of complete cases:", (complete_cases / total_rows * 100), "%\n")
  cat("Missing data by column (%):\n")
  print(missing_percentage)
  
  # Optionally, return a list with complete cases and missing percentages
  return(list(complete_cases = complete_cases, missing_percentage = missing_percentage))
}

# Missingness across depression definitions
missing_info_GPNoDep            = calculate_missingness(neuro_GPNoDep)
missing_info_GPpsy              = calculate_missingness(neuro_GPpsy)
missing_info_psypsy             = calculate_missingness(neuro_Psypsy)
missing_info_self_rep_dep       = calculate_missingness(neuro_self_rep_dep)
missing_info_DepAll             = calculate_missingness(neuro_DepAll)
missing_info_ICD10Dep           = calculate_missingness(neuro_ICD10Dep)
missing_info_ICD10Dep_exclpsych = calculate_missingness(neuro_ICD10Dep_exclpsych)
missing_info_lifetime_mdd       = calculate_missingness(neuro_lifetime_mdd)
missing_info_mdd_recurr         = calculate_missingness(neuro_mdd_recurr)

# Subset to complete cases and report case/control counts (last column must be binary 0/1)
#'
subset_complete_cases = function(data) {
  # Explicitly convert the last column to a vector
  last_column = as.vector(data[[ncol(data)]])
  
  # Check if the last column is numeric (integer or double) and binary (0 or 1)
  if (!(is.integer(last_column) || is.numeric(last_column)) || !all.equal(unique(last_column), c(0, 1))) {
    stop("The last column must be numeric (integer or double) with 0s and 1s only.")
  }
  
  # Subset to complete cases
  complete_data = data[complete.cases(data), ]
  
  # Count the number of controls (0) and cases (1)
  num_controls = sum(complete_data[[ncol(complete_data)]] == 0)
  num_cases    = sum(complete_data[[ncol(complete_data)]] == 1)
  
  # Print the counts
  #cat("\nSubsetted Dataframe:\n")
  cat("Total Subject Count:", nrow(complete_data), "\n")
  cat("Number of controls:", num_controls, "\n")
  cat("Number of cases:", num_cases, "\n")
  
  # Return a list containing the subsetted data frame and counts
  return(list(completed_df = complete_data, num_controls = num_controls, num_cases = num_cases))
}

# Complete cases for each depression definition
completed_DepAll     = subset_complete_cases(neuro_DepAll)
completed_GPNoDep    = subset_complete_cases(neuro_GPNoDep)
completed_Psypsy     = subset_complete_cases(neuro_Psypsy)
completed_SelfRepDep = subset_complete_cases(neuro_self_rep_dep)
completed_GPpsy      = subset_complete_cases(neuro_GPpsy)
completed_ICD10Dep   = subset_complete_cases(neuro_ICD10Dep)
completed_ICD10Dep_exclpsych = subset_complete_cases(neuro_ICD10Dep_exclpsych)
completed_LifetimeMDD        = subset_complete_cases(neuro_lifetime_mdd)
completed_MDDRecurr          = subset_complete_cases(neuro_mdd_recurr)

save(completed_DepAll, 
     completed_GPNoDep,
     completed_Psypsy,
     completed_SelfRepDep,
     completed_GPpsy,
     completed_ICD10Dep,
     completed_ICD10Dep_exclpsych,
     completed_LifetimeMDD,
     completed_MDDRecurr,
     file = here("output", "S1_neuroimaging_objects", "pre_z_score_neuroimaging_objects.RData")
)

# Standardizing features to N(0, 1) using z-scores

# Z-score standardize columns with 'bl_' prefix
standardize_bl_columns = function(df) {
  bl_columns = grep("^bl_", names(df), value = TRUE) # Identifying 'bl_' columns
  df[, bl_columns] = lapply(df[, bl_columns, with = FALSE], scale) # Applying Z-score standardization
  return(df)
}

original_dataframes = list(
  completed_GPNoDep = completed_GPNoDep$completed_df,
  completed_GPpsy      = completed_GPpsy$completed_df,
  completed_Psypsy  = completed_Psypsy$completed_df,
  completed_SelfRepDep = completed_SelfRepDep$completed_df,
  completed_DepAll  = completed_DepAll$completed_df,
  completed_ICD10Dep   = completed_ICD10Dep$completed_df,
  completed_ICD10Dep_exclpsych = completed_ICD10Dep_exclpsych$completed_df,
  completed_LifetimeMDD        = completed_LifetimeMDD$completed_df,
  completed_MDDRecurr          = completed_MDDRecurr$completed_df
)

for (name in names(original_dataframes)) {
  df_standardized = standardize_bl_columns(original_dataframes[[name]])
  assign(paste0("z_score_", name), df_standardized)
}

save(
  z_score_completed_GPNoDep,
  z_score_completed_GPpsy,
  z_score_completed_Psypsy,
  z_score_completed_SelfRepDep,
  z_score_completed_DepAll, 
  z_score_completed_ICD10Dep,
  z_score_completed_ICD10Dep_exclpsych,
  z_score_completed_LifetimeMDD,
  z_score_completed_MDDRecurr,
  file = here("output", "S1_neuroimaging_objects", "z_score_neuroimaging_objects.RData")
)

# Separate bl_-prefixed features from covariates and outcome into a named list
separate_features_and_covariates = function(dataframe, covariates, outcome_col) {
  baseline_pattern = "^bl_"
  baseline_cols = grep(baseline_pattern, names(dataframe), value = TRUE)
  exclude_cols = c(covariates, paste0("PC", 1:20), outcome_col, "subject_num")
  final_features = setdiff(baseline_cols, exclude_cols)

  features_df   = dataframe[, ..final_features, with = FALSE]
  covariates_df = dataframe[, ..covariates, with = FALSE]
  outcome_df    = dataframe[, ..outcome_col, with = FALSE]
  subject_num_df = dataframe[, .(subject_num = subject_num)]
  
  return(list(features = features_df, covariates = covariates_df, outcome = outcome_df, subject_num = subject_num_df))
}

# Define covariates
# age, sex, assessment center, estimated total intracranial volume, and head motion
covariates = c("bl_31", "bl_54", "bl_21022", "bl_26521", paste0("PC", 1:20))

# Separate features and covariates for each dataset
separated_GPNoDep    = separate_features_and_covariates(z_score_completed_GPNoDep, covariates, "GPNoDep")
separated_GPpsy      = separate_features_and_covariates(z_score_completed_GPpsy, covariates, "GPpsy")
separated_Psypsy     = separate_features_and_covariates(z_score_completed_Psypsy, covariates, "Psypsy")
separated_SelfRepDep = separate_features_and_covariates(z_score_completed_SelfRepDep, covariates, "SelfRepDep")
separated_DepAll     = separate_features_and_covariates(z_score_completed_DepAll, covariates, "DepAll")
separated_ICD10Dep   = separate_features_and_covariates(z_score_completed_ICD10Dep, covariates, "ICD10Dep")
separated_ICD10Dep_exclpsych = separate_features_and_covariates(z_score_completed_ICD10Dep_exclpsych, covariates, "ICD10Dep.exclpsych")
separated_LifetimeMDD        = separate_features_and_covariates(z_score_completed_LifetimeMDD, covariates, "LifetimeMDD")
separated_MDDRecurr          = separate_features_and_covariates(z_score_completed_MDDRecurr, covariates, "MDDRecur")

save(separated_GPNoDep,
     separated_GPpsy,
     separated_Psypsy,
     separated_SelfRepDep,
     separated_DepAll, 
     separated_ICD10Dep,
     separated_ICD10Dep_exclpsych,
     separated_LifetimeMDD,
     separated_MDDRecurr,
     file = here("output", "S1_neuroimaging_objects", "separated_feature_and_covariate_objects.RData")
)

# Gathering measurement IDs
t1_global_sMRI      = c("25010", "25009", "25008", "25007", "25006", "25005", "25002", "25001", "25004", "25003")
t1_subcortical_sMRI = c("26564", "26595", "26563", "26594", "26559", "26590", "26562", "26593", "26561", "26592", 
                        "26560", "26591", "26558", "26589")

surface_area_ids    = c("26721", "26822", "26722", "26823", "26723", "26824", "26724", "26825", "26725", "26826", "26726",
                        "26827", "26752", "26853", "26727", "26828", "26728", "26829", "26729", "26830", "26754", "26855", 
                        "26730", "26831", "26731", "26832", "26732", "26833", "26733", "26834", "26734", "26835", "26735", 
                        "26836", "26737", "26838", "26736", "26837", "26738", "26839", "26739", "26840", "26740", "26841", 
                        "26741", "26842", "26742", "26843", "26743", "26844", "26744", "26845", "26745", "26846", "26746", 
                        "26847", "26747", "26848", "26748", "26849", "26749", "26850", "26750", "26851", "26751", "26852", 
                        "26753", "26854")

thickness_ids       = c("26755", "26856", "26756", "26857", "26757", "26858", "26758", "26859", "26759", "26860", "26760", 
                        "26861", "26786", "26887", "26761", "26862", "26762", "26863", "26763", "26864", "26788", "26889", 
                        "26764", "26865", "26765", "26866", "26766", "26867", "26767", "26868", "26768", "26869", "26769", 
                        "26870", "26771", "26872", "26770", "26871", "26772", "26873", "26773", "26874", "26774", "26875", 
                        "26775", "26876", "26776", "26877", "26777", "26878", "26778", "26879", "26779", "26880", "26780", 
                        "26881", "26781", "26882", "26782", "26883", "26783", "26884", "26784", "26885", "26785", "26886", 
                        "26787", "26888")

fa_measures         = c("25079", "25078", "25073", "25072", "25059", "25071", "25070", "25091", "25090", "25093", "25092", 
                        "25063", "25062", "25089", "25088", "25095", "25094", "25061", "25058", "25067", "25066", "25065", 
                        "25064", "25056", "25057", "25083", "25082", "25075", "25074", "25085", "25084", "25077", "25076", 
                        "25087", "25086", "25060", "25069", "25068", "25081", "25080", "25099", "25098", "25097", "25096", 
                        "25103", "25102", "25101", "25100")


md_measures         = c("25127", "25126", "25121", "25120", "25107", "25119", "25118", "25139", "25138", "25141", "25140", 
                        "25111", "25110", "25137", "25136", "25143", "25142", "25109", "25106", "25115", "25114", "25113", 
                        "25112", "25104", "25105", "25131", "25130", "25123", "25122", "25133", "25132", "25125", "25124", 
                        "25135", "25134", "25108", "25117", "25116", "25129", "25128", "25147", "25146", "25145", "25144", 
                        "25151", "25150", "25149", "25148"
)



# Measurement types to analyze
measurement_types = c("global_volume", "subcortical_volume", "surface_area", "thickness", "FA", "MD")

# Return bl_-prefixed measurement IDs for a given neuroimaging type
get_measurement_ids = function(measurement_type) {
  if (measurement_type == "global_volume") {
    return(paste0("bl_", t1_global_sMRI))
  } else if (measurement_type == "subcortical_volume") {
    return(paste0("bl_", t1_subcortical_sMRI))
  } else if (measurement_type == "surface_area") {
    return(paste0("bl_", surface_area_ids))
  } else if (measurement_type == "thickness") {
    return(paste0("bl_", thickness_ids))
  } else if (measurement_type == "FA") {
    return(paste0("bl_", fa_measures))
  } else if (measurement_type == "MD") {
    return(paste0("bl_", md_measures))
  } else {
    stop("Invalid measurement type")
  }
}
# Subset data to diagnosis column, covariates, and measurement-type-specific features
subset_data_by_type = function(neuro_data, measurement_types, diagnosis_col, covariates) {
  if (!is.data.table(neuro_data)) {
    neuro_data = as.data.table(neuro_data)
  }
  measurement_ids = get_measurement_ids(measurement_types)
  required_columns = c(diagnosis_col, covariates, measurement_ids)
  return(neuro_data[, ..required_columns])
}

# Calculating logistic regression models

# availableCores() automatically respects cgroups, SLURM, and local CPU limits
n_workers = parallelly::availableCores()
plan(multisession, workers = n_workers - 1)

# run a small task on all workers to force-load the libraries
# so serialization warnings don't trigger later.
invisible(future.apply::future_lapply(seq_len(n_workers), function(x) {
  library(future.apply)
  library(data.table)
  library(dplyr)
}, future.seed = TRUE))

# Bootstrap CIs for logistic regression coefficients (parallelized in chunks)
bootstrap_ci_logic = function(model, data, feature, diagnosis_col, covariate, bootstrap_iterations, pb = FALSE, cores_for_bootstrap = 12) {
  
  formula_str = paste(diagnosis_col, "~", feature, "+", paste(covariate, collapse = " + "))
  
  chunk_size = max(2, floor(bootstrap_iterations / cores_for_bootstrap))

  process_chunk = function(chunk_id, total_chunks, iterations, data, formula_str) {
    start_iteration = (chunk_id - 1) * iterations + 1
    end_iteration = min(chunk_id * iterations, bootstrap_iterations)
    resampled_log_odds = vector("numeric", end_iteration - start_iteration + 1)
    
    for (i in start_iteration:end_iteration) {
      resampled_data = data[sample(nrow(data), replace = TRUE), ]
      tryCatch({
        resampled_model = glm(as.formula(formula_str), data = resampled_data, family = "binomial")
        
        # This prevents NaN errors for negative coefficients.
        resampled_log_odds[i - start_iteration + 1] = coef(resampled_model)[feature]
        
      }, error = function(e) {
        resampled_log_odds[i - start_iteration + 1] = NA
      })
    }
    return(resampled_log_odds)
  }
  
  results = future_lapply(1:cores_for_bootstrap, function(chunk_id) {
    process_chunk(chunk_id, cores_for_bootstrap, chunk_size, data, formula_str)
  }, future.seed = TRUE, future.packages = c("future.apply"))

  bootstrapped_log_odds = unlist(results)
  valid_log_odds = bootstrapped_log_odds[!is.na(bootstrapped_log_odds)]
  
  if (length(valid_log_odds) > 0) {
    lower_log_ci = quantile(valid_log_odds, probs = 0.025)
    upper_log_ci = quantile(valid_log_odds, probs = 0.975)
    return(exp(c(lower_log_ci, upper_log_ci)))
  } else {
    warning("Bootstrap CI calculation failed: No valid models from resampling.")
    return(NA_real_)
  }
}

# Profile likelihood CIs for logistic regression (robust to imbalanced/small samples)
profile_likelihood_ci = function(model, feature, conf_level = 0.95) {
  require(MASS) # for confint function
  
  coef_index = which(names(coef(model)) == feature)
  log_ci = confint(model, parm = coef_index, level = conf_level)
  return(exp(log_ci))
}

# Stratified bootstrap CIs preserving case/control ratio in each resample
stratified_bootstrap_ci_logic = function(model, data, feature, diagnosis_col, covariate, bootstrap_iterations, pb = FALSE, cores_for_bootstrap = NULL) {
  formula_str = paste(diagnosis_col, "~", feature, "+", paste(covariate, collapse = " + "))
  
  # Default to sequential if inside an existing parallel worker or if not specified
  if (is.null(cores_for_bootstrap)) {
    cores_for_bootstrap = 1
  }
  
  # Define chunk size
  chunk_size = max(2, floor(bootstrap_iterations / cores_for_bootstrap))
  
  # Function to process each chunk with stratified sampling
  process_chunk = function(chunk_id, total_chunks, iterations, data, formula_str) {
    start_iteration = (chunk_id - 1) * iterations + 1
    end_iteration = min(chunk_id * iterations, bootstrap_iterations)
    resampled_coefs = vector("numeric", end_iteration - start_iteration + 1)
    
    # Split data once per chunk to save overhead
    cases = data[data[[diagnosis_col]] == 1, ]
    controls = data[data[[diagnosis_col]] == 0, ]
    
    for (i in start_iteration:end_iteration) {
      # Stratified sampling with replacement
      resampled_data = rbind(
        cases[sample(nrow(cases), replace = TRUE), ],
        controls[sample(nrow(controls), replace = TRUE), ]
      )
      
      tryCatch({
        resampled_model = glm(as.formula(formula_str), data = resampled_data, family = "binomial")
        # Extract coefficient for the feature (log-odds scale)
        resampled_coefs[i - start_iteration + 1] = coef(resampled_model)[feature]
      }, error = function(e) {
        resampled_coefs[i - start_iteration + 1] = NA
      })
    }
    return(resampled_coefs)
  }
  
  # Parallel processing of chunks using the global plan
  results = future_lapply(1:cores_for_bootstrap, function(chunk_id) {
    process_chunk(chunk_id, cores_for_bootstrap, chunk_size, data, formula_str)
  }, future.seed = TRUE, future.packages = c("future.apply"))
  
  # Combine and process the results
  bootstrapped_coefs = unlist(results)
  valid_coefs = bootstrapped_coefs[!is.na(bootstrapped_coefs)]
  
  if (length(valid_coefs) > 0) {
    # Calculate quantiles on the log-odds scale
    lower_ci = quantile(valid_coefs, probs = 0.025)
    upper_ci = quantile(valid_coefs, probs = 0.975)
    # Exponentiate to return Odds Ratios
    return(exp(c(lower_ci, upper_ci)))
  } else {
    warning("Stratified Bootstrap CI calculation failed: No valid models from resampling.")
    return(NA_real_)
  }
}

#' Run Logistic Regression and Calculate Odds Ratios and Confidence Intervals
#'
#' @param data A data frame containing the dataset for the logistic regression.
#' @param feature The name of the feature (predictor variable) to be used in the model.
#' @param diagnosis_col The name of the diagnosis column (response variable) in the data.
#' @param covariate Vector of covariate names to be included in the model.
#' @param bootstrap_iterations The number of bootstrap iterations for CIs 
#' @param show_progress Logical, indicating if progress should be shown.
#' @param total_features Total number of features being analyzed, used to adjust
#'        the number of cores for parallel processing.
#'         For univariate analysis (default scenario), keeping total_features at 1 is ideal. 
#'         For multivariate or batch analyses, adjust total_features to match the number of features processed together.
#' @return A list containing the logistic regression model, p-value of the feature, 
#'         and a data frame with the feature, odds ratio, and confidence intervals.
run_logistic_regression = function(data, feature, diagnosis_col, covariate, 
                                   bootstrap_iterations = NULL, ci_method = "bootstrap", 
                                   show_progress = FALSE, total_features = 1) {
  # Check if feature exists in the dataset
  if (!feature %in% names(data)) {
    stop(paste("Feature", feature, "not found in the data."))
  }
  
  # Create the formula for the logistic regression model
  formula_str = paste(diagnosis_col, "~", feature, "+", paste(covariate, collapse = " + "))
  model = glm(as.formula(formula_str), data = data, family = "binomial")
  
  # Check if feature exists in the model coefficients
  if (!feature %in% names(coef(model))) {
    stop(paste("Feature", feature, "not found in the model coefficients."))
  }
  
  # Extract p-values and odds ratios from the model
  p_values = summary(model)$coefficients[, "Pr(>|z|)"]
  ors = exp(coef(model))
  
  # Calculate confidence intervals based on the specified method
  if (ci_method == "bootstrap" && !is.null(bootstrap_iterations)) {
    # Dynamically determine the number of cores based on total features
    cores_for_bootstrap = max(1, floor(75 / total_features))
    cis = bootstrap_ci_logic(model, data, feature, diagnosis_col, 
                             covariate, bootstrap_iterations, 
                             pb = show_progress, cores_for_bootstrap)
  } else {
    # Profile likelihood CI calculation
    cis = profile_likelihood_ci(model, feature)
  }
  
  # Create a dataframe with the results
  ors_df = data.frame(feature = feature, odds_ratio = ors[feature], 
                      lower_CI = cis[1], upper_CI = cis[2])
  
  # Return the model, p-value, and odds ratios with confidence intervals
  return(list(model = model, p_value = p_values[feature], ors_df = ors_df))
}

#' Run Logistic Regression for All Features
#'
#' This function applies logistic regression to each feature in a dataset, 
#' adjusting for specified covariates. 
#' It calculates p-values, applies Bonferroni and FDR corrections for multiple testing, and 
#' calculates odds ratios (ORs) and confidence intervals (CIs) for each feature.
#' 
#' The function allows the choice between bootstrap and profile likelihood methods 
#' for calculating confidence intervals.
#'
#' @param neuro_data A data frame containing the neuroimaging data.
#' @param features A character vector of feature names to include in the models.
#' @param diagnosis_col A character string specifying the name of the diagnosis column.
#' @param covariates A character vector of covariate names to adjust for in the models.
#' @param bootstrap_iterations (Optional) Number of bootstrap iterations for CI calculation 
#'        if bootstrap method is used. Default is NULL.
#'        
#' @param ci_method Method for CI calculation: "bootstrap" or "profile_likelihood". 
#'        Default is "bootstrap".
#'        
#' @param conf_level Confidence level for the intervals (used in profile likelihood method). 
#'        Default is 0.95.
#'        
#' @param num_cores The number of cores to use for parallel processing. 
#'        If NULL, half of the available cores are used. Default is NULL.
#'        
#' @return A list with one elements: 'results_df' 
#'         'results_df' contains a data frame with one row per feature, 
#'         including p-values, corrected p-values, odds ratios 
#'         and  confidence intervals for each feature.
#' @export
run_all_features = function(neuro_data, features, diagnosis_col, covariates, 
                            bootstrap_iterations = NULL, ci_method = "bootstrap", 
                            conf_level = 0.95, num_cores = NULL) {
  
  # Parallel processing of logistic regression for each feature
  results = future_lapply(features, function(feature) {
    tryCatch({
      # Calculate total_features for dynamic core allocation inside run_logistic_regression
      total_features = length(features)
      
      run_logistic_regression(neuro_data, feature, diagnosis_col, covariates, 
                              bootstrap_iterations, ci_method, TRUE, 
                              total_features)
    }, error = function(e) {
      message(paste("Error processing feature", feature, ":", e$message))
      return(NULL)
    })
    # Explicitly load all relevant packages on the workers
  }, future.seed = TRUE, future.packages = c("dplyr", "data.table", "future.apply", "MASS"))
  
  # Filter out NULL and error results
  valid_results = Filter(function(x) !is.null(x) && !is.list(x$error), results)
  if (length(valid_results) == 0) {
    stop("All feature analyses failed.")
  }
  
  # Combine results into a single data frame
  results_df = do.call(dplyr::bind_rows, lapply(valid_results, `[[`, "ors_df"))
  
  # Add p-values and adjust for multiple comparisons
  p_values = unlist(lapply(valid_results, `[[`, "p_value"))
  results_df$p_value = p_values
  results_df$bonferroni_p_value = p.adjust(p_values, method = "bonferroni")
  results_df$fdr_p_value = p.adjust(p_values, method = "fdr")
  
  # Determine significance
  results_df$significant = p_values < 0.05
  results_df$significant_bonferroni = results_df$bonferroni_p_value < 0.05
  results_df$significant_fdr = results_df$fdr_p_value < 0.05
  
  return(list(results_df = results_df))
}

# Helper function to filter results for a subset of features
filter_results_for_subset = function(all_results, subset_features) {
  # Assuming all_results is a data frame or list with feature names as one of the columns/indices
  return(all_results[all_results$feature %in% subset_features, ])
}

#' Run Logistic Regression Analysis for a Single Dataset
#'
#' Performs logistic regression on a dataset, analyzing the overall feature set 
#' and subsets based on neuroimaging types. It supports choosing between bootstrap 
#' and profile likelihood methods for CI calculation.
#'
#' The function first runs logistic regression for all features, then filters 
#' results for specific measurement types.
#'
#' @param features_df A data frame representing features in the neuroimaging data.
#' @param covariates_df A data frame representing covariates.
#' @param outcome_df A data frame representing the outcome variable.
#' @param ci_method Method for CI calculation: "bootstrap" or "profile_likelihood". 
#'        Default is "bootstrap".
#' @param bootstrap_iterations (Optional) Number of bootstrap iterations for CI calculation 
#'        if bootstrap method is used. Default is NULL. Set to NULL if profile likelihood is used.
#' @param conf_level Confidence level for the intervals (used in profile likelihood method). 
#'        Default is 0.95.
#' @return A list containing logistic regression results, including a key "All_Features" 
#'         for overall results and additional keys for each measurement type 
#'         (e.g., "volume", "surface_area", "thickness", "FA"), containing respective 
#'         analysis results for those feature subsets.
#' @export
#' @examples
#' # Running analysis for a specific definition of depression using bootstrap method
#'   results = run_analysis_for_a_dataset(features_df, covariates_df, outcome_df, 
#'                                        "bootstrap", 500)
#'
#' # Running analysis using profile likelihood method
#'   results = run_analysis_for_a_dataset(features_df, covariates_df, outcome_df, 
#'                                        "profile_likelihood", NULL)
run_analysis_for_a_dataset = function(features_df, covariates_df, outcome_df, 
                                      ci_method = "bootstrap", 
                                      bootstrap_iterations = NULL, 
                                      conf_level = 0.95) {
  results = list()
  
  # Combine features and covariates into a single dataset for analysis
  full_data = as.data.frame(cbind(features_df, covariates_df, outcome_df))
  
  # Determine the name of the outcome column
  outcome_col_name = names(outcome_df)
  
  # Extract all feature names
  all_features = setdiff(names(features_df), outcome_col_name)
  
  # Check if bootstrap iterations are provided when the method is bootstrap
  if(ci_method == "bootstrap" && is.null(bootstrap_iterations)) {
    stop("Bootstrap iterations must be provided for bootstrap method.")
  }
  
  # Run logistic regression for all features
  all_features_results = run_all_features(full_data, all_features, 
                                          outcome_col_name, 
                                          names(covariates_df), 
                                          bootstrap_iterations, 
                                          ci_method, 
                                          conf_level)
  
  # Store results for all features
  results[["All_Features"]] = list(
    results_df = all_features_results$results_df
  )
  
  # Iterate over each measurement type and filter results
  for (measurement_type in c("global_volume", "subcortical_volume", "surface_area", "thickness", "FA", "MD")) {
    # Retrieve measurement IDs
    measurement_ids = get_measurement_ids(measurement_type)
    subset_features = intersect(all_features, measurement_ids)
    
    # Filter the all_features_results$results_df for the subset of features
    subset_results = filter_results_for_subset(all_features_results$results_df, 
                                               subset_features)
    
    # Store results for the subset of features in a structured manner
    results[[measurement_type]] = list(
      results_df = subset_results
    )
  }
  
  return(results)
}

# Running Logistic Regression across depression definitions

# GPNoDep
# Individuals who sought help for nerves, anxiety, tension or depression from a GP (fields 2090/2100), 
# but did not meet criteria for depression based on the touchscreen questionnaire (fields 4631/4598). 
LR_GPNoDep  = run_analysis_for_a_dataset(separated_GPNoDep$features, 
                                         separated_GPNoDep$covariates, 
                                         separated_GPNoDep$outcome, 
                                         ci_method = "profile_likelihood",
                                         conf_level = 0.95,
                                         bootstrap_iterations = NULL)
saveRDS(LR_GPNoDep, file = here("output", "S1_neuroimaging_objects", "LR_GPNoDep_PL_95.RDS"))

# GPpsy
# Individuals who sought help for nerves, anxiety, tension, or depression from a GP (fields 2090/2100).
LR_GPpsy               = run_analysis_for_a_dataset(separated_GPpsy$features, 
                                                    separated_GPpsy$covariates, 
                                                    separated_GPpsy$outcome, 
                                                    ci_method = "profile_likelihood",
                                                    conf_level = 0.95,
                                                    bootstrap_iterations = NULL)

saveRDS(LR_GPpsy, file = here("output", "S1_neuroimaging_objects", "LR_GPpsy_PL_95.RDS"))

# Psypsy 
# Individuals who sought help for nerves, anxiety, tension, or depression from a psychiatrist (fields 2090/2100).
LR_Psypsy              = run_analysis_for_a_dataset(separated_Psypsy$features, 
                                                    separated_Psypsy$covariates, 
                                                    separated_Psypsy$outcome, 
                                                    ci_method = "profile_likelihood",
                                                    conf_level = 0.95,
                                                    bootstrap_iterations = NULL)
saveRDS(LR_Psypsy, file = here("output", "S1_neuroimaging_objects", "LR_Psypsy_PL_95.RDS"))

# SelfRepDep
# Individuals who self-reported having depression (field 20002) during a verbal interview with a trained nurse. 
LR_SelfRepDep          = run_analysis_for_a_dataset(separated_SelfRepDep$features, 
                                                    separated_SelfRepDep$covariates, 
                                                    separated_SelfRepDep$outcome, 
                                                    ci_method = "profile_likelihood",
                                                    conf_level = 0.95,
                                                    bootstrap_iterations = NULL)
saveRDS(LR_SelfRepDep, file = here("output", "S1_neuroimaging_objects", "LR_SelfRepDep_PL_95.RDS"))

# DepAll
# Individuals who reported at least one of the two cardinal symptoms of depression (low mood or apathy, fields 4631/4598) 
# for at least two weeks on the Touchscreen questionnaire and sought help from either a GP or psychiatrist. 
# This provides a more comprehensive, though still minimal, phenotyping approach. 
LR_DepAll              = run_analysis_for_a_dataset(separated_DepAll$features, 
                                                    separated_DepAll$covariates, 
                                                    separated_DepAll$outcome, 
                                                    ci_method = "profile_likelihood",
                                                    conf_level = 0.95,
                                                    bootstrap_iterations = NULL)

saveRDS(LR_DepAll, file = here("output", "S1_neuroimaging_objects", "LR_DepAll_PL_95.RDS"))

# ICD10Dep
# This phenotype includes individuals with primary (data field 41202) or secondary (data field 41204) 
# ICD-10 diagnoses of mood or affective disorders (F32-F39) from linked hospital admission records in the UK Biobank. 
# The specific diagnoses include: depressive episodes (F32), recurrent depressive disorder (F33), persistent mood disorders (F34), 
# other mood disorders (F38), and unspecified mood disorders (F39). 
# This definition aligns with the "ICD10-coded depression" classification used in Howard et al. 2018: https://doi.org/10.1038/s41467-018-03819-3
LR_ICD10Dep            = run_analysis_for_a_dataset(separated_ICD10Dep$features, 
                                                    separated_ICD10Dep$covariates, 
                                                    separated_ICD10Dep$outcome, 
                                                    ci_method = "profile_likelihood",
                                                    conf_level = 0.95,
                                                    bootstrap_iterations = NULL)
saveRDS(LR_ICD10Dep, file = here("output", "S1_neuroimaging_objects", "LR_ICD10Dep_PL_95.RDS"))

# ICD10Dep_exclpsych
# This phenotype includes individuals with primary or secondary ICD-10 diagnoses of mood or affective disorders (F32-F39), 
# specifically: depressive episodes (F32), recurrent depressive disorder (F33), persistent mood disorders (F34), other mood disorders (F38), 
# and unspecified mood disorders (F39). 
# Both cases and controls exclude individuals with diagnoses of schizophrenia, bipolar disorder, and other psychiatric conditions. 
LR_ICD10Dep_exclpsych  = run_analysis_for_a_dataset(separated_ICD10Dep_exclpsych$features, 
                                                    separated_ICD10Dep_exclpsych$covariates, 
                                                    separated_ICD10Dep_exclpsych$outcome,
                                                    ci_method = "profile_likelihood",
                                                    conf_level = 0.95,
                                                    bootstrap_iterations = NULL)

saveRDS(LR_ICD10Dep_exclpsych, file = here("output", "S1_neuroimaging_objects", "LR_ICD10Dep_exclpsych_PL_95.RDS"))

# LifetimeMDD 
# Individuals who met the DSM-V criteria for lifetime Major Depressive Disorder (MDD) based on responses to the 
# Composite International Diagnostic Interview Short Form (CIDI-SF) in the UK Biobank’s follow-up online questionnaire. 
# This uses full DSM-5 criteria for diagnosing lifetime MDD, providing a stricter phenotyping definition. 
LR_LifetimeMDD         = run_analysis_for_a_dataset(separated_LifetimeMDD$features, 
                                                    separated_LifetimeMDD$covariates, 
                                                    separated_LifetimeMDD$outcome,
                                                    ci_method = "profile_likelihood",
                                                    conf_level = 0.95,
                                                    bootstrap_iterations = NULL)

saveRDS(LR_LifetimeMDD, file = here("output", "S1_neuroimaging_objects", "LR_LifetimeMDD_PL_95.RDS"))

# LR_MDDRecur
# Individuals from the LifetimeMDD group who reported having two or more depressive episodes (field 20442), indicating 
# recurrent MDD.
LR_MDDRecur            = run_analysis_for_a_dataset(separated_MDDRecurr$features, 
                                                    separated_MDDRecurr$covariates, 
                                                    separated_MDDRecurr$outcome, 
                                                    ci_method = "profile_likelihood",
                                                    conf_level = 0.95,
                                                    bootstrap_iterations = NULL)

saveRDS(LR_MDDRecur, file = here("output", "S1_neuroimaging_objects", "LR_MDDRecur_PL_95.RDS"))

# Functions to Calculate Overlap in Subject Numbers Across Labels

# Storing depression datasets as a list 
depression_datasets = list(
  GPNoDep            = z_score_completed_GPNoDep,
  Psypsy             = z_score_completed_Psypsy,
  SelfRepDep         = z_score_completed_SelfRepDep,
  GPpsy              = z_score_completed_GPpsy,
  DepAll             = z_score_completed_DepAll,
  ICD10Dep           = z_score_completed_ICD10Dep,
  ICD10Dep_exclpsych = z_score_completed_ICD10Dep_exclpsych,
  LifetimeMDD        = z_score_completed_LifetimeMDD,
  MDDRecur           = z_score_completed_MDDRecurr
)


#' Calculate Overlap of Subject Numbers Between Depression Labels
#'
#' @description
#' Calculates the percentage overlap of subject numbers between pairs of depression labels.
#'
#' @param depression_datasets A list of dataframes, each representing a different depression label with 'subject_num' and diagnosis columns.
#' @return A matrix or dataframe showing the percentage overlap of subject numbers between each pair of labels.
#'
#' @export
remove_overlapping_controls_with_summary = function(depression_datasets, diagnosis_column_names) {
  
  # Initialize a list for storing subjects to be removed from each dataset
  subject_removal_list = rep(list(), length(depression_datasets))
  
  # Prepare a data frame to store initial and final case/control counts for each dataset
  summary_data = data.frame(
    dataset_name = names(depression_datasets),
    initial_cases = integer(length(depression_datasets)),
    initial_controls = integer(length(depression_datasets)),
    final_cases = integer(length(depression_datasets)),
    final_controls = integer(length(depression_datasets)),
    stringsAsFactors = FALSE
  )
  
  # Capturing initial case and control counts for each dataset
  for (i in seq_along(depression_datasets)) {
    dataset = depression_datasets[[i]]
    diagnosis_col = diagnosis_column_names[[i]]
    summary_data$initial_cases[i] = sum(dataset[[diagnosis_col]] == 1)
    summary_data$initial_controls[i] = sum(dataset[[diagnosis_col]] == 0)
  }
  
  # Iterating through each dataset to identify overlapping controls
  for (i in seq_along(depression_datasets)) {
    dataset1 = depression_datasets[[i]]
    diagnosis_col1 = diagnosis_column_names[[i]]
    subject_removal_list[[i]] = integer(0)  # Start with an empty list for the i-th dataset
    
    for (j in seq_along(depression_datasets)) {
      if (i != j) {
        dataset2 = depression_datasets[[j]]
        diagnosis_col2 = diagnosis_column_names[[j]]
        
        # Identifying subjects that are controls in dataset1 but cases in dataset2
        controls_in_dataset1 = dataset1$subject_num[dataset1[[diagnosis_col1]] == 0]
        cases_in_dataset2 = dataset2$subject_num[dataset2[[diagnosis_col2]] == 1]
        overlapping_subjects = intersect(controls_in_dataset1, cases_in_dataset2)
        
        # Append the overlapping subjects to the removal list for the i-th dataset
        subject_removal_list[[i]] = unique(c(subject_removal_list[[i]], overlapping_subjects))
      }
    }
    
    # Removing the identified overlapping subjects from the i-th dataset
    dataset = depression_datasets[[i]]
    subjects_to_remove = subject_removal_list[[i]]
    cleaned_dataset = dataset[!dataset$subject_num %in% subjects_to_remove, ]
    depression_datasets[[i]] = cleaned_dataset
    
    # Updating summary data with the counts after removal
    summary_data$final_cases[i]    = sum(cleaned_dataset[[diagnosis_col1]] == 1)
    summary_data$final_controls[i] = sum(cleaned_dataset[[diagnosis_col1]] == 0)
  }
  
  # Return the cleaned datasets and the summary data
  return(list(cleaned_datasets = depression_datasets, summary = summary_data))
}

diagnosis_column_names = c("GPNoDep", "Psypsy", "SelfRepDep", "GPpsy", "DepAll", "ICD10Dep", "ICD10Dep.exclpsych", "LifetimeMDD", "MDDRecur")
non_overlap_controls_df  = remove_overlapping_controls_with_summary(depression_datasets, diagnosis_column_names)

save(non_overlap_controls_df, file = here("output", "S1_neuroimaging_objects", "non_overlapping_control_df.RData"))

# Separating features and covariates for logistic regression analysis
# Each dataset corresponds to a different depression definition. 
# The 'separate_features_and_covariates' function isolates features, covariates, and outcomes
# for each depression dataset, preparing them for logistic regression analysis.
GPNoDep_nonoverlap_controls            = separate_features_and_covariates(non_overlap_controls_df$cleaned_datasets$GPNoDep, covariates, "GPNoDep")
Psypsy_nonoverlap_controls             = separate_features_and_covariates(non_overlap_controls_df$cleaned_datasets$Psypsy, covariates, "Psypsy")
SelfRepDep_nonoverlap_controls         = separate_features_and_covariates(non_overlap_controls_df$cleaned_datasets$SelfRepDep, covariates, "SelfRepDep")
GPpsy_nonoverlap_controls              = separate_features_and_covariates(non_overlap_controls_df$cleaned_datasets$GPpsy, covariates, "GPpsy")
ICD10Dep_nonoverlap_controls           = separate_features_and_covariates(non_overlap_controls_df$cleaned_datasets$ICD10Dep, covariates, "ICD10Dep")
ICD10Dep_nonoverlap_exclpsych_controls = separate_features_and_covariates(non_overlap_controls_df$cleaned_datasets$ICD10Dep_exclpsych, covariates, "ICD10Dep.exclpsych")
DepAll_nonoverlap_controls             = separate_features_and_covariates(non_overlap_controls_df$cleaned_datasets$DepAll, covariates, "DepAll")
LifetimeMDD_nonoverlap_controls        = separate_features_and_covariates(non_overlap_controls_df$cleaned_datasets$LifetimeMDD, covariates, "LifetimeMDD")
MDDRecur_nonoverlap_controls           = separate_features_and_covariates(non_overlap_controls_df$cleaned_datasets$MDDRecur, covariates, "MDDRecur")


# Function to subset features for both surface area and thickness, and save to CSV
save_combined_features_and_outcome = function(dataset, surface_area_ids, thickness_ids, outcome_name) {
  tryCatch({
    # Prefix the feature IDs with "bl_" and combine them
    combined_feature_ids = c(paste0("bl_", surface_area_ids), paste0("bl_", thickness_ids))
    
    # Subset the combined features
    if (!all(combined_feature_ids %in% names(dataset$features))) {
      stop("Not all feature IDs are found in the dataset.")
    }
    combined_feature_subset = dataset$features[, ..combined_feature_ids, with = FALSE]
    
    # Check if the outcome column exists
    if (!outcome_name %in% names(dataset$outcome)) {
      stop("Outcome column not found in the dataset.")
    }
    
    # Bind the outcome to the combined feature subset
    combined_data = cbind(combined_feature_subset, dataset$outcome[, ..outcome_name, with = FALSE])
    
    # Save to CSV
    file_path = paste0(outcome_name, "_surface_area_thickness.csv")
    fwrite(combined_data, file = file_path, row.names = FALSE)
    cat("Saved dataset to", file_path, "\n")
  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
  })
}

# Example usage for one of the datasets
save_combined_features_and_outcome(GPNoDep_nonoverlap_controls, surface_area_ids, thickness_ids, "GPNoDep")
save_combined_features_and_outcome(Psypsy_nonoverlap_controls, surface_area_ids, thickness_ids, "Psypsy")
save_combined_features_and_outcome(SelfRepDep_nonoverlap_controls, surface_area_ids, thickness_ids, "SelfRepDep")
save_combined_features_and_outcome(GPpsy_nonoverlap_controls, surface_area_ids, thickness_ids, "GPpsy")
save_combined_features_and_outcome(ICD10Dep_nonoverlap_controls, surface_area_ids, thickness_ids, "ICD10Dep")
save_combined_features_and_outcome(ICD10Dep_nonoverlap_exclpsych_controls, surface_area_ids, thickness_ids, "ICD10Dep.exclpsych")
save_combined_features_and_outcome(DepAll_nonoverlap_controls, surface_area_ids, thickness_ids, "DepAll")
save_combined_features_and_outcome(LifetimeMDD_nonoverlap_controls, surface_area_ids, thickness_ids, "LifetimeMDD")
save_combined_features_and_outcome(MDDRecur_nonoverlap_controls, surface_area_ids, thickness_ids, "MDDRecur")

# Running Logistic Regression with Non-overlapping controls
# This section conducts logistic regression analyses for each depression definition
# using the constant control groups. 
# It applies the 'profile_likelihood' method for CI calculation.
# This method is chosen for its robustness in smaller samples and/or imbalanced datasets 
# Note that run_analysis_for_a_dataset can also be supplied with ci_method = "bootstrapping" 
LR_GPNoDep_nonoverlap_controls   = run_analysis_for_a_dataset(GPNoDep_nonoverlap_controls$features, 
                                                              GPNoDep_nonoverlap_controls$covariates, 
                                                              GPNoDep_nonoverlap_controls$outcome, 
                                                              ci_method = "profile_likelihood",
                                                              conf_level = 0.95,
                                                              bootstrap_iterations = NULL)

LR_Psypsy_nonoverlap_controls   = run_analysis_for_a_dataset(Psypsy_nonoverlap_controls$features, 
                                                             Psypsy_nonoverlap_controls$covariates, 
                                                             Psypsy_nonoverlap_controls$outcome, 
                                                             ci_method = "profile_likelihood",
                                                             conf_level = 0.95,
                                                             bootstrap_iterations = NULL)

LR_SelfRepDep_nonoverlap_controls   = run_analysis_for_a_dataset(SelfRepDep_nonoverlap_controls$features, 
                                                                 SelfRepDep_nonoverlap_controls$covariates, 
                                                                 SelfRepDep_nonoverlap_controls$outcome, 
                                                                 ci_method = "profile_likelihood",
                                                                 conf_level = 0.95,
                                                                 bootstrap_iterations = NULL)

LR_GPpsy_nonoverlap_controls   = run_analysis_for_a_dataset(GPpsy_nonoverlap_controls$features, 
                                                            GPpsy_nonoverlap_controls$covariates, 
                                                            GPpsy_nonoverlap_controls$outcome, 
                                                            ci_method = "profile_likelihood",
                                                            conf_level = 0.95,
                                                            bootstrap_iterations = NULL)

LR_ICD10Dep_nonoverlap_controls   = run_analysis_for_a_dataset(ICD10Dep_nonoverlap_controls$features, 
                                                               ICD10Dep_nonoverlap_controls$covariates, 
                                                               ICD10Dep_nonoverlap_controls$outcome, 
                                                               ci_method = "profile_likelihood",
                                                               conf_level = 0.95,
                                                               bootstrap_iterations = NULL)

LR_ICD10Dep_exclpsych_nonoverlap_controls  = run_analysis_for_a_dataset(ICD10Dep_nonoverlap_exclpsych_controls$features, 
                                                                        ICD10Dep_nonoverlap_exclpsych_controls$covariates, 
                                                                        ICD10Dep_nonoverlap_exclpsych_controls$outcome, 
                                                                        ci_method = "profile_likelihood",
                                                                        conf_level = 0.95,
                                                                        bootstrap_iterations = NULL)

LR_DepAll_nonoverlap_controls  = run_analysis_for_a_dataset(DepAll_nonoverlap_controls$features, 
                                                            DepAll_nonoverlap_controls$covariates, 
                                                            DepAll_nonoverlap_controls$outcome, 
                                                            ci_method = "profile_likelihood",
                                                            conf_level = 0.95,
                                                            bootstrap_iterations = NULL)

LR_LifetimeMDD_nonoverlap_controls  = run_analysis_for_a_dataset(LifetimeMDD_nonoverlap_controls$features, 
                                                                 LifetimeMDD_nonoverlap_controls$covariates, 
                                                                 LifetimeMDD_nonoverlap_controls$outcome, 
                                                                 ci_method = "profile_likelihood",
                                                                 conf_level = 0.95,
                                                                 bootstrap_iterations = NULL)

LR_MDDRecur_nonoverlap_controls  = run_analysis_for_a_dataset(MDDRecur_nonoverlap_controls$features, 
                                                              MDDRecur_nonoverlap_controls$covariates, 
                                                              MDDRecur_nonoverlap_controls$outcome, 
                                                              ci_method = "profile_likelihood",
                                                              conf_level = 0.95,
                                                              bootstrap_iterations = NULL)

# Summary statistics for the LR results 

# Create a list of all logistic regression analysis objects
lr_objects = list(
  GPNoDep            = LR_GPNoDep_nonoverlap_controls,
  GPpsy              = LR_GPpsy_nonoverlap_controls,
  Psypsy             = LR_Psypsy_nonoverlap_controls,
  SelfRepDep         = LR_SelfRepDep_nonoverlap_controls,
  DepAll             = LR_DepAll_nonoverlap_controls,
  ICD10Dep           = LR_ICD10Dep_nonoverlap_controls,
  ICD10Dep_exclpsych = LR_ICD10Dep_exclpsych_nonoverlap_controls,
  LifetimeMDD        = LR_LifetimeMDD_nonoverlap_controls,
  MDDRecur           = LR_MDDRecur_nonoverlap_controls
)

# Function to summarize statistics from a given results data frame
summarize_statistics = function(results_df) {
  list(
    mean_or = mean(results_df$odds_ratio),      # Calculate mean of odds ratios
    median_or = median(results_df$odds_ratio),  # Calculate median of odds ratios
    min_or    = min(results_df$odds_ratio),     # Find minimum odds ratio
    max_or    = max(results_df$odds_ratio),     # Find maximum odds ratio
    sd_or     = sd(results_df$odds_ratio),      # Calculate standard deviation of odds ratios
    significant_count  = sum(results_df$significant),   # Count number of significant results
    mean_fdr_p_value   = mean(results_df$fdr_p_value),  # Calculate mean of FDR p-values
    median_fdr_p_value = median(results_df$fdr_p_value) # Calculate median of FDR p-values
  )
}

# Function to process each logistic regression object and summarize its statistics
process_lr_object = function(lr_object) {
  lapply(lr_object, function(category) summarize_statistics(category$results_df))
}

# Apply the function to each logistic regression object and store the results
all_summaries = lapply(lr_objects, process_lr_object)
print(all_summaries)

# Function to format the summaries into a data frame
format_summaries = function(summaries, condition) {
  do.call(rbind, lapply(summaries, function(cat_sum) {
    data.frame(
      Condition = condition,  # Condition for the summary
      Category  = names(cat_sum),  # Category of the summary
      Mean_OR   = cat_sum$mean_or,  # Mean odds ratio
      Median_OR = cat_sum$median_or,  # Median odds ratio
      Min_OR    = cat_sum$min_or,  # Minimum odds ratio
      Max_OR    = cat_sum$max_or,  # Maximum odds ratio
      SD_OR     = cat_sum$sd_or,  # Standard deviation of odds ratios
      Significant_Count  = cat_sum$significant_count,  # Count of significant results
      Mean_FDR_P_Value   = cat_sum$mean_fdr_p_value,  # Mean FDR p-value
      Median_FDR_P_Value = cat_sum$median_fdr_p_value  # Median FDR p-value
    )
  }))
}

# Apply the formatting function to each set of summaries and combine into a single data frame
final_summaries = do.call(rbind, mapply(format_summaries, all_summaries, names(all_summaries), SIMPLIFY = FALSE))

# View the final combined summaries
print(final_summaries)

save(final_summaries, file = here("output", "S1_neuroimaging_objects", "summary_stats_effect_sizes_across_phenotypes.RData"))

# Function to save all objects with a specific suffix to an RData file
save_objects_with_suffix = function(suffix, save_path) {
  # Retrieve all objects in the environment with the specified suffix
  objects_to_save = ls(envir = .GlobalEnv, pattern = paste0(suffix, "$"))
  
  # Check if there are objects to save
  if (length(objects_to_save) == 0) {
    warning("No objects found with the specified suffix.")
    return(invisible())
  }
  
  # Create a list of objects to save
  objects_list = mget(objects_to_save, envir = .GlobalEnv)
  
  # Save the list of objects to an RData file
  save(list = objects_to_save, file = save_path, envir = .GlobalEnv)
  
  # Provide feedback
  message("Saved ", length(objects_to_save), " objects with suffix '", suffix, "' to '", save_path, "'.")
}

# Identifying save path and calling the function
save_path = "/dw459/fuss/aim1_2024/molecular_psychiatry_scripts/preprocessing_objects/S1_neuroimaging_objects/sensitivity_analysis/sensitivity_nonoverlap_controls_objects.RData"
save_objects_with_suffix("_nonoverlap_controls", save_path)

# Plotting significant features across depression definitions 
# Using lr_objects list as input

# Define a function to create plots for each feature type
create_feature_plot = function(lr_objects, feature_type) {
  # Initialize an empty data frame to store the counts
  feature_counts = data.frame(Definition = character(), Count = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each depression definition and count the significant features
  for (definition_name in names(lr_objects)) {
    definition_data = lr_objects[[definition_name]][[feature_type]]$results_df
    significant_count = sum(definition_data$significant_fdr, na.rm = TRUE)
    feature_counts = rbind(feature_counts, data.frame(Definition = definition_name, Count = significant_count))
  }
  
  # Create the plot
  ggplot(feature_counts, aes(x = Definition, y = Count, fill = Definition)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.5, size = 3.5) + # Adjust 'vjust' for spacing and 'size' for text size
    theme_minimal() +
    labs(title = paste("Count of Significant Features for", feature_type),
         x = "Definition of Depression",
         y = "Count of Significant Features") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text(face = "bold"),
          legend.position = "none")  # Remove legend if not necessary
}

# Feature types to iterate over
feature_types = c("All_Features", "global_volume", "subcortical_volume", "surface_area", "thickness", "FA", "MD")

# List to store all the plots
all_plots = list()

# Iterate over feature types and create plots
for (feature_type in feature_types) {
  all_plots[[feature_type]] = create_feature_plot(lr_objects, feature_type)
}

# Saving all_plots object
save(all_plots, 
     file = here("output", "results_plots.RData"))

# Access the plots using all_plots$feature_type, e.g., all_plots$All_Features
