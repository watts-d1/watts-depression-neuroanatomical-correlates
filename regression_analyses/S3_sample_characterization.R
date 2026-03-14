#!/usr/bin/env Rscript
# S3_sample_characterization.R
# UKB Sample Characterization Table for Molecular Psychiatry 
#
# Purpose: For each of 9 depression phenotypes, extract demographics, severity
# proxies, chronicity, treatment exposure, comorbidities, and MHQ completion
# rates for cases and controls separately.
#
# Dependencies: data.table, MASS (for confint)
# Input:  csv_data (UKB phenotype file), z_score_completed_* objects
# Output: Formatted characterization table (CSV + console output)
#
# Design notes:
# 1. MHQ timing exclusion: Participants whose MRI preceded MHQ are excluded
#    from ALL MHQ-derived summary statistics (Section 4.2 of context file).
# 2. UKB negative sentinel values (-1, -3, -818) are recoded to NA, NEVER
#    treated as numeric data.
# 3. Illness duration and episode count are reported for CASES ONLY.
# Author: Devon Watts

# --- User-defined paths (modify for your environment) ---
UKBB_DATA_PATH   = "/path/to/ukbiobank/ukb674571.csv"
ZSCORE_DATA_PATH = "/path/to/preprocessing_objects/z_score_neuroimaging_objects.RData"

library(data.table)
library(here)

# 0. CONFIGURATION

# --- File paths ---
csv_file_path     = UKBB_DATA_PATH
z_score_path      = ZSCORE_DATA_PATH
output_dir        = here("output", "characterization_results")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Antidepressant coding4 codes (from UKB coding dictionary) ---
# These are the numeric codes stored in field 20003 for the target medications.
# Source: coding4.tsv matched against drug names specified in data dictionary.
# https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=4
antidepressant_codes = c(
  1140879540,  # fluoxetine
  1140867876,  # prozac 20mg capsule
  1140867878,  # sertraline
  1141180212,  # escitalopram
  1141200564,  # duloxetine
  1141201834,  # cymbalta 30mg gastro-resistant capsule
  1140921600,  # citalopram
  1140916282,  # venlafaxine
  1141176854,  # bupropion
  1141152732   # mirtazapine
)

# 1. LOAD DATA

cat("Loading z_score objects...\n")
load(z_score_path)

cat("Loading UKB CSV (this may take several minutes)...\n")
csv_data = fread(csv_file_path, showProgress = TRUE, nThread = 20)

# 2. PHENOTYPE REGISTRY
# Maps phenotype names to their z_score object and outcome column name.
# Case = 1, Control = 0 in the outcome column.

phenotype_registry = list(
  GPNoDep            = list(data = z_score_completed_GPNoDep,            col = "GPNoDep",            label = "GP (No Depression Screen)"),
  GPPsy              = list(data = z_score_completed_GPpsy,              col = "GPpsy",              label = "GP Help-Seeking"),
  PsyPsy             = list(data = z_score_completed_Psypsy,             col = "Psypsy",             label = "Psychiatrist Help-Seeking"),
  SelfRepDep         = list(data = z_score_completed_SelfRepDep,         col = "SelfRepDep",         label = "Self-Reported Depression"),
  DepAll             = list(data = z_score_completed_DepAll,             col = "DepAll",             label = "Cardinal Symptom + Help-Seeking"),
  ICD10Dep           = list(data = z_score_completed_ICD10Dep,           col = "ICD10Dep",           label = "ICD-10 Depression"),
  ICD10Dep_exclpsych = list(data = z_score_completed_ICD10Dep_exclpsych, col = "ICD10Dep.exclpsych", label = "ICD-10 (Excl. Psych Comorbidity)"),
  LifetimeMDD        = list(data = z_score_completed_LifetimeMDD,        col = "LifetimeMDD",        label = "Lifetime MDD (CIDI-SF)"),
  MDDRecur           = list(data = z_score_completed_MDDRecurr,          col = "MDDRecur",           label = "Recurrent MDD")
)

# Expected sample sizes from manuscript (Table in context file Section 3)
expected_n = list(
  GPNoDep            = list(cases = 683,  controls = 4789),
  GPPsy              = list(cases = 9529, controls = 20477),
  PsyPsy             = list(cases = 2706, controls = 27353),
  SelfRepDep         = list(cases = 1944, controls = 20421),
  DepAll             = list(cases = 1994, controls = 5491),
  ICD10Dep           = list(cases = 432,  controls = 15812),
  ICD10Dep_exclpsych = list(cases = 392,  controls = 15432),
  LifetimeMDD        = list(cases = 3019, controls = 9458),
  MDDRecur           = list(cases = 1869, controls = 9105)
)

# 3. VERIFICATION: CHECK CASE/CONTROL COUNTS

cat("\n=== Verification: Case/Control Counts ===\n")
for (pheno_name in names(phenotype_registry)) {
  reg = phenotype_registry[[pheno_name]]
  dt  = as.data.table(reg$data)
  
  obs_cases    = sum(dt[[reg$col]] == 1, na.rm = TRUE)
  obs_controls = sum(dt[[reg$col]] == 0, na.rm = TRUE)
  exp          = expected_n[[pheno_name]]
  
  case_match    = obs_cases == exp$cases
  control_match = obs_controls == exp$controls
  
  status = ifelse(case_match & control_match, "OK", "MISMATCH")
  cat(sprintf("  %s: Cases=%d (expected %d, %s), Controls=%d (expected %d, %s) [%s]\n",
              pheno_name, obs_cases, exp$cases, ifelse(case_match, "match", "MISMATCH"),
              obs_controls, exp$controls, ifelse(control_match, "match", "MISMATCH"),
              status))
  
  if (!case_match | !control_match) {
    cat(sprintf("    [ERROR] Sample size mismatch for %s. Verify phenotype object.\n", pheno_name))
  }
}

# 4. EXTRACT SUBJECT IDS PER PHENOTYPE
# We need subject_num from each z_score object to merge with csv_data.

get_subject_ids = function(pheno_name) {
  reg = phenotype_registry[[pheno_name]]
  dt  = as.data.table(reg$data)
  
  cases    = dt[dt[[reg$col]] == 1, .(subject_num)]
  controls = dt[dt[[reg$col]] == 0, .(subject_num)]
  
  return(list(cases = cases$subject_num, controls = controls$subject_num))
}

# 5. MHQ TIMING EXCLUSION
# Exclude participants whose brain MRI (field 53-2.0) preceded MHQ completion
# (field 20400-0.0). This exclusion applies to ALL MHQ-derived variables.

cat("\n=== MHQ Timing Exclusion ===\n")

# Identify the column names in csv_data
mhq_date_col     = "20400-0.0"
imaging_date_col  = "53-2.0"

# Determine subject ID column
id_col = ifelse("subject_num" %in% names(csv_data), "subject_num", "eid")

# Build the MHQ timing exclusion set
mhq_timing = csv_data[, c(id_col, mhq_date_col, imaging_date_col), with = FALSE]
setnames(mhq_timing, c("subject_num", "mhq_date", "imaging_date"))
mhq_timing[, subject_num := as.integer(subject_num)]

# Convert to Date if stored as character
if (is.character(mhq_timing$mhq_date)) {
  mhq_timing[, mhq_date := as.Date(mhq_date)]
}
if (is.character(mhq_timing$imaging_date)) {
  mhq_timing[, imaging_date := as.Date(imaging_date)]
}

# Participants to EXCLUDE: imaging occurred before MHQ
mri_before_mhq = mhq_timing[!is.na(mhq_date) & !is.na(imaging_date) & imaging_date < mhq_date, subject_num]
# Participants with valid MHQ timing (MHQ before MRI)
mhq_valid      = mhq_timing[!is.na(mhq_date) & !is.na(imaging_date) & mhq_date <= imaging_date, subject_num]

cat(sprintf("  MHQ before MRI (valid): n = %d\n", length(mhq_valid)))
cat(sprintf("  MRI before MHQ (exclude): n = %d\n", length(mri_before_mhq)))

# Participants who completed MHQ at all (regardless of timing)
mhq_completers = mhq_timing[!is.na(mhq_date), subject_num]
cat(sprintf("  Total MHQ completers: n = %d\n", length(mhq_completers)))

# 6. EXTRACT NEW VARIABLES FROM csv_data

cat("\n=== Extracting characterization variables from csv_data ===\n")

# --- 6a. Demographics ---
# Age and sex in the z_score objects (bl_54, bl_31) are z-score normalized
# and cannot be used for descriptive statistics. Extract raw values from csv_data.
# Age at imaging visit: field 21003-2.0 (instance 2 = imaging visit)
# Sex: field 31-0.0 (0 = Female, 1 = Male)
age_raw_field = "21003-2.0"
sex_raw_field = "31-0.0"

# --- 6b. Medication (field 20003, all array columns at instance 0) ---

# --- Antidepressant extraction from field 20003 (all array columns) ---
# Field 20003 stores medication codes across multiple array indices within
# instance 0 (20003-0.0 through 20003-0.N). Each array index represents a
# different medication the participant reported, NOT a different visit.
# We search ALL array columns to avoid undercounting participants who listed
# an antidepressant as their 2nd, 3rd, etc. medication.

extract_antidepressant = function(csv_data, id_col) {
  # Find all columns matching 20003-0.* (instance 0, all array indices)
  med_cols = grep("^20003-0\\.", names(csv_data), value = TRUE)
  cat(sprintf("  Found %d medication columns for field 20003 instance 0\n", length(med_cols)))
  
  if (length(med_cols) == 0) {
    stop("No columns matching '20003-0.*' found in csv_data. Check column naming.")
  }
  
  med_dt = csv_data[, c(id_col, med_cols), with = FALSE]
  setnames(med_dt, old = id_col, new = "subject_num")
  med_dt[, subject_num := as.integer(subject_num)]
  
  # Check each column for antidepressant codes
  med_dt[, any_antidepressant := 0L]
  for (col in med_cols) {
    med_dt[get(col) %in% antidepressant_codes, any_antidepressant := 1L]
  }
  
  return(med_dt[, .(subject_num, any_antidepressant)])
}

medication_dt = extract_antidepressant(csv_data, id_col)

cat(sprintf("  Antidepressant prevalence in full UKB imaging sample: %d / %d (%.1f%%)\n",
            sum(medication_dt$any_antidepressant == 1, na.rm = TRUE),
            nrow(medication_dt),
            100 * mean(medication_dt$any_antidepressant == 1, na.rm = TRUE)))

# --- 6c. MHQ-derived severity proxies ---
# PHQ-9 style items: 20510, 20514, 20507, 20508, 20519, 20513, 20517, 20511, 20518
phq_fields = c("20510-0.0", "20514-0.0", "20507-0.0", "20508-0.0",
               "20519-0.0", "20513-0.0", "20517-0.0", "20511-0.0", "20518-0.0")

# Suicidal ideation item (field 20513) also serves as a standalone severity indicator
suicidal_ideation_field = "20513-0.0"

# Duration of worst episode (field 20438)
worst_episode_duration_field = "20438-0.0"

# --- 6d. Chronicity proxies (CASES ONLY) ---
age_first_episode_field  = "20433-0.0"
num_episodes_field       = "20442-0.0"

# --- 6e. Comorbidity proxies (MHQ-based, per author preference) ---
anxiety_proxy_field      = "20428-0.0"    # Professional informed about anxiety
alcohol_addiction_field   = "20406-0.0"    # Ever addicted to alcohol
personality_disorder_field = "130932-0.0" # Date F60 first reported
bipolar_field            = "20126-0.0"    # Bipolar and major depression status

# --- 6f. Demographics (raw, from csv_data; z_score objects have normalized values) ---

# Build a single characterization data.table from csv_data
all_fields = c(age_raw_field, sex_raw_field,
               phq_fields, suicidal_ideation_field, worst_episode_duration_field,
               age_first_episode_field, num_episodes_field,
               anxiety_proxy_field, alcohol_addiction_field,
               personality_disorder_field, bipolar_field)

# Remove duplicates (suicidal_ideation_field is also in phq_fields)
all_fields = unique(all_fields)

# Check which fields exist in csv_data
missing_fields = setdiff(all_fields, names(csv_data))
if (length(missing_fields) > 0) {
  cat(sprintf("  WARNING: Fields not found in csv_data: %s\n", paste(missing_fields, collapse = ", ")))
  cat("  Check column naming convention. UKB may use different formats.\n")
  # Remove missing fields from extraction
  all_fields = intersect(all_fields, names(csv_data))
}

char_dt = csv_data[, c(id_col, all_fields), with = FALSE]
setnames(char_dt, old = id_col, new = "subject_num")
char_dt[, subject_num := as.integer(subject_num)]

# 7. RECODE VARIABLES

cat("\n=== Recoding variables ===\n")

# --- 7a. Raw demographics ---
# Age at imaging (field 21003-2.0): continuous, no sentinel values expected
if (age_raw_field %in% names(char_dt)) {
  char_dt[, age_at_imaging := as.numeric(get(age_raw_field))]
  cat(sprintf("  Age at imaging: n valid = %d, range = %.0f-%.0f\n",
              sum(!is.na(char_dt$age_at_imaging)),
              min(char_dt$age_at_imaging, na.rm = TRUE),
              max(char_dt$age_at_imaging, na.rm = TRUE)))
} else {
  cat(sprintf("  [AUTHOR NOTE] Field %s not found. Check column name for age at imaging.\n", age_raw_field))
}

# Sex (field 31-0.0): 0 = Female, 1 = Male
if (sex_raw_field %in% names(char_dt)) {
  char_dt[, sex_raw := as.integer(get(sex_raw_field))]
  cat(sprintf("  Sex: n Female = %d, n Male = %d\n",
              sum(char_dt$sex_raw == 0, na.rm = TRUE),
              sum(char_dt$sex_raw == 1, na.rm = TRUE)))
} else {
  cat(sprintf("  [AUTHOR NOTE] Field %s not found. Check column name for sex.\n", sex_raw_field))
}

# --- 7b. PHQ-9 sum score ---
# UKB coding: 1 = Not at all, 2 = Several days, 3 = More than half the days,
#             4 = Nearly every day, -818 = Prefer not to answer
# Standard PHQ-9: 0 = Not at all, 1 = Several days, 2 = More than half, 3 = Nearly every day
# Recode: subtract 1 from UKB values; set -818 to NA.

phq_field_names = phq_fields  # column names in char_dt

for (field in phq_field_names) {
  if (field %in% names(char_dt)) {
    char_dt[, (field) := {
      x = get(field)
      x = as.numeric(x)
      x[x < 0] = NA  # Remove all negative sentinel values (-818, etc.)
      x = x - 1      # Recode: UKB 1->0, 2->1, 3->2, 4->3
      # Sanity check: values should now be 0, 1, 2, 3 or NA
      x[!x %in% c(0, 1, 2, 3)] = NA
      x
    }]
  }
}

# Compute PHQ-9 mean score (range 0-3).
# Each item is recoded 0-3, then we take the mean of all non-missing items.
# This produces a score on the same 0-3 scale regardless of whether a
# participant answered 7, 8, or 9 items, avoiding listwise deletion for
# participants missing 1-2 items. Participants with fewer than 7 valid
# items are set to NA to avoid unstable estimates from sparse data.
if (all(phq_field_names %in% names(char_dt))) {
  char_dt[, phq9_n_valid := rowSums(!is.na(.SD)), .SDcols = phq_field_names]
  char_dt[, phq9_mean := ifelse(phq9_n_valid >= 7,
                                rowMeans(.SD, na.rm = TRUE),
                                NA_real_), .SDcols = phq_field_names]
  cat(sprintf("  PHQ-9 mean score (0-3) computed: n valid = %d (of which %d had 7-8 items, %d had all 9), n excluded (<7 items) = %d\n",
              sum(!is.na(char_dt$phq9_mean)),
              sum(char_dt$phq9_n_valid %in% c(7, 8), na.rm = TRUE),
              sum(char_dt$phq9_n_valid == 9, na.rm = TRUE),
              sum(is.na(char_dt$phq9_mean))))
} else {
  cat("  Error: Not all PHQ fields found. PHQ-9 mean not computed.\n")
}

# --- 7b. Suicidal ideation (binary) ---
# Already recoded as part of PHQ items above (0 = not at all, 1+ = any endorsement)
# Create a separate binary indicator: any endorsement (>=1) vs. none (0)
si_col = "20513-0.0"
if (si_col %in% names(char_dt)) {
  char_dt[, suicidal_ideation_binary := as.integer(get(si_col) >= 1)]
}

# --- 7c. Duration of worst episode (field 20438) ---
# UKB coding: continuous (months), with negative sentinels
worst_dur_col = worst_episode_duration_field
if (worst_dur_col %in% names(char_dt)) {
  char_dt[, (worst_dur_col) := {
    x = as.numeric(get(worst_dur_col))
    x[x < 0] = NA  # Remove sentinels
    x
  }]
}

# --- 7d. Age at first episode (field 20433) ---
# Negative values are sentinels (-1 = Do not know, -3 = Prefer not to answer)
age_first_col = age_first_episode_field
if (age_first_col %in% names(char_dt)) {
  char_dt[, (age_first_col) := {
    x = as.numeric(get(age_first_col))
    x[x < 0] = NA  # Remove sentinels
    x
  }]
}

# --- 7e. Number of episodes (field 20442) ---
# Negative values are sentinels. Positive integers only.
num_ep_col = num_episodes_field
if (num_ep_col %in% names(char_dt)) {
  char_dt[, (num_ep_col) := {
    x = as.numeric(get(num_ep_col))
    x[x < 0] = NA  # Remove sentinels; do NOT code as 0
    x
  }]
}

# --- 7f. Anxiety proxy (field 20428) ---
# 1 = Yes, 0 = No, negative = sentinel
anx_col = anxiety_proxy_field
if (anx_col %in% names(char_dt)) {
  char_dt[, (anx_col) := {
    x = as.numeric(get(anx_col))
    x[x < 0] = NA
    as.integer(x == 1)
  }]
}

# --- 7g. Alcohol addiction (field 20406) ---
# 1 = Yes, 0 = No, negative = sentinel
alc_col = alcohol_addiction_field
if (alc_col %in% names(char_dt)) {
  char_dt[, (alc_col) := {
    x = as.numeric(get(alc_col))
    x[x < 0] = NA
    as.integer(x == 1)
  }]
}

# --- 7h. Personality disorder (field 130932) ---
# Date-based: any non-NA date = 1 (has diagnosis), NA = 0
pd_col = personality_disorder_field
if (pd_col %in% names(char_dt)) {
  char_dt[, hx_personality_disorder := as.integer(!is.na(get(pd_col)) & get(pd_col) != "")]
}

# --- 7i. Bipolar disorder (field 20126) ---
# UKB field 20126 "Bipolar and major depression status" (derived):
#   1 = Bipolar I, 2 = Bipolar II, other values = non-bipolar categories
# Create binary: 1 if Bipolar I or II, 0 otherwise.
bp_col = bipolar_field
if (bp_col %in% names(char_dt)) {
  # Diagnostic: print raw value distribution to verify coding
  cat(sprintf("  Bipolar field (%s) raw distribution:\n", bp_col))
  print(table(char_dt[[bp_col]], useNA = "always"))
  
  char_dt[, bipolar_disorder_history := {
    x = as.numeric(get(bp_col))
    as.integer(x %in% c(1, 2))
  }]
  n_bp = sum(char_dt$bipolar_disorder_history == 1, na.rm = TRUE)
  cat(sprintf("  Bipolar positive (value 1 or 2): n = %d\n", n_bp))
  
  if (n_bp == 0) {
    cat("  Note: Bipolar indicator is 0 for all participants. Check:\n")
    cat("    (a) Does field 20126-0.0 exist with this exact name in csv_data?\n")
    cat("    (b) Are values 1/2 the correct codes for Bipolar I/II?\n")
    cat("    (c) If genuinely zero, this covariate adds nothing to sensitivity models.\n")
  }
}

# Merge medication data
char_dt = merge(char_dt, medication_dt, by = "subject_num", all.x = TRUE)

# 8. HELPER FUNCTIONS FOR SUMMARY STATISTICS

#' Compute mean (SD) with n available
#' @param x Numeric vector
#' @return Character string: "mean (SD) [n=X]"
mean_sd_n = function(x) {
  x = x[!is.na(x)]
  n = length(x)
  if (n == 0) return("NA [n=0]")
  sprintf("%.1f (%.1f) [n=%d]", mean(x), sd(x), n)
}

#' Compute n (%) with n available for binary variable
#' @param x Binary vector (0/1, may contain NA)
#' @return Character string: "n (pct%) [n_avail=X]"
n_pct = function(x) {
  x = x[!is.na(x)]
  n_avail = length(x)
  if (n_avail == 0) return("NA [n=0]")
  n_pos = sum(x == 1)
  sprintf("%d (%.1f%%) [n=%d]", n_pos, 100 * n_pos / n_avail, n_avail)
}

#' Compute median [IQR] with n available
#' @param x Numeric vector
#' @return Character string: "median [Q1-Q3] [n=X]"
median_iqr_n = function(x) {
  x = x[!is.na(x)]
  n = length(x)
  if (n == 0) return("NA [n=0]")
  q = quantile(x, probs = c(0.25, 0.50, 0.75))
  sprintf("%.1f [%.1f-%.1f] [n=%d]", q[2], q[1], q[3], n)
}

#' Flag variables with >50% missingness
#' @param x Vector
#' @param total_n Total sample size
#' @return Logical
high_missingness = function(x, total_n) {
  n_miss = sum(is.na(x))
  return(n_miss / total_n > 0.50)
}

# 9. COMPUTE CHARACTERIZATION TABLE

cat("\n=== Computing characterization table ===\n")

# Store results
results_list = list()

for (pheno_name in names(phenotype_registry)) {
  reg = phenotype_registry[[pheno_name]]
  dt  = as.data.table(reg$data)
  
  cat(sprintf("\nProcessing %s...\n", pheno_name))
  
  # Get subject IDs for cases and controls
  case_ids    = dt[dt[[reg$col]] == 1, subject_num]
  control_ids = dt[dt[[reg$col]] == 0, subject_num]
  
  # --- Demographics from csv_data (raw, not z-scored) ---
  # Age at imaging = field 21003-2.0, Sex = field 31-0.0 (0 = Female, 1 = Male)
  
  # Merge characterization variables (includes raw age and sex from csv_data)
  cases_char    = merge(data.table(subject_num = case_ids),    char_dt, by = "subject_num", all.x = TRUE)
  controls_char = merge(data.table(subject_num = control_ids), char_dt, by = "subject_num", all.x = TRUE)
  
  # --- MHQ completion and timing ---
  cases_mhq_completed   = sum(case_ids %in% mhq_completers)
  controls_mhq_completed = sum(control_ids %in% mhq_completers)
  
  cases_mhq_valid   = sum(case_ids %in% mhq_valid)
  controls_mhq_valid = sum(control_ids %in% mhq_valid)
  
  cases_mhq_excluded   = sum(case_ids %in% mri_before_mhq)
  controls_mhq_excluded = sum(control_ids %in% mri_before_mhq)
  
  n_cases    = length(case_ids)
  n_controls = length(control_ids)
  
  # --- Apply MHQ timing exclusion for MHQ-derived variables ---
  # For MHQ variables, restrict to participants with valid timing
  cases_char_mhq    = cases_char[subject_num %in% mhq_valid]
  controls_char_mhq = controls_char[subject_num %in% mhq_valid]
  
  # --- Build result row ---
  row = data.table(
    phenotype = pheno_name,
    category  = ifelse(pheno_name %in% c("ICD10Dep", "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur"),
                       "Deep", "Shallow"),
    
    # Sample sizes
    n_cases    = n_cases,
    n_controls = n_controls,
    
    # Demographics (raw values from csv_data, not z-scored)
    cases_age          = mean_sd_n(cases_char$age_at_imaging),
    controls_age       = mean_sd_n(controls_char$age_at_imaging),
    cases_pct_female   = sprintf("%.1f%% [n=%d]", 
                                 100 * sum(cases_char$sex_raw == 0, na.rm = TRUE) / 
                                   sum(!is.na(cases_char$sex_raw)), 
                                 sum(!is.na(cases_char$sex_raw))),
    controls_pct_female = sprintf("%.1f%% [n=%d]",
                                  100 * sum(controls_char$sex_raw == 0, na.rm = TRUE) / 
                                    sum(!is.na(controls_char$sex_raw)),
                                  sum(!is.na(controls_char$sex_raw))),
    
    # MHQ completion
    cases_mhq_completion   = sprintf("%d (%.1f%%) [n=%d]", cases_mhq_completed, 100 * cases_mhq_completed / n_cases, n_cases),
    controls_mhq_completion = sprintf("%d (%.1f%%) [n=%d]", controls_mhq_completed, 100 * controls_mhq_completed / n_controls, n_controls),
    cases_mhq_valid_timing   = sprintf("%d (%.1f%%) [n=%d]", cases_mhq_valid, 100 * cases_mhq_valid / n_cases, n_cases),
    controls_mhq_valid_timing = sprintf("%d (%.1f%%) [n=%d]", controls_mhq_valid, 100 * controls_mhq_valid / n_controls, n_controls),
    
    # Antidepressant use (from field 20003; NOT MHQ-dependent)
    cases_antidepressant   = n_pct(cases_char$any_antidepressant),
    controls_antidepressant = n_pct(controls_char$any_antidepressant),
    
    # --- MHQ-derived variables (apply timing exclusion) ---
    # PHQ-9 mean (0-3 scale)
    cases_phq9         = mean_sd_n(cases_char_mhq$phq9_mean),
    controls_phq9      = mean_sd_n(controls_char_mhq$phq9_mean),
    
    # Suicidal ideation
    cases_si           = n_pct(cases_char_mhq$suicidal_ideation_binary),
    controls_si        = n_pct(controls_char_mhq$suicidal_ideation_binary),
    
    # Duration of worst episode (cases only; meaningless for controls)
    cases_worst_episode_duration = mean_sd_n(cases_char_mhq[[worst_dur_col]]),
    
    # Anxiety proxy
    cases_anxiety      = n_pct(cases_char_mhq[[anx_col]]),
    controls_anxiety   = n_pct(controls_char_mhq[[anx_col]]),
    
    # Alcohol addiction
    cases_alcohol      = n_pct(cases_char_mhq[[alc_col]]),
    controls_alcohol   = n_pct(controls_char_mhq[[alc_col]]),
    
    # Personality disorder (HES-derived, not MHQ; no timing exclusion needed)
    cases_pd           = n_pct(cases_char$hx_personality_disorder),
    controls_pd        = n_pct(controls_char$hx_personality_disorder),
    
    # Bipolar disorder history
    cases_bipolar      = n_pct(cases_char$bipolar_disorder_history),
    controls_bipolar   = n_pct(controls_char$bipolar_disorder_history),
    
    # --- Chronicity (CASES ONLY, MHQ timing exclusion applied) ---
    cases_age_first_episode = mean_sd_n(cases_char_mhq[[age_first_col]]),
    cases_num_episodes      = median_iqr_n(cases_char_mhq[[num_ep_col]])
  )
  
  # --- Flag high-missingness variables ---
  mhq_vars_to_check = list(
    phq9_mean = cases_char_mhq$phq9_mean,
    anxiety_proxy = cases_char_mhq[[anx_col]],
    alcohol_addiction = cases_char_mhq[[alc_col]],
    age_first_episode = cases_char_mhq[[age_first_col]],
    num_episodes = cases_char_mhq[[num_ep_col]]
  )
  
  for (var_name in names(mhq_vars_to_check)) {
    if (high_missingness(mhq_vars_to_check[[var_name]], nrow(cases_char_mhq))) {
      cat(sprintf("  >50%% missingness for '%s' among %s cases (after MHQ timing exclusion).\n",
                  var_name, pheno_name))
    }
  }
  
  # Check if MHQ completion among cases is <50% (flag for sensitivity model eligibility)
  if (cases_mhq_valid / n_cases < 0.50) {
    cat(sprintf("  MHQ completion with valid timing is <50%% for %s cases (%d/%d = %.1f%%). MHQ-dependent sensitivity models should NOT be run for this phenotype.\n",
                pheno_name, cases_mhq_valid, n_cases, 100 * cases_mhq_valid / n_cases))
  }
  
  results_list[[pheno_name]] = row
}

# 10. COMBINE AND OUTPUT

results_table = rbindlist(results_list, fill = TRUE)

# Save as CSV
output_file = file.path(output_dir, "sample_characterization_table.csv")
fwrite(results_table, output_file)
cat(sprintf("\n=== Results saved to %s ===\n", output_file))

# Also save as RDS for downstream use
saveRDS(results_table, file = file.path(output_dir, "sample_characterization_table.RDS"))

# 11. PRINT FORMATTED TABLE

for (pheno_name in names(phenotype_registry)) {
  row = results_table[phenotype == pheno_name]
  cat(sprintf("--- %s (%s) ---\n", pheno_name, row$category))
  cat(sprintf("  N: Cases=%d, Controls=%d\n", row$n_cases, row$n_controls))
  cat(sprintf("  Age: Cases=%s, Controls=%s\n", row$cases_age, row$controls_age))
  cat(sprintf("  %% Female: Cases=%s, Controls=%s\n", row$cases_pct_female, row$controls_pct_female))
  cat(sprintf("  MHQ Completed: Cases=%s, Controls=%s\n", row$cases_mhq_completion, row$controls_mhq_completion))
  cat(sprintf("  MHQ Valid Timing: Cases=%s, Controls=%s\n", row$cases_mhq_valid_timing, row$controls_mhq_valid_timing))
  cat(sprintf("  Any Antidepressant: Cases=%s, Controls=%s\n", row$cases_antidepressant, row$controls_antidepressant))
  cat(sprintf("  PHQ-9 Sum (MHQ valid): Cases=%s, Controls=%s\n", row$cases_phq9, row$controls_phq9))
  cat(sprintf("  Suicidal Ideation: Cases=%s, Controls=%s\n", row$cases_si, row$controls_si))
  cat(sprintf("  Worst Episode Duration (cases): %s\n", row$cases_worst_episode_duration))
  cat(sprintf("  Anxiety Proxy: Cases=%s, Controls=%s\n", row$cases_anxiety, row$controls_anxiety))
  cat(sprintf("  Alcohol Addiction: Cases=%s, Controls=%s\n", row$cases_alcohol, row$controls_alcohol))
  cat(sprintf("  Personality Disorder: Cases=%s, Controls=%s\n", row$cases_pd, row$controls_pd))
  cat(sprintf("  Bipolar History: Cases=%s, Controls=%s\n", row$cases_bipolar, row$controls_bipolar))
  cat(sprintf("  Age at First Episode (cases, MHQ valid): %s\n", row$cases_age_first_episode))
  cat(sprintf("  Number of Episodes (cases, MHQ valid): %s\n", row$cases_num_episodes))
  cat("\n")
}

# 12. EXPORT MHQ COMPLETION RATES FOR SENSITIVITY MODEL ELIGIBILITY
# This table determines which phenotypes are eligible for MHQ-dependent models.

mhq_eligibility = results_table[, .(
  phenotype,
  n_cases,
  cases_mhq_valid = as.integer(sub("^(\\d+).*", "\\1", cases_mhq_valid_timing)),
  pct_mhq_valid   = as.numeric(sub(".*\\((.*)%\\)", "\\1", cases_mhq_valid_timing))
)]
mhq_eligibility[, eligible_mhq_models := pct_mhq_valid >= 50]

cat("\n=== MHQ Model Eligibility (>=50% completion with valid timing) ===\n")
print(mhq_eligibility)

saveRDS(mhq_eligibility, file = file.path(output_dir, "mhq_model_eligibility.RDS"))


# 13. EXPORT FORMATTED EXCEL WORKBOOK (Supplementary Table)
# One tab per phenotype, styled with headers, footnotes, and consistent formatting.

cat("\n=== Generating Excel workbook ===\n")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}
library(openxlsx)

wb = createWorkbook()

# Styles
header_style = createStyle(
  fontName = "Times New Roman", fontSize = 10, fontColour = "#FFFFFF",
  fgFill = "#4472C4", halign = "center", textDecoration = "bold",
  border = "Bottom", borderColour = "#2F5496"
)
header_left = createStyle(
  fontName = "Times New Roman", fontSize = 10, fontColour = "#FFFFFF",
  fgFill = "#4472C4", halign = "left", textDecoration = "bold",
  border = "Bottom", borderColour = "#2F5496"
)
title_style = createStyle(
  fontName = "Times New Roman", fontSize = 12, textDecoration = "bold", halign = "center"
)
subtitle_style = createStyle(
  fontName = "Times New Roman", fontSize = 10, textDecoration = "italic", halign = "center"
)
data_style = createStyle(fontName = "Times New Roman", fontSize = 10, halign = "center")
data_left  = createStyle(fontName = "Times New Roman", fontSize = 10, halign = "left")
note_style = createStyle(fontName = "Times New Roman", fontSize = 9, textDecoration = "italic")
row_border = createStyle(border = "Bottom", borderColour = "#D9D9D9", borderStyle = "thin")

# Variable labels and mapping to results_table columns
# Each entry: (label, cases_col, controls_col, cases_only)
var_map = list(
  list("Age at imaging, years",                   "cases_age",                    "controls_age",       FALSE),
  list("Female, %",                               "cases_pct_female",             "controls_pct_female", FALSE),
  list("Antidepressant use, n (%)",                "cases_antidepressant",         "controls_antidepressant", FALSE),
  list("MHQ completed, n (%)",                     "cases_mhq_completion",         "controls_mhq_completion", FALSE),
  list("MHQ with valid timing, n (%)",             "cases_mhq_valid_timing",       "controls_mhq_valid_timing", FALSE),
  list("PHQ-9 mean score (0-3)\u2020",             "cases_phq9",                   "controls_phq9",      FALSE),
  list("Suicidal ideation, n (%)\u2020",           "cases_si",                     "controls_si",        FALSE),
  list("Worst episode duration\u2020",             "cases_worst_episode_duration",  NA_character_,        TRUE),
  list("Anxiety (professional informed)\u2020",    "cases_anxiety",                "controls_anxiety",   FALSE),
  list("Alcohol addiction\u2020",                  "cases_alcohol",                "controls_alcohol",   FALSE),
  list("Personality disorder (HES)",               "cases_pd",                     "controls_pd",        FALSE),
  list("Bipolar disorder history",                 "cases_bipolar",                "controls_bipolar",   FALSE),
  list("Age at first episode, years\u2020\u2021",  "cases_age_first_episode",       NA_character_,        TRUE),
  list("Number of episodes, median [IQR]\u2020\u2021", "cases_num_episodes",        NA_character_,        TRUE)
)

# Helper: split "value [n=X]" or "value (X%) [n=X]" into value and n
split_val_n = function(x) {
  if (is.na(x) || x == "" || grepl("^NA", x)) return(list(val = "\u2014", n = "\u2014"))
  # Pattern: "value [n=X]" 
  m = regmatches(x, regexec("^(.+?)\\s*\\[n=(\\d+)\\]$", x))[[1]]
  if (length(m) == 3) return(list(val = trimws(m[2]), n = m[3]))
  # Pattern: "N (X%) [n=Y]" 
  m2 = regmatches(x, regexec("^(.+\\))\\s*\\[n=(\\d+)\\]$", x))[[1]]
  if (length(m2) == 3) return(list(val = trimws(m2[2]), n = m2[3]))
  # Pattern: "N (X%)" without [n=]
  m3 = regmatches(x, regexec("^(\\d+)\\s*\\((.+)\\)$", x))[[1]]
  if (length(m3) == 3) return(list(val = x, n = "\u2014"))
  return(list(val = x, n = "\u2014"))
}

for (pheno_name in names(phenotype_registry)) {
  reg = phenotype_registry[[pheno_name]]
  row = results_table[phenotype == pheno_name]
  
  addWorksheet(wb, sheetName = pheno_name)
  
  # Title
  mergeCells(wb, pheno_name, cols = 1:5, rows = 1)
  writeData(wb, pheno_name, x = pheno_name, startCol = 1, startRow = 1)
  addStyle(wb, pheno_name, title_style, rows = 1, cols = 1:5)
  
  # Subtitle
  mergeCells(wb, pheno_name, cols = 1:5, rows = 2)
  subtitle_text = sprintf("%s phenotype | Cases: n=%s | Controls: n=%s",
                          row$category,
                          format(row$n_cases, big.mark = ","),
                          format(row$n_controls, big.mark = ","))
  writeData(wb, pheno_name, x = subtitle_text, startCol = 1, startRow = 2)
  addStyle(wb, pheno_name, subtitle_style, rows = 2, cols = 1:5)
  
  # Headers
  headers = c("Variable", "Cases", "n", "Controls", "n")
  writeData(wb, pheno_name, x = t(headers), startCol = 1, startRow = 4, colNames = FALSE)
  addStyle(wb, pheno_name, header_left,  rows = 4, cols = 1)
  addStyle(wb, pheno_name, header_style, rows = 4, cols = 2:5)
  
  # Data rows
  for (i in seq_along(var_map)) {
    vm       = var_map[[i]]
    label    = vm[[1]]
    cases_col = vm[[2]]
    ctrl_col  = vm[[3]]
    cases_only = vm[[4]]
    
    r = 4 + i  # Excel row
    
    # Variable label
    writeData(wb, pheno_name, x = label, startCol = 1, startRow = r)
    addStyle(wb, pheno_name, data_left, rows = r, cols = 1)
    
    # Cases value + n
    cases_raw = as.character(row[[cases_col]])
    cs = split_val_n(cases_raw)
    writeData(wb, pheno_name, x = cs$val, startCol = 2, startRow = r)
    writeData(wb, pheno_name, x = cs$n,   startCol = 3, startRow = r)
    addStyle(wb, pheno_name, data_style, rows = r, cols = 2:3)
    
    # Controls value + n (or em dash if cases-only)
    if (cases_only || is.na(ctrl_col)) {
      writeData(wb, pheno_name, x = "\u2014", startCol = 4, startRow = r)
      writeData(wb, pheno_name, x = "\u2014", startCol = 5, startRow = r)
    } else {
      ctrl_raw = as.character(row[[ctrl_col]])
      ct = split_val_n(ctrl_raw)
      writeData(wb, pheno_name, x = ct$val, startCol = 4, startRow = r)
      writeData(wb, pheno_name, x = ct$n,   startCol = 5, startRow = r)
    }
    addStyle(wb, pheno_name, data_style, rows = r, cols = 4:5)
    
    # Row border
    addStyle(wb, pheno_name, row_border, rows = r, cols = 1:5, stack = TRUE)
  }
  
  # Footnotes
  fn_start = 4 + length(var_map) + 2
  footnotes = c(
    "Values are mean (SD) unless otherwise specified.",
    "\u2020 Restricted to participants with valid MHQ timing (MHQ completed before imaging visit).",
    "\u2021 Cases only. Chronicity variables are not applicable to controls.",
    "PHQ-9 mean score: mean of available items (minimum 7/9 required), on the 0-3 per-item scale.",
    "Worst episode duration: UKB field 20438, categorical 1-6 scale (1=<2 weeks, 6=>2 years).",
  )
  for (j in seq_along(footnotes)) {
    writeData(wb, pheno_name, x = footnotes[j], startCol = 1, startRow = fn_start + j - 1)
    addStyle(wb, pheno_name, note_style, rows = fn_start + j - 1, cols = 1)
  }
  
  # Column widths
  setColWidths(wb, pheno_name, cols = 1, widths = 42)
  setColWidths(wb, pheno_name, cols = 2, widths = 22)
  setColWidths(wb, pheno_name, cols = 3, widths = 10)
  setColWidths(wb, pheno_name, cols = 4, widths = 22)
  setColWidths(wb, pheno_name, cols = 5, widths = 10)
}

xlsx_path = file.path(output_dir, "Supplementary_Table_Clinical_Characterization.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(sprintf("Saved Excel workbook: %s\n", xlsx_path))
cat("Done. Proceed to S4\n")