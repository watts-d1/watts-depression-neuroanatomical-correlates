#!/usr/bin/env Rscript
# S4_sensitivity_analyses.R
# Sensitivity Analyses
#
# Case-control sensitivity models (each adds one covariate to the primary):
#   Primary (already run): IDP ~ depression + age + sex + center + eTIV + head_motion + PCs
#   Model 1 (+ antidepressant):  Primary + any_antidepressant (binary, field 20003)
#   Model 2 (+ anxiety):         Primary + anxiety_proxy (binary, field 20428, MHQ-derived)
#   Model 3 (+ PHQ-9 severity):  Primary + phq9_mean (continuous 0-3, MHQ-derived)
#   Model 4 (+ psychotherapy):   Primary + any_psychotherapy (binary, field 20547, MHQ-derived)
#
# Phenotype eligibility:
#   Model 1: All 9 phenotypes (antidepressant from assessment visit, no MHQ dependency)
#   Model 2: All 9 phenotypes (MHQ timing exclusion applied)
#   Model 3: LifetimeMDD and MDDRecur only (other phenotypes have <50% MHQ completion)
#   Model 4: All 9 phenotypes (MHQ timing exclusion applied)
#
# Within-cases models (episode count, onset age) are in S4b_sensitivity_within_cases.R.
#
# Design notes:
# 1. Antidepressant use and psychotherapy are plausible MEDIATORS (depression ->
#    treatment -> brain), not confounders. Adjustment attenuates the total effect.
#    Both are presented as sensitivity analyses, NOT "correct" models.
# 2. MHQ timing exclusion (n=8,245) applied before Models 2, 3, and 4.
# 3. Complete-case restriction: each model drops participants missing the new
#    covariate(s). Sample sizes will be smaller than primary analyses.
# 4. FDR correction applied within each phenotype, matching primary analysis.
# Author: Devon Watts

# --- User-defined paths (modify for your environment) ---
UKBB_DATA_PATH = "/path/to/ukbiobank/ukb674571.csv"  # UK Biobank phenotype CSV
S1_FUNCTIONS_PATH = "/path/to/S1_functions.R"  # Utility functions from S1 preprocessing
S1_Z_SCORE_PATH = "/path/to/z_score_neuroimaging_objects.RData"  # S1 output: z-scored neuroimaging data
S1_PRIMARY_RESULTS_DIR = "/path/to/sup_table2_results"  # S1 output: primary results directory

library(data.table)
library(future)
library(future.apply)
library(MASS)
library(openxlsx)
library(here)

# 0. Source the primary analysis functions

source(S1_FUNCTIONS_PATH)

# 1. Configuration

csv_file_path  = UKBB_DATA_PATH
z_score_path   = S1_Z_SCORE_PATH
output_dir     = here("output", "sensitivity_results")

# Primary ORs and FDR p-values needed for sensitivity vs primary comparison
primary_results_dir = S1_PRIMARY_RESULTS_DIR

primary_files = list(
  FA                 = file.path(primary_results_dir, "Supplementary Table S2a - Odds ratios across imaging-derived-phenotypes_FA.xlsx"),
  MD                 = file.path(primary_results_dir, "Supplementary Table S2b - Odds ratios across imaging-derived-phenotypes_MD.xlsx"),
  SA                 = file.path(primary_results_dir, "Supplementary Table S2c - Odds ratios across imaging-derived-phenotypes_SA.xlsx"),
  CT                 = file.path(primary_results_dir, "Supplementary Table S2d - Odds ratios across imaging-derived-phenotypes_CT.xlsx"),
  global_volume      = file.path(primary_results_dir, "Supplementary Table S2e - Odds ratios across imaging-derived-phenotypes_global_volume.xlsx"),
  subcortical_volume = file.path(primary_results_dir, "Supplementary Table S2f - Odds ratios across imaging-derived-phenotypes_subcortical_volume.xlsx")
)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Antidepressant codes (from coding4.tsv) ---
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

# Measurement ID vectors (needed for subsetting)
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
                        "25151", "25150", "25149", "25148")

# 2. Load data

cat("Loading z_score objects...\n")
load(z_score_path)

cat("Loading UKB CSV...\n")
csv_data = fread(csv_file_path, showProgress = TRUE, nThread = 20)

# 3. PHENOTYPE REGISTRY

phenotype_registry = list(
  GPNoDep            = list(data = z_score_completed_GPNoDep,            col = "GPNoDep"),
  GPPsy              = list(data = z_score_completed_GPpsy,              col = "GPpsy"),
  PsyPsy             = list(data = z_score_completed_Psypsy,             col = "Psypsy"),
  SelfRepDep         = list(data = z_score_completed_SelfRepDep,         col = "SelfRepDep"),
  DepAll             = list(data = z_score_completed_DepAll,             col = "DepAll"),
  ICD10Dep           = list(data = z_score_completed_ICD10Dep,           col = "ICD10Dep"),
  ICD10Dep_exclpsych = list(data = z_score_completed_ICD10Dep_exclpsych, col = "ICD10Dep.exclpsych"),
  LifetimeMDD        = list(data = z_score_completed_LifetimeMDD,        col = "LifetimeMDD"),
  MDDRecur           = list(data = z_score_completed_MDDRecurr,          col = "MDDRecur")
)

# 4. EXTRACT AND RECODE NEW COVARIATES FROM csv_data

cat("\n=== Extracting new covariates ===\n")

id_col = ifelse("subject_num" %in% names(csv_data), "subject_num", "eid")

# --- 4a. MHQ timing exclusion ---
mhq_date_col    = "20400-0.0"
imaging_date_col = "53-2.0"

mhq_timing = csv_data[, c(id_col, mhq_date_col, imaging_date_col), with = FALSE]
setnames(mhq_timing, c("subject_num", "mhq_date", "imaging_date"))
mhq_timing[, subject_num := as.integer(subject_num)]

if (is.character(mhq_timing$mhq_date))    mhq_timing[, mhq_date := as.Date(mhq_date)]
if (is.character(mhq_timing$imaging_date)) mhq_timing[, imaging_date := as.Date(imaging_date)]

mri_before_mhq_ids = mhq_timing[!is.na(mhq_date) & !is.na(imaging_date) & imaging_date < mhq_date, subject_num]
mhq_valid_ids      = mhq_timing[!is.na(mhq_date) & !is.na(imaging_date) & mhq_date <= imaging_date, subject_num]

cat(sprintf("  MHQ valid timing: n = %d\n", length(mhq_valid_ids)))
cat(sprintf("  MRI before MHQ (excluded): n = %d\n", length(mri_before_mhq_ids)))

# --- MHQ-imaging interval summary by phenotype ---
# Computed on the analysis sample for each phenotype, not the full csv_data.

mhq_timing[, mhq_imaging_interval_days := as.numeric(imaging_date - mhq_date)]

interval_summary = rbindlist(lapply(names(phenotype_registry), function(pheno_name) {
  
  reg = phenotype_registry[[pheno_name]]
  dt  = as.data.table(reg$data)
  
  # Merge interval column into phenotype sample
  dt = merge(dt, mhq_timing[, .(subject_num, mhq_imaging_interval_days)],
             by = "subject_num", all.x = TRUE)
  
  x = dt$mhq_imaging_interval_days
  x = x[!is.na(x)]
  
  if (length(x) == 0) {
    return(data.table(
      phenotype = pheno_name,
      n         = 0L,
      mean_days = NA_real_, sd_days   = NA_real_,
      median_days = NA_real_, q1_days = NA_real_, q3_days = NA_real_,
      mean_years = NA_real_, median_years = NA_real_
    ))
  }
  
  data.table(
    phenotype    = pheno_name,
    n            = length(x),
    mean_days    = mean(x),
    sd_days      = sd(x),
    median_days  = median(x),
    q1_days      = quantile(x, 0.25),
    q3_days      = quantile(x, 0.75),
    mean_years   = mean(x)   / 365.25,
    median_years = median(x) / 365.25
  )
}))

cat("\nMHQ-to-imaging interval (days) by phenotype:\n")
print(interval_summary[, .(
  phenotype, n,
  mean   = round(mean_days,   1),
  sd     = round(sd_days,     1),
  median = round(median_days, 1),
  Q1     = round(q1_days,     1),
  Q3     = round(q3_days,     1)
)])

cat("\nMHQ-to-imaging interval (years) by phenotype:\n")
print(interval_summary[, .(
  phenotype, n,
  mean   = round(mean_days   / 365.25, 2),
  sd     = round(sd_days     / 365.25, 2),
  median = round(median_days / 365.25, 2),
  Q1     = round(q1_days     / 365.25, 2),
  Q3     = round(q3_days     / 365.25, 2)
)])

saveRDS(interval_summary, file.path(output_dir, "mhq_imaging_interval_by_phenotype.RDS"))

# --- Initial assessment-to-imaging interval summary by phenotype ---
# Computed on the analysis sample for each phenotype, not the full csv_data.
initial_assessment_date_col = "53-0.0"
ia_timing = csv_data[, c(id_col, initial_assessment_date_col, imaging_date_col), with = FALSE]
setnames(ia_timing, c("subject_num", "initial_assessment_date", "imaging_date"))
ia_timing[, subject_num := as.integer(subject_num)]
if (is.character(ia_timing$initial_assessment_date)) ia_timing[, initial_assessment_date := as.Date(initial_assessment_date)]
if (is.character(ia_timing$imaging_date))            ia_timing[, imaging_date            := as.Date(imaging_date)]
ia_timing[, ia_imaging_interval_days := as.numeric(imaging_date - initial_assessment_date)]
assessment_imaging_interval_summary = rbindlist(lapply(names(phenotype_registry), function(pheno_name) {
  
  reg = phenotype_registry[[pheno_name]]
  dt  = as.data.table(reg$data)
  
  # Merge interval column into phenotype sample
  dt = merge(dt, ia_timing[, .(subject_num, ia_imaging_interval_days)],
             by = "subject_num", all.x = TRUE)
  
  x = dt$ia_imaging_interval_days
  x = x[!is.na(x)]
  
  if (length(x) == 0) {
    return(data.table(
      phenotype   = pheno_name,
      n           = 0L,
      mean_days   = NA_real_, sd_days      = NA_real_,
      median_days = NA_real_, q1_days      = NA_real_, q3_days = NA_real_,
      mean_years  = NA_real_, median_years = NA_real_
    ))
  }
  
  data.table(
    phenotype    = pheno_name,
    n            = length(x),
    mean_days    = mean(x),
    sd_days      = sd(x),
    median_days  = median(x),
    q1_days      = quantile(x, 0.25),
    q3_days      = quantile(x, 0.75),
    mean_years   = mean(x)   / 365.25,
    median_years = median(x) / 365.25
  )
}))
cat("\nInitial assessment-to-imaging interval (days) by phenotype:\n")
print(assessment_imaging_interval_summary[, .(
  phenotype, n,
  mean   = round(mean_days,   1),
  sd     = round(sd_days,     1),
  median = round(median_days, 1),
  Q1     = round(q1_days,     1),
  Q3     = round(q3_days,     1)
)])
cat("\nInitial assessment-to-imaging interval (years) by phenotype:\n")
print(assessment_imaging_interval_summary[, .(
  phenotype, n,
  mean   = round(mean_days   / 365.25, 2),
  sd     = round(sd_days     / 365.25, 2),
  median = round(median_days / 365.25, 2),
  Q1     = round(q1_days     / 365.25, 2),
  Q3     = round(q3_days     / 365.25, 2)
)])
saveRDS(assessment_imaging_interval_summary, file.path(output_dir, "assessment_imaging_interval_by_phenotype.RDS"))

# --- MHQ completion rate by phenotype ---
# Denominator = total n in the phenotype analysis sample (before any MHQ exclusion).
# Numerator   = participants with a non-missing mhq_date in that sample.
mhq_completion = rbindlist(lapply(names(phenotype_registry), function(pheno_name) {
  
  reg = phenotype_registry[[pheno_name]]
  dt  = as.data.table(reg$data)
  dt  = merge(dt, mhq_timing[, .(subject_num, mhq_date)],
              by = "subject_num", all.x = TRUE)
  
  n_total     = nrow(dt)
  n_mhq       = sum(!is.na(dt$mhq_date))
  n_valid_timing = sum(dt$subject_num %in% mhq_valid_ids)
  
  data.table(
    phenotype          = pheno_name,
    n_total            = n_total,
    n_mhq_completed    = n_mhq,
    pct_mhq_completed  = 100 * n_mhq / n_total,
    n_valid_timing     = n_valid_timing,
    pct_valid_timing   = 100 * n_valid_timing / n_total
  )
}))

cat("\nMHQ completion rate by phenotype:\n")
print(mhq_completion[, .(
  phenotype,
  n_total,
  n_mhq_completed,
  pct_completed  = round(pct_mhq_completed, 1),
  n_valid_timing,
  pct_valid      = round(pct_valid_timing, 1)
)])
saveRDS(mhq_completion, file.path(output_dir, "mhq_completion_by_phenotype.RDS"))

# --- 4b. Antidepressant use (field 20003, all array columns at instance 0) ---
# Searches 20003-0.0 through 20003-0.N for any antidepressant code.
# Each array index = a different medication reported, not a different visit.
med_cols = grep("^20003-0\\.", names(csv_data), value = TRUE)
cat(sprintf("  Found %d medication columns for field 20003\n", length(med_cols)))

med_dt = csv_data[, c(id_col, med_cols), with = FALSE]
setnames(med_dt, old = id_col, new = "subject_num")
med_dt[, subject_num := as.integer(subject_num)]
med_dt[, any_antidepressant := 0L]
for (col in med_cols) {
  med_dt[get(col) %in% antidepressant_codes, any_antidepressant := 1L]
}
med_dt = med_dt[, .(subject_num, any_antidepressant)]

cat(sprintf("  Antidepressant prevalence: %d (%.1f%%)\n",
            sum(med_dt$any_antidepressant == 1),
            100 * mean(med_dt$any_antidepressant == 1)))

# --- 4c. Comorbidity proxies ---
# Anxiety proxy (field 20428): 1 = Yes, 0 = No, negative = sentinel
#
# EXCLUDED COVARIATES:
# Alcohol addiction (field 20406): Conditional MHQ item asked only to participants
#   who screened positive on a prior alcohol question. Item-level response rates
#   were 1-5% of the imaging subsample across phenotypes, precluding its use as
#   a covariate. Prevalence is reported descriptively in Supplementary Table Clinical Characterization.
# Personality disorder (field 130932): 0% prevalence across all phenotypes.
# Bipolar disorder (field 20126): 0% prevalence across all phenotypes.

comorbidity_fields = c("20428-0.0")
available_comorbidity = intersect(comorbidity_fields, names(csv_data))

if (length(available_comorbidity) == 0) {
  cat("  Warning: Anxiety proxy field (20428-0.0) not found in csv_data.\n")
}
comor_dt = csv_data[, c(id_col, available_comorbidity), with = FALSE]
setnames(comor_dt, old = id_col, new = "subject_num")
comor_dt[, subject_num := as.integer(subject_num)]

# Anxiety proxy
if ("20428-0.0" %in% names(comor_dt)) {
  comor_dt[, anxiety_proxy := {
    x = as.numeric(`20428-0.0`)
    x[x < 0] = NA
    as.integer(x == 1)
  }]
} else {
  comor_dt[, anxiety_proxy := NA_integer_]
}

comor_dt = comor_dt[, .(subject_num, anxiety_proxy)]

# --- 4d. PHQ-9 sum score ---
phq_fields = c("20510-0.0", "20514-0.0", "20507-0.0", "20508-0.0",
               "20519-0.0", "20513-0.0", "20517-0.0", "20511-0.0", "20518-0.0")

available_phq = intersect(phq_fields, names(csv_data))
if (length(available_phq) == 9) {
  phq_dt = csv_data[, c(id_col, phq_fields), with = FALSE]
  setnames(phq_dt, old = id_col, new = "subject_num")
  phq_dt[, subject_num := as.integer(subject_num)]
  
  for (field in phq_fields) {
    phq_dt[, (field) := {
      x = as.numeric(get(field))
      x[x < 0] = NA
      x = x - 1  # Recode UKB 1-4 to PHQ 0-3
      x[!x %in% c(0, 1, 2, 3)] = NA
      x
    }]
  }
  
  phq_dt[, phq9_n_valid := rowSums(!is.na(.SD)), .SDcols = phq_fields]
  phq_dt[, phq9_mean := ifelse(phq9_n_valid >= 7,
                               rowMeans(.SD, na.rm = TRUE),
                               NA_real_), .SDcols = phq_fields]
  phq_dt = phq_dt[, .(subject_num, phq9_mean)]
  
  cat(sprintf("  PHQ-9 mean (0-3 scale): n valid = %d\n", sum(!is.na(phq_dt$phq9_mean))))
} else {
  cat("  Note: Not all PHQ fields found. Model 3 cannot be run.\n")
  phq_dt = data.table(subject_num = integer(0), phq9_mean = numeric(0))
}

# PHQ-9 per-item and total score descriptives by phenotype (imaging subsample)
#
# Extends the phq_dt construction to retain all 9 recoded items (0-3 scale)
# in phq_items_dt before collapsing to phq9_mean for Model 3.
#
# phq9_mean: row mean of available items (>=7/9 required), 0-3 scale.
# Mean rather than sum accommodates partial completion without listwise deletion.
#
# MRI-before-MHQ exclusion applied (mri_before_mhq_ids) for temporal validity.
# Control means for non-MHQ phenotypes (GP/ICD10) reflect a selected MHQ-
# completing minority and should be interpreted descriptively only.
#
# Outputs:
#   phq_items_dt       -- item-level data for imaging subsample descriptives
#   phq_dt             -- slim (subject_num + phq9_mean) for Model 3 downstream
#   phq9_summary       -- mean/SD/median PHQ-9 by phenotype x case/control
#   phq_wide           -- per-item means wide format for supplementary table
# Add row-level PHQ-9 mean (>=7 items required) to phq_items_dt
phq_items_dt[, phq9_n_valid := rowSums(!is.na(.SD)), .SDcols = phq_items]
phq_items_dt[, phq9_mean := fifelse(
  phq9_n_valid >= 7,
  rowMeans(.SD, na.rm = TRUE),
  NA_real_
), .SDcols = phq_items]

# Summarize by phenotype x case/control in imaging subsample
phq9_summary = rbindlist(lapply(names(phenotype_registry), function(pheno_name) {
  
  reg = phenotype_registry[[pheno_name]]
  dt  = as.data.table(reg$data)
  outcome_col = reg$col
  
  dt = merge(dt[, .(subject_num, outcome = get(outcome_col))],
             phq_items_dt[, .(subject_num, phq9_mean, phq9_n_valid)],
             by = "subject_num", all.x = TRUE)
  dt = dt[!subject_num %in% mri_before_mhq_ids]
  
  rbindlist(lapply(c(1L, 0L), function(grp) {
    sub = dt[outcome == grp & !is.na(phq9_mean)]
    data.table(
      phenotype   = pheno_name,
      group       = fifelse(grp == 1L, "case", "control"),
      n_total     = nrow(dt[outcome == grp]),
      n_phq_valid = nrow(sub),
      pct_valid   = round(100 * nrow(sub) / nrow(dt[outcome == grp]), 1),
      mean_phq9   = round(mean(sub$phq9_mean, na.rm = TRUE), 3),
      sd_phq9     = round(sd(sub$phq9_mean, na.rm = TRUE), 3),
      median_phq9 = round(median(sub$phq9_mean, na.rm = TRUE), 3),
      q25_phq9    = round(quantile(sub$phq9_mean, 0.25, na.rm = TRUE), 3),
      q75_phq9    = round(quantile(sub$phq9_mean, 0.75, na.rm = TRUE), 3)
    )
  }))
}))

# Wide format for clean reporting
phq9_wide = dcast(
  phq9_summary,
  phenotype ~ group,
  value.var = c("n_total", "n_phq_valid", "pct_valid", "mean_phq9", "sd_phq9")
)

print(phq9_summary[order(phenotype, group)])
# --- 4e. Psychotherapy (field 20547: Activities undertaken to treat depression) ---
# MHQ-derived field. Categorical (multiple), array indices 0-1.
# Coding 1406: 1 = Talking therapies (psychotherapy, counselling, group therapy, CBT)
#              3 = Other therapeutic activities (mindfulness, yoga, art classes)
#             -818 = Prefer not to answer
# Binary: 1 if value == 1 in any array position, 0 if valid response but != 1, NA otherwise.

psychotherapy_cols = grep("^20547-0\\.", names(csv_data), value = TRUE)
cat(sprintf("  Found %d columns for field 20547 (psychotherapy)\n", length(psychotherapy_cols)))

if (length(psychotherapy_cols) > 0) {
  psych_dt = csv_data[, c(id_col, psychotherapy_cols), with = FALSE]
  setnames(psych_dt, old = id_col, new = "subject_num")
  psych_dt[, subject_num := as.integer(subject_num)]
  
  # Check for talking therapies (value == 1) across all array positions
  psych_dt[, any_psychotherapy := NA_integer_]
  for (col in psychotherapy_cols) {
    psych_dt[, (col) := {
      x = as.numeric(get(col))
      x[x < 0] = NA  # Set sentinels (-818) to NA
      x
    }]
  }
  
  # Has valid response in any column (not all NA)
  psych_dt[, has_valid := rowSums(!is.na(.SD)) > 0, .SDcols = psychotherapy_cols]
  
  # Check if value 1 (talking therapies) appears in any column
  psych_dt[, endorsed_therapy := rowSums(.SD == 1, na.rm = TRUE) > 0, .SDcols = psychotherapy_cols]
  
  # Code: 1 if endorsed talking therapies, 0 if valid response but did not, NA otherwise
  psych_dt[has_valid == TRUE, any_psychotherapy := as.integer(endorsed_therapy)]
  
  psych_dt = psych_dt[, .(subject_num, any_psychotherapy)]
  
  cat(sprintf("  Psychotherapy (talking therapies): %d endorsed (%.1f%% of valid responses)\n",
              sum(psych_dt$any_psychotherapy == 1, na.rm = TRUE),
              100 * mean(psych_dt$any_psychotherapy == 1, na.rm = TRUE)))
  cat(sprintf("  Missing/sentinel: %d\n", sum(is.na(psych_dt$any_psychotherapy))))
} else {
  cat("  Note: Field 20547 not found in csv_data. Psychotherapy model cannot be run.\n")
  psych_dt = data.table(subject_num = integer(0), any_psychotherapy = integer(0))
}

# 5. DEFINE SENSITIVITY MODEL SPECIFICATIONS

# Base covariates (from primary analysis)
base_covariates = c("bl_31", "bl_54", "bl_21022", "bl_26521", "bl_24419",
                    paste0("PC", 1:20))

# Model specifications
model_specs = list(
  
  Model_1 = list(
    label       = "Model 1: Primary + Antidepressant Use",
    description = "Adjusts for any antidepressant medication use (binary, from field 20003).
Note: Antidepressant use is likely a mediator on the causal path
depression -> medication -> brain changes, rather than a confounder. Adjusting
for a mediator attenuates the total effect. This model tests what the
depression-brain association looks like after removing treatment-related
variance, but the attenuated estimate is NOT more 'correct' than the
primary model. Presented as sensitivity analysis only.",
    new_covariates = "any_antidepressant",
    merge_dt       = med_dt,
    requires_mhq_exclusion = FALSE,
    phenotypes = c("GPNoDep", "GPPsy", "PsyPsy", "SelfRepDep", "DepAll",
                   "ICD10Dep", "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur")
  ),
  
  Model_2 = list(
    label       = "Model 2: Primary + Anxiety",
    description = "Adjusts for anxiety proxy (professional informed, field 20428).
MHQ-derived -> MHQ timing exclusion applied. This is a help-seeking proxy
(whether a professional was informed about anxiety or nerves), not a
structured anxiety diagnosis.
Alcohol addiction (field 20406) was assessed but is a conditional MHQ item
with item-level response rates of 1-5% in the imaging subsample, precluding
its use as a covariate. All three are reported descriptively
in Supplementary Table S3.",
    new_covariates = "anxiety_proxy",
    merge_dt       = comor_dt,
    requires_mhq_exclusion = TRUE,
    phenotypes = c("GPNoDep", "GPPsy", "PsyPsy", "SelfRepDep", "DepAll",
                   "ICD10Dep", "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur")
  ),
  
  Model_3 = list(
    label       = "Model 3: Primary + PHQ-9 Severity",
    description = "Adjusts for current depressive symptom severity (PHQ-9 mean, 0-3 scale).
Mean of available items (minimum 7/9 required), which accommodates partial
missingness without listwise deletion.
Restricted to LifetimeMDD and MDDRecur only because:
(a) These phenotypes are defined from the MHQ, so nearly all cases have PHQ data,
    minimizing selection bias from complete-case restriction.
(b) Shallow phenotypes have low MHQ completion, creating severe selection bias.
(c) ICD10Dep cases may have been diagnosed years before the MHQ; PHQ is temporally
    disconnected from the phenotype-defining event.
This model tests whether adjusting for current severity further attenuates
depression-brain associations, informing the UKB-ENIGMA comparison.",
    new_covariates = "phq9_mean",
    merge_dt       = phq_dt,
    requires_mhq_exclusion = TRUE,
    phenotypes = c("LifetimeMDD", "MDDRecur")
  ),
  
  Model_4 = list(
    label       = "Model 4: Primary + Psychotherapy",
    description = "Adjusts for lifetime psychotherapy use (binary, from field 20547).
Field 20547 asks about 'Activities undertaken to treat depression' and is
MHQ-derived, so MHQ timing exclusion is applied. Coded 1 if participant
endorsed 'Talking therapies, such as psychotherapy, counselling, group
therapy or CBT' (coding value 1), 0 if a valid response was given but
talking therapies were not endorsed. Sentinel values (-818 = 'Prefer not
to answer') set to NA.
Like antidepressant use, psychotherapy is a plausible mediator rather than
a confounder (depression precedes treatment-seeking). The same mediator
framing applies: adjustment attenuates the total effect.",
    new_covariates = "any_psychotherapy",
    merge_dt       = psych_dt,
    requires_mhq_exclusion = TRUE,
    phenotypes = c("GPNoDep", "GPPsy", "PsyPsy", "SelfRepDep", "DepAll",
                   "ICD10Dep", "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur")
  )
)

# 6. PARALLEL EXECUTION SETUP

# Use sequential execution for profile likelihood CI (consistent with symptom
# dimension script rationale: GLM + profiling is fast per-feature at these
# sample sizes, and parallel overhead exceeds savings).
plan(sequential)

# 7. CORE SENSITIVITY ANALYSIS FUNCTION

#' Run a single sensitivity model for one phenotype
#'
#' Merges new covariates into the z_score object, applies MHQ timing exclusion
#' if required, separates features/covariates/outcome, and runs logistic
#' regression for all features using profile likelihood CIs.
#'
#' @param pheno_name Character: phenotype name (key into phenotype_registry)
#' @param model_spec List: model specification from model_specs
#' @param base_covariates Character vector of baseline covariate column names
#' @param mhq_valid_ids Integer vector of subject_nums with valid MHQ timing
#' @param mri_before_mhq_ids Integer vector of subject_nums to exclude
#' @return List with results_df and metadata
run_sensitivity_model = function(pheno_name, model_spec, base_covariates,
                                 mhq_valid_ids, mri_before_mhq_ids) {
  
  reg = phenotype_registry[[pheno_name]]
  dt  = as.data.table(reg$data)
  outcome_col = reg$col
  
  n_original = nrow(dt)
  cat(sprintf("\n  %s: Starting with n = %d\n", pheno_name, n_original))
  
  # Initialize sample_info (will be returned even on skip/error)
  sample_info = list(
    n_original  = n_original,
    n_post_mhq  = n_original,  # updated below if MHQ exclusion applied
    n_post_cc   = NA_integer_,
    n_cases     = NA_integer_,
    n_controls  = NA_integer_,
    pct_dropped = NA_real_,
    status      = "completed"
  )
  
  # --- Step 1: Merge new covariates ---
  dt = merge(dt, model_spec$merge_dt, by = "subject_num", all.x = TRUE)
  
  # --- Step 2: Apply MHQ timing exclusion if required ---
  if (model_spec$requires_mhq_exclusion) {
    # Exclude participants whose MRI preceded MHQ
    n_before_exclusion = nrow(dt)
    dt = dt[!subject_num %in% mri_before_mhq_ids]
    # Also exclude those who never completed MHQ (they have no MHQ-derived data)
    # For MHQ-dependent covariates, these will be NA anyway and dropped at complete-case,
    # but excluding explicitly makes the sample size reporting clearer.
    n_after_exclusion = nrow(dt)
    n_excluded = n_before_exclusion - n_after_exclusion
    cat(sprintf("    MHQ timing exclusion: removed %d, remaining n = %d\n",
                n_excluded, n_after_exclusion))
    sample_info$n_post_mhq = n_after_exclusion
  }
  
  # --- Step 3: Define full covariate set ---
  all_covariates = c(base_covariates, model_spec$new_covariates)
  
  # --- Step 4: Complete-case restriction on new covariates ---
  n_before_cc = nrow(dt)
  for (cov in model_spec$new_covariates) {
    if (cov %in% names(dt)) {
      n_missing = sum(is.na(dt[[cov]]))
      if (n_missing > 0) {
        cat(sprintf("    Covariate '%s': %d missing (%.1f%%)\n",
                    cov, n_missing, 100 * n_missing / nrow(dt)))
      }
    }
  }
  
  # Drop rows missing any new covariate
  # NOTE: Must extract column names to a local variable first.
  # data.table's .. prefix does not parse compound expressions like
  # model_spec$new_covariates correctly when length > 1.
  new_cov_cols = model_spec$new_covariates
  complete_mask = complete.cases(dt[, ..new_cov_cols])
  dt = dt[complete_mask]
  n_after_cc = nrow(dt)
  n_dropped_cc = n_before_cc - n_after_cc
  pct_dropped = 100 * n_dropped_cc / n_before_cc
  
  cat(sprintf("    Complete-case restriction: dropped %d (%.1f%%), remaining n = %d\n",
              n_dropped_cc, pct_dropped, n_after_cc))
  
  if (pct_dropped > 90) {
    cat(sprintf("    [SKIPPING] >90%% of sample dropped for %s in %s. Results would be unreliable.\n",
                pheno_name, model_spec$label))
    sample_info$n_post_cc   = n_after_cc
    sample_info$pct_dropped = pct_dropped
    sample_info$status      = "skipped_missingness"
    return(list(result = NULL, sample_info = sample_info))
  } else if (pct_dropped > 20) {
    cat(sprintf("   Complete-case restriction dropped >20%% of sample for %s in %s.\n",
                pheno_name, model_spec$label))
    cat("    Consider whether this introduces selection bias.\n")
  }
  
  # Report case/control split after restrictions
  n_cases_final    = sum(dt[[outcome_col]] == 1)
  n_controls_final = sum(dt[[outcome_col]] == 0)
  cat(sprintf("    Final: Cases = %d, Controls = %d\n", n_cases_final, n_controls_final))
  
  # --- Step 4b: Check covariate validity by case/control status ---
  # Flags covariates that are conditional on case status (e.g., field 20547
  # is only asked of participants who endorsed depression in the MHQ).
  # If controls have near-zero valid data, the covariate cannot adjust a
  # case-control comparison and the model is uninterpretable.
  for (cov in model_spec$new_covariates) {
    if (cov %in% names(dt)) {
      n_valid_cases    = sum(!is.na(dt[dt[[outcome_col]] == 1, get(cov)]))
      n_valid_controls = sum(!is.na(dt[dt[[outcome_col]] == 0, get(cov)]))
      pct_controls_valid = 100 * n_valid_controls / n_controls_final
      cat(sprintf("    Covariate '%s' validity: Cases = %d/%d valid, Controls = %d/%d valid (%.1f%%)\n",
                  cov, n_valid_cases, n_cases_final,
                  n_valid_controls, n_controls_final, pct_controls_valid))
      if (pct_controls_valid < 5) {
        cat(sprintf("    [Warning] <5%% of controls have valid '%s' data.\n", cov))
        cat("    This covariate may be conditional on case status (asked only of\n")
        cat("    participants who endorsed depression). The model compares treated\n")
        cat("    vs. untreated cases, not adjusted case-control. Results should be\n")
        cat("    reported descriptively, not as a sensitivity analysis.\n")
        sample_info$status = "flagged_conditional_covariate"
      }
    }
  }
  
  sample_info$n_post_cc   = n_after_cc
  sample_info$n_cases     = n_cases_final
  sample_info$n_controls  = n_controls_final
  sample_info$pct_dropped = pct_dropped
  
  if (n_cases_final < 50) {
    cat(sprintf("   Fewer than 50 cases for %s. Results may be unstable.\n", pheno_name))
    sample_info$status = "skipped_low_n"
    return(list(result = NULL, sample_info = sample_info))
  }
  
  # --- Step 5: Separate features, covariates, outcome ---
  # Use the same logic as the primary pipeline
  bl_cols   = grep("^bl_", names(dt), value = TRUE)
  feat_cols = setdiff(bl_cols, all_covariates)
  
  features_df   = dt[, ..feat_cols]
  covariates_df = dt[, ..all_covariates]
  outcome_df    = dt[, ..outcome_col]
  
  # --- Step 6: Check for covariate-phenotype definition overlap ---
  # Note:: confirm this covariate is independent of the case definition]
  
  # Model 2 specific checks:
  if ("anxiety_proxy" %in% model_spec$new_covariates) {
    # Anxiety proxy (professional informed about anxiety) partly overlaps with
    # GP/psychiatrist consultation phenotypes. The question asks about help-seeking
    # for anxiety, while GPPsy/PsyPsy ask about help-seeking for nerves/anxiety/
    # tension/depression. This is not perfectly independent.
    if (pheno_name %in% c("GPPsy", "PsyPsy", "GPNoDep")) {
      cat(sprintf("    Note: anxiety_proxy (field 20428) asks about professional contact for anxiety. %s is defined by GP/psychiatrist contact for nerves/anxiety/tension/depression. These share help-seeking behavior as a common cause. The covariate is not fully independent of the phenotype definition.\n", pheno_name))
    }
  }
  
  # --- Step 7: Run logistic regression for all features ---
  cat(sprintf("    Running logistic regression for %d features...\n", ncol(features_df)))
  
  # Capture warnings per model to track separation and convergence issues
  model_warnings = character(0)
  result = tryCatch(
    withCallingHandlers(
      run_analysis_for_a_dataset(
        features_df          = features_df,
        covariates_df        = covariates_df,
        outcome_df           = outcome_df,
        ci_method            = "profile_likelihood",
        conf_level           = 0.95,
        bootstrap_iterations = NULL
      ),
      warning = function(w) {
        model_warnings <<- c(model_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      cat(sprintf("    ERROR: %s\n", e$message))
      return(NULL)
    }
  )
  
  # Report warnings summary for this model x phenotype
  if (length(model_warnings) > 0) {
    warn_table = table(model_warnings)
    cat(sprintf("    Warnings captured: %d total\n", length(model_warnings)))
    for (wname in names(warn_table)) {
      cat(sprintf("      '%s': %d occurrences\n",
                  substr(wname, 1, 80), warn_table[wname]))
    }
    # Flag separation specifically
    n_separation = sum(grepl("fitted probabilities numerically 0 or 1", model_warnings))
    if (n_separation > 0) {
      cat(sprintf("    Note: %d features had quasi-complete separation in %s / %s.\n",
                  n_separation, model_spec$label, pheno_name))
      cat("    This may indicate the covariate(s) nearly perfectly predict the outcome\n")
      cat("    for some feature values. Results for these features should be interpreted cautiously.\n")
    }
  } else {
    cat("    No warnings.\n")
  }
  
  if (is.null(result)) {
    sample_info$status = "skipped_error"
    return(list(result = NULL, sample_info = sample_info))
  }
  
  # --- Step 8: Annotate results with metadata ---
  for (mtype in names(result)) {
    if (!is.null(result[[mtype]]$results_df) && nrow(result[[mtype]]$results_df) > 0) {
      result[[mtype]]$results_df$phenotype         = pheno_name
      result[[mtype]]$results_df$model              = model_spec$label
      result[[mtype]]$results_df$measurement_type   = mtype
      result[[mtype]]$results_df$n_cases_final      = n_cases_final
      result[[mtype]]$results_df$n_controls_final   = n_controls_final
      result[[mtype]]$results_df$n_total_final      = n_after_cc
      result[[mtype]]$results_df$pct_dropped_cc     = pct_dropped
    }
  }
  
  return(list(result = result, sample_info = sample_info))
}

# 8. RUN ALL SENSITIVITY MODELS

all_sensitivity_results = list()

# Global sample size tracker: one row per model x phenotype
sample_size_tracker = data.table(
  model       = character(0),
  model_label = character(0),
  phenotype   = character(0),
  n_original  = integer(0),
  n_post_mhq  = integer(0),
  n_post_cc   = integer(0),
  n_cases     = integer(0),
  n_controls  = integer(0),
  pct_dropped = numeric(0),
  status      = character(0)  # "completed", "skipped_missingness", "skipped_error", "skipped_low_n"
)

for (model_name in names(model_specs)) {
  spec = model_specs[[model_name]]
  
  cat(sprintf("\n\n{'='*60}\n"))
  cat(sprintf("=== %s ===\n", spec$label))
  cat(sprintf("Description: %s\n", gsub("\n", "\n  ", spec$description)))
  cat(sprintf("Phenotypes: %s\n", paste(spec$phenotypes, collapse = ", ")))
  cat(sprintf("New covariates: %s\n", paste(spec$new_covariates, collapse = ", ")))
  cat(sprintf("MHQ timing exclusion: %s\n", spec$requires_mhq_exclusion))
  cat(sprintf("{'='*60}\n"))
  
  model_results = list()
  
  for (pheno_name in spec$phenotypes) {
    run_output = run_sensitivity_model(
      pheno_name       = pheno_name,
      model_spec       = spec,
      base_covariates  = base_covariates,
      mhq_valid_ids    = mhq_valid_ids,
      mri_before_mhq_ids = mri_before_mhq_ids
    )
    
    # Populate sample size tracker
    si = run_output$sample_info
    sample_size_tracker = rbindlist(list(
      sample_size_tracker,
      data.table(
        model       = model_name,
        model_label = spec$label,
        phenotype   = pheno_name,
        n_original  = si$n_original,
        n_post_mhq  = si$n_post_mhq,
        n_post_cc   = si$n_post_cc,
        n_cases     = si$n_cases,
        n_controls  = si$n_controls,
        pct_dropped = si$pct_dropped,
        status      = si$status
      )
    ), fill = TRUE)
    
    result = run_output$result
    if (!is.null(result)) {
      model_results[[pheno_name]] = result
      
      # Save individual result
      saveRDS(result, file = file.path(output_dir,
                                       sprintf("sensitivity_%s_%s.RDS", model_name, pheno_name)))
      cat(sprintf("    Saved: sensitivity_%s_%s.RDS\n", model_name, pheno_name))
    }
  }
  
  all_sensitivity_results[[model_name]] = model_results
}

saveRDS(all_sensitivity_results, file.path(output_dir, "all_sensitivity_results.RDS"))

# 9. COMPILE RESULTS ACROSS EACH IDP

cat("\n\n=== Compiling sensitivity results ===\n")

# --- 9a. Modality classifier: map bl_XXXXX features to FA/MD/SA/CT ---
# These field IDs match the primary analysis supplementary tables (S2a-S2d).
# Field IDs already specified above (i.e. line 82-112)
# Build the lookup: bl_XXXXX -> modality
classify_feature_modality = function(feature_name) {
  # Extract numeric ID from bl_XXXXX
  fid = sub("^bl_", "", feature_name)
  
  if (fid %in% fa_measures) return("FA")
  if (fid %in% md_measures) return("MD")
  if (fid %in% surface_area_ids) return("SA")
  if (fid %in% thickness_ids) return("CT")
  if (fid %in% t1_global_sMRI) return("global_volume")
  if (fid %in% t1_subcortical_sMRI) return("subcortical_volume")
  return("Unknown")
}

# Vectorized version
assign_modality = function(features) {
  sapply(features, classify_feature_modality, USE.NAMES = FALSE)
}

# --- 9b. Compile with modality assignment ---

compile_sensitivity_results = function(model_results) {
  rows = list()
  
  for (pheno_name in names(model_results)) {
    result = model_results[[pheno_name]]
    df = result[["All_Features"]]$results_df
    if (!is.null(df) && nrow(df) > 0) {
      rows[[length(rows) + 1]] = df
    }
  }
  
  if (length(rows) == 0) return(data.table())
  compiled = rbindlist(rows, fill = TRUE)
  
  # Assign modality to each feature
  compiled[, modality := assign_modality(feature)]
  
  # Report classification
  mod_counts = compiled[, .N, by = modality]
  cat("  Modality assignment:\n")
  for (i in seq_len(nrow(mod_counts))) {
    cat(sprintf("    %s: %d rows\n", mod_counts$modality[i], mod_counts$N[i]))
  }
  if (any(mod_counts$modality == "Unknown")) {
    unknown_feats = unique(compiled[modality == "Unknown", feature])
    cat(sprintf("  Warning: %d features could not be classified: %s\n",
                length(unknown_feats), paste(head(unknown_feats, 5), collapse = ", ")))
    cat("  Note:: Check whether t1_global or t1_subcortical should map to a different modality.]\n")
  }
  
  return(compiled)
}

# Compile each model's results
for (model_name in names(all_sensitivity_results)) {
  cat(sprintf("\n  Compiling %s...\n", model_name))
  compiled = compile_sensitivity_results(all_sensitivity_results[[model_name]])
  
  if (nrow(compiled) > 0) {
    output_file = file.path(output_dir, sprintf("%s_compiled.RDS", model_name))
    saveRDS(compiled, output_file)
    cat(sprintf("  %s: %d rows compiled and saved.\n", model_name, nrow(compiled)))
    
    # Also save per-modality CSVs for inspection
    for (mod in unique(compiled$modality)) {
      mod_subset = compiled[modality == mod]
      mod_file = file.path(output_dir, sprintf("%s_%s_compiled.csv", model_name, mod))
      fwrite(mod_subset, mod_file)
      cat(sprintf("    Saved %s subset: %d rows -> %s\n", mod, nrow(mod_subset), basename(mod_file)))
    }
  }
}

# 9c. LOAD AND PARSE PRIMARY RESULTS FROM EXCEL (all 6 modalities)
#
# The primary supplementary tables S2a-S2f use a stacked multi-phenotype
# format within a single sheet. Each phenotype block is preceded by a header
# row where column 1 = "feature". The phenotype name appears either as the
# Excel column header (first block only) or in the row immediately above
# subsequent header rows. Data rows have bl_XXXXX in column 1.
#
# Column layout (consistent across S2a-S2f):
#   1: feature (bl_XXXXX)
#   2: region_name
#   3: odds_ratio
#   4: lower_CI
#   5: upper_CI
#   6: p_value (uncorrected)
#   7: bonferroni p-value
#   8: fdr_p_value

cat("\n\n=== Loading primary results from supplementary tables S2a-S2f ===\n")

library(openxlsx)

# Canonical phenotype name mapping (applied to all parsed results)
pheno_canonical = c(
  "gpnodep"            = "GPNoDep",
  "gppsy"              = "GPPsy",
  "psypsy"             = "PsyPsy",
  "selfrepdep"         = "SelfRepDep",
  "depall"             = "DepAll",
  "icd10dep"           = "ICD10Dep",
  "icd10dep_exclpsych" = "ICD10Dep_exclpsych",
  "lifetimemdd"        = "LifetimeMDD",
  "mddrecur"           = "MDDRecur"
)

# Valid field IDs per modality (from Section 1 measurement ID vectors).
# Used to filter out misassigned feature IDs in primary Excel files
# (e.g., S2b contains FA field IDs in some phenotype blocks).
valid_features_by_modality = list(
  FA                 = paste0("bl_", fa_measures),
  MD                 = paste0("bl_", md_measures),
  SA                 = paste0("bl_", surface_area_ids),
  CT                 = paste0("bl_", thickness_ids),
  global_volume      = paste0("bl_", t1_global_sMRI),
  subcortical_volume = paste0("bl_", t1_subcortical_sMRI)
)

#' Parse a stacked multi-phenotype supplementary table Excel file
#'
#' Reads one of the S2a-S2f Excel files and returns a data.table with
#' columns: feature, region_name, phenotype, odds_ratio, lower_CI, upper_CI,
#' p_value, fdr_p_value, modality. All phenotype names are normalized to
#' canonical form. Feature IDs are filtered to the valid set for the
#' specified modality to handle misassigned field IDs in the source Excel.
#'
#' @param filepath Path to the .xlsx file
#' @param modality_label Character string: FA, MD, SA, CT, global_volume,
#'   or subcortical_volume
#' @return data.table with one row per feature x phenotype
parse_primary_xlsx = function(filepath, modality_label) {
  
  if (!file.exists(filepath)) {
    cat(sprintf("  Warning: File not found: %s\n", filepath))
    return(data.table())
  }
  
  # --- Read twice: once with colNames to recover first phenotype name,
  #     once without for uniform row-level parsing ---
  raw_named = as.data.table(read.xlsx(filepath, sheet = 1, colNames = TRUE))
  first_pheno_raw = names(raw_named)[1]
  
  raw = as.data.table(read.xlsx(filepath, sheet = 1, colNames = FALSE))
  
  # --- Identify structural rows ---
  is_header = which(raw$X1 == "feature")
  is_data   = which(grepl("^bl_", raw$X1))
  
  if (length(is_header) == 0) {
    cat(sprintf("  Warning: No header rows found in %s\n", basename(filepath)))
    return(data.table())
  }
  
  # --- Assign each header row to a phenotype name ---
  pheno_names = character(length(is_header))
  pheno_names[1] = first_pheno_raw
  
  if (length(is_header) > 1) {
    for (i in 2:length(is_header)) {
      h = is_header[i]
      candidate = NA_character_
      for (offset in 1:3) {
        row_idx = h - offset
        if (row_idx < 1) break
        val = raw$X1[row_idx]
        if (!is.na(val) && val != "" && val != "feature" && !grepl("^bl_", val)) {
          candidate = val
          break
        }
      }
      if (is.na(candidate)) {
        cat(sprintf("  Warning: Could not identify phenotype for header at row %d in %s\n",
                    h, basename(filepath)))
        candidate = sprintf("Unknown_%d", i)
      }
      pheno_names[i] = candidate
    }
  }
  
  # --- Extract data rows for each phenotype block ---
  results = vector("list", length(is_header))
  
  for (i in seq_along(is_header)) {
    h_start = is_header[i] + 1
    if (i < length(is_header)) {
      h_end = is_header[i + 1] - 1
    } else {
      h_end = nrow(raw)
    }
    
    block_idx = is_data[is_data >= h_start & is_data <= h_end]
    if (length(block_idx) == 0) next
    
    block = raw[block_idx]
    
    # Extract bl_XXXXX portion only (handles any appended text)
    feature_ids = sub("^(bl_\\d+).*$", "\\1", block$X1)
    
    results[[i]] = data.table(
      feature     = feature_ids,
      region_name = block$X2,
      phenotype   = pheno_names[i],
      odds_ratio  = as.numeric(block$X3),
      lower_CI    = as.numeric(block$X4),
      upper_CI    = as.numeric(block$X5),
      p_value     = as.numeric(block$X6),
      fdr_p_value = as.numeric(block$X8),
      modality    = modality_label
    )
  }
  
  out = rbindlist(results[lengths(results) > 0])
  
  # --- Filter to valid field IDs for this modality ---
  # S2b (MD) is known to contain FA field IDs in some phenotype blocks.
  # This drops any row whose feature ID does not belong to the modality.
  valid_ids = valid_features_by_modality[[modality_label]]
  if (!is.null(valid_ids)) {
    n_before = nrow(out)
    invalid = out[!feature %in% valid_ids]
    out = out[feature %in% valid_ids]
    n_dropped = n_before - nrow(out)
    if (n_dropped > 0) {
      invalid_ids = unique(invalid$feature)
      cat(sprintf("  %s: dropped %d rows with %d invalid feature IDs (not in expected set)\n",
                  modality_label, n_dropped, length(invalid_ids)))
      cat(sprintf("    Examples: %s\n", paste(head(sort(invalid_ids), 5), collapse = ", ")))
    }
  }
  
  # Normalize phenotype names
  out[, phenotype := pheno_canonical[tolower(phenotype)]]
  
  n_na = sum(is.na(out$phenotype))
  if (n_na > 0) {
    cat(sprintf("  Warning: %d rows with unresolved phenotype names in %s\n",
                n_na, basename(filepath)))
  }
  
  cat(sprintf("  %s: %d rows, %d unique features, %d phenotypes, %d FDR-sig (p < 0.05)\n",
              modality_label, nrow(out), length(unique(out$feature)),
              length(unique(out$phenotype)),
              sum(out$fdr_p_value < 0.05, na.rm = TRUE)))
  
  return(out)
}

# --- Parse all 6 supplementary tables ---
# primary_files is defined in Section 1 and contains paths for:
#   FA, MD, SA, CT, global_volume, subcortical_volume

primary_compiled = rbindlist(
  lapply(names(primary_files), function(mod) {
    parse_primary_xlsx(primary_files[[mod]], mod)
  })
)

cat(sprintf("\nprimary_compiled: %d total rows, %d unique features, %d FDR-sig\n",
            nrow(primary_compiled),
            length(unique(primary_compiled$feature)),
            sum(primary_compiled$fdr_p_value < 0.05, na.rm = TRUE)))
cat(sprintf("Modalities: %s\n",
            paste(sort(unique(primary_compiled$modality)), collapse = ", ")))
cat(sprintf("Phenotypes: %s\n",
            paste(sort(unique(primary_compiled$phenotype)), collapse = ", ")))

# Validate expected feature counts per modality
cat("\nFeature counts per modality:\n")
mod_counts = primary_compiled[, .(n_features = length(unique(feature)),
                                  n_rows = .N,
                                  n_phenotypes = length(unique(phenotype)),
                                  n_sig = sum(fdr_p_value < 0.05, na.rm = TRUE)),
                              by = modality][order(modality)]
print(mod_counts)

expected_n = c(FA = 48, MD = 48, SA = 68, CT = 68, global_volume = 10, subcortical_volume = 14)
for (mod in names(expected_n)) {
  actual = mod_counts[modality == mod, n_features]
  if (length(actual) == 0 || actual != expected_n[mod]) {
    cat(sprintf("  *** VALIDATION FAIL: %s has %s features (expected %d) ***\n",
                mod, ifelse(length(actual) == 0, "0", as.character(actual)), expected_n[mod]))
  }
}
cat(sprintf("Total unique features: %d (expected 256)\n",
            length(unique(primary_compiled$feature))))

# Build feature-to-region lookup (used for table generation in Section 14)
feature_region_lookup = unique(primary_compiled[, .(feature, region_name, modality)])
feature_region_lookup = feature_region_lookup[!duplicated(feature)]

# Save compiled primary results for future use
saveRDS(primary_compiled, file.path(output_dir, "primary_results_compiled.RDS"))
saveRDS(feature_region_lookup, file.path(output_dir, "feature_region_lookup.RDS"))
cat(sprintf("Saved: %s\n", file.path(output_dir, "primary_results_compiled.RDS")))


# 9d. MERGE DIAGNOSTIC (verify primary vs sensitivity alignment)

cat("Diagnostic: Primary Analysis vs Sensitivity Analysis Results\n")

first_model = names(all_sensitivity_results)[1]
sens_compiled = readRDS(file.path(output_dir, sprintf("%s_compiled.RDS", first_model)))

# Feature overlap
p_features = unique(primary_compiled$feature)
s_features = unique(sens_compiled$feature)
shared_features = intersect(p_features, s_features)
p_only = setdiff(p_features, s_features)
s_only = setdiff(s_features, p_features)

cat(sprintf("Features: %d primary, %d sensitivity, %d shared\n",
            length(p_features), length(s_features), length(shared_features)))

if (length(p_only) > 0) {
  cat(sprintf("  In primary only (%d): %s\n", length(p_only),
              paste(head(sort(p_only), 10), collapse = ", ")))
  # Classify by modality
  p_only_mods = primary_compiled[feature %in% p_only, .(n = length(unique(feature))), by = modality]
  cat("  By modality:\n")
  for (i in seq_len(nrow(p_only_mods))) {
    cat(sprintf("    %s: %d features\n", p_only_mods$modality[i], p_only_mods$n[i]))
  }
}
if (length(s_only) > 0) {
  cat(sprintf("  In sensitivity only (%d): %s\n", length(s_only),
              paste(head(sort(s_only), 10), collapse = ", ")))
}

if (length(shared_features) == 0) {
  cat("*** Zero feature overlap. Check feature naming convention. ***\n")
  cat(sprintf("  Primary format (first 3): %s\n", paste(head(p_features, 3), collapse = ", ")))
  cat(sprintf("  Sensitivity format (first 3): %s\n", paste(head(s_features, 3), collapse = ", ")))
  stop("Cannot proceed with comparison: zero feature overlap.")
}

# Phenotype overlap
p_phenos = sort(unique(primary_compiled$phenotype))
s_phenos = sort(unique(sens_compiled$phenotype))
shared_phenos = intersect(p_phenos, s_phenos)

cat(sprintf("\nPhenotypes: %d primary, %d sensitivity, %d shared\n",
            length(p_phenos), length(s_phenos), length(shared_phenos)))
if (length(p_phenos) != length(s_phenos) || length(shared_phenos) != length(p_phenos)) {
  cat(sprintf("  Primary: %s\n", paste(p_phenos, collapse = ", ")))
  cat(sprintf("  Sensitivity: %s\n", paste(s_phenos, collapse = ", ")))
  cat(sprintf("  Missing from sensitivity: %s\n",
              paste(setdiff(p_phenos, s_phenos), collapse = ", ")))
}

# Modality overlap
p_mods = sort(unique(primary_compiled$modality))
s_mods = sort(unique(sens_compiled$modality))

cat(sprintf("\nModalities: primary (%s) / sensitivity (%s)\n",
            paste(p_mods, collapse = ", "), paste(s_mods, collapse = ", ")))
cat(sprintf("Shared: %d / %d\n",
            length(intersect(p_mods, s_mods)), length(union(p_mods, s_mods))))

# Primary FDR-significant by modality (these are what the comparison tests)
cat("\nPrimary FDR-significant associations (the denominator for comparison):\n")
print(primary_compiled[fdr_p_value < 0.05,
                       .(n_sig = .N,
                         n_features = length(unique(feature)),
                         n_phenotypes = length(unique(phenotype))),
                       by = modality][order(-n_sig)])

cat(sprintf("\nTotal primary FDR-significant: %d\n",
            sum(primary_compiled$fdr_p_value < 0.05, na.rm = TRUE)))
cat("\nMerge should proceed correctly.\n")


# 10. COMPARISON WITH PRIMARY RESULTS

cat("\n\n=== Comparing sensitivity results to primary analysis ===\n")

#' Compare sensitivity results to primary results
#'
#' For every feature x phenotype that was FDR-significant in the primary
#' analysis, reports: primary OR, sensitivity OR, whether it remains
#' significant, direction preservation, and attenuation.
#'
#' @param primary_df data.table with columns: feature, phenotype, odds_ratio,
#'   fdr_p_value, modality. Phenotype names must be canonical.
#' @param sensitivity_df data.table with columns: feature, phenotype,
#'   odds_ratio, fdr_p_value, model, modality. Phenotype names must be
#'   canonical (assigned by compile_sensitivity_results via assign_modality).
#' @return data.table with one row per primary FDR-significant association,
#'   annotated with sensitivity results and comparison metrics.
compare_to_primary = function(primary_df, sensitivity_df) {
  
  # Filter primary to FDR-significant only
  primary_sig = copy(primary_df[fdr_p_value < 0.05])
  
  if (nrow(primary_sig) == 0) {
    cat("  No FDR-significant primary results to compare.\n")
    return(data.table())
  }
  
  # Normalize phenotype names for merge (handles any residual case issues)
  primary_sig[, pheno_key := tolower(phenotype)]
  sens_copy = copy(sensitivity_df)
  sens_copy[, pheno_key := tolower(phenotype)]
  
  # Merge on feature + phenotype
  comparison = merge(
    primary_sig[, .(feature, pheno_key, phenotype, modality,
                    primary_OR = odds_ratio, primary_fdr_p = fdr_p_value)],
    sens_copy[, .(feature, pheno_key, model,
                  sensitivity_OR = odds_ratio, sensitivity_fdr_p = fdr_p_value)],
    by = c("feature", "pheno_key"),
    all.x = TRUE
  )
  
  # Report merge
  n_matched   = sum(!is.na(comparison$sensitivity_OR))
  n_unmatched = sum(is.na(comparison$sensitivity_OR))
  cat(sprintf("  Merge: %d matched, %d unmatched (primary features not in sensitivity)\n",
              n_matched, n_unmatched))
  
  if (n_unmatched > 0) {
    unmatched_feats = unique(comparison[is.na(sensitivity_OR), feature])
    unmatched_mods  = unique(comparison[is.na(sensitivity_OR), modality])
    cat(sprintf("  Unmatched features (%d unique): %s\n",
                length(unmatched_feats),
                paste(head(sort(unmatched_feats), 8), collapse = ", ")))
    cat(sprintf("  Unmatched modalities: %s\n", paste(unmatched_mods, collapse = ", ")))
  }
  
  # Compute comparison metrics (matched rows only)
  comparison[!is.na(sensitivity_OR), `:=`(
    remains_significant = sensitivity_fdr_p < 0.05,
    OR_change           = sensitivity_OR - primary_OR,
    OR_pct_change       = fifelse(primary_OR == 1, NA_real_,
                                  100 * (sensitivity_OR - primary_OR) / (primary_OR - 1)),
    direction_preserved = sign(log(sensitivity_OR)) == sign(log(primary_OR)),
    attenuated          = abs(log(sensitivity_OR)) < abs(log(primary_OR))
  )]
  
  # Normalize phenotype display names (use canonical from primary)
  comparison[, phenotype := pheno_canonical[pheno_key]]
  comparison[, pheno_key := NULL]
  
  return(comparison)
}


# --- Run comparison for each model ---

for (model_name in names(all_sensitivity_results)) {
  sens_path = file.path(output_dir, sprintf("%s_compiled.RDS", model_name))
  if (!file.exists(sens_path)) {
    cat(sprintf("\n  %s: compiled file not found, skipping.\n", model_name))
    next
  }
  sensitivity_compiled = readRDS(sens_path)
  
  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("=== %s ===\n", model_name))
  cat(sprintf("%s\n", paste(rep("=", 60), collapse = "")))
  
  comparison = compare_to_primary(primary_compiled, sensitivity_compiled)
  if (nrow(comparison) == 0) next
  
  matched = comparison[!is.na(sensitivity_OR)]
  n_primary_sig = nrow(matched)
  
  if (n_primary_sig == 0) {
    cat("  No matched features to compare.\n")
    next
  }
  
  n_still_sig  = sum(matched$remains_significant, na.rm = TRUE)
  n_attenuated = sum(matched$attenuated, na.rm = TRUE)
  n_direction  = sum(matched$direction_preserved, na.rm = TRUE)
  n_lost       = n_primary_sig - n_still_sig
  
  # --- Overall summary ---
  cat(sprintf("\n  OVERALL (matched features only):\n"))
  cat(sprintf("  Primary FDR-significant associations matched: %d\n", n_primary_sig))
  cat(sprintf("  Remain FDR-significant after adjustment: %d (%.1f%%)\n",
              n_still_sig, 100 * n_still_sig / n_primary_sig))
  cat(sprintf("  Attenuated (OR moved toward 1): %d (%.1f%%)\n",
              n_attenuated, 100 * n_attenuated / n_primary_sig))
  cat(sprintf("  Direction preserved: %d (%.1f%%)\n",
              n_direction, 100 * n_direction / n_primary_sig))
  cat(sprintf("  No longer significant: %d (%.1f%%)\n",
              n_lost, 100 * n_lost / n_primary_sig))
  
  # Median absolute OR change
  or_changes = matched[, abs(OR_change)]
  cat(sprintf("  Median |OR change|: %.4f (IQR: %.4f-%.4f)\n",
              median(or_changes, na.rm = TRUE),
              quantile(or_changes, 0.25, na.rm = TRUE),
              quantile(or_changes, 0.75, na.rm = TRUE)))
  
  # Direction flips summary (if any)
  n_flipped = n_primary_sig - n_direction
  if (n_flipped > 0) {
    flipped = matched[direction_preserved == FALSE]
    flip_by_mod   = flipped[, .N, by = modality][order(-N)]
    flip_by_pheno = flipped[, .N, by = phenotype][order(-N)]
    cat(sprintf("\n  DIRECTION FLIPS (%d associations):\n", n_flipped))
    cat("    By modality: ")
    cat(paste(sprintf("%s=%d", flip_by_mod$modality, flip_by_mod$N), collapse = ", "))
    cat("\n    By phenotype: ")
    cat(paste(sprintf("%s=%d", flip_by_pheno$phenotype, flip_by_pheno$N), collapse = ", "))
    cat("\n")
  }
  
  # --- Breakdown by each IDP ---
  cat("\n  Breakdown by each IDP:\n")
  for (mod in sort(unique(na.omit(matched$modality)))) {
    mod_matched = matched[modality == mod]
    n_mod = nrow(mod_matched)
    if (n_mod == 0) next
    n_sig_mod  = sum(mod_matched$remains_significant, na.rm = TRUE)
    n_dir_mod  = sum(mod_matched$direction_preserved, na.rm = TRUE)
    n_att_mod  = sum(mod_matched$attenuated, na.rm = TRUE)
    cat(sprintf("    %s: %d matched, %d (%.0f%%) remain sig, %d (%.0f%%) direction preserved, %d (%.0f%%) attenuated\n",
                mod, n_mod, n_sig_mod, 100 * n_sig_mod / n_mod,
                n_dir_mod, 100 * n_dir_mod / n_mod,
                n_att_mod, 100 * n_att_mod / n_mod))
  }
  
  # Report IDPs with 0 primary-sig (not in the comparison at all)
  all_mods_in_primary = unique(primary_compiled[fdr_p_value < 0.05, modality])
  mods_in_comparison  = unique(matched$modality)
  missing_mods = setdiff(all_mods_in_primary, mods_in_comparison)
  if (length(missing_mods) > 0) {
    cat(sprintf("    (No matched features for: %s)\n",
                paste(missing_mods, collapse = ", ")))
  }
  
  # --- Per-phenotype breakdown ---
  cat("\n  Per-phenotype breakdown:\n")
  for (pheno in sort(unique(comparison$phenotype))) {
    ph_matched = matched[phenotype == pheno]
    n_ph = nrow(ph_matched)
    if (n_ph == 0) {
      cat(sprintf("    %s: 0 matched features\n", pheno))
      next
    }
    n_sig_ph = sum(ph_matched$remains_significant, na.rm = TRUE)
    n_dir_ph = sum(ph_matched$direction_preserved, na.rm = TRUE)
    n_att_ph = sum(ph_matched$attenuated, na.rm = TRUE)
    cat(sprintf("    %s: %d matched, %d (%.0f%%) remain sig, %d (%.0f%%) direction preserved, %d (%.0f%%) attenuated\n",
                pheno, n_ph, n_sig_ph, 100 * n_sig_ph / n_ph,
                n_dir_ph, 100 * n_dir_ph / n_ph,
                n_att_ph, 100 * n_att_ph / n_ph))
  }
  
  # Phenotypes with 0 primary-sig (not in comparison)
  all_phenos_in_primary = unique(primary_compiled[fdr_p_value < 0.05, phenotype])
  phenos_in_comparison  = unique(comparison$phenotype)
  missing_phenos = setdiff(all_phenos_in_primary, phenos_in_comparison)
  if (length(missing_phenos) > 0) {
    cat(sprintf("    (No primary FDR-sig associations for: %s)\n",
                paste(missing_phenos, collapse = ", ")))
  }
  
  # --- Save ---
  saveRDS(comparison, file.path(output_dir, sprintf("comparison_%s.RDS", model_name)))
  fwrite(comparison, file.path(output_dir, sprintf("comparison_%s.csv", model_name)))
  cat(sprintf("\n  Saved: comparison_%s.RDS / .csv\n", model_name))
}

comp_m1 = readRDS(file.path(output_dir, "comparison_Model_1.RDS"))
matched_m1 = subset(comp_m1, !is.na(sensitivity_OR))
or_changes = abs(matched_m1$OR_change)
cat(sprintf("Median |OR change|: %.4f (IQR: %.4f-%.4f)\n",
            median(or_changes, na.rm = TRUE),
            quantile(or_changes, 0.25, na.rm = TRUE),
            quantile(or_changes, 0.75, na.rm = TRUE)))

# 11. MODEL DOCUMENTATION AND COMPARISON SUMMARY

for (model_name in names(model_specs)) {
  spec = model_specs[[model_name]]
  cat(sprintf("--- %s (%s) ---\n", model_name, spec$label))
  
  # Model specification
  cat(sprintf("  Covariates added: %s\n", paste(spec$new_covariates, collapse = ", ")))
  cat(sprintf("  MHQ timing exclusion: %s\n",
              ifelse(spec$requires_mhq_exclusion, "Yes", "No")))
  cat(sprintf("  FDR correction: Within each phenotype, matching primary analysis\n"))
  
  # Phenotypes included/excluded
  all_phenos = names(phenotype_registry)
  included = spec$phenotypes
  excluded = setdiff(all_phenos, included)
  cat(sprintf("  Phenotypes included (%d): %s\n",
              length(included), paste(included, collapse = ", ")))
  if (length(excluded) > 0) {
    cat(sprintf("  Phenotypes excluded (%d): %s\n",
                length(excluded), paste(excluded, collapse = ", ")))
    if (model_name == "Model_3") {
      cat("    Reason: PHQ-9 restricted to MHQ-derived phenotypes (LifetimeMDD, MDDRecur)\n")
      cat("    Shallow phenotypes have <50% MHQ completion with valid timing.\n")
    }
  }
  
  # Comparison summary (if comparison file exists)
  comp_path = file.path(output_dir, sprintf("comparison_%s.RDS", model_name))
  if (file.exists(comp_path)) {
    comp = readRDS(comp_path)
    matched = comp[!is.na(sensitivity_OR)]
    n_total = nrow(matched)
    if (n_total > 0) {
      n_sig = sum(matched$remains_significant, na.rm = TRUE)
      n_dir = sum(matched$direction_preserved, na.rm = TRUE)
      n_att = sum(matched$attenuated, na.rm = TRUE)
      cat(sprintf("  Comparison: %d matched, %d (%.1f%%) remain sig, %d (%.1f%%) direction preserved, %d (%.1f%%) attenuated\n",
                  n_total, n_sig, 100 * n_sig / n_total,
                  n_dir, 100 * n_dir / n_total,
                  n_att, 100 * n_att / n_total))
      # Flag direction flips if any
      n_flipped = n_total - n_dir
      if (n_flipped > 0) {
        flip_summary = matched[direction_preserved == FALSE, .N, by = .(modality, phenotype)]
        cat(sprintf("  Direction flips (%d): %s\n", n_flipped,
                    paste(sprintf("%s x %s (n=%d)",
                                  flip_summary$modality,
                                  flip_summary$phenotype,
                                  flip_summary$N),
                          collapse = "; ")))
      }
    }
  } else {
    cat("  Comparison: not yet run\n")
  }
  
  # Interpretation notes
  if (model_name == "Model_1") {
    cat("  Note: Antidepressant use is a plausible mediator. Adjustment attenuates\n")
    cat("    total effect of depression on brain structure; this is expected.\n")
  }
  if (model_name == "Model_2") {
    cat("  Note: Anxiety proxy (field 20428) overlaps conceptually with GPPsy,\n")
    cat("    PsyPsy, and GPNoDep definitions. MHQ non-completion causes severe\n")
    cat("    sample loss for shallow phenotypes.\n")
  }
  if (model_name == "Model_3") {
    cat("  Note: PHQ-9 adjustment tests whether associations survive control for\n")
    cat("    current symptom severity. Complete attenuation is consistent with\n")
    cat("    primary associations partly reflecting symptom burden at MHQ.\n")
  }
  
  cat("\n")
}

# 12. EXPORT SAMPLE SIZE TABLE (Excel, one tab per phenotype)

cat("\n=== Generating sample size summary workbook ===\n")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}
library(openxlsx)

# Print the full tracker for the log
cat("\nSample size tracker (all model x phenotype combinations):\n")
print(sample_size_tracker)

# Save tracker as RDS and CSV for programmatic use
saveRDS(sample_size_tracker, file = file.path(output_dir, "sensitivity_sample_sizes.RDS"))
fwrite(sample_size_tracker, file = file.path(output_dir, "sensitivity_sample_sizes.csv"))

# --- Build Excel workbook: one tab per phenotype ---

wb_ss = createWorkbook()

# Styles
ss_header = createStyle(
  fontName = "Times New Roman", fontSize = 12, fontColour = "#FFFFFF",
  fgFill = "#4472C4", halign = "center", textDecoration = "bold",
  border = "Bottom", borderColour = "#2F5496"
)
ss_header_left = createStyle(
  fontName = "Times New Roman", fontSize = 12, fontColour = "#FFFFFF",
  fgFill = "#4472C4", halign = "left", textDecoration = "bold",
  border = "Bottom", borderColour = "#2F5496"
)
ss_title = createStyle(
  fontName = "Times New Roman", fontSize = 12, textDecoration = "bold", halign = "center"
)
ss_data    = createStyle(fontName = "Times New Roman", fontSize = 12, halign = "center")
ss_data_l  = createStyle(fontName = "Times New Roman", fontSize = 12, halign = "left")
ss_note    = createStyle(fontName = "Times New Roman", fontSize = 9, textDecoration = "italic")
ss_skip    = createStyle(fontName = "Times New Roman", fontSize = 10, halign = "center",
                         fontColour = "#C00000")
ss_border  = createStyle(border = "Bottom", borderColour = "#D9D9D9", borderStyle = "thin")

all_phenotypes = names(phenotype_registry)
all_models     = names(model_specs)

for (pheno_name in all_phenotypes) {
  
  addWorksheet(wb_ss, sheetName = pheno_name)
  
  # Title
  mergeCells(wb_ss, pheno_name, cols = 1:7, rows = 1)
  writeData(wb_ss, pheno_name,
            x = sprintf("Sensitivity Analysis Sample Sizes: %s", pheno_name),
            startCol = 1, startRow = 1)
  addStyle(wb_ss, pheno_name, ss_title, rows = 1, cols = 1:7)
  
  # Headers
  headers = c("Model", "Covariates Added", "N (Original)", "N (Post-MHQ)",
              "N (Final)", "Cases", "Controls")
  writeData(wb_ss, pheno_name, x = t(headers), startCol = 1, startRow = 3, colNames = FALSE)
  addStyle(wb_ss, pheno_name, ss_header_left, rows = 3, cols = 1:2)
  addStyle(wb_ss, pheno_name, ss_header, rows = 3, cols = 3:7)
  
  current_row = 4
  
  for (model_name in all_models) {
    spec = model_specs[[model_name]]
    
    # Check if this phenotype is in this model's list
    if (!pheno_name %in% spec$phenotypes) {
      # Model not applicable
      writeData(wb_ss, pheno_name, x = spec$label, startCol = 1, startRow = current_row)
      addStyle(wb_ss, pheno_name, ss_data_l, rows = current_row, cols = 1)
      writeData(wb_ss, pheno_name, x = paste(spec$new_covariates, collapse = " + "),
                startCol = 2, startRow = current_row)
      addStyle(wb_ss, pheno_name, ss_data_l, rows = current_row, cols = 2)
      
      mergeCells(wb_ss, pheno_name, cols = 3:7, rows = current_row)
      writeData(wb_ss, pheno_name,
                x = sprintf("Not applicable (restricted to %s)",
                            paste(spec$phenotypes, collapse = ", ")),
                startCol = 3, startRow = current_row)
      addStyle(wb_ss, pheno_name, ss_skip, rows = current_row, cols = 3:7)
      addStyle(wb_ss, pheno_name, ss_border, rows = current_row, cols = 1:7, stack = TRUE)
      current_row = current_row + 1
      next
    }
    
    # Look up tracker row
    tracker_row = sample_size_tracker[model == model_name & phenotype == pheno_name]
    
    if (nrow(tracker_row) == 0) {
      writeData(wb_ss, pheno_name, x = spec$label, startCol = 1, startRow = current_row)
      addStyle(wb_ss, pheno_name, ss_data_l, rows = current_row, cols = 1)
      mergeCells(wb_ss, pheno_name, cols = 2:7, rows = current_row)
      writeData(wb_ss, pheno_name, x = "No data (model may not have run)",
                startCol = 2, startRow = current_row)
      addStyle(wb_ss, pheno_name, ss_skip, rows = current_row, cols = 2:7)
      addStyle(wb_ss, pheno_name, ss_border, rows = current_row, cols = 1:7, stack = TRUE)
      current_row = current_row + 1
      next
    }
    
    tr = tracker_row[1]
    
    # Model label
    writeData(wb_ss, pheno_name, x = tr$model_label, startCol = 1, startRow = current_row)
    addStyle(wb_ss, pheno_name, ss_data_l, rows = current_row, cols = 1)
    
    # Covariates
    writeData(wb_ss, pheno_name, x = paste(spec$new_covariates, collapse = " + "),
              startCol = 2, startRow = current_row)
    addStyle(wb_ss, pheno_name, ss_data_l, rows = current_row, cols = 2)
    
    # N original
    writeData(wb_ss, pheno_name, x = format(tr$n_original, big.mark = ","),
              startCol = 3, startRow = current_row)
    
    # N post-MHQ
    post_mhq_val = ifelse(spec$requires_mhq_exclusion,
                          format(tr$n_post_mhq, big.mark = ","),
                          "\u2014")
    writeData(wb_ss, pheno_name, x = post_mhq_val, startCol = 4, startRow = current_row)
    
    # N final, cases, controls (or skipped)
    if (tr$status == "completed") {
      writeData(wb_ss, pheno_name, x = format(tr$n_post_cc, big.mark = ","),
                startCol = 5, startRow = current_row)
      writeData(wb_ss, pheno_name, x = format(tr$n_cases, big.mark = ","),
                startCol = 6, startRow = current_row)
      writeData(wb_ss, pheno_name, x = format(tr$n_controls, big.mark = ","),
                startCol = 7, startRow = current_row)
      addStyle(wb_ss, pheno_name, ss_data, rows = current_row, cols = 3:7)
    } else {
      skip_label = switch(tr$status,
                          skipped_missingness = sprintf("Skipped (%.0f%% missing)", tr$pct_dropped),
                          skipped_error       = "Skipped (error)",
                          skipped_low_n       = sprintf("Skipped (<50 cases; n=%s)",
                                                        format(tr$n_cases, big.mark = ",")),
                          sprintf("Skipped (%s)", tr$status)
      )
      if (!is.na(tr$n_post_cc)) {
        writeData(wb_ss, pheno_name, x = format(tr$n_post_cc, big.mark = ","),
                  startCol = 5, startRow = current_row)
      }
      if (!is.na(tr$n_cases)) {
        writeData(wb_ss, pheno_name, x = format(tr$n_cases, big.mark = ","),
                  startCol = 6, startRow = current_row)
      }
      if (!is.na(tr$n_controls)) {
        writeData(wb_ss, pheno_name, x = format(tr$n_controls, big.mark = ","),
                  startCol = 7, startRow = current_row)
      }
      addStyle(wb_ss, pheno_name, ss_skip, rows = current_row, cols = 3:7)
      
      # Skip note row
      current_row = current_row + 1
      mergeCells(wb_ss, pheno_name, cols = 1:7, rows = current_row)
      writeData(wb_ss, pheno_name, x = skip_label, startCol = 1, startRow = current_row)
      addStyle(wb_ss, pheno_name, ss_note, rows = current_row, cols = 1)
    }
    
    addStyle(wb_ss, pheno_name, ss_border, rows = current_row, cols = 1:7, stack = TRUE)
    current_row = current_row + 1
  }
  
  # Footnotes
  fn_start = current_row + 1
  footnotes = c(
    "N (Original): Sample size in the z_score_completed object for this phenotype.",
    "N (Post-MHQ): After excluding participants whose MRI preceded the MHQ (\u2014 = no MHQ exclusion for this model).",
    "N (Final): After complete-case restriction on the new covariate(s).",
    "Model 3 is restricted to LifetimeMDD and MDDRecur only (see Supplementary Methods).",
    "Skipped models indicate >90% sample loss, <50 cases, or model fitting errors."
  )
  for (j in seq_along(footnotes)) {
    writeData(wb_ss, pheno_name, x = footnotes[j], startCol = 1, startRow = fn_start + j - 1)
    addStyle(wb_ss, pheno_name, ss_note, rows = fn_start + j - 1, cols = 1)
  }
  
  # Column widths
  setColWidths(wb_ss, pheno_name, cols = 1, widths = 40)
  setColWidths(wb_ss, pheno_name, cols = 2, widths = 30)
  setColWidths(wb_ss, pheno_name, cols = 3:7, widths = 14)
}

xlsx_path = file.path(output_dir, "Sensitivity_Sample_Sizes_New_Covariates.xlsx")
saveWorkbook(wb_ss, xlsx_path, overwrite = TRUE)
cat(sprintf("Saved sample size workbook: %s\n", xlsx_path))


# 14. SUPPLEMENTARY TABLE GENERATION (Models 1-4)
#
# Generates one .xlsx per model with tabs by modality and a Sample_Sizes tab.
# Within each modality tab, phenotypes are stacked with blank-row separators.
# Columns include sensitivity OR, CIs, p-values, primary OR for comparison,
# and direction_preserved / remains_significant flags.
#
# Output:
#   Supplementary_Table_Da_Antidepressant.xlsx   (Model 1)
#   Supplementary_Table_Db_Anxiety.xlsx          (Model 2)
#   Supplementary_Table_Dc_PHQ9.xlsx             (Model 3)
#   Supplementary_Table_Dd_Psychotherapy.xlsx    (Model 4)
#
# Within-cases tables (episode count, MDD age of onset) are generated in  
# S4b_sensitivity_within_cases.R, shown in Supplementary Tables S4e-S4f.

cat("\n\n=== Section 14: Generating Supplementary Tables Da-Dd ===\n")

phenotype_order = c("GPNoDep", "GPPsy", "PsyPsy", "SelfRepDep", "DepAll",
                    "ICD10Dep", "ICD10Dep_exclpsych", "LifetimeMDD", "MDDRecur")

modality_order = c("FA", "MD", "SA", "CT", "global_volume", "subcortical_volume")

modality_tab_names = c(
  FA                 = "FA",
  MD                 = "MD",
  SA                 = "SA",
  CT                 = "CT",
  global_volume      = "Global_Volume",
  subcortical_volume = "Subcortical_Volume"
)

# --- Styles ---

st14_header = createStyle(
  fontName = "Times New Roman", fontSize = 11, textDecoration = "bold",
  fgFill = "#4472C4", fontColour = "#FFFFFF",
  halign = "center", border = "TopBottom", borderColour = "#2F5496"
)

st14_pheno_label = createStyle(
  fontName = "Times New Roman", fontSize = 11, textDecoration = "bold",
  fgFill = "#D9E2F3", border = "Bottom", borderColour = "#4472C4"
)

st14_sig = createStyle(
  fontName = "Times New Roman", fontSize = 10, textDecoration = "bold"
)

# --- Load lookup if needed ---
if (!exists("feature_region_lookup")) {
  feature_region_lookup = readRDS(file.path(output_dir, "feature_region_lookup.RDS"))
}

# --- Helper: write one modality tab ---

write_modality_tab = function(wb, tab_name, sens_dt, prim_dt,
                              mod_name, phenotypes, lookup) {
  
  addWorksheet(wb, tab_name)
  
  col_headers = c("Feature", "Region", "Sensitivity OR", "Lower CI", "Upper CI",
                  "p-value", "FDR p-value", "Primary OR",
                  "Direction Preserved", "Remains Significant")
  
  current_row = 1
  
  for (pheno in phenotypes) {
    sens_sub = sens_dt[modality == mod_name & phenotype == pheno]
    if (nrow(sens_sub) == 0) next
    
    prim_sub = prim_dt[modality == mod_name & phenotype == pheno,
                       .(feature, primary_OR = odds_ratio, primary_fdr_p = fdr_p_value)]
    
    merged = merge(sens_sub, prim_sub, by = "feature", all.x = TRUE)
    merged = merge(merged, lookup[, .(feature, region_name)], by = "feature", all.x = TRUE)
    
    merged[, direction_preserved := fifelse(
      is.na(primary_OR), NA,
      sign(log(odds_ratio)) == sign(log(primary_OR))
    )]
    merged[, remains_sig := fdr_p_value < 0.05]
    setorder(merged, feature)
    
    n_cases    = merged$n_cases_final[1]
    n_controls = merged$n_controls_final[1]
    n_total    = merged$n_total_final[1]
    
    # Phenotype label row
    label_text = sprintf("%s (n_cases=%s, n_controls=%s, n_total=%s)",
                         pheno,
                         format(n_cases, big.mark = ","),
                         format(n_controls, big.mark = ","),
                         format(n_total, big.mark = ","))
    writeData(wb, tab_name, x = label_text,
              startRow = current_row, startCol = 1)
    mergeCells(wb, tab_name, cols = 1:length(col_headers), rows = current_row)
    addStyle(wb, tab_name, st14_pheno_label, rows = current_row,
             cols = 1:length(col_headers), gridExpand = TRUE)
    current_row = current_row + 1
    
    # Column headers
    writeData(wb, tab_name, x = as.data.frame(t(col_headers)),
              startRow = current_row, startCol = 1, colNames = FALSE)
    addStyle(wb, tab_name, st14_header, rows = current_row,
             cols = 1:length(col_headers), gridExpand = TRUE)
    current_row = current_row + 1
    
    # Data rows
    out_df = data.frame(
      feature     = merged$feature,
      region      = fifelse(is.na(merged$region_name), "", merged$region_name),
      sens_or     = merged$odds_ratio,
      lower_ci    = merged$lower_CI,
      upper_ci    = merged$upper_CI,
      p_val       = merged$p_value,
      fdr_p       = merged$fdr_p_value,
      primary_or  = merged$primary_OR,
      dir_pres    = fifelse(is.na(merged$direction_preserved), "",
                            fifelse(merged$direction_preserved, "TRUE", "FALSE")),
      remains_sig = fifelse(merged$remains_sig, "TRUE", "FALSE"),
      stringsAsFactors = FALSE
    )
    writeData(wb, tab_name, x = out_df,
              startRow = current_row, startCol = 1, colNames = FALSE)
    
    # Highlight FDR-significant rows
    sig_rows = which(merged$fdr_p_value < 0.05)
    if (length(sig_rows) > 0) {
      addStyle(wb, tab_name, st14_sig,
               rows = current_row + sig_rows - 1,
               cols = 3:7, gridExpand = TRUE, stack = TRUE)
    }
    
    current_row = current_row + nrow(merged) + 1  # +1 for blank separator
  }
  
  # Column widths
  setColWidths(wb, tab_name, cols = 1, widths = 14)
  setColWidths(wb, tab_name, cols = 2, widths = 55)
  setColWidths(wb, tab_name, cols = 3:8, widths = 14)
  setColWidths(wb, tab_name, cols = 9:10, widths = 18)
}

# --- Helper: write Sample_Sizes tab ---

write_ss_tab = function(wb, tracker_sub, model_label) {
  addWorksheet(wb, "Sample_Sizes")
  
  writeData(wb, "Sample_Sizes",
            x = sprintf("Sample sizes: %s", model_label),
            startRow = 1, startCol = 1)
  addStyle(wb, "Sample_Sizes",
           createStyle(fontSize = 12, textDecoration = "bold"),
           rows = 1, cols = 1)
  
  ss_cols = c("Phenotype", "N Original", "N Post-MHQ", "N Post-CC",
              "N Cases", "N Controls", "% Dropped", "Status")
  writeData(wb, "Sample_Sizes", x = as.data.frame(t(ss_cols)),
            startRow = 3, startCol = 1, colNames = FALSE)
  addStyle(wb, "Sample_Sizes", st14_header, rows = 3,
           cols = 1:length(ss_cols), gridExpand = TRUE)
  
  for (i in seq_len(nrow(tracker_sub))) {
    r = tracker_sub[i]
    row_data = data.frame(
      pheno  = r$phenotype,
      n_orig = format(r$n_original, big.mark = ","),
      n_mhq  = fifelse(is.na(r$n_post_mhq), "\u2014",
                       format(r$n_post_mhq, big.mark = ",")),
      n_cc   = fifelse(is.na(r$n_post_cc), "\u2014",
                       format(r$n_post_cc, big.mark = ",")),
      n_case = fifelse(is.na(r$n_cases), "\u2014",
                       format(r$n_cases, big.mark = ",")),
      n_ctrl = fifelse(is.na(r$n_controls), "\u2014",
                       format(r$n_controls, big.mark = ",")),
      pct    = fifelse(is.na(r$pct_dropped), "\u2014",
                       sprintf("%.1f%%", r$pct_dropped)),
      status = r$status,
      stringsAsFactors = FALSE
    )
    writeData(wb, "Sample_Sizes", x = row_data,
              startRow = 3 + i, startCol = 1, colNames = FALSE)
  }
  
  setColWidths(wb, "Sample_Sizes", cols = 1, widths = 22)
  setColWidths(wb, "Sample_Sizes", cols = 2:7, widths = 14)
  setColWidths(wb, "Sample_Sizes", cols = 8, widths = 20)
}


# --- Generate workbooks ---

model_table_map = list(
  list(model_key = "Model_1",
       filename  = "Supplementary_Table_Da_Antidepressant.xlsx",
       label     = "Model 1: Primary + Antidepressant Use"),
  list(model_key = "Model_2",
       filename  = "Supplementary_Table_Db_Anxiety.xlsx",
       label     = "Model 2: Primary + Anxiety Proxy"),
  list(model_key = "Model_3",
       filename  = "Supplementary_Table_Dc_PHQ9.xlsx",
       label     = "Model 3: Primary + PHQ-9 Severity"),
  list(model_key = "Model_4",
       filename  = "Supplementary_Table_Dd_Psychotherapy.xlsx",
       label     = "Model 4: Primary + Psychotherapy")
)

for (tbl in model_table_map) {
  
  cat(sprintf("\n--- Generating %s ---\n", tbl$filename))
  
  sens_path = file.path(output_dir, sprintf("%s_compiled.RDS", tbl$model_key))
  if (!file.exists(sens_path)) {
    cat(sprintf("  SKIPPED: %s not found\n", sens_path))
    next
  }
  sens_compiled = readRDS(sens_path)
  model_phenos  = intersect(phenotype_order, unique(sens_compiled$phenotype))
  
  wb = createWorkbook()
  
  for (mod in modality_order) {
    tab_name = modality_tab_names[mod]
    if (nrow(sens_compiled[modality == mod]) == 0) next
    
    write_modality_tab(
      wb         = wb,
      tab_name   = tab_name,
      sens_dt    = sens_compiled,
      prim_dt    = primary_compiled,
      mod_name   = mod,
      phenotypes = model_phenos,
      lookup     = feature_region_lookup
    )
    cat(sprintf("  %s: written\n", tab_name))
  }
  
  tracker_sub = sample_size_tracker[model == tbl$model_key]
  if (nrow(tracker_sub) > 0) {
    write_ss_tab(wb, tracker_sub, tbl$label)
    cat("  Sample_Sizes: written\n")
  }
  
  out_path = file.path(output_dir, tbl$filename)
  saveWorkbook(wb, out_path, overwrite = TRUE)
  cat(sprintf("  Saved: %s\n", out_path))
}

cat("Done. Proceed to S4b\n")
