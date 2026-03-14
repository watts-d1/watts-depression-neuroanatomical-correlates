# Supplementary Analysis: Symptom Dimension x Brain Structure
# Within-Cases Analysis for Molecular Psychiatry Revision
#
# PURPOSE:"Can alterations within specific white matter tracts 
#   be mapped onto discrete clinical dimensions of depression?"
#
# DESIGN:
#   Within-cases logistic regression. Among individuals meeting criteria 
#   for each depression phenotype, test whether each IDP (z-scored) 
#   predicts endorsement of specific symptom dimensions.
#
#   Outcome: Symptom dimension (0/1) — NOT depression case status
#   Predictors: Each IDP, one at a time, adjusting for covariates
#   Sample: Cases only (phenotype == 1) from each definition
#
# CONSTRUCTS (from MHQ Category 138, historical/lifetime items):
#   1. Cognitive slowing:  Field 20435 — "Difficulty concentrating during worst depression"
#   2. Anhedonia:          Field 20441 — "Ever had prolonged loss of interest in normal activities"
#   3. Core dysphoria:     Field 20446 — "Ever had prolonged feelings of sadness or depression"
#
# PREREQUISITES:
#   - csv_data loaded (UK Biobank main CSV via data.table::fread)
#   - z_score objects loaded from completed_ICD_objects_z_score.RData
#   - Packages: data.table, dplyr, future, future.apply, parallelly, MASS
#   - Source the function definitions below OR ensure they are in the environment
#     from S1_neuroimaging_feature_extraction.R
#
# NOTES:
#   - These MHQ fields are conditional: only asked of participants who endorsed 
#     depression in the online questionnaire. Therefore this analysis is inherently 
#     within-cases. Participants who did not complete the MHQ or did not endorse 
#     depression screening items will have NA for these fields.
#   - We use instance 0.0 for MHQ fields (online follow-up, no instance indexing 
#     by visit type for these fields — they are singular instance).
#   - Results are run for ALL IDPs but the manuscript will focus reporting on 
#     white matter tracts (FA, MD) as requested by the reviewer.

# --- User-defined paths (modify for your environment) ---
UKBB_DATA_PATH = "/path/to/ukbiobank/ukb674571.csv"  # UK Biobank phenotype CSV
S1_FUNCTIONS_PATH = "/path/to/S1_functions.R"  # Utility functions from S1 preprocessing
S1_Z_SCORE_PATH = "/path/to/z_score_neuroimaging_objects.RData"  # S1 output: z-scored neuroimaging data

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

packages = c('dplyr', 'data.table', 'magrittr', 'rsample', 'readr', 'stringr', 
             'car', 'future', 'future.apply', 'parallel', 'furrr', 'openxlsx')

package_install(packages)
library(here)

# 0. Sourcing S1_functions.R
source(S1_FUNCTIONS_PATH)

# 0. Configuration

# Paths — adjust as needed
z_score_path   = S1_Z_SCORE_PATH
main_script    = here("regression_analyses", "S1_neuroimaging_feature_extraction.R")
output_dir     = here("output", "symptom_dimension_results")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1. Load z_score objects
load(z_score_path)

# 2. Source functions from existing pipeline
# All UKB fields: https://biobank.ctsu.ox.ac.uk/showcase/label.cgi?id=100003
csv_file_path = UKBB_DATA_PATH
csv_data = fread(csv_file_path, showProgress = TRUE, nThread = 20)

# 2b. Measurement ID vectors (needed for results filtering)
#     Copy from main script if not already in environment

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
# 3. Extract symptom dimension constructs from csv_data

# --- 3a. Primary constructs (historical/lifetime items) ---
primary_fields = c("20435-0.0",   # Concentration difficulties: Difficulty concentrating during worst depression
                   "20441-0.0",   # Anhedonia: Ever had prolonged loss of interest in normal activities
                   "20450-0.0")   # Negative self-evaluation: Feelings of worthlessness during worst depression

# --- 3b. Sensitivity constructs (current/recent items, PHQ-style) ---
sensitivity_fields = c("20508-0.0",   # Concentration difficulties: Recent trouble concentrating on things
                       "20514-0.0",   # Anhedonia: Recent lack of interest or pleasure in doing things
                       "20507-0.0")   # Negative self-evaluation: Recent feelings of inadequacy

# Note: No current/recent item available for core dysphoria;
# sensitivity analysis covers 3 of 4 primary constructs.
all_construct_fields = c(primary_fields, sensitivity_fields)

# Verify fields exist in csv_data
missing_fields = setdiff(all_construct_fields, names(csv_data))
if (length(missing_fields) > 0) {
  stop(paste("Fields not found in csv_data:", paste(missing_fields, collapse = ", "),
             "\nCheck field naming convention (e.g., try without '-0.0' suffix or with 'p' prefix)."))
}

# Extract constructs with subject identifier
id_col = ifelse("subject_num" %in% names(csv_data), "subject_num", "eid")

construct_data = csv_data[, c(id_col, all_construct_fields), with = FALSE]

# Standardize column names
setnames(construct_data, 
         old = c(id_col, primary_fields, sensitivity_fields),
         new = c("subject_num", 
                 "cognitive_slowing_raw", "anhedonia_raw", "dysphoria_raw", "neg_self_eval_raw",
                 "cognitive_slowing_recent_raw", "anhedonia_recent_raw", "neg_self_eval_recent_raw"))

construct_data[, subject_num := as.integer(subject_num)]

cat("\n=== Raw construct distributions (before recoding) ===\n")
cat("\n--- Primary (historical/lifetime) ---\n")
cat("\nCognitive slowing (20435):\n");        print(table(construct_data$cognitive_slowing_raw, useNA = "always"))
cat("\nAnhedonia (20441):\n");                 print(table(construct_data$anhedonia_raw, useNA = "always"))
cat("\nCore dysphoria (20446):\n");            print(table(construct_data$dysphoria_raw, useNA = "always"))
cat("\nNegative self-evaluation (20450):\n");  print(table(construct_data$neg_self_eval_raw, useNA = "always"))
cat("\n--- Sensitivity (current/recent) ---\n")
cat("\nCognitive slowing recent (20508):\n");        print(table(construct_data$cognitive_slowing_recent_raw, useNA = "always"))
cat("\nAnhedonia recent (20514):\n");                 print(table(construct_data$anhedonia_recent_raw, useNA = "always"))
cat("\nNegative self-evaluation recent (20507):\n");  print(table(construct_data$neg_self_eval_recent_raw, useNA = "always"))

# 4. Recode constructs to binary
#
# --- Primary constructs (historical/lifetime) ---
# UK Biobank MHQ coding:
#   20435: 1 = Yes, 0 = No, -121 = Do not know, -818 = Prefer not to answer
#   20441: 1 = Yes, 0 = No, -818 = Prefer not to answer
#   20446: 1 = Yes, 0 = No, -818 = Prefer not to answer
#   20450: 1 = Yes, 0 = No, -121 = Do not know, -818 = Prefer not to answer

recode_binary = function(x) {
  x = as.numeric(x)
  x[x < 0] = NA
  x[!x %in% c(0, 1)] = NA
  return(as.integer(x))
}

recode_phq_binary = function(x) {
  # Dichotomize PHQ-style items: "Not at all" (1) -> 0, any endorsement (2/3/4) -> 1
  x = as.numeric(x)
  x[x < 0] = NA            # Remove PNA
  out = rep(NA_integer_, length(x))
  out[x == 1] = 0L          # Not at all
  out[x %in% c(2, 3, 4)] = 1L  # Any endorsement
  return(out)
}

# Primary constructs
construct_data[, cognitive_slowing := recode_binary(cognitive_slowing_raw)]
construct_data[, anhedonia         := recode_binary(anhedonia_raw)]
construct_data[, dysphoria         := recode_binary(dysphoria_raw)]
construct_data[, neg_self_eval     := recode_binary(neg_self_eval_raw)]

# Drop raw columns
raw_cols = grep("_raw$", names(construct_data), value = TRUE)
construct_data[, (raw_cols) := NULL]

cat("\n=== Recoded construct distributions ===\n")
cat("\n--- Primary (historical/lifetime): Yes / No / NA ---\n")
cat("\nCognitive slowing (20435):\n");        print(table(construct_data$cognitive_slowing, useNA = "always"))
cat("\nAnhedonia (20441):\n");                 print(table(construct_data$anhedonia, useNA = "always"))
cat("\nCore dysphoria (20446):\n");            print(table(construct_data$dysphoria, useNA = "always"))
cat("\nNegative self-evaluation (20450):\n");  print(table(construct_data$neg_self_eval, useNA = "always"))

# 5. Define phenotype registry

phenotype_registry = list(
  GPNoDep            = list(data = z_score_completed_GPNoDep,            col = "GPNoDep"),
  GPpsy              = list(data = z_score_completed_GPpsy,              col = "GPpsy"),
  Psypsy             = list(data = z_score_completed_Psypsy,             col = "Psypsy"),
  SelfRepDep         = list(data = z_score_completed_SelfRepDep,         col = "SelfRepDep"),
  DepAll             = list(data = z_score_completed_DepAll,             col = "DepAll"),
  ICD10Dep           = list(data = z_score_completed_ICD10Dep,           col = "ICD10Dep"),
  ICD10Dep_exclpsych = list(data = z_score_completed_ICD10Dep_exclpsych, col = "ICD10Dep.exclpsych"),
  LifetimeMDD        = list(data = z_score_completed_LifetimeMDD,        col = "LifetimeMDD"),
  MDDRecurr          = list(data = z_score_completed_MDDRecurr,          col = "MDDRecur")
)

# Construct names for primary analysis (historical/lifetime items)
construct_names_primary = c("cognitive_slowing", 
                            "anhedonia", 
                            "neg_self_eval")

# Construct names for sensitivity analysis (current/recent items)
# No current dysphoria item available
construct_names_sensitivity = c("cognitive_slowing_recent", "anhedonia_recent", "neg_self_eval_recent")

# 6. Build within-cases datasets
covariates = c("bl_31",    # Sex
               "bl_54",    # Age
               "bl_21022", # Assessment Center
               "bl_26521", # Estimated Total Intracranial Volume 
               "bl_24419", # Head Motion 
               paste0("PC", 1:20))

build_within_cases_data = function(z_score_dt, phenotype_col, construct_data, covariates, construct_names) {
  
  # Step 1: Subset to cases
  cases_dt = z_score_dt[z_score_dt[[phenotype_col]] == 1, ]
  cat(sprintf("  Cases for %s: %d\n", phenotype_col, nrow(cases_dt)))
  
  # Step 2: Merge construct columns
  cases_dt = merge(cases_dt, construct_data, by = "subject_num", all.x = TRUE)
  
  # Step 3: For each construct, build analysis-ready datasets
  construct_datasets = list()
  
  for (construct in construct_names) {
    
    bl_cols   = grep("^bl_", names(cases_dt), value = TRUE)
    feat_cols = setdiff(bl_cols, covariates)
    keep_cols = c(feat_cols, covariates, construct)
    
    analysis_dt = cases_dt[, ..keep_cols]
    analysis_dt = analysis_dt[complete.cases(analysis_dt), ]
    
    n_yes = sum(analysis_dt[[construct]] == 1)
    n_no  = sum(analysis_dt[[construct]] == 0)
    
    cat(sprintf("    Construct '%s': n=%d (Yes=%d, No=%d)\n", 
                construct, nrow(analysis_dt), n_yes, n_no))
    
    if (n_yes < 10 || n_no < 10) {
      cat(sprintf("    WARNING: Insufficient cases for '%s' — skipping.\n", construct))
      construct_datasets[[construct]] = NULL
      next
    }
    
    features_df   = analysis_dt[, ..feat_cols]
    covariates_df = analysis_dt[, ..covariates]
    outcome_df    = analysis_dt[, ..construct]
    
    construct_datasets[[construct]] = list(
      features   = features_df,
      covariates = covariates_df,
      outcome    = outcome_df,
      n_total    = nrow(analysis_dt),
      n_yes      = n_yes,
      n_no       = n_no
    )
  }
  
  return(construct_datasets)
}

# Build all within-cases datasets for PRIMARY analysis
cat("\n=== Building within-cases datasets (PRIMARY — historical/lifetime) ===\n\n")

within_cases_data = list()

for (pheno_name in names(phenotype_registry)) {
  cat(sprintf("Processing phenotype: %s\n", pheno_name))
  reg      = phenotype_registry[[pheno_name]]
  datasets = build_within_cases_data(reg$data, reg$col, construct_data, covariates, construct_names_primary)
  within_cases_data[[pheno_name]] = datasets
  cat("\n")
}

# 7. Run logistic regression for each phenotype x construct
#
# Uses run_analysis_for_a_dataset from the main pipeline.
# CI method: profile_likelihood (consistent with main analysis).
#
# Sequential execution: At n~1,500 per construct, each GLM + profile
# likelihood CI completes in ~1-2 sec. 
# Full run of 256 features takes ~5-10 min per construct — 
# parallel overhead exceeds compute savings at this sample size, 
# and avoids FutureInterruptError cascades from
# workers hanging on near-separation cases during profiling.
# The reason future_lapply was used in the main analysis is that sample sizes were 10,000-30,000, 
# where each GLM is meaningfully expensive. 
plan(sequential)

#' Run Within-Cases Symptom Dimension Analysis for a Single Phenotype
#'
#' Iterates over all three constructs (cognitive_slowing, anhedonia, dysphoria)
#' for a given phenotype, runs logistic regression via the existing pipeline,
#' and saves results incrementally.
#'
#' @param pheno_name Character string identifying the phenotype (e.g., "GPNoDep")
#' @param within_cases_data The list built in Section 6
#' @param construct_names Character vector of construct names to iterate over
#' @param output_dir Directory for saving .RDS files
#' @return A named list of results, one entry per construct
run_phenotype_constructs = function(pheno_name, within_cases_data, construct_names, output_dir) {
  
  cat(sprintf("\n=== %s: Running within-cases logistic regressions ===\n", pheno_name))
  pheno_results = list()
  
  for (construct in construct_names) {
    dataset = within_cases_data[[pheno_name]][[construct]]
    
    if (is.null(dataset)) {
      cat(sprintf("  Skipping %s x %s (insufficient data)\n", pheno_name, construct))
      next
    }
    
    cat(sprintf("  Running %s x %s (n=%d, Yes=%d, No=%d)...\n",
                pheno_name, construct, dataset$n_total, dataset$n_yes, dataset$n_no))
    
    tryCatch({
      result = run_analysis_for_a_dataset(
        features_df          = dataset$features,
        covariates_df        = dataset$covariates,
        outcome_df           = dataset$outcome,
        ci_method            = "profile_likelihood",
        conf_level           = 0.95,
        bootstrap_iterations = NULL
      )
      
      # Append metadata for downstream merging and plotting
      for (mtype in names(result)) {
        if (!is.null(result[[mtype]]$results_df) && nrow(result[[mtype]]$results_df) > 0) {
          result[[mtype]]$results_df$phenotype        = pheno_name
          result[[mtype]]$results_df$construct         = construct
          result[[mtype]]$results_df$measurement_type  = mtype
          result[[mtype]]$results_df$n_cases_yes       = dataset$n_yes
          result[[mtype]]$results_df$n_cases_no        = dataset$n_no
        }
      }
      
      pheno_results[[construct]] = result
      
      saveRDS(result, file = file.path(output_dir,
                                       sprintf("LR_within_%s_%s_PL_95.RDS", pheno_name, construct)))
      
      cat(sprintf("    Done. Saved.\n"))
      
    }, error = function(e) {
      cat(sprintf("    ERROR in %s - %s: %s\n", pheno_name, construct, e$message))
    })
  }
  
  return(pheno_results)
}

# Run each phenotype separately

results_GPNoDep            = run_phenotype_constructs("GPNoDep",            within_cases_data, construct_names_primary, output_dir)
results_GPpsy              = run_phenotype_constructs("GPpsy",              within_cases_data, construct_names_primary, output_dir)
results_Psypsy             = run_phenotype_constructs("Psypsy",             within_cases_data, construct_names_primary, output_dir)
results_SelfRepDep         = run_phenotype_constructs("SelfRepDep",         within_cases_data, construct_names_primary, output_dir)
results_DepAll             = run_phenotype_constructs("DepAll",             within_cases_data, construct_names_primary, output_dir)
#results_ICD10Dep           = run_phenotype_constructs("ICD10Dep",           within_cases_data, construct_names_primary, output_dir)
#results_ICD10Dep_exclpsych = run_phenotype_constructs("ICD10Dep_exclpsych", within_cases_data, construct_names_primary, output_dir)
results_LifetimeMDD        = run_phenotype_constructs("LifetimeMDD",        within_cases_data, construct_names_primary, output_dir)
results_MDDRecurr          = run_phenotype_constructs("MDDRecurr",          within_cases_data, construct_names_primary, output_dir)

# Collect into a single list for downstream compilation (Sections 8-10)
all_results = list(
  GPNoDep             = results_GPNoDep,
  GPpsy               = results_GPpsy,
  Psypsy              = results_Psypsy,
  SelfRepDep          = results_SelfRepDep,
  DepAll              = results_DepAll,
  #ICD10Dep           = results_ICD10Dep,
  #ICD10Dep_exclpsych = results_ICD10Dep_exclpsych,
  LifetimeMDD         = results_LifetimeMDD,
  MDDRecurr           = results_MDDRecurr
)

# 8. Region name lookup and results compilation

cat("\n=== Compiling summary tables ===\n")

# White matter tract region names keyed by UKB feature ID.
# FA: Fields 25056–25103; MD: Fields 25104–25151.
region_name_lookup = c(
  # --- Fractional Anisotropy (FA) ---
  "bl_25056" = "Mean FA in middle cerebellar peduncle on FA skeleton",
  "bl_25057" = "Mean FA in pontine crossing tract on FA skeleton",
  "bl_25058" = "Mean FA in genu of corpus callosum on FA skeleton",
  "bl_25059" = "Mean FA in body of corpus callosum on FA skeleton",
  "bl_25060" = "Mean FA in splenium of corpus callosum on FA skeleton",
  "bl_25061" = "Mean FA in fornix on FA skeleton",
  "bl_25062" = "Mean FA in corticospinal tract on FA skeleton (right)",
  "bl_25063" = "Mean FA in corticospinal tract on FA skeleton (left)",
  "bl_25064" = "Mean FA in medial lemniscus on FA skeleton (right)",
  "bl_25065" = "Mean FA in medial lemniscus on FA skeleton (left)",
  "bl_25066" = "Mean FA in inferior cerebellar peduncle on FA skeleton (right)",
  "bl_25067" = "Mean FA in inferior cerebellar peduncle on FA skeleton (left)",
  "bl_25068" = "Mean FA in superior cerebellar peduncle on FA skeleton (right)",
  "bl_25069" = "Mean FA in superior cerebellar peduncle on FA skeleton (left)",
  "bl_25070" = "Mean FA in cerebral peduncle on FA skeleton (right)",
  "bl_25071" = "Mean FA in cerebral peduncle on FA skeleton (left)",
  "bl_25072" = "Mean FA in anterior limb of internal capsule on FA skeleton (right)",
  "bl_25073" = "Mean FA in anterior limb of internal capsule on FA skeleton (left)",
  "bl_25074" = "Mean FA in posterior limb of internal capsule on FA skeleton (right)",
  "bl_25075" = "Mean FA in posterior limb of internal capsule on FA skeleton (left)",
  "bl_25076" = "Mean FA in retrolenticular part of internal capsule on FA skeleton (right)",
  "bl_25077" = "Mean FA in retrolenticular part of internal capsule on FA skeleton (left)",
  "bl_25078" = "Mean FA in anterior corona radiata on FA skeleton (right)",
  "bl_25079" = "Mean FA in anterior corona radiata on FA skeleton (left)",
  "bl_25080" = "Mean FA in superior corona radiata on FA skeleton (right)",
  "bl_25081" = "Mean FA in superior corona radiata on FA skeleton (left)",
  "bl_25082" = "Mean FA in posterior corona radiata on FA skeleton (right)",
  "bl_25083" = "Mean FA in posterior corona radiata on FA skeleton (left)",
  "bl_25084" = "Mean FA in posterior thalamic radiation on FA skeleton (right)",
  "bl_25085" = "Mean FA in posterior thalamic radiation on FA skeleton (left)",
  "bl_25086" = "Mean FA in sagittal stratum on FA skeleton (right)",
  "bl_25087" = "Mean FA in sagittal stratum on FA skeleton (left)",
  "bl_25088" = "Mean FA in external capsule on FA skeleton (right)",
  "bl_25089" = "Mean FA in external capsule on FA skeleton (left)",
  "bl_25090" = "Mean FA in cingulum cingulate gyrus on FA skeleton (right)",
  "bl_25091" = "Mean FA in cingulum cingulate gyrus on FA skeleton (left)",
  "bl_25092" = "Mean FA in cingulum hippocampus on FA skeleton (right)",
  "bl_25093" = "Mean FA in cingulum hippocampus on FA skeleton (left)",
  "bl_25094" = "Mean FA in fornix cres+stria terminalis on FA skeleton (right)",
  "bl_25095" = "Mean FA in fornix cres+stria terminalis on FA skeleton (left)",
  "bl_25096" = "Mean FA in superior longitudinal fasciculus on FA skeleton (right)",
  "bl_25097" = "Mean FA in superior longitudinal fasciculus on FA skeleton (left)",
  "bl_25098" = "Mean FA in superior fronto-occipital fasciculus on FA skeleton (right)",
  "bl_25099" = "Mean FA in superior fronto-occipital fasciculus on FA skeleton (left)",
  "bl_25100" = "Mean FA in uncinate fasciculus on FA skeleton (right)",
  "bl_25101" = "Mean FA in uncinate fasciculus on FA skeleton (left)",
  "bl_25102" = "Mean FA in tapetum on FA skeleton (right)",
  "bl_25103" = "Mean FA in tapetum on FA skeleton (left)",
  
  # --- Mean Diffusivity (MD) ---
  "bl_25104" = "Mean MD in middle cerebellar peduncle on FA skeleton",
  "bl_25105" = "Mean MD in pontine crossing tract on FA skeleton",
  "bl_25106" = "Mean MD in genu of corpus callosum on FA skeleton",
  "bl_25107" = "Mean MD in body of corpus callosum on FA skeleton",
  "bl_25108" = "Mean MD in splenium of corpus callosum on FA skeleton",
  "bl_25109" = "Mean MD in fornix on FA skeleton",
  "bl_25110" = "Mean MD in corticospinal tract on FA skeleton (right)",
  "bl_25111" = "Mean MD in corticospinal tract on FA skeleton (left)",
  "bl_25112" = "Mean MD in medial lemniscus on FA skeleton (right)",
  "bl_25113" = "Mean MD in medial lemniscus on FA skeleton (left)",
  "bl_25114" = "Mean MD in inferior cerebellar peduncle on FA skeleton (right)",
  "bl_25115" = "Mean MD in inferior cerebellar peduncle on FA skeleton (left)",
  "bl_25116" = "Mean MD in superior cerebellar peduncle on FA skeleton (right)",
  "bl_25117" = "Mean MD in superior cerebellar peduncle on FA skeleton (left)",
  "bl_25118" = "Mean MD in cerebral peduncle on FA skeleton (right)",
  "bl_25119" = "Mean MD in cerebral peduncle on FA skeleton (left)",
  "bl_25120" = "Mean MD in anterior limb of internal capsule on FA skeleton (right)",
  "bl_25121" = "Mean MD in anterior limb of internal capsule on FA skeleton (left)",
  "bl_25122" = "Mean MD in posterior limb of internal capsule on FA skeleton (right)",
  "bl_25123" = "Mean MD in posterior limb of internal capsule on FA skeleton (left)",
  "bl_25124" = "Mean MD in retrolenticular part of internal capsule on FA skeleton (right)",
  "bl_25125" = "Mean MD in retrolenticular part of internal capsule on FA skeleton (left)",
  "bl_25126" = "Mean MD in anterior corona radiata on FA skeleton (right)",
  "bl_25127" = "Mean MD in anterior corona radiata on FA skeleton (left)",
  "bl_25128" = "Mean MD in superior corona radiata on FA skeleton (right)",
  "bl_25129" = "Mean MD in superior corona radiata on FA skeleton (left)",
  "bl_25130" = "Mean MD in posterior corona radiata on FA skeleton (right)",
  "bl_25131" = "Mean MD in posterior corona radiata on FA skeleton (left)",
  "bl_25132" = "Mean MD in posterior thalamic radiation on FA skeleton (right)",
  "bl_25133" = "Mean MD in posterior thalamic radiation on FA skeleton (left)",
  "bl_25134" = "Mean MD in sagittal stratum on FA skeleton (right)",
  "bl_25135" = "Mean MD in sagittal stratum on FA skeleton (left)",
  "bl_25136" = "Mean MD in external capsule on FA skeleton (right)",
  "bl_25137" = "Mean MD in external capsule on FA skeleton (left)",
  "bl_25138" = "Mean MD in cingulum cingulate gyrus on FA skeleton (right)",
  "bl_25139" = "Mean MD in cingulum cingulate gyrus on FA skeleton (left)",
  "bl_25140" = "Mean MD in cingulum hippocampus on FA skeleton (right)",
  "bl_25141" = "Mean MD in cingulum hippocampus on FA skeleton (left)",
  "bl_25142" = "Mean MD in fornix cres+stria terminalis on FA skeleton (right)",
  "bl_25143" = "Mean MD in fornix cres+stria terminalis on FA skeleton (left)",
  "bl_25144" = "Mean MD in superior longitudinal fasciculus on FA skeleton (right)",
  "bl_25145" = "Mean MD in superior longitudinal fasciculus on FA skeleton (left)",
  "bl_25146" = "Mean MD in superior fronto-occipital fasciculus on FA skeleton (right)",
  "bl_25147" = "Mean MD in superior fronto-occipital fasciculus on FA skeleton (left)",
  "bl_25148" = "Mean MD in uncinate fasciculus on FA skeleton (right)",
  "bl_25149" = "Mean MD in uncinate fasciculus on FA skeleton (left)",
  "bl_25150" = "Mean MD in tapetum on FA skeleton (right)",
  "bl_25151" = "Mean MD in tapetum on FA skeleton (left)",
  
  # --- Surface Area (SA) ---
  "bl_26721" = "Area of TotalSurface (left hemisphere)",
  "bl_26722" = "Area of bankssts (left hemisphere)",
  "bl_26723" = "Area of caudalanteriorcingulate (left hemisphere)",
  "bl_26724" = "Area of caudalmiddlefrontal (left hemisphere)",
  "bl_26725" = "Area of cuneus (left hemisphere)",
  "bl_26726" = "Area of entorhinal (left hemisphere)",
  "bl_26727" = "Area of fusiform (left hemisphere)",
  "bl_26728" = "Area of inferiorparietal (left hemisphere)",
  "bl_26729" = "Area of inferiortemporal (left hemisphere)",
  "bl_26730" = "Area of isthmuscingulate (left hemisphere)",
  "bl_26731" = "Area of lateraloccipital (left hemisphere)",
  "bl_26732" = "Area of lateralorbitofrontal (left hemisphere)",
  "bl_26733" = "Area of lingual (left hemisphere)",
  "bl_26734" = "Area of medialorbitofrontal (left hemisphere)",
  "bl_26735" = "Area of middletemporal (left hemisphere)",
  "bl_26736" = "Area of parahippocampal (left hemisphere)",
  "bl_26737" = "Area of paracentral (left hemisphere)",
  "bl_26738" = "Area of parsopercularis (left hemisphere)",
  "bl_26739" = "Area of parsorbitalis (left hemisphere)",
  "bl_26740" = "Area of parstriangularis (left hemisphere)",
  "bl_26741" = "Area of pericalcarine (left hemisphere)",
  "bl_26742" = "Area of postcentral (left hemisphere)",
  "bl_26743" = "Area of posteriorcingulate (left hemisphere)",
  "bl_26744" = "Area of precentral (left hemisphere)",
  "bl_26745" = "Area of precuneus (left hemisphere)",
  "bl_26746" = "Area of rostralanteriorcingulate (left hemisphere)",
  "bl_26747" = "Area of rostralmiddlefrontal (left hemisphere)",
  "bl_26748" = "Area of superiorfrontal (left hemisphere)",
  "bl_26749" = "Area of superiorparietal (left hemisphere)",
  "bl_26750" = "Area of superiortemporal (left hemisphere)",
  "bl_26751" = "Area of supramarginal (left hemisphere)",
  "bl_26752" = "Area of frontalpole (left hemisphere)",
  "bl_26753" = "Area of transversetemporal (left hemisphere)",
  "bl_26754" = "Area of insula (left hemisphere)",
  "bl_26822" = "Area of TotalSurface (right hemisphere)",
  "bl_26823" = "Area of bankssts (right hemisphere)",
  "bl_26824" = "Area of caudalanteriorcingulate (right hemisphere)",
  "bl_26825" = "Area of caudalmiddlefrontal (right hemisphere)",
  "bl_26826" = "Area of cuneus (right hemisphere)",
  "bl_26827" = "Area of entorhinal (right hemisphere)",
  "bl_26828" = "Area of fusiform (right hemisphere)",
  "bl_26829" = "Area of inferiorparietal (right hemisphere)",
  "bl_26830" = "Area of inferiortemporal (right hemisphere)",
  "bl_26831" = "Area of isthmuscingulate (right hemisphere)",
  "bl_26832" = "Area of lateraloccipital (right hemisphere)",
  "bl_26833" = "Area of lateralorbitofrontal (right hemisphere)",
  "bl_26834" = "Area of lingual (right hemisphere)",
  "bl_26835" = "Area of medialorbitofrontal (right hemisphere)",
  "bl_26836" = "Area of middletemporal (right hemisphere)",
  "bl_26837" = "Area of parahippocampal (right hemisphere)",
  "bl_26838" = "Area of paracentral (right hemisphere)",
  "bl_26839" = "Area of parsopercularis (right hemisphere)",
  "bl_26840" = "Area of parsorbitalis (right hemisphere)",
  "bl_26841" = "Area of parstriangularis (right hemisphere)",
  "bl_26842" = "Area of pericalcarine (right hemisphere)",
  "bl_26843" = "Area of postcentral (right hemisphere)",
  "bl_26844" = "Area of posteriorcingulate (right hemisphere)",
  "bl_26845" = "Area of precentral (right hemisphere)",
  "bl_26846" = "Area of precuneus (right hemisphere)",
  "bl_26847" = "Area of rostralanteriorcingulate (right hemisphere)",
  "bl_26848" = "Area of rostralmiddlefrontal (right hemisphere)",
  "bl_26849" = "Area of superiorfrontal (right hemisphere)",
  "bl_26850" = "Area of superiorparietal (right hemisphere)",
  "bl_26851" = "Area of superiortemporal (right hemisphere)",
  "bl_26852" = "Area of supramarginal (right hemisphere)",
  "bl_26853" = "Area of frontalpole (right hemisphere)",
  "bl_26854" = "Area of transversetemporal (right hemisphere)",
  "bl_26855" = "Area of insula (right hemisphere)",
  
  # --- Cortical Thickness (CT) ---
  "bl_26755" = "Global Mean thickness (left hemisphere)",
  "bl_26756" = "Mean thickness of bankssts (left hemisphere)",
  "bl_26757" = "Mean thickness of caudalanteriorcingulate (left hemisphere)",
  "bl_26758" = "Mean thickness of caudalmiddlefrontal (left hemisphere)",
  "bl_26759" = "Mean thickness of cuneus (left hemisphere)",
  "bl_26760" = "Mean thickness of entorhinal (left hemisphere)",
  "bl_26761" = "Mean thickness of fusiform (left hemisphere)",
  "bl_26762" = "Mean thickness of inferiorparietal (left hemisphere)",
  "bl_26763" = "Mean thickness of inferiortemporal (left hemisphere)",
  "bl_26764" = "Mean thickness of isthmuscingulate (left hemisphere)",
  "bl_26765" = "Mean thickness of lateraloccipital (left hemisphere)",
  "bl_26766" = "Mean thickness of lateralorbitofrontal (left hemisphere)",
  "bl_26767" = "Mean thickness of lingual (left hemisphere)",
  "bl_26768" = "Mean thickness of medialorbitofrontal (left hemisphere)",
  "bl_26769" = "Mean thickness of middletemporal (left hemisphere)",
  "bl_26770" = "Mean thickness of parahippocampal (left hemisphere)",
  "bl_26771" = "Mean thickness of paracentral (left hemisphere)",
  "bl_26772" = "Mean thickness of parsopercularis (left hemisphere)",
  "bl_26773" = "Mean thickness of parsorbitalis (left hemisphere)",
  "bl_26774" = "Mean thickness of parstriangularis (left hemisphere)",
  "bl_26775" = "Mean thickness of pericalcarine (left hemisphere)",
  "bl_26776" = "Mean thickness of postcentral (left hemisphere)",
  "bl_26777" = "Mean thickness of posteriorcingulate (left hemisphere)",
  "bl_26778" = "Mean thickness of precentral (left hemisphere)",
  "bl_26779" = "Mean thickness of precuneus (left hemisphere)",
  "bl_26780" = "Mean thickness of rostralanteriorcingulate (left hemisphere)",
  "bl_26781" = "Mean thickness of rostralmiddlefrontal (left hemisphere)",
  "bl_26782" = "Mean thickness of superiorfrontal (left hemisphere)",
  "bl_26783" = "Mean thickness of superiorparietal (left hemisphere)",
  "bl_26784" = "Mean thickness of superiortemporal (left hemisphere)",
  "bl_26785" = "Mean thickness of supramarginal (left hemisphere)",
  "bl_26786" = "Mean thickness of frontalpole (left hemisphere)",
  "bl_26787" = "Mean thickness of transversetemporal (left hemisphere)",
  "bl_26788" = "Mean thickness of insula (left hemisphere)",
  "bl_26856" = "Global Mean thickness (right hemisphere)",
  "bl_26857" = "Mean thickness of bankssts (right hemisphere)",
  "bl_26858" = "Mean thickness of caudalanteriorcingulate (right hemisphere)",
  "bl_26859" = "Mean thickness of caudalmiddlefrontal (right hemisphere)",
  "bl_26860" = "Mean thickness of cuneus (right hemisphere)",
  "bl_26861" = "Mean thickness of entorhinal (right hemisphere)",
  "bl_26862" = "Mean thickness of fusiform (right hemisphere)",
  "bl_26863" = "Mean thickness of inferiorparietal (right hemisphere)",
  "bl_26864" = "Mean thickness of inferiortemporal (right hemisphere)",
  "bl_26865" = "Mean thickness of isthmuscingulate (right hemisphere)",
  "bl_26866" = "Mean thickness of lateraloccipital (right hemisphere)",
  "bl_26867" = "Mean thickness of lateralorbitofrontal (right hemisphere)",
  "bl_26868" = "Mean thickness of lingual (right hemisphere)",
  "bl_26869" = "Mean thickness of medialorbitofrontal (right hemisphere)",
  "bl_26870" = "Mean thickness of middletemporal (right hemisphere)",
  "bl_26871" = "Mean thickness of parahippocampal (right hemisphere)",
  "bl_26872" = "Mean thickness of paracentral (right hemisphere)",
  "bl_26873" = "Mean thickness of parsopercularis (right hemisphere)",
  "bl_26874" = "Mean thickness of parsorbitalis (right hemisphere)",
  "bl_26875" = "Mean thickness of parstriangularis (right hemisphere)",
  "bl_26876" = "Mean thickness of pericalcarine (right hemisphere)",
  "bl_26877" = "Mean thickness of postcentral (right hemisphere)",
  "bl_26878" = "Mean thickness of posteriorcingulate (right hemisphere)",
  "bl_26879" = "Mean thickness of precentral (right hemisphere)",
  "bl_26880" = "Mean thickness of precuneus (right hemisphere)",
  "bl_26881" = "Mean thickness of rostralanteriorcingulate (right hemisphere)",
  "bl_26882" = "Mean thickness of rostralmiddlefrontal (right hemisphere)",
  "bl_26883" = "Mean thickness of superiorfrontal (right hemisphere)",
  "bl_26884" = "Mean thickness of superiorparietal (right hemisphere)",
  "bl_26885" = "Mean thickness of superiortemporal (right hemisphere)",
  "bl_26886" = "Mean thickness of supramarginal (right hemisphere)",
  "bl_26887" = "Mean thickness of frontalpole (right hemisphere)",
  "bl_26888" = "Mean thickness of transversetemporal (right hemisphere)",
  "bl_26889" = "Mean thickness of insula (right hemisphere)",
  
  # --- Global Volume ---
  "bl_25001" = "Volume of peripheral cortical grey matter (normalised for head size)",
  "bl_25002" = "Volume of peripheral cortical grey matter",
  "bl_25003" = "Volume of ventricular cerebrospinal fluid (normalised for head size)",
  "bl_25004" = "Volume of ventricular cerebrospinal fluid",
  "bl_25005" = "Volume of grey matter (normalised for head size)",
  "bl_25006" = "Volume of grey matter",
  "bl_25007" = "Volume of white matter (normalised for head size)",
  "bl_25008" = "Volume of white matter",
  "bl_25009" = "Volume of brain, grey+white matter (normalised for head size)",
  "bl_25010" = "Volume of brain, grey+white matter",
  
  # --- Subcortical Volume ---
  "bl_26558" = "Volume of Thalamus-Proper (left hemisphere)",
  "bl_26559" = "Volume of Caudate (left hemisphere)",
  "bl_26560" = "Volume of Putamen (left hemisphere)",
  "bl_26561" = "Volume of Pallidum (left hemisphere)",
  "bl_26562" = "Volume of Hippocampus (left hemisphere)",
  "bl_26563" = "Volume of Amygdala (left hemisphere)",
  "bl_26564" = "Volume of Accumbens-area (left hemisphere)",
  "bl_26589" = "Volume of Thalamus-Proper (right hemisphere)",
  "bl_26590" = "Volume of Caudate (right hemisphere)",
  "bl_26591" = "Volume of Putamen (right hemisphere)",
  "bl_26592" = "Volume of Pallidum (right hemisphere)",
  "bl_26593" = "Volume of Hippocampus (right hemisphere)",
  "bl_26594" = "Volume of Amygdala (right hemisphere)",
  "bl_26595" = "Volume of Accumbens-area (right hemisphere)"
)

# 8b. IDP category classification by feature ID
# Since measurement_type in full results is "All_Features", we classify
# IDP category from the feature ID using the vectors defined in Section 2b.

# Canonical feature ordering within each IDP category
# (matches ordering in Supplementary Tables X, S2c-S2f)
fa_feature_order = paste0("bl_", c("25056", "25057", "25058", "25059", "25060", "25061",
                                   "25062", "25063", "25064", "25065", "25066", "25067", "25068", "25069",
                                   "25070", "25071", "25072", "25073", "25074", "25075", "25076", "25077",
                                   "25078", "25079", "25080", "25081", "25082", "25083", "25084", "25085",
                                   "25086", "25087", "25088", "25089", "25090", "25091", "25092", "25093",
                                   "25094", "25095", "25096", "25097", "25098", "25099", "25100", "25101",
                                   "25102", "25103"))

md_feature_order = paste0("bl_", c("25104", "25105", "25106", "25107", "25108", "25109",
                                   "25110", "25111", "25112", "25113", "25114", "25115", "25116", "25117",
                                   "25118", "25119", "25120", "25121", "25122", "25123", "25124", "25125",
                                   "25126", "25127", "25128", "25129", "25130", "25131", "25132", "25133",
                                   "25134", "25135", "25136", "25137", "25138", "25139", "25140", "25141",
                                   "25142", "25143", "25144", "25145", "25146", "25147", "25148", "25149",
                                   "25150", "25151"))

sa_feature_order = paste0("bl_", c("26721", "26722", "26723", "26724", "26725", "26726",
                                   "26727", "26728", "26729", "26730", "26731", "26732", "26733", "26734",
                                   "26735", "26736", "26737", "26738", "26739", "26740", "26741", "26742",
                                   "26743", "26744", "26745", "26746", "26747", "26748", "26749", "26750",
                                   "26751", "26752", "26753", "26754",
                                   "26822", "26823", "26824", "26825", "26826", "26827", "26828", "26829",
                                   "26830", "26831", "26832", "26833", "26834", "26835", "26836", "26837",
                                   "26838", "26839", "26840", "26841", "26842", "26843", "26844", "26845",
                                   "26846", "26847", "26848", "26849", "26850", "26851", "26852", "26853",
                                   "26854", "26855"))

ct_feature_order = paste0("bl_", c("26755", "26756", "26757", "26758", "26759", "26760",
                                   "26761", "26762", "26763", "26764", "26765", "26766", "26767", "26768",
                                   "26769", "26770", "26771", "26772", "26773", "26774", "26775", "26776",
                                   "26777", "26778", "26779", "26780", "26781", "26782", "26783", "26784",
                                   "26785", "26786", "26787", "26788",
                                   "26856", "26857", "26858", "26859", "26860", "26861", "26862", "26863",
                                   "26864", "26865", "26866", "26867", "26868", "26869", "26870", "26871",
                                   "26872", "26873", "26874", "26875", "26876", "26877", "26878", "26879",
                                   "26880", "26881", "26882", "26883", "26884", "26885", "26886", "26887",
                                   "26888", "26889"))

global_vol_feature_order = paste0("bl_", c("25001", "25002", "25003", "25004", "25005",
                                           "25006", "25007", "25008", "25009", "25010"))

subcort_vol_feature_order = paste0("bl_", c("26558", "26559", "26560", "26561", "26562",
                                            "26563", "26564", "26589", "26590", "26591", "26592", "26593", "26594", "26595"))

idp_category_registry = list(
  FA                  = list(label = "Fractional Anisotropy (FA)",     features = fa_feature_order),
  MD                  = list(label = "Mean Diffusivity (MD)",          features = md_feature_order),
  surface_area        = list(label = "Surface Area (SA)",              features = sa_feature_order),
  thickness           = list(label = "Cortical Thickness (CT)",        features = ct_feature_order),
  global_volume       = list(label = "Global Volume",                  features = global_vol_feature_order),
  subcortical_volume  = list(label = "Subcortical Volume",             features = subcort_vol_feature_order)
)

#' Classify a feature ID into its IDP category
#' @param feature Character vector of feature IDs (e.g., "bl_25056")
#' @return Character vector of IDP category names
classify_idp_category = function(feature) {
  category = rep(NA_character_, length(feature))
  for (cat_name in names(idp_category_registry)) {
    idx = feature %in% idp_category_registry[[cat_name]]$features
    category[idx] = cat_name
  }
  return(category)
}

#' Compile results across phenotypes and constructs into a single data.frame,
#' optionally filtering to specific measurement types (e.g., FA/MD only).
#' Adds region_name by matching the feature column against region_name_lookup.
compile_results = function(all_results, measurement_filter = NULL) {
  rows = list()
  
  for (pheno_name in names(all_results)) {
    for (construct in names(all_results[[pheno_name]])) {
      result = all_results[[pheno_name]][[construct]]
      
      if (!is.null(measurement_filter)) {
        # Pull from IDP-specific keys (e.g., "FA", "MD")
        mtypes = intersect(names(result), measurement_filter)
      } else {
        # Pull from "All_Features" only to avoid duplication
        mtypes = "All_Features"
      }
      
      for (mtype in mtypes) {
        df = result[[mtype]]$results_df
        if (!is.null(df) && nrow(df) > 0) {
          rows[[length(rows) + 1]] = df
        }
      }
    }
  }
  
  if (length(rows) == 0) return(data.frame())
  
  compiled = do.call(rbind, rows)
  
  # Map feature IDs to human-readable region names
  compiled$region_name = region_name_lookup[compiled$feature]
  compiled$region_name[is.na(compiled$region_name)] = compiled$feature[is.na(compiled$region_name)]
  
  # Sort by feature ID within each phenotype x construct for consistent ordering
  compiled = compiled[order(compiled$phenotype, compiled$construct, 
                            compiled$measurement_type, compiled$feature), ]
  
  return(compiled)
}

# Compile primary results: all IDPs and WM-only (FA + MD)
full_results_df = compile_results(all_results)
wm_results_df   = compile_results(all_results, measurement_filter = c("FA", "MD"))

# Save full results as .RDS files
if (nrow(full_results_df) > 0) {
  saveRDS(full_results_df, file = file.path(output_dir, "all_results_compiled.RDS"))
  cat(sprintf("  Full results: %d rows saved as RDS.\n", nrow(full_results_df)))
}
if (nrow(wm_results_df) > 0) {
  saveRDS(wm_results_df, file = file.path(output_dir, "wm_results_compiled.RDS"))
  cat(sprintf("  WM-only results: %d rows saved as RDS.\n", nrow(wm_results_df)))
}


# 9. Write formatted .xlsx workbook (one tab per phenotype)
#
# Each tab contains all WM results for that phenotype, organized
# by feature ID with region_name, across all estimable constructs.

#' Write compiled results to a formatted .xlsx workbook with one tab
#' per phenotype. Columns are ordered for readability: region identifiers
#' first, then construct, sample sizes, effect estimates, and significance.
#'
#' @param results_df  Compiled data.frame from compile_results()
#' @param file_path   Output .xlsx path
#' @param table_title Title string written to row 1 of each tab
write_results_workbook = function(results_df, file_path, table_title) {
  
  wb = createWorkbook()
  
  # Styles
  header_style = createStyle(
    textDecoration = "bold", halign = "center", valign = "center",
    wrapText = TRUE, border = "Bottom", fontSize = 11, fontName = "Arial"
  )
  cell_style = createStyle(halign = "center", valign = "center", fontName = "Arial", fontSize = 11)
  region_style = createStyle(halign = "left", valign = "center", fontName = "Arial", fontSize = 11)
  title_style = createStyle(textDecoration = "bold", fontSize = 12, fontName = "Arial")
  sig_style = createStyle(
    halign = "center", valign = "center", fontName = "Arial", fontSize = 11,
    fontColour = "#006400", textDecoration = "bold"
  )
  
  # Column order for the output table
  display_cols = c("feature", "region_name", "measurement_type", "construct",
                   "n_cases_yes", "n_cases_no", "odds_ratio", "lower_CI", "upper_CI",
                   "p_value", "fdr_p_value", "significant_fdr")
  
  # Readable column headers
  display_headers = c("Feature", "Region", "Measure", "Construct",
                      "N (Yes)", "N (No)", "OR", "CI Lower", "CI Upper",
                      "P-value (unadjusted)", "P-value (FDR-corrected)", "FDR Significant")
  
  for (pheno in unique(results_df$phenotype)) {
    pheno_df = results_df[results_df$phenotype == pheno, ]
    
    # Subset to available display columns (some may be named differently)
    available_cols = intersect(display_cols, names(pheno_df))
    pheno_df = pheno_df[, available_cols, drop = FALSE]
    
    # Sort by measurement type (FA before MD) then feature ID
    pheno_df = pheno_df[order(pheno_df$measurement_type, pheno_df$feature), ]
    rownames(pheno_df) = NULL
    
    # Truncate sheet name to 31 chars (Excel limit)
    sheet_name = substr(pheno, 1, 31)
    addWorksheet(wb, sheet_name)
    
    # Row 1: title
    writeData(wb, sheet_name, x = paste0(table_title, ": ", pheno), startRow = 1, startCol = 1)
    addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)
    mergeCells(wb, sheet_name, cols = 1:length(available_cols), rows = 1)
    
    # Row 3: column headers
    headers = display_headers[match(available_cols, display_cols)]
    writeData(wb, sheet_name, x = t(headers), startRow = 3, startCol = 1, colNames = FALSE)
    addStyle(wb, sheet_name, header_style, rows = 3, cols = 1:length(headers), gridExpand = TRUE)
    setRowHeights(wb, sheet_name, rows = 3, heights = 30)
    
    # Rows 4+: data
    for (i in 1:nrow(pheno_df)) {
      row_num = 3 + i
      
      for (j in 1:ncol(pheno_df)) {
        val = pheno_df[i, j]
        
        # Round numeric columns to 4 decimal places for readability
        if (is.numeric(val)) val = round(val, 4)
        
        writeData(wb, sheet_name, x = val, startRow = row_num, startCol = j)
        
        # Left-align the region_name column; center everything else
        col_name = available_cols[j]
        if (col_name == "region_name") {
          addStyle(wb, sheet_name, region_style, rows = row_num, cols = j)
        } else if (col_name == "significant_fdr" && isTRUE(val)) {
          # Bold green for FDR-significant rows
          addStyle(wb, sheet_name, sig_style, rows = row_num, cols = j)
        } else {
          addStyle(wb, sheet_name, cell_style, rows = row_num, cols = j)
        }
      }
    }
    
    # Column widths: wider for region name, narrower for numeric columns
    setColWidths(wb, sheet_name, cols = 1, widths = 14)                          # Feature ID
    setColWidths(wb, sheet_name, cols = which(available_cols == "region_name"), widths = 48)
    setColWidths(wb, sheet_name, cols = which(available_cols == "construct"), widths = 22)
    numeric_cols = which(available_cols %in% c("measurement_type", "n_cases_yes", "n_cases_no",
                                               "odds_ratio", "ci_lower", "ci_upper",
                                               "p_value", "fdr_adjusted_p", "significant_fdr"))
    setColWidths(wb, sheet_name, cols = numeric_cols, widths = 14)
    
    cat(sprintf("  Tab '%s': %d rows written.\n", sheet_name, nrow(pheno_df)))
  }
  
  saveWorkbook(wb, file_path, overwrite = TRUE)
  cat(sprintf("  Workbook saved: %s\n", file_path))
}

#' Write compiled results to one CSV per IDP category, with section headers
#' separating phenotypes. Feature ordering within each phenotype section
#' matches region_name_lookup order.
#'
#' Output files: one CSV each for FA, MD, SA, CT, Global Volume, Subcortical Volume.
#'
#' Phenotype order within each CSV:
#'   GPNoDep, GPpsy, Psypsy, SelfRepDep, DepAll, LifetimeMDD, MDDRecurr
#'
#' Region order within each phenotype section:
#'   Matches the sequential order in region_name_lookup
#'
#' Output columns (matching Table X):
#'   Feature, Region, Measure, Construct, N (Yes), N (No), OR,
#'   CI Lower, CI Upper, P-value (unadjusted), P-value (FDR-corrected),
#'   FDR Significant
#'
#' @param results_df  Compiled data.frame from compile_results()
#' @param output_dir  Directory for CSV files
#' @param file_prefix Prefix for CSV filenames (IDP category name is appended)
write_results_csvs_by_idp = function(results_df, output_dir, file_prefix = "Supplementary_Table") {
  
  if (nrow(results_df) == 0) {
    cat("  No results to write.\n")
    return(invisible(NULL))
  }
  
  # Classify IDP category from feature IDs
  results_df$idp_category = classify_idp_category(results_df$feature)
  
  # Re-map region_name using the extended lookup
  mapped = region_name_lookup[results_df$feature]
  results_df$region_name = ifelse(is.na(mapped), results_df$feature, mapped)
  
  # Feature ordering from region_name_lookup (sequential position)
  feature_rank = setNames(seq_along(names(region_name_lookup)), names(region_name_lookup))
  
  # Phenotype display order
  phenotype_order = c("GPNoDep", "GPpsy", "Psypsy", "SelfRepDep",
                      "DepAll", "LifetimeMDD", "MDDRecurr")
  
  # IDP category order and file-safe names
  idp_tab_order = c("FA", "MD", "surface_area", "thickness",
                    "global_volume", "subcortical_volume")
  idp_file_names = c(
    FA                 = "FA",
    MD                 = "MD",
    surface_area       = "Surface_Area",
    thickness          = "Cortical_Thickness",
    global_volume      = "Global_Volume",
    subcortical_volume = "Subcortical_Volume"
  )
  
  # Column header row (matches Table X)
  header_row = c("Feature", "Region", "Measure", "Construct",
                 "N (Yes)", "N (No)", "OR", "CI Lower", "CI Upper",
                 "P-value (unadjusted)", "P-value (FDR-corrected)", "FDR Significant")
  n_cols = length(header_row)
  
  for (idp_cat in idp_tab_order) {
    cat_df = results_df[results_df$idp_category == idp_cat, ]
    if (nrow(cat_df) == 0) next
    
    # Filter to phenotypes in order
    cat_df = cat_df[cat_df$phenotype %in% phenotype_order, ]
    if (nrow(cat_df) == 0) next
    
    cat_label = idp_category_registry[[idp_cat]]$label
    
    # Build output lines
    output_lines = list()
    
    # Title row
    output_lines[[1]] = c(
      paste0("Within-Cases Associations (", cat_label, ")"),
      rep("", n_cols - 1)
    )
    
    # Blank row
    output_lines[[2]] = rep("", n_cols)
    
    line_idx = 3
    
    for (pheno in phenotype_order) {
      pheno_df = cat_df[cat_df$phenotype == pheno, ]
      if (nrow(pheno_df) == 0) next
      
      # Sort by region_name_lookup order, then construct alphabetically
      pheno_df$feature_rank = feature_rank[pheno_df$feature]
      pheno_df = pheno_df[order(pheno_df$feature_rank, pheno_df$construct), ]
      
      # Section header row: phenotype name
      output_lines[[line_idx]] = c(pheno, rep("", n_cols - 1))
      line_idx = line_idx + 1
      
      # Column headers
      output_lines[[line_idx]] = header_row
      line_idx = line_idx + 1
      
      # Data rows
      for (i in 1:nrow(pheno_df)) {
        row = pheno_df[i, ]
        output_lines[[line_idx]] = c(
          row$feature,
          row$region_name,
          idp_cat,
          row$construct,
          as.character(row$n_cases_yes),
          as.character(row$n_cases_no),
          as.character(round(row$odds_ratio, 4)),
          as.character(round(row$lower_CI, 4)),
          as.character(round(row$upper_CI, 4)),
          as.character(row$p_value),
          as.character(row$fdr_p_value),
          as.character(row$significant_fdr)
        )
        line_idx = line_idx + 1
      }
      
      # Blank row between phenotype sections
      output_lines[[line_idx]] = rep("", n_cols)
      line_idx = line_idx + 1
    }
    
    # Convert to matrix and write CSV
    output_mat = do.call(rbind, output_lines)
    file_path  = file.path(output_dir, sprintf("%s_%s.csv", file_prefix, idp_file_names[idp_cat]))
    write.table(output_mat, file = file_path, sep = ",", row.names = FALSE, col.names = FALSE,
                quote = TRUE, na = "")
    
    cat(sprintf("  %s: %d data rows written to %s\n", cat_label, nrow(cat_df), basename(file_path)))
  }
  
  cat(sprintf("  All CSVs saved to: %s\n", output_dir))
}


# Write primary WM results workbook
if (nrow(wm_results_df) > 0) {
  write_results_workbook(
    wm_results_df,
    file.path(output_dir, "Supplementary_Table_X_WM_Primary.xlsx"),
    "Within-Cases WM Associations: Historical/Lifetime Symptom Dimensions"
  )
}

if (nrow(full_results_df) > 0) {
  write_results_csvs_by_idp(
    full_results_df,
    output_dir,
    file_prefix = "Supplementary_Table_X_Full_Primary_CHeck"
  )
}

# 9b. Print FDR-significant WM summary to console

cat("\n=== Summary: FDR-significant WM associations by construct ===\n\n")

if (nrow(wm_results_df) > 0) {
  sig_wm = wm_results_df[wm_results_df$significant_fdr == TRUE, ]
  
  if (nrow(sig_wm) > 0) {
    summary_table = sig_wm %>%
      group_by(phenotype, construct, measurement_type) %>%
      summarise(
        n_sig_features = n(),
        median_OR      = median(odds_ratio, na.rm = TRUE),
        range_OR       = paste0(round(min(odds_ratio, na.rm = TRUE), 3), "-", 
                                round(max(odds_ratio, na.rm = TRUE), 3)),
        .groups = "drop"
      ) %>%
      arrange(construct, phenotype, measurement_type)
    
    print(as.data.frame(summary_table), row.names = FALSE)
  } else {
    cat("No FDR-significant white matter associations found.\n")
  }
}

# 9c. Print FDR-significant summary across ALL IDPs to console
cat("\n=== Summary: FDR-significant associations across ALL IDPs by construct ===\n\n")

if (nrow(full_results_df) > 0) {
  
  # Derive IDP category from feature ID 
  full_results_df$idp_category = classify_idp_category(full_results_df$feature)
  
  sig_all = full_results_df[full_results_df$significant_fdr == TRUE, ]
  
  if (nrow(sig_all) > 0) {
    
    idp_order = c("FA", "MD", "surface_area", "thickness", "global_volume", "subcortical_volume")
    
    summary_table_all = sig_all %>%
      group_by(phenotype, construct, idp_category) %>%
      summarise(
        n_sig_features = n(),
        median_OR      = median(odds_ratio, na.rm = TRUE),
        range_OR       = paste0(round(min(odds_ratio, na.rm = TRUE), 3), "-",
                                round(max(odds_ratio, na.rm = TRUE), 3)),
        .groups = "drop"
      ) %>%
      arrange(idp_category, construct, phenotype)
    
    # Print per IDP category for readability
    for (cat in idp_order) {
      cat_subset = as.data.frame(summary_table_all[summary_table_all$idp_category == cat, ])
      if (nrow(cat_subset) == 0) next
      cat_label = idp_category_registry[[cat]]$label
      cat(sprintf("--- %s ---\n", cat_label))
      print(cat_subset, row.names = FALSE)
      cat("\n")
    }
    
    # Catch any unclassified features (e.g. bl_24419)
    unclassified = as.data.frame(summary_table_all[is.na(summary_table_all$idp_category), ])
    if (nrow(unclassified) > 0) {
      cat("--- Unclassified features ---\n")
      print(unclassified, row.names = FALSE)
      cat("\n")
    }
    
    # Overall tally by IDP category
    cat("--- Overall tally by IDP category ---\n")
    tally = sig_all %>%
      group_by(idp_category) %>%
      summarise(
        n_sig_total  = n(),
        n_phenotypes = n_distinct(phenotype),
        n_constructs = n_distinct(construct),
        median_OR    = round(median(odds_ratio, na.rm = TRUE), 4),
        .groups = "drop"
      ) %>%
      arrange(idp_category)
    print(as.data.frame(tally), row.names = FALSE)
    
    # --------------------------------------------------------
    # Detailed listing: each FDR-significant association
    # --------------------------------------------------------
    cat("\n--- FDR-significant associations by phenotype and region ---\n\n")
    
    detail_cols = c("phenotype", "construct", "idp_category", "region_name",
                    "n_cases_yes", "n_cases_no",
                    "odds_ratio", "lower_CI", "upper_CI",
                    "p_value", "fdr_p_value")
    
    sig_detail = sig_all[, detail_cols]
    sig_detail$idp_category = factor(sig_detail$idp_category, levels = idp_order)
    sig_detail = sig_detail[order(sig_detail$phenotype,
                                  sig_detail$idp_category,
                                  sig_detail$construct,
                                  sig_detail$region_name), ]
    
    for (pheno in sort(unique(sig_detail$phenotype))) {
      cat(sprintf("=== %s ===\n", pheno))
      pheno_rows = sig_detail[sig_detail$phenotype == pheno, ]
      
      for (cat in levels(pheno_rows$idp_category)) {
        cat_rows = pheno_rows[!is.na(pheno_rows$idp_category) & pheno_rows$idp_category == cat, ]
        if (nrow(cat_rows) == 0) next
        cat_label = idp_category_registry[[as.character(cat)]]$label
        cat(sprintf("  [ %s ]\n", cat_label))
        
        for (i in seq_len(nrow(cat_rows))) {
          r = cat_rows[i, ]
          cat(sprintf(
            "    %-55s | construct: %-20s | OR: %5.3f [%5.3f, %5.3f] | p=%.4f | FDR-p=%.4f\n",
            r$region_name, r$construct,
            r$odds_ratio, r$lower_CI, r$upper_CI,
            r$p_value, r$fdr_p_value
          ))
        }
        cat("\n")
      }
      
      # Unclassified features for this phenotype
      unclass_rows = pheno_rows[is.na(pheno_rows$idp_category), ]
      if (nrow(unclass_rows) > 0) {
        cat("  [ Unclassified ]\n")
        for (i in seq_len(nrow(unclass_rows))) {
          r = unclass_rows[i, ]
          cat(sprintf(
            "    %-55s | construct: %-20s | OR: %5.3f [%5.3f, %5.3f] | p=%.4f | FDR-p=%.4f\n",
            r$region_name, r$construct,
            r$odds_ratio, r$lower_CI, r$upper_CI,
            r$p_value, r$fdr_p_value
          ))
        }
        cat("\n")
      }
    }
    
  } else {
    cat("No FDR-significant associations found across any IDP category.\n")
  }
}

# 10. Sample size summary across all analyses

cat("\n=== Sample sizes per phenotype x construct ===\n\n")

for (pheno_name in names(within_cases_data)) {
  for (construct in construct_names_primary) {
    dataset = within_cases_data[[pheno_name]][[construct]]
    if (!is.null(dataset)) {
      cat(sprintf("  %-20s x %-20s: n=%5d (Yes=%5d, No=%5d)\n",
                  pheno_name, construct, dataset$n_total, dataset$n_yes, dataset$n_no))
    } else {
      cat(sprintf("  %-20s x %-20s: SKIPPED (insufficient data)\n",
                  pheno_name, construct))
    }
  }
}

cat("\n=== Analysis complete ===\n")
cat(sprintf("Results saved to: %s\n", output_dir))

# 11. Supplementary Table Z: Sample Characteristics for
#     Within-Cases Symptom Dimension Analyses
#
# 1. Primary (Historical/Lifetime): 4 constructs x 9 phenotypes
#
# Cells show: N (% endorsing) or "NE" for phenotypes excluded
# due to quasi-complete separation.

# Phenotypes excluded due to quasi-complete separation
excluded_primary = list(
  ICD10Dep           = c("cognitive_slowing", "anhedonia", "dysphoria", "neg_self_eval"),
  ICD10Dep_exclpsych = c("cognitive_slowing", "anhedonia", "dysphoria", "neg_self_eval"),
  MDDRecurr          = c("dysphoria")
)

# Readable labels
phenotype_labels = c(
  GPNoDep            = "GPNoDep",
  GPpsy              = "GPpsy",
  Psypsy             = "Psypsy",
  SelfRepDep         = "SelfRepDep",
  DepAll             = "DepAll",
  ICD10Dep           = "ICD10Dep",
  ICD10Dep_exclpsych = "ICD10Dep_exclpsych",
  LifetimeMDD        = "LifetimeMDD",
  MDDRecurr          = "MDDRecurr"
)

construct_labels_primary = c(
  cognitive_slowing = "Cognitive difficulty\n(Field 20435)",
  anhedonia         = "Anhedonia\n(Field 20441)",
  dysphoria         = "Core dysphoria\n(Field 20446)",
  neg_self_eval     = "Negative self-evaluation\n(Field 20450)"
)


build_sample_table = function(within_cases_data, construct_names, construct_labels, 
                              phenotype_labels, excluded = list()) {
  
  pheno_names = names(phenotype_labels)
  
  mat = matrix(NA_character_, nrow = length(pheno_names), ncol = 1 + length(construct_names))
  rownames(mat) = phenotype_labels
  colnames(mat) = c("Total Cases (with MHQ)", construct_labels[construct_names])
  
  for (i in seq_along(pheno_names)) {
    pheno = pheno_names[i]
    
    total_shown = FALSE
    
    for (construct in construct_names) {
      dataset = within_cases_data[[pheno]][[construct]]
      col_idx = which(construct_names == construct) + 1
      
      is_ne = pheno %in% names(excluded) && construct %in% excluded[[pheno]]
      
      if (is.null(dataset)) {
        mat[i, col_idx] = if (is_ne) "NE (insufficient data)" else "Insufficient data"
        next
      }
      
      if (!total_shown) {
        mat[i, 1] = as.character(dataset$n_total)
        total_shown = TRUE
      }
      
      pct = round(100 * dataset$n_yes / dataset$n_total, 1)
      cell_text = sprintf("%d / %d (%s%%)", dataset$n_yes, dataset$n_total, format(pct, nsmall = 1))
      
      # For NE cells, append footnote marker rather than replacing content
      mat[i, col_idx] = if (is_ne) paste0(cell_text, " \u2020") else cell_text
    }
  }
  
  return(as.data.frame(mat, stringsAsFactors = FALSE))
}

# Build both tables
primary_table = build_sample_table(
  within_cases_data, construct_names_primary, construct_labels_primary,
  phenotype_labels, excluded = excluded_primary
)


# Create workbook
wb = createWorkbook()

# --- Shared styles ---
header_style  = createStyle(textDecoration = "bold", halign = "center", valign = "center",
                            wrapText = TRUE, border = "Bottom", fontSize = 12, fontName = "Times New Roman")
row_label_style = createStyle(textDecoration = "bold", valign = "center", fontName = "Times New Roman", fontSize = 12)
cell_style    = createStyle(halign = "center", valign = "center", fontName = "Times New Roman", fontSize = 12)
ne_style      = createStyle(halign = "center", valign = "center", fontName = "Times New Roman", fontSize = 12,
                            fontColour = "#808080", textDecoration = "italic")
note_style    = createStyle(fontName = "Arial", fontSize = 10, textDecoration = "italic", wrapText = TRUE)

write_sheet = function(wb, sheet_name, df, title, footnote) {
  addWorksheet(wb, sheet_name)
  
  # Styles
  title_style     = createStyle(textDecoration = "bold", fontSize = 12, fontName = "Arial")
  header_style    = createStyle(textDecoration = "bold", halign = "center", valign = "center",
                                wrapText = TRUE, border = "Bottom", fontSize = 12, fontName = "Times New Roman")
  row_label_style = createStyle(textDecoration = "bold", valign = "center", 
                                fontName = "Times New Roman", fontSize = 12)
  cell_style      = createStyle(halign = "center", valign = "center", 
                                fontName = "Times New Roman", fontSize = 12)
  ne_style        = createStyle(halign = "center", valign = "center", 
                                fontName = "Times New Roman", fontSize = 12,
                                fontColour = "#8B0000", textDecoration = "italic")
  note_style      = createStyle(fontName = "Arial", fontSize = 10, 
                                textDecoration = "italic", wrapText = TRUE)
  
  # Title
  writeData(wb, sheet_name, x = title, startRow = 1, startCol = 1)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)
  mergeCells(wb, sheet_name, cols = 1:(ncol(df) + 1), rows = 1)
  
  # Column headers
  headers = c("Depression Phenotype", colnames(df))
  writeData(wb, sheet_name, x = t(headers), startRow = 3, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet_name, header_style, rows = 3, cols = 1:length(headers), gridExpand = TRUE)
  setRowHeights(wb, sheet_name, rows = 3, heights = 45)
  
  # Data rows
  for (i in 1:nrow(df)) {
    row_num = 3 + i
    writeData(wb, sheet_name, x = rownames(df)[i], startRow = row_num, startCol = 1)
    addStyle(wb, sheet_name, row_label_style, rows = row_num, cols = 1)
    
    for (j in 1:ncol(df)) {
      val = df[i, j]
      writeData(wb, sheet_name, x = val, startRow = row_num, startCol = j + 1)
      # Apply NE style (dark red italic) to cells with dagger marker
      if (!is.na(val) && grepl("\u2020", val)) {
        addStyle(wb, sheet_name, ne_style, rows = row_num, cols = j + 1)
      } else {
        addStyle(wb, sheet_name, cell_style, rows = row_num, cols = j + 1)
      }
    }
  }
  
  # Column widths
  setColWidths(wb, sheet_name, cols = 1, widths = 30)
  setColWidths(wb, sheet_name, cols = 2:(ncol(df) + 1), widths = 25)
  
  # Footnote
  footnote_row = 3 + nrow(df) + 2
  writeData(wb, sheet_name, x = footnote, startRow = footnote_row, startCol = 1)
  addStyle(wb, sheet_name, note_style, rows = footnote_row, cols = 1)
  mergeCells(wb, sheet_name, cols = 1:(ncol(df) + 1), rows = footnote_row)
}

# Updated footnotes referencing the dagger
primary_footnote = paste0(
  "Values shown as: N endorsing / N total (% endorsing). ",
  "\u2020 Not estimable (NE): model excluded due to quasi-complete separation driven by ",
  "near-ceiling endorsement rates relative to the number of model parameters; ",
  "raw endorsement figures are shown for completeness. ",
  "Historical/lifetime items from the Mental Health Questionnaire (MHQ, Category 138). ",
  "Cognitive difficulty: Field 20435; Anhedonia: Field 20441; Core dysphoria: Field 20446; ",
  "Negative self-evaluation: Field 20450."
)


write_sheet(wb, "Primary (Lifetime)", primary_table,
            "Within-Cases Sample Characteristics: Historical/Lifetime Symptom Dimensions",
            primary_footnote)

# Save
table_path = file.path(output_dir, "Supplementary_Table_S3_Clinical_Characterization_of_UK_Biobank_Depression_Phenotypes.xlsx")
saveWorkbook(wb, table_path, overwrite = TRUE)
cat(sprintf("\nSupplementary Table Z saved to: %s\n", table_path))

# END
