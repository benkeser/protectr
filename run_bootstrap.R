# --------------------------------------------------------------------------------
#    Main R script to run PROTECT analysis for bootstrap estimates on HPC Cluster
# 0. Load configuration
# 1. Load data
# 2. Subset people
# 3. Make boot data
# 4. Fit propensity models
# 5. Cloning process
# 6. MSM models
# --------------------------------------------------------------------------------

# For initial debugging scratch file
options(echo = TRUE)

# increase size allowed future (default 500mb)
options(future.globals.maxSize = 5 * 1024^3) # 5 GB
options(future.globals.onReference = "ignore")

# Path to installed packages on cluster
.libPaths(c("~/Rlibs", "/apps/R/4.4.0/lib64/R/site/library", .libPaths()))

here::i_am("run_bootstrap.R")

library(tidyverse)
library(fastverse)
library(future.apply)
library(progressr)

source(here::here("code/create_weekly_records.R"))
source(here::here("code/fit_propensity_models.R"))
source(here::here("code/compute_cf_init_dist.R"))
source(here::here("code/create_cloned_data_set.R"))
source(here::here("code/fit_msm.R"))
source(here::here("code/compute_cuminc.R"))

# Set seed, read in config file
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
setting <- Sys.getenv("SETTING")

set.seed(seed)
config <- config::get(file = "config.yml", config = setting)

# config args passed into bootstrap functions
grace_pd_wks <- config$grace_pd_wks
denom_model_formula <- config$propensity_formulas$denom_model_formula
num_model_formula <- config$propensity_formulas$num_model_formula
right_cens_model_formula <- config$propensity_formulas$right_cens_model_formula
msm_formulas_tb <- config$msm_formulas
msm_formulas_death <- config$msm_formulas
msm_formulas_death_for_tb <- config$msm_formulas
admin_cens_wks <- config$admin_cens_wk
skip_cens_for_tb <- config$skip_cens_for_tb

# 1. Load weekly records data --------------------------------------------------
if (file.exists(here::here(paste0("data/", setting, "_weekly_records_data.rds")))) {
  weekly_records_data <- setDT(readRDS(
    here::here(paste0("data/", setting, "_weekly_records_data.rds"))
  ))
} else {
  stop(paste0("Weekly records data for setting ", setting, " cannot be found"))
}

# 2. Adjust weekly records data based on config file ---------------------------
if (!is.na(config$exclusion_period)) {
  weekly_records_data <- weekly_records_data[
    is.na(tb_diagnosis_date) | (tb_diagnosis_date - enroll_date > config$exclusion_period)
  ]
}

if (!is.na(config$exclusion_date)) {
  exclusion_date <- as.Date(config$exclusion_date)
  weekly_records_data <- weekly_records_data[
    enroll_date > exclusion_date
  ]
}

if (!is.na(config$admin_cens_wk)) {
  # should be <= admin_cens_wk used in call to create_weekly_records
  weekly_records_data <- weekly_records_data[
    wk < config$admin_cens_wk
  ]
  
  weekly_records_data$admin_cens_wk <- config$admin_cens_wk
}

# 3. Create bootstrap sample ---------------------------------------------------
sampled_ids <- sample(unique(weekly_records_data$id), replace = TRUE)

weekly_records_data_bootstrap_by_id <- lapply(sampled_ids, function(this_id) {
  return(weekly_records_data[id == this_id])
})

weekly_records_data_bootstrap <- rbindlist(
  weekly_records_data_bootstrap_by_id,
  idcol = "new_id"
)
setnames(weekly_records_data_bootstrap, old = "id", new = "orig_id")
setnames(weekly_records_data_bootstrap, old = "new_id", new = "id")

# 4. Fit propensity models -----------------------------------------------------

propensity_output <- fit_propensity_models(
  weekly_records_data_bootstrap,
  grace_pd_wks = grace_pd_wks,
  denom_model_formula = denom_model_formula,
  num_model_formula = num_model_formula,
  right_cens_model_formula = right_cens_model_formula,
  return_models = FALSE,
  skip_cens_for_tb = skip_cens_for_tb
)

# skipped this in the main function bc can be recreated later but left in here 
# to avoid breaking existing boot functions
cf_init_dist <- compute_cf_init_dist(
  propensity_output = propensity_output
)

# 5. Clone data ----------------------------------------------------------------
cloned_data_sets <- create_cloned_data_set(
  propensity_output = propensity_output
)

# 6. Fit MSMs ------------------------------------------------------------------
msm_fits_tb <- sapply(msm_formulas_tb,
                      FUN = fit_msm,
                      cloned_data_set = cloned_data_sets$tb,
                      return_msm_model = TRUE,
                      return_msm_vcov = FALSE,
                      simplify = FALSE,
                      USE.NAMES = TRUE
)

msm_fits_death <- sapply(msm_formulas_death,
                         FUN = fit_msm,
                         cloned_data_set = cloned_data_sets$death,
                         return_msm_model = TRUE,
                         return_msm_vcov = FALSE,
                         simplify = FALSE,
                         USE.NAMES = TRUE
)

msm_fits_death_for_tb <- sapply(msm_formulas_death_for_tb,
                                FUN = fit_msm,
                                cloned_data_set = cloned_data_sets$death_for_tb,
                                return_msm_model = TRUE,
                                return_msm_vcov = FALSE,
                                simplify = FALSE,
                                USE.NAMES = TRUE
)

#  Get cumulative incidence (doing later for point ests)
cuminc <- mapply(
  FUN = compute_cuminc,
  msm_fit_tb = msm_fits_tb,
  msm_fit_death = msm_fits_death,
  msm_fit_death_for_tb = msm_fits_death_for_tb,
  MoreArgs = list(max_wk = 52 * 2),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)

# !!!!!!! replace msm_models with NULL to avoid ballooning object sizes !!!!!!
for (i in seq_len(length(msm_fits_tb))) {
  msm_fits_tb[[i]]$msm_model <- NULL
}
for (i in seq_len(length(msm_fits_death))) {
  msm_fits_death[[i]]$msm_model <- NULL
}
for (i in seq_len(length(msm_fits_death_for_tb))) {
  msm_fits_death_for_tb[[i]]$msm_model <- NULL
}

# Make & save output object for given seed
out <- list(
  cf_init_dist = list(cf_init_dist),
  msm_fits_tb = msm_fits_tb,
  msm_fits_death = msm_fits_death,
  msm_fits_death_for_tb = msm_fits_death_for_tb,
  cuminc = cuminc,
  error = FALSE
)

saveRDS(out, here::here(paste0("/projects/dbenkes/allison/protectr/boot_res/", setting, "_seed_", seed, ".Rds")))
