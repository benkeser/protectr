# --------------------------------------------------------------------------------
#             Main R script to run PROTECT analysis on HPC Cluster
# 0. Set desired settings in config.yml
# 1. Create data
# 2. Subset people
# 3. Fit propensity models
# 4. Cloning process
# 5. MSM models
# 6. Bootstrap
# --------------------------------------------------------------------------------

# For initial debugging scratch file
options(echo = TRUE)

# increase size allowed future (default 500mb)
options(future.globals.maxSize = 5 * 1024^3) # 5 GB
options(future.globals.onReference = NULL)

# Path to installed packages on cluster
.libPaths(c("~/Rlibs", "/apps/R/4.4.0/lib64/R/site/library", .libPaths()))
# .libPaths("/apps/R/4.4.0/lib64/R/site/library")

.libPaths()

here::i_am("run_simulation.R")

library(tidyverse)
library(fastverse)
library(future.apply)
library(progressr)

set.seed(404)

# Source files needed
source(here::here("code/create_weekly_records.R"))
source(here::here("code/fit_propensity_models.R"))
source(here::here("code/compute_cf_init_dist.R"))
source(here::here("code/create_cloned_data_set.R"))
source(here::here("code/fit_msm.R"))
source(here::here("code/compute_cuminc.R"))
source(here::here("code/bootstrap.R"))

# Setup parallelization

# NOTE on the cluster this detects total cores on the node, not necessarily the ones that have been allocated by the scheduler
# Opt for parallelly package instead, request appropriate number of nodes allocated in bash script
# ncores <- parallel::detectCores()
ncores <- parallelly::availableCores()
#ncores_for_future <- max(ncores - 1, 1)
#future::plan("multicore", workers = ncores)
future::plan("multisession", workers = ncores)

# 0. Get settings from config file ----------------------------------------------------------
setting <- Sys.getenv("SETTING")
config <- config::get(file = "config.yml", config = setting)

# Create dir to save results if does not exist
if (!file.exists(here::here(paste0("results/", setting)))) {
  dir.create(here::here(paste0("results/", setting)), recursive = TRUE)
}

# 1. Create or load data --------------------------------------------------------------------

# if data/weekly_records_data_SETTING.csv exists, load; otherwise, make weekly records data
if (file.exists(here::here(paste0("data/", setting, "_weekly_records_data.rds")))) {
  weekly_records_data <- setDT(readRDS(
    here::here(paste0("data/", setting, "_weekly_records_data.rds"))
  ))
} else {
  # Read in raw data
  dat <- readRDS(here::here(paste0("data/", config$weekly_records$data_name)))

  # temp to track time for creating weekly records
  system.time({
    # Create weekly record data
    weekly_records_data <- create_weekly_record_data(
      dat = dat,
      baseline_covariates = config$weekly_records$baseline_covariates,
      k = config$weekly_records$k,
      admin_cens_wks = config$weekly_records$admin_cens_wks
    )

    saveRDS(
      weekly_records_data,
      here::here(paste0("data/", setting, "_weekly_records_data.rds"))
    )
  })
}

# 2. Subset data ----------------------------------------------------------------------------

# Potential adjustments:
# (1) excluding participants with TB diagnosis date within some specified period after enrollment into care
# (2) subsetting data based on year of enrollment into care
# (3) shortening the administrative censoring follow-up period

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

# 3. Fitting propensity scores --------------------------------------------------------------


system.time({
  propensity_output <- fit_propensity_models(
    weekly_records_data = weekly_records_data,
    grace_pd_wks = config$grace_pd_wks,
    denom_model_formula = config$propensity_formulas$denom_model_formula,
    num_model_formula = config$propensity_formulas$num_model_formula,
    right_cens_model_formula = config$propensity_formulas$right_cens_model_formula
  )
})

# save later because need weekly records data with weights

# 4. Cloning procedure -----------------------------------------------------------------------

cloned_data_sets <- create_cloned_data_set(
  propensity_output = propensity_output
)

# Null out parts of propensity_output no longer needed & save
propensity_output$weekly_records_data <- NULL
propensity_output[["models"]][["denom_model"]]$qr <- NULL
propensity_output[["models"]][["num_model"]]$qr <- NULL
propensity_output[["models"]][["cens_model_tb"]]$qr <- NULL
propensity_output[["models"]][["cens_model_death"]]$qr <- NULL

saveRDS(
  propensity_output,
  here::here(paste0("results/", setting, "/propensity_output_", setting, ".rds"))
)

# 5. Fit marginal structural models (MSMs) ----------------------------------------------------

# iterate through all formulas specified in msm_formulas in configuration file
msm_formula_list <- vector("list", length = length(config$msm_formulas))

for (i in 1:length(config$msm_formulas)) {
  msm_formula <- config$msm_formulas[[i]]

  msm_fit_tb <- fit_msm(
    cloned_data_set = cloned_data_sets$tb,
    msm_formula = msm_formula,
    return_msm_model = TRUE,
    gee = FALSE,
    return_msm_vcov = FALSE
  )

  msm_fit_death <- fit_msm(
    cloned_data_set = cloned_data_sets$death,
    msm_formula = msm_formula,
    return_msm_model = TRUE,
    gee = FALSE,
    return_msm_vcov = FALSE
  )

  msm_fit_death_for_tb <- fit_msm(
    cloned_data_set = cloned_data_sets$death_for_tb,
    msm_formula = msm_formula,
    return_msm_model = TRUE,
    gee = FALSE,
    return_msm_vcov = FALSE
  )

  fit_models <- list(
    msm_fit_tb = msm_fit_tb,
    msm_fit_death = msm_fit_death,
    msm_fit_death_for_tb = msm_fit_death_for_tb
  )

  # Save formulas in list
  msm_formula_list[[i]] <- fit_models
}

# Combine MSMs and null out parts no longer needed
for (i in seq_along(msm_formula_list)) {
  msm_formula_list[[i]]$msm_fit_tb$msm_model$qr <- NULL
  msm_formula_list[[i]]$msm_fit_death$msm_model$qr <- NULL
  msm_formula_list[[i]]$msm_fit_death_for_tb$msm_model$qr <- NULL
}

# Save MSMs
saveRDS(
  msm_formula_list,
  here::here(paste0("results/", setting, "/msm_formula_list_", setting, ".rds"))
)

# 6. Bootstrap ----------------------------------------------------------------

system.time({
  bootstrap_results <- run_bootstrap(
    nboot = config$n_boot,
    weekly_records_data = weekly_records_data,
    grace_pd_wks = config$grace_pd_wks,
    denom_model_formula = config$propensity_formulas$denom_model_formula,
    num_model_formula = config$propensity_formulas$num_model_formula,
    right_cens_model_formula = config$propensity_formulas$right_cens_model_formula,
    admin_cens_wks = config$admin_cens_wk,
    msm_formulas_tb = config$msm_formulas,
    msm_formulas_death = config$msm_formulas,
    msm_formulas_death_for_tb = config$msm_formulas
  )
})

saveRDS(
  bootstrap_results,
  here::here(paste0("results/", setting, "/bootstrap_results_", setting, ".rds"))
)

# Get bootstrap CI
bootstrap_ci <- get_bootstrap_ci(bootstrap_results)

saveRDS(
  bootstrap_ci,
  here::here(paste0("results/", setting, "/bootstrap_ci_", setting, ".rds"))
)
