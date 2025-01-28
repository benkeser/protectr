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
.libPaths("~/Rlibs")

here::i_am("run_simulation.R")

library(fastverse)
library(tidyverse)
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
# Opt for parallelly package instead, request appropriate number of nodes in bash script
# ncores <- parallel::detectCores() 
ncores <- parallelly::availableCores()
ncores_for_future <- max(ncores - 1, 1) 
future::plan('multisession', workers = ncores_for_future)

# 0. Get settings from config file ----------------------------------------------------------
setting <- Sys.getenv("SETTING")
config <- config::get(file = "config.yml", config = setting)

# 1. Create or load data --------------------------------------------------------------------

#if data/weekly_records_data_SETTING.csv exists, load; otherwise, make weekly records data
if(file.exists(here::here(paste0("data/", setting, "_weekly_records_data.rds")))){
    weekly_records_data <- setDT(readRDS(
      here::here(paste0("data/", setting, "_weekly_records_data.rds"))
    ))
} else{
  # Read in raw data
  dat <- readRDS(here::here(paste0("data/", config$weekly_records$data_name)))
  
  # Create weekly record data
  weekly_records_data <- create_weekly_record_data(
      dat = dat,
      baseline_covariates = config$weekly_records$baseline_covariates,
      k = config$weekly_records$k,
      admin_cens_wks = config$weekly_records$admin_cens_wks
  )

  saveRDS(weekly_records_data, 
      here::here(paste0("data/weekly_records_data_", setting, ".rds"))
  )
}

# 2. Subset data ----------------------------------------------------------------------------

# Potential adjustments:
# (1) excluding participants with TB diagnosis date within some specified period after enrollment into care
# (2) subsetting data based on year of enrollment into care
# (3) shortening the administrative censoring follow-up period

if(!is.na(config$exclusion_period)){
  
  # QUESTION - tb_diagnosis_date seems to be missing from haiti but it's in uganda and zimbabwe
  # is it supposed to be in all files? could it be named something else
    weekly_records_data_exclude <- weekly_records_data[
        is.na(tb_diagnosis_date) | (tb_diagnosis_date - enroll_date > config$exclusion_period)
    ]
}

if(!is.na(config$exclusion_date)){
  
  # QUESTION same as above with exclusion date haiti
    exclusion_date <- as.Date(config$exclusion_date)
    weekly_records_data_exclude <- weekly_records_data[
	    enroll_date > exclusion_date
    ]
}

if(!is.na(config$new_admin_cens_wk)){
    # should be <= admin_cens_wk used in call to create_weekly_records
    weekly_records_data_exclude <- weekly_records_data[
        wk < config$new_admin_cens_wk
    ]
    # confirm this is correct
    weekly_records_data_exclude$admin_cens_wk <- config$new_admin_cens_wk
}

# also example in Kenya code using exclude_weeks?? 
# not sure if that's doing the same thing with dif var names/if one of the files is more updated than the other

# 3. Fitting propensity scores --------------------------------------------------------------

propensity_output <- fit_propensity_models(
	weekly_records_data = weekly_records_data,
	grace_pd_wks = config$grace_pd_wks,
	denom_model_formula = config$propensity_formulas$num_and_denom_model_formula,
	num_model_formula = config$propensity_formulas$num_and_denom_model_formula,
	right_cens_model_formula = config$propensity_formulas$right_cens_model_formula
)

# i think we probably wanted to save/visualize some of these too?
# TBD once seth shares script

# 4. Cloning procedure -----------------------------------------------------------------------

cloned_data_sets <- create_cloned_data_set(
	propensity_output = propensity_output
)

# 5. Fit marginal structural models (MSMs) ----------------------------------------------------

# iterate through all formulas specified in msm_formulas in configuration file

msm_formula_list <- list(length = length(config$msm_formulas))

for(i in 1:length(config$msm_formulas)){

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
        msm_fit_tb = msm_fit_tb ,
        msm_fit_death = msm_fit_death,
        msm_fit_death_for_tb = msm_fit_death_for_tb
    )

    # Save formulas in list
    msm_formula_list[[i]] <- fit_models

}

# visualize output with cuminc? or do later

# 6. Bootstrap ----------------------------------------------------------------

future::plan('multisession', workers = ncores - 1)

bootstrap_results <- run_bootstrap(
  nboot = config$nboot,
  weekly_records_data = weekly_records_data,
  grace_pd_wks = config$grace_pd_wks,
  denom_model_formula = config$propensity_formulas$num_and_denom_model_formula,
  num_model_formula = config$propensity_formulas$num_and_denom_model_formula,
  right_cens_model_formula = config$propensity_formulas$right_cens_model_formula,  
  admin_cens_wks = config$new_admin_cens_wk * 2,  # Check if this is supposed to be double?? 
  msm_formulas_tb = config$msm_formulas, # Check if these are allowed to vary or should all be same
  msm_formulas_death = config$msm_formulas,
  msm_formulas_death_for_tb = config$msm_formulas
)

# save bootstrap results
saveRDS(
	bootstrap_results,
	here::here(paste0("data/bootstrap_results_", setting, ".rds"))
)

# Additional processing in example analysis?? then save?
