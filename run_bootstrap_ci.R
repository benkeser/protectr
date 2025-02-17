# ------------------------------------------------------------------------
# Script to aggregate bootstrap results and produce confidence interval
# ------------------------------------------------------------------------

here::i_am("run_bootstrap_ci.R")

# For initial debugging scratch file
options(echo = TRUE)

# Path to installed packages on cluster
.libPaths(c("~/Rlibs", "/apps/R/4.4.0/lib64/R/site/library", .libPaths()))

source(here::here("code/bootstrap.R"))

cargs <- commandArgs(TRUE)
setting <- cargs[[1]]

config <- config::get(file = "config.yml", config = setting)

# Path to individual bootstrap results
dir <- "/projects/dbenkes/allison/protectr/boot_res/"

# pull all results for given setting
pattern <- paste0("^", setting, "_seed_[0-9]+\\.Rds$")
all_files <- list.files(dir, pattern = pattern, full.names = TRUE)

# get idxs of correct order 1-1000 vs alphabetical and sort
numeric_order <- as.numeric(gsub(".*_([0-9]+)\\.Rds$", "\\1", all_files))
all_files_sorted <- all_files[order(numeric_order)]

boot_results <- lapply(all_files_sorted, readRDS)

# transform to matrix format compatible with existing functions
boot_results_transform <- t(do.call(rbind, boot_results))

# save all in one boot res in results dir
saveRDS(
  boot_results_transform,
  here::here(paste0("results/", setting, "/bootstrap_results_", setting, ".rds"))
)

# Get bootstrap CI
bootstrap_ci <- get_bootstrap_ci(boot_results_transform)

saveRDS(
  bootstrap_ci,
  here::here(paste0("results/", setting, "/bootstrap_ci_", setting, ".rds"))
)
