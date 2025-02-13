#' Helper function to run one iteration of bootstrap
#'
#' The function samples with replacement from \code{weekly_records_data}
#' and performs a complete analysis of the bootstrap data.
#'
#' The analysis includes estimating the propensity models, creating cloned
#' data sets, and fitting (possibly multiple) MSMs. The output of these MSMs
#' is returned.
#'
#' The analysis additionally summarizes the output of the MSMs to compute
#' cumulative incidence that can be used for plotting functions.
#'
#' The SAP for the primary PROTECT analysis specifies that two MSMs will be
#' fit for each outcome: a simple MSM (similar to a Cox-model) that will be
#' used to compute a log-rank-like test statistic and a more flexible MSM that
#' will be used to visualize curves in a more nonparametric fashion.
#'
#' @param weekly_records_data a \code{data.table} of the original
#'   weekly records for each client
#' @param grace_pd_wks a \code{numeric} specifying the length of
#'   the grace period
#' @param denom_model_formula a \code{string} specifying the
#'   right-hand side of the formula for the TPT initiation model
#'   used in the denominator of the weights calculation
#' @param num_model_formula a \code{string} specifying the
#'   right-hand side of the formula for the TPT intiation model
#'   used in the numerator of the weights calculation
#' @param right_cens_model_formula a \code{string} specifying the
#'   right-hand side of the formula for the right-censoring model
#'   used in the weights calculation
#' @param msm_formulas_tb a \code{vector} or \code{strings} specifying
#'   the right-hand side of the formula for the MSM for the TB outcome
#' @param msm_formulas_death a \code{vector} or \code{strings} specifying
#'   the right-hand side of the formula for the MSM for the mortality outcome
#' @param msm_formulas_death a \code{vector} or \code{strings} specifying
#'   the right-hand side of the formula for the MSM for the TB outcome
#' @param admin_cens_wks a \code{numeric} indicating the length of follow
#'   up prior to administrative censoring
#' @param ... other options (not currently used)
#'
#' @return a list with named elements:
#'   cf_init_dist = a \code{data.frame} with columns wk, haz, pmf, cdf
#'     corresponding respectively to the week, hazard for TPT initiation,
#'     PMF for TPT initiation, and CDF for TPT initiation
#'   msm_fits_tb = a \code{list} with elements named according the formula
#'     used to fit the MSM. Each element of the list includes the output of
#'     \code{fit_msm}
#'   msm_fits_death = same as above, but for death outcome
#'   msm_fits_death_for_tb = same as above, but for the death outcome where
#'     TB cases are treated as being censored
#'   cuminc = a list with named elements, as above. each element contains
#'     a data.frame with information on the cumulative incidence implied by
#'     each of the MSM fits

try_one_bootstrap <- function(
    boot_num, # added for sapply
    weekly_records_data,
    grace_pd_wks = 8,
    denom_model_formula = "1",
    num_model_formula = denom_model_formula,
    right_cens_model_formula = "1",
    msm_formulas_tb = "splines::ns(wk, 3) + z",
    msm_formulas_death = "splines::ns(wk, 3) + z",
    msm_formulas_death_for_tb = msm_formulas_death,
    admin_cens_wks = 52 * 2) {
    tryCatch(
        {
            return(
                do_one_bootstrap(
                    weekly_records_data = weekly_records_data,
                    grace_pd_wks = grace_pd_wks,
                    denom_model_formula = denom_model_formula,
                    num_model_formula = num_model_formula,
                    right_cens_model_formula = right_cens_model_formula,
                    msm_formulas_tb = msm_formulas_tb,
                    msm_formulas_death = msm_formulas_death,
                    msm_formulas_death_for_tb = msm_formulas_death_for_tb,
                    admin_cens_wks = admin_cens_wks
                )
            )
        },
        error = function(e) {
            return(
                list(
                    cf_init_dist = NULL,
                    msm_fits_tb = NULL,
                    msm_fits_death = NULL,
                    msm_fits_death_for_tb = NULL,
                    cuminc = NULL,
                    error = TRUE ,
                    error_msg = e
                )
            )
        }
    )
}

do_one_bootstrap <- function(
    weekly_records_data,
    grace_pd_wks = 8,
    denom_model_formula = "1",
    num_model_formula = denom_model_formula,
    right_cens_model_formula = "1",
    msm_formulas_tb = "splines::ns(wk, 3) + z",
    msm_formulas_death = "splines::ns(wk, 3) + z",
    msm_formulas_death_for_tb = msm_formulas_death,
    admin_cens_wks = 52 * 2,
    ...) {
    
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

    propensity_output <- fit_propensity_models(
        weekly_records_data_bootstrap,
        grace_pd_wks = grace_pd_wks,
        denom_model_formula = denom_model_formula,
        num_model_formula = num_model_formula,
        right_cens_model_formula = right_cens_model_formula,
        return_models = FALSE
    )

    cf_init_dist <- compute_cf_init_dist(
        propensity_output = propensity_output
    )

    cloned_data_sets <- create_cloned_data_set(
        propensity_output = propensity_output
    )

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

    out <- list(
        cf_init_dist = list(cf_init_dist),
        msm_fits_tb = msm_fits_tb,
        msm_fits_death = msm_fits_death,
        msm_fits_death_for_tb = msm_fits_death_for_tb,
        cuminc = cuminc,
        error = FALSE
    )

    return(out)
}


run_bootstrap <- function(
    nboot = 1e3,
    weekly_records_data,
    grace_pd_wks = 8,
    denom_model_formula = "1",
    num_model_formula = denom_model_formula,
    right_cens_model_formula = "1",
    msm_formulas_tb = "splines::ns(wk, 3) + z",
    msm_formulas_death = "splines::ns(wk, 3) + z",
    msm_formulas_death_for_tb = msm_formulas_death,
    admin_cens_wks = 52 * 2) {
    
	# bootstrap_results <-
    #     future.apply::future_replicate(
    #         nboot, try_one_bootstrap(
    #             weekly_records_data = weekly_records_data,
    #             grace_pd_wks = grace_pd_wks,
    #             denom_model_formula = denom_model_formula,
    #             num_model_formula = num_model_formula,
    #             right_cens_model_formula = right_cens_model_formula,
    #             msm_formulas_tb = msm_formulas_tb,
    #             msm_formulas_death = msm_formulas_death,
    #             msm_formulas_death_for_tb = msm_formulas_death_for_tb,
    #             admin_cens_wks = admin_cens_wks
    #         ),
	# 		future.globals = c("weekly_records_data", 
	# 		"fit_propensity_models",
	# 		"create_cloned_data_set",
	# 		"fit_msm",
	# 		"compute_cf_init_dist", 
	# 		"compute_cuminc")
    #     )
  
  
  bootstrap_results <-
    future.apply::future_sapply(
      1:nboot, try_one_bootstrap,
      weekly_records_data = weekly_records_data,
      grace_pd_wks = grace_pd_wks,
      denom_model_formula = denom_model_formula,
      num_model_formula = num_model_formula,
      right_cens_model_formula = right_cens_model_formula,
      msm_formulas_tb = msm_formulas_tb,
      msm_formulas_death = msm_formulas_death,
      msm_formulas_death_for_tb = msm_formulas_death_for_tb,
      admin_cens_wks = admin_cens_wks,
      future.globals = c("weekly_records_data", 
                         "fit_propensity_models",
                         "create_cloned_data_set",
                         "fit_msm",
                         "compute_cf_init_dist", 
                         "compute_cuminc")
    )

    return(bootstrap_results)
}

get_bootstrap_ci <- function(bootstrap_results) {
    # confidence intervals for cf_init_dist
    cf_init_dist_ci <- get_fn_by_wk_boot_result(
        bootstrap_results, "cf_init_dist"
    )

    # bootstrap results for msm_fits_tb
    msm_fits_tb_ci <- get_msm_boot_result(
        bootstrap_results, "msm_fits_tb"
    )

    # bootstrap results for msm_fits_death
    msm_fits_death_ci <- get_msm_boot_result(
        bootstrap_results, "msm_fits_death"
    )

    # bootstrap results for msm_fits_death_for_tb
    msm_fits_death_for_tb_ci <- get_msm_boot_result(
        bootstrap_results, "msm_fits_death_for_tb"
    )

    # bootstrap results for cuminc
    n_models <- length(bootstrap_results["cuminc", 1])
    cuminc_ci <- get_fn_by_wk_boot_result(
        bootstrap_results = bootstrap_results,
        fn_by_wk_name = "cuminc"
    )

    out <- list(
        cf_init_dist_ci = cf_init_dist_ci,
        msm_fits_tb_ci = msm_fits_tb_ci,
        msm_fits_death_ci = msm_fits_death_ci,
        msm_fits_death_for_tb_ci = msm_fits_death_for_tb_ci,
        cuminc_ci = cuminc_ci
    )

    return(out)
}

get_fn_by_wk_boot_result <- function(bootstrap_results,
                                     fn_by_wk_name = "cf_init_dist",
                                     subset_idx = NULL) {
    all_fn_by_wk_list <- bootstrap_results[fn_by_wk_name, ]
    n_models <- length(all_fn_by_wk_list[[1]])
    out <- vector(mode = "list", length = n_models)

    for (i in seq_len(n_models)) {
        if (fn_by_wk_name == "cf_init_dist") {
            # this should always be length 1
            n_subgrp <- 1
        } else {
            # if effect hetero this will have length > 1
            # else length == 1
            n_subgrp <- length(all_fn_by_wk_list[[1]][[i]])
            names(out)[i] <- paste0(names(all_fn_by_wk_list[[1]][i]))
        }

        out_j <- vector(mode = "list", length = n_subgrp)
        for (j in seq_len(n_subgrp)) {
            if (fn_by_wk_name == "cf_init_dist") {
                all_fn_by_wk_df_list <- list(lapply(all_fn_by_wk_list, "[[", i))
                all_fn_by_wk_df <- all_fn_by_wk_df_list[[j]]
            } else {
                all_fn_by_wk_df <- lapply(all_fn_by_wk_list, function(x) {
                    x[[i]][[j]]
                })
            }


            if (!is.null(subset_idx)) {
                all_fn_by_wk_df <- lapply(
                    all_fn_by_wk_df, "[[", subset_idx
                )
            }

            all_fn_by_wk_df_no_wk <- lapply(
                all_fn_by_wk_df, function(x) {
                    x[, -1] # remove week column
                }
            )

            colnames_fn_by_wk <- colnames(all_fn_by_wk_df_no_wk[[1]])

            fn_by_wk_by_col <- sapply(colnames_fn_by_wk, function(x) {
                lapply(all_fn_by_wk_df_no_wk, "[[", x)
            }, simplify = FALSE, USE.NAMES = TRUE)

            fn_by_wk_by_col_all_boot <- lapply(fn_by_wk_by_col, Reduce, f = cbind)

            fn_by_wk_ci_by_col <- lapply(fn_by_wk_by_col_all_boot, function(x) {
                t(apply(x, 1, quantile, p = c(0.025, 0.975), na.rm = TRUE))
            })
            for (x in names(fn_by_wk_ci_by_col)) {
                colnames(fn_by_wk_ci_by_col[[x]]) <- paste0(x, c("_cil", "_ciu"))
            }
            fn_by_wk_ci_matrix <- Reduce(cbind, fn_by_wk_ci_by_col)
            fn_by_wk_ci <- data.frame(
                wk = all_fn_by_wk_df[[1]]$wk,
                fn_by_wk_ci_matrix
            )
            out_j[[j]] <- fn_by_wk_ci
            if (n_subgrp == 1) {
                names(out_j)[j] <- paste0(names(all_fn_by_wk_list[[1]][i]))
            } else {
                names(out_j)[j] <- paste0(names(all_fn_by_wk_list[[1]][[i]][j]))
            }
        }
        out[[i]] <- out_j
    }
    return(out)
}

get_msm_boot_result <- function(bootstrap_results,
                                msm_fits_name = "msm_fits_tb") {
    all_z_coef <- sapply(
        bootstrap_results[msm_fits_name, ],
        function(y) {
            unlist(lapply(y, function(yy) {
                if (length(yy$msm_coef) == 5) { # this is 5 because models without effect heterogeneity will have 5 terms (intercept, 3 wk splines, z)
                    yy$msm_coef["z"]
                } else {
                    c(yy$msm_coef["z"], (yy$msm_coef["z"] + yy$msm_coef[str_detect(names(yy$msm_coef), "z:")]))
                }
            }))
        },
        simplify = FALSE, USE.NAMES = TRUE
    )

    all_z_coef_matrix <- Reduce(rbind, all_z_coef)
    all_z_coef_ci <- t(apply(all_z_coef_matrix, 2, quantile, p = c(0.025, 0.975)))
    all_z_coef_sd <- t(apply(all_z_coef_matrix, 2, sd))
    msm_fits_ci <- data.frame(
        se = c(all_z_coef_sd),
        ci_l = all_z_coef_ci[, 1],
        ci_u = all_z_coef_ci[, 2]
    )
    row.names(msm_fits_ci) <- NULL
    msm_fits_ci_split <- split(msm_fits_ci, seq_len(dim(msm_fits_ci)[1]))
    names(msm_fits_ci_split) <- names(all_z_coef[[1]])
    return(msm_fits_ci_split)
}
