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

do_one_bootstrap <- function(
  weekly_records_data,
  grace_pd_wks,
  denom_model_formula,
  num_model_formula,
  right_cens_model_formula,   
  msm_formulas_tb = "splines::ns(wk, 3) + z",
  msm_formulas_death = "splines::ns(wk, 3) + z",
  msm_formulas_death_for_tb = msm_formulas_death,
  admin_cens_wks = 52 * 2
  ...                   
){

	sampled_ids <- sample(unique(weekly_records_data$id), replace = TRUE)
	
	weekly_records_data_bootstrap_by_id <- lapply(sampled_ids, function(this_id) {
  	return(weekly_records_data[id == this_id])
	})

	weekly_records_data_bootstrap <- rbindlist(
		weekly_records_data_bootstrap_by_id , idcol = "new_id"
	)
	setnames(weekly_records_data_bootstrap, old = "id", new = "orig_id")
	setnames(weekly_records_data_bootstrap, old = "new_id", new = "id")

	propensity_model_output <- fit_propensity_models(
		weekly_records_data,
		grace_pd_wks = grace_pd_wks,
		denom_model_formula = denom_model_formula,
		num_model_formula = num_model_formula,
		right_cens_model_formula = right_cens_model_formula,
		return_models = FALSE
	)
	
	cf_init_dist <- compute_cf_init_dist(
		weekly_records_data = propensity_model_output$weekly_records_data,
		grace_pd_wks = grace_pd_wks
	)

	cloned_data_sets <- create_cloned_data_set(
		weekly_records_data = propensity_model_output$weekly_records_data,
		grace_pd_wks = grace_pd_wks
	)

	msm_fits_tb <- sapply(msm_formulas_tb, 
		FUN = fit_msm,
		cloned_data_set = cloned_data_sets$tb, 
		return_msm_model = FALSE,
		return_msm_vcov = FALSE,
		simplify = FALSE,
		USE.NAMES = TRUE
	)

	msm_fits_death <- sapply(msm_formulas_death, 
		FUN = fit_msm,
		cloned_data_set = cloned_data_sets$death, 
		return_msm_model = FALSE,
		return_msm_vcov = FALSE,
		simplify = FALSE,
		USE.NAMES = TRUE
	)

	msm_fits_death_for_tb <- sapply(msm_formulas_death_for_tb,
		FUN = fit_msm, 
		cloned_data_set = cloned_data_sets$death_for_tb, 
		return_msm_model = FALSE,
		return_msm_vcov = FALSE,
		simplify = FALSE,
		USE.NAMES = TRUE
	)

	cuminc <- mapply(
		FUN = compute_cuminc,
		msm_fit_tb = msm_fits_tb,
		msm_fit_death = msm_fits_death,
		msm_fit_death_for_tb = msm_fits_death_for_tb,
		MoreArgs = list(max_wk = max_wk),
		simplify = FALSE,
		USE.NAMES = TRUE
	)

	out <- list(
		cf_init_dist = cf_init_dist,
		msm_fits_tb = msm_fits_tb,
		msm_fits_death = msm_fits_death,
		msm_fits_death_for_tb = msm_fits_death_for_tb,
		cuminc = cuminc
	)

	return(out)
}