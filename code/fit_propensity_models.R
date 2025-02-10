#' Fit propensity models and add weights to data set
#' 
#' This function estimates TPT initiation probabilities and
#' right-censoring probabilities that are used to create weights
#' needed for the marginal structural model.
#' 
#' @param weekly_records_data a \code{data.table} as outputted by
#'   \code{create_weekly_records} function
#' @param grace_pd_wks a \code{numeric} indicating the length of 
#'   the grace period
#' @param denom_model_formula a \code{string} with the right-hand side
#'   of the model formula used to estimate TPT initiation probabilities
#'   used in the denominator of the probability weights
#' @param num_model_formula a \code{string} with the right-hand side
#'   of the model formula used to estimate TPT initiation probabilities
#'   used in the numerator of the probability weights. Unless there is a
#'   good reason to, this should use the same model formula as the 
#'   denominator model
#' @param right_cens_model_formula a \code{string} with the right-hand side
#'   of the model formula used to estimate right-censoring probabilities
#' @param return_models a \code{boolean} indicating whether to return estimated
#'   model fits as part of the output
#' 
#' @return a named \code{list} with elements: \code{weekly_records_data} = same
#'   as the input data but with estimated probabilities and weights added; 
#'   \code{models} = a \code{list} of model fits if \code{return_models};
#'   \code{grace_pd_wks} = passed through from the input

fit_propensity_models <- function(
	weekly_records_data,
	grace_pd_wks = 8,
	denom_model_formula = "1",
	num_model_formula = denom_model_formula,
	right_cens_model_formula = "1",
	return_models = TRUE,
	...
){
	
	#----------------------------------------------
	# DENOMINATOR 
	#----------------------------------------------

	denom_model_data <- weekly_records_data[
	  wk <= tpt_start_wk &
	  wk <= tb_wk &
	  wk <= death_wk &
	  wk <= right_cens_wk_tb
	]

	# create an outcome variable column
	denom_model_data[, tpt_outcome := as.numeric(wk == tpt_start_wk)]

	# regress the outcome variable against covariates
	denom_model_form <- paste0("tpt_outcome ~ ", denom_model_formula)
	denom_model <- glm(
		denom_model_form,
		family = stats::binomial(),
		data = denom_model_data
	)
	weekly_records_data[, denom_model_prediction := predict(denom_model, newdata = weekly_records_data, type = "response")]

	# TPT weights
	# weights for TB outcome
	weekly_records_data[, prob_wt_denom_tpt_tb := -99999]

	# people who initiate TPT in the grace period
	weekly_records_data[(wk <  tpt_start_wk) & (tpt_start_wk <= grace_pd_wks), prob_wt_denom_tpt_tb := 1 - denom_model_prediction]
	weekly_records_data[(wk == tpt_start_wk) & (tpt_start_wk <= grace_pd_wks), prob_wt_denom_tpt_tb := denom_model_prediction]
	weekly_records_data[(wk >  tpt_start_wk) & (tpt_start_wk <= grace_pd_wks), prob_wt_denom_tpt_tb := 1]

	# people who initiate TPT outside the grace period or never
	weekly_records_data[(wk <= grace_pd_wks) & (tpt_start_wk > grace_pd_wks), prob_wt_denom_tpt_tb := 1 - denom_model_prediction]
	weekly_records_data[(wk >  grace_pd_wks) & (tpt_start_wk > grace_pd_wks), prob_wt_denom_tpt_tb := 99999]
  
  # wks after case of TB should not be included
	weekly_records_data[(wk > tb_wk), prob_wt_denom_tpt_tb := 99999]
  
	# weights for death outcome
	weekly_records_data[, prob_wt_denom_tpt_death := -99999]
	# people who initiate TPT in the grace period
	weekly_records_data[(wk <  tpt_start_wk) & (tpt_start_wk <= grace_pd_wks), prob_wt_denom_tpt_death := 1 - denom_model_prediction]
	weekly_records_data[(wk == tpt_start_wk) & (tpt_start_wk <= grace_pd_wks), prob_wt_denom_tpt_death := denom_model_prediction]
	weekly_records_data[(wk >  tpt_start_wk) & (tpt_start_wk <= grace_pd_wks), prob_wt_denom_tpt_death := 1]

	# people who initiate TPT outside the grace period
	weekly_records_data[(wk <= grace_pd_wks) & (tpt_start_wk > grace_pd_wks) & (tpt_start_wk < 99999), prob_wt_denom_tpt_death := 1 - denom_model_prediction]
	weekly_records_data[(wk >  grace_pd_wks) & (tpt_start_wk > grace_pd_wks) & (tpt_start_wk < 99999), prob_wt_denom_tpt_death := 99999]
  
  # people who get TB in the grace period and so never initiate TPT
  weekly_records_data[(wk <= tb_wk) & (tb_wk <= grace_pd_wks), prob_wt_denom_tpt_death := 1 - denom_model_prediction]
  weekly_records_data[(wk >  tb_wk) & (tb_wk <= grace_pd_wks), prob_wt_denom_tpt_death := 1]
  
  # people who either (get TB outside the grace period or never get TB) and never initiate TPT
  weekly_records_data[(wk <= grace_pd_wks) & (tb_wk > grace_pd_wks) & (tpt_start_wk == 99999), prob_wt_denom_tpt_death := 1 - denom_model_prediction]
  weekly_records_data[(wk >  grace_pd_wks) & (tb_wk > grace_pd_wks) & (tpt_start_wk == 99999), prob_wt_denom_tpt_death := 99999]
  

  # control weights
  weekly_records_data[, prob_wt_denom_cntrl_tb := 99999]
	
	# these weeks should not be included in any z = 0 data set, so fill in 
	# with 99999 to trigger a flag extreme weights if these rows are accidentally included
	weekly_records_data[(wk >= tpt_start_wk), prob_wt_denom_cntrl_tb := 99999]
	weekly_records_data[(wk > tb_wk), prob_wt_denom_cntrl_tb := 99999]
	
	# for any week prior to tpt start, the weight is 1 - prediction for those not on TPT 
	weekly_records_data[(wk < tpt_start_wk & wk <= tb_wk), prob_wt_denom_cntrl_tb := 1 - denom_model_prediction]


	weekly_records_data[, prob_wt_denom_cntrl_death := 99999]
	
	# same logic as for cntrl_tb above
	weekly_records_data[(wk >= tpt_start_wk), prob_wt_denom_cntrl_death := 99999]
	weekly_records_data[(wk > death_wk), prob_wt_denom_cntrl_death := 99999]
	
	weekly_records_data[(wk < tpt_start_wk & wk <= death_wk), prob_wt_denom_cntrl_death := 1 - denom_model_prediction]
	weekly_records_data[(wk > tb_wk), prob_wt_denom_cntrl_death := 1]
	

	#----------------------------------------------
	# NUMERATOR 
	#----------------------------------------------
	num_model_data <- denom_model_data[tpt_start_wk < grace_pd_wks]

	num_model_form <- paste0("tpt_outcome ~ ", num_model_formula)
	num_model <- glm(
		num_model_form,
		family = stats::binomial(),
		data = num_model_data
	)

	weekly_records_data[, num_model_prediction := predict(num_model, newdata = weekly_records_data, type = "response")]

	# if not initiated prior to the end of the grace period, prob of
	# initiating in last week of grace period = 1
	weekly_records_data[wk == grace_pd_wks, num_model_prediction := 1]

	weekly_records_data[, prob_wt_num_tpt_tb := 99999]
	weekly_records_data[(wk > grace_pd_wks), prob_wt_num_tpt_tb := 1]
	weekly_records_data[(wk > tpt_start_wk), prob_wt_num_tpt_tb := 1]

	# none of these weeks should receive any weight
	weekly_records_data[(wk > tb_wk), prob_wt_num_tpt_tb := 0]
	weekly_records_data[(wk <  tpt_start_wk & wk <= tb_wk & wk <= grace_pd_wks), prob_wt_num_tpt_tb := 1 - num_model_prediction]
	weekly_records_data[(wk == tpt_start_wk & wk <= tb_wk & wk <= grace_pd_wks), prob_wt_num_tpt_tb := num_model_prediction]


	weekly_records_data[, prob_wt_num_tpt_death := 99999]
	weekly_records_data[(wk > grace_pd_wks), prob_wt_num_tpt_death := 1]
	weekly_records_data[(wk > tpt_start_wk), prob_wt_num_tpt_death := 1]

	weekly_records_data[(wk > death_wk), prob_wt_num_tpt_death := 0]
	weekly_records_data[(wk <  tpt_start_wk & wk <= death_wk & wk <= grace_pd_wks), prob_wt_num_tpt_death := 1 - num_model_prediction]
	weekly_records_data[(wk == tpt_start_wk & wk <= death_wk & wk <= grace_pd_wks), prob_wt_num_tpt_death := num_model_prediction]
	# special care for tb cases in the grace period -- these are still "compatible" with
	# our interventional TPT initiation distribution
	weekly_records_data[(tb_wk <= grace_pd_wks & wk > tb_wk & wk <= death_wk), prob_wt_num_tpt_death := 1]



	weekly_records_data[, prob_wt_num_cntrl_tb := NA]
	weekly_records_data[(wk > grace_pd_wks), prob_wt_num_cntrl_tb := 1]
	weekly_records_data[(wk > tpt_start_wk), prob_wt_num_cntrl_tb := 0]
	weekly_records_data[(wk > tb_wk), prob_wt_num_cntrl_tb := 0]
	weekly_records_data[(wk < tpt_start_wk & wk <= tb_wk), prob_wt_num_cntrl_tb := 1]
	weekly_records_data[(wk == tpt_start_wk & wk <= tb_wk), prob_wt_num_cntrl_tb := 0]

	weekly_records_data[, prob_wt_num_cntrl_death := NA]
	weekly_records_data[(wk > grace_pd_wks), prob_wt_num_cntrl_death := 1]
	weekly_records_data[(wk > tpt_start_wk), prob_wt_num_cntrl_death := 0]
	weekly_records_data[(wk > death_wk), prob_wt_num_cntrl_death := 0]
	weekly_records_data[(wk < tpt_start_wk & wk <= death_wk), prob_wt_num_cntrl_death := 1]
	weekly_records_data[(wk == tpt_start_wk & wk <= death_wk), prob_wt_num_cntrl_death := 0]


  # RIGHT-CENSORING weights for tb
	cens_model_data_tb <- weekly_records_data[wk < admin_cens_wk]

	cens_model_data_tb[, cens_outcome := (wk == right_cens_wk_tb) & ((tb_wk == 99999) | (death_wk == 99999))]

	# regress the outcome variable against covariates
	cens_model_form <- paste0("cens_outcome ~ ", right_cens_model_formula)
	cens_model_tb <- glm(
		cens_model_form,
		family = stats::binomial(),
		data = cens_model_data_tb
	)

	fitted_values <- predict(cens_model_tb, newdata = weekly_records_data, type = "response")
	uncens_probs <- c(1, 1 - fitted_values[1:(length(fitted_values)-1)])
	uncens_probs[weekly_records_data$wk == 1] <- 1

	weekly_records_data[, prob_wt_cens_tb := uncens_probs]


  # RIGHT-CENSORING weights for death
	cens_model_data_death <- weekly_records_data[wk < admin_cens_wk]

	cens_model_data_death[, cens_outcome := (wk == last_visit_wk) & (death_wk == 99999)]

	# regress the outcome variable against covariates
	cens_model_form <- paste0("cens_outcome ~ ", right_cens_model_formula)
	cens_model_death <- glm(
		cens_model_form,
		family = stats::binomial(),
		data = cens_model_data_death
	)

	fitted_values <- predict(cens_model_death, newdata = weekly_records_data, type = "response")
	uncens_probs <- c(1, 1 - fitted_values[1:(length(fitted_values)-1)])
	uncens_probs[weekly_records_data$wk == 1] <- 1

	weekly_records_data[, prob_wt_cens_death := uncens_probs]


	# 3) Turn calculated columns into appropriate weights to be used in the MSM
	
	# NEW split beforehand to use less memory + only required columns (thanks gpt)
	
	# give idx to rejoin later
	weekly_records_data$idx <- 1:nrow(weekly_records_data)
	
	sub_weekly <- weekly_records_data[, .(id, 
	                                      idx,
	                                      prob_wt_num_tpt_tb,
	                                      prob_wt_denom_tpt_tb,
	                                      prob_wt_cens_tb,
	                                      prob_wt_num_tpt_death,
	                                      prob_wt_denom_tpt_death,
	                                      prob_wt_cens_death,
	                                      prob_wt_denom_cntrl_tb,
	                                      prob_wt_denom_cntrl_death)]
	
	
	weekly_records_data_this_id <- split(sub_weekly, sub_weekly$id)
	
	wts_by_id <- future.apply::future_lapply(
	  weekly_records_data_this_id, 
	  function(data_chunk) {
  	  # Weight calculations
  	  data_chunk[, wt_tpt_tb := cumprod(prob_wt_num_tpt_tb / (prob_wt_denom_tpt_tb * prob_wt_cens_tb))]
  	  data_chunk[, wt_tpt_death := cumprod(prob_wt_num_tpt_death / (prob_wt_denom_tpt_death * prob_wt_cens_death))]
  	  data_chunk[, wt_cntrl_tb := cumprod(1 / (prob_wt_denom_cntrl_tb * prob_wt_cens_tb))]
  	  data_chunk[, wt_cntrl_death := cumprod(1 / (prob_wt_denom_cntrl_death * prob_wt_cens_death))]
  	  
  	  return(data_chunk[, .(id, idx, wt_tpt_tb, wt_tpt_death, wt_cntrl_tb, wt_cntrl_death)])
	  },
	  future.envir = baseenv(),
	  future.packages = c("data.table")
	)
	
	weekly_records_data_wts <- rbindlist(wts_by_id)
	
	# rejoin weights data with original weekly records data (needed later)
	weekly_records_data <- merge(
	  weekly_records_data, 
	  weekly_records_data_wts, 
	  by = c("idx", "id"), 
	  all.x = TRUE
	)

	out <- list()
	out$weekly_records_data <- weekly_records_data_wts
	out$models <- list(
		denom_model = NULL,
		num_model = NULL,
		cens_model_tb = NULL,
		cens_model_death = NULL
	)
	out$grace_pd_wks <- grace_pd_wks

	if(return_models){
		out$models$denom_model <- strip_glm(denom_model)
		out$models$num_model <- strip_glm(num_model)
		out$models$cens_model_tb <- strip_glm(cens_model_tb)
		out$models$cens_model_death <- strip_glm(cens_model_death)
	}

	return(out)
}


strip_glm <- function(model) {
  model$env <- NULL
  model$y = c()
  model$model = c()

  model$residuals = c()
  model$fitted.values = c()
  model$effects = c()
  model$linear.predictors = c()
  model$weights = c()
  model$prior.weights = c()
  model$data = c()

  model$family$variance = c()
  model$family$dev.resids = c()
  model$family$aic = c()
  model$family$validmu = c()
  model$family$simulate = c()

  attr(model$terms,".Environment") = c()
  attr(model$formula,".Environment") = c()

  return(model)
}