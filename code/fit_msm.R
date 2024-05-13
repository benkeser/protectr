fit_msm <- function(
	cloned_data_set, 
	msm_formula = "splines::ns(wk, 3) + z",
	gee = FALSE,
	return_msm_model = FALSE,
	return_msm_vcov = FALSE,
	...
){
	msm_form <- paste0("dN ~ ", msm_formula)
	if(!gee){
		msm_fit <- glm(
			msm_form,
			data = cloned_data_set,
			family = stats::binomial(),
			weights = wt_k
		)
		msm_fit <- strip_glm(msm_fit)
	}else{
		msm_fit <- geepack::geeglm(
			msm_form,
			data = clone_data_set,
			id = id,
			family = stats::binomial(),
			weights = wt_k,
			corstr = "independence"
		)
		# TODO: strip model?
	}

	msm_coef <- msm_fit$coefficients

	if(return_vcov){
		msm_vcov <- vcov(msm_fit)
	}else{
		msm_vcov <- NULL
	}

	out <- list(
		msm_coef = msm_coef,
		msm_vcov = msm_vcov,
		msm_model = NULL
	)
	if(return_msm_model){
		out$msm_model <- msm_fit
	}
	
	return(out)
}