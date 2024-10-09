fit_msm <- function(
	cloned_data_set,
	msm_formula = "splines::ns(wk, 3) + z",
	outcome = "dN", # could add dN_mdr_tb, dN_not_mdr_tb
	gee = FALSE,
	return_msm_model = FALSE,
	return_msm_vcov = FALSE,
	...
){
	# if(is.null(msm_formula)){
	# 	msm_form <- paste0(outcome, " ~ ", baseline_haz_model, " + z")
	# 	if(effect_hetero_variable != "none"){
	# 		msm_form <- paste0(msm_form, "*", effect_hetero_variable)
	# 	}
	# }else{
	#   msm_form <- paste0(outcome, " ~ ", msm_formula)
	# }
	
	# look in msm_formula for z* or *z
	# split msm_formula on * and take the element that is not "z"
	# > str_split(msm_formula, "\\*")[[1]]
	# [1] "blah" "z"   
		
  assertthat::assert_that(startsWith(msm_formula, "splines::ns(wk, 3) + z") == TRUE)
  if(!is.na(str_split(msm_formula, "\\*")[[1]][2])){
    effect_hetero_variable <- str_split(msm_formula, "\\*")[[1]][2]
  }else{
    effect_hetero_variable <- NULL
  }
  
  msm_form <- paste0(outcome, " ~ ", msm_formula)
	
	# if(is.null(msm_formula)){
	# 	msm_form <- paste0("dN ~ ", baseline_haz_model, " + z")
	# 	if(effect_hetero_variable != "none"){
	# 		msm_form <- paste0(msm_form, "*", effect_hetero_variable)
	# 	}
	# }else{
	# 	msm_form <- paste0("dN ~ ", msm_formula)
	# }

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
			as.formula(msm_form),
			data = cloned_data_set,
			id = id,
			family = stats::binomial(),
			weights = wt_k,
			corstr = "independence"
		)
		# TODO: strip model?
	}

	msm_coef <- msm_fit$coefficients

	if(return_msm_vcov){
		msm_vcov <- vcov(msm_fit)
	}else{
		msm_vcov <- NULL
	}

	if(!is.null(effect_hetero_variable)){
		effect_hetero_variable_vals <- sort(
			unique(cloned_data_set[[effect_hetero_variable]])
		)
	}else{
		effect_hetero_variable_vals <- NULL
	}

	out <- list(
		msm_coef = msm_coef,
		msm_vcov = msm_vcov,
		msm_model = NULL,
		effect_hetero_variable = effect_hetero_variable,
		effect_hetero_variable_vals = effect_hetero_variable_vals
	)
	if(return_msm_model){
		out$msm_model <- msm_fit
	}
	
	return(out)
}
