compute_cuminc_ <- function(
	df_z1, df_z0, max_wk,
	msm_fit_tb,
	msm_fit_mdr_tb = NULL,
	msm_fit_not_mdr_tb = NULL,
	msm_fit_death,
	msm_fit_death_for_tb
){

	if(!is.null(msm_fit_tb)){
		haz_tb_z1 <- predict(msm_fit_tb$msm_model, newdata = df_z1, type = "response")
		haz_tb_z0 <- predict(msm_fit_tb$msm_model, newdata = df_z0, type = "response")
		haz_mdr_tb_z1 <- haz_mdr_tb_z0 <- haz_not_mdr_tb_z0 <- haz_not_mdr_tb_z1 <- NA
	}else{
		haz_mdr_tb_z1 <- predict(msm_fit_mdr_tb$msm_model, newdata = df_z1, type = "response")
		haz_not_mdr_tb_z1 <- predict(msm_fit_not_mdr_tb$msm_model, newdata = df_z1, type = "response")
		haz_mdr_tb_z0 <- predict(msm_fit_mdr_tb$msm_model, newdata = df_z0, type = "response")
		haz_not_mdr_tb_z0 <- predict(msm_fit_not_mdr_tb$msm_model, newdata = df_z0, type = "response")
		haz_tb_z1 <- haz_mdr_tb_z1 + haz_not_mdr_tb_z1
		haz_tb_z0 <- haz_mdr_tb_z0 + haz_not_mdr_tb_z0
	}

	haz_death_z1 <- predict(msm_fit_death$msm_model, newdata = df_z1, type = "response")
	haz_death_z0 <- predict(msm_fit_death$msm_model, newdata = df_z0, type = "response")

	haz_death_for_tb_z1 <- predict(msm_fit_death_for_tb$msm_model, newdata = df_z1, type = "response")
	haz_death_for_tb_z0 <- predict(msm_fit_death_for_tb$msm_model, newdata = df_z0, type = "response")
	
	# allcause = haz_death_for_tb + haz_mdr_tb + haz_non_mdr_tb
	allcause_haz_tb_z1 <- haz_death_for_tb_z1 + haz_tb_z1
	allcause_haz_tb_z0 <- haz_death_for_tb_z0 + haz_tb_z0
	
	allcause_surv_z1 <- c(1, cumprod(1 - allcause_haz_tb_z1[-max_wk]))
	allcause_surv_z0 <- c(1, cumprod(1 - allcause_haz_tb_z0[-max_wk]))

	cif_tb_z1 <- cumsum(haz_tb_z1 * allcause_surv_z1)
	cif_tb_z0 <- cumsum(haz_tb_z0 * allcause_surv_z0)
	if(!is.null(msm_fit_mdr_tb)){
		cif_mdr_tb_z1 <- cumsum(haz_mdr_tb_z1 * allcause_surv_z1)
		cif_mdr_tb_z0 <- cumsum(haz_mdr_tb_z0 * allcause_surv_z0)
		cif_not_mdr_tb_z1 <- cumsum(haz_not_mdr_tb_z1 * allcause_surv_z1)
		cif_not_mdr_tb_z0 <- cumsum(haz_not_mdr_tb_z0 * allcause_surv_z0)
	}else{
		cif_mdr_tb_z0 <- cif_mdr_tb_z1 <- cif_not_mdr_tb_z0 <- cif_not_mdr_tb_z1 <- NA
	}

	cif_death_z1 <- 1 - cumprod(1 - haz_death_z1)
	cif_death_z0 <- 1 - cumprod(1 - haz_death_z0)

	# add mdr to output, as appropriate
	out <- data.frame(
		wk = seq_len(max_wk),
		haz_tb_z1 = haz_tb_z1,
		haz_tb_z0 = haz_tb_z0,
		haz_mdr_tb_z1 = haz_mdr_tb_z1,
		haz_mdr_tb_z0 = haz_mdr_tb_z0,
		haz_not_mdr_tb_z1 = haz_not_mdr_tb_z1,
		haz_not_mdr_tb_z0 = haz_not_mdr_tb_z0,
		haz_death_z1 = haz_death_z1,
		haz_death_z0 = haz_death_z0,
		cif_tb_z1 = cif_tb_z1,	                  
		cif_mdr_tb_z1 = cif_mdr_tb_z1,	                  
		cif_not_mdr_tb_z1 = cif_not_mdr_tb_z1,	                  
		cif_tb_z0 = cif_tb_z0,
		cif_mdr_tb_z0 = cif_mdr_tb_z0,
		cif_not_mdr_tb_z0 = cif_not_mdr_tb_z0,
		cif_death_z1 = cif_death_z1,
		cif_death_z0 = cif_death_z0                  
	)
	return(out)
}

compute_cuminc <- function(
	msm_fit_tb,
	msm_fit_mdr_tb = NULL,
	msm_fit_not_mdr_tb = NULL,
	msm_fit_death,
	msm_fit_death_for_tb,
	max_wk
){

	if(!is.null(msm_fit_tb)){
		effect_hetero_variable_present <- (msm_fit_tb$effect_hetero_variable != "none")
		effect_hetero_variable <- msm_fit_tb$effect_hetero_variable
	}else{
		effect_hetero_variable_present <- (msm_fit_mdr_tb$effect_hetero_variable != "none")
		effect_hetero_variable <- msm_fit_mdr_tb$effect_hetero_variable
	}

	if(!effect_hetero_variable_present){
		df_z1 <- data.frame(z = 1, wk = seq_len(max_wk))
		df_z0 <- data.frame(z = 0, wk = seq_len(max_wk))
		out <- list(compute_cuminc_(df_z1 = df_z1, df_z0 = df_z0, max_wk = max_wk, 
			msm_fit_tb = msm_fit_tb,
			msm_fit_mdr_tb = msm_fit_mdr_tb,
			msm_fit_not_mdr_tb = msm_fit_not_mdr_tb,
			msm_fit_death = msm_fit_death,
			msm_fit_death_for_tb = msm_fit_death_for_tb)
		)
	}else{
		out <- sapply(msm_fit_tb$effect_hetero_variable_vals, simplify = FALSE, function(val){
			df_z1 <- data.frame(z = 1, wk = seq_len(max_wk), tmp = val)
			df_z0 <- data.frame(z = 0, wk = seq_len(max_wk), tmp = val)
			colnames(df_z1)[3] <- colnames(df_z0)[3] <- effect_hetero_variable
			compute_cuminc_(df_z1, df_z0, max_wk)
		}, msm_fit_tb = msm_fit_tb,
			msm_fit_mdr_tb = msm_fit_mdr_tb,
			msm_fit_not_mdr_tb = msm_fit_not_mdr_tb,
			msm_fit_death = msm_fit_death,
			msm_fit_death_for_tb = msm_fit_death_for_tb)
	}
	return(out)
}
