compute_cuminc_ <- function(
	df_z1, df_z0, max_wk
){
	haz_tb_z1 <- predict(msm_fit_tb$msm_model, newdata = df_z1, type = "response")
	haz_tb_z0 <- predict(msm_fit_tb$msm_model, newdata = df_z0, type = "response")
	
	haz_death_z1 <- predict(msm_fit_death$msm_model, newdata = df_z1, type = "response")
	haz_death_z0 <- predict(msm_fit_death$msm_model, newdata = df_z0, type = "response")

	haz_death_for_tb_z1 <- predict(msm_fit_death_for_tb$msm_model, newdata = df_z1, type = "response")
	haz_death_for_tb_z0 <- predict(msm_fit_death_for_tb$msm_model, newdata = df_z0, type = "response")
	
	allcause_haz_tb_z1 <- haz_death_for_tb_z1 + haz_tb_z1
	allcause_haz_tb_z0 <- haz_death_for_tb_z0 + haz_tb_z0
	
	allcause_surv_z1 <- c(1, cumprod(1 - allcause_haz_tb_z1[-max_wk]))
	allcause_surv_z0 <- c(1, cumprod(1 - allcause_haz_tb_z0[-max_wk]))

	cif_tb_z1 <- cumsum(haz_tb_z1 * allcause_surv_z1)
	cif_tb_z0 <- cumsum(haz_tb_z0 * allcause_surv_z0)

	cif_death_z1 <- 1 - cumprod(1 - haz_death_z1)
	cif_death_z0 <- 1 - cumprod(1 - haz_death_z0)

	out <- data.frame(
		wk = seq_len(max_wk),
		haz_tb_z1 = haz_tb_z1,
		haz_tb_z0 = haz_tb_z0,
		haz_death_z1 = haz_death_z1,
		haz_death_z0 = haz_death_z0,
		cif_tb_z1 = cif_tb_z1,	                  
		cif_tb_z0 = cif_tb_z0,
		cif_death_z1 = cif_death_z1,
		cif_death_z0 = cif_death_z0                  
	)
	return(out)
}

compute_cuminc <- function(
	msm_fit_tb,
	msm_fit_death,
	msm_fit_death_for_tb,
	max_wk
){

	effect_hetero_variable_present <- (msm_fit_tb$effect_hetero_variable != "none")

	if(!effect_hetero_variable_present){
		df_z1 <- data.frame(z = 1, wk = seq_len(max_wk))
		df_z0 <- data.frame(z = 0, wk = seq_len(max_wk))
		out <- list(compute_cuminc_(df_z1, df_z0, max_wk))
	}else{
		out <- sapply(msm_fit_tb$effect_hetero_variable_vals, function(val){
			df_z1 <- data.frame(z = 1, wk = seq_len(max_wk), tmp = val)
			df_z0 <- data.frame(z = 0, wk = seq_len(max_wk), tmp = val)
			colnames(df_z1)[3] <- colnames(df_z0)[3] <- msm_fit_tb$effect_hetero_variable
			compute_cuminc_(df_z1, df_z0, max_wk)
		})
	}
	return(out)
}
