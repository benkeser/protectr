compute_cf_init_dist <- function(propensity_output){
	
	assertthat::assert_that("num_model_prediction" %in% colnames(propensity_output$weekly_records_data))

	haz <- rep(NA, propensity_output$grace_pd_wks)
	for(k in seq_len(propensity_output$grace_pd_wks - 1)){
		haz[k] <- mean(
			propensity_output$weekly_records_data[wk == k, num_model_prediction], na.rm = TRUE
		)
	}
	# everyone not started by final week, starts in final week
	haz[propensity_output$grace_pd_wks] <- 1

	Sbar <- c(1, cumprod(1 - haz[-propensity_output$grace_pd_wks]))

	pmf <- haz * Sbar
	cdf <- cumsum(pmf)

	out <- data.frame(
    wk = seq_len(propensity_output$grace_pd_wks),
    haz = haz,
    pmf = pmf,
    cdf = cdf
	)
	return(out)
}