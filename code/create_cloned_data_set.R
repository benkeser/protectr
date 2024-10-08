#' Create cloned data sets used to fit marginal structural models
#' 
#' Given weekly data with probability weights added in, this function
#' creates "cloned" data sets that are used to fit MSMs.
#' 
#' Up to three data sets are created for the TB outcome, death outcome,
#' and death outcome where TB events are treated as censored (needed for
#' calculation of the cumulative incidence of TB treating death as a 
#' competing risk).
#' 
#' @param propensity_output
#' @param endpoint a \code{vector} of strings indicating which of the three
#'   outcomes to create a cloned data set for
#' 
#' @return a \code{list} with named elements for each of the possible \code{endpoint}
#'   values. If an endpoint is not included, then the element is \code{NULL} in the output

create_cloned_data_set <- function(
	propensity_output,
	endpoint = c("tb", "death", "death_for_tb"),
	mdr = "no",
	...
){

	assertthat::assert_that(all(endpoint %in% c("tb", "death", "death_for_tb")))

	if("tb" %in% endpoint){
		weekly_records_tb_z1 <- propensity_output$weekly_records_data[,
			.SD[(wk <= pmin(tb_wk, right_cens_wk_tb) & tpt_start_wk <= propensity_output$grace_pd_wks) |
				  ((tpt_start_wk > propensity_output$grace_pd_wks) & (wk <= propensity_output$grace_pd_wks) & (wk <= pmin(tb_wk, right_cens_wk_tb)))],
			by = id
		]
	  weekly_records_tb_z1[, z := 1]
	}
	
	if("death" %in% endpoint){
		weekly_records_death_z1 <- propensity_output$weekly_records_data[
		  (tpt_start_wk <= propensity_output$grace_pd_wks) |
	    ((tpt_start_wk > propensity_output$grace_pd_wks) & (wk <= propensity_output$grace_pd_wks))
	  ]
	  weekly_records_death_z1[, z := 1]
	}

	if("death_for_tb" %in% endpoint){
		weekly_records_death_for_tb_z1 <- propensity_output$weekly_records_data[,
		  .SD[((tpt_start_wk <= propensity_output$grace_pd_wks) & (wk <= pmin(death_wk, right_cens_wk_tb)) & (wk <= pmin(tb_wk, right_cens_wk_tb))) |
	    		((tpt_start_wk > propensity_output$grace_pd_wks) & (wk <= propensity_output$grace_pd_wks) & (wk <= pmin(death_wk, right_cens_wk_tb)) & (wk <= pmin(tb_wk, right_cens_wk_tb)))],
	    by = id
	  ]
	  weekly_records_death_for_tb_z1[, z := 1]
	}

	# for the "half" of the data set where z = 0
	if("tb" %in% endpoint){
		weekly_records_tb_z0 <- propensity_output$weekly_records_data[,
		  .SD[((tpt_start_wk < 99999) & (wk <= tpt_start_wk) & (wk <= pmin(tb_wk, right_cens_wk_tb))) |
		  	  ((tpt_start_wk == 99999) & (wk <= pmin(tb_wk, right_cens_wk_tb)))],
		  by = id
		]
		weekly_records_tb_z0[, z := 0]
	}

	if("death" %in% endpoint){
		weekly_records_death_z0 <- propensity_output$weekly_records_data[
		  ((tpt_start_wk < 99999) & (wk <= tpt_start_wk)) |
		  (tpt_start_wk == 99999)
		]
		weekly_records_death_z0[, z := 0]
	}

	if("death_for_tb" %in% endpoint){
		weekly_records_death_for_tb_z0 <- propensity_output$weekly_records_data[,
		  .SD[((tpt_start_wk < 99999) & (wk <= tpt_start_wk)) & (wk <= min(death_wk, right_cens_wk_tb)) & (wk <= min(tb_wk, right_cens_wk_tb)) |
		  		((tpt_start_wk == 99999) & (wk <= min(death_wk, right_cens_wk_tb)) & (wk <= min(tb_wk, right_cens_wk_tb)))],
		  by = id
		]
		weekly_records_death_for_tb_z0[, z := 0]
	}

	# combining the z = 1 and z = 0 data
	if("tb" %in% endpoint){
		cloned_data_tb <- rbindlist(list(weekly_records_tb_z1, weekly_records_tb_z0))

		if(mdr == "no"){
			cloned_data_tb[, dN := (wk == tb_wk)]
		}else{
			cloned_data_tb[, dN_mdr := (wk == tb_wk) & (mdr == 1) ]
			cloned_data_tb[, dN_not_mdr := (wk == tb_wk) & (mdr == 0)]
		}

		cloned_data_tb[, wt_k := ifelse(z == 1, wt_tpt_tb, wt_cntrl_tb)]
	}else{
		cloned_data_tb <- NULL
	}

	# creating the death outcome
	if("death" %in% endpoint){
		cloned_data_death <- rbindlist(list(weekly_records_death_z1, weekly_records_death_z0))
		cloned_data_death[, dN := (wk == death_wk)]
		cloned_data_death[, wt_k := ifelse(z == 1, wt_tpt_death, wt_cntrl_death)]
	}else{
		cloned_data_death <- NULL
	}

	if("death_for_tb" %in% endpoint){
		cloned_data_death_for_tb <- rbindlist(list(weekly_records_death_for_tb_z1, weekly_records_death_for_tb_z0))
		cloned_data_death_for_tb[, dN := (wk == death_wk)]
		cloned_data_death_for_tb[, wt_k := ifelse(z == 1, wt_tpt_death, wt_cntrl_death)]
	}else{
		cloned_data_death_for_tb <- NULL
	}

	out <- list(
		tb = cloned_data_tb,
		death = cloned_data_death,
		death_for_tb = cloned_data_death_for_tb
	)

	return(out)
}
