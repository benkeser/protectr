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

	  # truncating weights
	  weekly_records_tb_z1 <- truncate_wts(weekly_records_tb_z1, tb = "yes", tpt = "yes")
	}
	
	if("death" %in% endpoint){
		weekly_records_death_z1 <- propensity_output$weekly_records_data[,
			.SD[((((death_wk == 99999) & (wk <= last_visit_wk)) | ((death_wk < 99999) & (wk <= death_wk))) & tpt_start_wk <= propensity_output$grace_pd_wks) |
				  ((tpt_start_wk > propensity_output$grace_pd_wks) & (wk <= propensity_output$grace_pd_wks) & (((death_wk == 99999) & (wk <= last_visit_wk)) | ((death_wk < 99999) & (wk <= death_wk))))],
			by = id
	  ]
	  weekly_records_death_z1[, z := 1]

	  # truncating weights
	  weekly_records_death_z1 <- truncate_wts(weekly_records_death_z1, tb = "no", tpt = "yes")
	}

	if("death_for_tb" %in% endpoint){
		weekly_records_death_for_tb_z1 <- propensity_output$weekly_records_data[,
		  .SD[((tpt_start_wk <= propensity_output$grace_pd_wks) & (((death_wk == 99999) & (wk <= last_visit_wk)) | ((death_wk < 99999) & (wk <= death_wk))) & (wk <= pmin(tb_wk, right_cens_wk_tb))) |
	    		((tpt_start_wk > propensity_output$grace_pd_wks) & (wk <= propensity_output$grace_pd_wks) & (((death_wk == 99999) & (wk <= last_visit_wk)) | ((death_wk < 99999) & (wk <= death_wk))) & (wk <= pmin(tb_wk, right_cens_wk_tb)))],
	    by = id
	  ]
	  weekly_records_death_for_tb_z1[, z := 1]

	  # truncating weights
	  weekly_records_death_for_tb_z1 <- truncate_wts(weekly_records_death_for_tb_z1, tb = "no", tpt = "yes")
	}

	# for the "half" of the data set where z = 0
	if("tb" %in% endpoint){
		weekly_records_tb_z0 <- propensity_output$weekly_records_data[,
		  .SD[((tpt_start_wk < 99999) & (wk <= tpt_start_wk) & (wk <= pmin(tb_wk, right_cens_wk_tb))) |
		  	  ((tpt_start_wk == 99999) & (wk <= pmin(tb_wk, right_cens_wk_tb)))],
		  by = id
		]
		weekly_records_tb_z0[, z := 0]

		# truncating weights
	  	weekly_records_tb_z0 <- truncate_wts(weekly_records_tb_z0, tb = "yes", tpt = "no")
	}

	if("death" %in% endpoint){
		weekly_records_death_z0 <- propensity_output$weekly_records_data[,
		  .SD[((tpt_start_wk < 99999) & (wk <= tpt_start_wk) & (((death_wk == 99999) & (wk <= last_visit_wk)) | ((death_wk < 99999) & (wk <= death_wk)))) |
		  	  ((tpt_start_wk == 99999) & (((death_wk == 99999) & (wk <= last_visit_wk)) | ((death_wk < 99999) & (wk <= death_wk))))],
		  by = id
		]
		weekly_records_death_z0[, z := 0]

		# truncating weights
	  	weekly_records_death_z0 <- truncate_wts(weekly_records_death_z0, tb = "no", tpt = "no")
	}

	if("death_for_tb" %in% endpoint){
		weekly_records_death_for_tb_z0 <- propensity_output$weekly_records_data[,
		  .SD[((tpt_start_wk < 99999) & (wk <= tpt_start_wk)) & (((death_wk == 99999) & (wk <= last_visit_wk)) | ((death_wk < 99999) & (wk <= death_wk))) & (wk <= min(tb_wk, right_cens_wk_tb)) |
		  		((tpt_start_wk == 99999) & (((death_wk == 99999) & (wk <= last_visit_wk)) | ((death_wk < 99999) & (wk <= death_wk))) & (wk <= min(tb_wk, right_cens_wk_tb)))],
		  by = id
		]
		weekly_records_death_for_tb_z0[, z := 0]

		# truncating weights
	  	weekly_records_death_for_tb_z0 <- truncate_wts(weekly_records_death_for_tb_z0, tb = "no", tpt = "no")
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

truncate_wts <- function(weekly_records_data, tb, tpt){
  # Convert to data.table if not already
  setDT(weekly_records_data)
  
  # Summarizing the data table - for each ID, counting weeks, censoring, etc.
  if (tb == "yes") {
    wk_summary <- weekly_records_data[, .(wks = .N, 
                                          cens_wk = min(first(tb_wk), first(death_wk), 104)),
                                      by = id]
  } else {
    wk_summary <- weekly_records_data[, .(wks = .N, 
                                          cens_wk = min(first(death_wk), 104)),
                                      by = id]
  }
  
  # Calculate how many weeks to add
  wk_summary[, wks_to_add := cens_wk - wks]
  
  # Filter rows where wks_to_add > 0
  wk_summary <- wk_summary[wks_to_add > 0]
  
  # Creating new rows
  rows_add <- rbindlist(lapply(1:nrow(wk_summary), function(i) {
    this_id <- wk_summary$id[i]
    wks_to_add <- wk_summary$wks_to_add[i]
    last_wk <- wk_summary$wks[i]

    wk <- seq(from = last_wk + 1, by = 1, length.out = wks_to_add)
    id <- rep(this_id, wks_to_add)

    if(tb == "yes"){
      if(tpt == "yes"){
        wt_tpt_tb <- rep(0, wks_to_add)
        data.table(id = id, wk = wk, wt_tpt_tb = wt_tpt_tb)
      } else {
        wt_cntrl_tb <- rep(0, wks_to_add)
        data.table(id = id, wk = wk, wt_cntrl_tb = wt_cntrl_tb)
      }
    } else {
      if(tpt == "yes"){
        wt_tpt_death <- rep(0, wks_to_add)
        data.table(id = id, wk = wk, wt_tpt_death = wt_tpt_death)
      } else {
        wt_cntrl_death <- rep(0, wks_to_add)
        data.table(id = id, wk = wk, wt_cntrl_death = wt_cntrl_death)
      }
    }
  }))
  rm(wk_summary)
  
  # Combine with original data and compute truncation limits
  if (tpt == "yes") {
    if (tb == "yes") {
      weekly_wts <- weekly_records_data[, .(id, wk, wt_tpt_tb)]
    } else {
      weekly_wts <- weekly_records_data[, .(id, wk, wt_tpt_death)]
    }
  } else {
    if (tb == "yes") {
      weekly_wts <- weekly_records_data[, .(id, wk, wt_cntrl_tb)]
    } else {
      weekly_wts <- weekly_records_data[, .(id, wk, wt_cntrl_death)]
    }
  }
  
  # Combine original data with new rows
  weekly_wts <- rbindlist(list(weekly_wts, rows_add))
  
  # Calculate 99th percentile
  weekly_wts <- weekly_wts[, .(wt_trunc = quantile(get(names(weekly_wts)[3]), 0.99, na.rm = TRUE)), 
                           by = wk]
  # if(tpt == "yes"){
  #   if(tb == "yes"){
  #     weekly_wts[, wt_tpt_tb := NULL]
  #   } else {
  #     weekly_wts[, wt_tpt_death := NULL]
  #   }
  # } else {
  #   if(tb == "yes"){
  #     weekly_wts[, wt_cntrl_tb := NULL]
  #   } else {
  #     weekly_wts[, wt_cntrl_death := NULL]
  #   }
  # }
  
  # Merge with original data
  weekly_records_data <- merge(weekly_records_data, weekly_wts, by = "wk", all.x = TRUE)
  
  # Truncate based on the quantiles
  if (tpt == "yes") {
    if (tb == "yes") {
      weekly_records_data[, wt_tpt_tb := pmin(wt_tpt_tb, wt_trunc)]
    } else {
      weekly_records_data[, wt_tpt_death := pmin(wt_tpt_death, wt_trunc)]
    }
  } else {
    if (tb == "yes") {
      weekly_records_data[, wt_cntrl_tb := pmin(wt_cntrl_tb, wt_trunc)]
    } else {
      weekly_records_data[, wt_cntrl_death := pmin(wt_cntrl_death, wt_trunc)]
    }
  }
  
  # Remove unnecessary column
  weekly_records_data[, wt_trunc := NULL]

  return(weekly_records_data)
}