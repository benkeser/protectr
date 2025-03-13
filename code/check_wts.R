#' Function to add 0's back onto weekly records data
#' 
#' @param  cloned_data_sets \code{list}. Cloned data set object
#'   created by \code{create_cloned_data_set} function
check_wts <- function(cloned_data_sets){
	tb_z1 <- add_zeroes(cloned_data_sets = cloned_data_sets, condition = "tb_z1")
	tb_z0 <- add_zeroes(cloned_data_sets = cloned_data_sets, condition = "tb_z0")
	death_z1 <- add_zeroes(cloned_data_sets = cloned_data_sets, condition = "death_z1")
	death_z0 <- add_zeroes(cloned_data_sets = cloned_data_sets, condition = "death_z0")

	tb_z1_mean <- tb_z1 %>% group_by(wk) %>% summarize(mean_wt_k = mean(wt_k, na.rm = TRUE))
	tb_z0_mean <- tb_z0 %>% group_by(wk) %>% summarize(mean_wt_k = mean(wt_k, na.rm = TRUE))
	death_z1_mean <- death_z1 %>% group_by(wk) %>% summarize(mean_wt_k = mean(wt_k, na.rm = TRUE))
	death_z0_mean <- death_z0 %>% group_by(wk) %>% summarize(mean_wt_k = mean(wt_k, na.rm = TRUE))

	out <- list()
	out$weights <- list(
		tb_z1 <- tb_z1,
		tb_z0 <- tb_z0,
		death_z1 <- death_z1,
		death_z0 <- death_z0
    )
	out$weekly_means <- list(
		tb_z1_mean <- tb_z1_mean,
		tb_z0_mean <- tb_z0_mean,
		death_z1_mean <- death_z1_mean,
		death_z0_mean <- death_z0_mean
    )

	return(out)
}

#' Function to add 0's back onto weekly records data
add_zeroes <- function(cloned_data_sets, condition){
  
  # Filter based on condition and z value
  if(condition == "tb_z1"){
    tb_data <- as.data.table(cloned_data_sets$tb)
    data <- tb_data[z == 1, .(id, wk, tb_wk, death_wk, wt_k)]
    rm(tb_data)
  } else if(condition == "tb_z0"){
    tb_data <- as.data.table(cloned_data_sets$tb)
    data <- tb_data[z == 0, .(id, wk, tb_wk, death_wk, wt_k)]
    rm(tb_data)
  } else if(condition == "death_z1"){
    death_data <- as.data.table(cloned_data_sets$death)
    data <- death_data[z == 1, .(id, wk, tb_wk, death_wk, wt_k)]
    rm(death_data)
  } else {
    death_data <- as.data.table(cloned_data_sets$death)
    data <- death_data[z == 0, .(id, wk, tb_wk, death_wk, wt_k)]
    rm(death_data)
  }

  # Create the summary data.table
  if (condition %in% c("tb_z1", "tb_z0")) {
    data_summary <- data[, .(wks = .N, cens_wk = min(c(first(tb_wk), first(death_wk), 103))), by = id]
  } else {
    data_summary <- data[, .(wks = .N, cens_wk = min(c(first(death_wk), 103))), by = id]
  }
  
  # Compute the number of weeks to add
  data_summary[, wks_to_add := cens_wk - wks]
  
  # Filter rows where wks_to_add > 0
  data_summary <- data_summary[wks_to_add > 0]
  
  # Create the additional data.table to add
  data_add <- data_summary[, {
    # For each id, replicate values and create new rows
    .(id = rep(id, wks_to_add),
      wk = seq(from = wks + 1, by = 1, length.out = wks_to_add),
      tb_wk = NA_real_,
      death_wk = NA_real_,
      wt_k = 0)
  }, by = id]
  data_add[,2] <- NULL
  
  # Combine original data with the new data
  data_all <- rbindlist(list(data, data_add))
  
  # Order by id and wk
  setorder(data_all, id, wk)
  
  return(data_all)
}