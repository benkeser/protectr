#' Function to add 0's back onto weekly records data
#' 
#' @param  cloned_data_sets \code{list}. Cloned data set object
#'   created by \code{create_cloned_data_set} function
#' @param condition a \code{character} in one of "tb_z1", "tb_z0",
#'   "death_z1", or "death_z0" to indicate which weights are desired

add_zeroes <- function(cloned_data_sets, condition){
  if(condition == "tb_z1"){
    data <- cloned_data_sets$tb %>% filter(z == 1) %>% select(id, wk, tb_wk, death_wk, prob_wt_cens_tb)
  } else if(condition == "tb_z0"){
    data <- cloned_data_sets$tb %>% filter(z == 0) %>% select(id, wk, tb_wk, death_wk, wt_k)
  } else if(condition == "death_z1"){
    data <- cloned_data_sets$death %>% filter(z == 1) %>% select(id, wk, tb_wk, death_wk, prob_wt_cens_death)
  } else {
    data <- cloned_data_sets$death %>% filter(z == 0) %>% select(id, wk, tb_wk, death_wk, wt_k)
  }
  data_summary <- data %>% group_by(id) %>%
    summarize(wks = length(id), cens_wk = min(first(tb_wk), first(death_wk), 103)) %>%
    mutate(wks_to_add = cens_wk - wks) %>%
    filter(wks_to_add > 0)
  iterations <- nrow(data_summary)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  if(condition == "tb_z1"){
    data_add <- foreach(i = 1:nrow(data_summary), .combine = rbind, .options.snow = opts) %dopar% {
      this_id <- data_summary$id[i]
      id <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], this_id)
      wk <- seq(from = data_summary$wks[data_summary$id %in% this_id] + 1, by = 1, length.out = data_summary$wks_to_add[data_summary$id %in% this_id])
      tb_wk <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], NA)
      death_wk <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], NA)
      prob_wt_cens_tb <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], 0)
      data.frame(id, wk, tb_wk, death_wk, prob_wt_cens_tb)
    }
  } else if(condition == "tb_z0"){
    data_add <- foreach(i = 1:nrow(data_summary), .combine = rbind, .options.snow = opts) %dopar% {
      this_id <- data_summary$id[i]
      id <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], this_id)
      wk <- seq(from = data_summary$wks[data_summary$id %in% this_id] + 1, by = 1, length.out = data_summary$wks_to_add[data_summary$id %in% this_id])
      tb_wk <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], NA)
      death_wk <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], NA)
      wt_k <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], 0)
      data.frame(id, wk, tb_wk, death_wk, wt_k)
    }
  } else if(condition == "death_z1"){
    data_add <- foreach(i = 1:nrow(data_summary), .combine = rbind, .options.snow = opts) %dopar% {
      this_id <- data_summary$id[i]
      id <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], this_id)
      wk <- seq(from = data_summary$wks[data_summary$id %in% this_id] + 1, by = 1, length.out = data_summary$wks_to_add[data_summary$id %in% this_id])
      tb_wk <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], NA)
      death_wk <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], NA)
      prob_wt_cens_death <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], 0)
      data.frame(id, wk, tb_wk, death_wk, prob_wt_cens_death)
    }
  } else {
    data_add <- foreach(i = 1:nrow(data_summary), .combine = rbind, .options.snow = opts) %dopar% {
      this_id <- data_summary$id[i]
      id <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], this_id)
      wk <- seq(from = data_summary$wks[data_summary$id %in% this_id] + 1, by = 1, length.out = data_summary$wks_to_add[data_summary$id %in% this_id])
      tb_wk <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], NA)
      death_wk <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], NA)
      wt_k <- replicate(data_summary$wks_to_add[data_summary$id %in% this_id], 0)
      data.frame(id, wk, tb_wk, death_wk, wt_k)
    }
  }
  data_all <- rbind(data, data_add)
  data_all <- data_all[order(data_all$id, data_all$wk),]
  
  return(data_all)
}