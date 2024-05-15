#' Function to create a data set with weekly records
#' 
#' @param dat a \code{data.table}. See function for variable name
#'   requirements
#' @param baseline_covariates a \code{string} of column names for 
#'   baseline covariates
#' @param k a \code{numeric} indicating the number of past encounters
#'   to summarize
#' @param verbose not currently in use
#' @param admin_cens_wks a \code{numeric} indicating the number of 
#'   weeks prior to administrative censoring
#' @param ... not currently in use
#' @return a \code{data.table}

create_weekly_record_data <- function(
  dat,
  baseline_covariates,
  k = 3,
  verbose = FALSE, 
  admin_cens_wks = 52 * 2, # two years
  ...
){

  assertthat::assert_that("id" %in% colnames(dat))
  assertthat::assert_that("visit_date" %in% colnames(dat))
  assertthat::assert_that("right_cens_date_tb" %in% colnames(dat))
  assertthat::assert_that("last_visit_date" %in% colnames(dat))

  assertthat::assert_that("enroll_date" %in% colnames(dat))
  assertthat::assert_that("art_start_date" %in% colnames(dat))
  assertthat::assert_that("tpt_start_date" %in% colnames(dat))
  assertthat::assert_that("tb_diagnosis_date" %in% colnames(dat))
  assertthat::assert_that("death_date" %in% colnames(dat))
  
  assertthat::assert_that("art_adherence" %in% colnames(dat))
  # because we fill in with 0's later, make sure adherence is coded as numeric
  assertthat::assert_that(class(dat[, art_adherence]) == "numeric")

  assertthat::assert_that("tb_sx" %in% colnames(dat))
  # because we fill in with 0's later, make sure symptoms are coded as numeric
  assertthat::assert_that(class(dat[, tb_sx]) == "numeric")

  assertthat::assert_that(all(baseline_covariates %in% colnames(dat)))

  n_id <- length(unique(dat[,id]))

  obs_by_week <- future.apply::future_lapply(dat[, unique(id)], function(this_id, k) {
    
    dat_id <- dat[id == this_id]
    dat_id[, admin_cens_date := enroll_date + admin_cens_wks * 7]

    fup_start <- as_date(dat_id$enroll_date[1])
    
    if(!is.na(dat_id$death_date)){
      fup_end <- min(
        dat_id$admin_cens_date, dat_id$death_date
      )
    }else{
      fup_end <- min(
        dat_id$admin_cens_date, dat_id$last_visit_date
      )
    }

    visit_weeks_all <- seq(fup_start, fup_end + 7, by = "weeks")
    n_weeks <- length(visit_weeks_all)

    weekly_records <- data.table(
      wk = 1:(length(visit_weeks_all) - 1),
      week_start = visit_weeks_all[-length(visit_weeks_all)],
      week_end = visit_weeks_all[-1]
    )

    uniq_art_init_date <- unique(dat_id$art_start_date)
    assertthat::assert_that(length(uniq_art_init_date) == 1)
    art_init_date <- as_date(uniq_art_init_date)
    if (!is.na(art_init_date)) {
      weekly_records[,
        art_init := art_init_date %within% interval(fup_start, week_end)
      ]
    }else{
      weekly_records[, art_init := FALSE]
    }

    uniq_tpt_start_date <- unique(dat_id$tpt_start_date)
    assertthat::assert_that(length(uniq_tpt_start_date) == 1)
    tpt_init_date <- as_date(uniq_tpt_start_date)
    if (!is.na(tpt_init_date)) {
      weekly_records[,
        tpt_init := tpt_init_date %within% interval(fup_start, week_end)
      ]
    } else {
      weekly_records[, tpt_init := FALSE]
    }

    cd4_count_this_id <- unique(dat_id$cd4_count)
    cd4_count_date_this_id <- unique(dat_id$cd4_count_date)
    assertthat::assert_that(length(cd4_count_this_id) == 1)

    # find k closest clinic visits to week j 
    # only include visits in week j if tpt_init_date not in 
    # interval(week_start, week_end)
    for (j in seq_len(n_weeks-1)) {

      if(!is.na(tpt_init_date)){
        tpt_init_wk_j <- tpt_init_date %within% interval(weekly_records[j, week_start], weekly_records[j, week_end])
      }else{
        tpt_init_wk_j <- FALSE
      }
      
      if(tpt_init_wk_j){
        visits_before_date <- tpt_init_date
      }else{
        visits_before_date <- weekly_records$week_start[j]
      }

      all_prior_visit_dates <- dat_id[visit_date < visits_before_date, visit_date]
      if(length(all_prior_visit_dates) > 0){
        k_most_recent_visit_dates <- dat_id$visit_date[order(visits_before_date - all_prior_visit_dates)[1:k]]
        k_most_recent_visits <- dat_id[visit_date %in% k_most_recent_visit_dates]
        n_visits <- nrow(k_most_recent_visits) # n_visits <= k
      }else{
        k_most_recent_visits <- NULL
        n_visits <- 0
      }
      
      cd4_count_variables <- summarize_cd4_count(
        cd4_count_date = cd4_count_this_id, 
        cd4_count = cd4_count_this_id,
        visits_before_date = visits_before_date                    
      )

      # variables indicating whether a visit is present
      have_visit_variables <- "have_visit_k"
      if(k > 1){
        have_visit_variables <- c(have_visit_variables, paste0("have_visit_kminus", 1:(k-1)))
      }
      visit_counter <- 0
      for(have_visit_variable in have_visit_variables){
        visit_counter <- visit_counter + 1
        weekly_records[j, (have_visit_variable) := as.numeric(n_visits >= visit_counter)]
      }

      # variables summarizing ART adherence
      art_adherence_variables <- summarize_art_adherence(k_most_recent_visits)
      for(art_adherence_variable in names(art_adherence_variables)){
        weekly_records[j, (art_adherence_variable) := art_adherence_variables[[art_adherence_variable]]]
      }
      
      # variables summarizing tb symptom screen results adherence
      tb_symptoms_variables <- summarize_tb_symptoms(k_most_recent_visits)
      for(tb_symptoms_variable in names(tb_symptoms_variables)){
        weekly_records[j, (tb_symptoms_variable) := tb_symptoms_variables[[tb_symptoms_variable]]]
      }

      visit_times_variables <- summarize_visit_times(
        all_prior_visit_dates = all_prior_visit_dates,
        visits_before_date = visits_before_date,
        k_most_recent_visits = k_most_recent_visits
      )
      for(visit_times_variable in names(visit_times_variables)){
        weekly_records[j, (visit_times_variable) := visit_times_variables[[visit_times_variable]]]
      }
    }


    enroll_date_this_id <- dat_id[1, enroll_date]

    tpt_start_date_this_id <- dat_id[1, tpt_start_date]
    tpt_start_wk_this_id <- get_week_of_event(tpt_start_date_this_id, enroll_date_this_id)
    weekly_records[, tpt_start_wk := tpt_start_wk_this_id]

    tb_diagnosis_date_this_id <- dat_id[1, tb_diagnosis_date]
    tb_wk_this_id <- get_week_of_event(tb_diagnosis_date_this_id, enroll_date_this_id)
    weekly_records[, tb_wk := tb_wk_this_id]

    death_date_this_id <- dat_id[1, death_date]
    death_wk_this_id <- get_week_of_event(death_date_this_id, enroll_date_this_id)
    weekly_records[, death_wk := death_wk_this_id]

    right_cens_date_this_id <- dat_id[1, right_cens_date_tb]
    right_cens_wk_this_id <- get_week_of_event(right_cens_date_this_id, enroll_date_this_id)
    weekly_records[, right_cens_wk_tb := right_cens_wk_this_id]
    
    last_visit_date_this_id <- dat_id[1, last_visit_date]
    last_visit_wk_this_id <- get_week_of_event(last_visit_date_this_id, enroll_date_this_id)
    weekly_records[, last_visit_wk := last_visit_wk_this_id]

    admin_cens_date_this_id <- dat_id[1, admin_cens_date]
    admin_cens_wk_this_id <- get_week_of_event(admin_cens_date_this_id, enroll_date_this_id)
    weekly_records[, admin_cens_wk := admin_cens_wk_this_id]
    
    # output DT of weekly records for given study unit
    weekly_records[, id := this_id]
    setcolorder(weekly_records, "id")
    return(weekly_records)
  }, k = k) # end lapply

  weekly_records_allids <- rbindlist(obs_by_week)

  # magic gpt code to get baseline only data
  baseline_dat <- dat[, lapply(.SD, first), by = .(id), .SDcols = baseline_covariates]

  weekly_data <- weekly_records_allids[baseline_dat, on = .(id)]

  return(weekly_data)
}

#' @param k_most_recent_visits A \code{data.table} of the most recent k encounters
#' with relevant ART adherence included. If \code{k_most_recent_visits} is \code{NULL}
#' then it is assumed that there were no encounters prior to a particular week.

summarize_art_adherence <- function(
  k_most_recent_visits
){
  out <- list(
    art_adherence_measured_visit_k = 0, 
    art_adherence_measured_visit_kminus1 = 0,
    art_adherence_measured_visit_kminus2 = 0,
    art_adherence_visit_k = 0,
    art_adherence_visit_kminus1 = 0,
    art_adherence_visit_kminus2 = 0
  )
  if(!is.null(k_most_recent_visits)){
    n_visits <- nrow(k_most_recent_visits)
    out$art_adherence_measured_visit_k <- as.numeric(!is.na(
      k_most_recent_visits[visit_date == max(visit_date), art_adherence]
    ))
    if(n_visits > 1){
      out$art_adherence_measured_visit_kminus1 <- as.numeric(!is.na(
        k_most_recent_visits[visit_date == visit_date[order(visit_date)[n_visits - 1]], art_adherence]
      ))
    }else{
      out$art_adherence_measured_visit_kminus1 <- 0
    }
    if(n_visits > 2){
      out$art_adherence_measured_visit_kminus2 <- as.numeric(!is.na(
        k_most_recent_visits[visit_date == visit_date[order(visit_date)[n_visits - 2]], art_adherence]
      ))
    }else{
      out$art_adherence_measured_visit_kminus2 <- 0
    }
  }

  if(!is.null(k_most_recent_visits)){
    n_visits <- nrow(k_most_recent_visits)
    out$art_adherence_visit_k <- zero_if_missing(
      k_most_recent_visits[visit_date == max(visit_date), art_adherence]
    )
    if(n_visits > 1){
      out$art_adherence_visit_kminus1 <- zero_if_missing(
        k_most_recent_visits[visit_date == visit_date[order(visit_date)[n_visits - 1]], art_adherence]
      )
    }else{
      out$art_adherence_visit_kminus1 <- 0
    }
    if(n_visits > 2){
      out$art_adherence_visit_kminus2 <- zero_if_missing(
        k_most_recent_visits[visit_date == visit_date[order(visit_date)[n_visits - 2]], art_adherence]
      )
    }else{
      out$art_adherence_visit_kminus2 <- 0
    }
  }
  return(out)
}

summarize_cd4_count <- function(
  cd4_count_date,
  cd4_count,
  visits_before_date                                
){
  out <- list(
    have_cd4 = 0,
    cd4_count = 0,
    days_since_cd4_count = 0                  
  )
  if((cd4_count_date < visits_before_date) & 
     !(is.na(cd4_count))){
    out$have_cd4 <- 1
    out$cd4_cound <- cd4_count
    out$days_since_cd4_count <- round(visits_before_date - cd4_count_date)
  }
  return(out)
}

summarize_visit_times <- function(
  all_prior_visit_dates,
  visits_before_date,
  k_most_recent_visits
){
  out <- list(
    number_encounters_last_2_weeks = 0,
    number_encounters_last_4_weeks = 0,
    days_between_week_date_and_visit_k = 0,
    days_between_week_date_and_visit_kminus1 = 0,
    days_between_week_date_and_visit_kminus2 = 0
  )
  if(length(all_prior_visit_dates) > 0){
    out$number_encounters_last_2_weeks <- sum(all_prior_visit_dates > (visits_before_date - 2 * 7))
    out$number_encounters_last_4_weeks <- sum(all_prior_visit_dates > (visits_before_date - 4 * 7))
  }

  if(!is.null(k_most_recent_visits)){
    n_visits <- nrow(k_most_recent_visits)
    out$days_between_week_date_and_visit_k <- round(as.numeric(
      visits_before_date - max(k_most_recent_visits[,visit_date])
    ))
    if(n_visits > 1){
      out$days_between_week_date_and_visit_kminus1 <- round(as.numeric(
        visits_before_date - k_most_recent_visits[visit_date == visit_date[order(visit_date)[n_visits - 1]], visit_date]
      ))
    }else{
      out$days_between_week_date_and_visit_kminus1 <- 0
    }
    if(n_visits > 2){
      out$days_between_week_date_and_visit_kminus2 <- round(as.numeric(
        visits_before_date - k_most_recent_visits[visit_date == visit_date[order(visit_date)[n_visits - 2]], visit_date]
      ))
    }else{
      out$days_between_week_date_and_visit_kminus2 <- 0
    }
  }
  return(out)
}
    


#' @param k_most_recent_visits A \code{data.table} of the most recent k encounters
#' with relevant ART adherence included. If \code{k_most_recent_visits} is \code{NULL}
#' then it is assumed that there were no encounters prior to a particular week.

summarize_tb_symptoms <- function(
  k_most_recent_visits
){
  out <- list(
    tb_symptoms_measured_visit_k = 0, 
    tb_symptoms_measured_visit_kminus1 = 0,
    tb_symptoms_measured_visit_kminus2 = 0,
    tb_symptoms_visit_k = 0,
    tb_symptoms_visit_kminus1 = 0,
    tb_symptoms_visit_kminus2 = 0
  )
  if(!is.null(k_most_recent_visits)){
    n_visits <- nrow(k_most_recent_visits)
    out$tb_symptoms_measured_visit_k <- as.numeric(!is.na(
      k_most_recent_visits[visit_date == max(visit_date), tb_sx]
    ))
    if(n_visits > 1){
      out$tb_symptoms_measured_visit_kminus1 <- as.numeric(!is.na(
        k_most_recent_visits[visit_date == visit_date[order(visit_date)][2], tb_sx]
      ))
    }else{
      out$tb_symptoms_measured_visit_kminus1 <- 0
    }
    if(n_visits > 2){
      out$tb_symptoms_measured_visit_kminus2 <- as.numeric(!is.na(
        k_most_recent_visits[visit_date == visit_date[order(visit_date)][3], tb_sx]
      ))
    }else{
      out$tb_symptoms_measured_visit_kminus2 <- 0
    }
  }

  if(!is.null(k_most_recent_visits)){
    n_visits <- nrow(k_most_recent_visits)
    out$tb_symptoms_visit_k <- zero_if_missing(
      k_most_recent_visits[visit_date == max(visit_date), tb_sx]
    )
    if(n_visits > 1){
      out$tb_symptoms_visit_kminus1 <- zero_if_missing(
        k_most_recent_visits[visit_date == visit_date[order(visit_date)[n_visits - 1]], tb_sx]
      )
    }else{
      out$tb_symptoms_visit_kminus1 <- 0
    }
    if(n_visits > 2){
      out$tb_symptoms_visit_kminus2 <- zero_if_missing(
        k_most_recent_visits[visit_date == visit_date[order(visit_date)[n_visits - 2]], tb_sx]
      )
    }else{
      out$tb_symptoms_visit_kminus2 <- 0
    }
  }
  return(out)
}

#' Return value if observed and replace with 0 otherwise
#' 
#' This function is used while summarizing encounter-level data
#' 
#' @param x A \code{numeric}
#' @return A \code{numeric} with NAs replaced with zeros

zero_if_missing <- function(x){
  return(ifelse(is.na(x), 0, x))
}

#' Convert event dates to weeks and replace missing values with
#' an arbitrary numeric
#' 
#' @param event_date a \code{Date} object
#' @param enroll_date a \code{Date} object
#' @param missing_value a \code{numeric} value that is used to fill
#'   in missing values
get_week_of_event <- function(event_date, enroll_date, missing_value = 99999){
  if(!is.na(event_date)){
    event_wk <- 1 + floor((event_date - enroll_date) / 7)
    class(event_wk) <- "numeric"
  }else{
    event_wk <- missing_value
  }
  return(event_wk)
}