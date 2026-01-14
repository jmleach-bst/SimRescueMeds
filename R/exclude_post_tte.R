#' Exclude Longitudinal Measurements Occuring Post-Rescue Therapy
#'
#' Excludes longitudinal measurements that occur after rescue therapy (time-to-event)
#' and merges the longitudinal and time-to-event data.
#'
#' @param original_data A data frame containing longitudinal measurements in long format.
#' @param tte_data A data frame containing time-to-event measurements with 1 row per subject.
#' @param id_var Character. The name of the subject identifier variable. The variable name
#' should be the same for both \code{original_data} and \code{tte_data} data frames.
#' @param y_var Character. The name of the longitudinal outcome measurement variable.
#' @param time_var Character. The name of the time variable in the longitudinal data.
#' @param tte_var Character. The variable name for the column containing the time-to-event measurements.
#'
#' @details
#' For the purposes of analysis with linear mixed models or joint longitudinal and
#' time-to-event models, longitudinal measurements that occur after the event, e.g.,
#' rescue therapy, are excluded from the analysis. This function removes the rows of
#' longitudinal data that occur post-event, merges the longitudinal and time-to-event
#' data, and produces 2 new variables, \code{last_time_preRT} and \code{last_y_preRT},
#' that indicate the time the last pre-event longitudinal measurement occurred and the
#' measured value at that time, respectively.
#'
#'
#' @examples
#' set.seed(3002242)
#'
#' xmat_t <- build_design_matrix(
#'     N_t = 1,
#'     N_c = 0,
#'     K = 5,
#'     time_start = 0,
#'     time_scale = 1/4
#'   )
#'
#' df <- simulate_lmm_rct(
#'   N_t = 50,
#'   N_c = 50,
#'   xmat_t = xmat_t,
#'   xmat_c = build_design_matrix(
#'     N_t = 0,
#'     N_c = 1,
#'     K = 5,
#'     time_start = 0,
#'     time_scale = 1/4
#'   ),
#'   betas = c(
#'     beta0 = 8,
#'     beta1 = 0,
#'     beta2 = -1,
#'     beta3 = -1
#'   ),
#'   Sigma = build_lmm_cov(
#'     zmat = xmat_t[, c("intercept", "time")],
#'     re_sigma = c(1.15, 1.05),
#'     re_corr_mat = matrix(
#'       c(1, 0,
#'         0, 1),
#'       byrow = TRUE,
#'       nrow = 2
#'       ),
#'     ws_sigma = 1.25,
#'     ws_corr_mat = diag(5),
#'     print_intermediate = FALSE
#'     )
#' )
#'
#' set.seed(707582)
#'
#' tte_ex  <- sim_hazard_thomadakis_df(
#'   data = df,
#'   theta1 = 1/3,
#'   theta2 = 0,
#'   lambda = 0.05,
#'   print_intermediate = FALSE
#' )
#'
#' post_tte_ex <- exclude_post_tte(
#'   original_data = df,
#'   tte_data = tte_ex
#' )
#'
#' head(post_tte_ex, 25)
#'
#'
#' @export
exclude_post_tte <- function(
    original_data,
    tte_data,
    id_var = "id",
    y_var = "y",
    time_var = "time",
    tte_var = "eventtime"
) {
  merged_data <- merge(
    x = original_data,
    y = tte_data,
    all.x = TRUE,
    by = id_var
  )
  merged_data <- merged_data[order(merged_data[[id_var]]), ]
  merged_data_ex <- merged_data[
    merged_data[[time_var]] <= merged_data[[tte_var]],
  ]
  last_measurement <- do.call(
    'rbind',
    by(data = merged_data_ex,
       INDICES = merged_data_ex[[id_var]],
       FUN = function(x) x[which.max(x[[time_var]]), ])
  )[, c(id_var, time_var, y_var)]
  colnames(last_measurement) <- c(id_var, "last_time_preRT", "last_y_preRT")
  merged_data_ex <- merge(
    x = merged_data_ex,
    y = last_measurement,
    all.x = TRUE,
    by = id_var
  )
  merged_data_ex <- merged_data_ex[order(merged_data_ex[[id_var]]), ]
  rownames(merged_data_ex) <- NULL
  return(merged_data_ex)
}
