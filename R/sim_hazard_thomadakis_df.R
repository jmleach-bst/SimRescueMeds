#' Simulate Time-to-Rescue-Therapy for a Data Set
#'
#' Use proportional hazards models to simulate time-to-rescue therapy
#' that is dependent on longitudinal trajectories.
#'
#' @inheritParams sim_hazard_thomadakis
#'
#' @param data A data frame of longitudinal measurements in long format, i.e.
#' 1 row per observation, not 1 row per subject.
#' @param id_var Character. The name of the subject identifier variable.
#' @param y_var Character. The name of the longitudinal outcome measurement variable.
#' @param time_var Character. The name of the time variable.
#' @param event_factor Logical. When \code{TRUE} creates a factor variable for the event
#' indicator, which may be beneficial if no events occur.
#' @param print_intermediate Logical. When \code{TRUE} prints intermediate results.
#'
#' @importFrom utils tail
#'
#' @return A data frame.
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
#' sim_hazard_thomadakis_df(
#'   data = head(df, 15),
#'   theta1 = 1/3,
#'   theta2 = 0,
#'   lambda = 0.05,
#'   print_intermediate = TRUE
#' )
#'
#' sim_hazard_thomadakis_df(
#'   data = head(df, 15),
#'   theta1 = 0,
#'   theta2 = 1/3,
#'   lambda = 0.05,
#'   print_intermediate = TRUE,
#'   dist = "weibull"
#' )
#'
#'
#' @details
#' This function simulates time-to-rescue medication based on the approach of Thomadakis (2019).
#' A general formulation for simulating from missing-at-random (MAR) and missing-not-at-random
#' (MNAR) data is to allow the hazard of rescue therapy within some interval \eqn{t_{i,j} \le t < t_{i,j+1}}
#' to be a function of the last observed outcome measurement, \eqn{y_{i,j}}, and the next measurement,
#' \eqn{y_{i,j+1}}, the latter of which is unobserved in the study since it occurs after rescue therapy.
#' If the event-time distribution is exponential, then we have the following hazard function, but it is
#' straightforward to extend to Weibull, for which exponential is a special case.
#'
#' \deqn{
#' h_i(t) = \lambda\exp\left(y_{i,j}\theta_1 + y_{i,j+1}\theta_2\right), \quad t_{i,j} \le t < t_{i,j+1}.
#' }
#'
#' This requires simulating times-to-rescue therapy within each interval between measurements until either
#' one of those intervals produces rescue therapy or the study ends. If \eqn{\theta_2 \ne 0}, then there is
#' MNAR missingness. If \eqn{\theta_1 = \theta_2 = 0}, then MCAR. If \eqn{\theta_1 \ne 0} and \eqn{\theta_2 = 0},
#' then MAR.
#'
#' @references
#' \insertRef{Thomadakis2019}{SimRescueMeds}
#'
#' @export
sim_hazard_thomadakis_df <- function(
    data,
    id_var = "id",
    y_var = "y",
    time_var = "time",
    theta1 = 0,
    theta2 = 0,
    dist = "exponential",
    lambda = 1,
    gamma = 1,
    maxt = 1,
    print_intermediate = FALSE,
    event_factor = TRUE # create factor variable for event indicator (helps if no events occur)
) {
  unique_ids <- unique(data[[id_var]])
  n <- length(unique_ids)
  tte_list <- vector(mode = "list", length = n)
  for (j in seq_len(n)) {
    id_j <- unique_ids[j]
    data_j <- data[data[[id_var]] == id_j,]
    tte_j <- sim_hazard_thomadakis(
      y = data_j[[y_var]],
      time = data_j[[time_var]],
      theta1 = theta1,
      theta2 = theta2,
      dist = dist,
      lambda = lambda,
      gamma = gamma
    )
    tte_list[[j]] <- cbind(id = id_j, tail(tte_j, n = 1))
    if (print_intermediate == TRUE) {
      cat("Here are all draws for individual", j, "\n")
      print(tte_j)
      cat("Here is final time and event indicator for individual ", j, "\n")
      print(tte_list[[j]])
    }
  }
  tte_df <- as.data.frame(
    do.call('rbind', tte_list)
  )
  rownames(tte_df) <- NULL
  if (event_factor == TRUE) {
    tte_df$event_factor <- factor(
      tte_df$event,
      levels = c(0, 1)
    )
  }
  return(tte_df)
}
