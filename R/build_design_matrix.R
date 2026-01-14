#' Design Matrices for Longitudinal RCT
#'
#' Builds a design matrix for a 2-group RCT assuming longitudinal follow-up.
#'
#' @param N_t Numeric. Number of patients that receive active treatment.
#' @param N_c Numeric. Number of patients that receive control treatment.
#' @param K Numeric. The number of longitudinal measurements.
#' @param time_start Numeric. Defines the starting time of follow-up.
#' @param time_scale Numeric. A value multiplied by the time vector to scale times of measurement.
#' @param interaction Logical. When \code{TRUE} includes an interaction term between treatment and time.
#' @param num_binary_covariates Numeric. The number of binary covariates to draw from Bernoulli distribution.
#' @param success_probs Numeric vector of success probabilities for binary covariates.
#' @param num_continuous_covariates Numeric. The number of continuous covariates to draw from a standard
#' Normal distribution.
#' @param return_type Character. Incicates whether to return a \code{"matrix"} or \code{"data.frame"}.
#'
#' @importFrom stats rbinom
#' @importFrom stats rnorm
#'
#' @return A matrix or data frame.
#'
#' @note
#' This function produces a design matrix that can be used to generate longitudinal outcome data. That is,
#' it is a first step in simulating longitudinal RCT data. If no covariates are used then it may be
#' produced a single time and does not need to be recreated anew for each data set. If covariates are included
#' the intention of the simulation comes into play, i.e., generating data anew for each data set implies that
#' the simulations sample new individuals, whereas generating once assumes we generate new outcomes for the
#' same individuals.
#'
#' @examples
#' # fixed effects design for treated
#' xmat_trt <- build_design_matrix(
#'   N_t = 1,
#'   N_c = 0,
#'   K = 5,
#'   time_start = 0,
#'   time_scale = 1/4
#' )
#'
#' # fixed effects design for control
#' xmat_ctrl <- build_design_matrix(
#'   N_t = 0,
#'   N_c = 1,
#'   K = 5,
#'   time_start = 0,
#'   time_scale = 1/4
#' )
#'
#'
#' @export
build_design_matrix <- function(
    N_t,
    N_c,
    K,
    time_start = 0,
    time_scale = 1,
    interaction = TRUE,
    num_binary_covariates = 0,
    success_probs = 0.5,
    num_continuous_covariates = 0,
    return_type = "matrix"
){
  N <- N_t + N_c
  subj_id <- rep(1:N, each = K)
  x <- rep(
    c(rep(1, N_t),
      rep(0, N_c)),
    each = K
  )
  time <- rep(time_start:(time_start + (K - 1)), N)
  if (time_scale != 1) {
    time <- time*time_scale
  }
  if (num_binary_covariates < 0 | num_continuous_covariates < 0) {
    warning(
      "Require # of covariates >= 0. No covariates added. \n")
  }
  if (interaction == TRUE) {
    xmat <- cbind(
      intercept = 1, x = x, time = time, x_by_time = x*time
    )
  } else {
    xmat <- cbind(
      intercept = 1, x = x, time = time
    )
  }
  if (num_binary_covariates > 0) {
    if (!(length(success_probs) %in% c(1, num_binary_covariates))) {
      stop("Require # of success probs == 1 or == num_binary_covariates.")
    }
    binary_list <- list()
    if (num_binary_covariates > 1 & length(success_probs) == 1) {
      success_probs <- rep(success_probs, num_binary_covariates)
    }
    for (i in 1:num_binary_covariates) {
      binary_list[[i]] <- rep(
        rbinom(n = N_t + N_c, size = 1, prob = success_probs[i]), each = K
      )
    }
    binary_mat <- do.call(cbind, binary_list)
    colnames(binary_mat) <- paste0("BC_", 1:num_binary_covariates)
    xmat <- cbind(xmat, binary_mat)
  }
  if (num_continuous_covariates > 0) {
    continuous_list <- list()
    for (i in 1:num_continuous_covariates) {
      continuous_list[[i]] <- rep(
        rnorm(n = N_t + N_c), each = K,
      )
    }
    continuous_mat <- do.call(cbind, continuous_list)
    colnames(continuous_mat) <- paste0("CC_", 1:num_continuous_covariates)
    xmat <- cbind(xmat, continuous_mat)
  }
  if (return_type == "matrix") {
    return(xmat)
  } else {
    return(
      as.data.frame(
        cbind(subj_id = rep(1:(N_t + N_c), each = K), xmat)
      )
    )
  }
}
