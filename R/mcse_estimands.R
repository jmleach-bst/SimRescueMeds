#' Calculate MCSE Relevant for Estimands
#'
#' Combines estimation of MCSE for bias, EmpSE, ModSE, and MSE into a single function.
#'
#' @param data A data frame where each row contains results from a unique simulated data set.
#' @param true_value Numeric. The true value of the estimand.
#' @param estimand_name Character. The variable name for the column containing the estimates
#' of the estimand from each simulated data set.
#' @param se_name Character. The variable name for the column containing the estimates
#' of the SE of the estimator from each simulated data set.
#'
#' @details
#' Combines results from  \link[SimRescueMeds]{mcse_bias}, \link[SimRescueMeds]{mcse_empse},
#' \link[SimRescueMeds]{mcse_modse},  \link[SimRescueMeds]{mcse_mse} into a single data frame.
#'
#' @return A data frame.
#'
#' @examples
#' set.seed(283964)
#' smry_ex <- data.frame(
#'   beta_hat = rnorm(100),
#'   se_hat = rgamma(100, shape = 1)
#' )
#'
#' mcse_estimands(
#'   data = smry_ex,
#'   true_value = 1,
#'   estimand_name = "beta_hat",
#'   se_name = "se_hat"
#' )
#'
#' mcse_estimands(
#'   data = smry_ex,
#'   true_value = 0,
#'   estimand_name = "beta_hat",
#'   se_name = "se_hat"
#' )
#'
#' @references
#' \insertRef{Morris2019}{SimRescueMeds}
#'
#' @export
mcse_estimands <- function(
    data,
    true_value = 0,
    estimand_name = "x_diff",
    se_name = "se_diff"
) {
  M <- nrow(data)
  measure_obs <- data[[estimand_name]]
  se_obs <- data[[se_name]]
  mean_value <- mean(measure_obs, na.rm = TRUE)

  bias_df <- mcse_bias(
    measure_obs = measure_obs,
    true_value = true_value,
    mean_value = mean_value,
    M = M
  )

  empse_df <- mcse_empse(
    measure_obs = measure_obs,
    mean_value = mean_value,
    M = M
  )

  modse_df <- mcse_modse(
    se_obs = se_obs,
    M = M
  )

  mse_df <- mcse_mse(
    measure_obs = measure_obs,
    true_value = true_value,
    M = M
  )

  return(do.call('rbind', list(bias_df, empse_df, modse_df, mse_df)))
}
