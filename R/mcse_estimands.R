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
#' @param jackknife Logical vector that determines whether or not jackknife estimator is
#' used for MCSE. Default is \code{FALSE} for all. Note that you need to correctly name
#' the elements of the vector to avoid errors. See examples.
#' @param include_bias_percent Logical. Default is \code{FALSE}. If \code{true_value = 0},
#' percent-bias is undefined and will not be returned, although a \code{warning} will
#' be thrown.
#'
#' @inheritParams mcse_bias_percent
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
#'   true_value = 0.5,
#'   estimand_name = "beta_hat",
#'   se_name = "se_hat",
#'   include_bias_percent = TRUE,
#'   report_proportion = TRUE
#' )
#'
#' mcse_estimands(
#'   data = smry_ex,
#'   true_value = 0.5,
#'   estimand_name = "beta_hat",
#'   se_name = "se_hat",
#'   include_bias_percent = TRUE,
#'   report_proportion = FALSE
#' )
#'
#' ## How to specify jackknife, e.g., for ModSE
#' mcse_estimands(
#'   data = smry_ex,
#'   true_value = 1,
#'   estimand_name = "beta_hat",
#'   se_name = "se_hat",
#'   jackknife = c(
#'     bias = FALSE,
#'     empse = FALSE,
#'     modse = TRUE,
#'     mse = FALSE
#'     )
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
#' \insertRef{Koehler2009}{SimRescueMeds}
#'
#' \insertRef{Morris2019}{SimRescueMeds}
#'
#' @export
mcse_estimands <- function(
    data,
    true_value = 0,
    estimand_name = "x_diff",
    se_name = "se_diff",
    jackknife = c(
      bias = FALSE,
      empse = FALSE,
      modse = FALSE,
      mse = FALSE
    ),
    include_bias_percent = FALSE,
    report_proportion = FALSE
) {
  M <- nrow(data)
  measure_obs <- data[[estimand_name]]
  se_obs <- data[[se_name]]
  mean_value <- mean(measure_obs, na.rm = TRUE)

  if (include_bias_percent == TRUE) {
    if (true_value == 0) {
      warning("When true_value = 0, percent-bias is undefined")
      bias_percent_df <- NULL
    } else {
      bias_percent_df <- mcse_bias_percent(
        measure_obs = measure_obs,
        true_value = true_value,
        mean_value = mean_value,
        M = M,
        report_proportion = report_proportion
      )
    }
  } else {
    bias_percent_df <- NULL
  }

  bias_df <- mcse_bias(
    measure_obs = measure_obs,
    true_value = true_value,
    mean_value = mean_value,
    M = M,
    jackknife = jackknife["bias"]
  )

  empse_df <- mcse_empse(
    measure_obs = measure_obs,
    mean_value = mean_value,
    M = M,
    jackknife = jackknife["empse"]
  )

  modse_df <- mcse_modse(
    se_obs = se_obs,
    M = M,
    jackknife = jackknife["modse"]
  )

  mse_df <- mcse_mse(
    measure_obs = measure_obs,
    true_value = true_value,
    M = M,
    jackknife = jackknife["mse"]
  )

  return(do.call('rbind', list(bias_percent_df, bias_df, empse_df, modse_df, mse_df)))
}
