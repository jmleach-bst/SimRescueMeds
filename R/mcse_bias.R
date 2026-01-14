#' Calculate MCSE for Bias
#'
#' Calculate MCSE for bias based on Table 6 of Morris (2019).
#'
#' @param measure_obs A numeric vector where each element is an estimated value from
#' a simulated data set.
#' @param mean_value Numeric. The mean value of the estimates. If \code{NULL}, then calculated
#' internally.
#' @param true_value Numeric. The true value of the estimand.
#' @param M Numeric. The number of simulated data sets. If \code{NULL}, then calculated
#' internally.
#'
#' @details
#' This function calculates the Monte Carlo standard error (MCSE) for the bias of an estimator
#' for an estimand from a set of \eqn{M} simulated data sets. If \eqn{\theta} is the true value
#' of the estimand, \eqn{\bar{\theta}} is the mean estimate from the simulation study, and
#' \eqn{\hat{\theta}_m} is the estimate in data set \eqn{m}, then estimated bias is:
#' \deqn{
#'    \widehat{Bias} = \frac{1}{M}\sum_{m=1}^{M}\hat{\theta}_i - \theta
#' }
#' The MCSE for \eqn{\widehat{Bias}} is given by:
#' \deqn{
#'    MCSE_{bias} = \sqrt{\frac{1}{M(M-1)}\sum_{m = 1}^{M}(\hat{\theta}_m-\bar{\theta})^2}
#' }
#'
#' @examples
#' set.seed(283964)
#' smry_ex <- data.frame(
#'   beta_hat = rnorm(100),
#'   se_hat = rgamma(100, shape = 1)
#' )
#'
#' mcse_bias(
#'   measure_obs = smry_ex$beta_hat,
#'   mean_value = NULL,
#'   true_value = 1
#' )
#'
#' mcse_bias(
#'   measure_obs = smry_ex$beta_hat,
#'   mean_value = NULL,
#'   true_value = 0
#' )
#'
#' @references
#' \insertRef{Morris2019}{SimRescueMeds}
#'
#' @export
mcse_bias <- function(
    measure_obs,
    true_value,
    mean_value = NULL,
    M = NULL
) {
  if (is.null(M) == TRUE) {
    M <- length(measure_obs)
  }
  if (is.null(mean_value) == TRUE) {
    mean_value <- mean(measure_obs)
  }
  bias <- (1 / M) * sum(measure_obs - true_value)
  mult <- 1 / (M * (M - 1))
  sum_msq <- sum((measure_obs - mean_value)^2)
  bias_mcse <- sqrt(mult * sum_msq)
  return(
    data.frame(
      measure = "Bias",
      estimate = bias,
      mcse = bias_mcse
    )
  )
}

