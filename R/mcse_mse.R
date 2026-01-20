#' Calculate MCSE for Mean Squared Error
#'
#' Calculate MCSE for mean squared error of an estimator based on Table 6 of Morris (2019).
#'
#' @param measure_obs A numeric vector where each element is an estimated value from
#' a simulated data set.
#' @param true_value Numeric. The true value of the estimand.
#' @param M Numeric. The number of simulated data sets. If \code{NULL}, then calculated
#' internally.
#' @param jackknife Logical. When \code{TRUE}, calculates the jackknife estimate of MCSE based
#' on Koehler (2009).
#'
#' @details
#' This function calculates the Monte Carlo standard error (MCSE) for the mean squared error (MSE)
#' of an estimator for an estimand from a set of \eqn{M} simulated data sets. If \eqn{\theta} is
#' the true value of the estimand and \eqn{\hat{\theta}_m} is the estimate in data set \eqn{m},
#' then estimate of MSE is:
#' \deqn{
#'    \widehat{MSE} = \frac{1}{M}\sum_{m = 1}^{M}(\hat{\theta}_m - \theta)^2
#' }
#' The MCSE of the MSE is:
#' \deqn{
#'    MCSE_{MSE} = \sqrt{\frac{\sum_{m=1}^{M}[(\hat{\theta}_m -\theta)^2 - \widehat{MSE}]^2}{M(M-1)}}
#' }
#'
#' @examples
#' set.seed(283964)
#' smry_ex <- data.frame(
#'   beta_hat = rnorm(100),
#'   se_hat = rgamma(100, shape = 1)
#' )
#'
#' mcse_mse(
#'   measure_obs = smry_ex$beta_hat,
#'   true_value = 1
#' )
#'
#' mcse_mse(
#'   measure_obs = smry_ex$beta_hat,
#'   true_value = 1,
#'   jackknife = TRUE
#' )
#'
#' mcse_mse(
#'   measure_obs = smry_ex$beta_hat,
#'   true_value = 0
#' )
#'
#' mcse_mse(
#'   measure_obs = smry_ex$beta_hat,
#'   true_value = 0,
#'   jackknife = TRUE
#' )
#'
#' @references
#' \insertRef{Koehler2009}{SimRescueMeds}
#'
#' \insertRef{Morris2019}{SimRescueMeds}
#'
#' @export
mcse_mse <- function(
    measure_obs,
    true_value,
    M = NULL,
    jackknife = FALSE
) {
  if (is.null(M) == TRUE) {
    M <- length(measure_obs)
  }
  mse <- (1 / M) * sum((measure_obs - true_value)^2)
  num_mse <- sum(((measure_obs - true_value)^2 - mse)^2)
  if (jackknife == FALSE) {
    mse_mcse <- sqrt(num_mse / (M * (M - 1)))
  } else {
    mse_mcse <- mcse_jackknife(
      measure_obs = (measure_obs - true_value)^2
    )
  }

  return(
    data.frame(
      measure = "MSE",
      estimate = mse,
      mcse = mse_mcse
    )
  )
}
