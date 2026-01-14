#' Calculate MCSE for Empirical SE
#'
#' Calculate MCSE for empirical standard error of an estimator based on Table 6 of Morris (2019).
#'
#' @param measure_obs A numeric vector where each element is an estimated value from
#' a simulated data set.
#' @param mean_value The mean value of the estimates. If \code{NULL}, then calculated
#' internally.
#' @param M Numeric. The number of simulated data sets. If \code{NULL}, then calculated
#' internally.
#'
#' @details
#' This function calculates the Monte Carlo standard error (MCSE) for the standard error estimated
#' of an estimator of an estimand \eqn{\theta} from a set of \eqn{M} simulated data sets. If \eqn{\bar{\theta}}
#' is the mean estimate from the simulation study and \eqn{\hat{\theta}_m} is the estimate of is the
#' estimate in the data set \eqn{m}, then empirical standard error (EmpSE) of the estimator is given by:
#' \deqn{
#'    \widehat{EmpSE} = \sqrt{\frac{1}{M-1}\sum_{m = 1}^{M}(\hat{\theta}_m-\bar{\theta})^2}
#' }
#' Then the MCSE for of EmpSE is given by:
#'
#' \deqn{
#'    MCSE_{EmpSE} = \frac{\widehat{EmpSE}}{\sqrt{2(M - 1)}}
#' }
#'
#' @examples
#' set.seed(283964)
#' smry_ex <- data.frame(
#'   beta_hat = rnorm(100),
#'   se_hat = rgamma(100, shape = 1)
#' )
#'
#' mcse_empse(
#'   measure_obs = smry_ex$beta_hat
#' )
#'
#' @references
#' \insertRef{Morris2019}{SimRescueMeds}
#'
#' @export
mcse_empse <- function(
    measure_obs,
    mean_value = NULL,
    M = NULL
) {
  if (is.null(M) == TRUE) {
    M <- length(measure_obs)
  }
  if (is.null(mean_value) == TRUE) {
    mean_value <- mean(measure_obs)
  }
  sum_EmpSE <- sum((measure_obs - mean_value)^2)
  EmpSE <- sqrt((1 / (M - 1)) * sum_EmpSE)
  EmpSE_mcse <- EmpSE / sqrt(2 * (M - 1))
  return(
    data.frame(
      measure = "EmpSE",
      estimate = EmpSE,
      mcse = EmpSE_mcse
    )
  )
}
