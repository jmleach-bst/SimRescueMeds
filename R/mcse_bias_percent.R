#' Calculate MCSE for Percent-Bias
#'
#' Calculate MCSE for percent-bias based on a jackknife estimator.
#'
#' @param measure_obs A numeric vector where each element is an estimated value from
#' a simulated data set.
#' @param mean_value Numeric. The mean value of the estimates. If \code{NULL}, then calculated
#' internally.
#' @param true_value Numeric. The true value of the estimand.
#' @param M Numeric. The number of simulated data sets. If \code{NULL}, then calculated
#' internally.
#' @param report_proportion Logical. It may be preferable to report bias proportion, which is
#' done if \code{report_proportion = TRUE}.
#'
#' @details
#' This function calculates the Monte Carlo standard error (MCSE) for the percent bias of an estimator
#' for an estimand from a set of \eqn{M} simulated data sets. If \eqn{\theta} is the true value
#' of the estimand, \eqn{\bar{\theta}} is the mean estimate from the simulation study, and
#' \eqn{\hat{\theta}_m} is the estimate in data set \eqn{m}, then estimated bias is:
#' \deqn{
#'    \widehat{Bias} = \frac{1}{M}\sum_{m=1}^{M}\hat{\theta}_i - \theta
#' }
#' The percent-bias for simulated data set \eqn{i} and estimator \eqn{\hat{\theta}} is
#' \deqn{
#'  \widehat{Bias}_{percent} = 100 \times \frac{\hat{\theta}_i - \theta}{\theta}
#' }
#' There is no formula for percent-bias in Morris (2019), but Koehler (2009) describes a
#' jackknife estimator that is also applied more broadly in Kelter (2024), which we can use.
#' See \link[SimRescueMeds]{mcse_jackknife} for details.
#'
#' @examples
#' set.seed(283964)
#' smry_ex <- data.frame(
#'   beta_hat = rnorm(1000, 1, 1),
#'   se_hat = rgamma(1000, shape = 1)
#' )
#'
#' mcse_bias_percent(
#'   measure_obs = smry_ex$beta_hat,
#'   mean_value = NULL,
#'   true_value = 1
#' )
#'
#' mcse_bias_percent(
#'   measure_obs = smry_ex$beta_hat,
#'   mean_value = NULL,
#'   true_value = 1,
#'   report_proportion = TRUE
#' )
#'
#' @references
#' \insertRef{Koehler2009}{SimRescueMeds}
#'
#' \insertRef{Morris2019}{SimRescueMeds}
#'
#' \insertRef{Kelter2024}{SimRescueMeds}
#'
#' @export
mcse_bias_percent <- function(
    measure_obs,
    true_value,
    mean_value = NULL,
    M = NULL,
    report_proportion = FALSE
) {
  if (is.null(M) == TRUE) {
    M <- length(measure_obs)
  }
  if (is.null(mean_value) == TRUE) {
    mean_value <- mean(measure_obs)
  }
  if (true_value == 0) {
    stop("Percent bias is undefined when true_value = 0.")
  }
  if (report_proportion == FALSE) {
    pbias_m <- 100*((measure_obs - true_value) / true_value)
    measure_label <- "Percent-Bias"
  } else {
    pbias_m <- (measure_obs - true_value) / true_value
    measure_label <- "Proportion-Bias"
  }

  pbias <- mean(pbias_m, na.rm = TRUE)
  pbias_mcse <- mcse_jackknife(measure_obs = pbias_m)
  return(
    data.frame(
      measure = measure_label,
      estimate = pbias,
      mcse = pbias_mcse
    )
  )
}

