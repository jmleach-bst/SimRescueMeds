#' Calculate MCSE for Coverage
#'
#' Calculate MCSE for coverage of an interval estimator based on Table 6 of Morris (2019).
#'
#' @param data A data frame where each row contains results from a unique simulated data set.
#' @param true_value Numeric. The true value of the estimand.
#' @param lcl_name Character. The variable name for the column containing the estimates
#' of the lower confidence limit from each simulated data set.
#' @param ucl_name Character. The variable name for the column containing the estimates
#' of the upper confidence limit from each simulated data set.
#' @param measure_label Character. What value should be assigned to column \code{measure}?
#' Most users will not need to change the default of \code{"Coverage"}.
#' @param jackknife Logical. When \code{TRUE}, calculates the jackknife estimate of MCSE based
#' on Koehler (2009).
#'
#' @details
#' If If \eqn{\theta} is the true value of the estimand, \eqn{\hat{\theta}_{L,m}} is the lower
#' confidence limit in data set \eqn{m}, \eqn{\hat{\theta}_{U,m}} is the upper confidence limit
#' in data set \eqn{m}, then the estimated coverage, \eqn{\widehat{Cover}} is the proportion of
#' simulated data sets in which \eqn{\hat{\theta}_{L,m} \le \theta \le \hat{\theta}_{U,m}}.
#' The MCSE of coverage is then
#' \deqn{
#'    MCSE_{Cover} = \sqrt{\frac{\widehat{Cover} \times (1 - \widehat{Cover})}{M}}
#' }
#'
#' @examples
#' set.seed(283964)
#' smry_ex <- data.frame(
#'   beta_hat = rnorm(100),
#'   se_hat = rgamma(100, shape = 1)
#' )
#'
#' smry_ex$lcl <- smry_ex$beta_hat - 1.96
#' smry_ex$ucl <- smry_ex$beta_hat + 1.96
#'
#' mcse_coverage(
#'   data = smry_ex,
#'   true_value = 0,
#'   lcl_name = "lcl",
#'   ucl_name = "ucl"
#' )
#'
#' mcse_coverage(
#'   data = smry_ex,
#'   true_value = 1,
#'   lcl_name = "lcl",
#'   ucl_name = "ucl"
#' )
#'
#' mcse_coverage(
#'   data = smry_ex,
#'   true_value = 1,
#'   lcl_name = "lcl",
#'   ucl_name = "ucl",
#'   jackknife = TRUE
#' )
#'
#' @references
#' \insertRef{Koehler2009}{SimRescueMeds}
#'
#' \insertRef{Morris2019}{SimRescueMeds}
#'
#' @export
mcse_coverage <- function(
    data,
    true_value = 0,
    lcl_name,
    ucl_name,
    measure_label = "Coverage",
    jackknife = FALSE
) {
  M <- nrow(data)
  lower_ci <- data[[lcl_name]]
  upper_ci <- data[[ucl_name]]
  covered <- ifelse(true_value >= lower_ci & true_value <= upper_ci, 1, 0)
  coverage <- mean(covered, na.rm = TRUE)
  if (jackknife == FALSE) {
    coverage_mcse <- sqrt((coverage * (1 - coverage)) / M)
  } else {
    coverage_mcse <- mcse_jackknife(measure_obs = covered)
  }

  return(
    data.frame(
      measure = measure_label,
      estimate = coverage,
      mcse = coverage_mcse
    )
  )
}
