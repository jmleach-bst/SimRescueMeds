#' Calculate Jackknife MCSE
#'
#' Calculate MCSE based on a jackknife estimator.
#'
#' @param measure_obs A numeric vector where each element is an estimated value from
#' a simulated data set.
#' @param user_function In some cases it will necessary to define a function to take over
#' \code{measure_obs}. For example, if we need to calculate the empirical SE of treatment
#' effect estimates, then this is usually calculated as the standard deviation (SD) of the
#' estimates. For a jackknife MCSE of empirical SE, we need to recalculate the SD \eqn{M}
#' times, each time excluding one of the treatment effect estimates.
#'
#' @details
#' This function calculates the Monte Carlo standard error (MCSE) using a jackknife estimator
#' described initially in Koehler (2009) and applied more generally in Appendix A.3 of Kelter (2024).
#' For \eqn{M} simulated data sets, let \eqn{\theta} be the true value of an estimand, \eqn{\hat{\theta}}
#' the estimate (here the mean) over all simulated data sets, and \eqn{\hat{\theta}_{(-m)}} is a
#' recalculation of \eqn{\hat{\theta}} with data set \eqn{m} left out. The MCSE jackknife estimate
#' is given by
#' \deqn{
#'  MCSE_{JK} = \sqrt{\frac{M - 1}{M}\sum_{m = 1}^{M}(\hat{\theta}_{(-m)} - \bar{\theta}_{JK})^2}
#' }
#' where \eqn{\bar{\theta}_{JK}} is the mean of the \eqn{\hat{\theta}_{(-m)}}:
#' \deqn{
#'  \bar{\theta}_{JK} = \frac{1}{M} \sum_{m = 1}^{M} \hat{\theta}_{(-m)}
#' }
#' In many cases, all that is required to calculate the MCSE in this way is the vector of
#' observed estimates, which are specified using \code{measure_obs}. In most cases the user
#' will not need to call this function directly, as it is incorporated into relevant functions
#' as an option, or in some cases, may be the only option if no formula is available for MCSE.
#'
#' @references
#' \insertRef{Koehler2009}{SimRescueMeds}
#'
#' \insertRef{Kelter2024}{SimRescueMeds}
#'
#' @export
mcse_jackknife <- function(
    measure_obs,
    user_function = NULL
) {
  if (!is.numeric(measure_obs)) {
    stop("Requires numeric vector input.")
  }
  M <- length(measure_obs)
  measure_obs_jk <- c()
  for (m in 1:M) {
    if (is.null(user_function)) {
      measure_obs_jk[m] <- mean(measure_obs[-m], na.rm = TRUE)
    } else {
      measure_obs_jk[m] <- user_function(x = measure_obs[-m])
    }

  }
  mult_jk <- (M - 1)/M
  mean_jk <- mean(measure_obs_jk, na.rm = TRUE)
  sum_jk <- sum((measure_obs_jk - mean_jk)^2)
  return(sqrt(mult_jk*sum_jk))
}
