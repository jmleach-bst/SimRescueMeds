#' Calculate MCSE for Rejection Rate
#'
#' Calculate MCSE for rejection rate of a hypothesis test based on Table 6 of Morris (2019).
#'
#' @param data A numeric vector of either p-values or binary such that rejected = 1 or not rejected = 0.
#' @param type Character. \code{"pvalues"} if \code{data} is a vector of p-values and \code{"rejection"}
#' if a binary vector.
#' @param alpha Numeric. The desired significance level, which should be greater than 0 and less than 1.
#' @param measure_label Character. What value should be assigned to column \code{measure}?
#' Most users will not need to change the default of \code{"Rejection"}.
#' @param jackknife Logical. When \code{TRUE}, calculates the jackknife estimate of MCSE based
#' on Koehler (2009).
#'
#' @details
#' If If \eqn{\theta} is the true value of the estimand and \eqn{p_m} is a p-value for the hypothesis
#' test evaluated for simulated data set \eqn{m}, and \eqn{\alpha} is the desired significance level,
#' then the estimated rejection rate, \eqn{\widehat{Reject}} is the proportion of simulated data sets
#' in which \eqn{p_m < \alpha}. The MCSE of coverage is then
#' \deqn{
#'    MCSE_{Reject} = \sqrt{\frac{\widehat{Reject} \times (1 - \widehat{Reject})}{M}}
#' }
#'
#' @examples
#' set.seed(8972364)
#' mcse_rejection(
#'   data = runif(100),
#'   type = "pvalues"
#' )
#'
#' rb2 <- rbinom(100, size = 1, prob = 0.15)
#'
#' mcse_rejection(
#'   data = rb2,
#'   type = "rejection"
#' )
#'
#' mcse_rejection(
#'   data = rb2,
#'   type = "rejection",
#'   jackknife = TRUE
#' )
#'
#' @references
#' \insertRef{Koehler2009}{SimRescueMeds}
#'
#' \insertRef{Morris2019}{SimRescueMeds}
#'
#' @export
mcse_rejection <- function(
    data,
    type = "pvalues", # alternatively "rejection"
    alpha = 0.05,
    measure_label = "Rejection",
    jackknife = FALSE
    ) {
      M <- length(data)
      if (type == "pvalues") {
        reject <- ifelse(data < alpha, 1, 0)
      } else {
        reject <- data
      }
      reject_rate <- mean(reject, na.rm = TRUE)
      if (jackknife == FALSE) {
        reject_rate_mcse <- sqrt((reject_rate * (1 - reject_rate)) / M)
      } else {
        reject_rate_mcse <- mcse_jackknife(reject)
      }

      return(
        data.frame(
          measure = measure_label,
          estimate = reject_rate,
          mcse = reject_rate_mcse
        )
      )
    }
