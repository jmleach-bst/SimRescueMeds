#' Calculate MCSE for Model SE
#'
#' Calculate MCSE for model-based standard error of an estimator based on Table 6 of Morris (2019).
#'
#' @param se_obs A numeric vector where each element is an estimated model-based SE value from
#' a simulated data set.
#' @param M Numeric. The number of simulated data sets. If \code{NULL}, then calculated
#' internally.
#' @param jackknife Logical. When \code{TRUE}, calculates the jackknife estimate of MCSE based
#' on Koehler (2009).
#'
#' @importFrom stats var
#'
#' @details
#' This function calculates the Monte Carlo standard error (MCSE) for the model-based standard
#' error based on \eqn{M} simulated data sets. The average model-based SE is the average SE
#' estimate for the estimator:
#' \deqn{
#'    \widehat{ModSE} = \sqrt{\frac{1}{M}\sum_{m = 1}^{M}\widehat{var}(\hat{\theta}_i)}
#' }
#' Then the MCSE for \eqn{\widehat{ModSE}} is
#' \deqn{
#'    MCSE_{ModSE} = \sqrt{\frac{\widehat{var}[\widehat{var}(\hat{\theta})]}{4M\widehat{ModSE}^2}}
#' }
#'
#' @examples
#' set.seed(283964)
#' smry_ex <- data.frame(
#'   beta_hat = rnorm(100),
#'   se_hat = rgamma(100, shape = 1)
#' )
#'
#' mcse_modse(
#'   se_obs = smry_ex$se_hat
#' )
#'
#' mcse_modse(
#'   se_obs = smry_ex$se_hat,
#'   jackknife = TRUE
#' )
#'
#' @references
#' \insertRef{Koehler2009}{SimRescueMeds}
#'
#' \insertRef{Morris2019}{SimRescueMeds}
#'
#' @export
mcse_modse <- function(
    se_obs,
    M = NULL,
    jackknife = FALSE
) {
  if (is.null(M) == TRUE) {
    M <- length(se_obs)
  }
  var_obs <- se_obs^2
  ModSE <- sqrt(mean(var_obs))
  if (jackknife == FALSE) {
    vhvhht <- var(var_obs)
    ModSE_mcse <- sqrt(vhvhht / (4*M*ModSE^2))
  } else {
    ModSE_mcse <- mcse_jackknife(
      measure_obs = var_obs,
      user_function = function(x) sqrt(mean(x))
    )
  }

  return(
    data.frame(
      measure = "ModSE",
      estimate = ModSE,
      mcse = ModSE_mcse
    )
  )
}
