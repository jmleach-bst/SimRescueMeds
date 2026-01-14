#' Extract Time-to-Event Summaries for Bayesian Fitted Joint Longitudinal and Time-to-Event Models
#'
#' Extract treatment effect summaries from the time-to-event sub-model for a single data set
#' analyzed with joint longitudinal and time-to-event models using \link[JMbayes2]{jm}.
#'
#' @param JMbayes2_object An An object of class \code{"jointModel"} or \code{"summary.jointModel"}, produced by fitting
#' a joint longitudinal and time-to-event mode \link[JMbayes2]{jm} to analyze a data set.
#' @param JMbayes2_object_summary Logical. When \code{TRUE} it is expected that \code{summary()} has alreadly been
#' applied to the results of fitting \link[JMbayes2]{jm}.
#' @param alpha Numeric. The level specified for hypothesis testing and confidence intervals.
#' @param true_value Numeric. The true value of association parameter of interest. This may be
#' the association between the time-varying covariate from the longitudinal sub-model or perhaps
#' treatment assignment.
#' @param predictor_name Character. The name of the treatment effect in the time-to-event submodel.
#'
#' @importFrom methods is
#'
#' @details
#' This function is narrowly designed to specifically extract (1) a treatment
#' effect estimate in the time-to-event sub-model, (2) standard errors for the
#' treatment effect, (3) confidence intervals for the treatment effect, and (4) a
#' p-value for the treatment effect. Based on the level specified in \code{alpha}
#' argument, whether the null hypothesis is rejected and whether the CI contains the
#' "true" value are also included. We additionally output \code{Rhat} for diagnostics.
#'
#' The idea is to, e.g., use \link[base]{lapply} or \link[purrr]{map} with this function
#' to apply it to a list of \link[JMbayes2]{jm} objects to ultimately produce a data frame
#' where each row contains the results from a single simulated data set.
#'
#' @export
output_JMbayes2_rct_surv <- function(
    JMbayes2_object,
    JMbayes2_object_summary = FALSE,
    predictor_name = "value(y)",
    alpha = 0.05,
    true_value
) {
  if (JMbayes2_object_summary == FALSE) {
    if (is(JMbayes2_object, "jm") == FALSE) {
      stop('Expecting object of class "jm".')
    }
    smry <- summary(JMbayes2_object)
  } else {
    if (is(JMbayes2_object, "summary.jm") == FALSE) {
      stop('Expecting object of class "summary.jm".')
    }
    smry <- JMbayes2_object
  }
  df <- smry$Survival[predictor_name,]
  colnames(df) <- c("effect", "se_effect", "lcl", "ucl", "p", "Rhat")
  df$rejected = ifelse(df$p < alpha, 1, 0)
  df$covered = ifelse(df$lcl <= true_value & true_value <= df$ucl, 1, 0)
  rownames(df) <- NULL
  df <- df[, c("effect", "se_effect", "lcl", "ucl", "p", "rejected", "covered", "Rhat")]
  return(df)
}
