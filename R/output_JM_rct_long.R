#' Extract Longitudinal Summaries for Fitted Joint Longitudinal and Time-to-Event Models
#'
#' Extract treatment effect summaries from the longitudinal sub-model for a single data set
#' analyzed with joint longitudinal and time-to-event models using \link[JM]{jointModel}.
#'
#' @param JM_object An An object of class \code{"jointModel"} or \code{"summary.jointModel"}, produced by fitting
#' a joint longitudinal and time-to-event mode \link[JM]{jointModel} to analyze a data set.
#' @param JM_object_summary Logical. When \code{TRUE} it is expected that \code{summary()} has alreadly been
#' applied to the results of fitting \link[JM]{jointModel}.
#' @param alpha Numeric. The level specified for hypothesis testing and confidence intervals.
#' @param true_beta Numeric. The true value of the treatment effect, usually the interaction
#' between time and treatment assignment.
#' @param predictor_name Character. The name of the treatment effect as specified in the
#' \code{formula} argument in \link[nlme]{lme}.
#'
#' @importFrom stats qnorm
#' @importFrom stats confint
#' @importFrom methods is
#'
#' @details
#' This function is narrowly designed to specifically extract (1) a treatment
#' effect estimate, usually the interaction between time and treatment group in
#' longitudinal sub-model, (2) standard errors for the treatment effect,
#' (3) confidence intervals for the treatment effect, and (4) a p-value for the
#' treatment effect. Based on the level specified in \code{alpha} argument, whether
#' the null hypothesis is rejected and whether the CI contains the "true" value are
#' also included. The idea is to, e.g., use \link[base]{lapply} or \link[purrr]{map}
#' with this function to apply it to a list of \link[JM]{jointModel} objects and produce a
#' data frame where each row contains the results from a single simulated data set.
#'
#' @export
output_JM_rct_long <- function(
    JM_object,
    JM_object_summary = FALSE,
    predictor_name = "x:time",
    alpha = 0.05,
    true_beta
) {
  if (JM_object_summary == FALSE) {
    if (is(JM_object, "jointModel") == FALSE) {
      stop('Expecting object of class "jointModel".')
    }
    ci_effect <- confint(
      object = JM_object,
      level = 1 - alpha
    )[paste0("Y.", predictor_name),]
    df <- data.frame(
      effect = ci_effect["est."],
      se_effect = summary(JM_object)$`CoefTable-Long`[predictor_name, "Std.Err"],
      lcl = ci_effect["2.5 %"],
      ucl = ci_effect["97.5 %"],
      p = summary(JM_object)$`CoefTable-Long`[predictor_name, "p-value"]
    )
    df$rejected = ifelse(df$p < alpha, 1, 0)
    df$covered = ifelse(df$lcl <= true_beta & true_beta <= df$ucl, 1, 0)
    rownames(df) <- NULL
    return(df)
  } else {
    if (is(JM_object, "summary.jointModel") == FALSE) {
      stop('Expecting object of class "summary.jointModel".')
    }
    df <- data.frame(
      effect = JM_object$`CoefTable-Long`[predictor_name, "Value"],
      se_effect = JM_object$`CoefTable-Long`[predictor_name, "Std.Err"]
    )
    df$lcl <- df$effect - qnorm(1 - alpha/2)*df$se_effect
    df$ucl <- df$effect + qnorm(1 - alpha/2)*df$se_effect
    df$p <- JM_object$`CoefTable-Long`[predictor_name, "p-value"]
    df$rejected = ifelse(df$p < alpha, 1, 0)
    df$covered = ifelse(df$lcl <= true_beta & true_beta <= df$ucl, 1, 0)
    rownames(df) <- NULL
    return(df)
  }

}
