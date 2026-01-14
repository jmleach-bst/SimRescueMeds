#' Extract Summaries for Fitted Linear Models
#'
#' Extract treatment effect summaries for a single data set analyzed with \link[stats]{lm}.
#'
#' @param lm_object An object of class \code{lm}, produce by fitting a linear model using \link[stats]{lm}
#' to analyze a data set.
#' @param alpha Numeric. The level specified for hypothesis testing and confidence intervals.
#' @param true_beta Numeric. The true value of the treatment effect, usually the interaction
#' between time and treatment assignment.
#' @param predictor_name Character. The name of the treatment effect as specified in the
#' \code{formula} argument in \link[stats]{lm}.
#'
#' @importFrom methods is
#'
#' @details
#' This function is narrowly designed to specifically extract (1) a treatment
#' effect estimate, (2) standard errors for the treatment effect, (3) confidence
#' intervals for the treatment effect, and (4) a p-value for the treatment effect.
#' Based on the level specified in \code{alpha} argument, whether the null hypothesis is
#' rejected and whether the CI contains the "true" value are also included. The idea is
#' to, e.g., use \link[base]{lapply} or \link[purrr]{map} with this function to
#' apply it to a list of \link[stats]{lm} objects and produce a data frame where
#' each row contains the results from a single simulated data set.
#'
#' @return A data frame.
#'
#' @export
output_lm_rct <- function(
    lm_object,
    predictor_name = "x",
    alpha = 0.05,
    true_beta
) {
  if (is(lm_object, "lm") == FALSE) {
    stop('Expect object of class "lm".')
  }
  ci_effect <- confint.default(
    object = lm_object,
    parm = predictor_name,
    level = 1 - alpha
  )
  df <- data.frame(
    effect = lm_object$coefficients[predictor_name],
    se_effect = summary(lm_object)$coefficients[predictor_name, "Std. Error"],
    lcl = ci_effect[1],
    ucl = ci_effect[2],
    p = summary(lm_object)$coefficients[predictor_name, "Pr(>|t|)"]
  )
  df$rejected = ifelse(df$p < alpha, 1, 0)
  df$covered = ifelse(df$lcl <= true_beta & true_beta <= df$ucl, 1, 0)
  rownames(df) <- NULL
  return(df)
}
