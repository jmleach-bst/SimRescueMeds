#' Extract Summaries for Fitted Linear Mixed Models
#'
#' Extract treatment effect summaries for a single data set analyzed with \link[nlme]{lme}.
#'
#' @param lme_object An object of class \code{lme}, produce by fitting a linear mixed
#' model using \link[nlme]{lme} to analyze a data set.
#' @param alpha Numeric. The level specified for hypothesis testing and confidence intervals.
#' @param true_beta Numeric. The true value of the treatment effect, usually the interaction
#' between time and treatment assignment.
#' @param predictor_name Character. The name of the treatment effect as specified in the
#' \code{formula} argument in \link[nlme]{lme}.
#'
#' @importFrom stats confint.default
#' @importFrom methods is
#'
#' @details
#' This function is narrowly designed to specifically extract (1) a treatment
#' effect estimate, usually the interaction between time and treatment group,
#' (2) standard errors for the treatment effect, (3) confidence intervals for
#' the treatment effect, and (4) a p-value for the treatment effect. Based on
#' the level specified in \code{alpha} argument, whether the null hypothesis is rejected
#' and whether the CI contains the "true" value are also included. The idea is
#' to, e.g., use \link[base]{lapply} or \link[purrr]{map} with this function to
#' apply it to a list of \link[nlme]{lme} objects and produce a data frame where
#' each row contains the results from a single simulated data set.
#'
#' @return A data frame.
#'
#' @examples
#' xmat_t <- build_design_matrix(
#'     N_t = 1,
#'     N_c = 0,
#'     K = 5,
#'     time_start = 0,
#'     time_scale = 1/4
#'     )
#' xmat_c <- build_design_matrix(
#'     N_t = 0,
#'     N_c = 1,
#'     K = 5,
#'     time_start = 0,
#'     time_scale = 1/4
#'     )
#'
#' set.seed(2834628)
#'
#' df <- simulate_lmm_rct(
#'   N_t = 50,
#'   N_c = 50,
#'   xmat_t = xmat_t,
#'   xmat_c = xmat_c,
#'   betas = c(
#'     beta0 = 8,
#'     beta1 = 0,
#'     beta2 = -1,
#'     beta3 = -1
#'     ),
#'   Sigma = build_lmm_cov(
#'     zmat = xmat_t[, c("intercept", "time")],
#'     re_sigma = c(1.15, 1.05),
#'     re_corr_mat = matrix(
#'       c(1, 0,
#'         0, 1),
#'       byrow = TRUE,
#'       nrow = 2
#'       ),
#'     ws_sigma = 1.25,
#'     ws_corr_mat = diag(5),
#'     print_intermediate = FALSE
#'     )
#' )
#'
#' lme_df <- nlme::lme(
#'     data = df,
#'     fixed = y ~ x + time + x:time,
#'     random = ~ 1 + time | id
#' )
#'
#' output_lmm_rct(
#'   lme_object = lme_df,
#'   true_beta = -1
#' )
#'
#' @export
output_lmm_rct <- function(
    lme_object,
    alpha = 0.05,
    true_beta,
    predictor_name = "x:time"
) {
  if (is(lme_object, "lme") == FALSE) {
    stop('Expect object of class "lme".')
  }
  ci_effect <- nlme::intervals(lme_object,
                               level = 1 - alpha,
                               which = "fixed")$fixed[predictor_name, c("lower", "upper")]
  df <- data.frame(
    effect = lme_object$coefficients$fixed[predictor_name],
    se_effect = summary(lme_object)$tTable[predictor_name, "Std.Error"],
    lcl = ci_effect[1],
    ucl = ci_effect[2],
    p = summary(lme_object)$tTable[predictor_name, "p-value"]
  )
  df$rejected = ifelse(df$p < alpha, 1, 0)
  df$covered = ifelse(df$lcl <= true_beta & true_beta <= df$ucl, 1, 0)
  rownames(df) <- NULL
  return(df)
}
