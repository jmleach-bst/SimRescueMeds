#' Simulate Longitudinal RCT Data Set via Linear Mixed Model
#'
#' Simulate continuous-valued longitudinal data for RCT using linear mixed model
#' and subject-specific MVN formulation.
#'
#' @inheritParams build_design_matrix
#' @inheritParams simulate_subject_MVN
#' @param xmat_t A design matrix for those on active treatment.
#' @param xmat_c A design matrix for those on control treatment.
#'
#' @return A data frame containing simulated data for a single
#' data set.
#'
#' @details
#' This function simply uses \link[SimRescueMeds]{simulate_subject_MVN} and
#' \link[base]{replicate} to generate an entire data set. Each treatment group
#' assumes the same design matrix for all subjects and the dimension of the design
#' matrices for each group should be equal. That is, they should be followed for
#' the same number of time points and have the same covariates (usually omitted).
#' In short, this function is narrowly written for longitudinal RCT-style data.
#' While subject-specific MVN generation may be slightly more involved to program
#' than whole-data-set MVN generation, the former is likely to be more computationally
#' efficient, especially as the sample size increases.
#'
#' @examples
#'
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
#' simulate_lmm_rct(
#'   N_t = 3,
#'   N_c = 4,
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
#'
#' @export
simulate_lmm_rct <- function(
    N_t,
    N_c,
    xmat_t,
    xmat_c,
    betas,
    Sigma
) {
  if (all(dim(xmat_t) != dim(xmat_c))) {
    stop("Design matrices should all have same dimension")
  }
  K <- nrow(xmat_t)
  df_t <- replicate(
    n = N_t,
    expr = simulate_subject_MVN(
      xmat = xmat_t,
      betas = betas,
      Sigma = Sigma
    ),
    simplify = FALSE
  )
  df_t <- do.call('rbind', df_t)
  df_c <- replicate(
    n = N_c,
    expr = simulate_subject_MVN(
      xmat = xmat_c,
      betas = betas,
      Sigma = Sigma
    ),
    simplify = FALSE
  )
  df_c <- do.call('rbind', df_c)
  df <- rbind(df_t, df_c)
  df$id <- factor(rep(1:(N_t + N_c), each = K))
  return(df)
}
