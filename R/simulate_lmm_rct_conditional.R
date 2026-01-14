#' Simulate Longitudinal RCT Data Set via Linear Mixed Model
#'
#' Simulate continuous-valued longitudinal data for RCT using linear mixed model
#' and subject-specific MVN formulation, while saving generated random effects
#' and within-subject errors.
#'
#' @inheritParams build_design_matrix
#' @inheritParams simulate_lmm_rct
#' @inheritParams simulate_subject_MVN_conditional
#' @param zmat_t The random effects design matrix for the active treatment group.
#' @param zmat_c The random effects design matrix for the control treatment group.
#'
#' @details
#' This function uses \link[SimRescueMeds]{simulate_subject_MVN_conditional} to generate
#' entire data sets. It is similar to \link[SimRescueMeds]{simulate_lmm_rct}, except that
#' it expliticly produces and saves random effects and within-subject errors.
#'
#' @examples
#' set.seed(78634)
#' x_ex <- build_design_matrix(
#'   N_t = 1,
#'   N_c = 0,
#'   K = 5,
#'   time_start = 0,
#'   time_scale = 1/4
#' )
#' z_ex <- x_ex[, c("intercept", "time")]
#'
#' cov_ex <- build_lmm_cov(
#'   zmat = z_ex,
#'   re_corr_mat = matrix(
#'     c(1, -0.15,
#'       -0.15, 1),
#'     byrow = TRUE, nrow = 2
#'   ),
#'   re_sigma = c(2, 3),
#'   ws_corr_mat = diag(1, 5),
#'   ws_sigma = 4,
#'   partition = TRUE,
#'   print_intermediate = FALSE
#' )
#'
#' x_ex0 <- build_design_matrix(
#'   N_t = 0,
#'   N_c = 1,
#'   K = 5,
#'   time_start = 0,
#'   time_scale = 1/4
#' )
#' z_ex0 <- x_ex[, c("intercept", "time")]
#'
#' simulate_lmm_rct_conditional(
#'   N_t = 2,
#'   N_c = 1,
#'   xmat_t = x_ex,
#'   xmat_c = x_ex0,
#'   zmat_t = z_ex,
#'   zmat_c = z_ex0,
#'   betas = c(
#'     beta0 = 8,
#'     beta1 = 0,
#'     beta2 = -1,
#'     beta3 = -1
#'   ),
#'   re_cov_mat = cov_ex[["re_cov_mat"]],
#'   ws_cov_mat = cov_ex[["ws_cov_mat"]],
#'   re_names = c("b0", "b2")
#' )
#'
#'
#' @export
simulate_lmm_rct_conditional <- function(
    N_t,
    N_c,
    xmat_t,
    xmat_c,
    betas,
    zmat_t,
    zmat_c,
    re_cov_mat = 1,
    ws_cov_mat,
    re_names = NULL
) {
  if (all(dim(xmat_t) != dim(xmat_c))) {
    stop("Design matrices should all have same dimension")
  }
  K <- nrow(xmat_t)
  df_t <- replicate(
    n = N_t,
    expr = simulate_subject_MVN_conditional(
      xmat = xmat_t,
      zmat = zmat_t,
      betas = betas,
      re_cov_mat = re_cov_mat,
      ws_cov_mat = ws_cov_mat,
      re_names = re_names
    ),
    simplify = FALSE
  )
  df_t <- do.call('rbind', df_t)
  df_c <- replicate(
    n = N_c,
    expr = simulate_subject_MVN_conditional(
      xmat = xmat_c,
      zmat = zmat_c,
      betas = betas,
      re_cov_mat = re_cov_mat,
      ws_cov_mat = ws_cov_mat,
      re_names = re_names
    ),
    simplify = FALSE
  )
  df_c <- do.call('rbind', df_c)
  df <- rbind(df_t, df_c)
  df$id <- factor(rep(1:(N_t + N_c), each = K))
  return(df)
}
