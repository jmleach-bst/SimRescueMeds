#' Simulate Single-Subject Outcomes via Linear Mixed Model
#'
#' Use \link[MASS]{mvrnorm} to draw single-subject data from MVN based on LMM.
#'
#' @param xmat A fixed effects design matrix.
#' @param betas A numeric vector of fixed effects parameters whose length should equal \code{ncol(xmat)}.
#' @param Sigma A covariance matrix. Recommended to generate from \link[SimRescueMeds]{build_lmm_cov}.
#'
#' @details We use \link[MASS]{mvrnorm} to draw from a multivariate Normal distribution. The model used is
#' \eqn{\mathbf{Y}_i \sim MVN(\mathbf{X}_i\mathbf{\beta}, \mathbf{\Sigma}_i)} where \eqn{\mathbf{X}_i}
#' is the subject-specific fixed effects design matrix, \eqn{\mathbf{\beta}} is the fixed effects parameter vector,
#' \eqn{\mathbf{\Sigma}_i = \mathbf{Z}_i\mathbf{\Sigma}_b\mathbf{Z}_i^\top + \mathbf{\Sigma}_\epsilon},
#' \eqn{\mathbf{Z}_i} is the subject-specific random effects design matrix, \eqn{\mathbf{\Sigma}_b}is the random
#' effects covariance matrix, \eqn{\mathbf{\Sigma}_\epsilon} is the within-subject random error covariance matrix.
#'
#' @importFrom MASS mvrnorm
#'
#' @return A data frame containing both the generated outcomes and the design matrix from \code{xmat},
#' although, we drop the vector of ones for the intercept in the output.
#'
#' @examples
#' xmat_trt <- build_design_matrix(
#'   N_t = 1,
#'   N_c = 0,
#'   K = 5,
#'   time_start = 0,
#'   time_scale = 1/4
#' )
#'
#' Sigma_i <- build_lmm_cov(
#' zmat = xmat_trt[, c("intercept", "time")],
#' re_sigma = c(1.15, 1.05),
#' re_corr_mat = matrix(
#'   c(1, 0,
#'     0, 1),
#'   byrow = TRUE,
#'   nrow = 2
#' ),
#' ws_sigma = 1.25,
#' ws_corr_mat = diag(5),
#' print_intermediate = FALSE
#' )
#'
#' set.seed(32857)
#' df_i <- simulate_subject_MVN(xmat = xmat_trt,
#'                             betas = c(
#'                               beta0 = 8,
#'                               beta1 = 0,
#'                               beta2 = -1,
#'                               beta3 = -1
#'                              ),
#'                             Sigma = Sigma_i
#'                             )
#' print(df_i)
#'
#' @export
simulate_subject_MVN <- function(
    xmat,
    betas,
    Sigma
) {
  y <- MASS::mvrnorm(
    n = 1,
    mu = xmat %*% betas,
    Sigma = Sigma
  )
  df <- cbind(
    y = y,
    xmat[, -1] # remove column of 1's from analysis data set
  )
  return(as.data.frame(df))
}
