#' Simulate Single-Subject Outcomes via Linear Mixed Model
#'
#' Draw single-subject data from MVN based on LMM and save random effects.
#'
#' @inheritParams simulate_subject_MVN
#' @param zmat A random effect design matrx.
#' @param re_cov_mat A random effects covariance matrix (scalar if only intercept).
#' This can be obtained from \link[SimRescueMeds]{build_lmm_cov} with \code{partition = TRUE}.
#' @param ws_cov_mat A within-subject covariance matrix, usually diagonal. This can be
#' obtained from \link[SimRescueMeds]{build_lmm_cov} with \code{partition = TRUE}.
#' @param re_names An optional vector for variable names for random effects.
#' @param print_intermediate Logical. When \code{TRUE} prints intermediate results.
#'
#' @details
#' While \link[MASS]{mvrnorm} is efficient for drawing single-subject outcomes,
#' it does not explicitly produce or save random effects. For some approaches to
#' simulating time-to-event data that depends on the longitudinal trajector,
#' the random effects may be necessary. Here we use \link[MASS]{mvrnorm} to draw
#' the random effects and within-subject errors separately.
#'
#' @return A data frame.
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
#'
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
#' simulate_subject_MVN_conditional(
#'   xmat = x_ex,
#'   zmat = z_ex,
#'   betas = c(
#'     beta0 = 8,
#'     beta1 = 0,
#'     beta2 = -1,
#'     beta3 = -1
#'   ),
#'   re_cov_mat = cov_ex[["re_cov_mat"]],
#'   ws_cov_mat = cov_ex[["ws_cov_mat"]],
#'   print_intermediate = TRUE
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
#' simulate_subject_MVN_conditional(
#'   xmat = x_ex0,
#'   zmat = z_ex0,
#'   betas = c(
#'     beta0 = 8,
#'     beta1 = 0,
#'     beta2 = -1,
#'     beta3 = -1
#'   ),
#'   re_cov_mat = cov_ex[["re_cov_mat"]],
#'   ws_cov_mat = cov_ex[["ws_cov_mat"]],
#'   re_names = c("b0", "b2"),
#'   print_intermediate = TRUE
#' )
#'
#' @export
simulate_subject_MVN_conditional <- function(
    xmat,
    betas,
    zmat,
    re_cov_mat,
    ws_cov_mat,
    re_names = NULL,
    print_intermediate = FALSE
) {
  # random effects covariance
  if (is.matrix(re_cov_mat) == TRUE) {
    b <- MASS::mvrnorm(
      n = 1,
      mu = rep(0, ncol(zmat)),
      Sigma = re_cov_mat
    )
  } else {
    b <- rnorm(1, 0, sd = re_cov_mat)
  }
  if (print_intermediate) {
    cat("Random effects = [", b, "] \n")
  }
  # within-subjects covariance
  e <- MASS::mvrnorm(
    n = 1,
    mu = rep(0, nrow(xmat)),
    Sigma = ws_cov_mat
  )
  if (print_intermediate) {
    cat("Within-subject errors = [", e, "] \n")
  }

  mt <- xmat %*% betas + zmat %*% b
  y <- mt + e

  b_mat <- matrix(
    b,
    byrow = TRUE,
    nrow = nrow(xmat),
    ncol = length(b)
  )
  if (is.null(re_names) == FALSE) {
    colnames(b_mat) <- re_names
  } else {
    colnames(b_mat) <- paste0("b", 0:(nrow(re_cov_mat) - 1))
  }

  df <- as.data.frame(
    cbind(
      y = as.vector(y),
      mt = as.vector(mt),
      xmat[, -1], # remove column of 1's from analysis data set
      b_mat,
      e
    )
  )

  return(df)
}
