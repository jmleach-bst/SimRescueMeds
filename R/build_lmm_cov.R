#' Build Linear Mixed Model Covariance Matrix
#'
#' Constructs the covariance matrix for linear mixed model based on random effects
#' and within-subject covariance specifications.
#'
#' @param zmat Random effects design matrix.
#' @param re_corr_mat Random effects correlation matrix (1 if random intercept only).
#' @param re_sigma Random effects standard deviation vector.
#' @param ws_corr_mat Within-subject correlation matrix, usually identity.
#' @param ws_sigma Within-subject standard deviation vector, usually a single scalar.
#' @param partition Logical. when \code{TRUE} returns list containing each covariance matrix.
#' @param print_intermediate Logical. When \code{TRUE} prints intermediate output.
#'
#' @return A matrix or list of matrices.
#'
#' @importFrom matrixcalc is.square.matrix
#' @importFrom matrixcalc is.positive.definite
#'
#' @details
#' This function builds a subject-specific covariance matrix as
#' \deqn{\mathbf{\Sigma}_i = \mathbf{Z}_i\mathbf{\Sigma}_b\mathbf{Z}_i^\top + \mathbf{\Sigma}_e}
#' where \eqn{\mathbf{Z}_i} is the subject-specific random effects matrix, \eqn{\mathbf{\Sigma}_b} is the
#' random-effects covariance, i.e., constructed from \code{re_corr_mat} and \code{re_sigma},
#' and \eqn{\mathbf{\Sigma}_e} is the within-subject covariance matrix, which is constructed from
#' \code{ws_corr_mat} and \code{ws_sigma}. If \code{partition = FALSE}, then only the final matrix
#' \eqn{\mathbf{Z}_i} is returned, but if \code{partition = TRUE}, then \eqn{\mathbf{Z}_b}
#' and \eqn{\mathbf{Z}_e} are returned separately in a list.
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
#' print(Sigma_i)
#'
#' Sigma_partition <- build_lmm_cov(
#'   zmat = xmat_trt[, c("intercept", "time")],
#'   re_sigma = c(1.15, 1.05),
#'   re_corr_mat = matrix(
#'     c(1, 0,
#'       0, 1),
#'     byrow = TRUE,
#'     nrow = 2
#'   ),
#'   ws_sigma = 1.25,
#'   ws_corr_mat = diag(5),
#'   partition = TRUE,
#'   print_intermediate = FALSE
#' )
#'
#' print(Sigma_partition)
#'
#' @export
build_lmm_cov <- function(
    zmat, # random effects design matrix
    re_corr_mat = 1, # random effect correlation matrix (1 if random intercept)
    re_sigma, # random effects SD vector
    ws_corr_mat, # within-subject correlation matrix, usually identity
    ws_sigma, # within-subject SD vector, usually a single scalar
    partition = FALSE, # when TRUE returns list containing each covariance matrix
    print_intermediate = FALSE
) {
  # random effects covariance
  if (is.matrix(re_corr_mat) == TRUE) {
    if (matrixcalc::is.square.matrix(re_corr_mat) == FALSE) {
      stop("Random effects correlation matrix is not square")
    }
  }
  if (is.matrix(re_corr_mat) == TRUE) {
    if (dim(re_corr_mat)[1] != length(re_sigma)) {
      stop("Random effects variances and correlation matrix not compatible (check dimensions)")
    }
  }
  if (length(re_sigma) == 1) {
    re_cov_mat <- re_sigma^2
    zSzt <- (zmat %*% t(zmat)) * re_sigma^2
  } else {
    re_cov_mat <- diag(re_sigma) %*% re_corr_mat %*% diag(re_sigma)
    if (matrixcalc::is.positive.definite(re_cov_mat) == FALSE) {
      stop("Random effects covariance matrix is not positive definite")
    }
    zSzt <- zmat %*% re_cov_mat %*% t(zmat)
  }
  if (print_intermediate == TRUE) {
    cat("Random effects covariance matrix \n")
    print(re_cov_mat)
    cat("Z %*% Sigma_b %*% t(Z) \n")
    print(zSzt)
  }

  # within-subjects covariance
  if (is.matrix(ws_corr_mat) == TRUE) {
    if(matrixcalc::is.square.matrix(ws_corr_mat) == FALSE) {
      stop("Within subjects correlation matrix is not square")
    }
  }
  if (length(ws_sigma) == 1) {
    ws_cov_mat <- ws_corr_mat * ws_sigma^2
  } else {
    ws_cov_mat <- diag(ws_sigma) %*% ws_corr_mat %*% diag(ws_sigma)
  }
  if (matrixcalc::is.positive.definite(ws_cov_mat) == FALSE) {
    stop("Within subjects covariance matrix is not positive definite")
  }
  if (print_intermediate == TRUE) {
    cat("Within-subjects covariance matrix \n")
    print(ws_cov_mat)
  }

  if (partition == FALSE) {
    # combined covariance
    cov_mat <- zSzt + ws_cov_mat
    return(cov_mat)
  } else {
    return(list(re_cov_mat = re_cov_mat,
                ws_cov_mat = ws_cov_mat))
  }
}
