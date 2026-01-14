#' Estimate MCSE for Competing Statistical Approaches
#'
#' Applies \link[SimRescueMeds]{mcse_estimands}, \link[SimRescueMeds]{mcse_coverage}, and \link[SimRescueMeds]{mcse_rejection}
#' to a list of data frames containing summaries of analysis results on the same set of simulated data sets.
#'
#' @param smry_list A list of data frames, each of which contains summaries from applying some statistical method
#' to the same set of simulated data sets.
#' @param methods_labels A character vector of labels for the different analysis approaches.
#' @param true_value Numeric. The true value of the estimand.
#' @param alpha Numeric. The desired significance level, which should be greater than 0 and less than 1.
#' @param type_rejection Character. \code{"pvalues"} if rejection data is a vector of p-values and \code{"rejection"}
#' if a binary vector.
#' @param rejection_name Character. The name of the column containing p-values or binary rejections.
#' @param remove_measure_labels Logical. If \code{TRUE}, removes column of measure labels for all but 1st data set.
#' This can make it easier to combine data frames into a single data frame for general summaries.
#' @param return_list Logical. When \code{TRUE}, return a list instead of a data frame.
#'
#' @inheritParams mcse_estimands
#' @inheritParams mcse_coverage
#'
#' @examples
#' set.seed(895642)
#' df1 <- data.frame(
#'   effect = rnorm(100),
#'   se_effect = rgamma(100, shape = 1)
#' )
#'
#' df1$lcl <- df1$effect - 1.96*df1$se_effect
#' df1$ucl <- df1$effect + 1.96*df1$se_effect
#' df1$reject <- ifelse(df1$lcl < 0 & df1$ucl > 0, 1, 0)
#'
#' df2 <- data.frame(
#'   effect = rnorm(100),
#'   se_effect = rgamma(100, shape = 1)
#' )
#'
#' df2$lcl <- df2$effect - 1.96*df2$se_effect
#' df2$ucl <- df2$effect + 1.96*df2$se_effect
#' df2$reject <- ifelse(df2$lcl < 0 & df2$ucl > 0, 1, 0)
#'
#' mcse_ex <- mcse_many(
#'   smry_list = list(df1 = df1, df2 = df2),
#'   methods_labels = c("Spiderman", "Batman"),
#'   true_value = 0,
#'   type_rejection = "rejection",
#'   rejection_name = "reject"
#' )
#'
#' colnames(mcse_ex) <- c("Measure",
#'                        "Spiderman.estimate",
#'                        "Spiderman.mcse",
#'                        "Batman.estimate",
#'                        "Batman.mcse")
#'
#' mcse_ex
#'
#' @references
#' \insertRef{Morris2019}{SimRescueMeds}
#'
#' @export
mcse_many <- function(
    smry_list,
    methods_labels = c("CFB", "ANCOVA", "LMM", "JM", "JMbayes2"),
    true_value = -1,
    alpha = 0.05,
    estimand_name = "effect",
    se_name = "se_effect",
    lcl_name = "lcl",
    ucl_name = "ucl",
    rejection_name = "p",
    type_rejection = "pvalues",
    remove_measure_labels = TRUE,
    return_list = FALSE
) {
  mcse_list <- vector(mode = "list", length = length(smry_list))
  for (i in seq_len(length(smry_list))) {
    if (i == 1 | remove_measure_labels == FALSE) {
      mcse_list[[i]] <- rbind(
        mcse_estimands(data = smry_list[[i]],
                       true_value = true_value,
                       estimand_name = estimand_name,
                       se_name = se_name),
        mcse_coverage(data = smry_list[[i]],
                      true_value = true_value,
                      lcl_name = lcl_name,
                      ucl_name = ucl_name),
        mcse_rejection(data = smry_list[[i]][[rejection_name]],
                       type = type_rejection,
                       alpha = alpha)
      )
    } else {
      mcse_list[[i]] <- rbind(
        mcse_estimands(data = smry_list[[i]],
                       true_value = true_value,
                       estimand_name = estimand_name,
                       se_name = se_name),
        mcse_coverage(data = smry_list[[i]],
                      true_value = true_value,
                      lcl_name = lcl_name,
                      ucl_name = ucl_name),
        mcse_rejection(data = smry_list[[i]][[rejection_name]],
                       type = type_rejection,
                       alpha = alpha)
      )[, -1]
    }
  }
  names(mcse_list) <- methods_labels
  if (return_list == FALSE) {
    mcse_df <- do.call(
      "cbind",
      mcse_list
    )
    return(mcse_df)
  } else {
    return(mcse_list)
  }
}
