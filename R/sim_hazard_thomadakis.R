#' Simulate Time-to-Rescue-Therapy for Single Subject
#'
#' Use proportional hazards models to simulate time-to-rescue therapy
#' that is dependent on longitudinal trajectories.
#'
#' @param y A vector of longitudinal outcome measurements for a single subject.
#' @param time A vector of time points corresponding to each longitudinal measurement.
#' @param theta1 Numeric. An association parameter for \eqn{y_{i, j}}.
#' @param theta2 Numeric. An association parameter for \eqn{y_{i,j+1}}. If \code{theta2} is
#' non-zero, then this implies at least some missingness is not-at-random (MNAR).
#' @param dist Character. Either \code{"exponential"} or \code{"weibull"} distributions for the
#' survival function.
#' @param lambda Numeric. The rate parameter for \link[stats]{rexp} or scale for \link[stats]{rweibull}.
#' @param gamma Numeric. The shape parameter for \link[stats]{rweibull}.
#' @param maxt Numeric. The last possible follow-up time, after which observations
#' are censored. In the present context, the time the study ends.
#'
#' @importFrom stats rexp
#' @importFrom stats rweibull
#'
#' @return A data frame.
#'
#' @importFrom Rdpack reprompt
#'
#' @details
#' This function simulates time-to-rescue medication based on the approach of Thomadakis (2019).
#' A general formulation for simulating from missing-at-random (MAR) and missing-not-at-random
#' (MNAR) data is to allow the hazard of rescue therapy within some interval \eqn{t_{i,j} \le t < t_{i,j+1}}
#' to be a function of the last observed outcome measurement, \eqn{y_{i,j}}, and the next measurement,
#' \eqn{y_{i,j+1}}, the latter of which is unobserved in the study since it occurs after rescue therapy.
#' If the event-time distribution is exponential, then we have the following hazard function, but it is
#' straightforward to extend to Weibull, for which exponential is a special case.
#'
#' \deqn{
#' h_i(t) = \lambda\exp\left(y_{i,j}\theta_1 + y_{i,j+1}\theta_2\right), \quad t_{i,j} \le t < t_{i,j+1}.
#' }
#'
#' This requires simulating times-to-rescue therapy within each interval between measurements until either
#' one of those intervals produces rescue therapy or the study ends. If \eqn{\theta_2 \ne 0}, then there is
#' MNAR missingness. If \eqn{\theta_1 = \theta_2 = 0}, then MCAR. If \eqn{\theta_1 \ne 0} and \eqn{\theta_2 = 0},
#' then MAR.
#'
#' @examples
#' set.seed(65378)
#' sim_hazard_thomadakis(
#'   y = rnorm(n = 5),
#'   time = c(0, 3, 6, 9, 12)/12
#'   )
#' sim_hazard_thomadakis(
#'   y = rnorm(n = 5),
#'   time = c(0, 3, 6, 9, 12)/12
#'   )
#' set.seed(15378)
#' sim_hazard_thomadakis(
#'   y = rnorm(n = 5),
#'   time = c(0, 3, 6, 9, 12)/12
#'   )
#'
#' set.seed(28394)
#' sim_hazard_thomadakis(
#'   y = rnorm(n = 5),
#'   time = c(0, 3, 6, 9, 12)/12,
#'   dist = "weibull",
#'   gamma = 2
#'   )
#'
#' set.seed(28235)
#' sim_hazard_thomadakis(
#'   y = rnorm(n = 5),
#'   time = c(0, 3, 6, 9, 12)/12,
#'   dist = "weibull",
#'   gamma = 2
#'   )
#'
#' @references
#' \insertRef{Thomadakis2019}{SimRescueMeds}
#'
#' @export
sim_hazard_thomadakis <- function(
    y,
    time,
    theta1 = 0,
    theta2 = 0,
    dist = "exponential",
    lambda = 1,
    gamma = 1,
    maxt = 1
) {
  if (lambda <= 0 | gamma <= 0) {
    stop("Require lambda > 0 and gamma > 0")
  }
  if (length(y) != length(time)) {
    stop("Require observation and time vectors to have equal length")
  }
  if (is.numeric(y) == FALSE | is.numeric(time) == FALSE) {
    stop("Require observation and time vectors to be numeric")
  }
  ni <- length(y)
  tte <- c()
  event <- c()
  eventtime <- c()
  for (i in seq_len(ni - 1)) {
    if (dist == "exponential") {
      tte[i] <- time[i] + rexp(n = 1,
                               rate = lambda*exp(y[i]*theta1 + y[i + 1]*theta2))
    } else {
      tte[i] <- time[i] + rweibull(n = 1,
                                   shape = gamma,
                                   scale = (1/(lambda*exp(y[i]*theta1 + y[i + 1]*theta2)))^(1/gamma))
    }

    if (tte[i] < time[i + 1]) {
      event[i] <- 1
      eventtime[i] <- tte[i]
      break
    } else {
      event[i] <- 0
      eventtime[i] <- maxt
    }
  }
  return(cbind(tte, eventtime, event))
}
