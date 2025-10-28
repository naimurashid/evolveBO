
#' Global and local sensitivity via emulator
#' @export
sa_sobol <- function(surr, bounds, outcome="EN", n_mc=1000) {
  message("Sobol indices placeholder.")
  data.frame(parameter=names(bounds), S_first=runif(length(bounds)))
}

#' Local gradients of emulator mean
#' @export
sa_gradients <- function(surr, theta, outcome="EN", eps=1e-4) {
  message("Gradient sensitivity placeholder.")
  data.frame(parameter=names(theta), grad=rnorm(length(theta)))
}

#' Covariance of effects
#' @export
cov_effects <- function(surr, bounds, outcome="EN", n_mc=500, eps=1e-4) {
  message("Covariance structure placeholder.")
  matrix(rnorm(length(bounds)^2), ncol=length(bounds))
}
