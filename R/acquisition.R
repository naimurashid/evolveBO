#' Expected constrained improvement acquisition function
#'
#' @param unit_x list of candidate points on the unit hypercube.
#' @param surrogates named list of fitted surrogate models (as returned by [fit_surrogates()]).
#' @param constraint_tbl tibble from [parse_constraints()].
#' @param objective name of the objective metric (character scalar).
#' @param best_feasible best observed value of the objective among feasible points
#'   (numeric scalar or `Inf` if none observed).
#'
#' @return numeric vector of acquisition scores for each candidate.
#' @export
acq_eci <- function(unit_x,
                    surrogates,
                    constraint_tbl,
                    objective,
                    best_feasible) {
  if (!objective %in% names(surrogates)) {
    stop("Objective surrogate not available.", call. = FALSE)
  }
  pred <- predict_surrogates(surrogates, unit_x)

  obj_pred <- pred[[objective]]
  mu <- obj_pred$mean
  sd <- obj_pred$sd

  metric_names <- names(pred)
  prob_feas <- purrr::map_dbl(seq_along(mu), function(i) {
    mean_vec <- purrr::map_dbl(pred, ~ .x$mean[[i]]) |>
      purrr::set_names(metric_names)
    sd_vec <- purrr::map_dbl(pred, ~ .x$sd[[i]]) |>
      purrr::set_names(metric_names)
    prob_feasibility(mean_vec, sd_vec, constraint_tbl)
  })

  ei <- compute_ei(mu, sd, best_feasible)
  ei * prob_feas
}

#' Quasi Expected Hypervolume Improvement (constraint-aware)
#'
#' @description
#' \strong{Note:} This function is currently a placeholder that calls
#' [acq_eci()]. A full qEHVI implementation with proper hypervolume computation
#' is planned for a future release.
#'
#' For batch size greater than one, the function uses a sequential (Kriging
#' believer) approximation via constrained expected improvement.
#'
#' @inheritParams acq_eci
#' @return numeric vector of acquisition scores.
#' @keywords internal
acq_qehvi <- function(unit_x,
                      surrogates,
                      constraint_tbl,
                      objective,
                      best_feasible) {
  warning("acq_qehvi is currently an alias for acq_eci. Full qEHVI implementation is planned for a future release.",
          call. = FALSE)
  acq_eci(unit_x, surrogates, constraint_tbl, objective, best_feasible)
}

#' @keywords internal
compute_ei <- function(mu, sd, best_feasible) {
  ei <- numeric(length(mu))
  has_feasible <- is.finite(best_feasible)
  if (!has_feasible) {
    return(pmax(sd, 0))
  }
  improvement <- best_feasible - mu
  positive <- sd > 0
  z <- improvement[positive] / sd[positive]
  phi <- stats::dnorm(z)
  Phi <- stats::pnorm(z)
  ei[positive] <- improvement[positive] * Phi + sd[positive] * phi
  ei[!positive] <- pmax(improvement[!positive], 0)
  ei
}
