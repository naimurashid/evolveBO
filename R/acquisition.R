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
  mu_obj <- obj_pred$mean
  sd_obj <- obj_pred$sd

  metric_names <- names(pred)
  prob_feas <- purrr::map_dbl(seq_along(mu_obj), function(i) {
    mean_vec <- purrr::map_dbl(pred, ~ .x$mean[[i]]) |>
      purrr::set_names(metric_names)
    sd_vec <- purrr::map_dbl(pred, ~ .x$sd[[i]]) |>
      purrr::set_names(metric_names)
    prob_feasibility(mean_vec, sd_vec, constraint_tbl)
  })

  has_feasible <- is.finite(best_feasible)

  if (!has_feasible) {
    # No feasible solution yet
    # Strategy: Minimize expected constraint violation weighted by feasibility probability
    violations <- compute_expected_violation(pred, constraint_tbl, metric_names)

    # Normalize violations to [0, 1] range
    max_viol <- max(violations, na.rm = TRUE)
    if (is.finite(max_viol) && max_viol > 0) {
      normalized_viol <- violations / max_viol
    } else {
      normalized_viol <- violations
    }

    # Acquisition = (1 - violation) + exploration bonus, weighted by feasibility
    # Scale to be comparable to EI values
    acq <- (1 - normalized_viol + 0.3 * sd_obj) * (0.3 + 0.7 * prob_feas)
    return(acq)
  }

  # Standard ECI: EI × P(feasible)
  ei <- compute_ei(mu_obj, sd_obj, best_feasible)
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
    # No feasible solution yet - use exploration weighted by uncertainty
    # Add small epsilon to prevent division by zero
    return(pmax(sd + 1e-10, 0))
  }
  improvement <- best_feasible - mu
  positive <- sd > 1e-10  # Add epsilon for numerical stability
  z <- improvement[positive] / (sd[positive] + 1e-10)
  phi <- stats::dnorm(z)
  Phi <- stats::pnorm(z)
  ei[positive] <- improvement[positive] * Phi + (sd[positive] + 1e-10) * phi
  ei[!positive] <- pmax(improvement[!positive], 0)
  ei
}

#' Compute expected constraint violation
#'
#' For each candidate point, computes the expected magnitude of constraint
#' violations under the GP posterior. Used to guide search toward feasible
#' region when no feasible solutions have been found yet.
#'
#' @param pred list of predictions from predict_surrogates (contains mean and sd for each metric)
#' @param constraint_tbl tibble from parse_constraints
#' @param metric_names character vector of metric names
#'
#' @return numeric vector of expected violations for each candidate
#' @keywords internal
compute_expected_violation <- function(pred, constraint_tbl, metric_names) {
  n_candidates <- length(pred[[1]]$mean)
  violations <- numeric(n_candidates)

  for (i in seq_len(n_candidates)) {
    total_viol <- 0

    for (j in seq_len(nrow(constraint_tbl))) {
      metric <- constraint_tbl$metric[j]
      direction <- constraint_tbl$direction[j]
      threshold <- constraint_tbl$threshold[j]

      if (!metric %in% metric_names) next

      mu_val <- pred[[metric]]$mean[[i]]
      sd_val <- pred[[metric]]$sd[[i]]

      if (direction == "ge") {
        # Constraint: metric >= threshold
        # Violation when metric < threshold
        # E[max(0, threshold - metric)]
        z <- (threshold - mu_val) / (sd_val + 1e-10)
        prob_violate <- stats::pnorm(z)
        # Expected shortfall (expected violation magnitude)
        expected_shortfall <- (sd_val + 1e-10) * stats::dnorm(z) / (prob_violate + 1e-10)
        total_viol <- total_viol + prob_violate * expected_shortfall

      } else {  # "le"
        # Constraint: metric <= threshold
        # Violation when metric > threshold
        # E[max(0, metric - threshold)]
        z <- (mu_val - threshold) / (sd_val + 1e-10)
        prob_violate <- stats::pnorm(z)
        # Expected excess (expected violation magnitude)
        expected_excess <- (sd_val + 1e-10) * stats::dnorm(z) / (prob_violate + 1e-10)
        total_viol <- total_viol + prob_violate * expected_excess
      }
    }

    violations[i] <- total_viol
  }

  violations
}

#' Select diverse batch using local penalization
#'
#' Implements the local penalization strategy of González et al. (2016).
#' Iteratively selects points by penalizing acquisition near previously
#' selected points, ensuring spatial diversity in batch selection.
#'
#' @param candidates list of candidate points (unit scale)
#' @param acq_scores numeric vector of acquisition values for each candidate
#' @param q batch size (number of points to select)
#' @param lipschitz Lipschitz constant for penalization. Higher values enforce
#'   greater diversity. Default estimated from typical BO landscapes.
#'
#' @return integer vector of indices of selected candidates
#' @keywords internal
#'
#' @references
#' González, J., Dai, Z., Hennig, P., & Lawrence, N. (2016).
#' Batch Bayesian Optimization via Local Penalization. AISTATS.
select_batch_local_penalization <- function(candidates, acq_scores, q,
                                             lipschitz = 10) {
  n_candidates <- length(acq_scores)
  if (q >= n_candidates) {
    return(seq_len(n_candidates))
  }

  selected_indices <- integer(q)
  penalized_scores <- acq_scores
  candidates_matrix <- do.call(rbind, candidates)

  for (i in seq_len(q)) {
    # Select point with highest penalized acquisition
    best_idx <- which.max(penalized_scores)
    selected_indices[i] <- best_idx

    if (i < q) {
      # Penalize acquisition near selected point
      selected_point <- candidates_matrix[best_idx, , drop = FALSE]

      # Compute distances to selected point
      distances <- compute_distances(candidates_matrix, selected_point)

      # Penalization function: max(0, L * r - acq_best)
      # where r is distance, L is Lipschitz constant
      # Points within distance acq_best/L of selected point get penalized
      penalty <- pmax(0, lipschitz * distances - penalized_scores[best_idx])

      # Apply penalty
      penalized_scores <- penalized_scores - penalty

      # Ensure we don't select same point again
      penalized_scores[best_idx] <- -Inf
    }
  }

  selected_indices
}

#' Compute Euclidean distances from points to reference
#'
#' @param points n × d matrix of points
#' @param reference 1 × d matrix (or vector) of reference point
#'
#' @return numeric vector of n distances
#' @keywords internal
compute_distances <- function(points, reference) {
  # Ensure reference is a matrix
  if (is.vector(reference)) {
    reference <- matrix(reference, nrow = 1)
  }

  # Compute differences
  diff <- sweep(points, 2, reference[1, ], "-")

  # Euclidean distance
  sqrt(rowSums(diff^2))
}

#' Estimate Lipschitz constant from GP lengthscales
#'
#' Uses GP lengthscales to estimate a reasonable Lipschitz constant
#' for local penalization. Conservative estimate to ensure diversity.
#'
#' @param surrogates list of GP models
#' @param objective name of objective metric
#'
#' @return numeric scalar, Lipschitz constant estimate
#' @keywords internal
estimate_lipschitz <- function(surrogates, objective) {
  model <- surrogates[[objective]]

  # Try to extract lengthscales from DiceKriging model
  if (inherits(model, "km")) {
    tryCatch({
      lengthscales <- model@covariance@range.val
      # Lipschitz constant inversely related to lengthscale
      # Use minimum lengthscale (most sensitive direction)
      L <- 1 / min(lengthscales)
      # Conservative factor: multiply by 2
      return(L * 2)
    }, error = function(e) {
      # If extraction fails, use default
      return(10)
    })
  }

  # Default for other model types
  10
}
