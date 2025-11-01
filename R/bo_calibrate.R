#' Calibrate adaptive trial designs with constrained Bayesian optimisation
#'
#' Implements the methodology described in Section 2 of the manuscript:
#' heteroskedastic Gaussian process emulation, expected constrained improvement,
#' and multi-fidelity simulation budgeting. The function orchestrates initial
#' design selection, surrogate fitting, sequential acquisition, and diagnostic
#' bookkeeping required for downstream benchmarking and reporting.
#'
#' @param sim_fun callback with signature
#'   `function(theta, fidelity = c("low","med","high"), ...)` returning a named
#'   numeric vector of operating characteristics (e.g., `power`, `type1`, `EN`,
#'   `ET`). The function should attach attributes `variance` (named numeric
#'   vector of Monte Carlo variances) and `n_rep` (number of simulator
#'   replications) for best results. If `variance` is not provided, a binomial
#'   approximation is used for metrics in [0,1], and a small nugget is used for
#'   other metrics.
#' @param bounds named list of length-two numeric vectors specifying lower and
#'   upper limits for each design parameter.
#' @param objective name of the operating characteristic to minimise (character).
#' @param constraints named list mapping metric -> `c(direction, threshold)` with
#'   `direction` in `{"ge","le"}`.
#' @param n_init number of initial design evaluations generated via Latin
#'   hypercube sampling.
#' @param q batch size for each BO iteration (number of new evaluations per
#'   acquisition round).
#' @param budget total number of simulator evaluations (initial + BO iterations).
#' @param seed base RNG seed for reproducibility.
#' @param fidelity_levels named numeric vector giving the number of simulator
#'   replications associated with each fidelity level (defaults to low/med/high).
#' @param acquisition acquisition rule to use (currently only `"eci"` is
#'   implemented).
#' @param candidate_pool number of random candidate points assessed per
#'   acquisition step.
#' @param covtype covariance kernel passed to [fit_surrogates()].
#' @param integer_params optional character vector of parameters that should be
#'   rounded to the nearest integer prior to simulation.
#' @param progress logical; if `TRUE` (default) messages are emitted at key
#'   milestones.
#' @param ... additional arguments forwarded to `sim_fun`.
#'
#' @return An object of class `evolveBO_fit` containing the optimisation history,
#'   best design, fitted surrogates, policy configuration, and posterior draws
#'   supporting sensitivity diagnostics.
#' @export
bo_calibrate <- function(sim_fun,
                         bounds,
                         objective,
                         constraints,
                         n_init = 40,
                         q = 8,
                         budget = 150,
                         seed = 2025,
                         fidelity_levels = c(low = 200, med = 1000, high = 10000),
                         acquisition = c("eci"),
                         candidate_pool = 2000,
                         covtype = "matern5_2",
                         integer_params = NULL,
                         progress = TRUE,
                         ...) {
  acquisition <- match.arg(acquisition)
  validate_bounds(bounds)
  constraint_tbl <- parse_constraints(constraints)

  # Validate numeric parameters
  if (n_init <= 0) {
    stop("`n_init` must be positive.", call. = FALSE)
  }
  if (q <= 0) {
    stop("`q` (batch size) must be positive.", call. = FALSE)
  }
  if (budget <= 0) {
    stop("`budget` must be positive.", call. = FALSE)
  }
  if (candidate_pool <= 0) {
    stop("`candidate_pool` must be positive.", call. = FALSE)
  }
  if (budget < n_init) {
    warning("`budget` is less than `n_init`; only initial design will be evaluated.",
            call. = FALSE)
  }

  rng_seed <- if (is.null(seed)) sample.int(1e6, 1) else seed
  set.seed(rng_seed)

  fidelity_levels <- validate_fidelity_levels(fidelity_levels)
  primary_fidelity <- names(fidelity_levels)[1]

  if (progress) {
    message("Initialising evolveBO calibration with seed = ", rng_seed)
  }

  history <- initialise_history()

  n_init <- min(n_init, budget)
  initial_design <- lhs_design(n_init, bounds, seed = rng_seed)

  eval_counter <- 0L
  iter_counter <- 0L

  for (theta in initial_design) {
    theta <- coerce_theta_types(theta, integer_params)
    eval_counter <- eval_counter + 1L
    history <- record_evaluation(
      history = history,
      sim_fun = sim_fun,
      theta = theta,
      bounds = bounds,
      objective = objective,
      constraint_tbl = constraint_tbl,
      fidelity = primary_fidelity,
      fidelity_levels = fidelity_levels,
      eval_id = eval_counter,
      iter = 0L,
      seed = rng_seed + eval_counter,
      progress = progress,
      ...
    )
  }

  while (eval_counter < budget) {
    iter_counter <- iter_counter + 1L
    if (progress) {
      message(sprintf("Iteration %d: fitting surrogates on %d evaluations.",
                      iter_counter, nrow(history)))
    }

    surrogates <- fit_surrogates(history, objective, constraint_tbl, covtype = covtype)
    best_feasible_value <- best_feasible_objective(history, objective)
    candidate_pool_size <- max(candidate_pool, 50 * length(bounds))
    unit_candidates <- lhs_candidate_pool(candidate_pool_size, bounds)

    acquisition_scores <- evaluate_acquisition(
      acquisition = acquisition,
      unit_candidates = unit_candidates,
      surrogates = surrogates,
      constraint_tbl = constraint_tbl,
      objective = objective,
      best_feasible = best_feasible_value
    )

    order_idx <- order(acquisition_scores, decreasing = TRUE)
    n_new <- min(q, budget - eval_counter)
    selected_idx <- order_idx[seq_len(n_new)]
    selected_candidates <- unit_candidates[selected_idx]

    for (i in seq_len(n_new)) {
      eval_counter <- eval_counter + 1L
      chosen_unit <- selected_candidates[[i]]
      theta <- scale_from_unit(chosen_unit, bounds)
      theta <- coerce_theta_types(theta, integer_params)

      prob_feas <- estimate_candidate_feasibility(
        surrogates = surrogates,
        unit_x = list(chosen_unit),
        constraint_tbl = constraint_tbl
      )

      # Get CV (coefficient of variation) from surrogate for objective
      pred_obj <- predict_surrogates(surrogates, list(chosen_unit))[[objective]]
      cv_estimate <- pred_obj$sd[[1]] / max(abs(pred_obj$mean[[1]]), 1e-6)

      fidelity <- select_fidelity_staged(prob_feas, cv_estimate,
                                         iter_counter, fidelity_levels)

      history <- record_evaluation(
        history = history,
        sim_fun = sim_fun,
        theta = theta,
        bounds = bounds,
        objective = objective,
        constraint_tbl = constraint_tbl,
        fidelity = fidelity,
        fidelity_levels = fidelity_levels,
        eval_id = eval_counter,
        iter = iter_counter,
        seed = rng_seed + eval_counter,
        progress = progress,
        ...
      )
    }
  }

  final_surrogates <- fit_surrogates(history, objective, constraint_tbl, covtype = covtype)
  best_idx <- best_feasible_index(history, objective)
  best_theta <- history$theta[[best_idx]]

  diagnostics <- generate_diagnostics(final_surrogates, history, objective)

  structure(
    list(
      history = history,
      best_theta = best_theta,
      surrogates = final_surrogates,
      policies = list(
        acquisition = acquisition,
        q = q,
        budget = budget,
        seed = rng_seed,
        fidelity_levels = fidelity_levels,
        candidate_pool = candidate_pool,
        objective = objective
      ),
      diagnostics = diagnostics
    ),
    class = c("evolveBO_fit", "list")
  )
}

#' @keywords internal
initialise_history <- function() {
  tibble::tibble(
    iter = integer(),
    eval_id = integer(),
    theta = list(),
    unit_x = list(),
    theta_id = character(),
    fidelity = character(),
    n_rep = integer(),
    metrics = list(),
    variance = list(),
    objective = numeric(),
    feasible = logical()
  )
}

#' @keywords internal
record_evaluation <- function(history,
                              sim_fun,
                              theta,
                              bounds,
                              objective,
                              constraint_tbl,
                              fidelity,
                              fidelity_levels,
                              eval_id,
                              iter,
                              seed,
                              progress,
                              ...) {
  unit_theta <- scale_to_unit(theta, bounds)
  theta_id <- theta_to_id(unit_theta)
  result <- invoke_simulator(
    sim_fun = sim_fun,
    theta = theta,
    fidelity = fidelity,
    n_rep = fidelity_levels[[fidelity]],
    seed = seed,
    ...
  )
  metrics <- result$metrics
  variance <- result$variance
  n_rep <- result$n_rep
  feasible <- is_feasible(metrics, constraint_tbl)
  objective_value <- metrics[[objective]]
  if (is.null(objective_value)) {
    stop(sprintf("Objective '%s' not returned by simulator.", objective), call. = FALSE)
  }

  if (progress) {
    msg <- sprintf("  Eval %03d (iter %02d, %s fidelity): %s = %.4f%s",
                   eval_id, iter, fidelity, objective, objective_value,
                   if (feasible) " [feasible]" else "")
    message(msg)
  }

  dplyr::bind_rows(
    history,
    tibble::tibble(
      iter = iter,
      eval_id = eval_id,
      theta = list(theta),
      unit_x = list(unit_theta),
      theta_id = theta_id,
      fidelity = fidelity,
      n_rep = n_rep,
      metrics = list(metrics),
      variance = list(variance),
      objective = objective_value,
      feasible = feasible
    )
  )
}

#' @keywords internal
invoke_simulator <- function(sim_fun, theta, fidelity, n_rep, seed, ...) {
  res <- sim_fun(theta, fidelity = fidelity, seed = seed, ...)
  variance <- attr(res, "variance", exact = TRUE)
  rep_attr <- attr(res, "n_rep", exact = TRUE)
  if (is.list(res) && !is.null(res$metrics)) {
    metrics <- res$metrics
    if (is.null(variance) && !is.null(res$variance)) {
      variance <- res$variance
    }
    if (is.null(rep_attr) && !is.null(res$n_rep)) {
      rep_attr <- res$n_rep
    }
  } else {
    metrics <- res
  }
  metrics <- metrics |> as.list() |> purrr::imap_dbl(function(value, name) as.numeric(value))
  if (is.null(variance)) {
    variance <- default_variance_estimator(metrics, n_rep)
  } else {
    variance <- variance |> as.list() |> purrr::imap_dbl(function(value, name) as.numeric(value))
  }
  if (is.null(rep_attr) || is.na(rep_attr)) {
    rep_attr <- n_rep
  }
  list(
    metrics = metrics,
    variance = variance,
    n_rep = rep_attr
  )
}

#' Default variance estimator for simulator outputs
#'
#' Provides a fallback variance estimate when the simulator does not return
#' variance information. Uses binomial variance formula for proportion-type
#' metrics (values in [0,1]) and returns NA for other metrics.
#'
#' @param metrics named numeric vector of operating characteristics.
#' @param n_rep number of simulator replications.
#'
#' @details
#' For metrics that appear to be proportions (e.g., power, type I error),
#' this function estimates variance using the binomial formula: p(1-p)/n.
#' For continuous metrics outside [0,1] (e.g., expected sample size, trial
#' duration), it returns NA, which triggers the use of a small homoskedastic
#' nugget in the Gaussian process surrogate.
#'
#' \strong{Recommendation:} For best results, simulators should return variance
#' estimates as an attribute, especially for non-proportion metrics where the
#' binomial approximation does not apply.
#'
#' @return Named numeric vector of variance estimates (may contain NA values).
#' @keywords internal
default_variance_estimator <- function(metrics, n_rep) {
  if (is.null(n_rep) || is.na(n_rep) || n_rep <= 0) {
    return(rep(NA_real_, length(metrics)) |> purrr::set_names(names(metrics)))
  }
  purrr::imap_dbl(metrics, function(value, name) {
    # Use binomial variance for proportion-type metrics (values in [0,1])
    if (!is.na(value) && value >= 0 && value <= 1) {
      max(value * (1 - value) / n_rep, 1e-6)
    } else {
      # For continuous metrics outside [0,1], return NA
      # The GP fitting will use a small nugget for these cases
      NA_real_
    }
  })
}

#' @keywords internal
validate_fidelity_levels <- function(fidelity_levels) {
  if (is.null(names(fidelity_levels))) {
    stop("`fidelity_levels` must be a named numeric vector.", call. = FALSE)
  }
  purrr::iwalk(fidelity_levels, function(value, name) {
    if (!is.finite(value) || value <= 0) {
      stop("Fidelity replication counts must be positive.", call. = FALSE)
    }
  })
  fidelity_levels
}

#' @keywords internal
lhs_candidate_pool <- function(n, bounds) {
  lhs <- lhs::randomLHS(n, length(bounds))
  purrr::map(seq_len(n), function(i) {
    values <- lhs[i, , drop = TRUE]
    names(values) <- names(bounds)
    values
  })
}

#' @keywords internal
best_feasible_objective <- function(history, objective) {
  idx <- which(history$feasible)
  if (length(idx) == 0L) {
    return(Inf)
  }
  min(history$objective[idx], na.rm = TRUE)
}

#' @keywords internal
best_feasible_index <- function(history, objective) {
  feasible_idx <- which(history$feasible)
  if (length(feasible_idx) == 0L) {
    warning("No feasible solutions found; returning best infeasible point.",
            call. = FALSE)
    return(which.min(history$objective))
  }
  feasible_idx[which.min(history$objective[feasible_idx])]
}

#' @keywords internal
estimate_candidate_feasibility <- function(surrogates, unit_x, constraint_tbl) {
  pred <- predict_surrogates(surrogates, unit_x)
  names <- names(pred)
  mean_vec <- purrr::map_dbl(pred, ~ .x$mean[[1]]) |> purrr::set_names(names)
  sd_vec <- purrr::map_dbl(pred, ~ .x$sd[[1]]) |> purrr::set_names(names)
  prob_feasibility(mean_vec, sd_vec, constraint_tbl)
}

#' Select fidelity using staged multi-fidelity strategy
#'
#' Implements the staged approach described in manuscript Section 2.4:
#' - Iterations 1-30: uniform low fidelity for exploration
#' - Iterations 31-100: adaptive selection based on feasibility and CV
#' - Iterations 101+: high fidelity near optimum and boundaries
#'
#' @param prob_feasible probability of constraint satisfaction
#' @param cv_estimate coefficient of variation from surrogate
#' @param iter current iteration number
#' @param fidelity_levels named vector of fidelity levels
#' @keywords internal
select_fidelity_staged <- function(prob_feasible, cv_estimate, iter, fidelity_levels) {
  if (length(fidelity_levels) == 1L) {
    return(names(fidelity_levels))
  }

  # Stage 1: Initial exploration (iterations 1-30) - uniform low fidelity
  if (iter <= 30) {
    return(names(fidelity_levels)[1])  # "low"
  }

  # Stage 3: Final refinement (iterations 101+) - high fidelity near optimum
  if (iter > 100 && prob_feasible > 0.6 && "high" %in% names(fidelity_levels)) {
    return("high")
  }

  # Stage 2: Focused search (iterations 31-100) - adaptive via MFKG heuristic
  # Promote fidelity if:
  # 1. High uncertainty (CV > 0.18), AND
  # 2. Near feasibility boundary (0.2 < prob_feasible < 0.8)
  if (cv_estimate > 0.18 && prob_feasible >= 0.2 && prob_feasible <= 0.8) {
    if ("high" %in% names(fidelity_levels)) {
      return("high")
    } else if ("med" %in% names(fidelity_levels)) {
      return("med")
    }
  }

  # Otherwise use medium fidelity for promising regions
  if (prob_feasible >= 0.4 && "med" %in% names(fidelity_levels)) {
    return("med")
  }

  # Default to low fidelity
  names(fidelity_levels)[1]
}

#' Legacy fidelity selection (simple threshold-based)
#' @keywords internal
select_fidelity <- function(prob_feasible, fidelity_levels) {
  if (length(fidelity_levels) == 1L) {
    return(names(fidelity_levels))
  }
  if (prob_feasible >= 0.75 && "high" %in% names(fidelity_levels)) {
    "high"
  } else if (prob_feasible >= 0.4 && "med" %in% names(fidelity_levels)) {
    "med"
  } else {
    names(fidelity_levels)[1]
  }
}

#' @keywords internal
generate_diagnostics <- function(surrogates, history, objective, n_draws = 1000) {
  unit_best <- history$unit_x[[best_feasible_index(history, objective)]]
  draws <- purrr::imap(surrogates, function(model, metric) {
    newdata <- matrix(unname(unlist(unit_best)), nrow = 1)
    colnames(newdata) <- names(unit_best)
    pred <- DiceKriging::predict.km(model,
                                    newdata = newdata,
                                    type = "UK",
                                    se.compute = TRUE,
                                    cov.compute = TRUE)
    mean <- as.numeric(pred$mean)
    sd <- max(as.numeric(pred$sd), 1e-8)
    stats::rnorm(n_draws, mean, sd)
  })
  list(
    posterior_draws = draws,
    history_summary = history |> dplyr::select(eval_id, objective, feasible)
  )
}

#' @keywords internal
evaluate_acquisition <- function(acquisition,
                                 unit_candidates,
                                 surrogates,
                                 constraint_tbl,
                                 objective,
                                 best_feasible) {
  if (identical(acquisition, "eci")) {
    return(acq_eci(
      unit_x = unit_candidates,
      surrogates = surrogates,
      constraint_tbl = constraint_tbl,
      objective = objective,
      best_feasible = best_feasible
    ))
  }
  if (identical(acquisition, "qehvi")) {
    return(acq_qehvi(
      unit_x = unit_candidates,
      surrogates = surrogates,
      constraint_tbl = constraint_tbl,
      objective = objective,
      best_feasible = best_feasible
    ))
  }
  stop(sprintf("Acquisition '%s' is not supported.", acquisition), call. = FALSE)
}
