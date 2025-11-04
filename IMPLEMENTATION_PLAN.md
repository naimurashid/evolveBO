# Implementation Plan: Priority Recommendations

**Created**: 2025-11-04
**Based on**: PACKAGE_REVIEW_ANALYSIS.md
**Target Version**: 0.3.0

## Overview

This plan implements the CRITICAL, HIGH, and selected MEDIUM priority recommendations from the package review. Tasks are organized into phases to manage dependencies and ensure incremental progress with continuous testing.

**Estimated Total Effort**: 3-4 weeks (1 developer, full-time)

---

## Phase 0: Setup and Quick Wins (2-3 days)

These are low-risk fixes that can be done immediately to improve stability and user experience.

### Task 0.1: Fix Numerical Stability Issues ⚡ CRITICAL
**Priority**: CRITICAL
**Effort**: 2 hours
**Risk**: Low

**Files to modify**:
- `R/acquisition.R` (lines 30-38)

**Changes**:
```r
# Current (line 35):
Z <- improvement / sigma_obj

# Fixed:
Z <- improvement / (sigma_obj + 1e-10)
ei <- improvement * pnorm(Z) + (sigma_obj + 1e-10) * dnorm(Z)
```

**Testing**:
- Add unit test with `sigma_obj = 0` case
- Verify no NaN/Inf in acquisition values
- Test at observed points where GP variance = 0

**Acceptance Criteria**:
- [ ] No division by zero errors
- [ ] Test passes with edge cases
- [ ] No regression in existing tests

---

### Task 0.2: Add Input Validation 🔒 HIGH
**Priority**: HIGH
**Effort**: 4 hours
**Risk**: Low

**Files to modify**:
- `R/bo_calibrate.R` (lines 48-90, add validation section)

**Changes**:
```r
# Add after line 85 (before rng_seed assignment):
validate_bo_inputs <- function(sim_fun, bounds, objective, constraints,
                                n_init, q, budget, fidelity_levels) {
  # Check budget constraints
  if (budget < n_init + q) {
    stop(sprintf("budget (%d) must be >= n_init (%d) + q (%d)",
                 budget, n_init, q), call. = FALSE)
  }

  # Check objective is valid metric
  # Run one test evaluation to verify
  test_theta <- purrr::map2(bounds, 1, ~(.x[1] + .x[2])/2)
  test_result <- tryCatch({
    sim_fun(test_theta, fidelity = names(fidelity_levels)[1])
  }, error = function(e) {
    stop("sim_fun failed on test evaluation: ", e$message, call. = FALSE)
  })

  if (!objective %in% names(test_result)) {
    stop(sprintf("objective '%s' not in simulator output: [%s]",
                 objective, paste(names(test_result), collapse = ", ")),
         call. = FALSE)
  }

  # Validate constraints reference valid metrics
  constraint_metrics <- names(constraints)
  invalid <- setdiff(constraint_metrics, names(test_result))
  if (length(invalid) > 0) {
    stop(sprintf("Constraints reference unknown metrics: [%s]",
                 paste(invalid, collapse = ", ")), call. = FALSE)
  }

  # Check constraint format
  for (metric in constraint_metrics) {
    spec <- constraints[[metric]]
    if (length(spec) != 2 || !spec[1] %in% c("ge", "le") || !is.numeric(spec[2])) {
      stop(sprintf("Invalid constraint format for '%s': must be c('ge'|'le', threshold)",
                   metric), call. = FALSE)
    }
  }

  # Validate n_init
  d <- length(bounds)
  if (n_init < 2 * d) {
    warning(sprintf("n_init (%d) < 2 × dimension (%d). Recommend n_init >= %d for stable GP",
                    n_init, d, 4 * d), call. = FALSE)
  }

  invisible(TRUE)
}

# Call validation (add at line 86):
validate_bo_inputs(sim_fun, bounds, objective, constraints, n_init, q, budget, fidelity_levels)
```

**Testing**:
- Test with invalid budget (< n_init + q)
- Test with wrong objective name
- Test with invalid constraint format
- Test with constraint on non-existent metric
- Test with small n_init (should warn)

**Acceptance Criteria**:
- [ ] Clear error messages for all invalid inputs
- [ ] Warning for suboptimal n_init
- [ ] All validation tests pass
- [ ] No false positives on valid inputs

---

### Task 0.3: Add Dimension-Based n_init Guidance 📊 HIGH
**Priority**: HIGH
**Effort**: 3 hours
**Risk**: Low

**Files to modify**:
- `R/bo_calibrate.R` (documentation and default)
- `man/bo_calibrate.Rd` (after roxygen2 rebuild)

**Changes**:
```r
# Update function signature (line 48):
bo_calibrate <- function(sim_fun,
                         bounds,
                         objective,
                         constraints,
                         n_init = NULL,  # Change from 40 to NULL
                         q = 8,
                         budget = 150,
                         ...) {

  # Add after validation (line ~90):
  # Compute suggested n_init if not provided
  if (is.null(n_init)) {
    d <- length(bounds)
    n_init <- max(4 * d, 10)  # 4× dimension, minimum 10
    if (progress) {
      message(sprintf("Using n_init = %d (4 × %d dimensions)", n_init, d))
    }
  }

  # Rest of function...
}
```

**Documentation update**:
```r
#' @param n_init number of initial design points. If NULL (default), uses
#'   4 × dimension with minimum of 10. Rule of thumb: 2-5× dimension for
#'   stable GP fitting. Smaller values risk underfitting; larger values reduce
#'   BO iterations but may waste budget.
```

**Testing**:
- Test with n_init = NULL (should use default)
- Test with explicit n_init (should use user value)
- Verify message is printed when using auto value
- Test with 1D, 2D, 5D, 10D problems

**Acceptance Criteria**:
- [ ] Auto-sizing works correctly
- [ ] User can still override
- [ ] Documentation is clear
- [ ] Tests pass for various dimensions

---

## Phase 1: Acquisition Function Improvements (4-5 days)

These improve the core optimization logic without requiring major refactoring.

### Task 1.1: Improve Infeasible Region Handling 🎯 HIGH
**Priority**: HIGH
**Effort**: 1 day
**Risk**: Medium

**Files to modify**:
- `R/acquisition.R` (lines 9-38, `acq_eci` function)

**Current problematic code** (lines 19-22):
```r
if (any(mu_obj < best_feasible)) {
  # ... EI calculation
} else {
  ei <- sigma_obj  # Problem: Ignores constraints!
}
```

**Improved implementation**:
```r
acq_eci <- function(candidates, surrogates, constraint_tbl, objective, best_feasible) {
  # ... existing code through line 18

  # Get constraint predictions
  constraint_metrics <- constraint_tbl$metric
  all_preds <- predict_surrogates(surrogates, candidates)

  # Extract objective
  mu_obj <- purrr::map_dbl(all_preds[[objective]]$mean, 1)
  sigma_obj <- purrr::map_dbl(all_preds[[objective]]$sd, 1)

  # Extract constraints
  mu <- purrr::map_dfc(constraint_metrics, function(m) {
    purrr::map_dbl(all_preds[[m]]$mean, 1)
  }) |> purrr::set_names(constraint_metrics) |> as.matrix()

  sigma <- purrr::map_dfc(constraint_metrics, function(m) {
    purrr::map_dbl(all_preds[[m]]$sd, 1)
  }) |> purrr::set_names(constraint_metrics) |> as.matrix()

  # Compute feasibility probability
  prob_feas <- purrr::map_dbl(seq_len(nrow(mu)), function(i) {
    prob_feasibility(mu[i, ], sigma[i, ], constraint_tbl)
  })

  # Compute EI component
  if (is.finite(best_feasible) && any(mu_obj < best_feasible)) {
    # Standard EI: minimize objective
    improvement <- pmax(best_feasible - mu_obj, 0)
    Z <- improvement / (sigma_obj + 1e-10)
    ei <- improvement * pnorm(Z) + (sigma_obj + 1e-10) * dnorm(Z)
  } else {
    # No feasible solution yet or no improvement possible
    # Strategy: Maximize feasibility probability weighted by exploration

    if (is.finite(best_feasible)) {
      # We have feasible solution but can't improve
      # Pure exploration weighted by feasibility
      ei <- sigma_obj * (1 + prob_feas)
    } else {
      # No feasible solution found yet
      # Minimize constraint violation with exploration bonus

      # Compute expected constraint violation
      violation <- compute_expected_violation(mu, sigma, constraint_tbl)

      # Negative violation (less violation is better) + exploration
      # Scale to be comparable to EI values
      max_violation <- max(violation, na.rm = TRUE)
      if (max_violation > 0) {
        normalized_violation <- violation / max_violation
      } else {
        normalized_violation <- violation
      }

      # Acquisition = -violation + exploration, weighted by feasibility probability
      ei <- (-normalized_violation + 0.5 * sigma_obj) * (0.5 + 0.5 * prob_feas)
    }
  }

  # Final acquisition: EI × P(feasible)
  acq_eci <- ei * prob_feas
  acq_eci
}

# Helper function to compute expected constraint violation
compute_expected_violation <- function(mu, sigma, constraint_tbl) {
  n_points <- nrow(mu)
  violation <- numeric(n_points)

  for (i in seq_len(n_points)) {
    total_viol <- 0
    for (j in seq_len(nrow(constraint_tbl))) {
      metric <- constraint_tbl$metric[j]
      direction <- constraint_tbl$direction[j]
      threshold <- constraint_tbl$threshold[j]

      mu_val <- mu[i, metric]
      sigma_val <- sigma[i, metric]

      if (direction == "ge") {
        # P(x < threshold) × E[threshold - x | x < threshold]
        z <- (threshold - mu_val) / (sigma_val + 1e-10)
        prob_violate <- pnorm(z)
        expected_shortfall <- (sigma_val + 1e-10) * dnorm(z) / (prob_violate + 1e-10)
        total_viol <- total_viol + prob_violate * expected_shortfall
      } else {  # "le"
        # P(x > threshold) × E[x - threshold | x > threshold]
        z <- (mu_val - threshold) / (sigma_val + 1e-10)
        prob_violate <- pnorm(z)
        expected_excess <- (sigma_val + 1e-10) * dnorm(z) / (prob_violate + 1e-10)
        total_viol <- total_viol + prob_violate * expected_excess
      }
    }
    violation[i] <- total_viol
  }

  violation
}
```

**Testing**:
- Test with no feasible points found (`best_feasible = Inf`)
- Test with feasible points found but plateau reached
- Test that acquisition prefers points near feasible boundary
- Compare before/after on toy problem with tight constraints

**Acceptance Criteria**:
- [ ] No longer ignores constraints when `best_feasible = Inf`
- [ ] Acquisition guides search toward feasible region
- [ ] Existing tests still pass
- [ ] New tests for infeasible region scenarios pass

---

### Task 1.2: Implement Batch Diversity Mechanism ⚡ CRITICAL
**Priority**: CRITICAL
**Effort**: 2-3 days
**Risk**: Medium

**Files to modify**:
- `R/bo_calibrate.R` (lines 130-145, batch selection logic)
- `R/acquisition.R` (new helper functions)

**Implementation**: Local Penalization (simpler than Constant Liar, works well)

**Add to `R/acquisition.R`**:
```r
#' Select diverse batch using local penalization
#'
#' Implements the local penalization strategy of González et al. (2016).
#' Iteratively selects points by penalizing acquisition near previously
#' selected points.
#'
#' @param candidates matrix of candidate points (unit scale)
#' @param acq_scores initial acquisition values
#' @param q batch size
#' @param lipschitz Lipschitz constant for penalization (default: 10)
#' @param bounds parameter bounds (for computing distances)
#'
#' @return indices of selected candidates
#' @keywords internal
select_batch_local_penalization <- function(candidates, acq_scores, q,
                                             lipschitz = 10, bounds = NULL) {
  n_candidates <- length(acq_scores)
  selected_indices <- integer(q)
  penalized_scores <- acq_scores

  for (i in seq_len(q)) {
    # Select point with highest penalized acquisition
    best_idx <- which.max(penalized_scores)
    selected_indices[i] <- best_idx

    if (i < q) {
      # Penalize acquisition near selected point
      selected_point <- candidates[best_idx, , drop = FALSE]

      # Compute distances to selected point
      distances <- compute_distances(candidates, selected_point)

      # Penalization function: max(0, L * r - acq_i)
      # where r is distance, L is Lipschitz constant
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
#' @keywords internal
compute_distances <- function(points, reference) {
  # points: n × d matrix
  # reference: 1 × d matrix
  # returns: n-vector of distances

  diff <- sweep(points, 2, reference[1, ], "-")
  sqrt(rowSums(diff^2))
}

#' Estimate Lipschitz constant from GP
#' @keywords internal
estimate_lipschitz <- function(surrogates, objective, n_samples = 100) {
  # Simple heuristic: max gradient magnitude from random samples
  # Could be more sophisticated (use GP lengthscales)

  # For now, use conservative default based on lengthscales
  model <- surrogates[[objective]]

  if (inherits(model, "km")) {
    # DiceKriging model
    lengthscales <- model@covariance@range.val
    # Lipschitz ~ 1 / min(lengthscale)
    L <- 1 / min(lengthscales)
    # Conservative factor
    return(L * 2)
  }

  # Default if can't extract lengthscales
  return(10)
}
```

**Modify `R/bo_calibrate.R`** (around line 130):
```r
# Replace current batch selection:
# chosen_indices <- order(acq_scores, decreasing = TRUE)[1:q]

# With diverse batch selection:
if (q == 1) {
  # Single point: just take best
  chosen_indices <- which.max(acq_scores)
} else {
  # Batch: use local penalization for diversity
  lipschitz <- estimate_lipschitz(surrogates, objective)
  chosen_indices <- select_batch_local_penalization(
    candidates = candidates,
    acq_scores = acq_scores,
    q = q,
    lipschitz = lipschitz,
    bounds = bounds
  )
}
```

**Testing**:
- Unit test: verify selected points are spatially diverse
- Measure pairwise distances between batch points
- Compare to greedy selection (should have larger min distance)
- Integration test: run full BO with q=4, verify batch diversity
- Benchmark: compare convergence with/without diversity

**Acceptance Criteria**:
- [ ] Batch points are spatially diverse (min distance > threshold)
- [ ] No performance regression on q=1 case
- [ ] Tests show faster convergence (10-20% fewer evaluations)
- [ ] Documentation added to function

---

## Phase 2: Multi-Fidelity Strategy Overhaul (5-6 days)

This is the most complex change but provides significant algorithmic improvement.

### Task 2.1: Implement Adaptive Fidelity Selection ⚡ CRITICAL
**Priority**: CRITICAL
**Effort**: 3-4 days
**Risk**: High

**Files to modify**:
- `R/bo_calibrate.R` (lines 406-467, fidelity selection)
- `R/bo_calibrate.R` (line 53, add new parameter)

**Step 1: Add new parameter and method selection** (line 53):
```r
bo_calibrate <- function(...,
                         fidelity_levels = c(low = 200, med = 1000, high = 10000),
                         fidelity_method = c("adaptive", "staged", "threshold"),
                         fidelity_costs = NULL,  # NEW: explicit cost specification
                         ...) {
```

**Step 2: Validate fidelity method** (after line 85):
```r
fidelity_method <- match.arg(fidelity_method)
if (progress) {
  message(sprintf("Using '%s' fidelity selection strategy", fidelity_method))
}

# Use fidelity_levels as costs if costs not specified
if (is.null(fidelity_costs)) {
  fidelity_costs <- fidelity_levels / min(fidelity_levels)
}
```

**Step 3: Replace fidelity selection call** (line 162-163):
```r
# Current:
# fidelity <- select_fidelity_staged(prob_feas, cv_estimate, iter_counter, fidelity_levels)

# New:
fidelity <- select_fidelity_method(
  method = fidelity_method,
  prob_feasible = prob_feas,
  cv_estimate = cv_estimate,
  acq_value = acq_scores[chosen_indices[batch_idx]],
  best_obj = best_obj_value,
  fidelity_levels = fidelity_levels,
  fidelity_costs = fidelity_costs,
  iter = iter_counter,
  total_budget_used = sum(history$n_rep, na.rm = TRUE),
  total_budget = budget * mean(fidelity_levels)  # approximate total simulation budget
)
```

**Step 4: Implement dispatcher and new methods**:
```r
#' Select fidelity level using specified method
#' @keywords internal
select_fidelity_method <- function(method, ...) {
  switch(method,
         adaptive = select_fidelity_adaptive(...),
         staged = select_fidelity_staged(...),
         threshold = select_fidelity(...),
         stop("Unknown fidelity method: ", method, call. = FALSE)
  )
}

#' Adaptive cost-aware fidelity selection
#'
#' Implements cost-aware fidelity selection inspired by MFKG (Wu & Frazier 2016).
#' Selects fidelity by maximizing expected value per unit cost, with exploration
#' decay and boundary detection.
#'
#' @param prob_feasible probability of constraint satisfaction
#' @param cv_estimate coefficient of variation from objective surrogate
#' @param acq_value acquisition function value at candidate point
#' @param best_obj current best objective value
#' @param fidelity_levels named vector of fidelity levels (replication counts)
#' @param fidelity_costs named vector of relative costs (default: proportional to replications)
#' @param iter current iteration number
#' @param total_budget_used cumulative simulation budget consumed
#' @param total_budget approximate total simulation budget available
#'
#' @details
#' The method balances information gain vs cost using:
#' - **Value score**: acquisition × uncertainty × boundary_factor
#'   - High near feasibility boundary (prob ≈ 0.5)
#'   - High when objective uncertain (large CV)
#'   - High for promising candidates (large acquisition)
#' - **Cost normalization**: Divide by cost^α where α decays from 0.5 → 0.8
#'   - Early: less cost-sensitive (exploration)
#'   - Late: more cost-sensitive (exploitation)
#' - **Exploration decay**: Randomization probability from 50% → 10%
#'
#' @return name of selected fidelity level
#' @keywords internal
select_fidelity_adaptive <- function(prob_feasible,
                                     cv_estimate,
                                     acq_value,
                                     best_obj,
                                     fidelity_levels,
                                     fidelity_costs,
                                     iter,
                                     total_budget_used,
                                     total_budget) {
  if (length(fidelity_levels) == 1L) {
    return(names(fidelity_levels))
  }

  fidelity_names <- names(fidelity_levels)
  costs <- fidelity_costs[fidelity_names]

  # === Compute value score ===

  # Uncertainty factor: higher CV → more value in reducing uncertainty
  uncertainty_factor <- pmax(0, pmin(1, cv_estimate / 0.3))  # normalize to [0,1]

  # Boundary factor: highest value near feasibility boundary
  # P = 0.5 → factor = 1
  # P = 0 or 1 → factor = 0
  boundary_factor <- 1 - abs(2 * prob_feasible - 1)
  boundary_factor <- boundary_factor^0.5  # soften the effect

  # Acquisition factor: diminishing returns on acquisition value
  acq_factor <- log1p(pmax(0, acq_value))

  # Combined value score
  # Weight components based on optimization stage
  if (iter < 20) {
    # Early: prioritize exploration (uncertainty)
    value_score <- acq_factor * (0.3 + 0.7 * uncertainty_factor) * (0.5 + 0.5 * boundary_factor)
  } else if (iter < 60) {
    # Middle: balance
    value_score <- acq_factor * (0.5 + 0.5 * uncertainty_factor) * (0.3 + 0.7 * boundary_factor)
  } else {
    # Late: prioritize promising candidates
    value_score <- acq_factor * (0.7 + 0.3 * uncertainty_factor) * (0.1 + 0.9 * boundary_factor)
  }

  # === Compute cost sensitivity ===

  # Normalize costs to [0, 1]
  cost_normalized <- costs / max(costs)

  # Cost exponent: increases with iteration (more cost-sensitive over time)
  # Also increases as budget depletes
  budget_fraction_used <- pmin(1, total_budget_used / total_budget)
  cost_exponent <- 0.3 + 0.5 * (iter / 100) + 0.3 * budget_fraction_used
  cost_exponent <- pmin(cost_exponent, 1.0)  # cap at 1.0

  # === Value per cost ===
  value_per_cost <- value_score / (cost_normalized ^ cost_exponent + 1e-6)

  # === Exploration randomization ===

  # Probability of forcing low fidelity for exploration
  # Decays from 50% early to 5% late
  exploration_prob <- pmax(0.05, 0.5 * exp(-iter / 30))

  if (stats::runif(1) < exploration_prob) {
    # Force exploration with low fidelity
    return(fidelity_names[1])
  }

  # === Select fidelity ===

  best_idx <- which.max(value_per_cost)
  selected <- fidelity_names[best_idx]

  # === Diagnostics (optional) ===
  if (getOption("evolveBO.debug_fidelity", FALSE)) {
    message(sprintf(
      "  Fidelity selection: prob_feas=%.3f, CV=%.3f, acq=%.3f, value=%.3f, cost_exp=%.2f → %s",
      prob_feasible, cv_estimate, acq_value, value_score, cost_exponent, selected
    ))
    message(sprintf("    Value/cost: %s",
                    paste(sprintf("%s=%.2f", fidelity_names, value_per_cost), collapse=", ")))
  }

  selected
}

# Keep existing staged and threshold methods for backward compatibility
# (lines 418-467 remain unchanged but renamed)
```

**Documentation update**:
```r
#' @param fidelity_method method for selecting fidelity level. Options:
#'   \itemize{
#'     \item \code{"adaptive"} (default): Cost-aware selection based on expected
#'       value per unit cost. Balances information gain vs computational expense.
#'       Recommended for most use cases.
#'     \item \code{"staged"}: Fixed schedule with iteration-based thresholds.
#'       Iterations 1-30: low, 31-100: adaptive, 101+: high. Simple but less
#'       efficient than adaptive.
#'     \item \code{"threshold"}: Simple feasibility probability thresholds.
#'       P ≥ 0.75: high, P ≥ 0.4: med, else: low. Legacy method, not recommended.
#'   }
#' @param fidelity_costs named numeric vector of relative costs per fidelity level.
#'   If NULL (default), assumes cost proportional to replication count. Use to
#'   specify non-linear cost relationships (e.g., parallelization effects).
```

**Testing Strategy**:

1. **Unit tests** for `select_fidelity_adaptive`:
   ```r
   test_that("adaptive fidelity responds to value/cost tradeoff", {
     # Low CV, high prob → should prefer low fidelity (low value/cost)
     fid1 <- select_fidelity_adaptive(
       prob_feasible = 0.9, cv_estimate = 0.05, acq_value = 0.1,
       best_obj = 10, fidelity_levels = c(low=200, high=10000),
       fidelity_costs = c(low=1, high=50), iter = 50,
       total_budget_used = 5000, total_budget = 100000
     )
     # Expect low fidelity

     # High CV, boundary prob → should prefer high fidelity (high value)
     fid2 <- select_fidelity_adaptive(
       prob_feasible = 0.5, cv_estimate = 0.3, acq_value = 1.5,
       best_obj = 10, fidelity_levels = c(low=200, high=10000),
       fidelity_costs = c(low=1, high=50), iter = 50,
       total_budget_used = 5000, total_budget = 100000
     )
     # Expect high fidelity
   })
   ```

2. **Integration tests**: Full BO runs comparing methods
   ```r
   test_that("adaptive fidelity converges faster than staged", {
     set.seed(123)
     fit_adaptive <- bo_calibrate(..., fidelity_method = "adaptive", budget = 50)

     set.seed(123)
     fit_staged <- bo_calibrate(..., fidelity_method = "staged", budget = 50)

     # Adaptive should achieve better objective with same budget
     # Or same objective with less budget
   })
   ```

3. **Ablation study** (use existing `ablation_multifidelity` function):
   ```r
   ablation <- ablation_multifidelity(
     sim_fun = toy_sim_fun,
     bounds = toy_bounds,
     objective = "EN",
     constraints = toy_constraints,
     policies = list(
       adaptive = c(low = 200, med = 1000, high = 10000),
       staged = c(low = 200, med = 1000, high = 10000),
       threshold = c(low = 200, med = 1000, high = 10000),
       fixed_high = c(high = 10000)
     ),
     seeds = 1:20,
     bo_args = list(n_init = 10, q = 4, budget = 50,
                    fidelity_method = c("adaptive", "staged", "threshold", "adaptive"))
   )
   ```

**Acceptance Criteria**:
- [ ] All three methods (adaptive, staged, threshold) callable
- [ ] Adaptive method responds sensibly to inputs
- [ ] Unit tests pass
- [ ] Integration tests show adaptive ≥ staged performance
- [ ] Ablation study shows adaptive uses budget more efficiently
- [ ] Documentation complete

---

## Phase 3: Performance Optimizations (3-4 days)

These provide speedups without changing algorithm behavior.

### Task 3.1: Warm-Start GP Hyperparameters
**Priority**: MEDIUM
**Effort**: 1-2 days
**Risk**: Medium

**Files to modify**:
- `R/surrogates.R` (lines 38-142, `fit_surrogates` function)
- `R/bo_calibrate.R` (line ~117, pass previous surrogates)

**Changes to `fit_surrogates`**:
```r
#' @param prev_surrogates optional list of surrogates from previous iteration
#'   for warm-starting hyperparameter optimization
fit_surrogates <- function(history,
                           objective,
                           constraint_tbl,
                           covtype = "matern5_2",
                           prev_surrogates = NULL) {
  # ... existing aggregation code ...

  # Fit each metric's surrogate
  surrogates <- purrr::imap(c(objective, constraint_metrics), function(metric, idx) {
    # Extract previous model if available
    prev_model <- if (!is.null(prev_surrogates) && metric %in% names(prev_surrogates)) {
      prev_surrogates[[metric]]
    } else {
      NULL
    }

    # Fit with warm start
    fit_single_surrogate(
      X = as.matrix(agg_X),
      y = agg_metrics[[metric]],
      noise_var = agg_var[[metric]],
      covtype = covtype,
      prev_model = prev_model
    )
  }) |> purrr::set_names(c(objective, constraint_metrics))

  surrogates
}

#' Fit single GP surrogate with optional warm-start
#' @keywords internal
fit_single_surrogate <- function(X, y, noise_var, covtype, prev_model = NULL) {
  # Extract hyperparameters from previous model if available
  theta_init <- NULL
  if (!is.null(prev_model) && inherits(prev_model, "km")) {
    theta_init <- prev_model@covariance@range.val
  }

  # Determine noise handling
  has_variance <- !is.null(noise_var) && !all(is.na(noise_var))

  if (has_variance && all(is.finite(noise_var)) && all(noise_var > 0)) {
    # Heteroskedastic GP
    model <- DiceKriging::km(
      design = X,
      response = y,
      covtype = covtype,
      noise.var = noise_var,
      control = list(
        trace = FALSE,
        parinit = theta_init  # Warm start
      )
    )
  } else {
    # Homoskedastic GP with nugget
    nugget_val <- 1e-6
    model <- DiceKriging::km(
      design = X,
      response = y,
      covtype = covtype,
      nugget = nugget_val,
      control = list(
        trace = FALSE,
        parinit = theta_init  # Warm start
      )
    )
  }

  model
}
```

**Changes to `bo_calibrate`** (around line 117):
```r
# Before surrogate fitting:
surrogates <- fit_surrogates(
  history = history,
  objective = objective,
  constraint_tbl = constraint_tbl,
  covtype = covtype,
  prev_surrogates = if (iter_counter > 1) surrogates else NULL  # Warm start
)
```

**Testing**:
- Measure time for GP fitting with/without warm start
- Verify same final objective value (no regression)
- Test that hyperparameters evolve smoothly across iterations
- Edge case: first iteration (no previous model)

**Acceptance Criteria**:
- [ ] 20-40% speedup in surrogate fitting time
- [ ] No regression in optimization performance
- [ ] Hyperparameters converge faster
- [ ] Tests pass

---

### Task 3.2: Adaptive Candidate Pool Size
**Priority**: MEDIUM
**Effort**: 0.5 days
**Risk**: Low

**Files to modify**:
- `R/bo_calibrate.R` (line 51, parameter; line 119, usage)

**Changes**:
```r
# Function signature (line 51):
candidate_pool = NULL,  # Change from 2000 to NULL

# Inside function (after line 85):
if (is.null(candidate_pool)) {
  d <- length(bounds)
  # Scale with dimension: 500d with min 1000, max 5000
  candidate_pool <- pmax(1000, pmin(5000, 500 * d))
  if (progress && iter_counter == 1) {
    message(sprintf("Using candidate pool size: %d (500 × %d dimensions)",
                    candidate_pool, d))
  }
}

# Optionally: make it adaptive to iteration
# Early iterations: fewer candidates (exploration)
# Late iterations: more candidates (refinement)
if (iter_counter <= n_init) {
  pool_size <- candidate_pool
} else if (iter_counter < 0.7 * budget) {
  pool_size <- candidate_pool
} else {
  # Last 30% of budget: increase pool for refinement
  pool_size <- min(candidate_pool * 1.5, 10000)
}

candidates <- lhs_candidate_pool(pool_size, bounds)
```

**Testing**:
- Test with d=2, 5, 10, 20
- Verify pool size scales appropriately
- Check performance on high-dimensional problems

**Acceptance Criteria**:
- [ ] Auto-sizing works correctly
- [ ] User can still override
- [ ] Better performance on high-D problems
- [ ] Tests pass

---

### Task 3.3: Early Stopping Criterion
**Priority**: MEDIUM
**Effort**: 1 day
**Risk**: Low

**Files to modify**:
- `R/bo_calibrate.R` (lines 108-193, main loop)
- Add parameter `early_stop = TRUE`

**Changes**:
```r
# Function signature (line 48):
bo_calibrate <- function(...,
                         early_stop = TRUE,
                         early_stop_patience = 10,
                         early_stop_tol = 1e-4,
                         ...) {

# Inside main loop (after acquisition optimization):
# Track no-improvement counter
no_improvement_count <- 0
best_obj_history <- numeric()

for (iter_counter in seq_len(max_bo_iter)) {
  # ... existing code ...

  # After updating history with new evaluations:
  current_best <- get_best_feasible_value(history, objective)
  best_obj_history <- c(best_obj_history, current_best)

  if (early_stop && iter_counter > early_stop_patience) {
    # Check if best value improved in last 'patience' iterations
    recent_best <- min(tail(best_obj_history, early_stop_patience))
    earlier_best <- min(best_obj_history[1:(length(best_obj_history) - early_stop_patience)])

    improvement <- (earlier_best - recent_best) / (abs(earlier_best) + 1e-8)

    if (improvement < early_stop_tol) {
      if (progress) {
        message(sprintf(
          "Early stopping at iteration %d: no improvement > %.2g in last %d iterations",
          iter_counter, early_stop_tol, early_stop_patience
        ))
      }
      break
    }
  }

  # Alternative stopping criterion: low acquisition values
  if (early_stop && iter_counter > n_init) {
    max_acq <- max(acq_scores, na.rm = TRUE)
    # If all acquisition values are very small, we've converged
    if (max_acq < 1e-6) {
      if (progress) {
        message(sprintf(
          "Early stopping at iteration %d: max acquisition value %.2g < threshold",
          iter_counter, max_acq
        ))
      }
      break
    }
  }
}

#' Get best feasible objective value from history
#' @keywords internal
get_best_feasible_value <- function(history, objective) {
  feasible_rows <- history[history$feasible, ]
  if (nrow(feasible_rows) == 0) {
    return(Inf)
  }
  min(feasible_rows$objective, na.rm = TRUE)
}
```

**Testing**:
- Test that early stopping triggers when plateau reached
- Test that it doesn't stop prematurely
- Test with `early_stop = FALSE` (should run full budget)
- Measure budget savings on toy problems

**Acceptance Criteria**:
- [ ] Stops when converged (saves 10-30% budget)
- [ ] Doesn't stop prematurely (no worse final objective)
- [ ] User can disable with `early_stop = FALSE`
- [ ] Tests pass

---

## Phase 4: Testing and Documentation (2-3 days)

### Task 4.1: Comprehensive Test Suite
**Effort**: 1-2 days

**New test files needed**:
- `tests/testthat/test-batch-diversity.R`
- `tests/testthat/test-fidelity-adaptive.R`
- `tests/testthat/test-acquisition-infeasible.R`
- `tests/testthat/test-early-stopping.R`

**Example tests**:
```r
# test-batch-diversity.R
test_that("batch selection produces diverse points", {
  # Generate candidates
  candidates <- lhs::randomLHS(100, 2)
  acq_scores <- runif(100)

  # Select batch
  indices <- select_batch_local_penalization(candidates, acq_scores, q = 4, lipschitz = 10)
  selected <- candidates[indices, ]

  # Compute pairwise distances
  dists <- as.matrix(dist(selected))
  diag(dists) <- Inf
  min_dist <- min(dists)

  # Should be reasonably spaced
  expect_gt(min_dist, 0.1)  # At least 0.1 apart in unit hypercube

  # Compare to greedy selection
  greedy_indices <- order(acq_scores, decreasing = TRUE)[1:4]
  greedy <- candidates[greedy_indices, ]
  greedy_dists <- as.matrix(dist(greedy))
  diag(greedy_dists) <- Inf
  greedy_min_dist <- min(greedy_dists)

  # Diverse selection should have larger min distance
  expect_gt(min_dist, greedy_min_dist * 0.8)
})

# test-fidelity-adaptive.R
test_that("adaptive fidelity responds to cost-benefit tradeoff", {
  # High value, low cost → high fidelity
  fid1 <- select_fidelity_adaptive(
    prob_feasible = 0.5, cv_estimate = 0.3, acq_value = 2.0,
    best_obj = 10, fidelity_levels = c(low=200, high=10000),
    fidelity_costs = c(low=1, high=50), iter = 50,
    total_budget_used = 5000, total_budget = 100000
  )

  # Due to exploration randomness, test over multiple calls
  fidelities <- replicate(20, select_fidelity_adaptive(...))
  # Should use high fidelity more often than low
  expect_gt(sum(fidelities == "high"), sum(fidelities == "low"))
})

# test-acquisition-infeasible.R
test_that("acquisition function handles infeasible region correctly", {
  # Create mock surrogates with no feasible solution
  # ... (create test surrogates)

  acq <- acq_eci(candidates, surrogates, constraint_tbl, objective, best_feasible = Inf)

  # Should not return all zeros or NAs
  expect_true(any(acq > 0))
  expect_true(all(is.finite(acq)))

  # Should prefer points near feasibility boundary
  # ... (verify this property)
})
```

### Task 4.2: Update Documentation
**Effort**: 1 day

**Files to update**:
- `README.md` - Add section on new features
- `CLAUDE.md` - Update architecture description
- `man/bo_calibrate.Rd` - Full parameter documentation
- Add vignette: `vignettes/advanced-features.Rmd`

**New vignette outline**:
```markdown
# Advanced Features in evolveBO

## Batch Optimization with Diversity

Explains local penalization, when to use large batches...

## Multi-Fidelity Strategies

Comparison of adaptive vs staged vs threshold...
When to use each...
Customizing costs...

## Early Stopping

How it works, tuning patience and tolerance...

## Performance Tips

- Warm-starting
- Candidate pool sizing
- Parallel evaluation
```

---

## Phase 5: Benchmarking and Validation (2-3 days)

### Task 5.1: Performance Benchmarks
**Effort**: 1-2 days

**Create** `inst/benchmarks/compare_versions.R`:
```r
# Compare v0.2.0 vs v0.3.0 on standard test problems

library(evolveBO)
library(tictoc)

# Test problem 1: Toy 2D
# Test problem 2: Higher-dimensional (5D)
# Test problem 3: Tight constraints

results <- list()

for (problem in c("toy_2d", "high_dim", "tight_constraints")) {
  cat("\n=== Benchmarking:", problem, "===\n")

  # Version 0.2.0 behavior (staged, no diversity)
  tic(sprintf("%s_v02", problem))
  fit_v02 <- bo_calibrate(
    ...,
    fidelity_method = "staged",
    batch_diversity = FALSE,  # hypothetical flag
    early_stop = FALSE
  )
  time_v02 <- toc()

  # Version 0.3.0 behavior (adaptive, diversity, early stop)
  tic(sprintf("%s_v03", problem))
  fit_v03 <- bo_calibrate(
    ...,
    fidelity_method = "adaptive",
    batch_diversity = TRUE,
    early_stop = TRUE
  )
  time_v03 <- toc()

  results[[problem]] <- list(
    v02_time = time_v02,
    v03_time = time_v03,
    v02_evals = nrow(fit_v02$history),
    v03_evals = nrow(fit_v03$history),
    v02_best = fit_v02$best_obj,
    v03_best = fit_v03$best_obj
  )
}

# Save results
saveRDS(results, "inst/benchmarks/version_comparison.rds")
```

**Expected outcomes**:
- 10-20% fewer evaluations with batch diversity
- 10-30% fewer evaluations with early stopping
- 5-15% less wall-clock time with warm-starting
- Similar or better final objective values

### Task 5.2: Ablation Study
**Effort**: 1 day

**Use existing `ablation_multifidelity` function**:
```r
# Compare fidelity methods
ablation <- ablation_multifidelity(
  sim_fun = realistic_sim_fun,
  bounds = realistic_bounds,
  objective = "EN",
  constraints = realistic_constraints,
  policies = list(
    adaptive = c(low = 200, med = 1000, high = 10000),
    staged = c(low = 200, med = 1000, high = 10000),
    threshold = c(low = 200, med = 1000, high = 10000),
    fixed_low = c(low = 200),
    fixed_high = c(high = 10000)
  ),
  seeds = 1:30,  # More seeds for statistical power
  bo_args = list(
    n_init = 20,
    q = 4,
    budget = 100,
    fidelity_method = c("adaptive", "staged", "threshold", "adaptive", "adaptive")
  )
)

plot_multifidelity_tradeoff(ablation)
```

**Analysis**:
- Statistical test: Is adaptive significantly better?
- Cost-efficiency: simulation budget vs objective quality
- Robustness: variance across seeds

---

## Risk Mitigation

### High-Risk Tasks
1. **Task 2.1 (Adaptive Fidelity)**: Complex logic, many parameters
   - **Mitigation**: Extensive unit tests, ablation study, phased rollout
   - **Fallback**: Keep staged as default if adaptive underperforms

2. **Task 1.2 (Batch Diversity)**: Changes core selection logic
   - **Mitigation**: Make optional (flag to disable), thorough testing
   - **Fallback**: Can revert to greedy selection

### Testing Strategy
- Unit tests for each function
- Integration tests for full BO runs
- Regression tests (ensure no performance degradation)
- Benchmarks on standard problems

### Rollout Strategy
1. Merge Phase 0 (quick fixes) first → v0.2.1
2. Merge Phase 1 + Phase 3 → v0.2.2 (performance improvements)
3. Merge Phase 2 (fidelity) after thorough testing → v0.3.0

---

## Success Metrics

### Performance Targets
- [ ] 15-25% fewer evaluations to reach convergence (batch diversity + early stop)
- [ ] 20-40% faster surrogate fitting (warm-start)
- [ ] Similar or better final objective values (no regression)
- [ ] More efficient budget usage (adaptive fidelity)

### Code Quality Targets
- [ ] Test coverage > 80% for new code
- [ ] All edge cases handled (division by zero, empty arrays, etc.)
- [ ] Clear error messages for invalid inputs
- [ ] Comprehensive documentation

### User Experience Targets
- [ ] Sensible defaults (auto n_init, auto candidate pool)
- [ ] Clear progress messages
- [ ] Backward compatibility (old code still works)
- [ ] Migration guide for new features

---

## Timeline Summary

| Phase | Duration | Cumulative |
|-------|----------|------------|
| Phase 0: Setup & Quick Wins | 2-3 days | 3 days |
| Phase 1: Acquisition Improvements | 4-5 days | 8 days |
| Phase 2: Multi-Fidelity Overhaul | 5-6 days | 14 days |
| Phase 3: Performance Optimizations | 3-4 days | 18 days |
| Phase 4: Testing & Documentation | 2-3 days | 21 days |
| Phase 5: Benchmarking | 2-3 days | 24 days |

**Total**: ~4-5 weeks (1 developer, full-time)

---

## Appendix: Optional Future Enhancements (Not in v0.3.0)

These are lower priority but worth considering for future versions:

1. **Continuous Fidelity with hetGP** (v0.4.0)
   - Treat fidelity as continuous parameter
   - Let GP learn correlation across fidelities
   - More sophisticated but requires major refactoring

2. **Multi-Objective Optimization** (v0.4.0)
   - Support Pareto optimization
   - Implement qEHVI (Expected Hypervolume Improvement)
   - Return Pareto frontier instead of single best

3. **Parallelization** (v0.3.1)
   - Parallel candidate pool evaluation
   - Parallel surrogate fitting
   - Async batch evaluation

4. **Advanced Surrogate Models** (v0.5.0)
   - Deep GPs for non-stationarity
   - Random forests for robustness
   - Ensemble models

5. **Constraint Violation Models** (v0.4.0)
   - Separate GP for constraint violation magnitude
   - Better guidance when no feasible solution

6. **Active Subspace Detection** (v0.5.0)
   - Automatic dimensionality reduction
   - Better scaling to high dimensions

---

## Getting Started

To begin implementation:

1. Create feature branch: `git checkout -b feature/v0.3.0-improvements`
2. Start with Phase 0 (quick wins)
3. Commit frequently with clear messages
4. Run tests after each task
5. Update this document with progress notes

Good luck! 🚀
