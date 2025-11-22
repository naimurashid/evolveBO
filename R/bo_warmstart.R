# bo_warmstart.R
# Warm-starting and parameter refinement for Bayesian Optimization
#
# This module provides practical warm-start capabilities for iterative BO calibration:
# 1. Save/load BO state for resumption
# 2. Refine parameter bounds around best designs
# 3. Fix unimportant parameters identified by analysis
#
# These functions work with existing evolveBO::bo_calibrate() output and can be
# integrated into the evolveBO package later.

#' Save BO State for Warm-Starting
#'
#' Saves a BO fit object with all information needed to resume optimization.
#' Includes history, surrogates (if available), and configuration.
#'
#' @param fit Output from evolveBO::bo_calibrate()
#' @param path File path to save (default: adds .rds extension if missing)
#' @param compress Compression level (default: "xz" for maximum compression)
#'
#' @return Invisibly returns the file path
#' @export
#'
#' @examples
#' \dontrun{
#' fit1 <- evolveBO::bo_calibrate(...)
#' save_bo_state(fit1, "results/phase1_fit.rds")
#' }
save_bo_state <- function(fit, path, compress = "xz") {
  # Ensure .rds extension
  if (!grepl("\\.rds$", path)) {
    path <- paste0(path, ".rds")
  }

  # Create directory if needed
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  # Save with compression
  saveRDS(fit, file = path, compress = compress)

  cat(sprintf("✓ BO state saved to: %s\n", path))
  cat(sprintf("  - Evaluations: %d\n", nrow(fit$history)))
  objective_name <- fit$policies$objective

  # Handle case where there are no feasible designs
  feasible_vals <- fit$history[[objective_name]][fit$history$feasible]
  if (length(feasible_vals) > 0 && any(!is.na(feasible_vals))) {
    cat(sprintf("  - Best %s: %.4f\n", objective_name, min(feasible_vals, na.rm = TRUE)))
  } else {
    cat(sprintf("  - Best %s: (no feasible designs)\n", objective_name))
  }

  invisible(path)
}


#' Load BO State
#'
#' Loads a previously saved BO fit object.
#'
#' @param path File path to load
#'
#' @return BO fit object
#' @export
load_bo_state <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("BO state file not found: %s", path))
  }

  fit <- readRDS(path)

  cat(sprintf("✓ BO state loaded from: %s\n", path))
  cat(sprintf("  - Evaluations: %d\n", nrow(fit$history)))
  objective_name <- fit$policies$objective

  # Handle case where there are no feasible designs
  feasible_vals <- fit$history[[objective_name]][fit$history$feasible]
  if (length(feasible_vals) > 0 && any(!is.na(feasible_vals))) {
    cat(sprintf("  - Best %s: %.4f\n", objective_name, min(feasible_vals, na.rm = TRUE)))
  } else {
    cat(sprintf("  - Best %s: (no feasible designs)\n", objective_name))
  }

  return(fit)
}


#' Refine Parameter Bounds Around Best Design
#'
#' Shrinks parameter bounds to focus on region near the best feasible design.
#' Useful for sequential refinement: coarse exploration → focused optimization.
#'
#' @param fit Output from evolveBO::bo_calibrate()
#' @param shrink_factor Fraction to shrink bounds by (0.5 = 50% of original range)
#' @param respect_original_bounds Ensure refined bounds stay within original limits
#' @param min_range Minimum allowed range for each parameter (prevents collapse)
#'
#' @return List of refined bounds
#' @export
#'
#' @examples
#' \dontrun{
#' fit1 <- evolveBO::bo_calibrate(...)
#' bounds_refined <- refine_bounds(fit1, shrink_factor = 0.4)
#'
#' fit2 <- evolveBO::bo_calibrate(
#'   sim_fun = wrapper,
#'   bounds = bounds_refined,  # Use refined bounds
#'   n_init = 20,              # Fewer initial samples needed
#'   budget = 100
#' )
#' }
refine_bounds <- function(fit,
                          shrink_factor = 0.5,
                          respect_original_bounds = TRUE,
                          min_range = 0.05) {

  # Extract best feasible design
  history <- fit$history
  feasible <- history$feasible

  if (!any(feasible)) {
    warning("No feasible designs found. Returning original bounds.")
    return(fit$bounds)
  }

  # Get best feasible design
  feasible_history <- history[feasible, ]
  objective_name <- fit$policies$objective
  best_idx <- which.min(feasible_history[[objective_name]])
  best_design <- feasible_history[best_idx, ]

  # Extract original bounds
  original_bounds <- fit$bounds
  param_names <- names(original_bounds)

  refined_bounds <- list()

  cat("\nRefining bounds around best design (shrink factor:", shrink_factor, ")\n")
  cat("─────────────────────────────────────────────────────────\n")

  for (param in param_names) {
    # Original bounds
    orig_lb <- original_bounds[[param]][1]
    orig_ub <- original_bounds[[param]][2]
    orig_range <- orig_ub - orig_lb

    # Best value
    best_val <- best_design[[param]]

    # New range (shrunk)
    new_range <- max(orig_range * shrink_factor, min_range)

    # New bounds centered on best design
    new_lb <- best_val - new_range / 2
    new_ub <- best_val + new_range / 2

    # Respect original bounds
    if (respect_original_bounds) {
      new_lb <- max(orig_lb, new_lb)
      new_ub <- min(orig_ub, new_ub)

      # If we hit a boundary, extend the other side
      actual_range <- new_ub - new_lb
      if (actual_range < new_range * 0.8) {  # Lost >20% of desired range
        if (new_lb == orig_lb) {
          # Hit lower bound, try to extend upper
          new_ub <- min(orig_ub, new_lb + new_range)
        } else if (new_ub == orig_ub) {
          # Hit upper bound, try to extend lower
          new_lb <- max(orig_lb, new_ub - new_range)
        }
      }
    }

    refined_bounds[[param]] <- c(new_lb, new_ub)

    # Report
    pct_shrink <- 100 * (1 - (new_ub - new_lb) / orig_range)
    cat(sprintf("  %s: [%.3f, %.3f] → [%.3f, %.3f] (−%.0f%%, center: %.3f)\n",
                param, orig_lb, orig_ub, new_lb, new_ub, pct_shrink, best_val))
  }

  cat("─────────────────────────────────────────────────────────\n\n")

  return(refined_bounds)
}


#' Create Wrapper with Fixed Parameters
#'
#' Wraps a simulator function to fix certain parameters at specific values.
#' Useful after identifying unimportant parameters via ARD analysis.
#'
#' @param sim_fun Original simulator function
#' @param fixed_params Named vector of parameters to fix
#' @param original_bounds Original bounds (for validation)
#'
#' @return Wrapped simulator function with reduced dimensionality
#' @export
#'
#' @examples
#' \dontrun{
#' # After ARD analysis identifies that margin has low importance
#' wrapper_reduced <- fix_parameters(
#'   wrapper_original,
#'   fixed_params = c(margin = 0.1, beat = 3),
#'   original_bounds
#' )
#'
#' # Remove fixed parameters from bounds
#' bounds_reduced <- original_bounds[!names(original_bounds) %in% names(fixed_params)]
#'
#' # Run BO with reduced parameter space
#' fit <- evolveBO::bo_calibrate(
#'   sim_fun = wrapper_reduced,
#'   bounds = bounds_reduced,
#'   ...
#' )
#' }
fix_parameters <- function(sim_fun, fixed_params, original_bounds) {

  cat("\nCreating wrapper with fixed parameters:\n")
  for (param in names(fixed_params)) {
    cat(sprintf("  %s = %.3f\n", param, fixed_params[[param]]))
  }
  cat("\n")

  # Create wrapper function
  wrapper <- function(theta, fidelity = c("low", "med", "high"), seed = 2025, ...) {
    # Inject fixed parameters
    theta_full <- c(as.list(theta), as.list(fixed_params))

    # Call original simulator
    result <- sim_fun(theta = theta_full, fidelity = fidelity, seed = seed, ...)

    return(result)
  }

  # Store metadata
  attr(wrapper, "fixed_params") <- fixed_params
  attr(wrapper, "original_sim_fun") <- sim_fun
  attr(wrapper, "reduced_params") <- names(original_bounds)[!names(original_bounds) %in% names(fixed_params)]

  return(wrapper)
}


#' Remove Fixed Parameters from Bounds
#'
#' Helper to remove fixed parameters from bounds list.
#'
#' @param bounds Original bounds
#' @param fixed_params Named vector of fixed parameters
#'
#' @return Reduced bounds list
#' @export
remove_fixed_from_bounds <- function(bounds, fixed_params) {
  bounds[!names(bounds) %in% names(fixed_params)]
}


#' Sequential BO Refinement
#'
#' Convenience function to run coarse → refined BO in one call.
#'
#' @param sim_fun Simulator function
#' @param bounds Initial (wide) bounds
#' @param objective Objective to minimize
#' @param constraints Constraint list
#' @param coarse_budget Budget for coarse phase
#' @param refine_budget Budget for refined phase
#' @param shrink_factor Bound shrink factor (default: 0.5)
#' @param save_intermediate Whether to save intermediate results
#' @param output_dir Directory for intermediate saves
#' @param ... Additional args to bo_calibrate
#'
#' @return List with coarse_fit and refined_fit
#' @export
#'
#' @examples
#' \dontrun{
#' result <- sequential_refinement(
#'   sim_fun = wrapper,
#'   bounds = bounds_wide,
#'   objective = "EN",
#'   constraints = list(power = c("ge", 0.80), type1 = c("le", 0.10)),
#'   coarse_budget = 100,
#'   refine_budget = 50,
#'   save_intermediate = TRUE,
#'   output_dir = "results/sequential"
#' )
#'
#' # Access results
#' coarse_best <- result$coarse_fit$history[result$coarse_fit$history$feasible, ][1, ]
#' refined_best <- result$refined_fit$history[result$refined_fit$history$feasible, ][1, ]
#' }
sequential_refinement <- function(sim_fun,
                                   bounds,
                                   objective,
                                   constraints,
                                   coarse_budget = 100,
                                   refine_budget = 50,
                                   shrink_factor = 0.5,
                                   save_intermediate = TRUE,
                                   output_dir = "results/sequential",
                                   ...) {

  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Sequential BO Refinement\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  # Phase 1: Coarse exploration
  cat("PHASE 1: Coarse Exploration\n")
  cat("─────────────────────────────────────────────────────────\n")
  cat("  Budget:", coarse_budget, "evaluations\n")
  cat("  Bounds: WIDE (full parameter space)\n\n")

  coarse_fit <- evolveBO::bo_calibrate(
    sim_fun = sim_fun,
    bounds = bounds,
    objective = objective,
    constraints = constraints,
    budget = coarse_budget,
    ...
  )

  if (save_intermediate) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    save_bo_state(coarse_fit, file.path(output_dir, "coarse_fit.rds"))
  }

  # Phase 2: Refined optimization
  cat("\n")
  cat("PHASE 2: Refined Optimization\n")
  cat("─────────────────────────────────────────────────────────\n")

  bounds_refined <- refine_bounds(coarse_fit, shrink_factor = shrink_factor)

  cat("  Budget:", refine_budget, "evaluations\n")
  cat("  Bounds: NARROW (focused around best design)\n\n")

  refined_fit <- evolveBO::bo_calibrate(
    sim_fun = sim_fun,
    bounds = bounds_refined,
    objective = objective,
    constraints = constraints,
    n_init = max(10, length(bounds_refined) * 5),  # Fewer initial samples
    budget = refine_budget,
    ...
  )

  if (save_intermediate) {
    save_bo_state(refined_fit, file.path(output_dir, "refined_fit.rds"))
  }

  # Summary
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Summary\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat(sprintf("Total evaluations: %d (coarse) + %d (refined) = %d\n",
              nrow(coarse_fit$history), nrow(refined_fit$history),
              nrow(coarse_fit$history) + nrow(refined_fit$history)))

  coarse_best <- min(coarse_fit$history[[objective]][coarse_fit$history$feasible])
  refined_best <- min(refined_fit$history[[objective]][refined_fit$history$feasible])
  improvement <- 100 * (coarse_best - refined_best) / coarse_best

  cat(sprintf("Best %s (coarse): %.4f\n", objective, coarse_best))
  cat(sprintf("Best %s (refined): %.4f\n", objective, refined_best))
  cat(sprintf("Improvement: %.2f%%\n", improvement))
  cat("═══════════════════════════════════════════════════════════\n\n")

  return(list(
    coarse_fit = coarse_fit,
    refined_fit = refined_fit,
    bounds_refined = bounds_refined,
    improvement_pct = improvement
  ))
}
