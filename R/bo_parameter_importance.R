# bo_parameter_importance.R
# Parameter importance analysis for Bayesian optimization results
#
# This module provides functions to identify which parameters are key vs unimportant
# using ARD (Automatic Relevance Determination) from GP lengthscales, local sensitivity
# analysis, and global sensitivity (Sobol indices).

#' Extract lengthscales from GP surrogate models
#'
#' Extracts ARD (Automatic Relevance Determination) lengthscales from DiceKriging
#' or hetGP surrogate models. Lengthscales indicate parameter importance:
#' - Small lengthscale → parameter is important (function varies quickly)
#' - Large lengthscale → parameter is unimportant (function varies slowly)
#'
#' @param surrogates List of surrogate models from bo_calibrate result
#' @param metric Which metric to extract lengthscales from (default: objective metric)
#'
#' @return Named numeric vector of lengthscales for each parameter
#'
#' @export
extract_lengthscales <- function(surrogates, metric = NULL) {
  # If metric not specified, use first surrogate (usually objective)
  if (is.null(metric)) {
    metric <- names(surrogates)[1]
  }

  if (!metric %in% names(surrogates)) {
    stop(sprintf("Metric '%s' not found in surrogates. Available: %s",
                 metric, paste(names(surrogates), collapse = ", ")))
  }

  model <- surrogates[[metric]]

  # Extract lengthscales depending on model type
  if (inherits(model, "km")) {
    # DiceKriging model
    cov_params <- DiceKriging::coef(model)

    # Handle different covariance structures
    if (!is.null(cov_params$range)) {
      # Range parameters are lengthscales
      lengthscales <- cov_params$range
    } else if (!is.null(cov_params$theta)) {
      # Theta parameters (convert to lengthscales: l = 1/sqrt(theta))
      lengthscales <- 1 / sqrt(cov_params$theta)
    } else {
      stop("Could not extract lengthscales from DiceKriging model")
    }

  } else if (inherits(model, "homGP") || inherits(model, "hetGP")) {
    # hetGP model
    lengthscales <- 1 / sqrt(model$theta)

  } else {
    stop(sprintf("Unknown surrogate model type: %s. Expected 'km' or 'hetGP'.",
                 class(model)[1]))
  }

  # Add parameter names if available
  if (!is.null(names(lengthscales))) {
    return(lengthscales)
  } else if (!is.null(model@d)) {
    # DiceKriging uses @d for number of dimensions
    param_names <- paste0("x", seq_len(model@d))
    names(lengthscales) <- param_names
    return(lengthscales)
  } else {
    # Fallback to generic names
    param_names <- paste0("x", seq_along(lengthscales))
    names(lengthscales) <- param_names
    return(lengthscales)
  }
}


#' Analyze parameter importance using ARD lengthscales
#'
#' Computes parameter importance scores by inverting GP lengthscales.
#' Returns a ranked data frame of parameter importance.
#'
#' @param fit Result object from bo_calibrate()
#' @param metric Which metric to analyze (default: objective metric)
#'
#' @return Data frame with columns:
#'   - parameter: Parameter name
#'   - lengthscale: GP lengthscale
#'   - importance: Normalized importance score (higher = more important)
#'   - rank: Importance rank (1 = most important)
#'
#' @export
analyze_parameter_importance <- function(fit, metric = NULL) {
  if (is.null(fit$surrogates)) {
    stop("fit object does not contain surrogates. Did you run bo_calibrate()?")
  }

  # Extract lengthscales
  lengthscales <- extract_lengthscales(fit$surrogates, metric = metric)

  # Compute importance scores (inverse lengthscales)
  importance <- 1 / lengthscales
  importance_normalized <- importance / sum(importance)

  # Rank parameters
  param_ranking <- data.frame(
    parameter = names(lengthscales),
    lengthscale = lengthscales,
    importance = importance_normalized,
    rank = rank(-importance_normalized),
    stringsAsFactors = FALSE
  )

  # Sort by importance (descending)
  param_ranking <- param_ranking[order(-param_ranking$importance), ]
  rownames(param_ranking) <- NULL

  return(param_ranking)
}


#' Compute local sensitivity via numerical gradients
#'
#' Computes gradient-based sensitivity at the best design found by BO.
#' Uses central differences to numerically approximate gradients.
#'
#' @param fit Result object from bo_calibrate()
#' @param epsilon Step size for finite differences (default: 0.01 = 1% of range)
#' @param metric Which metric to analyze (default: objective metric)
#'
#' @return Data frame with columns:
#'   - parameter: Parameter name
#'   - gradient: Numerical gradient at best design
#'   - abs_sensitivity: Absolute gradient magnitude
#'   - normalized_sensitivity: Sensitivity normalized to sum to 1
#'
#' @export
compute_local_sensitivity <- function(fit, epsilon = 0.01, metric = NULL) {
  if (is.null(fit$best_theta)) {
    stop("fit object does not contain best_theta. Did you run bo_calibrate()?")
  }

  if (is.null(fit$surrogates)) {
    stop("fit object does not contain surrogates. Cannot predict gradients.")
  }

  # Determine which metric to analyze
  if (is.null(metric)) {
    metric <- names(fit$surrogates)[1]
  }

  surrogate <- fit$surrogates[[metric]]
  best_theta <- fit$best_theta

  # Get parameter names
  param_names <- names(best_theta)

  # Compute numerical gradients
  gradients <- sapply(param_names, function(param) {
    # Perturb parameter up and down
    theta_plus <- theta_minus <- best_theta
    theta_plus[[param]] <- best_theta[[param]] + epsilon
    theta_minus[[param]] <- best_theta[[param]] - epsilon

    # Convert to matrix for prediction
    x_plus <- matrix(unlist(theta_plus), nrow = 1)
    x_minus <- matrix(unlist(theta_minus), nrow = 1)

    # Predict using surrogate
    if (inherits(surrogate, "km")) {
      y_plus <- DiceKriging::predict(surrogate, newdata = x_plus, type = "UK")$mean
      y_minus <- DiceKriging::predict(surrogate, newdata = x_minus, type = "UK")$mean
    } else if (inherits(surrogate, "homGP") || inherits(surrogate, "hetGP")) {
      y_plus <- predict(surrogate, x = x_plus)$mean
      y_minus <- predict(surrogate, x = x_minus)$mean
    } else {
      stop(sprintf("Unknown surrogate type: %s", class(surrogate)[1]))
    }

    # Central difference
    (y_plus - y_minus) / (2 * epsilon)
  })

  # Absolute sensitivity
  abs_sensitivity <- abs(gradients)

  # Normalized sensitivity
  if (sum(abs_sensitivity) > 0) {
    normalized_sensitivity <- abs_sensitivity / sum(abs_sensitivity)
  } else {
    normalized_sensitivity <- rep(0, length(abs_sensitivity))
  }

  # Build result data frame
  result <- data.frame(
    parameter = param_names,
    gradient = gradients,
    abs_sensitivity = abs_sensitivity,
    normalized_sensitivity = normalized_sensitivity,
    stringsAsFactors = FALSE
  )

  # Sort by absolute sensitivity (descending)
  result <- result[order(-result$abs_sensitivity), ]
  rownames(result) <- NULL

  return(result)
}


#' Compute global sensitivity using Sobol indices
#'
#' Wrapper around compute_sobol_indices() with parameter importance interpretation.
#' Uses variance-based decomposition to identify main effects and interactions.
#'
#' @param fit Result object from bo_calibrate()
#' @param n_samples Number of samples for Sobol estimation (default: 1000)
#' @param metric Which metric to analyze (default: objective metric)
#'
#' @return Data frame with columns:
#'   - parameter: Parameter name
#'   - first_order: First-order Sobol index (main effect)
#'   - total_order: Total Sobol index (including interactions)
#'   - interaction_strength: Difference (total - first_order)
#'
#' @export
compute_global_sensitivity <- function(fit, n_samples = 1000, metric = NULL) {
  if (!requireNamespace("sensitivity", quietly = TRUE)) {
    warning("Package 'sensitivity' not available. Using local sensitivity only.")
    return(NULL)
  }

  # Call existing Sobol computation function
  sobol_result <- compute_sobol_indices(
    fit = fit,
    n_samples = n_samples,
    metric = metric
  )

  # Add interaction strength
  sobol_result$interaction_strength <- sobol_result$total_order - sobol_result$first_order

  return(sobol_result)
}


#' Generate parameter recommendations (REMOVE, FIX, or OPTIMIZE)
#'
#' Analyzes parameter importance using multiple methods and generates
#' actionable recommendations for dimensionality reduction.
#'
#' @param fit Result object from bo_calibrate()
#' @param threshold_irrelevant Importance below this → REMOVE (default: 0.01 = 1%)
#' @param threshold_unimportant Importance below this → FIX (default: 0.05 = 5%)
#' @param confidence Confidence level for consensus (default: 0.8)
#'
#' @return List with class "param_recommendations" containing:
#'   - analysis: Combined importance data frame
#'   - recommendations: List of recommendations per parameter
#'   - thresholds: Applied thresholds
#'
#' @export
generate_recommendations <- function(
  fit,
  threshold_irrelevant = 0.01,
  threshold_unimportant = 0.05,
  confidence = 0.8
) {
  # Combine ARD and local sensitivity
  ard_importance <- analyze_parameter_importance(fit)
  local_sens <- compute_local_sensitivity(fit)

  # Merge results
  combined <- merge(
    ard_importance,
    local_sens[, c("parameter", "normalized_sensitivity")],
    by = "parameter",
    all.x = TRUE
  )

  # Compute consensus importance (handle NAs from merge)
  combined$normalized_sensitivity[is.na(combined$normalized_sensitivity)] <- 0
  combined$consensus_importance <- (combined$importance + combined$normalized_sensitivity) / 2

  # Determine action for each parameter
  combined$action <- ifelse(
    !is.na(combined$consensus_importance) & combined$consensus_importance < threshold_irrelevant, "REMOVE",
    ifelse(!is.na(combined$consensus_importance) & combined$consensus_importance < threshold_unimportant, "FIX", "OPTIMIZE")
  )

  # Sort by consensus importance
  combined <- combined[order(-combined$consensus_importance), ]
  rownames(combined) <- NULL

  # Generate specific recommendations
  recommendations <- list()

  for (i in seq_len(nrow(combined))) {
    param <- combined$parameter[i]
    action <- combined$action[i]
    consensus_imp <- combined$consensus_importance[i]

    if (action == "FIX") {
      # Recommend value from best design
      recommended_value <- fit$best_theta[[param]]

      recommendations[[param]] <- list(
        action = "FIX",
        value = recommended_value,
        reason = sprintf(
          "Low importance (%.1f%%). Fix at best value: %.3f",
          consensus_imp * 100,
          recommended_value
        )
      )

    } else if (action == "REMOVE") {
      recommendations[[param]] <- list(
        action = "REMOVE",
        value = NULL,
        reason = sprintf(
          "Negligible importance (%.1f%%). Safe to remove.",
          consensus_imp * 100
        )
      )

    } else {  # OPTIMIZE
      recommendations[[param]] <- list(
        action = "OPTIMIZE",
        value = NULL,
        reason = sprintf(
          "Important parameter (%.1f%%). Continue optimizing.",
          consensus_imp * 100
        )
      )
    }
  }

  # Return structured result
  result <- structure(
    list(
      analysis = combined,
      recommendations = recommendations,
      thresholds = list(
        irrelevant = threshold_irrelevant,
        unimportant = threshold_unimportant
      )
    ),
    class = "param_recommendations"
  )

  return(result)
}


#' Print method for parameter recommendations
#'
#' @param x Object of class "param_recommendations"
#' @param ... Additional arguments (unused)
#'
#' @export
print.param_recommendations <- function(x, ...) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  Parameter Importance Recommendations\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

  # Group by action
  remove_params <- names(x$recommendations)[sapply(x$recommendations, function(r) r$action == "REMOVE")]
  fix_params <- names(x$recommendations)[sapply(x$recommendations, function(r) r$action == "FIX")]
  optimize_params <- names(x$recommendations)[sapply(x$recommendations, function(r) r$action == "OPTIMIZE")]

  if (length(remove_params) > 0) {
    cat("REMOVE (negligible importance):\n")
    for (param in remove_params) {
      cat(sprintf("  ✗ %s: %s\n", param, x$recommendations[[param]]$reason))
    }
    cat("\n")
  }

  if (length(fix_params) > 0) {
    cat("FIX (low importance):\n")
    for (param in fix_params) {
      cat(sprintf("  ⚙ %s: %s\n", param, x$recommendations[[param]]$reason))
    }
    cat("\n")
  }

  if (length(optimize_params) > 0) {
    cat("OPTIMIZE (key parameters):\n")
    for (param in optimize_params) {
      cat(sprintf("  ✓ %s: %s\n", param, x$recommendations[[param]]$reason))
    }
    cat("\n")
  }

  cat("───────────────────────────────────────────────────────────\n")
  cat(sprintf("Dimensionality reduction: %dD → %dD\n",
              length(x$recommendations),
              length(optimize_params)))
  cat("═══════════════════════════════════════════════════════════\n\n")

  invisible(x)
}
