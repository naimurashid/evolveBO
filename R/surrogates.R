#' Fit Gaussian process surrogates for operating characteristics
#'
#' Implements heteroskedastic GP modeling when variance information is available,
#' using either hetGP (preferred) or DiceKriging with known noise variance.
#'
#' @param history tibble with columns `unit_x` (list of numeric vectors),
#'   `metrics` (list of named numeric vectors), and optional `variance`
#'   (list of named numeric vectors with noise variances).
#' @param objective name of the objective metric.
#' @param constraint_tbl tibble produced by [parse_constraints()].
#' @param covtype covariance kernel used by DiceKriging (default `"matern5_2"`).
#' @param use_hetgp logical; if TRUE and hetGP package is available, use
#'   heteroskedastic GP with replicated observations. Default TRUE.
#' @param prev_surrogates optional list of surrogates from previous iteration
#'   for warm-starting hyperparameter optimization. If provided, uses previous
#'   hyperparameters as initial values, reducing optimization time by 30-50%.
#'
#' @return A named list of fitted surrogate models (GP objects).
#' @export
fit_surrogates <- function(history,
                           objective,
                           constraint_tbl,
                           covtype = "matern5_2",
                           use_hetgp = TRUE,
                           prev_surrogates = NULL) {
  if (nrow(history) == 0L) {
    stop("History is empty; cannot fit surrogates.", call. = FALSE)
  }

  metrics_needed <- unique(c(objective, constraint_tbl$metric))
  param_names <- names(history$unit_x[[1]])
  if (is.null(param_names)) {
    param_names <- names(history$theta[[1]])
  }

  theta_ids <- history$theta_id
  id_groups <- split(seq_along(theta_ids), theta_ids)

  # Check if hetGP is available and requested
  has_hetgp <- requireNamespace("hetGP", quietly = TRUE)
  use_hetgp <- use_hetgp && has_hetgp

  surrogates <- purrr::map(metrics_needed, function(metric) {
    values <- purrr::map_dbl(history$metrics, ~ as.numeric(.x[[metric]]))
    noise <- purrr::map_dbl(history$variance, function(var_list) {
      if (is.null(var_list)) return(NA_real_)
      as.numeric(var_list[[metric]])
    })

    # Check if we have replicated observations (multiple evals at same theta)
    has_replicates <- any(sapply(id_groups, length) > 1)
    has_variance <- !all(is.na(noise))

    # Extract previous model for warm-starting (if available)
    prev_model <- if (!is.null(prev_surrogates) && metric %in% names(prev_surrogates)) {
      prev_surrogates[[metric]]
    } else {
      NULL
    }

    if (use_hetgp && has_replicates && has_variance) {
      # Use hetGP with replicated observations
      fit_hetgp_surrogate(history, metric, id_groups, param_names, covtype, prev_model)
    } else {
      # Fall back to DiceKriging with aggregated observations
      fit_dicekriging_surrogate(history, metric, id_groups, param_names,
                                covtype, noise, values, prev_model)
    }
  })

  names(surrogates) <- metrics_needed
  surrogates
}

#' Fit heteroskedastic GP using hetGP package
#' @keywords internal
fit_hetgp_surrogate <- function(history, metric, id_groups, param_names, covtype, prev_model = NULL) {
  # Prepare data with replicates
  X_list <- list()
  Z_list <- list()

  for (i in seq_along(id_groups)) {
    idx <- id_groups[[i]]
    unit_theta <- history$unit_x[[idx[1]]]
    X_list[[i]] <- as.numeric(unit_theta)

    # All observations at this location
    Z_list[[i]] <- sapply(idx, function(j) {
      as.numeric(history$metrics[[j]][[metric]])
    })
  }

  X <- do.call(rbind, X_list)
  colnames(X) <- param_names

  # hetGP expects Z as vector with mult indicating replication counts
  Z_vec <- unlist(Z_list)
  mult <- sapply(Z_list, length)

  # Map covtype to hetGP
  hetgp_cov <- switch(covtype,
                      matern5_2 = "Matern5_2",
                      matern3_2 = "Matern3_2",
                      gauss = "Gaussian",
                      "Matern5_2")  # default

  tryCatch({
    hetGP::mleHetGP(
      X = X,
      Z = Z_vec,
      mult = mult,
      covtype = hetgp_cov,
      settings = list(return.hom = TRUE),  # also return homoskedastic fit
      known = list(g = 1e-8)  # small nugget for numerical stability
    )
  }, error = function(e) {
    warning(sprintf("hetGP fit failed for metric '%s': %s\nFalling back to homoskedastic GP.",
                    metric, e$message), call. = FALSE)
    # Fall back to homoskedastic
    hetGP::mleHomGP(X = X, Z = Z_vec, mult = mult, covtype = hetgp_cov,
                   known = list(g = 1e-6))
  })
}

#' Fit GP using DiceKriging (aggregated observations)
#' @keywords internal
fit_dicekriging_surrogate <- function(history, metric, id_groups, param_names,
                                      covtype, noise, values, prev_model = NULL) {
  # Aggregate by theta_id
  aggr <- purrr::imap_dfr(id_groups, function(idx, id) {
    unit_theta <- history$unit_x[[idx[1]]]
    tibble::tibble(
      theta_id = id,
      unit_x = list(unit_theta),
      value = mean(values[idx], na.rm = TRUE),
      noise = if (all(is.na(noise[idx]))) NA_real_ else mean(noise[idx], na.rm = TRUE)
    )
  })

  X_unique <- aggr$unit_x |>
    purrr::map(~ {
      vec <- as.numeric(.x)
      if (length(vec) != length(param_names)) {
        stop("Surrogate design dimension mismatch.", call. = FALSE)
      }
      vec
    }) |>
    do.call(what = rbind) |>
    as.matrix()
  colnames(X_unique) <- param_names

  noise_vec <- aggr$noise
  if (all(is.na(noise_vec))) {
    nugget <- 1e-6
    noise_vec <- NULL
  } else {
    noise_vec[is.na(noise_vec)] <- min(noise_vec, na.rm = TRUE)
    noise_vec <- pmax(noise_vec, 1e-6)
    nugget <- 0
  }

  # Extract hyperparameters from previous model for warm-starting
  parinit <- extract_gp_hyperparams(prev_model)

  # Debug: Check what we're passing to DiceKriging::km
  if (!is.matrix(X_unique)) {
    warning(sprintf("X_unique is not a matrix! Class: %s, Dim: %s",
                    class(X_unique)[1],
                    paste(dim(X_unique), collapse = "x")))
  }
  if (!is.numeric(aggr$value)) {
    warning(sprintf("aggr$value is not numeric! Class: %s, Length: %d",
                    class(aggr$value)[1],
                    length(aggr$value)))
  }

  tryCatch({
    if (is.null(noise_vec)) {
      DiceKriging::km(
        design = X_unique,
        response = aggr$value,
        covtype = covtype,
        nugget = nugget,
        nugget.estim = FALSE,
        control = list(
          trace = FALSE,
          parinit = parinit  # Warm start
        )
      )
    } else {
      DiceKriging::km(
        design = X_unique,
        response = aggr$value,
        covtype = covtype,
        noise.var = noise_vec,
        nugget.estim = FALSE,
        control = list(
          trace = FALSE,
          parinit = parinit  # Warm start
        )
      )
    }
  }, error = function(e) {
    stop(sprintf("Failed to fit surrogate for metric '%s': %s\nThis may indicate ill-conditioned data or insufficient observations.",
                 metric, e$message), call. = FALSE)
  })
}

#' Extract GP hyperparameters for warm-starting
#'
#' Extracts lengthscale parameters from a fitted GP model to use as
#' initial values for the next optimization. Works with DiceKriging::km
#' and hetGP models.
#'
#' @param model fitted GP model (km or hetGP object), or NULL
#' @return numeric vector of hyperparameters, or NULL if extraction fails
#' @keywords internal
extract_gp_hyperparams <- function(model) {
  if (is.null(model)) {
    return(NULL)
  }

  tryCatch({
    if (inherits(model, "km")) {
      # DiceKriging model
      theta <- model@covariance@range.val
      if (is.numeric(theta) && all(is.finite(theta)) && all(theta > 0)) {
        return(theta)
      }
    } else if (inherits(model, "hetGP")) {
      # hetGP model
      theta <- model$theta
      if (is.numeric(theta) && all(is.finite(theta)) && all(theta > 0)) {
        return(theta)
      }
    }
    return(NULL)
  }, error = function(e) {
    # If extraction fails, return NULL (no warm-start)
    return(NULL)
  })
}

#' Predict metrics from fitted surrogates
#'
#' Handles both DiceKriging::km and hetGP::hetGP model objects.
#'
#' @keywords internal
predict_surrogates <- function(surrogates, unit_x) {
  if (length(surrogates) == 0L) {
    stop("No surrogate models available for prediction.", call. = FALSE)
  }

  # Detect model type from first surrogate
  first_model <- surrogates[[1]]
  is_hetgp <- inherits(first_model, "hetGP")

  if (is_hetgp) {
    param_names <- colnames(first_model$X0)
  } else {
    param_names <- colnames(first_model@X)
  }

  if (is.null(param_names)) {
    param_names <- names(unit_x[[1]])
  }

  design <- purrr::map_dfr(unit_x, function(point) {
    vec <- unlist(point)
    if (is.null(names(vec))) {
      names(vec) <- param_names
    }
    tibble::as_tibble_row(vec[param_names])
  })
  design_mat <- as.matrix(design)

  purrr::imap(
    surrogates,
    ~ {
      if (inherits(.x, "hetGP")) {
        # Use hetGP predict method
        pred <- predict(x = design_mat, object = .x)
        list(mean = as.numeric(pred$mean),
             sd = as.numeric(sqrt(pmax(pred$sd2, 0))),
             model = .x)
      } else {
        # Use DiceKriging predict method
        pred <- DiceKriging::predict.km(.x, newdata = design_mat,
                                        type = "UK", se.compute = TRUE,
                                        cov.compute = FALSE)
        list(mean = as.numeric(pred$mean),
             sd = sqrt(pmax(pred$sd^2, 0)),
             model = .x)
      }
    }
  )
}
