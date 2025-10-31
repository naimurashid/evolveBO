#' Fit Gaussian process surrogates for operating characteristics
#'
#' @param history tibble with columns `unit_x` (list of numeric vectors),
#'   `metrics` (list of named numeric vectors), and optional `variance`
#'   (list of named numeric vectors with noise variances).
#' @param objective name of the objective metric.
#' @param constraint_tbl tibble produced by [parse_constraints()].
#' @param covtype covariance kernel used by DiceKriging (default `"matern5_2"`).
#' @export
fit_surrogates <- function(history,
                           objective,
                           constraint_tbl,
                           covtype = "matern5_2") {
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

  surrogates <- purrr::map(metrics_needed, function(metric) {
    values <- purrr::map_dbl(history$metrics, ~ as.numeric(.x[[metric]]))
    noise <- purrr::map_dbl(history$variance, function(var_list) {
      if (is.null(var_list)) return(NA_real_)
      as.numeric(var_list[[metric]])
    })
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
      do.call(what = rbind)
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

    tryCatch({
      if (is.null(noise_vec)) {
        DiceKriging::km(
          design = X_unique,
          response = aggr$value,
          covtype = covtype,
          nugget = nugget,
          nugget.estim = FALSE,
          control = list(trace = FALSE)
        )
      } else {
        DiceKriging::km(
          design = X_unique,
          response = aggr$value,
          covtype = covtype,
          noise.var = noise_vec,
          nugget.estim = FALSE,
          control = list(trace = FALSE)
        )
      }
    }, error = function(e) {
      stop(sprintf("Failed to fit surrogate for metric '%s': %s\nThis may indicate ill-conditioned data or insufficient observations.",
                   metric, e$message), call. = FALSE)
    })
  })
  names(surrogates) <- metrics_needed
  surrogates
}

#' Predict metrics from fitted surrogates
#' @keywords internal
predict_surrogates <- function(surrogates, unit_x) {
  if (length(surrogates) == 0L) {
    stop("No surrogate models available for prediction.", call. = FALSE)
  }
  param_names <- colnames(surrogates[[1]]@X)
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
      pred <- DiceKriging::predict.km(.x, newdata = design_mat, type = "UK", se.compute = TRUE, cov.compute = FALSE)
      list(mean = as.numeric(pred$mean), sd = sqrt(pmax(pred$sd^2, 0)), model = .x)
    }
  )
}
