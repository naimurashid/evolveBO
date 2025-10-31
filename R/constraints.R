#' Normalise user-supplied constraint specification
#' @keywords internal
parse_constraints <- function(constraints) {
  if (is.null(constraints) || length(constraints) == 0) {
    return(tibble::tibble(metric = character(), direction = character(), threshold = numeric()))
  }
  if (!is.list(constraints)) {
    stop("`constraints` must be a named list.", call. = FALSE)
  }
  tibble::tibble(
    metric = names(constraints),
    direction = purrr::map_chr(constraints, ~ {
      dir <- .x[[1]]
      if (!dir %in% c("ge", "le")) {
        stop("Constraint directions must be 'ge' or 'le'.", call. = FALSE)
      }
      dir
    }),
    threshold = purrr::map_dbl(constraints, ~ as.numeric(.x[[2]]))
  )
}

#' Deterministic feasibility check for observed metrics
#' @keywords internal
is_feasible <- function(metrics, constraint_tbl) {
  if (nrow(constraint_tbl) == 0L) {
    return(TRUE)
  }
  purrr::pmap_lgl(constraint_tbl, function(metric, direction, threshold) {
    value <- metrics[[metric]]
    if (is.null(value) || is.na(value)) {
      return(FALSE)
    }
    if (direction == "ge") {
      value >= threshold
    } else {
      value <= threshold
    }
  }) |>
    all()
}

#' Probability of feasibility under Gaussian predictive distribution
#' @keywords internal
prob_feasibility <- function(mean, sd, constraint_tbl) {
  if (nrow(constraint_tbl) == 0L) {
    return(1)
  }
  purrr::pmap_dbl(constraint_tbl, function(metric, direction, threshold) {
    mu <- mean[[metric]]
    sigma <- sd[[metric]]
    if (is.null(mu) || is.na(mu) || is.null(sigma) || is.na(sigma) || sigma <= 0) {
      return(0)
    }
    if (direction == "ge") {
      stats::pnorm((mu - threshold) / sigma, lower.tail = TRUE)
    } else {
      stats::pnorm((threshold - mu) / sigma, lower.tail = TRUE)
    }
  }) |>
    prod()
}
