#' Internal helpers for parameter transformations
#' @keywords internal
scale_to_unit <- function(theta, bounds) {
  # Use base R instead of purrr::imap_dbl to avoid "In index: X" errors
  theta_names <- names(theta)
  result <- vapply(theta_names, function(name) {
    value <- as.numeric(theta[[name]])
    range <- bounds[[name]]
    if (length(range) != 2L) {
      stop(sprintf("Bounds for '%s' must have length 2.", name), call. = FALSE)
    }
    (value - range[[1]]) / (range[[2]] - range[[1]])
  }, FUN.VALUE = numeric(1))
  names(result) <- theta_names
  result
}

#' @keywords internal
scale_from_unit <- function(unit_x, bounds) {
  # Use base R instead of purrr::imap to avoid "In index: X" errors
  param_names <- names(unit_x)
  result <- lapply(param_names, function(name) {
    value <- as.numeric(unit_x[[name]])
    range <- bounds[[name]]
    range[[1]] + value * (range[[2]] - range[[1]])
  })
  names(result) <- param_names
  result
}

#' Sample an initial Latin hypercube design
#' @keywords internal
lhs_design <- function(n, bounds, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  d <- length(bounds)
  lhs <- lhs::randomLHS(n, d)
  theta_names <- names(bounds)
  # Use base R lapply instead of purrr::map
  lapply(seq_len(n), function(i) {
    unit_list <- as.list(lhs[i, , drop = TRUE])
    names(unit_list) <- theta_names
    scaled <- scale_from_unit(unit_list, bounds)
    names(scaled) <- theta_names
    scaled
  })
}

#' Convert theta list to named numeric vector
#' @keywords internal
theta_list_to_vec <- function(theta) {
  # Use base R vapply instead of purrr::map_dbl
  result <- vapply(theta, as.numeric, FUN.VALUE = numeric(1))
  names(result) <- names(theta)
  result
}

#' Create deterministic identifier for theta in unit space
#' @keywords internal
theta_to_id <- function(unit_theta) {
  paste(format(round(as.numeric(unit_theta), digits = 6), nsmall = 6), collapse = "|")
}

#' Enforce integer parameters if requested
#' @keywords internal
coerce_theta_types <- function(theta, integer_params = NULL) {
  if (is.null(integer_params) || length(integer_params) == 0L) {
    return(theta)
  }
  # Use base R lapply instead of purrr::imap
  theta_names <- names(theta)
  result <- lapply(theta_names, function(name) {
    value <- theta[[name]]
    if (name %in% integer_params) {
      value <- round(as.numeric(value))
    }
    value
  })
  names(result) <- theta_names
  result
}

#' Ensure bounds are valid numeric intervals
#' @keywords internal
validate_bounds <- function(bounds) {
  if (!is.list(bounds)) {
    stop("`bounds` must be a named list.", call. = FALSE)
  }
  if (is.null(names(bounds)) || any(names(bounds) == "")) {
    stop("`bounds` must have names.", call. = FALSE)
  }
  # Use base R for loop instead of purrr::walk
  for (rng in bounds) {
    if (length(rng) != 2L) {
      stop("Each bound entry must have length 2.", call. = FALSE)
    }
    if (!is.numeric(rng)) {
      stop("Bounds must be numeric.", call. = FALSE)
    }
    if (rng[[1]] >= rng[[2]]) {
      stop("Lower bound must be < upper bound.", call. = FALSE)
    }
  }
  invisible(bounds)
}

#' Null-coalescing helper
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}
