#' Internal helpers for parameter transformations
#' @keywords internal
scale_to_unit <- function(theta, bounds) {
  purrr::imap_dbl(theta, function(value, name) {
    value <- as.numeric(value)
    range <- bounds[[name]]
    if (length(range) != 2L) {
      stop(sprintf("Bounds for '%s' must have length 2.", name), call. = FALSE)
    }
    (value - range[[1]]) / (range[[2]] - range[[1]])
  })
}

#' @keywords internal
scale_from_unit <- function(unit_x, bounds) {
  purrr::imap(unit_x, function(value, name) {
    value <- as.numeric(value)
    range <- bounds[[name]]
    range[[1]] + value * (range[[2]] - range[[1]])
  })
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
  purrr::map(
    seq_len(n),
    ~ scale_from_unit(as.list(lhs[.x, , drop = TRUE]), bounds) |>
      purrr::set_names(theta_names)
  )
}

#' Convert theta list to named numeric vector
#' @keywords internal
theta_list_to_vec <- function(theta) {
  purrr::map_dbl(theta, as.numeric)
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
  purrr::imap(theta, function(value, name) {
    if (name %in% integer_params) {
      value <- round(as.numeric(value))
    }
    value
  })
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
  purrr::walk(bounds, function(rng) {
    if (length(rng) != 2L) {
      stop("Each bound entry must have length 2.", call. = FALSE)
    }
    if (!is.numeric(rng)) {
      stop("Bounds must be numeric.", call. = FALSE)
    }
    if (rng[[1]] >= rng[[2]]) {
      stop("Lower bound must be < upper bound.", call. = FALSE)
    }
  })
  invisible(bounds)
}
#' Null-coalescing helper
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}
