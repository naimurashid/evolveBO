#' Fully Adaptive Fidelity Selection (no iteration cutoffs)
#'
#' This is an alternative to select_fidelity_staged() that makes fidelity
#' decisions based ONLY on the current state (uncertainty, feasibility),
#' not on iteration number.
#'
#' @param prob_feasible Probability that candidate satisfies all constraints
#' @param cv_estimate Coefficient of variation for objective prediction
#' @param fidelity_levels Named numeric vector of fidelity levels
#' @return Name of selected fidelity level
#' @keywords internal
select_fidelity_adaptive <- function(prob_feasible, cv_estimate, fidelity_levels) {
  if (length(fidelity_levels) == 1L) {
    return(names(fidelity_levels))
  }

  has_med <- "med" %in% names(fidelity_levels)
  has_high <- "high" %in% names(fidelity_levels)

  # Decision based ONLY on current state, no iteration checks!

  # High fidelity: High uncertainty near constraint boundary
  # This is where precision matters most
  if (cv_estimate > 0.18 && prob_feasible >= 0.2 && prob_feasible <= 0.8) {
    if (has_high) {
      return("high")
    } else if (has_med) {
      return("med")
    }
  }

  # High fidelity: Very promising candidate (likely near optimum)
  if (prob_feasible > 0.7 && cv_estimate > 0.10) {
    if (has_high) {
      return("high")
    } else if (has_med) {
      return("med")
    }
  }

  # Medium fidelity: Moderate uncertainty in promising regions
  if (prob_feasible >= 0.4 && cv_estimate > 0.05) {
    if (has_med) {
      return("med")
    }
  }

  # Medium fidelity: Near boundary with moderate uncertainty
  if (prob_feasible >= 0.15 && prob_feasible <= 0.85 && cv_estimate > 0.08) {
    if (has_med) {
      return("med")
    }
  }

  # Default: Low fidelity for exploration and clearly infeasible regions
  names(fidelity_levels)[1]  # "low"
}
