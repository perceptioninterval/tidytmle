#' Modified Treatment Policy Constructors for Continuous Exposures
#'
#' Define causal Modified Treatment Policies (MTPs) for continuous treatments
#' following the Longitudinal Modified Treatment Policy (LMTP) framework
#' (Díaz, Williams, Hoffman, & Schenck, 2021).
#'
#' @details
#' Unlike static binary interventions ("treat everyone vs. treat no one"),
#' continuous exposures often violate positivity if we try to force every
#' individual to the exact same dose. Modified Treatment Policies avoid
#' extrapolation by shifting each person's observed exposure realistically:
#' - `shift(delta)` defines an additive shift \eqn{d(a, w) = a + \delta}
#'   (e.g., increasing daily exercise by +15 minutes or reducing dose by -1 unit).
#' - `scale_tx(factor)` defines a multiplicative policy
#'   \eqn{d(a, w) = a \cdot \text{factor}} (e.g., reducing exposure by 10%).
#'
#' `tidytmle` estimates the post-intervention density ratio
#' \eqn{g_d(a \mid w) / g(a \mid w)} nonparametrically via stacked
#' cross-validated classification without parametric normality assumptions.
#'
#' @references
#' Díaz, I., Williams, N., Hoffman, K. L., & Schenck, E. J. (2021).
#' Nonparametric causal effects based on longitudinal modified treatment
#' policies. *Journal of the American Statistical Association*, 118(542),
#' 846–857.
#'
#' @param delta Numeric scalar specifying the additive shift applied to observed
#'   treatment values.
#' @param factor Strictly positive numeric scalar specifying the multiplicative
#'   scaling factor applied to observed treatment values.
#'
#' @return An S3 policy object inheriting from `"tmle_policy"` (`"tmle_shift"`
#'   or `"tmle_scale"`).
#'
#' @examples
#' spec_shift <- tmle_spec(
#'   outcome = y,
#'   treatment = dose,
#'   covariates = c(age, bmi),
#'   estimand = shift(1.0)
#' )
#'
#' spec_scale <- tmle_spec(
#'   outcome = y,
#'   treatment = dose,
#'   covariates = c(age, bmi),
#'   estimand = scale_tx(0.9)
#' )
#'
#' @export
shift <- function(delta) {
  if (!is.numeric(delta) || length(delta) != 1L || is.na(delta) || !is.finite(delta)) {
    cli::cli_abort("{.arg delta} must be a single finite numeric value.")
  }
  structure(
    list(
      type = "shift",
      delta = as.numeric(delta),
      label = sprintf("shift(%s)", format(delta, trim = TRUE))
    ),
    class = c("tmle_shift", "tmle_policy")
  )
}

#' @rdname shift
#' @export
scale_tx <- function(factor) {
  if (!is.numeric(factor) || length(factor) != 1L || is.na(factor) ||
      !is.finite(factor) || factor <= 0) {
    cli::cli_abort("{.arg factor} must be a single finite positive numeric value.")
  }
  structure(
    list(
      type = "scale",
      factor = as.numeric(factor),
      label = sprintf("scale_tx(%s)", format(factor, trim = TRUE))
    ),
    class = c("tmle_scale", "tmle_policy")
  )
}

#' @export
print.tmle_policy <- function(x, ...) {
  cat(x$label, "\n", sep = "")
  invisible(x)
}

#' Apply a Modified Treatment Policy to a treatment vector
#' @noRd
apply_policy <- function(policy, a) {
  if (inherits(policy, "tmle_shift")) {
    return(a + policy$delta)
  }
  if (inherits(policy, "tmle_scale")) {
    return(a * policy$factor)
  }
  cli::cli_abort("Unknown {.cls tmle_policy} object of class {.cls {class(policy)}}.")
}

#' Validate and normalize requested estimands
#' @noRd
normalize_estimand <- function(estimand) {
  valid_binary <- c("ate", "att", "atc", "rr", "or", "ey1", "ey0")

  if (inherits(estimand, "tmle_policy")) {
    return(list(
      kind = "policy",
      policies = list(estimand),
      names = estimand$label
    ))
  }

  if (is.list(estimand) && all(vapply(estimand, inherits, logical(1L), "tmle_policy"))) {
    if (length(estimand) == 0L) {
      cli::cli_abort("{.arg estimand} policy list cannot be empty.")
    }
    if (length(estimand) > 1L) {
      cli::cli_abort(c(
        "A single {.fn tmle_spec} fits one continuous Modified Treatment Policy at a time (received {length(estimand)} policies).",
        "i" = "Each continuous policy requires its own density-ratio model. Fit each policy separately and combine with {.code dplyr::bind_rows(tidy(fit1), tidy(fit2))}."
      ))
    }
    p1 <- estimand[[1L]]
    return(list(
      kind = "policy",
      policies = list(p1),
      names = p1$label
    ))
  }

  if (!is.character(estimand) || length(estimand) == 0L || anyNA(estimand)) {
    cli::cli_abort(c(
      "{.arg estimand} must be a non-empty character vector or a {.cls tmle_policy} object.",
      "i" = "Supported string estimands: {.val {valid_binary}}.",
      "i" = "Supported policy estimands: {.fn shift} and {.fn scale_tx}."
    ))
  }

  est_lower <- tolower(trimws(estimand))
  bad <- setdiff(est_lower, valid_binary)
  if (length(bad) > 0L) {
    cli::cli_abort(c(
      "Unsupported estimand{?s}: {.val {bad}}.",
      "i" = "Supported estimands are {.val {valid_binary}}, {.fn shift}, or {.fn scale_tx}."
    ))
  }

  est_unique <- unique(est_lower)
  # Ensure ATT/ATC aren't mixed ambiguously with incompatible policy modes;
  # multiple standard binary estimands (e.g. c("ate", "rr", "or")) are jointly supported!
  list(
    kind = "binary",
    names = est_unique
  )
}
