#' Modified Treatment Policy Constructors for Continuous Exposures
#'
#' Define causal Modified Treatment Policies (MTPs) for continuous treatments
#' following the Longitudinal Modified Treatment Policy (LMTP) framework.
#'
#' - `shift(delta)` defines a deterministic additive intervention
#'   \eqn{d(a, w) = a + \delta}.
#' - `scale_tx(factor)` defines a deterministic multiplicative intervention
#'   \eqn{d(a, w) = a \cdot \text{factor}}.
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
#' @keywords internal
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
#' @keywords internal
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
    labels <- vapply(estimand, function(p) p$label, character(1L))
    return(list(
      kind = "policy",
      policies = estimand,
      names = labels
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
