#' Extract Tidy Causal Inference Table from a TMLE Fit
#'
#' Returns a publication-ready inference tibble containing point estimates,
#' influence-function standard errors, test statistics, p-values, confidence
#' intervals, optional simultaneous uniform confidence bands, and optional
#' VanderWeele-Ding E-values.
#'
#' @param x A fitted `"tmle_fit"` object.
#' @param conf.int Logical; whether to include `conf.low` and `conf.high`
#'   columns (default `TRUE`).
#' @param conf.level Confidence level for intervals (default `0.95`).
#' @param simultaneous Logical; if `TRUE` and multiple estimands are present,
#'   computes simultaneous uniform confidence bands via a fast Rademacher wild
#'   multiplier bootstrap over the joint EIF matrix.
#' @param e_value Logical; if `TRUE`, appends `e_value` and `e_value_ci`
#'   sensitivity columns (VanderWeele & Ding, 2017) quantifying robustness to
#'   unmeasured confounding.
#' @param ... Currently unused.
#'
#' @return A [tibble::tibble()].
#' @method tidy tmle_fit
#' @export
tidy.tmle_fit <- function(x,
                          conf.int = TRUE,
                          conf.level = 0.95,
                          simultaneous = FALSE,
                          e_value = FALSE,
                          ...) {
  if (!is.numeric(conf.level) || length(conf.level) != 1L || conf.level <= 0 || conf.level >= 1) {
    cli::cli_abort("{.arg conf.level} must be strictly between 0 and 1.")
  }

  # Recompute inference on cached targeting outputs if custom conf.level, simultaneous, or e_value requested
  if (abs(conf.level - 0.95) > 1e-12 || isTRUE(simultaneous) || isTRUE(e_value)) {
    idx <- x$subgroup_idx
    cache <- x$cache
    y_sub <- cache$y_num[idx]
    a_sub <- cache$a_num[idx]
    d_sub <- cache$delta[idx]
    y_bounds_sub <- if (cache$outcome_type == "continuous") {
      range(y_sub[d_sub == 1L], na.rm = TRUE)
    } else {
      c(0, 1)
    }
    res <- run_tmle_targeting(
      y_num = y_sub,
      a_num = a_sub,
      delta = d_sub,
      q0_init = cache$q0_init[idx],
      q1_init = cache$q1_init[idx],
      qa_init = cache$qa_init[idx],
      g_hat = cache$g_hat[idx],
      c1_hat = cache$c1_hat[idx],
      c0_hat = cache$c0_hat[idx],
      weights = cache$weights[idx],
      estimand = x$spec$estimand,
      outcome_type = cache$outcome_type,
      y_bounds = y_bounds_sub,
      g_bounds_arg = x$spec$targeting$g_bounds,
      fluctuation = x$spec$targeting$fluctuation,
      cluster_id = if (!is.null(cache$id_vec)) cache$id_vec[idx] else NULL,
      conf.level = conf.level,
      simultaneous = simultaneous,
      e_value = e_value
    )
    tbl <- res$estimates
  } else {
    tbl <- x$targeting$estimates
  }

  if (!isTRUE(conf.int)) {
    tbl <- dplyr::select(tbl, -dplyr::any_of(c("conf.low", "conf.high")))
  }
  tbl
}

#' Extract 1-Row Causal Health Check from a TMLE Fit
#'
#' Summarizes sample sizes, cross-fitting metadata, targeting convergence
#' (`epsilon`, `mean_eif`), plug-in vs. TMLE shift, propensity positivity bounds
#' and truncation rate, Kish effective sample sizes (`ess_treated`,
#' `ess_control`), maximum single-observation EIF leverage, and Super Learner
#' nuisance performance.
#'
#' @details
#' ## How to Read the Causal Health Check
#' - **`mean_eif`**: The sample average of the Efficient Influence Function
#'   \eqn{\frac{1}{n}\sum D_i^*}. A healthy TMLE solves the EIF score equation,
#'   so `abs(mean_eif)` must be near zero (`< 1e-4`).
#' - **`epsilon`**: The magnitude of the targeting fluctuation step. Values
#'   close to `0` mean the initial ML outcome fit already had little residual
#'   confounding bias; larger magnitudes show how strongly TMLE adjusted the
#'   initial predictions to remove bias.
#' - **`plugin_vs_tmle_shift`**: Difference between the final targeted TMLE ATE
#'   and the untargeted G-computation plug-in ATE. Quantifies the exact causal
#'   bias correction performed by the targeting step.
#' - **`pct_truncated` & `ess_treated`/`ess_control`**: Share of propensity
#'   scores clipped to `[g_bound_low, g_bound_high]` and Kish's effective sample
#'   size under inverse-propensity weighting. High truncation (`> 5%`) or tiny
#'   ESS signals near-positivity violations (Petersen et al., 2012).
#' - **`max_eif_leverage`**: Largest single observation's share of total EIF
#'   variance. Values `> 0.10` flag an influential observation dominating the
#'   standard error.
#'
#' @references
#' Petersen, M. L., Porter, K. E., Gruber, S., Wang, Y., & van der Laan, M. J.
#' (2012). Diagnosing and responding to violations in the positivity assumption.
#' *Statistical Methods in Medical Research*, 21(1), 31–54.
#'
#' Schuler, M. S., & Rose, S. (2017). Targeted maximum likelihood estimation for
#' causal inference in observational studies. *American Journal of Epidemiology*,
#' 185(1), 65–73.
#'
#' @param x A fitted `"tmle_fit"` object.
#' @param ... Currently unused.
#'
#' @return A 1-row [tibble::tibble()].
#' @method glance tmle_fit
#' @export
glance.tmle_fit <- function(x, ...) {
  idx <- x$subgroup_idx
  t_res <- x$targeting
  cache <- x$cache

  n_sub <- length(idx)
  n_obs <- sum(cache$delta[idx] == 1L)
  n_clusters <- if (!is.null(cache$id_vec)) length(unique(cache$id_vec[idx])) else NA_integer_

  q_mets <- x$nuisance$Q$metrics
  q_sl <- dplyr::filter(q_mets, learner == "SuperLearner")
  q_metric_val <- if (cache$outcome_type == "binary") {
    mean(q_sl$.estimate[q_sl$.metric == "roc_auc"], na.rm = TRUE)
  } else {
    mean(q_sl$.estimate[q_sl$.metric == "rsq"], na.rm = TRUE)
  }

  g_mets <- x$nuisance$g$metrics
  g_sl <- dplyr::filter(g_mets, learner == "SuperLearner")
  g_auc <- mean(g_sl$.estimate[g_sl$.metric == "roc_auc"], na.rm = TRUE)
  g_ll <- mean(g_sl$.estimate[g_sl$.metric == "mn_log_loss"], na.rm = TRUE)

  tibble::tibble(
    n = n_sub,
    n_observed = n_obs,
    n_clusters = n_clusters,
    v_folds = cache$v_folds,
    cv_tmle = x$spec$targeting$cv_tmle,
    epsilon = t_res$epsilon,
    mean_eif = t_res$mean_eif,
    plugin_vs_tmle_shift = t_res$plugin_vs_tmle_shift,
    g_min = t_res$g_info$g_min,
    g_max = t_res$g_info$g_max,
    g_bound_low = t_res$g_info$bounds[1L],
    g_bound_high = t_res$g_info$bounds[2L],
    pct_truncated = t_res$g_info$pct_truncated,
    ess_treated = t_res$ess_treated,
    ess_control = t_res$ess_control,
    max_eif_leverage = max(t_res$eif_leverage),
    q_metric = q_metric_val,
    g_auc = g_auc,
    g_log_loss = g_ll
  )
}

#' Augment Cohort Data with TMLE Nuisance Predictions and EIF Diagnostics
#'
#' Appends initial and targeted counterfactual predictions (`.q0_init`,
#' `.q1_init`, `.q0_star`, `.q1_star`), raw and truncated propensity scores
#' (`.g_hat`, `.g_trunc`), clever covariates (`.clever_h`), targeted individual
#' treatment effects (`.ite_star`), Efficient Influence Function values
#' (`.eif`), and fractional EIF variance leverage (`.eif_leverage`).
#'
#' @details
#' ## Column Decoder Ring
#' - **`.q0_init`, `.q1_init`**: Initial cross-fitted outcome predictions
#'   \eqn{\bar{Q}^0(0, W)} and \eqn{\bar{Q}^0(1, W)} before targeting
#'   (G-computation predictions).
#' - **`.q0_star`, `.q1_star`**: Targeted TMLE counterfactual predictions
#'   \eqn{\bar{Q}^*(0, W)} and \eqn{\bar{Q}^*(1, W)} after fluctuating along
#'   the Clever Covariate.
#' - **`.g_hat`, `.g_trunc`**: Raw cross-fitted propensity scores
#'   \eqn{P(A = 1 \mid W)} and scores truncated to `g_bounds`.
#' - **`.clever_h`**: The Clever Covariate \eqn{H_n(A, W)} linking propensity
#'   and censoring weights to the parameter's Efficient Influence Function
#'   (Luque-Fernandez et al., 2018).
#' - **`.ite_star`**: Individual targeted treatment contrast
#'   \eqn{\bar{Q}^*(1, W_i) - \bar{Q}^*(0, W_i)}.
#' - **`.eif`, `.eif_leverage`**: Observation-level Efficient Influence Function
#'   contribution \eqn{D_i^*} and variance share \eqn{{D_i^*}^2 / \sum {D_j^*}^2}.
#'
#' @references
#' Luque-Fernandez, M. A., Schomaker, M., Rachet, B., & Schnitzer, M. E.
#' (2018). Targeted maximum likelihood estimation for a binary treatment: A
#' tutorial. *Statistics in Medicine*, 37(16), 2530–2546.
#'
#' @param x A fitted `"tmle_fit"` object.
#' @param ... Currently unused.
#'
#' @return A [tibble::tibble()] matching the input data (or selected subgroup).
#' @method augment tmle_fit
#' @export
augment.tmle_fit <- function(x, ...) {
  idx <- x$subgroup_idx
  t_res <- x$targeting
  df <- x$cache$data[idx, , drop = FALSE]

  dplyr::mutate(
    df,
    .q0_init = t_res$q0_init,
    .q1_init = t_res$q1_init,
    .q0_star = t_res$q0_star,
    .q1_star = t_res$q1_star,
    .g_hat = t_res$g_hat,
    .g_trunc = t_res$g_trunc,
    .clever_h = t_res$clever_h,
    .ite_star = t_res$q1_star - t_res$q0_star,
    .eif = t_res$eif,
    .eif_leverage = t_res$eif_leverage
  )
}

#' Collect Cross-Validated Performance Metrics for Nuisance Learners
#'
#' @param x A fitted `"tmle_fit"` object.
#' @param nuisance Which nuisance factor(s) to return: `"all"` (default), `"Q"`,
#'   `"g"`, or `"C"`.
#' @param summarize Logical; if `TRUE` (default), averages across folds and
#'   reports `mean`, `n`, and `std_err`. If `FALSE`, returns fold-level metrics.
#' @param ... Currently unused.
#'
#' @return A [tibble::tibble()] of `yardstick` CV metrics.
#' @method collect_metrics tmle_fit
#' @export
collect_metrics.tmle_fit <- function(x,
                                     nuisance = c("all", "Q", "g", "C"),
                                     summarize = TRUE,
                                     ...) {
  nuisance <- match.arg(nuisance)
  parts <- list(
    Q = x$nuisance$Q$metrics,
    g = x$nuisance$g$metrics,
    C = if (!is.null(x$nuisance$C)) x$nuisance$C$metrics else tibble::tibble()
  )

  raw <- if (nuisance == "all") {
    dplyr::bind_rows(parts)
  } else {
    parts[[nuisance]]
  }

  if (nrow(raw) == 0L || !isTRUE(summarize)) {
    return(raw)
  }

  raw |>
    dplyr::group_by(nuisance, learner, .metric, .estimator) |>
    dplyr::summarise(
      mean = mean(.estimate, na.rm = TRUE),
      n = sum(!is.na(.estimate)),
      std_err = stats::sd(.estimate, na.rm = TRUE) / sqrt(pmax(sum(!is.na(.estimate)), 1L)),
      .groups = "drop"
    )
}

#' Extract Super Learner Metalearner Weights Across Nuisance Models
#'
#' @param x A fitted `"tmle_fit"` object.
#' @param nuisance Which nuisance factor(s) to inspect: `"all"` (default), `"Q"`,
#'   `"g"`, or `"C"`.
#' @param by_fold Logical; if `FALSE` (default), returns the overall CV Super
#'   Learner ensemble weights. If `TRUE`, returns per-fold weights.
#'
#' @return A [tibble::tibble()] with columns `nuisance`, (optionally `fold`),
#'   `learner`, and `weight`.
#' @export
extract_sl_weights <- function(x,
                               nuisance = c("all", "Q", "g", "C"),
                               by_fold = FALSE) {
  if (!inherits(x, "tmle_fit")) {
    cli::cli_abort("{.arg x} must be a {.cls tmle_fit} object.")
  }
  nuisance <- match.arg(nuisance)
  factors <- if (nuisance == "all") c("Q", "g", "C") else nuisance

  purrr::map_dfr(factors, function(f) {
    sl_obj <- x$nuisance[[f]]
    if (is.null(sl_obj)) {
      return(tibble::tibble())
    }
    if (isTRUE(by_fold)) {
      return(sl_obj$fold_weights)
    }
    tibble::tibble(
      nuisance = f,
      learner = names(sl_obj$sl_weights),
      weight = as.numeric(sl_obj$sl_weights)
    )
  })
}

#' @export
print.tmle_spec <- function(x, ...) {
  lines <- cli::cli_format_method({
    cli::cli_h2("Targeted Maximum Likelihood Estimation (TMLE) Specification")
    cli::cli_bullets(c(
      "*" = "Outcome: {.field {rlang::as_label(x$roles$outcome)}}",
      "*" = "Treatment: {.field {rlang::as_label(x$roles$treatment)}}",
      "*" = "Covariates: {.field {if (rlang::quo_is_null(x$roles$covariates)) '<all remaining>' else rlang::as_label(x$roles$covariates)}}",
      "*" = "Estimand(s): {.val {x$estimand$names}}",
      "*" = "CV-TMLE: {.val {x$targeting$cv_tmle}} ({x$resamples$v_folds}-fold)",
      "*" = "Positivity bounds: {.val {paste(x$targeting$g_bounds, collapse = ', ')}} | Fluctuation: {.val {x$targeting$fluctuation}}"
    ))
  })
  cat(lines, sep = "\n")
  invisible(x)
}

#' @export
print.tmle_fit <- function(x, ...) {
  gl <- glance.tmle_fit(x)
  td <- x$targeting$estimates

  lines <- cli::cli_format_method({
    cli::cli_h2("TMLE Fit Summary")
    meta <- sprintf(
      "n = %d (observed = %d) | %d-fold %s | %s fluctuation | %.2fs",
      gl$n, gl$n_observed, gl$v_folds,
      if (gl$cv_tmle) "CV-TMLE" else "TMLE",
      x$spec$targeting$fluctuation,
      x$fit_time_sec
    )
    if (!is.na(gl$n_clusters)) {
      meta <- paste0(meta, sprintf(" | %d clusters", gl$n_clusters))
    }
    if (!is.null(x$subgroup_expr)) {
      meta <- paste0(meta, sprintf(" | subgroup: %s", x$subgroup_expr))
    }
    cli::cli_text("{.emph {meta}}")

    cli::cli_h3("Causal Estimates")
    for (i in seq_len(nrow(td))) {
      r <- td[i, ]
      p_str <- if (r$p.value < 1e-4) "< 0.0001" else sprintf("%.4f", r$p.value)
      cli::cli_text(
        "  {sprintf('%-6s', r$term)}  est = {sprintf('%8.4f', r$estimate)}  (SE = {sprintf('%.4f', r$std.error)})  95% CI [{sprintf('%.4f', r$conf.low)}, {sprintf('%.4f', r$conf.high)}]  p = {p_str}"
      )
    }

    cli::cli_h3("Nuisance Super Learners")
    format_sl_line <- function(sl_obj, label, metric_key, metric_label) {
      if (is.null(sl_obj)) return(NULL)
      m <- sl_obj$metrics |>
        dplyr::filter(learner == "SuperLearner", .metric == metric_key) |>
        dplyr::summarise(val = mean(.estimate, na.rm = TRUE)) |>
        dplyr::pull(val)
      w_sorted <- sort(sl_obj$sl_weights[sl_obj$sl_weights > 0.001], decreasing = TRUE)
      top_str <- paste(sprintf("%s %.2f", names(w_sorted), as.numeric(w_sorted)), collapse = ", ")
      cli::cli_text("  {label} Super Learner: {metric_label} = {sprintf('%.3f', m)} (weights: {top_str})")
    }

    if (x$cache$outcome_type == "binary") {
      format_sl_line(x$nuisance$Q, "Q", "roc_auc", "AUC")
    } else {
      format_sl_line(x$nuisance$Q, "Q", "rsq", "R2")
    }
    if (x$cache$treatment_type == "binary") {
      format_sl_line(x$nuisance$g, "g", "roc_auc", "AUC")
    } else {
      format_sl_line(x$nuisance$g, "g (ratio)", "mn_log_loss", "LogLoss")
    }
    if (!is.null(x$nuisance$C)) {
      format_sl_line(x$nuisance$C, "C", "roc_auc", "AUC")
    }

    cli::cli_h3("Diagnostics & Targeting Shift")
    cli::cli_text(
      "  Targeting shift (TMLE - plug-in): {sprintf('%+.4f', gl$plugin_vs_tmle_shift)} (epsilon = {sprintf('%.4f', gl$epsilon)})"
    )
    warned <- FALSE
    if (abs(gl$mean_eif) > 1e-4) {
      cli::cli_alert_danger(
        "EIF mean is non-zero ({sprintf('%.2e', gl$mean_eif)} > 1e-4); targeting equation may not have converged."
      )
      warned <- TRUE
    }
    if (gl$pct_truncated > 0.05) {
      cli::cli_alert_warning(
        "Positivity stress: {sprintf('%.1f%%', 100 * gl$pct_truncated)} of propensity scores hit g_bounds [{sprintf('%.3f', gl$g_bound_low)}, {sprintf('%.3f', gl$g_bound_high)}]."
      )
      warned <- TRUE
    }
    if (gl$max_eif_leverage > 0.10) {
      cli::cli_alert_warning(
        "High single-observation leverage: top observation explains {sprintf('%.1f%%', 100 * gl$max_eif_leverage)} of EIF variance."
      )
      warned <- TRUE
    }
    if (!warned) {
      cli::cli_alert_success(
        "All checks healthy: mean EIF = {sprintf('%.2e', gl$mean_eif)}, truncated = {sprintf('%.1f%%', 100 * gl$pct_truncated)}, max EIF leverage = {sprintf('%.1f%%', 100 * gl$max_eif_leverage)}."
      )
    }
  })

  cat(lines, sep = "\n")
  invisible(x)
}
