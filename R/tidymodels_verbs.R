#' Extract Tidy Causal Inference Table from a TMLE Fit
#'
#' Returns a publication-ready inference tibble containing point estimates,
#' influence-function standard errors, test statistics, p-values, confidence
#' intervals, optional simultaneous uniform confidence bands, and optional
#' VanderWeele-Ding E-values.
#'
#' @details
#' For additive contrasts (`"ate"`, `"att"`, `"atc"`, `"ey1"`, `"ey0"`, and
#' continuous Modified Treatment Policies), `std.error` is on the native outcome
#' scale and intervals are symmetric Wald bounds. For ratio estimands (`"rr"`
#' and `"or"`), `std.error` is reported on the **log scale** (`SE[log(RR)]` /
#' `SE[log(OR)]`) and confidence intervals are computed on the log scale and
#' exponentiated back so intervals remain strictly positive.
#'
#' @param x A fitted `"tmle_fit"` object.
#' @param conf.int Logical; whether to include `conf.low` and `conf.high`
#'   columns (default `TRUE`).
#' @param conf.level Confidence level for intervals; defaults to the level stored
#'   on `x` (`0.95` initially unless changed via [retarget()]).
#' @param simultaneous Logical; if `TRUE` and multiple estimands are present,
#'   computes simultaneous uniform confidence bands via a fast Rademacher wild
#'   multiplier bootstrap over the joint EIF matrix.
#' @param e_value Logical; if `TRUE`, appends `e_value` and `e_value_ci`
#'   sensitivity columns (VanderWeele & Ding, 2017) quantifying robustness to
#'   unmeasured confounding: exact on risk ratios (`"rr"`), square-root
#'   converted for common binary outcomes (`"or"` when outcome prevalence exceeds
#'   `15%`, exact otherwise), standardized-difference converted via
#'   `exp(0.91 * d)` on additive contrasts, and `NA` for non-contrast marginal
#'   means (`"ey1"`, `"ey0"`, uncontrasted continuous policies).
#' @param contrast Optional causal contrast override for continuous Modified
#'   Treatment Policies: `NULL` (default; uses the contrast declared inside
#'   [shift()] / [scale_tx()]), `"none"` (counterfactual policy mean), or
#'   `"observed"` (policy contrast minus observed outcome mean).
#' @param ... Currently unused.
#'
#' @return A [tibble::tibble()].
#' @method tidy tmle_fit
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, outcome_type = "binary", seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = c("ate", "rr", "or")) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' tidy(fit)
#' tidy(fit, simultaneous = TRUE, e_value = TRUE)
#'
#' @export
tidy.tmle_fit <- function(x,
                          conf.int = TRUE,
                          conf.level = NULL,
                          simultaneous = FALSE,
                          e_value = FALSE,
                          contrast = NULL,
                          ...) {
  stored_cl <- if (!is.null(x$targeting$conf_level)) x$targeting$conf_level else 0.95
  eff_cl <- if (is.null(conf.level)) stored_cl else conf.level
  if (!is.numeric(eff_cl) || length(eff_cl) != 1L || eff_cl <= 0 || eff_cl >= 1) {
    cli::cli_abort("{.arg conf.level} must be strictly between 0 and 1.")
  }

  est_spec <- x$spec$estimand
  contrast_override <- FALSE
  if (!is.null(contrast)) {
    contrast <- match.arg(contrast, c("none", "observed"))
    if (!identical(est_spec$kind, "policy")) {
      cli::cli_abort("{.arg contrast} in {.fn tidy} only applies to continuous Modified Treatment Policy fits.")
    }
    p_mod <- est_spec$policies[[1L]]
    p_mod$contrast <- contrast
    base_lbl <- if (!is.null(p_mod$base_label)) p_mod$base_label else p_mod$label
    p_mod$label <- if (identical(contrast, "observed")) paste0(base_lbl, " - observed") else base_lbl
    est_spec$policies[[1L]] <- p_mod
    est_spec$names <- p_mod$label
    contrast_override <- TRUE
  }

  # Recompute inference on cached targeting outputs if custom conf.level, simultaneous, e_value, or contrast requested
  if (abs(eff_cl - stored_cl) > 1e-12 || isTRUE(simultaneous) || isTRUE(e_value) || isTRUE(contrast_override)) {
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
      estimand = est_spec,
      outcome_type = cache$outcome_type,
      y_bounds = y_bounds_sub,
      g_bounds_arg = x$spec$targeting$g_bounds,
      fluctuation = x$spec$targeting$fluctuation,
      cluster_id = if (!is.null(cache$id_vec)) cache$id_vec[idx] else NULL,
      conf.level = eff_cl,
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
#' - **`epsilon`**: Scalar summary of the arm-specific targeting fluctuation
#'   parameters stored in `fit$targeting$epsilon_arms`. Values close to `0` mean
#'   the initial ML outcome fit already had little residual confounding bias;
#'   larger magnitudes show how strongly TMLE fluctuated the initial predictions
#'   to remove bias.
#' - **`plugin_vs_tmle_shift`**: Difference between the targeted estimate and
#'   the untargeted G-computation plug-in estimate for the primary estimand.
#'   Quantifies the exact causal bias calibration performed by the targeting step.
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
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' glance(fit)
#'
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
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' head(augment(fit)[, c(".q0_star", ".q1_star", ".g_trunc", ".clever_h", ".eif")])
#'
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
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' collect_metrics(fit)
#'
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
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' extract_sl_weights(fit)
#'
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
      "*" = "Outcome to measure: {.field {rlang::as_label(x$roles$outcome)}}",
      "*" = "Treatment / exposure: {.field {rlang::as_label(x$roles$treatment)}}",
      "*" = "Confounders adjusted for: {.field {if (rlang::quo_is_null(x$roles$covariates)) '<all remaining columns>' else rlang::as_label(x$roles$covariates)}}",
      "*" = "Causal question(s): {.val {x$estimand$names}}",
      "*" = "Cross-fitting: {.val {x$targeting$cv_tmle}} ({x$resamples$v_folds}-fold out-of-sample)",
      "*" = "Safety guardrails: bounds = {.val {paste(x$targeting$g_bounds, collapse = ', ')}} | calibration = {.val {x$targeting$fluctuation}}"
    ))
    cli::cli_text("{.emph Next: pass to fit(data = your_data) to estimate causal effects.}")
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
      "n = %d (%d observed) | %d-fold %s | %s calibration | %.2fs",
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

    ci_pct <- round(100 * (if (!is.null(x$targeting$conf_level)) x$targeting$conf_level else 0.95))
    cli::cli_h3("Causal Effect Estimates")
    for (i in seq_len(nrow(td))) {
      r <- td[i, ]
      p_str <- format_p_human(r$p.value)
      se_lab <- if (r$estimand %in% c("rr", "or")) "log-SE" else "SE"
      cli::cli_text(
        "  {sprintf('%-6s', r$term)}  est = {sprintf('%8.4f', r$estimate)}  ({se_lab} = {sprintf('%.4f', r$std.error)})  {ci_pct}% CI [{sprintf('%.4f', r$conf.low)}, {sprintf('%.4f', r$conf.high)}]  {p_str}"
      )
    }

    cli::cli_h3("Nuisance Models (Cross-Validated Super Learner)")
    format_sl_line <- function(sl_obj, label, metric_key, metric_label) {
      if (is.null(sl_obj)) return(NULL)
      m <- sl_obj$metrics |>
        dplyr::filter(learner == "SuperLearner", .metric == metric_key) |>
        dplyr::summarise(val = mean(.estimate, na.rm = TRUE)) |>
        dplyr::pull(val)
      w_sorted <- sort(sl_obj$sl_weights[sl_obj$sl_weights > 0.001], decreasing = TRUE)
      top_str <- paste(sprintf("%s %.0f%%", names(w_sorted), 100 * as.numeric(w_sorted)), collapse = ", ")
      cli::cli_text("  {label}: {metric_label} = {sprintf('%.3f', m)} (weights: {top_str})")
    }

    if (x$cache$outcome_type == "binary") {
      format_sl_line(x$nuisance$Q, "Outcome (Q)", "roc_auc", "AUC")
    } else {
      format_sl_line(x$nuisance$Q, "Outcome (Q)", "rsq", "R2")
    }
    if (x$cache$treatment_type == "binary") {
      format_sl_line(x$nuisance$g, "Treatment (g)", "roc_auc", "AUC")
    } else {
      format_sl_line(x$nuisance$g, "Policy weight (g)", "mn_log_loss", "LogLoss")
    }
    if (!is.null(x$nuisance$C)) {
      format_sl_line(x$nuisance$C, "Missingness (C)", "roc_auc", "AUC")
    }

    cli::cli_h3("Causal Health & Bias Adjustment")
    cli::cli_text(
      "  Bias calibration (TMLE - initial ML): {sprintf('%+.4f', gl$plugin_vs_tmle_shift)} (epsilon = {sprintf('%.4f', gl$epsilon)})"
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
        "Positivity stress: {sprintf('%.1f%%', 100 * gl$pct_truncated)} clipped at g_bounds [{sprintf('%.3f', gl$g_bound_low)}, {sprintf('%.3f', gl$g_bound_high)}]."
      )
      warned <- TRUE
    }
    if (gl$max_eif_leverage > 0.10) {
      cli::cli_alert_warning(
        "High single-row influence: top observation explains {sprintf('%.1f%%', 100 * gl$max_eif_leverage)} of uncertainty."
      )
      warned <- TRUE
    }
    if (!warned) {
      cli::cli_alert_success(
        "All core checks healthy: mean EIF = {sprintf('%.2e', gl$mean_eif)}, clipped = {sprintf('%.1f%%', 100 * gl$pct_truncated)}, top leverage = {sprintf('%.1f%%', 100 * gl$max_eif_leverage)}."
      )
    }
    cli::cli_text("{.emph Next: summary(fit) | check_tmle(fit) | report_tmle(fit) | autoplot(fit)}")
  })

  cat(lines, sep = "\n")
  invisible(x)
}

#' Plain-English Causal Interpretation & Audit for a TMLE Fit
#'
#' Translates a `"tmle_fit"` object into plain-English causal sentences,
#' explains how strongly the TMLE targeting step corrected residual bias from
#' the initial machine-learning fit, and runs a traffic-light diagnostic audit
#' with actionable next steps (`retarget()`, `autoplot()`, `augment()`).
#' `summary(fit)` delegates directly to `explain_tmle()`.
#'
#' @param x,object A fitted `"tmle_fit"` object.
#' @param conf.level Confidence level for interval statements; defaults to the
#'   level stored on `x` (`0.95` initially).
#' @param ... Currently unused.
#'
#' @return Invisibly returns a list with elements `estimates` (tibble) and
#'   `health` (1-row `glance()` tibble).
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, outcome_type = "binary", seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = c("ate", "rr")) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' summary(fit)
#'
#' @export
explain_tmle <- function(x, conf.level = NULL) {
  if (!inherits(x, "tmle_fit")) {
    cli::cli_abort("{.arg x} must be a {.cls tmle_fit} object.")
  }
  eff_cl <- if (is.null(conf.level)) {
    if (!is.null(x$targeting$conf_level)) x$targeting$conf_level else 0.95
  } else {
    conf.level
  }
  gl <- glance.tmle_fit(x)
  td <- tidy.tmle_fit(x, conf.level = eff_cl)
  pct_lbl <- sprintf("%.0f%%", 100 * eff_cl)

  describe_row <- function(r) {
    p_txt <- format_p_human(r$p.value)
    if (r$estimand %in% c("ate", "att", "atc")) {
      pop <- switch(
        r$estimand,
        ate = "across the full study population",
        att = "among treated individuals (ATT)",
        atc = "among untreated/comparison individuals (ATC)"
      )
      dir_txt <- if (r$estimate >= 0) "increases" else "decreases"
      if (identical(x$cache$outcome_type, "binary")) {
        sprintf(
          "%s (%s): Treatment %s the expected outcome probability by %.4f (%.1f percentage points; %s CI [%.4f, %.4f], %s).",
          r$term, pop, dir_txt, abs(r$estimate), 100 * abs(r$estimate), pct_lbl, r$conf.low, r$conf.high, p_txt
        )
      } else {
        sprintf(
          "%s (%s): Treatment %s the expected outcome by %.4f units (%s CI [%.4f, %.4f], %s).",
          r$term, pop, dir_txt, abs(r$estimate), pct_lbl, r$conf.low, r$conf.high, p_txt
        )
      }
    } else if (r$estimand == "rr") {
      pct_ch <- 100 * (r$estimate - 1)
      dir_txt <- if (r$estimate >= 1) "higher" else "lower"
      sprintf(
        "RR: Treatment corresponds to a %.2fx risk ratio (%.1f%% %s risk; %s CI [%.2f, %.2f], %s).",
        r$estimate, abs(pct_ch), dir_txt, pct_lbl, r$conf.low, r$conf.high, p_txt
      )
    } else if (r$estimand == "or") {
      sprintf(
        "OR: Treatment corresponds to an odds ratio of %.2f (%s CI [%.2f, %.2f], %s).",
        r$estimate, pct_lbl, r$conf.low, r$conf.high, p_txt
      )
    } else if (r$estimand %in% c("ey1", "ey0")) {
      arm_txt <- if (r$estimand == "ey1") "everyone received treatment" else "no one received treatment"
      sprintf(
        "%s: Counterfactual average if %s is %.4f (%s CI [%.4f, %.4f]).",
        r$term, arm_txt, r$estimate, pct_lbl, r$conf.low, r$conf.high
      )
    } else if (grepl(" - observed$", r$estimand)) {
      sprintf(
        "%s: Modified Treatment Policy shifted expected outcome relative to observed exposure by %.4f units (%s CI [%.4f, %.4f], %s).",
        r$term, r$estimate, pct_lbl, r$conf.low, r$conf.high, p_txt
      )
    } else {
      sprintf(
        "%s: Expected outcome under Modified Treatment Policy is %.4f (%s CI [%.4f, %.4f], %s).",
        r$term, r$estimate, pct_lbl, r$conf.low, r$conf.high, p_txt
      )
    }
  }

  bal_info <- if (identical(x$cache$treatment_type, "binary")) compute_covariate_balance(x) else NULL

  lines <- cli::cli_format_method({
    cli::cli_h2("Plain-English TMLE Summary")
    cli::cli_text(
      "Cohort: {gl$n} rows ({gl$n_observed} with observed outcomes) | {gl$v_folds}-fold {if (gl$cv_tmle) 'cross-fitted TMLE' else 'TMLE'}"
    )

    cli::cli_h3("1. What the Causal Estimates Mean")
    for (i in seq_len(nrow(td))) {
      cli::cli_alert_info("{describe_row(td[i, ])}")
    }

    cli::cli_h3("2. What the Targeting Step Did (\"Predict, Then Target\")")
    if (identical(x$spec$targeting$fluctuation, "none")) {
      cli::cli_alert_warning(
        "Targeting fluctuation is disabled (fluctuation = 'none'): this fit reports an untargeted ML plug-in estimate (mean EIF = {sprintf('%.2e', gl$mean_eif)})."
      )
    } else {
      shift_abs <- abs(gl$plugin_vs_tmle_shift)
      if (shift_abs < 1e-3) {
        cli::cli_text(
          "The initial ML predictions already captured nearly all measured confounding bias (calibration shift = {sprintf('%+.4f', gl$plugin_vs_tmle_shift)}, epsilon = {sprintf('%.4f', gl$epsilon)})."
        )
      } else {
        cli::cli_text(
          "Targeting nudged the initial ML plug-in predictions by {sprintf('%+.4f', gl$plugin_vs_tmle_shift)} (epsilon = {sprintf('%.4f', gl$epsilon)}) to zero out residual confounding bias."
        )
      }
    }

    cli::cli_h3("3. Can I Trust These Numbers? (Health Check & Next Steps)")
    if (abs(gl$mean_eif) <= 1e-4) {
      cli::cli_alert_success("Convergence: Bias-removal score equation solved (mean EIF = {sprintf('%.2e', gl$mean_eif)}).")
    } else {
      cli::cli_alert_danger("Convergence: mean EIF is {sprintf('%.2e', gl$mean_eif)} (> 1e-4); run check_tmle(fit) or inspect extreme weights.")
    }

    if (gl$pct_truncated <= 0.05) {
      cli::cli_alert_success(
        "Comparison overlap: {sprintf('%.1f%%', 100 * gl$pct_truncated)} clipped at [{sprintf('%.3f', gl$g_bound_low)}, {sprintf('%.3f', gl$g_bound_high)}] (effective n: {round(gl$ess_treated, 1)} treated, {round(gl$ess_control, 1)} control)."
      )
    } else {
      cli::cli_alert_warning(
        "Positivity stress: {sprintf('%.1f%%', 100 * gl$pct_truncated)} clipped at [{sprintf('%.3f', gl$g_bound_low)}, {sprintf('%.3f', gl$g_bound_high)}]. Test stability with {.code autoplot(fit, type = 'sensitivity')} or {.code retarget(fit, g_bounds = c(0.05, 0.95))}."
      )
    }

    if (!is.null(bal_info) && nrow(bal_info) > 0L) {
      max_raw <- max(bal_info$smd_raw, na.rm = TRUE)
      max_wt <- max(bal_info$smd_wt, na.rm = TRUE)
      if (max_wt <= 0.10) {
        bal_msg <- if (max_wt < max_raw) {
          "Confounder balance: Max standardized difference fell from {sprintf('%.3f', max_raw)} to {sprintf('%.3f', max_wt)} (well below 0.10)."
        } else {
          "Confounder balance: Max standardized difference remained below 0.10 ({sprintf('%.3f', max_raw)} raw vs. {sprintf('%.3f', max_wt)} weighted)."
        }
        cli::cli_alert_success(bal_msg)
      } else {
        cli::cli_alert_warning(
          "Confounder balance: Max weighted |SMD| is {sprintf('%.3f', max_wt)} (> 0.10). Inspect {.code autoplot(fit, type = 'balance')}."
        )
      }
    }

    if (gl$max_eif_leverage <= 0.10) {
      cli::cli_alert_success(
        "Uncertainty stability: No single person dominates SEs (top observation share = {sprintf('%.1f%%', 100 * gl$max_eif_leverage)})."
      )
    } else {
      cli::cli_alert_warning(
        "High leverage: Top observation accounts for {sprintf('%.1f%%', 100 * gl$max_eif_leverage)} of variance. Inspect {.code augment(fit) |> dplyr::slice_max(.eif_leverage, n = 5)}."
      )
    }

    cli::cli_text("{.emph Next steps: report_tmle(fit) for paper-ready prose | autoplot(fit, type = 'dashboard') for visual audit.}")
  })

  cat(lines, sep = "\n")
  invisible(list(estimates = td, health = gl))
}

#' @rdname explain_tmle
#' @method summary tmle_fit
#' @export
summary.tmle_fit <- function(object, conf.level = NULL, ...) {
  explain_tmle(object, conf.level = conf.level)
}
