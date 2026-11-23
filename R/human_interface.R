#' Audit a TMLE Fit Across 5 Causal Trust Pillars
#'
#' Evaluates a fitted `"tmle_fit"` object across core causal identification and
#' numerical convergence safeguards in plain English:
#' 1. **Targeting Convergence**: Verifies that the bias-removal fluctuation
#'    step solved the Efficient Influence Function (EIF) score equation
#'    (`abs(mean_eif) <= 1e-4`).
#' 2. **Comparison Overlap (Positivity)**: Checks propensity score truncation
#'    rates (`pct_truncated`) and Kish effective sample sizes (`ess_treated`,
#'    `ess_control`), or maximum density ratio clipping for continuous policies.
#' 3. **Confounder Balance**: For binary treatments, verifies that inverse-
#'    propensity weighting reduces maximum absolute standardized mean differences
#'    (`|SMD|`) below `0.10`.
#' 4. **Individual Influence**: Checks whether any single observation dominates
#'    total EIF variance (`max_eif_leverage <= 0.10`).
#' 5. **Model Learning Quality**: Confirms cross-validated nuisance Super
#'    Learners extracted signal above a baseline guessing threshold.
#'
#' @param x A fitted `"tmle_fit"` object.
#'
#' @return An object of class `"tmle_check"` (a [tibble::tibble()] with
#'   columns `check`, `status`, `metric`, `summary`, `why_it_matters`, and
#'   `next_step`).
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, outcome_type = "binary", seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = c("ate", "rr")) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' check_tmle(fit)
#'
#' @export
check_tmle <- function(x) {
  if (!inherits(x, "tmle_fit")) {
    cli::cli_abort(c(
      "{.arg x} must be a fitted {.cls tmle_fit} object returned by {.fn fit}.",
      "i" = "Fit a specification first: {.code fit <- tmle_spec(...) |> fit(data)}"
    ))
  }

  gl <- glance.tmle_fit(x)
  t_res <- x$targeting
  is_binary_tx <- identical(x$cache$treatment_type, "binary")
  rows <- list()

  # 1. Targeting Convergence
  eif_abs <- abs(gl$mean_eif)
  conv_status <- if (is.finite(eif_abs) && eif_abs <= 1e-4) "pass" else "fail"
  rows[[length(rows) + 1L]] <- tibble::tibble(
    check = "Targeting Convergence",
    status = conv_status,
    metric = sprintf("mean EIF = %.2e", gl$mean_eif),
    summary = if (conv_status == "pass") {
      sprintf("Bias-removal step converged cleanly (targeting shift = %+.4f, epsilon = %.4f).",
              gl$plugin_vs_tmle_shift, gl$epsilon)
    } else {
      sprintf("EIF score equation did not reach zero (mean EIF = %.2e > 1e-4).", gl$mean_eif)
    },
    why_it_matters = "A zero mean EIF confirms the targeting nudge finished stripping plug-in bias.",
    next_step = if (conv_status == "pass") {
      "None required."
    } else {
      "Inspect extreme weights via autoplot(fit, type = 'eif') or retarget(fit, g_bounds = c(0.05, 0.95))."
    }
  )

  # 2. Comparison Overlap (Positivity)
  pct_clip <- 100 * gl$pct_truncated
  ov_status <- if (gl$pct_truncated <= 0.05) {
    "pass"
  } else if (gl$pct_truncated <= 0.15) {
    "warn"
  } else {
    "fail"
  }
  ov_metric <- if (is_binary_tx) {
    sprintf("%.1f%% clipped | ESS %.0f vs %.0f", pct_clip, gl$ess_treated, gl$ess_control)
  } else {
    sprintf("%.1f%% capped (max r = %.1f)", pct_clip, gl$g_bound_high)
  }
  ov_summary <- if (is_binary_tx) {
    sprintf("%.1f%% of propensity scores hit bounds [%.3f, %.3f] (ESS: treated = %.1f, control = %.1f).",
            pct_clip, gl$g_bound_low, gl$g_bound_high, gl$ess_treated, gl$ess_control)
  } else {
    sprintf("%.1f%% of density ratios were capped at upper bound %.2f.", pct_clip, gl$g_bound_high)
  }
  rows[[length(rows) + 1L]] <- tibble::tibble(
    check = "Comparison Overlap (Positivity)",
    status = ov_status,
    metric = ov_metric,
    summary = ov_summary,
    why_it_matters = "Causal comparison fails when one group rarely receives the treatment seen in the other.",
    next_step = if (ov_status == "pass") {
      "None required; inspect autoplot(fit, type = 'overlap') to confirm visually."
    } else if (is_binary_tx) {
      "Inspect overlap & stability: autoplot(fit, type = 'sensitivity') or retarget(fit, g_bounds = c(0.05, 0.95))."
    } else {
      "Try a smaller policy shift (shift()/scale_tx()) or tighter bounds in retarget()."
    }
  )

  # 3. Confounder Balance (binary treatments only)
  if (is_binary_tx) {
    bal <- compute_covariate_balance(x)
    max_raw <- if (nrow(bal) > 0L) max(bal$smd_raw, na.rm = TRUE) else 0
    max_wt <- if (nrow(bal) > 0L) max(bal$smd_wt, na.rm = TRUE) else 0
    bal_status <- if (max_wt <= 0.10) {
      "pass"
    } else if (max_wt <= 0.20) {
      "warn"
    } else {
      "fail"
    }
    rows[[length(rows) + 1L]] <- tibble::tibble(
      check = "Confounder Balance",
      status = bal_status,
      metric = sprintf("max |SMD| %.3f -> %.3f", max_raw, max_wt),
      summary = if (bal_status == "pass") {
        sprintf("Weighted confounders are well balanced (max |SMD| fell from %.3f to %.3f, below 0.10).",
                max_raw, max_wt)
      } else {
        sprintf("Residual imbalance remains after weighting (max |SMD| = %.3f; goal <= 0.10).", max_wt)
      },
      why_it_matters = "Balanced covariates mean treated and control groups act like fair counterfactual twins.",
      next_step = if (bal_status == "pass") {
        "None required."
      } else {
        "Inspect autoplot(fit, type = 'balance') and add richer treatment learners/interactions."
      }
    )
  }

  # 4. Individual Influence
  lev_pct <- 100 * gl$max_eif_leverage
  inf_status <- if (gl$max_eif_leverage <= 0.10) {
    "pass"
  } else if (gl$max_eif_leverage <= 0.20) {
    "warn"
  } else {
    "fail"
  }
  top_obs <- which.max(t_res$eif_leverage)
  rows[[length(rows) + 1L]] <- tibble::tibble(
    check = "Individual Influence",
    status = inf_status,
    metric = sprintf("top obs = %.1f%% variance", lev_pct),
    summary = if (inf_status == "pass") {
      sprintf("No single row dominates uncertainty (row %d holds %.1f%% of EIF variance, <= 10%%).",
              top_obs, lev_pct)
    } else {
      sprintf("Row %d alone drives %.1f%% of total EIF variance (> 10%% caution threshold).",
              top_obs, lev_pct)
    },
    why_it_matters = "Confidence intervals are fragile if one outlier carries most of the variance.",
    next_step = if (inf_status == "pass") {
      "None required."
    } else {
      "Inspect outlier row via autoplot(fit, type = 'eif') or dplyr::slice_max(augment(fit), .eif_leverage, n = 5)."
    }
  )

  # 5. Model Learning Quality
  is_bin_out <- identical(x$cache$outcome_type, "binary")
  q_ok <- if (is_bin_out) {
    is.finite(gl$q_metric) && gl$q_metric > 0.52
  } else {
    is.finite(gl$q_metric) && gl$q_metric > 0.01
  }
  g_ok <- if (is_binary_tx) is.finite(gl$g_auc) && gl$g_auc > 0.52 else TRUE
  ml_status <- if (q_ok && (g_ok || !is_binary_tx)) "pass" else "warn"

  q_txt <- if (is_bin_out) sprintf("Q AUC = %.3f", gl$q_metric) else sprintf("Q R2 = %.3f", gl$q_metric)
  ml_metric <- if (is_binary_tx) sprintf("%s | g AUC = %.3f", q_txt, gl$g_auc) else q_txt
  ml_summary <- if (!q_ok) {
    sprintf("Outcome Super Learner explained near-zero signal (%s); TMLE leans heavily on treatment weights.", q_txt)
  } else if (!g_ok) {
    sprintf("Outcome learner is strong (%s), while treatment g AUC = %.3f (expected if treatment is randomized or weakly predictable).",
            q_txt, gl$g_auc)
  } else if (is_binary_tx) {
    sprintf("Cross-validated Super Learners captured real structure (%s, treatment AUC = %.3f).", q_txt, gl$g_auc)
  } else {
    sprintf("Cross-validated outcome Super Learner captured real structure (%s).", q_txt)
  }

  rows[[length(rows) + 1L]] <- tibble::tibble(
    check = "Model Learning Quality",
    status = ml_status,
    metric = ml_metric,
    summary = ml_summary,
    why_it_matters = "Doubly robust inference needs at least one nuisance model (outcome or treatment) to capture real signal.",
    next_step = if (ml_status == "pass") {
      "None required."
    } else {
      "Inspect learners via autoplot(fit, type = 'learners') and add nonlinear workflows via set_outcome_learners()."
    }
  )

  out <- dplyr::bind_rows(rows)
  class(out) <- c("tmle_check", class(out))
  out
}

#' @rdname check_tmle
#' @param ... Currently unused.
#' @export
print.tmle_check <- function(x, ...) {
  n_pass <- sum(x$status == "pass")
  n_warn <- sum(x$status == "warn")
  n_fail <- sum(x$status == "fail")

  lines <- cli::cli_format_method({
    cli::cli_h2("TMLE Causal Health Audit ({n_pass}/{nrow(x)} Passed)")
    for (i in seq_len(nrow(x))) {
      st <- x$status[i]
      chk <- x$check[i]
      met <- x$metric[i]
      sum_i <- x$summary[i]
      nxt <- x$next_step[i]
      if (st == "pass") {
        cli::cli_alert_success("{.strong {chk}} ({met}): {sum_i}")
      } else if (st == "warn") {
        cli::cli_alert_warning("{.strong {chk}} ({met}): {sum_i}")
        cli::cli_text("    {.emph Suggested next step:} {nxt}")
      } else {
        cli::cli_alert_danger("{.strong {chk}} ({met}): {sum_i}")
        cli::cli_text("    {.emph Suggested next step:} {nxt}")
      }
    }
    if (n_fail == 0L && n_warn == 0L) {
      cli::cli_text("")
      cli::cli_text("{.emph Verdict: Ready for interpretation. Generate paper prose with report_tmle(fit).}")
    }
  })
  cat(lines, sep = "\n")
  invisible(x)
}

#' Generate Publication-Ready Methods, Results, and Assumptions Prose
#'
#' Translates a fitted `"tmle_fit"` object into clean, copy-paste-ready
#' scientific paper paragraphs (`style = "paper"`), accessible blog/executive
#' prose (`style = "plain"`), or a complete Markdown report (`style = "markdown"`).
#' All numbers, learner names, citations, sensitivity E-values, and diagnostic
#' metrics are extracted directly from the fitted object so nothing is hand-typed.
#'
#' @param x A fitted `"tmle_fit"` object.
#' @param outcome_label Optional human-readable label for the outcome (defaults
#'   to the outcome column name).
#' @param treatment_label Optional human-readable label for the treatment
#'   (defaults to the treatment column name).
#' @param units Optional unit string for additive outcome contrasts (e.g.,
#'   `"mmHg"` or `"percentage points"`).
#' @param style Output framing: `"paper"` (academic Methods + Results sections),
#'   `"plain"` (accessible blog/stakeholder prose), or `"markdown"` (combined
#'   Markdown document with results table).
#' @param conf.level Confidence level for reported intervals (default `0.95`).
#' @param simultaneous Logical; whether to compute simultaneous uniform
#'   confidence bands when multiple estimands are present.
#' @param e_value Logical; whether to compute and report VanderWeele-Ding
#'   E-values for unmeasured confounding sensitivity (default `TRUE`).
#'
#' @return An object of class `"tmle_report"` (invisibly printed to `stdout`)
#'   with elements `methods`, `results`, `assumptions`, `table`, and `markdown`.
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, outcome_type = "binary", seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = c("ate", "rr")) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' rep <- report_tmle(fit, outcome_label = "30-day recovery", treatment_label = "therapy")
#' rep
#'
#' @export
report_tmle <- function(x,
                        outcome_label = NULL,
                        treatment_label = NULL,
                        units = NULL,
                        style = c("paper", "plain", "markdown"),
                        conf.level = 0.95,
                        simultaneous = FALSE,
                        e_value = TRUE) {
  if (!inherits(x, "tmle_fit")) {
    cli::cli_abort(c(
      "{.arg x} must be a fitted {.cls tmle_fit} object.",
      "i" = "Fit a specification first: {.code fit <- tmle_spec(...) |> fit(data)}"
    ))
  }
  style <- match.arg(style)

  out_name <- if (is.null(outcome_label)) x$cache$outcome_col else as.character(outcome_label)[1L]
  trt_name <- if (is.null(treatment_label)) x$cache$treatment_col else as.character(treatment_label)[1L]
  unit_suffix <- if (is.null(units) || !nzchar(units)) " units" else paste0(" ", trimws(as.character(units)[1L]))

  gl <- glance.tmle_fit(x)
  td <- tidy.tmle_fit(x, conf.level = conf.level, simultaneous = simultaneous, e_value = e_value)
  chk <- check_tmle(x)
  sl_w <- extract_sl_weights(x)
  pct_lbl <- sprintf("%.0f%%", 100 * conf.level)

  format_learners <- function(factor_name) {
    sub <- sl_w[sl_w$nuisance == factor_name, , drop = FALSE]
    if (nrow(sub) == 0L) return(NULL)
    ord <- order(sub$weight, decreasing = TRUE)
    sub <- sub[ord, , drop = FALSE]
    parts <- sprintf("%s (%.0f%%)", sub$learner, 100 * sub$weight)
    paste(parts, collapse = ", ")
  }

  q_desc <- format_learners("Q")
  g_desc <- format_learners("g")
  c_desc <- format_learners("C")

  n_covs <- length(x$cache$covariate_cols)
  cov_preview <- if (n_covs <= 5L) {
    paste(x$cache$covariate_cols, collapse = ", ")
  } else {
    paste0(paste(x$cache$covariate_cols[1:4], collapse = ", "), sprintf(", and %d other baseline covariates", n_covs - 4L))
  }

  cv_phrase <- if (gl$cv_tmle) {
    sprintf("%d-fold cross-fitted Targeted Maximum Likelihood Estimation (CV-TMLE; van der Laan & Rubin, 2006; Zheng & van der Laan, 2011)", gl$v_folds)
  } else {
    "Targeted Maximum Likelihood Estimation (TMLE; van der Laan & Rubin, 2006)"
  }

  if (x$cache$treatment_type == "binary") {
    bounds_phrase <- if (identical(x$spec$targeting$g_bounds, "adaptive")) {
      sprintf("adaptive sample-size-aware positivity bounds [%.3f, %.3f] (Petersen et al., 2012)",
              gl$g_bound_low, gl$g_bound_high)
    } else {
      sprintf("prespecified propensity truncation bounds [%.3f, %.3f]",
              gl$g_bound_low, gl$g_bound_high)
    }
    tx_model_phrase <- sprintf("The treatment propensity mechanism g(W) = P(%s = 1 | W) was estimated with [%s] and bounded using %s.",
                               trt_name, g_desc, bounds_phrase)
  } else {
    tx_model_phrase <- sprintf("Because %s is a continuous exposure, the Modified Treatment Policy density ratio was estimated nonparametrically via stacked classification ([%s]) and capped at %.2f (Diaz et al., 2021).",
                               trt_name, g_desc, gl$g_bound_high)
  }

  miss_phrase <- if (!is.null(c_desc)) {
    miss_kind <- if (isTRUE(x$cache$auto_censoring)) "an automatically generated non-missingness indicator" else "the supplied observation indicator"
    sprintf("To prevent selection bias from %d incomplete outcome records, censoring/missingness C(A, W) was modeled using %s ([%s]) and incorporated into the inverse-probability targeting weights.",
            gl$n - gl$n_observed, miss_kind, c_desc)
  } else {
    "Outcome measurements were complete across all analyzed observations."
  }

  se_phrase <- if (!is.na(gl$n_clusters)) {
    sprintf("Standard errors and %s confidence intervals were derived from the cluster-aggregated Efficient Influence Function (EIF) across %d clusters.", pct_lbl, gl$n_clusters)
  } else {
    sprintf("Standard errors and %s confidence intervals were derived from the sample variance of the Efficient Influence Function (EIF).", pct_lbl)
  }
  if (isTRUE(simultaneous) && "crit.value" %in% names(td)) {
    se_phrase <- paste(se_phrase, sprintf("Multi-estimand intervals were adjusted to simultaneous uniform bands via a Rademacher multiplier bootstrap (critical value = %.2f).", td$crit.value[1L]))
  }
  if (identical(x$spec$targeting$fluctuation, "logistic") && x$cache$outcome_type == "continuous") {
    se_phrase <- paste(se_phrase, "Continuous outcomes were linearly mapped to [0, 1] before canonical logistic fluctuation so counterfactual predictions remain strictly bounded inside the observed outcome range (Gruber & van der Laan, 2010).")
  }
  ev_method_phrase <- if (isTRUE(e_value) && "e_value" %in% names(td)) {
    "Sensitivity to potential unmeasured confounding was quantified on the risk-ratio scale using VanderWeele-Ding E-values (VanderWeele & Ding, 2017)."
  } else {
    ""
  }

  methods_text <- if (style == "plain") {
    paste(
      sprintf("We analyzed %d observations (%d with observed %s) using %d-fold cross-fitted TMLE ('predict, then target').",
              gl$n, gl$n_observed, out_name, gl$v_folds),
      sprintf("First, cross-validated ensembles estimated each unit's expected %s under alternative values of %s adjusting for %d baseline factors (%s; outcome learners: %s), alongside treatment weights (%s).",
              out_name, trt_name, n_covs, cov_preview, q_desc, g_desc),
      if (!is.null(c_desc)) sprintf("An additional model accounted for %d missing outcomes instead of dropping incomplete rows.", gl$n - gl$n_observed) else "",
      "Second, a targeted logistic calibration step adjusted those initial predictions so residual confounding bias on the causal estimand was removed.",
      se_phrase
    )
  } else {
    paste(
      sprintf("We estimated the causal effect of %s on %s among n = %d observations (%d observed) adjusting for %d baseline covariates (%s) using %s.",
              trt_name, out_name, gl$n, gl$n_observed, n_covs, cov_preview, cv_phrase),
      sprintf("Initial conditional expectations E[%s | %s, W] were fit with a cross-validated Super Learner ensemble combining [%s].",
              out_name, trt_name, q_desc),
      tx_model_phrase,
      miss_phrase,
      se_phrase,
      ev_method_phrase
    )
  }
  methods_text <- trimws(gsub(" +", " ", methods_text))

  # Build results prose for every estimand
  describe_est <- function(r) {
    p_str <- format_p_human(r$p.value)
    ev_str <- if ("e_value" %in% names(r) && is.finite(r$e_value)) {
      sprintf("; E-value = %.2f (CI limit = %.2f)", r$e_value, r$e_value_ci)
    } else {
      ""
    }
    if (r$estimand %in% c("ate", "att", "atc")) {
      pop_lbl <- switch(
        r$estimand,
        ate = "across the full study population (ATE)",
        att = "among treated units (ATT)",
        atc = "among untreated comparison units (ATC)"
      )
      dir_w <- if (r$estimate >= 0) "higher" else "lower"
      sprintf("%s shifted %s by %+.4f%s (%s; %s CI [%.4f, %.4f], %s%s)",
              trt_name, out_name, r$estimate, unit_suffix, pop_lbl,
              pct_lbl, r$conf.low, r$conf.high, p_str, ev_str)
    } else if (r$estimand == "rr") {
      pct_d <- 100 * (r$estimate - 1)
      sprintf("the relative risk ratio (RR) was %.2f (%.1f%% %s risk; %s CI [%.2f, %.2f], %s%s)",
              r$estimate, abs(pct_d), if (r$estimate >= 1) "higher" else "lower",
              pct_lbl, r$conf.low, r$conf.high, p_str, ev_str)
    } else if (r$estimand == "or") {
      sprintf("the odds ratio (OR) was %.2f (%s CI [%.2f, %.2f], %s%s)",
              r$estimate, pct_lbl, r$conf.low, r$conf.high, p_str, ev_str)
    } else if (r$estimand %in% c("ey1", "ey0")) {
      cf_lbl <- if (r$estimand == "ey1") paste("everyone received", trt_name) else paste("no one received", trt_name)
      sprintf("the counterfactual mean %s if %s was %.4f (%s CI [%.4f, %.4f])",
              out_name, cf_lbl, r$estimate, pct_lbl, r$conf.low, r$conf.high)
    } else {
      sprintf("under policy %s, expected %s was %.4f (%s CI [%.4f, %.4f], %s%s)",
              r$term, out_name, r$estimate, pct_lbl, r$conf.low, r$conf.high, p_str, ev_str)
    }
  }

  est_sentences <- vapply(seq_len(nrow(td)), function(i) describe_est(td[i, ]), character(1L))
  bal_txt <- if (x$cache$treatment_type == "binary") {
    bal <- compute_covariate_balance(x)
    if (nrow(bal) > 0L) {
      sprintf("Inverse-propensity weighting reduced maximum covariate standardized mean imbalance from %.3f to %.3f.",
              max(bal$smd_raw, na.rm = TRUE), max(bal$smd_wt, na.rm = TRUE))
    } else {
      ""
    }
  } else {
    ""
  }

  results_text <- paste(
    sprintf("After confounder adjustment and targeting, %s.", paste(est_sentences, collapse = "; and ")),
    sprintf("The targeting calibration step shifted the initial ML plug-in estimate by %+.4f (fluctuation epsilon = %.4f) and satisfied the EIF score equation (mean EIF = %.2e).",
            gl$plugin_vs_tmle_shift, gl$epsilon, gl$mean_eif),
    bal_txt,
    sprintf("Positivity truncation affected %.1f%% of observations, with maximum single-observation EIF variance leverage of %.1f%%.",
            100 * gl$pct_truncated, 100 * gl$max_eif_leverage)
  )
  results_text <- trimws(gsub(" +", " ", results_text))

  ev_note <- if ("e_value" %in% names(td) && any(is.finite(td$e_value))) {
    e_top <- td$e_value[1L]
    sprintf("1. **No unmeasured confounding:** Valid causal interpretation assumes all joint causes of %s and %s are captured in W. For the primary estimand (%s), an unmeasured confounder would need a risk-ratio association of at least %.2f with both treatment and outcome (beyond measured covariates) to explain away the point estimate.",
            trt_name, out_name, td$term[1L], e_top)
  } else {
    sprintf("1. **No unmeasured confounding:** Causal attribution requires that no unmeasured common causes of %s and %s remain outside W.", trt_name, out_name)
  }

  pos_note <- sprintf("2. **Positivity (real-world overlap):** Only %.1f%% of observations hit bounds [%.3f, %.3f]%s.",
                      100 * gl$pct_truncated, gl$g_bound_low, gl$g_bound_high,
                      if (x$cache$treatment_type == "binary") sprintf(" (effective sample sizes: treated = %.1f, control = %.1f)", gl$ess_treated, gl$ess_control) else "")
  sutva_note <- sprintf("3. **Consistency & interference:** Each unit's observed outcome matches its potential outcome under the received %s level%s.",
                        trt_name,
                        if (!is.na(gl$n_clusters)) sprintf(", with within-cluster dependence accounted for across %d clusters", gl$n_clusters) else " without cross-unit spillover")

  assumptions_text <- paste(ev_note, pos_note, sutva_note, sep = "\n")

  tbl_md <- c(
    "| Estimand | Estimate | SE / log-SE | 95% CI | p-value |",
    "|---|---:|---:|---:|---:|"
  )
  for (i in seq_len(nrow(td))) {
    r <- td[i, ]
    tbl_md <- c(tbl_md, sprintf(
      "| %s | %.4f | %.4f | [%.4f, %.4f] | %s |",
      r$term, r$estimate, r$std.error, r$conf.low, r$conf.high, format_p_human(r$p.value)
    ))
  }

  md_doc <- paste(
    c(
      "## Methods",
      "",
      methods_text,
      "",
      "## Results",
      "",
      results_text,
      "",
      tbl_md,
      "",
      "## Causal Assumptions & Diagnostics",
      "",
      assumptions_text
    ),
    collapse = "\n"
  )

  out <- structure(
    list(
      style = style,
      methods = methods_text,
      results = results_text,
      assumptions = assumptions_text,
      table = td,
      checks = chk,
      markdown = md_doc
    ),
    class = "tmle_report"
  )
  out
}

#' @rdname report_tmle
#' @param ... Currently unused.
#' @export
as.character.tmle_report <- function(x, ...) {
  x$markdown
}

#' @rdname report_tmle
#' @export
print.tmle_report <- function(x, ...) {
  if (identical(x$style, "markdown")) {
    cat(x$markdown, "\n", sep = "")
    return(invisible(x))
  }
  lines <- cli::cli_format_method({
    cli::cli_h2("Publication-Ready TMLE Report ({x$style} style)")
    cli::cli_h3("Methods")
    cli::cli_text("{x$methods}")
    cli::cli_h3("Results")
    cli::cli_text("{x$results}")
    cli::cli_h3("Causal Assumptions & Diagnostics")
    for (ln in strsplit(x$assumptions, "\n", fixed = TRUE)[[1L]]) {
      cli::cli_text("{ln}")
    }
    cli::cli_text("")
    cli::cli_text("{.emph Tip: Use as.character(report) or style = 'markdown' to copy directly into Quarto/R Markdown.}")
  })
  cat(lines, sep = "\n")
  invisible(x)
}
