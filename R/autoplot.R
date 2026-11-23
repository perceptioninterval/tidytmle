#' Diagnostic Plots for a TMLE Fit
#'
#' Produces `ggplot2` causal diagnostics for a fitted `"tmle_fit"` object,
#' turning every step of the TMLE estimator into an inspectable visual check:
#' - `"overlap"`: **Step 2 (Positivity / Overlap Check)**. Mirrored propensity
#'   score density by treatment arm with `g_bounds` truncation lines (or
#'   density-ratio distribution for LMTPs). Validates whether both treated and
#'   control units exist across confounder strata (Petersen et al., 2012).
#' - `"targeting"`: **Steps 3–4 (Initial G-Computation vs. Targeted TMLE)**.
#'   Compares initial ML counterfactual contrasts (\eqn{\bar{Q}^0(1, W) - \bar{Q}^0(0, W)})
#'   against fluctuated TMLE contrasts (\eqn{\bar{Q}^*(1, W) - \bar{Q}^*(0, W)}),
#'   shaded by Clever Covariate magnitude \eqn{|H_n(A, W)|}, illustrating how the
#'   fluctuation step removes residual confounding bias from the initial fit.
#' - `"learners"`: **Cross-Validated Super Learner Ensemble**. Compares CV risk
#'   and simplex metalearner weights across base workflows for \eqn{Q}, \eqn{g},
#'   and \eqn{C}.
#' - `"eif"`: **Step 5 (Efficient Influence Function Variance Leverage)**.
#'   Observation-level share of total EIF variance (\eqn{D_i^2 / \sum D_j^2})
#'   to catch single-row outliers driving standard errors.
#' - `"sensitivity"`: **Positivity Sensitivity Sweep**. Instantaneous
#'   `retarget()` sweep across truncation thresholds `g_grid` showing stability
#'   of point estimates, 95% confidence intervals, and sample truncation rates.
#' - `"balance"`: **Weighted Covariate Balance Love Plot**. Compares unadjusted
#'   vs. IPW/TMLE-weighted absolute standardized mean differences (`|SMD|`)
#'   against the standard `0.10` balance benchmark.
#'
#' @param object A fitted `"tmle_fit"` object.
#' @param type Diagnostic plot type: `"overlap"`, `"targeting"`, `"learners"`,
#'   `"eif"`, `"sensitivity"`, or `"balance"`.
#' @param g_grid Numeric vector of lower positivity thresholds in `(0, 0.5)`
#'   swept when `type = "sensitivity"`.
#' @param ... Currently unused.
#'
#' @references
#' Schuler, M. S., & Rose, S. (2017). Targeted maximum likelihood estimation for
#' causal inference in observational studies. *American Journal of Epidemiology*,
#' 185(1), 65–73.
#'
#' Luque-Fernandez, M. A., Schomaker, M., Rachet, B., & Schnitzer, M. E.
#' (2018). Targeted maximum likelihood estimation for a binary treatment: A
#' tutorial. *Statistics in Medicine*, 37(16), 2530–2546.
#'
#' Petersen, M. L., Porter, K. E., Gruber, S., Wang, Y., & van der Laan, M. J.
#' (2012). Diagnosing and responding to violations in the positivity assumption.
#' *Statistical Methods in Medical Research*, 21(1), 31–54.
#'
#' @return A `ggplot2` object.
#' @method autoplot tmle_fit
#' @export
autoplot.tmle_fit <- function(object,
                              type = c("overlap", "targeting", "learners", "eif", "sensitivity", "balance"),
                              g_grid = seq(0.01, 0.15, by = 0.01),
                              ...) {
  type <- match.arg(type)

  if (type == "overlap") {
    return(autoplot_overlap(object))
  }
  if (type == "targeting") {
    return(autoplot_targeting(object))
  }
  if (type == "learners") {
    return(autoplot_learners(object))
  }
  if (type == "eif") {
    return(autoplot_eif(object))
  }
  if (type == "sensitivity") {
    return(autoplot_sensitivity(object, g_grid = g_grid))
  }
  autoplot_balance(object)
}

autoplot_targeting <- function(object) {
  t_res <- object$targeting
  idx <- object$subgroup_idx
  w <- object$cache$weights[idx]
  w <- w / mean(w)

  if (object$cache$treatment_type == "binary") {
    x_init <- t_res$q1_init - t_res$q0_init
    y_star <- t_res$q1_star - t_res$q0_star
    plug_mean <- stats::weighted.mean(x_init, w = w)
    tmle_mean <- stats::weighted.mean(y_star, w = w)
    x_lab <- "Initial ML Individual Contrast Q0(1, W) - Q0(0, W)"
    y_lab <- "Targeted TMLE Contrast Q*(1, W) - Q*(0, W)"
    title_txt <- "TMLE Targeting Step: Initial G-Computation vs. Fluctuated Contrasts"
    sub_txt <- sprintf(
      "Sample-average plug-in ATE = %.4f -> Targeted TMLE ATE = %.4f (bias-correction shift = %+.4f)",
      plug_mean, tmle_mean, tmle_mean - plug_mean
    )
  } else {
    x_init <- t_res$q1_init
    y_star <- t_res$q1_star
    plug_mean <- stats::weighted.mean(x_init, w = w)
    tmle_mean <- stats::weighted.mean(y_star, w = w)
    x_lab <- "Initial ML Shifted Prediction Q0(d(A, W), W)"
    y_lab <- "Targeted TMLE Prediction Q*(d(A, W), W)"
    title_txt <- "TMLE Targeting Step: Initial vs. Fluctuated Policy Outcomes"
    sub_txt <- sprintf(
      "Plug-in policy mean = %.4f -> Targeted TMLE = %.4f (shift = %+.4f)",
      plug_mean, tmle_mean, tmle_mean - plug_mean
    )
  }

  df_t <- tibble::tibble(
    x_init = x_init,
    y_star = y_star,
    h_mag = abs(t_res$clever_h)
  )

  ggplot2::ggplot(df_t, ggplot2::aes(x = x_init, y = y_star, color = h_mag)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#737A82", linewidth = 0.7) +
    ggplot2::geom_point(alpha = 0.72, size = 1.9) +
    ggplot2::scale_color_gradient(low = "#3182ce", high = "#CA225E") +
    ggplot2::labs(
      title = title_txt,
      subtitle = sub_txt,
      x = x_lab,
      y = y_lab,
      color = "|Clever H|"
    ) +
    ggplot2::theme_minimal()
}

autoplot_overlap <- function(object) {
  idx <- object$subgroup_idx
  t_res <- object$targeting
  bounds <- t_res$g_info$bounds

  if (object$cache$treatment_type == "continuous") {
    df_p <- tibble::tibble(g = t_res$g_trunc)
    return(
      ggplot2::ggplot(df_p, ggplot2::aes(x = g)) +
        ggplot2::geom_histogram(bins = 30L, fill = "#2b6cb0", alpha = 0.8, color = "white") +
        ggplot2::geom_vline(xintercept = bounds[2L], linetype = "dashed", color = "#c53030", linewidth = 0.8) +
        ggplot2::labs(
          title = "Modified Treatment Policy Density-Ratio Overlap",
          subtitle = sprintf("Truncated at max ratio %.2f (%.1f%% clipped)", bounds[2L], 100 * t_res$g_info$pct_truncated),
          x = "Truncated Density Ratio r(A, W)",
          y = "Count"
        ) +
        ggplot2::theme_minimal()
    )
  }

  a_sub <- object$cache$a_num[idx]
  df_p <- tibble::tibble(
    g = t_res$g_hat,
    arm = factor(ifelse(a_sub == 1L, "Treated (A = 1)", "Control (A = 0)"),
                 levels = c("Treated (A = 1)", "Control (A = 0)"))
  )

  ggplot2::ggplot(df_p, ggplot2::aes(x = g)) +
    ggplot2::geom_density(
      data = dplyr::filter(df_p, arm == "Treated (A = 1)"),
      ggplot2::aes(y = ggplot2::after_stat(density), fill = arm),
      alpha = 0.65, color = "#1f4e78"
    ) +
    ggplot2::geom_density(
      data = dplyr::filter(df_p, arm == "Control (A = 0)"),
      ggplot2::aes(y = -ggplot2::after_stat(density), fill = arm),
      alpha = 0.65, color = "#9c3d10"
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = bounds, linetype = "dashed", color = "#c53030", linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = c("Treated (A = 1)" = "#3182ce", "Control (A = 0)" = "#ed8936")) +
    ggplot2::labs(
      title = "Propensity Score Overlap & Positivity Bounds",
      subtitle = sprintf("g_bounds = [%.3f, %.3f] | %.1f%% truncated", bounds[1L], bounds[2L], 100 * t_res$g_info$pct_truncated),
      x = "Estimated Propensity Score P(A = 1 | W)",
      y = "Mirrored Density",
      fill = "Treatment Arm"
    ) +
    ggplot2::theme_minimal()
}

autoplot_learners <- function(object) {
  w_tbl <- extract_sl_weights(object) |>
    dplyr::mutate(metric_type = "Super Learner Weight") |>
    dplyr::rename(value = weight)

  m_tbl <- collect_metrics.tmle_fit(object, summarize = TRUE)
  # Pick primary metric per nuisance
  m_primary <- m_tbl |>
    dplyr::filter(
      (nuisance == "Q" & .metric %in% c("rmse", "mn_log_loss")) |
        (nuisance %in% c("g", "C") & .metric == "mn_log_loss")
    ) |>
    dplyr::mutate(metric_type = paste0("CV ", .metric)) |>
    dplyr::select(nuisance, learner, value = mean, metric_type)

  plot_df <- dplyr::bind_rows(w_tbl, m_primary)
  ggplot2::ggplot(plot_df, ggplot2::aes(x = stats::reorder(learner, value), y = value, fill = nuisance)) +
    ggplot2::geom_col(show.legend = FALSE, alpha = 0.85) +
    ggplot2::coord_flip() +
    ggplot2::facet_grid(nuisance ~ metric_type, scales = "free") +
    ggplot2::labs(
      title = "Nuisance Super Learner Weights & Cross-Validated Risk",
      x = "Learner / Ensemble",
      y = NULL
    ) +
    ggplot2::theme_minimal()
}

autoplot_eif <- function(object) {
  eif <- object$targeting$eif
  lev <- object$targeting$eif_leverage
  df <- tibble::tibble(
    obs_id = seq_along(eif),
    eif = eif,
    leverage = lev,
    is_outlier = lev > 0.05
  )

  ggplot2::ggplot(df, ggplot2::aes(x = obs_id, y = leverage, color = is_outlier)) +
    ggplot2::geom_segment(ggplot2::aes(xend = obs_id, y = 0, yend = leverage), alpha = 0.65) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_hline(yintercept = 0.10, linetype = "dashed", color = "#c53030") +
    ggplot2::scale_color_manual(values = c("FALSE" = "#2b6cb0", "TRUE" = "#c53030"),
                                labels = c("<= 5%", "> 5% EIF variance")) +
    ggplot2::labs(
      title = "Observation-Level Efficient Influence Function (EIF) Variance Leverage",
      subtitle = sprintf("Max leverage = %.1f%% (dashed warning threshold = 10%%)", 100 * max(lev)),
      x = "Observation Index",
      y = "Share of Total EIF Variance (D_i^2 / sum D^2)",
      color = "Leverage Flag"
    ) +
    ggplot2::theme_minimal()
}

autoplot_sensitivity <- function(object, g_grid = seq(0.01, 0.15, by = 0.01)) {
  primary_est <- object$spec$estimand$names[1L]
  sweep_df <- purrr::map_dfr(g_grid, function(b) {
    rf <- retarget(object, g_bounds = c(b, 1 - b))
    td <- rf$targeting$estimates |>
      dplyr::filter(estimand == primary_est) |>
      dplyr::slice(1L)
    td$threshold <- b
    td$pct_truncated <- 100 * rf$targeting$g_info$pct_truncated
    td
  })

  ggplot2::ggplot(sweep_df, ggplot2::aes(x = threshold, y = estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high), fill = "#3182ce", alpha = 0.2) +
    ggplot2::geom_line(color = "#2b6cb0", linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(size = pct_truncated), color = "#1a365d") +
    ggplot2::labs(
      title = sprintf("Positivity Sensitivity Sweep via Instant retarget() (%s)", toupper(primary_est)),
      subtitle = "Estimates & 95% CIs across propensity score truncation thresholds",
      x = "Lower Propensity Bound (b, with upper = 1 - b)",
      y = sprintf("Estimated %s", sweep_df$term[1L]),
      size = "% Truncated"
    ) +
    ggplot2::theme_minimal()
}

autoplot_balance <- function(object) {
  if (object$cache$treatment_type != "binary") {
    cli::cli_abort("Covariate balance Love plot is only available for binary treatment estimands.")
  }
  idx <- object$subgroup_idx
  df <- object$cache$data[idx, object$cache$covariate_cols, drop = FALSE]
  a <- object$cache$a_num[idx]
  w_base <- object$cache$weights[idx]
  g_t <- object$targeting$g_trunc

  # IPW weights using truncated propensity
  w_ipw <- w_base * ifelse(a == 1L, 1 / g_t, 1 / (1 - g_t))

  # Expand numeric & factor covariates into numeric indicators cleanly
  mm <- stats::model.matrix(~ . - 1, data = df)
  calc_smd <- function(x, a_vec, wt) {
    m1 <- stats::weighted.mean(x[a_vec == 1L], wt[a_vec == 1L], na.rm = TRUE)
    m0 <- stats::weighted.mean(x[a_vec == 0L], wt[a_vec == 0L], na.rm = TRUE)
    v1 <- stats::var(x[a_vec == 1L], na.rm = TRUE)
    v0 <- stats::var(x[a_vec == 0L], na.rm = TRUE)
    pooled_sd <- sqrt(0.5 * (v1 + v0))
    if (!is.finite(pooled_sd) || pooled_sd <= 1e-12) return(0)
    abs(m1 - m0) / pooled_sd
  }

  cov_names <- colnames(mm)
  bal_df <- purrr::map_dfr(cov_names, function(cn) {
    xv <- mm[, cn]
    tibble::tibble(
      covariate = cn,
      group = c("Unadjusted", "TMLE / IPW Weighted"),
      smd = c(calc_smd(xv, a, w_base), calc_smd(xv, a, w_ipw))
    )
  })

  ggplot2::ggplot(bal_df, ggplot2::aes(x = smd, y = stats::reorder(covariate, smd), color = group, shape = group)) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_vline(xintercept = 0.10, linetype = "dashed", color = "#c53030") +
    ggplot2::scale_color_manual(values = c("Unadjusted" = "#dd6b20", "TMLE / IPW Weighted" = "#2b6cb0")) +
    ggplot2::labs(
      title = "Covariate Balance Love Plot (Unadjusted vs Weighted)",
      subtitle = "Dashed line marks standard |SMD| = 0.10 balance threshold",
      x = "Absolute Standardized Mean Difference |SMD|",
      y = NULL,
      color = "Sample",
      shape = "Sample"
    ) +
    ggplot2::theme_minimal()
}
