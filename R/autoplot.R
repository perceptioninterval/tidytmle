#' Internal Palette, Formatting, Balance, and Editorial Tufte Theme Helpers
#' @noRd
tmle_palette <- function() {
  list(
    crimson = "#CA225E",
    navy = "#1D4E89",
    amber = "#D97706",
    teal = "#0F766E",
    teal_bg = "#ECFDF5",
    red = "#B91C1C",
    red_bg = "#FEF2F2",
    ink = "#1E293B",
    muted = "#475569",
    light = "#94A3B8",
    grid = "#E2E8F0",
    axis = "#CBD5E1"
  )
}

fmt_pct <- function(x, digits = 0L) {
  out <- rep("", length(x))
  ok <- is.finite(x)
  out[ok] <- sprintf(paste0("%.", as.integer(digits), "f%%"), 100 * x[ok])
  out
}

format_p_human <- function(p) {
  vapply(p, function(val) {
    if (is.na(val) || !is.finite(val)) return("p = NA")
    if (val < 1e-4) return("p < 0.0001")
    sprintf("p = %.4f", val)
  }, character(1L))
}

theme_tidytmle <- function(base_size = 11,
                           grid = c("y", "x", "both", "none"),
                           compact = FALSE) {
  grid <- match.arg(grid)
  pal <- tmle_palette()

  grid_line <- ggplot2::element_line(color = pal$grid, linewidth = 0.35)
  blank <- ggplot2::element_blank()

  th <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold",
        color = pal$ink,
        size = if (compact) base_size * 0.98 else base_size * 1.15,
        margin = ggplot2::margin(b = if (compact) 2 else 4)
      ),
      plot.subtitle = ggplot2::element_text(
        color = pal$muted,
        size = if (compact) base_size * 0.80 else base_size * 0.92,
        margin = ggplot2::margin(b = if (compact) 6 else 10)
      ),
      plot.caption = ggplot2::element_text(
        color = "#64748B",
        size = if (compact) base_size * 0.72 else base_size * 0.80,
        hjust = 0,
        margin = ggplot2::margin(t = 6)
      ),
      axis.title = ggplot2::element_text(
        color = pal$ink,
        size = if (compact) base_size * 0.82 else base_size * 0.92
      ),
      axis.text = ggplot2::element_text(
        color = pal$muted,
        size = if (compact) base_size * 0.78 else base_size * 0.86
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        color = pal$ink,
        size = if (compact) base_size * 0.82 else base_size * 0.92
      ),
      panel.grid.minor = blank,
      panel.grid.major.x = if (grid %in% c("x", "both")) grid_line else blank,
      panel.grid.major.y = if (grid %in% c("y", "both")) grid_line else blank,
      axis.line = ggplot2::element_line(color = pal$axis, linewidth = 0.4),
      legend.title = ggplot2::element_text(
        face = "bold",
        color = pal$ink,
        size = if (compact) base_size * 0.78 else base_size * 0.86
      ),
      legend.text = ggplot2::element_text(
        color = pal$muted,
        size = if (compact) base_size * 0.75 else base_size * 0.84
      ),
      plot.margin = if (compact) {
        ggplot2::margin(8, 10, 6, 8)
      } else {
        ggplot2::margin(12, 16, 10, 12)
      }
    )
  th
}

compute_covariate_balance <- function(object) {
  empty_bal <- tibble::tibble(
    covariate = character(),
    smd_raw = numeric(),
    smd_wt = numeric(),
    improved = logical()
  )
  if (!identical(object$cache$treatment_type, "binary")) {
    return(empty_bal)
  }
  idx <- object$subgroup_idx
  df <- object$cache$data[idx, object$cache$covariate_cols, drop = FALSE]
  a <- object$cache$a_num[idx]
  w_base <- object$cache$weights[idx]
  g_t <- object$targeting$g_trunc
  w_ipw <- w_base * ifelse(a == 1L, 1 / g_t, 1 / (1 - g_t))

  # Expand column-by-column so NA covariates keep exact row alignment and
  # 1-level factors/characters in subgroups are skipped cleanly instead of erroring.
  cols_expanded <- list()
  for (cn in names(df)) {
    col <- df[[cn]]
    if (is.numeric(col) || is.integer(col) || is.logical(col)) {
      v <- as.numeric(col)
      if (length(unique(v[is.finite(v)])) >= 2L) {
        cols_expanded[[cn]] <- v
      }
    } else if (is.factor(col) || is.character(col)) {
      obs_vals <- as.character(col[!is.na(col)])
      lvls <- if (is.factor(col)) {
        intersect(levels(col), unique(obs_vals))
      } else {
        sort(unique(obs_vals))
      }
      if (length(lvls) >= 2L) {
        col_chr <- as.character(col)
        for (lv in lvls[-1L]) {
          cols_expanded[[paste0(cn, lv)]] <- ifelse(is.na(col_chr), NA_real_, as.numeric(col_chr == lv))
        }
      }
    }
  }

  if (length(cols_expanded) == 0L) {
    return(empty_bal)
  }

  calc_smd <- function(x, a_vec, wt) {
    ok <- is.finite(x) & !is.na(a_vec) & is.finite(wt)
    i1 <- ok & (a_vec == 1L)
    i0 <- ok & (a_vec == 0L)
    if (sum(i1) < 2L || sum(i0) < 2L) return(0)
    m1 <- stats::weighted.mean(x[i1], wt[i1])
    m0 <- stats::weighted.mean(x[i0], wt[i0])
    v1 <- stats::var(x[i1])
    v0 <- stats::var(x[i0])
    pooled_sd <- sqrt(0.5 * (v1 + v0))
    if (!is.finite(pooled_sd) || pooled_sd <= 1e-12) return(0)
    abs(m1 - m0) / pooled_sd
  }

  purrr::map_dfr(names(cols_expanded), function(cn) {
    xv <- cols_expanded[[cn]]
    s_raw <- calc_smd(xv, a, w_base)
    s_wt <- calc_smd(xv, a, w_ipw)
    tibble::tibble(
      covariate = cn,
      smd_raw = s_raw,
      smd_wt = s_wt,
      improved = s_wt <= s_raw + 1e-12
    )
  })
}

#' Human-Centered Causal Visual Diagnostics for a TMLE Fit
#'
#' Produces clean, Tufte-inspired `ggplot2` causal figures turning every step
#' of the TMLE estimator into an intuitive visual check:
#' - `"effect"`: **Causal Effect Forest Plot (Default)**. Publication-ready
#'   effect estimates and 95% confidence intervals with direct numeric labels, a
#'   `"No effect"` reference line, and (for binary ATE) a hollow marker showing
#'   the untargeted initial ML estimate before bias calibration.
#' - `"overlap"`: **Step 2 (Fair Comparison / Positivity Check)**. Mirrored
#'   propensity distributions directly labeled `"Treated"` and `"Control"`,
#'   with shaded tail clipping zones (or trimmed density-ratio distribution for
#'   continuous LMTPs).
#' - `"targeting"`: **Steps 3–4 (How Much Bias Calibration Changed Things)**.
#'   Initial ML counterfactual contrasts plotted against targeted TMLE contrasts,
#'   with direct `"No adjustment"` reference line and population shift marker.
#' - `"balance"`: **Confounder Balance Dumbbell Plot**. Before vs. after IPW
#'   imbalance (`|SMD|`) ordered from largest baseline imbalance to smallest,
#'   highlighting the balanced `< 0.10` zone without duplicate factor mirror rows.
#' - `"sensitivity"`: **Positivity Strictness Sweep**. Instant `retarget()`
#'   sensitivity curve across clipping thresholds showing stability of point
#'   estimates and 95% confidence intervals.
#' - `"learners"`: **Cross-Validated Model Weights & Risk**. Clean paired panel
#'   comparing Super Learner weights and cross-validated error across candidate
#'   workflows for Outcome (Q), Treatment (g), and Missingness (C).
#' - `"eif"`: **Step 5 (Outlier Influence Check)**. Observation-level share of
#'   total variance (`%`) flagging single rows exceeding the `10%` caution line.
#' - `"dashboard"`: **4-Panel Causal Audit Dashboard**. Uncluttered 2x2 view
#'   combining Overlap, Targeting, Sensitivity, and Balance (or Influence for
#'   continuous policies) via `patchwork`.
#'
#' @param object A fitted `"tmle_fit"` object.
#' @param type Diagnostic plot type: `"effect"` (default), `"overlap"`,
#'   `"targeting"`, `"balance"`, `"sensitivity"`, `"learners"`, `"eif"`, or
#'   `"dashboard"`.
#' @param g_grid Numeric vector of lower positivity thresholds in `(0, 0.5)`
#'   swept when `type = "sensitivity"` or `"dashboard"`.
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
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, outcome_type = "binary", seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = c("ate", "rr")) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' autoplot(fit)
#' autoplot(fit, type = "overlap")
#'
#' @export
autoplot.tmle_fit <- function(object,
                              type = c("effect", "overlap", "targeting", "balance", "sensitivity", "learners", "eif", "dashboard"),
                              g_grid = seq(0.01, 0.15, by = 0.01),
                              ...) {
  type <- match.arg(type)

  if (type == "effect") {
    return(autoplot_effect(object))
  }
  if (type == "overlap") {
    return(autoplot_overlap(object))
  }
  if (type == "targeting") {
    return(autoplot_targeting(object))
  }
  if (type == "balance") {
    return(autoplot_balance(object))
  }
  if (type == "sensitivity") {
    return(autoplot_sensitivity(object, g_grid = g_grid))
  }
  if (type == "learners") {
    return(autoplot_learners(object))
  }
  if (type == "eif") {
    return(autoplot_eif(object))
  }
  autoplot_dashboard(object, g_grid = g_grid)
}

autoplot_effect <- function(object, compact = FALSE) {
  pal <- tmle_palette()
  gl <- glance.tmle_fit(object)
  td <- object$targeting$estimates

  human_term <- function(est, trm) {
    switch(
      est,
      ate = "Average Effect (ATE)",
      att = "Effect Among Treated (ATT)",
      atc = "Effect Among Untreated (ATC)",
      rr  = "Risk Ratio (RR)",
      or  = "Odds Ratio (OR)",
      ey1 = "Mean If Everyone Treated",
      ey0 = "Mean If No One Treated",
      paste("Policy:", trm)
    )
  }

  is_ratio <- td$estimand %in% c("rr", "or")
  df_e <- td |>
    dplyr::mutate(
      is_obs_contrast = grepl(" - observed$", estimand),
      scale_group = ifelse(
        estimand %in% c("rr", "or"),
        "Ratio Scale (1 = no difference)",
        ifelse(
          estimand %in% c("ey1", "ey0") | (startsWith(estimand, "shift(") & !is_obs_contrast) | (startsWith(estimand, "scale_tx(") & !is_obs_contrast),
          "Counterfactual Means",
          "Difference Scale (0 = no difference)"
        )
      ),
      null_ref = ifelse(
        estimand %in% c("rr", "or"), 1,
        ifelse(scale_group == "Difference Scale (0 = no difference)", 0, NA_real_)
      ),
      label_human = vapply(seq_len(dplyr::n()), function(i) human_term(estimand[i], term[i]), character(1L)),
      direct_lbl = ifelse(
        estimand %in% c("rr", "or"),
        sprintf("%.2fx [%.2f, %.2f]", estimate, conf.low, conf.high),
        ifelse(
          scale_group == "Difference Scale (0 = no difference)",
          sprintf("%+.3f [%.3f, %.3f]", estimate, conf.low, conf.high),
          sprintf("%.3f [%.3f, %.3f]", estimate, conf.low, conf.high)
        )
      )
    )
  df_e$label_human <- factor(df_e$label_human, levels = rev(unique(df_e$label_human)))

  primary_est <- df_e$estimand[1L]
  show_plugin <- identical(object$cache$treatment_type, "binary") &&
    primary_est %in% c("ate", "att", "atc")
  df_plugin <- if (show_plugin) {
    p_row <- df_e[1L, , drop = FALSE]
    p_row$plugin_est <- p_row$estimate - gl$plugin_vs_tmle_shift
    p_row
  } else {
    NULL
  }

  ci_pct <- round(100 * (if (!is.null(object$targeting$conf_level)) object$targeting$conf_level else 0.95))

  n_scales <- length(unique(df_e$scale_group))
  ref_df <- df_e |>
    dplyr::filter(!is.na(null_ref)) |>
    dplyr::distinct(scale_group, null_ref)

  p <- ggplot2::ggplot(df_e, ggplot2::aes(x = estimate, y = label_human)) +
    ggplot2::geom_vline(
      data = ref_df,
      ggplot2::aes(xintercept = null_ref),
      linetype = "dashed", color = pal$light, linewidth = 0.65
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = conf.low, xend = conf.high, y = label_human, yend = label_human),
      linewidth = 1.15, color = pal$navy, lineend = "round"
    ) +
    ggplot2::geom_point(size = 3.4, color = pal$crimson) +
    ggplot2::geom_text(
      ggplot2::aes(x = conf.high, label = paste0("  ", direct_lbl)),
      hjust = 0, size = 3.25, fontface = "bold", color = pal$ink
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.36)))

  if (!is.null(df_plugin)) {
    p <- p + ggplot2::geom_point(
      data = df_plugin,
      ggplot2::aes(x = plugin_est, y = label_human),
      shape = 23, size = 2.8, fill = "white", color = pal$muted, stroke = 1.0
    )
  }

  if (n_scales > 1L) {
    p <- p + ggplot2::facet_wrap(~scale_group, scales = "free", ncol = 1L)
  }

  sub_txt <- sprintf(
    "n = %d (%d-fold %s) | Targeting shifted initial ML %s by %+.4f",
    gl$n, gl$v_folds, if (gl$cv_tmle) "CV-TMLE" else "TMLE", td$term[1L], gl$plugin_vs_tmle_shift
  )
  cap_txt <- if (show_plugin) {
    sprintf("Solid crimson dot = targeted TMLE estimate with %d%% CI; hollow diamond = initial untargeted ML plug-in.", ci_pct)
  } else {
    sprintf("Solid crimson dot = targeted TMLE estimate with %d%% confidence interval.", ci_pct)
  }

  p +
    ggplot2::labs(
      title = "Estimated Causal Effect & Uncertainty",
      subtitle = sub_txt,
      x = sprintf("Causal Estimate (%d%% Confidence Interval)", ci_pct),
      y = NULL,
      caption = cap_txt
    ) +
    theme_tidytmle(grid = "x", compact = compact) +
    ggplot2::coord_cartesian(clip = "off")
}

autoplot_dashboard <- function(object, g_grid = seq(0.01, 0.15, by = 0.01)) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg patchwork} is required for {.code autoplot(fit, type = \"dashboard\")}.",
      "i" = "Install it with {.code install.packages(\"patchwork\")} or plot individual types directly."
    ))
  }

  pal <- tmle_palette()
  gl <- glance.tmle_fit(object)
  td1 <- object$targeting$estimates[1L, ]
  est_headline <- if (td1$estimand %in% c("rr", "or")) {
    sprintf("%s = %.2fx [%.2f, %.2f]", td1$term, td1$estimate, td1$conf.low, td1$conf.high)
  } else {
    sprintf("%s = %+.3f [%.3f, %.3f]", td1$term, td1$estimate, td1$conf.low, td1$conf.high)
  }

  p1 <- autoplot_overlap(object, compact = TRUE) +
    ggplot2::labs(title = "1. Fair Comparison Overlap", subtitle = NULL, caption = NULL)
  p2 <- autoplot_targeting(object, compact = TRUE) +
    ggplot2::labs(title = "2. Targeting Bias Adjustment", subtitle = NULL, caption = NULL)
  p3 <- autoplot_sensitivity(object, g_grid = g_grid, compact = TRUE) +
    ggplot2::labs(title = "3. Stability Across Strictness Bounds", subtitle = NULL, caption = NULL)
  p4 <- if (object$cache$treatment_type == "binary") {
    autoplot_balance(object, compact = TRUE) +
      ggplot2::labs(title = "4. Confounder Balance (Before -> After)", subtitle = NULL, caption = NULL)
  } else {
    autoplot_eif(object, compact = TRUE) +
      ggplot2::labs(title = "4. Individual Outlier Influence", subtitle = NULL, caption = NULL)
  }

  patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2L) +
    patchwork::plot_annotation(
      title = "TMLE Causal Audit Dashboard",
      subtitle = sprintf(
        "%s  |  n = %d (%d-fold %s)  |  %.1f%% clipped  |  mean EIF = %.1e",
        est_headline, gl$n, gl$v_folds,
        if (gl$cv_tmle) "CV-TMLE" else "TMLE",
        100 * gl$pct_truncated, gl$mean_eif
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", color = pal$ink, size = 14),
        plot.subtitle = ggplot2::element_text(color = pal$muted, size = 10, margin = ggplot2::margin(b = 8))
      )
    )
}

autoplot_targeting <- function(object, compact = FALSE) {
  pal <- tmle_palette()
  t_res <- object$targeting
  idx <- object$subgroup_idx
  w <- object$cache$weights[idx]
  w <- w / mean(w)

  if (object$cache$treatment_type == "binary") {
    x_init <- t_res$q1_init - t_res$q0_init
    y_star <- t_res$q1_star - t_res$q0_star
    plug_mean <- stats::weighted.mean(x_init, w = w)
    tmle_mean <- stats::weighted.mean(y_star, w = w)
    x_lab <- if (compact) "Initial ML Effect" else "Initial ML Effect Contrast Q0(1, W) - Q0(0, W)"
    y_lab <- if (compact) "Targeted TMLE Effect" else "Targeted TMLE Contrast Q*(1, W) - Q*(0, W)"
    title_txt <- "How Targeting Nudge Removed Bias"
    sub_txt <- sprintf(
      "Initial ML avg = %.4f -> Targeted TMLE avg = %.4f (bias calibration shift = %+.4f)",
      plug_mean, tmle_mean, tmle_mean - plug_mean
    )
  } else {
    x_init <- t_res$q1_init
    y_star <- t_res$q1_star
    plug_mean <- stats::weighted.mean(x_init, w = w)
    tmle_mean <- stats::weighted.mean(y_star, w = w)
    x_lab <- if (compact) "Initial ML Policy" else "Initial ML Shifted Prediction Q0(d(A, W), W)"
    y_lab <- if (compact) "Targeted TMLE Policy" else "Targeted TMLE Prediction Q*(d(A, W), W)"
    title_txt <- "How Targeting Nudge Removed Policy Bias"
    sub_txt <- sprintf(
      "Initial avg = %.4f -> Targeted TMLE avg = %.4f (shift = %+.4f)",
      plug_mean, tmle_mean, tmle_mean - plug_mean
    )
  }

  df_t <- tibble::tibble(
    x_init = x_init,
    y_star = y_star,
    h_mag = abs(t_res$clever_h)
  )
  rng_x <- range(df_t$x_init, na.rm = TRUE)
  diag_x <- rng_x[1L] + 0.18 * diff(rng_x)

  ggplot2::ggplot(df_t, ggplot2::aes(x = x_init, y = y_star, color = h_mag)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = pal$light, linewidth = 0.7) +
    ggplot2::geom_point(alpha = 0.72, size = if (compact) 1.55 else 1.95) +
    ggplot2::annotate(
      "point", x = plug_mean, y = tmle_mean,
      shape = 23, size = if (compact) 3.0 else 3.6,
      fill = pal$crimson, color = "white", stroke = 0.9
    ) +
    ggplot2::annotate(
      "text", x = diag_x, y = diag_x,
      label = "No adjustment", angle = 0, vjust = 1.4, hjust = 0,
      size = if (compact) 2.7 else 3.1, color = pal$muted, fontface = "italic"
    ) +
    ggplot2::scale_color_gradientn(colors = c(pal$light, pal$navy, pal$crimson)) +
    ggplot2::labs(
      title = title_txt,
      subtitle = sub_txt,
      x = x_lab,
      y = y_lab,
      color = if (compact) "|H|" else "Calibration\nWeight |H|",
      caption = if (!compact) "Points off the dashed line moved during targeting; crimson diamond marks population averages." else NULL
    ) +
    theme_tidytmle(grid = "both", compact = compact) +
    ggplot2::theme(legend.position = if (compact) "right" else "right")
}

autoplot_overlap <- function(object, compact = FALSE) {
  pal <- tmle_palette()
  idx <- object$subgroup_idx
  t_res <- object$targeting
  bounds <- t_res$g_info$bounds

  if (object$cache$treatment_type == "continuous") {
    g_vals <- t_res$g_trunc
    pct_clip <- 100 * t_res$g_info$pct_truncated
    show_cap <- (max(g_vals, na.rm = TRUE) >= 0.75 * bounds[2L]) || (t_res$g_info$pct_truncated > 0)
    x_lim_hi <- if (show_cap) bounds[2L] * 1.05 else max(2.5, max(g_vals, na.rm = TRUE) * 1.15)
    df_p <- tibble::tibble(g = g_vals)
    sub_cont <- sprintf(
      "Median weight = %.2fx | Safety cap = %.1fx (%.1f%% clipped)",
      stats::median(g_vals, na.rm = TRUE), bounds[2L], pct_clip
    )
    p <- ggplot2::ggplot(df_p, ggplot2::aes(x = g)) +
      ggplot2::geom_histogram(bins = 28L, fill = pal$navy, alpha = 0.82, color = "white") +
      ggplot2::geom_vline(xintercept = 1, linetype = "dotted", color = pal$muted, linewidth = 0.7) +
      ggplot2::annotate(
        "text", x = 1, y = Inf, label = " 1x (no reweight)",
        vjust = 1.5, hjust = 0, size = if (compact) 2.7 else 3.1, color = pal$muted
      ) +
      ggplot2::coord_cartesian(xlim = c(0, x_lim_hi)) +
      ggplot2::labs(
        title = "Continuous Policy Weight Overlap",
        subtitle = sub_cont,
        x = if (compact) "Density Ratio Weight r(A, W)" else "Reweighting Factor for Shifted Dose r(A, W)",
        y = "Individuals"
      ) +
      theme_tidytmle(grid = "y", compact = compact)
    if (show_cap) {
      p <- p + ggplot2::geom_vline(xintercept = bounds[2L], linetype = "dashed", color = pal$red, linewidth = 0.8)
    }
    return(p)
  }

  a_sub <- object$cache$a_num[idx]
  n_tr <- sum(a_sub == 1L)
  n_ct <- sum(a_sub == 0L)
  df_p <- tibble::tibble(
    g = t_res$g_hat,
    arm = factor(ifelse(a_sub == 1L, "Treated", "Control"), levels = c("Treated", "Control"))
  )

  d_tr <- stats::density(df_p$g[df_p$arm == "Treated"], from = 0, to = 1, n = 256L)
  d_ct <- stats::density(df_p$g[df_p$arm == "Control"], from = 0, to = 1, n = 256L)
  y_top <- max(d_tr$y, 1) * 0.92
  y_bot <- -max(d_ct$y, 1) * 0.92
  x_lbl <- min(0.22, bounds[1L] + 0.025)

  p <- ggplot2::ggplot(df_p, ggplot2::aes(x = g)) +
    ggplot2::annotate("rect", xmin = 0, xmax = bounds[1L], ymin = -Inf, ymax = Inf, fill = pal$red_bg, alpha = 0.85) +
    ggplot2::annotate("rect", xmin = bounds[2L], xmax = 1, ymin = -Inf, ymax = Inf, fill = pal$red_bg, alpha = 0.85) +
    ggplot2::geom_density(
      data = dplyr::filter(df_p, arm == "Treated"),
      ggplot2::aes(y = ggplot2::after_stat(density)),
      fill = pal$navy, alpha = 0.72, color = pal$navy, linewidth = 0.65, show.legend = FALSE
    ) +
    ggplot2::geom_density(
      data = dplyr::filter(df_p, arm == "Control"),
      ggplot2::aes(y = -ggplot2::after_stat(density)),
      fill = pal$amber, alpha = 0.72, color = pal$amber, linewidth = 0.65, show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, color = pal$ink, linewidth = 0.45) +
    ggplot2::geom_vline(xintercept = bounds, linetype = "dashed", color = pal$red, linewidth = 0.75) +
    ggplot2::annotate(
      "text", x = x_lbl, y = y_top,
      label = sprintf("Treated (n = %d)", n_tr),
      hjust = 0, vjust = 1, fontface = "bold",
      size = if (compact) 3.0 else 3.6, color = pal$navy
    ) +
    ggplot2::annotate(
      "text", x = x_lbl, y = y_bot,
      label = sprintf("Control (n = %d)", n_ct),
      hjust = 0, vjust = 0, fontface = "bold",
      size = if (compact) 3.0 else 3.6, color = pal$amber
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = fmt_pct
    ) +
    ggplot2::labs(
      title = "Who Is Comparable? Propensity Overlap Check",
      subtitle = sprintf(
        "Safety bounds [%s, %s] | %.1f%% clipped | ESS: %.0f treated, %.0f control",
        fmt_pct(bounds[1L], 1L), fmt_pct(bounds[2L], 1L),
        100 * t_res$g_info$pct_truncated,
        t_res$ess_treated, t_res$ess_control
      ),
      x = if (compact) "Propensity Score P(Treated | W)" else "Estimated Chance of Receiving Treatment P(A = 1 | W)",
      y = "Relative Density (Treated vs. Control)",
      caption = if (!compact) "Mirrored curves in the middle confirm real-world overlap; shaded red tails outside dashed lines are clipped." else NULL
    ) +
    theme_tidytmle(grid = "x", compact = compact) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )
  p
}

autoplot_learners <- function(object, compact = FALSE) {
  pal <- tmle_palette()
  w_tbl <- extract_sl_weights(object)
  m_tbl <- collect_metrics.tmle_fit(object, summarize = TRUE)

  factors_present <- intersect(c("Q", "g", "C"), unique(m_tbl$nuisance))
  blocks <- list()

  for (f in factors_present) {
    m_f <- m_tbl |>
      dplyr::filter(nuisance == f)

    chosen_metric <- if (f == "Q") {
      if ("rmse" %in% m_f$.metric) "rmse" else "mn_log_loss"
    } else {
      if ("mn_log_loss" %in% m_f$.metric) "mn_log_loss" else unique(m_f$.metric)[1L]
    }
    metric_lbl <- switch(chosen_metric, rmse = "RMSE", mn_log_loss = "Log-Loss", toupper(chosen_metric))
    row_lbl <- switch(
      f,
      Q = sprintf("Outcome (Q)\nCV %s", metric_lbl),
      g = sprintf("Treatment (g)\nCV %s", metric_lbl),
      C = sprintf("Missingness (C)\nCV %s", metric_lbl)
    )

    err_df <- m_f |>
      dplyr::filter(.metric == chosen_metric) |>
      dplyr::transmute(
        nuisance_order = match(f, c("Q", "g", "C")),
        row_panel = row_lbl,
        learner = learner,
        panel_col = "Cross-Validated Error (Lower = Better)",
        value = mean,
        val_lbl = sprintf("%.3f", mean),
        is_sl = learner == "SuperLearner"
      )

    w_f <- w_tbl |>
      dplyr::filter(nuisance == f) |>
      dplyr::transmute(
        nuisance_order = match(f, c("Q", "g", "C")),
        row_panel = row_lbl,
        learner = learner,
        panel_col = "Ensemble Weight",
        value = weight,
        val_lbl = sprintf("%.0f%%", 100 * weight),
        is_sl = FALSE
      )

    blocks[[length(blocks) + 1L]] <- dplyr::bind_rows(w_f, err_df)
  }

  plot_df <- dplyr::bind_rows(blocks)
  row_levels <- plot_df |>
    dplyr::distinct(nuisance_order, row_panel) |>
    dplyr::arrange(nuisance_order) |>
    dplyr::pull(row_panel)

  plot_df$row_panel <- factor(plot_df$row_panel, levels = row_levels)
  plot_df$panel_col <- factor(
    plot_df$panel_col,
    levels = c("Ensemble Weight", "Cross-Validated Error (Lower = Better)")
  )
  plot_df$role_fill <- ifelse(plot_df$is_sl, "Super Learner", "Base Algorithm")

  ggplot2::ggplot(plot_df, ggplot2::aes(x = stats::reorder(learner, value), y = value, fill = role_fill)) +
    ggplot2::geom_col(alpha = 0.88, width = 0.65, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = val_lbl),
      hjust = -0.15, size = if (compact) 2.7 else 3.1, color = pal$ink
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::facet_grid(row_panel ~ panel_col, scales = "free") +
    ggplot2::scale_fill_manual(values = c("Super Learner" = pal$crimson, "Base Algorithm" = pal$navy)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.24))) +
    ggplot2::labs(
      title = "Super Learner Ensemble Weights & Cross-Validated Error",
      subtitle = "Ensemble combines diverse models in each column; SuperLearner row highlighted in crimson",
      x = NULL,
      y = NULL
    ) +
    theme_tidytmle(grid = "x", compact = compact) +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0, vjust = 0.5),
      plot.margin = ggplot2::margin(12, 18, 10, 12)
    )
}

autoplot_eif <- function(object, compact = FALSE) {
  pal <- tmle_palette()
  eif <- object$targeting$eif
  lev <- object$targeting$eif_leverage
  df <- tibble::tibble(
    obs_id = seq_along(eif),
    eif = eif,
    leverage = lev,
    flag = factor(ifelse(lev > 0.05, "Elevated (> 5%)", "Normal"), levels = c("Normal", "Elevated (> 5%)"))
  )
  top_i <- which.max(lev)
  top_v <- lev[top_i]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = obs_id, y = leverage, color = flag)) +
    ggplot2::geom_segment(ggplot2::aes(xend = obs_id, y = 0, yend = leverage), alpha = 0.55, show.legend = FALSE) +
    ggplot2::geom_point(size = if (compact) 1.35 else 1.75, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0.10, linetype = "dashed", color = pal$red, linewidth = 0.7) +
    ggplot2::annotate(
      "text", x = 1, y = 0.10, label = "10% caution threshold",
      vjust = -0.45, hjust = 0, size = if (compact) 2.6 else 3.0, color = pal$red
    ) +
    ggplot2::scale_color_manual(values = c("Normal" = pal$navy, "Elevated (> 5%)" = pal$crimson)) +
    ggplot2::scale_y_continuous(
      labels = fmt_pct,
      expand = ggplot2::expansion(mult = c(0, 0.16))
    ) +
    ggplot2::labs(
      title = "Does Any Single Row Dominate Uncertainty?",
      subtitle = sprintf("Max single-row leverage = %.1f%% (row %d) | Caution threshold = 10%%", 100 * top_v, top_i),
      x = "Observation Index",
      y = if (compact) "Variance Share" else "Share of Total Uncertainty (EIF Variance)"
    ) +
    theme_tidytmle(grid = "y", compact = compact)

  if (top_v > 0.05) {
    p <- p + ggplot2::annotate(
      "text", x = top_i, y = top_v,
      label = sprintf(" Row %d (%.1f%%)", top_i, 100 * top_v),
      hjust = if (top_i > 0.8 * length(lev)) 1.1 else -0.1,
      vjust = 0.4, size = if (compact) 2.6 else 3.0,
      fontface = "bold", color = pal$crimson
    )
  }
  p
}

autoplot_sensitivity <- function(object, g_grid = seq(0.01, 0.15, by = 0.01), compact = FALSE) {
  pal <- tmle_palette()
  primary_est <- object$spec$estimand$names[1L]
  cur_low <- object$targeting$g_info$bounds[1L]
  is_binary_tx <- identical(object$cache$treatment_type, "binary")

  sweep_df <- purrr::map_dfr(g_grid, function(b) {
    rf <- retarget(object, g_bounds = c(b, 1 - b))
    td <- rf$targeting$estimates |>
      dplyr::filter(estimand == primary_est) |>
      dplyr::slice(1L)
    td$threshold <- b
    td$pct_truncated <- 100 * rf$targeting$g_info$pct_truncated
    td
  })

  any_clipped <- any(sweep_df$pct_truncated > 0.05)
  x_lab <- if (is_binary_tx) {
    if (compact) "Min Propensity Cutoff (b)" else "Positivity Cutoff b (clips propensity to [b, 1 - b])"
  } else {
    "Lower Bound Scalar b"
  }

  p <- ggplot2::ggplot(sweep_df, ggplot2::aes(x = threshold, y = estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high), fill = pal$navy, alpha = 0.14) +
    ggplot2::geom_line(color = pal$navy, linewidth = 0.95) +
    ggplot2::geom_point(size = if (compact) 2.0 else 2.5, color = pal$crimson) +
    ggplot2::scale_x_continuous(labels = function(z) fmt_pct(z, 0L)) +
    ggplot2::labs(
      title = sprintf("Does Strictness Change the %s?", sweep_df$term[1L]),
      subtitle = "Estimates & 95% CIs recomputed instantaneously via retarget() without refitting ML",
      x = x_lab,
      y = sprintf("Estimated %s", sweep_df$term[1L]),
      caption = if (!compact && is_binary_tx) sprintf("Dotted vertical line shows active lower bound (%.1f%%).", 100 * cur_low) else NULL
    ) +
    theme_tidytmle(grid = "y", compact = compact)

  if (is_binary_tx && is.finite(cur_low) && cur_low >= min(g_grid) && cur_low <= max(g_grid)) {
    p <- p +
      ggplot2::geom_vline(xintercept = cur_low, linetype = "dotted", color = pal$muted, linewidth = 0.7) +
      ggplot2::annotate(
        "text", x = cur_low, y = max(sweep_df$conf.high, na.rm = TRUE),
        label = sprintf(" Active (%.1f%%)", 100 * cur_low),
        hjust = 0, vjust = 1.2, size = if (compact) 2.6 else 3.0, color = pal$muted
      )
  }

  if (any_clipped) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f%%", pct_truncated)),
      vjust = -1.1, size = if (compact) 2.3 else 2.7, color = pal$muted
    ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.14)))
  }
  p
}

autoplot_balance <- function(object, compact = FALSE) {
  if (!identical(object$cache$treatment_type, "binary")) {
    cli::cli_abort(c(
      "Covariate balance dumbbell plot is only available for binary treatment estimands.",
      "i" = "For continuous Modified Treatment Policies, inspect {.code autoplot(fit, type = \"overlap\")} or {.code type = \"eif\"}."
    ))
  }
  pal <- tmle_palette()
  bal_wide <- compute_covariate_balance(object)
  if (nrow(bal_wide) == 0L) {
    cli::cli_abort(c(
      "No varying baseline covariates remain in the analyzed sample or subgroup to plot balance.",
      "i" = "Reset the subgroup with {.code retarget(fit, subgroup = TRUE)} or inspect {.code autoplot(fit, type = \"overlap\")}."
    ))
  }
  max_raw <- max(bal_wide$smd_raw, na.rm = TRUE)
  max_wt <- max(bal_wide$smd_wt, na.rm = TRUE)
  bal_sub <- if (max_wt < max_raw) {
    sprintf("Max standardized difference dropped from %.3f (raw) to %.3f (weighted)", max_raw, max_wt)
  } else {
    sprintf("Max standardized difference was %.3f (raw) vs. %.3f (weighted)", max_raw, max_wt)
  }

  bal_wide$covariate <- stats::reorder(bal_wide$covariate, bal_wide$smd_raw)
  bal_long <- dplyr::bind_rows(
    tibble::tibble(covariate = bal_wide$covariate, group = "Before weighting", smd = bal_wide$smd_raw),
    tibble::tibble(covariate = bal_wide$covariate, group = "After weighting", smd = bal_wide$smd_wt)
  )
  bal_long$group <- factor(bal_long$group, levels = c("Before weighting", "After weighting"))

  ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = 0, xmax = 0.10, ymin = -Inf, ymax = Inf, fill = pal$teal_bg, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0.10, linetype = "dashed", color = pal$teal, linewidth = 0.7) +
    ggplot2::annotate(
      "text", x = 0.10, y = Inf, label = " Balanced (< 0.10)",
      vjust = 1.4, hjust = 0, size = if (compact) 2.6 else 3.0, color = pal$teal, fontface = "bold"
    ) +
    ggplot2::geom_segment(
      data = bal_wide,
      ggplot2::aes(x = smd_raw, xend = smd_wt, y = covariate, yend = covariate),
      color = pal$axis, linewidth = 1.0
    ) +
    ggplot2::geom_point(
      data = bal_long,
      ggplot2::aes(x = smd, y = covariate, color = group, shape = group),
      size = if (compact) 2.4 else 3.0
    ) +
    ggplot2::scale_color_manual(values = c("Before weighting" = pal$amber, "After weighting" = pal$navy)) +
    ggplot2::scale_shape_manual(values = c("Before weighting" = 17L, "After weighting" = 16L)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.14))) +
    ggplot2::labs(
      title = "Did Propensity Weighting Balance Confounders?",
      subtitle = bal_sub,
      x = if (compact) "Imbalance |SMD|" else "Absolute Standardized Mean Difference |SMD|",
      y = NULL,
      color = NULL,
      shape = NULL,
      caption = if (!compact) "Points inside the green shaded band (< 0.10) indicate well-balanced confounders across arms." else NULL
    ) +
    theme_tidytmle(grid = "x", compact = compact) +
    ggplot2::theme(legend.position = "top")
}
