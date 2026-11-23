#' Estimate Heterogeneous Treatment Effects (CATE) via DR-Learner
#'
#' Feeds the doubly-robust Efficient Influence Function pseudo-outcomes from a
#' fitted `"tmle_fit"` object into a second-stage `tidymodels` regression
#' workflow or Super Learner (Kennedy, 2023) to estimate Conditional Average
#' Treatment Effects \eqn{\tau(V) = E[Y(1) - Y(0) \mid V]}.
#'
#' @details
#' Because the individual treatment effect \eqn{Y_i(1) - Y_i(0)} is never
#' directly observed, Kennedy's Doubly Robust (DR) Learner constructs the
#' cross-fitted pseudo-outcome
#' \deqn{\varphi_i = H_i \left(Y_i - \bar{Q}^*(A_i, W_i)\right) + \left(\bar{Q}^*(1, W_i) - \bar{Q}^*(0, W_i)\right)}
#' whose conditional expectation given covariates \eqn{V} equals the true CATE
#' \eqn{\tau(V)} whenever *either* \eqn{\bar{Q}} or \eqn{g} is consistently
#' estimated. Regressing \eqn{\varphi_i} on \eqn{V_i} lets any `tidymodels`
#' regression workflow learn treatment effect heterogeneity while `tidy()`
#' reports Chernozhukov et al.'s (2018) Best Linear Projection (BLP) test for
#' significant effect modification.
#'
#' @references
#' Kennedy, E. H. (2023). Towards optimal doubly robust estimation of
#' heterogeneous causal effects. *Electronic Journal of Statistics*, 17(2),
#' 3008–3049.
#'
#' Chernozhukov, V., Demirer, M., Duflo, E., & Fernandez-Val, I. (2018).
#' Generic machine learning inference on heterogeneous treatment effects in
#' randomized experiments. *NBER Working Paper No. 24678*.
#'
#' @param fit A fitted `"tmle_fit"` object with binary treatment.
#' @param x A fitted `"tmle_cate"` object passed to `tidy()` or `augment()`.
#' @param ... Additional arguments (currently unused).
#' @param learners `"default"`, a `parsnip` `model_spec` (regression mode), a
#'   `workflow`, a `workflow_set`, or a list of specs/workflows.
#' @param covariates Optional character vector or tidyselect expression of
#'   effect modifiers \eqn{V \subseteq W}. Defaults to all baseline covariates
#'   from `fit`.
#' @param v_folds Number of cross-fitting folds for out-of-fold CATE evaluation
#'   and Best Linear Projection heterogeneity testing (default `5L`).
#'
#' @return An S3 object of class `"tmle_cate"` supporting `predict()`,
#'   `tidy()`, `augment()`, `autoplot()`, and [qini_auc()].
#'
#' @examples
#' df <- sim_causal_cohort(n = 200, outcome_type = "continuous", seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' cate <- fit_cate(fit, v_folds = 3L)
#' tidy(cate)
#' qini_auc(cate)
#'
#' @export
fit_cate <- function(fit,
                     learners = "default",
                     covariates = NULL,
                     v_folds = 5L) {
  if (!inherits(fit, "tmle_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls tmle_fit} object.")
  }
  if (fit$cache$treatment_type != "binary") {
    cli::cli_abort("DR-Learner CATE estimation requires a binary treatment {.cls tmle_fit}.")
  }

  idx <- fit$subgroup_idx
  aug <- augment.tmle_fit(fit)
  t_res <- fit$targeting

  # Doubly robust pseudo-outcome:
  # phi_i = H_i * (Y_i - Q*_i(A_i, W_i)) + (Q*_i(1, W_i) - Q*_i(0, W_i))
  y_obs <- fit$cache$y_num[idx]
  y_safe <- ifelse(fit$cache$delta[idx] == 1L, y_obs, 0)
  c_obs <- ifelse(fit$cache$a_num[idx] == 1L, t_res$c1_hat, t_res$c0_hat)
  c_obs <- pmax(c_obs, 0.025)
  h_dr <- ifelse(fit$cache$delta[idx] == 1L, t_res$clever_h, 0)
  phi <- h_dr * (y_safe - t_res$qa_star) + (t_res$q1_star - t_res$q0_star)

  cov_quo <- rlang::enquo(covariates)
  v_cols <- if (rlang::quo_is_null(cov_quo)) {
    fit$cache$covariate_cols
  } else {
    names(tidyselect::eval_select(cov_quo, data = aug))
  }
  if (length(v_cols) == 0L) {
    cli::cli_abort("At least 1 effect-modifier covariate must be selected for CATE estimation.")
  }

  cate_df <- aug[, v_cols, drop = FALSE]
  cate_df[[".phi"]] <- phi

  fml <- stats::reformulate(termlabels = v_cols, response = ".phi")
  rec_cate <- recipes::recipe(fml, data = cate_df) |>
    recipes::step_impute_median(recipes::all_numeric_predictors()) |>
    recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
    recipes::step_novel(recipes::all_nominal_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE) |>
    recipes::step_zv(recipes::all_predictors())

  spec_list <- if (identical(learners, "default") || is.null(learners)) {
    default_nuisance_learners(mode = "regression")
  } else if (inherits(learners, "model_spec")) {
    stats::setNames(list(parsnip::set_mode(learners, "regression")), class(learners)[1L])
  } else if (inherits(learners, "workflow")) {
    list(custom_cate = learners)
  } else if (inherits(learners, "workflow_set")) {
    ids <- learners$wflow_id
    stats::setNames(lapply(ids, function(id) workflowsets::extract_workflow(learners, id = id)), ids)
  } else if (is.list(learners)) {
    nms <- names(learners)
    if (is.null(nms) || any(nms == "")) nms <- paste0("cate_", seq_along(learners))
    stats::setNames(learners, nms)
  } else {
    cli::cli_abort("Unsupported {.arg learners} input for {.fn fit_cate}.")
  }

  wfs <- lapply(spec_list, function(item) {
    if (inherits(item, "model_spec")) {
      workflows::workflow() |>
        workflows::add_model(parsnip::set_mode(item, "regression")) |>
        workflows::add_recipe(rec_cate)
    } else {
      ms <- parsnip::set_mode(workflows::extract_spec_parsnip(item), "regression")
      workflows::update_model(item, ms) |>
        workflows::update_recipe(rec_cate)
    }
  })

  splits <- rsample::vfold_cv(cate_df, v = min(as.integer(v_folds), nrow(cate_df)))$splits
  sl_cate <- fit_cv_super_learner(
    workflows_list = wfs,
    splits = splits,
    data = cate_df,
    y_target = phi,
    mode = "regression",
    eval_datasets = list(),
    weights = fit$cache$weights[idx],
    metalearner = "nnls",
    cv_tmle = TRUE,
    save_nuisance = "full",
    parallel = FALSE,
    nuisance_name = "CATE"
  )

  # Also fit full-sample workflows for new-data prediction
  full_models <- lapply(wfs, function(wf) fit_single_workflow_safe(wf, cate_df))

  # Best Linear Projection heterogeneity test (Chernozhukov et al. 2018):
  # Regress phi on centered OOF CATE predictions
  tau_oof <- sl_cate$sl_preds
  tau_centered <- tau_oof - mean(tau_oof)
  blp_fit <- stats::lm(phi ~ tau_centered, weights = fit$cache$weights[idx])
  blp_tidy <- broom::tidy(blp_fit, conf.int = TRUE)
  blp_tidy$term <- c("ATE (BLP Intercept)", "Heterogeneity Loading (BLP Slope)")

  structure(
    list(
      data = aug,
      v_cols = v_cols,
      phi = phi,
      tau_oof = tau_oof,
      sl_weights = sl_cate$sl_weights,
      metrics = sl_cate$metrics,
      full_models = full_models,
      blp = blp_tidy
    ),
    class = "tmle_cate"
  )
}

#' Predict Conditional Average Treatment Effects on New Data
#'
#' @param object A fitted `"tmle_cate"` object.
#' @param new_data A `data.frame` containing the CATE effect modifiers.
#' @param ... Currently unused.
#'
#' @return A 1-column [tibble::tibble()] with column `.pred`.
#' @method predict tmle_cate
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, outcome_type = "continuous", seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' cate <- fit_cate(fit, v_folds = 3L)
#' predict(cate, new_data = head(df, 5L))
#'
#' @export
predict.tmle_cate <- function(object, new_data, ...) {
  w <- object$sl_weights
  active <- which(w > 0)
  pred_mat <- sapply(active, function(k) {
    predict_workflow_safe(object$full_models[[k]], new_df = new_data, mode = "regression")
  })
  if (!is.matrix(pred_mat)) {
    pred_mat <- matrix(pred_mat, nrow = nrow(new_data), ncol = length(active))
  }
  preds <- as.numeric(pred_mat %*% w[active])
  tibble::tibble(.pred = preds)
}

#' @rdname fit_cate
#' @method tidy tmle_cate
#' @export
tidy.tmle_cate <- function(x, ...) {
  x$blp
}

#' @rdname fit_cate
#' @method augment tmle_cate
#' @export
augment.tmle_cate <- function(x, ...) {
  dplyr::mutate(x$data, .cate = x$tau_oof, .eif_pseudo = x$phi)
}

#' @export
print.tmle_cate <- function(x, ...) {
  lines <- cli::cli_format_method({
    cli::cli_h2("Doubly Robust CATE (DR-Learner) Fit")
    cli::cli_text("Modifiers: {.field {x$v_cols}}")
    cli::cli_text("Out-of-fold CATE range: [{sprintf('%.3f', min(x$tau_oof))}, {sprintf('%.3f', max(x$tau_oof))}] (SD = {sprintf('%.3f', stats::sd(x$tau_oof))})")
  })
  cat(lines, sep = "\n")
  print(x$blp)
  invisible(x)
}

#' Compute Qini Area Under the Uplift Curve (AUUC / Qini)
#'
#' Evaluates treatment-prioritization performance of predicted CATE scores
#' against doubly-robust pseudo-outcomes (or observed outcomes).
#'
#' @param cate_fit A `"tmle_cate"` object (or numeric vector of CATE predictions).
#' @param pseudo_outcome Optional numeric vector of pseudo-outcomes (defaults
#'   to `cate_fit$phi` when `cate_fit` is a `"tmle_cate"` object).
#'
#' @return A 1-row [tibble::tibble()] with columns `.metric`, `.estimator`, and
#'   `.estimate` (`qini_auc` and `auuc`).
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, outcome_type = "continuous", seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' cate <- fit_cate(fit, v_folds = 3L)
#' qini_auc(cate)
#'
#' @export
qini_auc <- function(cate_fit, pseudo_outcome = NULL) {
  if (inherits(cate_fit, "tmle_cate")) {
    score <- cate_fit$tau_oof
    phi <- if (is.null(pseudo_outcome)) cate_fit$phi else pseudo_outcome
  } else {
    score <- as.numeric(cate_fit)
    phi <- as.numeric(pseudo_outcome)
  }
  n <- length(score)
  ord <- order(score, decreasing = TRUE)
  phi_ord <- phi[ord]
  frac <- seq_len(n) / n
  cum_gain <- cumsum(phi_ord) / n
  ate_overall <- mean(phi_ord)
  random_gain <- frac * ate_overall

  # Trapezoidal integration of cumulative uplift above random targeting
  df_q <- cum_gain - random_gain
  qini_val <- sum(0.5 * (df_q[-1L] + df_q[-n]) * diff(frac))
  auuc_val <- sum(0.5 * (cum_gain[-1L] + cum_gain[-n]) * diff(frac))

  tibble::tibble(
    .metric = c("qini_auc", "auuc"),
    .estimator = "standard",
    .estimate = c(qini_val, auuc_val)
  )
}

#' Diagnostic Plots for a Doubly Robust CATE (DR-Learner) Fit
#'
#' Visualizes treatment-effect heterogeneity and prioritization uplift for a
#' `"tmle_cate"` object returned by [fit_cate()]:
#' - `"qini"`: Plots the **Qini Cumulative Uplift Curve** ordering individuals
#'   from highest to lowest predicted benefit \eqn{\hat{\tau}(V)} against the
#'   random-targeting diagonal baseline, shading the incremental uplift gain
#'   (Qini AUC).
#' - `"cate"`: Plots the distribution of out-of-fold predicted CATEs
#'   \eqn{\hat{\tau}(V)} relative to the population Average Treatment Effect,
#'   annotating the Chernozhukov et al. (2018) heterogeneity p-value.
#'
#' @param object A fitted `"tmle_cate"` object.
#' @param type Plot type: `"qini"` (default) or `"cate"`.
#' @param ... Currently unused.
#'
#' @return A `ggplot2` object.
#' @method autoplot tmle_cate
#'
#' @examples
#' df <- sim_causal_cohort(n = 200, outcome_type = "continuous", seed = 42)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' cate <- fit_cate(fit, v_folds = 3L)
#' autoplot(cate, type = "qini")
#' autoplot(cate, type = "cate")
#'
#' @export
autoplot.tmle_cate <- function(object, type = c("qini", "cate"), ...) {
  type <- match.arg(type)
  pal <- tmle_palette()

  if (type == "qini") {
    score <- object$tau_oof
    phi <- object$phi
    n <- length(score)
    ord <- order(score, decreasing = TRUE)
    phi_ord <- phi[ord]
    frac <- c(0, seq_len(n) / n)
    cum_gain <- c(0, cumsum(phi_ord) / n)
    random_gain <- frac * mean(phi_ord)
    q_stats <- qini_auc(object)

    mid_i <- ceiling(0.5 * length(frac))
    tot_ben <- cum_gain[length(cum_gain)]
    half_share <- if (is.finite(tot_ben) && abs(tot_ben) > 1e-8) {
      100 * cum_gain[mid_i] / tot_ben
    } else {
      NA_real_
    }
    sub_qini <- if (is.finite(half_share)) {
      sprintf("Top 50%% prioritized captures %.0f%% of total benefit | Qini AUC = %.4f",
              half_share, q_stats$.estimate[1L])
    } else {
      sprintf("Qini AUC = %.4f | AUUC = %.4f (higher curve = stronger targeting)",
              q_stats$.estimate[1L], q_stats$.estimate[2L])
    }

    df_curve <- tibble::tibble(
      frac = frac,
      cum_gain = cum_gain,
      random_gain = random_gain
    )
    df_lines <- dplyr::bind_rows(
      tibble::tibble(frac = frac, gain = cum_gain, curve = "DR-Learner"),
      tibble::tibble(frac = frac, gain = random_gain, curve = "Random")
    )

    q_lbl_i <- ceiling(0.38 * length(frac))
    q_mid_x <- frac[q_lbl_i]
    q_mid_y <- cum_gain[q_lbl_i]
    r_mid_y <- random_gain[mid_i]

    return(
      ggplot2::ggplot() +
        ggplot2::geom_ribbon(
          data = df_curve,
          ggplot2::aes(x = frac, ymin = random_gain, ymax = cum_gain),
          fill = pal$crimson, alpha = 0.14
        ) +
        ggplot2::geom_line(
          data = df_lines,
          ggplot2::aes(x = frac, y = gain, color = curve, linetype = curve),
          linewidth = 0.95, show.legend = FALSE
        ) +
        ggplot2::annotate(
          "text", x = q_mid_x, y = q_mid_y,
          label = "Targeted by DR-Learner",
          vjust = -1.25, hjust = 0.72, fontface = "bold",
          size = 3.2, color = pal$crimson
        ) +
        ggplot2::annotate(
          "text", x = 0.60, y = r_mid_y,
          label = "Random order",
          vjust = 1.7, hjust = 0,
          size = 3.1, color = pal$muted
        ) +
        ggplot2::scale_color_manual(values = c("DR-Learner" = pal$crimson, "Random" = pal$muted)) +
        ggplot2::scale_linetype_manual(values = c("DR-Learner" = "solid", "Random" = "dashed")) +
        ggplot2::scale_x_continuous(labels = fmt_pct) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.14))) +
        ggplot2::labs(
          title = "Does Targeting by Predicted Benefit Outperform Random Care?",
          subtitle = sub_qini,
          x = "Share of Population Targeted (Highest to Lowest Predicted Benefit)",
          y = "Cumulative Incremental Outcome Gain",
          caption = "Shaded pink area above dashed line measures surplus benefit gained from personalized targeting."
        ) +
        theme_tidytmle(grid = "y") +
        ggplot2::coord_cartesian(clip = "off")
    )
  }

  ate_val <- mean(object$tau_oof)
  pct_pos <- mean(object$tau_oof > 0)
  het_p <- object$blp$p.value[2L]
  p_txt <- format_p_human(het_p)
  df_c <- tibble::tibble(
    cate_val = object$tau_oof,
    benefit_group = ifelse(object$tau_oof >= 0, "Predicted Benefit (>= 0)", "No / Negative Benefit (< 0)")
  )

  ggplot2::ggplot(df_c, ggplot2::aes(x = cate_val, fill = benefit_group)) +
    ggplot2::geom_histogram(bins = 28L, alpha = 0.84, color = "white", show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = 0, color = pal$ink, linewidth = 0.55) +
    ggplot2::geom_vline(xintercept = ate_val, linetype = "dashed", color = pal$navy, linewidth = 0.85) +
    ggplot2::annotate(
      "text", x = ate_val, y = Inf,
      label = sprintf(" Avg = %+.2f", ate_val),
      vjust = 1.15, hjust = 0, fontface = "bold",
      size = 3.2, color = pal$navy
    ) +
    ggplot2::scale_fill_manual(
      values = c("Predicted Benefit (>= 0)" = pal$crimson, "No / Negative Benefit (< 0)" = pal$amber)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.14))) +
    ggplot2::labs(
      title = "Who Benefits Most? Individualized Treatment Effects",
      subtitle = sprintf(
        "%.0f%% predicted to benefit (> 0) | SD = %.2f | Heterogeneity test: %s",
        100 * pct_pos, stats::sd(object$tau_oof), p_txt
      ),
      x = "Predicted Individual Effect tau(V)",
      y = "Number of People",
      caption = "Solid black line marks 0 (no effect); dashed navy line marks overall average effect."
    ) +
    theme_tidytmle(grid = "y")
}

#' Simulate a Synthetic Causal Cohort with Exact Ground-Truth Parameters
#'
#' Generates a realistic biomedical cohort featuring nonlinear confounding,
#' heterogeneous treatment effects, optional informative right-censoring,
#' optional cluster correlation, and optional positivity stress. Attaches
#' ground-truth causal estimands as `"truth"` attribute metadata for simulation
#' benchmarking and verification.
#'
#' @param n Number of individuals (default `500L`).
#' @param outcome_type `"continuous"` (default) or `"binary"`.
#' @param treatment_type `"binary"` (default) or `"continuous"`.
#' @param censoring_rate Approximate proportion of missing/censored outcomes
#'   in `[0, 0.5)` (default `0`, i.e., complete outcomes).
#' @param n_clusters Optional integer number of clinic clusters (default `NULL`
#'   for i.i.d. data).
#' @param positivity_stress Logical; if `TRUE`, creates strong propensity
#'   overlap violations for testing adaptive positivity bounds and bounded
#'   logistic fluctuation.
#' @param seed Optional integer random seed.
#'
#' @return A [tibble::tibble()] with columns `age`, `sex`, `lab_1`, `lab_2`,
#'   `lab_3`, `a`, `y`, and optionally `delta` and `clinic_id`, plus an attribute
#'   `"truth"` containing `true_ate`, `true_att`, `true_atc`, `true_rr`,
#'   `true_or`, `true_shift_1`, and `true_cate`.
#'
#' @examples
#' cohort <- sim_causal_cohort(n = 150, outcome_type = "binary", seed = 42)
#' attr(cohort, "truth")$true_ate
#'
#' @export
sim_causal_cohort <- function(n = 500L,
                              outcome_type = c("continuous", "binary"),
                              treatment_type = c("binary", "continuous"),
                              censoring_rate = 0,
                              n_clusters = NULL,
                              positivity_stress = FALSE,
                              seed = NULL) {
  outcome_type <- match.arg(outcome_type)
  treatment_type <- match.arg(treatment_type)

  if (!is.null(seed)) {
    has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (has_seed) {
      old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit({
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
          rm(".Random.seed", envir = globalenv())
        }
      }, add = TRUE)
    }
    set.seed(seed)
  }

  age <- stats::rnorm(n, mean = 55, sd = 12)
  sex <- factor(sample(c("F", "M"), size = n, replace = TRUE, prob = c(0.52, 0.48)))
  sex_num <- as.numeric(sex == "M")
  lab_1 <- stats::rnorm(n, mean = 0, sd = 1)
  lab_2 <- stats::rnorm(n, mean = 0.3 * lab_1, sd = 1)
  lab_3 <- stats::rbinom(n, size = 1L, prob = 0.4)

  z_age <- (age - 55) / 12
  cl_re <- rep(0, n)
  clinic_id <- NULL
  if (!is.null(n_clusters)) {
    clinic_id <- factor(sprintf("clinic_%02d", sample(seq_len(n_clusters), n, replace = TRUE)))
    cl_effects <- stats::rnorm(n_clusters, mean = 0, sd = 0.35)
    cl_re <- cl_effects[as.integer(clinic_id)]
  }

  # Treatment mechanism
  if (treatment_type == "binary") {
    scale_ps <- if (isTRUE(positivity_stress)) 3.2 else 0.75
    eta_g <- scale_ps * (0.4 * z_age - 0.5 * lab_1 + 0.35 * lab_2^2 - 0.3 * sex_num + 0.25 * lab_3 - 0.15)
    p_a <- stats::plogis(eta_g)
    a <- stats::rbinom(n, size = 1L, prob = p_a)
  } else {
    mu_a <- 1.5 + 0.6 * z_age - 0.5 * lab_1 + 0.3 * sex_num + 0.2 * lab_2
    a <- mu_a + stats::rnorm(n, mean = 0, sd = 0.9)
  }

  # Baseline outcome potential and heterogeneous treatment effect
  mu_0 <- 2.0 + 0.8 * z_age + 0.9 * lab_1 + 0.4 * (lab_2^2) - 0.5 * sex_num + 0.6 * lab_3 + cl_re
  tau_i <- 1.4 + 0.7 * z_age + 0.6 * (lab_1 > 0) - 0.4 * sex_num

  if (outcome_type == "continuous") {
    eps_y <- stats::rnorm(n, mean = 0, sd = 1.0)
    if (treatment_type == "binary") {
      y0_pot <- mu_0
      y1_pot <- mu_0 + tau_i
      y <- ifelse(a == 1L, y1_pot, y0_pot) + eps_y
      truth <- list(
        true_ey1 = mean(y1_pot),
        true_ey0 = mean(y0_pot),
        true_ate = mean(tau_i),
        true_att = mean(tau_i[a == 1L]),
        true_atc = mean(tau_i[a == 0L]),
        true_cate = tau_i
      )
    } else {
      y_obs_mean <- mu_0 + 0.8 * a + 0.25 * a * z_age
      y_shift_mean <- mu_0 + 0.8 * (a + 1) + 0.25 * (a + 1) * z_age
      y <- y_obs_mean + eps_y
      truth <- list(
        true_shift_1 = mean(y_shift_mean),
        true_scale_09 = mean(mu_0 + 0.8 * (0.9 * a) + 0.25 * (0.9 * a) * z_age)
      )
    }
  } else {
    eta_0 <- -0.8 + 0.5 * z_age + 0.6 * lab_1 + 0.3 * lab_2 - 0.3 * sex_num + 0.4 * cl_re
    tau_logit <- 0.85 + 0.35 * z_age
    if (treatment_type == "binary") {
      p0_pot <- stats::plogis(eta_0)
      p1_pot <- stats::plogis(eta_0 + tau_logit)
      p_obs <- ifelse(a == 1L, p1_pot, p0_pot)
      y <- stats::rbinom(n, size = 1L, prob = p_obs)
      ey1 <- mean(p1_pot)
      ey0 <- mean(p0_pot)
      truth <- list(
        true_ey1 = ey1,
        true_ey0 = ey0,
        true_ate = ey1 - ey0,
        true_att = mean((p1_pot - p0_pot)[a == 1L]),
        true_atc = mean((p1_pot - p0_pot)[a == 0L]),
        true_rr = ey1 / ey0,
        true_or = (ey1 / (1 - ey1)) / (ey0 / (1 - ey0)),
        true_cate = p1_pot - p0_pot
      )
    } else {
      p_obs <- stats::plogis(eta_0 + 0.45 * a)
      p_shift <- stats::plogis(eta_0 + 0.45 * (a + 1))
      y <- stats::rbinom(n, size = 1L, prob = p_obs)
      truth <- list(
        true_shift_1 = mean(p_shift)
      )
    }
  }

  out <- tibble::tibble(
    age = age,
    sex = sex,
    lab_1 = lab_1,
    lab_2 = lab_2,
    lab_3 = lab_3,
    a = a,
    y = y
  )

  if (!is.null(clinic_id)) {
    out$clinic_id <- clinic_id
  }

  if (censoring_rate > 0) {
    base_c <- stats::qlogis(1 - censoring_rate)
    p_obs_c <- stats::plogis(base_c - 0.5 * (if (treatment_type == "binary") a else scale(a)[, 1L]) - 0.4 * z_age + 0.3 * lab_1)
    delta <- stats::rbinom(n, size = 1L, prob = p_obs_c)
    out$y[delta == 0L] <- NA_real_
    out$delta <- delta
  }

  attr(out, "truth") <- truth
  out
}
