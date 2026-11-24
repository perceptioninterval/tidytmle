#' Fit a TMLE Specification on Data
#'
#' Executes cross-fitted nuisance Super Learning for outcome \eqn{Q}, treatment
#' mechanism \eqn{g} (or density ratio for Modified Treatment Policies), and
#' optional censoring/missingness mechanism \eqn{C}, caches out-of-fold nuisance
#' predictions, and executes the targeted fluctuation and EIF inference step.
#'
#' @param object A [tmle_spec()] object.
#' @param data A `data.frame` or `tibble` containing the cohort observations.
#' @param ... Currently unused; included for S3 consistency with `generics::fit`.
#'
#' @return An S3 object of class `"tmle_fit"`.
#' @method fit tmle_spec
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, outcome_type = "binary", seed = 1)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = c("ate", "rr")) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(data = df)
#' tidy(fit)
#'
#' @export
fit.tmle_spec <- function(object, data, ...) {
  t_start <- proc.time()[["elapsed"]]
  check_is_tmle_spec(object)

  # 1. Resolve causal roles and handle automatic informative missingness
  roles <- resolve_causal_roles(object, data)

  # 2. Construct default rsample folds if user didn't supply set_resamples()
  resamples <- if (is.null(object$resamples$obj)) {
    build_default_resamples(
      roles = roles,
      v_folds = object$resamples$v_folds,
      strata = object$resamples$strata
    )
  } else {
    object$resamples$obj
  }
  splits <- resamples$splits

  # 3. Build nuisance data and role-safe recipes (zero leakage)
  nuisance_data <- build_nuisance_data(roles)
  recipes_split <- split_causal_recipes(object$preprocessing, roles, nuisance_data)
  has_custom_rec <- !is.null(object$preprocessing)

  q_mode <- if (roles$outcome_type == "binary") "classification" else "regression"
  is_q_strat <- isTRUE(object$learners$Q$stratify) && identical(roles$treatment_type, "binary")
  q_role <- if (is_q_strat) "Q_strat" else "Q"
  q_default_rec <- if (is_q_strat) {
    adapt_recipe_for_role(object$preprocessing, role_type = "Q_strat", roles = roles, template_data = nuisance_data$Q)
  } else {
    recipes_split$Q
  }

  q_wfs <- normalize_learner_library(
    learners_cfg = object$learners$Q,
    mode = q_mode,
    role_type = q_role,
    roles = roles,
    template_data = nuisance_data$Q,
    default_recipe = q_default_rec,
    has_custom_recipe = has_custom_rec
  )

  obs_idx <- which(roles$delta == 1L)

  # 4. Cross-fitted Super Learning for Q (pooled S-learner or arm-stratified T-learner)
  if (is_q_strat) {
    obs_idx_1 <- intersect(obs_idx, which(roles$a_num == 1L))
    obs_idx_0 <- intersect(obs_idx, which(roles$a_num == 0L))
    sl_Q1 <- fit_cv_super_learner(
      workflows_list = q_wfs,
      splits = splits,
      data = nuisance_data$Q,
      y_target = roles$y_num,
      mode = q_mode,
      eval_datasets = list(all = nuisance_data$Q),
      train_subset_idx = obs_idx_1,
      weights = roles$weights,
      metalearner = object$learners$Q$metalearner,
      cv_tmle = object$targeting$cv_tmle,
      save_nuisance = object$targeting$save_nuisance,
      parallel = object$targeting$parallel,
      nuisance_name = "Q"
    )
    sl_Q0 <- fit_cv_super_learner(
      workflows_list = q_wfs,
      splits = splits,
      data = nuisance_data$Q,
      y_target = roles$y_num,
      mode = q_mode,
      eval_datasets = list(all = nuisance_data$Q),
      train_subset_idx = obs_idx_0,
      weights = roles$weights,
      metalearner = object$learners$Q$metalearner,
      cv_tmle = object$targeting$cv_tmle,
      save_nuisance = object$targeting$save_nuisance,
      parallel = object$targeting$parallel,
      nuisance_name = "Q"
    )
    q1_s <- sl_Q1$sl_eval_preds$all
    q0_s <- sl_Q0$sl_eval_preds$all
    qa_s <- ifelse(roles$a_num == 1L, q1_s, q0_s)
    oof_comb <- sl_Q0$oof_matrix
    i1_rows <- which(roles$a_num == 1L)
    oof_comb[i1_rows, ] <- sl_Q1$oof_matrix[i1_rows, , drop = FALSE]
    wt_comb <- 0.5 * (sl_Q1$sl_weights + sl_Q0$sl_weights)
    wt_comb <- wt_comb / sum(wt_comb)
    fw_comb <- dplyr::bind_rows(sl_Q1$fold_weights, sl_Q0$fold_weights) |>
      dplyr::group_by(nuisance, fold, learner) |>
      dplyr::summarise(weight = mean(weight), .groups = "drop")
    mets_comb <- compute_nuisance_metrics(
      oof_obs = oof_comb,
      sl_obs = qa_s,
      y_target = roles$y_num,
      fit_idx = obs_idx,
      fold_id_vec = sl_Q1$fold_id,
      mode = q_mode,
      nuisance_name = "Q"
    )
    sl_Q <- list(
      nuisance = "Q",
      mode = q_mode,
      sl_preds = qa_s,
      sl_eval_preds = list(q1 = q1_s, q0 = q0_s),
      oof_matrix = oof_comb,
      sl_weights = wt_comb,
      fold_weights = fw_comb,
      metrics = mets_comb,
      fold_id = sl_Q1$fold_id,
      models = c(sl_Q1$models, sl_Q0$models)
    )
  } else {
    # Build counterfactual evaluation datasets for Q (numeric treatment_col)
    if (roles$treatment_type == "binary") {
      df_a1 <- nuisance_data$Q
      df_a1[[roles$treatment_col]] <- 1
      df_a0 <- nuisance_data$Q
      df_a0[[roles$treatment_col]] <- 0
      eval_Q <- list(q1 = df_a1, q0 = df_a0)
    } else {
      policy <- object$estimand$policies[[1L]]
      df_shift <- nuisance_data$Q
      df_shift[[roles$treatment_col]] <- apply_policy(policy, roles$a_num)
      eval_Q <- list(qa = nuisance_data$Q, q_shift = df_shift)
    }

    sl_Q <- fit_cv_super_learner(
      workflows_list = q_wfs,
      splits = splits,
      data = nuisance_data$Q,
      y_target = roles$y_num,
      mode = q_mode,
      eval_datasets = eval_Q,
      train_subset_idx = obs_idx,
      weights = roles$weights,
      metalearner = object$learners$Q$metalearner,
      cv_tmle = object$targeting$cv_tmle,
      save_nuisance = object$targeting$save_nuisance,
      parallel = object$targeting$parallel,
      nuisance_name = "Q"
    )
  }

  # 5. Cross-fitted Super Learning for g (propensity or stacked density-ratio)
  if (roles$treatment_type == "binary") {
    g_wfs <- normalize_learner_library(
      learners_cfg = object$learners$g,
      mode = "classification",
      role_type = "g",
      roles = roles,
      template_data = nuisance_data$g,
      default_recipe = recipes_split$g,
      has_custom_recipe = has_custom_rec
    )
    sl_g <- fit_cv_super_learner(
      workflows_list = g_wfs,
      splits = splits,
      data = nuisance_data$g,
      y_target = roles$a_num,
      mode = "classification",
      eval_datasets = list(),
      train_subset_idx = NULL,
      weights = roles$weights,
      metalearner = object$learners$g$metalearner,
      cv_tmle = object$targeting$cv_tmle,
      save_nuisance = object$targeting$save_nuisance,
      parallel = object$targeting$parallel,
      nuisance_name = "g"
    )
  } else {
    policy <- object$estimand$policies[[1L]]
    dr_template <- nuisance_data$g
    dr_template[[".lambda"]] <- factor(0L, levels = c("0", "1"))
    rec_g_dr <- adapt_recipe_for_role(
      object$preprocessing,
      role_type = "g_dr",
      roles = roles,
      template_data = dr_template
    )
    g_wfs <- normalize_learner_library(
      learners_cfg = object$learners$g,
      mode = "classification",
      role_type = "g_dr",
      roles = roles,
      template_data = dr_template,
      default_recipe = rec_g_dr,
      has_custom_recipe = has_custom_rec
    )
    sl_g <- fit_cv_density_ratio_sl(
      workflows_list = g_wfs,
      splits = splits,
      data = nuisance_data$g,
      treatment_col = roles$treatment_col,
      policy = policy,
      weights = roles$weights,
      metalearner = object$learners$g$metalearner,
      cv_tmle = object$targeting$cv_tmle,
      save_nuisance = object$targeting$save_nuisance,
      parallel = object$targeting$parallel
    )
  }

  # 6. Optional Censoring / Missingness Super Learner C(A, W)
  sl_C <- NULL
  if (roles$has_censoring) {
    c_wfs <- normalize_learner_library(
      learners_cfg = object$learners$C,
      mode = "classification",
      role_type = "C",
      roles = roles,
      template_data = nuisance_data$C,
      default_recipe = recipes_split$C,
      has_custom_recipe = has_custom_rec
    )
    if (roles$treatment_type == "binary") {
      df_c1 <- nuisance_data$C
      df_c1[[roles$treatment_col]] <- 1
      df_c0 <- nuisance_data$C
      df_c0[[roles$treatment_col]] <- 0
      eval_C <- list(c1 = df_c1, c0 = df_c0)
    } else {
      eval_C <- list(c_obs = nuisance_data$C)
    }
    sl_C <- fit_cv_super_learner(
      workflows_list = c_wfs,
      splits = splits,
      data = nuisance_data$C,
      y_target = roles$delta,
      mode = "classification",
      eval_datasets = eval_C,
      train_subset_idx = NULL,
      weights = roles$weights,
      metalearner = object$learners$C$metalearner,
      cv_tmle = object$targeting$cv_tmle,
      save_nuisance = object$targeting$save_nuisance,
      parallel = object$targeting$parallel,
      nuisance_name = "C"
    )
    c1_hat <- if (roles$treatment_type == "binary") sl_C$sl_eval_preds$c1 else sl_C$sl_preds
    c0_hat <- if (roles$treatment_type == "binary") sl_C$sl_eval_preds$c0 else sl_C$sl_preds
  } else {
    c1_hat <- rep(1, roles$n)
    c0_hat <- rep(1, roles$n)
  }

  # Assemble cached out-of-fold nuisance predictions
  if (roles$treatment_type == "binary") {
    q1_init <- sl_Q$sl_eval_preds$q1
    q0_init <- sl_Q$sl_eval_preds$q0
    qa_init <- ifelse(roles$a_num == 1L, q1_init, q0_init)
  } else {
    qa_init <- sl_Q$sl_eval_preds$qa
    q1_init <- sl_Q$sl_eval_preds$q_shift
    q0_init <- qa_init
  }

  cache <- list(
    data = roles$data,
    y_num = roles$y_num,
    a_num = roles$a_num,
    delta = roles$delta,
    weights = roles$weights,
    weights_kind = roles$weights_kind,
    id_vec = roles$id_vec,
    outcome_col = roles$outcome_col,
    treatment_col = roles$treatment_col,
    covariate_cols = roles$covariate_cols,
    outcome_type = roles$outcome_type,
    treatment_type = roles$treatment_type,
    y_bounds = roles$y_bounds,
    has_censoring = roles$has_censoring,
    auto_censoring = roles$auto_censoring,
    q0_init = q0_init,
    q1_init = q1_init,
    qa_init = qa_init,
    g_hat = sl_g$sl_preds,
    c1_hat = c1_hat,
    c0_hat = c0_hat,
    fold_id = sl_Q$fold_id,
    v_folds = length(splits)
  )

  # 7. Execute targeting & EIF inference
  targeting_res <- run_tmle_targeting(
    y_num = cache$y_num,
    a_num = cache$a_num,
    delta = cache$delta,
    q0_init = cache$q0_init,
    q1_init = cache$q1_init,
    qa_init = cache$qa_init,
    g_hat = cache$g_hat,
    c1_hat = cache$c1_hat,
    c0_hat = cache$c0_hat,
    weights = cache$weights,
    estimand = object$estimand,
    outcome_type = cache$outcome_type,
    y_bounds = cache$y_bounds,
    g_bounds_arg = object$targeting$g_bounds,
    fluctuation = object$targeting$fluctuation,
    cluster_id = cache$id_vec,
    conf.level = 0.95
  )

  elapsed <- proc.time()[["elapsed"]] - t_start

  structure(
    list(
      spec = object,
      cache = cache,
      targeting = targeting_res,
      nuisance = list(Q = sl_Q, g = sl_g, C = sl_C),
      subgroup_idx = seq_len(roles$n),
      subgroup_expr = NULL,
      fit_time_sec = elapsed
    ),
    class = "tmle_fit"
  )
}

#' Instantaneously Retarget a Fitted TMLE Without Refitting ML Nuisance Models
#'
#' Re-runs only the positivity truncation, fluctuation submodel, and Efficient
#' Influence Function (EIF) inference using the cached out-of-fold nuisance
#' predictions stored inside a `tmle_fit` object. Because nuisance Super
#' Learning accounts for 95–99% of runtime while the targeting fluctuation
#' takes milliseconds, `retarget()` lets you run positivity sensitivity sweeps,
#' switch estimands, or inspect subgroup effects interactively with zero ML
#' refitting.
#'
#' @param object A fitted `"tmle_fit"` object.
#' @param g_bounds Optional updated positivity bounds (`"adaptive"`, scalar, or
#'   `c(low, high)`). Defaults to the bound currently used in `object`.
#' @param estimand Optional updated estimand(s) (e.g., `c("rr", "or")`,
#'   `"att"`, `"atc"`, `"ate"`). Defaults to the estimand in `object`.
#' @param subgroup Optional unquoted filtering expression (e.g. `age >= 65`)
#'   evaluated on the cached cohort data to perform instant subgroup TMLE
#'   re-targeting. When `NULL` (default), preserves the currently active
#'   subgroup; pass `subgroup = TRUE` to reset back to the full cohort.
#' @param fluctuation Optional updated fluctuation submodel (`"logistic"`,
#'   `"logistic_covariate"`, `"linear"`, `"none"`).
#' @param conf.level Confidence level for intervals; defaults to the level stored
#'   in `object` (`0.95` initially).
#'
#' @return An updated `"tmle_fit"` object.
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, outcome_type = "binary", seed = 1)
#' fit <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
#'   set_resamples(v_folds = 3L) |>
#'   fit(df)
#' fit_strict <- retarget(fit, g_bounds = c(0.05, 0.95))
#' fit_ratios <- retarget(fit, estimand = c("ate", "rr", "or"))
#' fit_older <- retarget(fit, subgroup = age >= 50)
#'
#' @export
retarget <- function(object,
                     g_bounds = NULL,
                     estimand = NULL,
                     subgroup = NULL,
                     fluctuation = NULL,
                     conf.level = NULL) {
  if (!inherits(object, "tmle_fit")) {
    cli::cli_abort("{.arg object} must be a {.cls tmle_fit} object.")
  }

  cache <- object$cache
  new_spec <- object$spec
  eff_conf <- if (is.null(conf.level)) {
    if (!is.null(object$targeting$conf_level)) object$targeting$conf_level else 0.95
  } else {
    if (!is.numeric(conf.level) || length(conf.level) != 1L || conf.level <= 0 || conf.level >= 1) {
      cli::cli_abort("{.arg conf.level} must be strictly between 0 and 1.")
    }
    conf.level
  }

  if (!is.null(g_bounds)) {
    validate_g_bounds_arg(g_bounds)
    new_spec$targeting$g_bounds <- g_bounds
  }
  if (!is.null(fluctuation)) {
    fluctuation <- match.arg(fluctuation, c("logistic", "logistic_covariate", "linear", "none"))
    new_spec$targeting$fluctuation <- fluctuation
  }
  if (!is.null(estimand)) {
    est_norm <- normalize_estimand(estimand)
    if (!identical(est_norm$kind, new_spec$estimand$kind)) {
      cli::cli_abort(
        "Cannot switch between binary point-treatment estimands and continuous Modified Treatment Policy estimands inside {.fn retarget}; re-fit with {.fn fit}."
      )
    }
    new_spec$estimand <- est_norm
    new_spec$raw_estimand <- estimand
  }

  sub_quo <- rlang::enquo(subgroup)
  if (!rlang::quo_is_null(sub_quo)) {
    sub_raw_expr <- rlang::quo_get_expr(sub_quo)
    if (isTRUE(sub_raw_expr)) {
      idx <- seq_len(nrow(cache$data))
      sub_expr <- NULL
    } else {
      mask <- rlang::eval_tidy(sub_quo, data = cache$data)
      if (!is.logical(mask) || length(mask) != nrow(cache$data)) {
        cli::cli_abort("{.arg subgroup} must evaluate to a logical vector of length {nrow(cache$data)}.")
      }
      mask[is.na(mask)] <- FALSE
      idx <- which(mask)
      if (length(idx) < 10L) {
        cli::cli_abort("Subgroup contains fewer than 10 observations ({length(idx)}).")
      }
      sub_expr <- rlang::quo_text(sub_quo)
    }
  } else {
    idx <- if (!is.null(object$subgroup_idx)) object$subgroup_idx else seq_len(nrow(cache$data))
    sub_expr <- object$subgroup_expr
  }

  y_sub <- cache$y_num[idx]
  a_sub <- cache$a_num[idx]
  d_sub <- cache$delta[idx]
  y_bounds_sub <- if (cache$outcome_type == "continuous") {
    range(y_sub[d_sub == 1L], na.rm = TRUE)
  } else {
    c(0, 1)
  }

  if (cache$treatment_type == "binary") {
    if (length(unique(a_sub)) < 2L) {
      cli::cli_abort("Selected subgroup does not contain both treated and control observations.")
    }
  }

  targeting_res <- run_tmle_targeting(
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
    estimand = new_spec$estimand,
    outcome_type = cache$outcome_type,
    y_bounds = y_bounds_sub,
    g_bounds_arg = new_spec$targeting$g_bounds,
    fluctuation = new_spec$targeting$fluctuation,
    cluster_id = if (!is.null(cache$id_vec)) cache$id_vec[idx] else NULL,
    conf.level = eff_conf
  )

  out <- object
  out$spec <- new_spec
  out$targeting <- targeting_res
  out$subgroup_idx <- idx
  out$subgroup_expr <- sub_expr
  out
}
