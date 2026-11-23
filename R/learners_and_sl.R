#' Build Default Nuisance Learner Library
#'
#' Creates a fast, diverse default `parsnip` learner library (`glm`, `glmnet`,
#' and tree/random forest) for regression or binary classification nuisances.
#'
#' @param mode `"regression"` or `"classification"`.
#'
#' @return A named list of `parsnip` `model_spec` objects.
#' @export
default_nuisance_learners <- function(mode = c("regression", "classification")) {
  mode <- match.arg(mode)

  if (mode == "regression") {
    learners <- list(
      glm = parsnip::linear_reg() |>
        parsnip::set_engine("glm")
    )
    if (requireNamespace("glmnet", quietly = TRUE)) {
      learners$glmnet <- parsnip::linear_reg(penalty = 0.01, mixture = 0.5) |>
        parsnip::set_engine("glmnet")
    }
    if (requireNamespace("ranger", quietly = TRUE)) {
      learners$ranger <- parsnip::rand_forest(trees = 150L, min_n = 10L) |>
        parsnip::set_engine("ranger", num.threads = 1L, seed = 1L) |>
        parsnip::set_mode("regression")
    } else if (requireNamespace("rpart", quietly = TRUE)) {
      learners$rpart <- parsnip::decision_tree(tree_depth = 5L, min_n = 10L) |>
        parsnip::set_engine("rpart") |>
        parsnip::set_mode("regression")
    }
    return(learners)
  }

  learners <- list(
    glm = parsnip::logistic_reg() |>
      parsnip::set_engine("glm")
  )
  if (requireNamespace("glmnet", quietly = TRUE)) {
    learners$glmnet <- parsnip::logistic_reg(penalty = 0.01, mixture = 0.5) |>
      parsnip::set_engine("glmnet")
  }
  if (requireNamespace("ranger", quietly = TRUE)) {
    learners$ranger <- parsnip::rand_forest(trees = 150L, min_n = 10L) |>
      parsnip::set_engine("ranger", num.threads = 1L, seed = 1L) |>
      parsnip::set_mode("classification")
  } else if (requireNamespace("rpart", quietly = TRUE)) {
    learners$rpart <- parsnip::decision_tree(tree_depth = 5L, min_n = 10L) |>
      parsnip::set_engine("rpart") |>
      parsnip::set_mode("classification")
  }
  learners
}

#' Normalize Polymorphic Learner Inputs into Role-Safe Workflows
#' @noRd
normalize_learner_library <- function(learners_cfg,
                                      mode,
                                      role_type,
                                      roles,
                                      template_data,
                                      default_recipe) {
  spec_input <- learners_cfg$spec
  if (identical(spec_input, "default") || is.null(spec_input)) {
    spec_list <- default_nuisance_learners(mode = mode)
  } else if (inherits(spec_input, "model_spec")) {
    nm <- paste0(class(spec_input)[1L], "_", spec_input$engine)
    spec_list <- stats::setNames(list(spec_input), nm)
  } else if (inherits(spec_input, "workflow")) {
    spec_list <- list(custom_wf = spec_input)
  } else if (inherits(spec_input, "workflow_set")) {
    ids <- spec_input$wflow_id
    spec_list <- stats::setNames(
      lapply(ids, function(id) workflowsets::extract_workflow(spec_input, id = id)),
      ids
    )
  } else if (is.list(spec_input)) {
    nms <- names(spec_input)
    if (is.null(nms) || any(nms == "")) {
      nms <- paste0("learner_", seq_along(spec_input))
    }
    spec_list <- stats::setNames(spec_input, nms)
  } else {
    cli::cli_abort("Unsupported learner specification of class {.cls {class(spec_input)}}.")
  }

  # Ensure mode compatibility and role-safe preprocessor for every workflow
  wf_list <- lapply(names(spec_list), function(nm) {
    item <- spec_list[[nm]]
    if (inherits(item, "model_spec")) {
      ms <- parsnip::set_mode(item, mode)
      wf <- workflows::workflow() |>
        workflows::add_model(ms) |>
        workflows::add_recipe(default_recipe)
    } else if (inherits(item, "workflow")) {
      wf <- item
      ms <- workflows::extract_spec_parsnip(wf)
      if (!identical(ms$mode, mode)) {
        ms <- parsnip::set_mode(ms, mode)
        wf <- workflows::update_model(wf, ms)
      }
      # Inspect preprocessor: adapt recipe or formula for the exact causal role
      has_rec <- tryCatch(!is.null(workflows::extract_preprocessor(wf)), error = function(e) FALSE)
      if (has_rec) {
        prep_obj <- workflows::extract_preprocessor(wf)
        if (inherits(prep_obj, "recipe")) {
          safe_rec <- adapt_recipe_for_role(prep_obj, role_type = role_type, roles = roles, template_data = template_data)
          wf <- workflows::update_recipe(wf, safe_rec)
        } else {
          # Formula or variables preprocessor: bind role-safe default_recipe instead if mismatched
          wf <- workflows::remove_formula(wf) |>
            workflows::add_recipe(default_recipe)
        }
      } else {
        wf <- workflows::add_recipe(wf, default_recipe)
      }
    } else {
      cli::cli_abort("All elements of learner list must be {.cls model_spec} or {.cls workflow} objects.")
    }

    # Attach hardhat case weights if present and supported by parsnip engine
    if (roles$weights_kind != "none" && ".tmle_wt" %in% names(template_data)) {
      ms_cur <- workflows::extract_spec_parsnip(wf)
      supports_wt <- isTRUE(tryCatch(
        parsnip::get_from_env(paste0(class(ms_cur)[1L], "_fit")) |>
          (\(tbl) {
            row <- tbl[tbl$engine == ms_cur$engine & tbl$mode == ms_cur$mode, , drop = FALSE]
            any(grepl("weights", as.character(row$value[[1L]]$interface))) || TRUE
          })(),
        error = function(e) FALSE
      ))
      # Safer: try add_case_weights; if engine rejects during fit, fallback cleanly in fit_single_workflow
      wf <- tryCatch(
        workflows::add_case_weights(wf, col = !!rlang::sym(".tmle_wt")),
        error = function(e) wf
      )
    }

    wf
  })
  names(wf_list) <- names(spec_list)
  wf_list
}

#' Fit a Single Tidymodels Workflow Safely (with Case-Weight Fallback)
#' @noRd
fit_single_workflow_safe <- function(wf, train_df) {
  fit_quietly <- function(workflow_obj, df) {
    withCallingHandlers(
      workflows::fit(workflow_obj, data = df),
      warning = function(w) {
        if (grepl("non-integer #successes", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  tryCatch(
    fit_quietly(wf, train_df),
    error = function(e) {
      # If engine doesn't accept case_weights, drop case_weights and retry
      has_wt <- ".tmle_wt" %in% names(train_df) &&
        !is.null(wf$pre$actions$case_weights)
      if (has_wt) {
        wf_nowt <- tryCatch(workflows::remove_case_weights(wf), error = function(e2) wf)
        return(tryCatch(fit_quietly(wf_nowt, train_df), error = function(e3) NULL))
      }
      NULL
    }
  )
}

#' Predict Numeric or Probability-of-1 Vector from Fitted Workflow
#' @noRd
predict_workflow_safe <- function(fitted_wf, new_df, mode) {
  if (is.null(fitted_wf)) {
    return(rep(NA_real_, nrow(new_df)))
  }
  if (mode == "regression") {
    preds <- tryCatch(
      stats::predict(fitted_wf, new_data = new_df, type = "numeric")$.pred,
      error = function(e) rep(NA_real_, nrow(new_df))
    )
    return(as.numeric(preds))
  }
  preds_tbl <- tryCatch(
    stats::predict(fitted_wf, new_data = new_df, type = "prob"),
    error = function(e) NULL
  )
  if (is.null(preds_tbl)) {
    return(rep(NA_real_, nrow(new_df)))
  }
  if (".pred_1" %in% names(preds_tbl)) {
    p <- as.numeric(preds_tbl$.pred_1)
  } else if (ncol(preds_tbl) >= 2L) {
    p <- as.numeric(preds_tbl[[2L]])
  } else {
    p <- as.numeric(preds_tbl[[1L]])
  }
  pmin(pmax(p, 1e-6), 1 - 1e-6)
}

#' Project a Numeric Vector onto the Probability Simplex
#' @noRd
project_simplex <- function(v) {
  n <- length(v)
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u) - 1
  ind <- seq_len(n)
  cond <- u - cssv / ind > 0
  rho <- max(ind[cond])
  theta <- cssv[rho] / rho
  pmax(v - theta, 0)
}

#' Solve Simplex-Constrained Super Learner Metalearner Weights
#'
#' Solves convex optimization over \eqn{\Delta^{K-1}} using FISTA accelerated
#' projected gradient descent so weights satisfy \eqn{\alpha_k \ge 0},
#' \eqn{\sum_k \alpha_k = 1}, and naturally admit sparse zeros.
#'
#' @noRd
solve_sl_weights <- function(oof_mat,
                             y_target,
                             mode,
                             weights = NULL,
                             metalearner = "auto") {
  K <- ncol(oof_mat)
  learner_names <- colnames(oof_mat)
  if (is.null(learner_names)) {
    learner_names <- paste0("learner_", seq_len(K))
  }

  valid_cols <- which(colSums(!is.na(oof_mat)) == nrow(oof_mat) &
                        apply(oof_mat, 2L, function(z) all(is.finite(z))))
  if (length(valid_cols) == 0L) {
    cli::cli_abort("All candidate nuisance learners failed to produce valid out-of-fold predictions.")
  }

  if (length(valid_cols) == 1L) {
    w_out <- stats::setNames(numeric(K), learner_names)
    w_out[valid_cols] <- 1
    return(w_out)
  }

  Z <- oof_mat[, valid_cols, drop = FALSE]
  Kv <- ncol(Z)
  n <- nrow(Z)
  wt <- if (is.null(weights)) rep(1 / n, n) else as.numeric(weights) / sum(weights)

  if (identical(metalearner, "auto")) {
    metalearner <- if (mode == "regression") "nnls" else "logloss"
  }

  # Compute per-learner CV loss for "best" or fallback
  learner_losses <- vapply(seq_len(Kv), function(k) {
    zk <- Z[, k]
    if (metalearner == "logloss") {
      zk <- pmin(pmax(zk, 1e-6), 1 - 1e-6)
      -sum(wt * (y_target * log(zk) + (1 - y_target) * log(1 - zk)))
    } else {
      sum(wt * (y_target - zk)^2)
    }
  }, numeric(1L))

  if (identical(metalearner, "best")) {
    best_idx <- valid_cols[which.min(learner_losses)]
    w_out <- stats::setNames(numeric(K), learner_names)
    w_out[best_idx] <- 1
    return(w_out)
  }

  # FISTA / Accelerated Projected Gradient on the probability simplex
  if (metalearner == "nnls") {
    # Weighted Gram matrix for 0.5 * sum(wt * (y - Z %*% alpha)^2)
    Zw <- Z * sqrt(wt)
    yw <- y_target * sqrt(wt)
    G <- crossprod(Zw)
    g0 <- crossprod(Zw, yw)
    L <- max(eigen(G, symmetric = TRUE, only.values = TRUE)$values, 1e-8)

    alpha <- rep(1 / Kv, Kv)
    z_var <- alpha
    t_k <- 1
    for (iter in seq_len(400L)) {
      grad <- as.numeric(G %*% z_var - g0)
      alpha_next <- project_simplex(z_var - grad / L)
      if (max(abs(alpha_next - alpha)) < 1e-9) {
        alpha <- alpha_next
        break
      }
      t_next <- 0.5 * (1 + sqrt(1 + 4 * t_k^2))
      z_var <- alpha_next + ((t_k - 1) / t_next) * (alpha_next - alpha)
      alpha <- alpha_next
      t_k <- t_next
    }
  } else {
    # Binomial log-loss minimization over convex combination of probabilities
    Zc <- pmin(pmax(Z, 1e-6), 1 - 1e-6)
    loss_fn <- function(a) {
      p <- pmin(pmax(as.numeric(Zc %*% a), 1e-6), 1 - 1e-6)
      -sum(wt * (y_target * log(p) + (1 - y_target) * log(1 - p)))
    }
    grad_fn <- function(a) {
      p <- pmin(pmax(as.numeric(Zc %*% a), 1e-6), 1 - 1e-6)
      r <- wt * ((1 - y_target) / (1 - p) - y_target / p)
      as.numeric(crossprod(Zc, r))
    }

    alpha <- rep(1 / Kv, Kv)
    step <- 0.25
    cur_loss <- loss_fn(alpha)
    for (iter in seq_len(350L)) {
      gr <- grad_fn(alpha)
      # Backtracking projected gradient
      step_cur <- step
      accepted <- FALSE
      delta_alpha <- 0
      for (bt in seq_len(25L)) {
        cand <- project_simplex(alpha - step_cur * gr)
        cand_loss <- loss_fn(cand)
        if (cand_loss <= cur_loss - 1e-4 * sum(gr * (alpha - cand)) + 1e-12) {
          delta_alpha <- max(abs(cand - alpha))
          alpha <- cand
          cur_loss <- cand_loss
          step <- min(step_cur * 1.5, 20.0)
          accepted <- TRUE
          break
        }
        step_cur <- step_cur * 0.5
      }
      if (!accepted || delta_alpha < 1e-8) {
        break
      }
    }
  }

  # Clean tiny numerical weights < 1e-6 and renormalize
  alpha[alpha < 1e-6] <- 0
  if (sum(alpha) == 0) {
    alpha[which.min(learner_losses)] <- 1
  } else {
    alpha <- alpha / sum(alpha)
  }

  w_out <- stats::setNames(numeric(K), learner_names)
  w_out[valid_cols] <- alpha
  w_out
}

#' Cross-Fit Super Learner Over rsample Folds and Candidate Workflows
#' @noRd
fit_cv_super_learner <- function(workflows_list,
                                 splits,
                                 data,
                                 y_target,
                                 mode,
                                 eval_datasets = list(),
                                 train_subset_idx = NULL,
                                 weights = NULL,
                                 metalearner = "auto",
                                 cv_tmle = TRUE,
                                 save_nuisance = "predictions",
                                 parallel = FALSE,
                                 nuisance_name = "Q") {
  n <- nrow(data)
  K <- length(workflows_list)
  learner_names <- names(workflows_list)
  V <- length(splits)

  # Helper to run a single (fold, learner) task (self-contained for mirai workers)
  run_fold_learner <- function(v_idx, k_idx) {
    sp <- splits[[v_idx]]
    wf <- workflows_list[[k_idx]]
    in_id <- if (cv_tmle) sp$in_id else seq_len(n)
    out_id <- if (cv_tmle) setdiff(seq_len(n), in_id) else seq_len(n)

    if (!is.null(train_subset_idx)) {
      in_id <- intersect(in_id, train_subset_idx)
    }

    train_df <- data[in_id, , drop = FALSE]
    valid_df <- data[out_id, , drop = FALSE]

    fitted_wf <- fit_fn(wf, train_df)
    p_valid <- pred_fn(fitted_wf, valid_df, mode = mode)

    p_eval <- lapply(eval_datasets, function(edf) {
      pred_fn(fitted_wf, edf[out_id, , drop = FALSE], mode = mode)
    })

    kept_wf <- NULL
    if (identical(save_nuisance, "full")) {
      kept_wf <- fitted_wf
    } else if (identical(save_nuisance, "butcher") && !is.null(fitted_wf) &&
               requireNamespace("butcher", quietly = TRUE)) {
      kept_wf <- butcher::butcher(fitted_wf)
    }

    list(
      v_idx = v_idx,
      k_idx = k_idx,
      out_id = out_id,
      p_valid = p_valid,
      p_eval = p_eval,
      fitted_wf = kept_wf
    )
  }

  fit_fn <- fit_single_workflow_safe
  pred_fn <- predict_workflow_safe
  tasks <- expand.grid(v_idx = seq_len(V), k_idx = seq_len(K))

  use_mirai <- isTRUE(parallel) &&
    requireNamespace("mirai", quietly = TRUE) &&
    isTRUE(tryCatch(mirai::daemons_set(), error = function(e) FALSE))

  if (use_mirai) {
    task_list <- lapply(seq_len(nrow(tasks)), function(r) {
      list(v_idx = tasks$v_idx[r], k_idx = tasks$k_idx[r])
    })
    shared_env <- list(
      splits = splits,
      workflows_list = workflows_list,
      cv_tmle = cv_tmle,
      n = n,
      train_subset_idx = train_subset_idx,
      data = data,
      eval_datasets = eval_datasets,
      mode = mode,
      save_nuisance = save_nuisance,
      fit_fn = fit_single_workflow_safe,
      pred_fn = predict_workflow_safe
    )
    results <- mirai::mirai_map(
      .x = task_list,
      .f = function(t, shared) {
        sp <- shared$splits[[t$v_idx]]
        wf <- shared$workflows_list[[t$k_idx]]
        in_id <- if (shared$cv_tmle) sp$in_id else seq_len(shared$n)
        out_id <- if (shared$cv_tmle) setdiff(seq_len(shared$n), in_id) else seq_len(shared$n)
        if (!is.null(shared$train_subset_idx)) {
          in_id <- intersect(in_id, shared$train_subset_idx)
        }
        train_df <- shared$data[in_id, , drop = FALSE]
        valid_df <- shared$data[out_id, , drop = FALSE]
        fitted_wf <- shared$fit_fn(wf, train_df)
        p_valid <- shared$pred_fn(fitted_wf, valid_df, mode = shared$mode)
        p_eval <- lapply(shared$eval_datasets, function(edf) {
          shared$pred_fn(fitted_wf, edf[out_id, , drop = FALSE], mode = shared$mode)
        })
        kept_wf <- NULL
        if (identical(shared$save_nuisance, "full")) {
          kept_wf <- fitted_wf
        } else if (identical(shared$save_nuisance, "butcher") && !is.null(fitted_wf) &&
                   requireNamespace("butcher", quietly = TRUE)) {
          kept_wf <- butcher::butcher(fitted_wf)
        }
        list(
          v_idx = t$v_idx,
          k_idx = t$k_idx,
          out_id = out_id,
          p_valid = p_valid,
          p_eval = p_eval,
          fitted_wf = kept_wf
        )
      },
      .args = list(shared = shared_env)
    )[]
    bad_m <- vapply(results, function(r) inherits(r, "miraiError") || inherits(r, "errorValue"), logical(1L))
    if (any(bad_m)) {
      cli::cli_abort("Parallel {.pkg mirai} worker failed: {as.character(results[[which(bad_m)[1L]]])}")
    }
  } else {
    results <- lapply(seq_len(nrow(tasks)), function(r) {
      run_fold_learner(tasks$v_idx[r], tasks$k_idx[r])
    })
  }

  # Assemble out-of-fold prediction matrices
  oof_obs <- matrix(NA_real_, nrow = n, ncol = K)
  colnames(oof_obs) <- learner_names

  eval_names <- names(eval_datasets)
  oof_eval <- stats::setNames(
    lapply(eval_names, function(...) {
      m <- matrix(NA_real_, nrow = n, ncol = K)
      colnames(m) <- learner_names
      m
    }),
    eval_names
  )

  fold_id_vec <- integer(n)
  kept_models <- vector("list", length(results))

  for (idx in seq_along(results)) {
    res <- results[[idx]]
    out_id <- res$out_id
    k_idx <- res$k_idx
    oof_obs[out_id, k_idx] <- res$p_valid
    if (k_idx == 1L) {
      fold_id_vec[out_id] <- res$v_idx
    }
    for (enm in eval_names) {
      oof_eval[[enm]][out_id, k_idx] <- res$p_eval[[enm]]
    }
    if (!is.null(res$fitted_wf)) {
      kept_models[[idx]] <- list(fold = res$v_idx, learner = learner_names[k_idx], model = res$fitted_wf)
    }
  }

  # Solve Super Learner metalearner weights on observed/training subset
  fit_idx <- if (is.null(train_subset_idx)) seq_len(n) else train_subset_idx
  wt_sub <- if (is.null(weights)) NULL else weights[fit_idx]

  sl_weights <- solve_sl_weights(
    oof_mat = oof_obs[fit_idx, , drop = FALSE],
    y_target = y_target[fit_idx],
    mode = mode,
    weights = wt_sub,
    metalearner = metalearner
  )

  # Also compute per-fold weights for diagnostics in extract_sl_weights()
  fold_weights <- purrr::map_dfr(seq_len(V), function(v) {
    v_fit_idx <- intersect(fit_idx, which(fold_id_vec == v))
    if (length(v_fit_idx) < 5L) {
      fw <- sl_weights
    } else {
      wt_v <- if (is.null(weights)) NULL else weights[v_fit_idx]
      fw <- tryCatch(
        solve_sl_weights(
          oof_mat = oof_obs[v_fit_idx, , drop = FALSE],
          y_target = y_target[v_fit_idx],
          mode = mode,
          weights = wt_v,
          metalearner = metalearner
        ),
        error = function(e) sl_weights
      )
    }
    tibble::tibble(
      nuisance = nuisance_name,
      fold = v,
      learner = learner_names,
      weight = as.numeric(fw)
    )
  })

  # Compute ensemble OOF predictions
  valid_k <- which(sl_weights > 0)
  combine_preds <- function(mat) {
    as.numeric(mat[, valid_k, drop = FALSE] %*% sl_weights[valid_k])
  }

  sl_obs <- combine_preds(oof_obs)
  sl_eval <- lapply(oof_eval, combine_preds)

  # Compute yardstick CV metrics for all base learners + Super Learner
  metrics_tbl <- compute_nuisance_metrics(
    oof_obs = oof_obs,
    sl_obs = sl_obs,
    y_target = y_target,
    fit_idx = fit_idx,
    fold_id_vec = fold_id_vec,
    mode = mode,
    nuisance_name = nuisance_name
  )

  list(
    nuisance = nuisance_name,
    mode = mode,
    sl_preds = sl_obs,
    sl_eval_preds = sl_eval,
    oof_matrix = oof_obs,
    sl_weights = sl_weights,
    fold_weights = fold_weights,
    metrics = metrics_tbl,
    fold_id = fold_id_vec,
    models = if (identical(save_nuisance, "predictions")) NULL else Filter(Negate(is.null), kept_models)
  )
}

#' Cross-Fit Stacked Classification Density-Ratio Super Learner for LMTP Policies
#' @noRd
fit_cv_density_ratio_sl <- function(workflows_list,
                                    splits,
                                    data,
                                    treatment_col,
                                    policy,
                                    weights = NULL,
                                    metalearner = "auto",
                                    cv_tmle = TRUE,
                                    save_nuisance = "predictions",
                                    parallel = FALSE) {
  n <- nrow(data)
  K <- length(workflows_list)
  learner_names <- names(workflows_list)
  V <- length(splits)

  a_obs <- data[[treatment_col]]
  a_shifted <- apply_policy(policy, a_obs)

  run_fold_learner <- function(v_idx, k_idx) {
    sp <- splits[[v_idx]]
    wf <- workflows_list[[k_idx]]
    in_id <- if (cv_tmle) sp$in_id else seq_len(n)
    out_id <- if (cv_tmle) setdiff(seq_len(n), in_id) else seq_len(n)

    tr_0 <- data[in_id, , drop = FALSE]
    tr_0[[".lambda"]] <- factor(0L, levels = c("0", "1"))
    tr_1 <- data[in_id, , drop = FALSE]
    tr_1[[treatment_col]] <- a_shifted[in_id]
    tr_1[[".lambda"]] <- factor(1L, levels = c("0", "1"))
    train_stacked <- dplyr::bind_rows(tr_0, tr_1)

    va_0 <- data[out_id, , drop = FALSE]
    va_0[[".lambda"]] <- factor(0L, levels = c("0", "1"))
    va_1 <- data[out_id, , drop = FALSE]
    va_1[[treatment_col]] <- a_shifted[out_id]
    va_1[[".lambda"]] <- factor(1L, levels = c("0", "1"))
    valid_stacked <- dplyr::bind_rows(va_0, va_1)

    fitted_wf <- fit_single_workflow_safe(wf, train_stacked)
    p_valid_stacked <- predict_workflow_safe(fitted_wf, valid_stacked, mode = "classification")
    p_obs_arm <- predict_workflow_safe(fitted_wf, va_0, mode = "classification")

    kept_wf <- NULL
    if (identical(save_nuisance, "full")) {
      kept_wf <- fitted_wf
    } else if (identical(save_nuisance, "butcher") && !is.null(fitted_wf) &&
               requireNamespace("butcher", quietly = TRUE)) {
      kept_wf <- butcher::butcher(fitted_wf)
    }

    list(
      v_idx = v_idx,
      k_idx = k_idx,
      out_id = out_id,
      p_valid_stacked = p_valid_stacked,
      y_valid_stacked = c(rep(0, length(out_id)), rep(1, length(out_id))),
      p_obs_arm = p_obs_arm,
      fitted_wf = kept_wf
    )
  }

  tasks <- expand.grid(v_idx = seq_len(V), k_idx = seq_len(K))
  use_mirai <- isTRUE(parallel) &&
    requireNamespace("mirai", quietly = TRUE) &&
    isTRUE(tryCatch(mirai::daemons_set(), error = function(e) FALSE))

  if (use_mirai) {
    task_list <- lapply(seq_len(nrow(tasks)), function(r) {
      list(v_idx = tasks$v_idx[r], k_idx = tasks$k_idx[r])
    })
    shared_env <- list(
      splits = splits,
      workflows_list = workflows_list,
      cv_tmle = cv_tmle,
      n = n,
      data = data,
      treatment_col = treatment_col,
      a_shifted = a_shifted,
      save_nuisance = save_nuisance,
      fit_fn = fit_single_workflow_safe,
      pred_fn = predict_workflow_safe
    )
    results <- mirai::mirai_map(
      .x = task_list,
      .f = function(t, shared) {
        sp <- shared$splits[[t$v_idx]]
        wf <- shared$workflows_list[[t$k_idx]]
        in_id <- if (shared$cv_tmle) sp$in_id else seq_len(shared$n)
        out_id <- if (shared$cv_tmle) setdiff(seq_len(shared$n), in_id) else seq_len(shared$n)

        tr_0 <- shared$data[in_id, , drop = FALSE]
        tr_0[[".lambda"]] <- factor(0L, levels = c("0", "1"))
        tr_1 <- shared$data[in_id, , drop = FALSE]
        tr_1[[shared$treatment_col]] <- shared$a_shifted[in_id]
        tr_1[[".lambda"]] <- factor(1L, levels = c("0", "1"))
        train_stacked <- dplyr::bind_rows(tr_0, tr_1)

        va_0 <- shared$data[out_id, , drop = FALSE]
        va_0[[".lambda"]] <- factor(0L, levels = c("0", "1"))
        va_1 <- shared$data[out_id, , drop = FALSE]
        va_1[[shared$treatment_col]] <- shared$a_shifted[out_id]
        va_1[[".lambda"]] <- factor(1L, levels = c("0", "1"))
        valid_stacked <- dplyr::bind_rows(va_0, va_1)

        fitted_wf <- shared$fit_fn(wf, train_stacked)
        p_valid_stacked <- shared$pred_fn(fitted_wf, valid_stacked, mode = "classification")
        p_obs_arm <- shared$pred_fn(fitted_wf, va_0, mode = "classification")

        kept_wf <- NULL
        if (identical(shared$save_nuisance, "full")) {
          kept_wf <- fitted_wf
        } else if (identical(shared$save_nuisance, "butcher") && !is.null(fitted_wf) &&
                   requireNamespace("butcher", quietly = TRUE)) {
          kept_wf <- butcher::butcher(fitted_wf)
        }

        list(
          v_idx = t$v_idx,
          k_idx = t$k_idx,
          out_id = out_id,
          p_valid_stacked = p_valid_stacked,
          y_valid_stacked = c(rep(0, length(out_id)), rep(1, length(out_id))),
          p_obs_arm = p_obs_arm,
          fitted_wf = kept_wf
        )
      },
      .args = list(shared = shared_env)
    )[]
    bad_m <- vapply(results, function(r) inherits(r, "miraiError") || inherits(r, "errorValue"), logical(1L))
    if (any(bad_m)) {
      cli::cli_abort("Parallel {.pkg mirai} worker failed: {as.character(results[[which(bad_m)[1L]]])}")
    }
  } else {
    results <- lapply(seq_len(nrow(tasks)), function(r) {
      run_fold_learner(tasks$v_idx[r], tasks$k_idx[r])
    })
  }

  oof_stacked <- matrix(NA_real_, nrow = 2L * n, ncol = K)
  colnames(oof_stacked) <- learner_names
  y_stacked <- c(rep(0, n), rep(1, n))

  oof_prob_lambda1_at_a <- matrix(NA_real_, nrow = n, ncol = K)
  colnames(oof_prob_lambda1_at_a) <- learner_names
  fold_id_vec <- integer(n)
  kept_models <- vector("list", length(results))

  for (idx in seq_along(results)) {
    res <- results[[idx]]
    out_id <- res$out_id
    k_idx <- res$k_idx
    m_out <- length(out_id)
    oof_stacked[out_id, k_idx] <- res$p_valid_stacked[seq_len(m_out)]
    oof_stacked[n + out_id, k_idx] <- res$p_valid_stacked[m_out + seq_len(m_out)]
    oof_prob_lambda1_at_a[out_id, k_idx] <- res$p_obs_arm
    if (k_idx == 1L) {
      fold_id_vec[out_id] <- res$v_idx
    }
    if (!is.null(res$fitted_wf)) {
      kept_models[[idx]] <- list(fold = res$v_idx, learner = learner_names[k_idx], model = res$fitted_wf)
    }
  }

  wt_stacked <- if (is.null(weights)) NULL else c(weights, weights)
  sl_weights <- solve_sl_weights(
    oof_mat = oof_stacked,
    y_target = y_stacked,
    mode = "classification",
    weights = wt_stacked,
    metalearner = metalearner
  )

  valid_k <- which(sl_weights > 0)
  p_lambda1 <- as.numeric(oof_prob_lambda1_at_a[, valid_k, drop = FALSE] %*% sl_weights[valid_k])
  p_lambda1_stacked <- as.numeric(oof_stacked[, valid_k, drop = FALSE] %*% sl_weights[valid_k])

  fold_weights <- tibble::tibble(
    nuisance = "g",
    fold = rep(seq_len(V), each = K),
    learner = rep(learner_names, times = V),
    weight = rep(as.numeric(sl_weights), times = V)
  )

  metrics_tbl <- compute_nuisance_metrics(
    oof_obs = oof_stacked,
    sl_obs = p_lambda1_stacked,
    y_target = y_stacked,
    fit_idx = seq_len(2L * n),
    fold_id_vec = c(fold_id_vec, fold_id_vec),
    mode = "classification",
    nuisance_name = "g"
  )

  list(
    nuisance = "g",
    mode = "classification",
    sl_preds = p_lambda1,
    sl_eval_preds = list(),
    oof_matrix = oof_prob_lambda1_at_a,
    sl_weights = sl_weights,
    fold_weights = fold_weights,
    metrics = metrics_tbl,
    fold_id = fold_id_vec,
    models = if (identical(save_nuisance, "predictions")) NULL else Filter(Negate(is.null), kept_models)
  )
}

#' Compute Yardstick Metrics for Base Learners and Super Learner
#' @noRd
compute_nuisance_metrics <- function(oof_obs,
                                     sl_obs,
                                     y_target,
                                     fit_idx,
                                     fold_id_vec,
                                     mode,
                                     nuisance_name) {
  learner_names <- c(colnames(oof_obs), "SuperLearner")
  pred_cols <- cbind(oof_obs, SuperLearner = sl_obs)
  y_sub <- y_target[fit_idx]
  folds_sub <- fold_id_vec[fit_idx]
  unique_folds <- sort(unique(folds_sub))

  rows <- purrr::map_dfr(seq_along(learner_names), function(j) {
    nm <- learner_names[j]
    p_all <- pred_cols[fit_idx, j]
    if (anyNA(p_all)) {
      return(tibble::tibble())
    }
    purrr::map_dfr(unique_folds, function(v) {
      idx_v <- which(folds_sub == v)
      yv <- y_sub[idx_v]
      pv <- p_all[idx_v]
      if (mode == "regression") {
        rmse_val <- yardstick::rmse_vec(truth = yv, estimate = pv)
        rsq_val <- if (stats::sd(yv) > 0 && stats::sd(pv) > 0) {
          yardstick::rsq_vec(truth = yv, estimate = pv)
        } else {
          NA_real_
        }
        tibble::tibble(
          nuisance = nuisance_name,
          learner = nm,
          fold = v,
          .metric = c("rmse", "rsq"),
          .estimator = "standard",
          .estimate = c(rmse_val, rsq_val)
        )
      } else {
        y_fac <- factor(as.character(as.integer(yv)), levels = c("0", "1"))
        pv_c <- pmin(pmax(pv, 1e-6), 1 - 1e-6)
        ll_val <- yardstick::mn_log_loss_vec(truth = y_fac, estimate = pv_c, event_level = "second")
        bc_val <- yardstick::brier_class_vec(truth = y_fac, estimate = pv_c, event_level = "second")
        auc_val <- if (length(unique(yv)) == 2L) {
          yardstick::roc_auc_vec(truth = y_fac, estimate = pv_c, event_level = "second")
        } else {
          NA_real_
        }
        tibble::tibble(
          nuisance = nuisance_name,
          learner = nm,
          fold = v,
          .metric = c("mn_log_loss", "roc_auc", "brier_class"),
          .estimator = "binary",
          .estimate = c(ll_val, auc_val, bc_val)
        )
      }
    })
  })
  rows
}
