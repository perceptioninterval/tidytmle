#' Resolve Causal Roles and Validate Input Data
#' @keywords internal
resolve_causal_roles <- function(spec, data) {
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a {.cls data.frame} or {.cls tibble}.")
  }
  data <- tibble::as_tibble(data)
  if (nrow(data) < 10L) {
    cli::cli_abort("{.arg data} must have at least 10 rows for TMLE estimation.")
  }

  eval_single_col <- function(quo, arg_name, required = TRUE) {
    if (rlang::quo_is_null(quo) || rlang::quo_is_missing(quo)) {
      if (required) {
        cli::cli_abort("{.arg {arg_name}} must be specified.")
      }
      return(NULL)
    }
    sel <- tidyselect::eval_select(quo, data = data)
    if (length(sel) == 0L) {
      if (required) {
        cli::cli_abort("{.arg {arg_name}} did not match any columns in {.arg data}.")
      }
      return(NULL)
    }
    if (length(sel) > 1L) {
      cli::cli_abort("{.arg {arg_name}} must select exactly 1 column, not {length(sel)} ({.val {names(sel)}}).")
    }
    names(sel)[[1L]]
  }

  outcome_col <- eval_single_col(spec$roles$outcome, "outcome", required = TRUE)
  treatment_col <- eval_single_col(spec$roles$treatment, "treatment", required = TRUE)
  censoring_col <- eval_single_col(spec$roles$censoring, "censoring", required = FALSE)
  id_col <- eval_single_col(spec$roles$id, "id", required = FALSE)
  weights_col <- eval_single_col(spec$roles$weights, "weights", required = FALSE)

  if (identical(outcome_col, treatment_col)) {
    cli::cli_abort("{.arg outcome} and {.arg treatment} cannot refer to the same column ({.val {outcome_col}}).")
  }

  special_cols <- unique(c(outcome_col, treatment_col, censoring_col, id_col, weights_col))

  if (rlang::quo_is_null(spec$roles$covariates) || rlang::quo_is_missing(spec$roles$covariates)) {
    covariate_cols <- setdiff(names(data), special_cols)
  } else {
    cov_sel <- tidyselect::eval_select(spec$roles$covariates, data = data)
    covariate_cols <- setdiff(names(cov_sel), special_cols)
  }

  if (length(covariate_cols) == 0L) {
    cli::cli_abort("At least 1 covariate column must be selected for {.arg covariates}.")
  }

  # Extract case weights
  if (!is.null(weights_col)) {
    w_raw <- data[[weights_col]]
    weights_kind <- if (hardhat::is_importance_weights(w_raw)) {
      "importance"
    } else if (hardhat::is_frequency_weights(w_raw)) {
      "frequency"
    } else {
      "numeric"
    }
    weights_vec <- as.numeric(w_raw)
    if (anyNA(weights_vec) || any(!is.finite(weights_vec)) || any(weights_vec < 0)) {
      cli::cli_abort("Case weights in {.field {weights_col}} must be finite and non-negative.")
    }
    if (sum(weights_vec) <= 0) {
      cli::cli_abort("Sum of case weights in {.field {weights_col}} must be strictly positive.")
    }
  } else {
    weights_kind <- "none"
    weights_vec <- rep(1, nrow(data))
  }

  # Extract cluster IDs
  id_vec <- NULL
  if (!is.null(id_col)) {
    id_vec <- data[[id_col]]
    if (anyNA(id_vec)) {
      cli::cli_abort("Cluster identifier {.field {id_col}} must not contain missing values.")
    }
    if (length(unique(id_vec)) < 2L) {
      cli::cli_abort("Cluster identifier {.field {id_col}} must contain at least 2 distinct clusters.")
    }
  }

  # Handle censoring & automatic informative missingness
  y_raw <- data[[outcome_col]]
  if (is.null(censoring_col)) {
    na_y <- is.na(y_raw)
    if (any(na_y)) {
      n_miss <- sum(na_y)
      pct_miss <- round(100 * n_miss / nrow(data), 1)
      cli::cli_alert_warning(
        "Detected {n_miss} missing outcome{?s} ({pct_miss}%); automatically introducing missingness indicator {.field .delta} and modeling C(A, W) = P(Delta = 1 | A, W) instead of listwise deletion."
      )
      delta_vec <- as.integer(!na_y)
      has_censoring <- TRUE
      auto_censoring <- TRUE
    } else {
      delta_vec <- rep(1L, nrow(data))
      has_censoring <- FALSE
      auto_censoring <- FALSE
    }
  } else {
    delta_vec <- coerce_binary_01(data[[censoring_col]], censoring_col)
    if (any(delta_vec == 1L & is.na(y_raw))) {
      cli::cli_abort(
        "Outcome {.field {outcome_col}} contains `NA`s on rows marked as observed ({.field {censoring_col}} == 1)."
      )
    }
    has_censoring <- any(delta_vec == 0L)
    auto_censoring <- FALSE
  }

  if (sum(delta_vec == 1L) < 5L) {
    cli::cli_abort("Too few uncensored/observed outcomes ({sum(delta_vec == 1L)}).")
  }

  # Resolve outcome vector and type
  y_obs <- y_raw[delta_vec == 1L]
  outcome_info <- parse_outcome_vector(y_raw, y_obs, delta_vec, outcome_col)

  # Resolve treatment vector and type
  trt_raw <- data[[treatment_col]]
  if (anyNA(trt_raw)) {
    cli::cli_abort("Treatment column {.field {treatment_col}} must not contain missing values.")
  }

  if (spec$estimand$kind == "binary") {
    a_num <- coerce_binary_01(trt_raw, treatment_col)
    n1 <- sum(a_num == 1L)
    n0 <- sum(a_num == 0L)
    if (n1 == 0L || n0 == 0L) {
      cli::cli_abort("Treatment {.field {treatment_col}} must contain both treated (1) and control (0) units.")
    }
    treatment_type <- "binary"
    a_factor <- factor(a_num, levels = c("0", "1"))
  } else {
    if (!is.numeric(trt_raw)) {
      cli::cli_abort(
        "Continuous treatment {.field {treatment_col}} must be numeric for Modified Treatment Policy estimands."
      )
    }
    a_num <- as.numeric(trt_raw)
    if (any(!is.finite(a_num))) {
      cli::cli_abort("Treatment {.field {treatment_col}} must contain only finite values.")
    }
    treatment_type <- "continuous"
    a_factor <- NULL
  }

  list(
    data = data,
    n = nrow(data),
    outcome_col = outcome_col,
    treatment_col = treatment_col,
    covariate_cols = covariate_cols,
    censoring_col = censoring_col,
    id_col = id_col,
    weights_col = weights_col,
    weights = weights_vec,
    weights_kind = weights_kind,
    id_vec = id_vec,
    delta = delta_vec,
    has_censoring = has_censoring,
    auto_censoring = auto_censoring,
    outcome_type = outcome_info$type,
    y_num = outcome_info$y_num,
    y_factor = outcome_info$y_factor,
    y_bounds = outcome_info$bounds,
    treatment_type = treatment_type,
    a_num = a_num,
    a_factor = a_factor,
    a_raw_template = trt_raw
  )
}

coerce_binary_01 <- function(x, col_name, call = rlang::caller_env()) {
  if (anyNA(x)) {
    cli::cli_abort("Binary column {.field {col_name}} must not contain missing values.", call = call)
  }
  if (is.logical(x)) {
    return(as.integer(x))
  }
  if (is.numeric(x) || is.integer(x)) {
    vals <- sort(unique(as.numeric(x)))
    if (!all(vals %in% c(0, 1))) {
      cli::cli_abort(
        "Numeric binary column {.field {col_name}} must only take values in {{0, 1}}, found {.val {vals}}.",
        call = call
      )
    }
    return(as.integer(x))
  }
  if (is.factor(x) || is.character(x)) {
    xf <- factor(x)
    lvls <- levels(xf)
    if (length(lvls) != 2L) {
      cli::cli_abort(
        "Categorical binary column {.field {col_name}} must have exactly 2 levels, not {length(lvls)}.",
        call = call
      )
    }
    if (identical(lvls, c("0", "1"))) {
      return(as.integer(as.character(xf)))
    }
    return(as.integer(xf) - 1L)
  }
  cli::cli_abort("Unsupported type for binary column {.field {col_name}}.", call = call)
}

parse_outcome_vector <- function(y_raw, y_obs, delta_vec, col_name, call = rlang::caller_env()) {
  if (is.logical(y_raw) || is.factor(y_raw) || is.character(y_raw)) {
    yf <- factor(y_obs)
    lvls <- levels(yf)
    if (length(lvls) != 2L) {
      cli::cli_abort(
        "Non-numeric outcome {.field {col_name}} must have exactly 2 levels on observed rows, found {length(lvls)}.",
        call = call
      )
    }
    y_num <- rep(NA_real_, length(y_raw))
    obs_idx <- which(delta_vec == 1L)
    y_num[obs_idx] <- if (is.logical(y_raw)) {
      as.numeric(y_raw[obs_idx])
    } else if (identical(lvls, c("0", "1"))) {
      as.numeric(as.character(yf))
    } else {
      as.numeric(yf) - 1
    }
    y_fac <- factor(ifelse(is.na(y_num), NA_character_, as.character(y_num)), levels = c("0", "1"))
    return(list(
      type = "binary",
      y_num = y_num,
      y_factor = y_fac,
      bounds = c(0, 1)
    ))
  }

  if (!is.numeric(y_raw)) {
    cli::cli_abort("Outcome column {.field {col_name}} must be numeric, logical, or a 2-level factor.", call = call)
  }

  y_num <- as.numeric(y_raw)
  y_num[delta_vec == 0L] <- NA_real_
  obs_vals <- y_num[delta_vec == 1L]
  if (any(!is.finite(obs_vals))) {
    cli::cli_abort("Outcome column {.field {col_name}} contains non-finite values.", call = call)
  }

  if (all(obs_vals %in% c(0, 1))) {
    y_fac <- factor(ifelse(is.na(y_num), NA_character_, as.character(as.integer(y_num))), levels = c("0", "1"))
    return(list(
      type = "binary",
      y_num = y_num,
      y_factor = y_fac,
      bounds = c(0, 1)
    ))
  }

  rng <- range(obs_vals)
  if (rng[1L] == rng[2L]) {
    cli::cli_abort("Continuous outcome {.field {col_name}} has zero variance.", call = call)
  }
  list(
    type = "continuous",
    y_num = y_num,
    y_factor = NULL,
    bounds = rng
  )
}

#' Map Counterfactual Values Back to Original Column Representation
#' @keywords internal
encode_treatment_column <- function(template, values_01_or_num) {
  if (is.factor(template)) {
    lvls <- levels(template)
    idx <- as.integer(values_01_or_num) + 1L
    return(factor(lvls[idx], levels = lvls))
  }
  if (is.logical(template)) {
    return(as.logical(values_01_or_num))
  }
  if (is.integer(template) && all(values_01_or_num %in% c(0, 1))) {
    return(as.integer(values_01_or_num))
  }
  as.numeric(values_01_or_num)
}

#' Build Clean Modeling Dataframes for Nuisance Factors
#' @keywords internal
build_nuisance_data <- function(roles) {
  d <- roles$data
  y_col <- roles$outcome_col
  a_col <- roles$treatment_col

  # Replace outcome with clean binary factor or numeric
  d_Q <- d
  if (roles$outcome_type == "binary") {
    d_Q[[y_col]] <- roles$y_factor
  } else {
    d_Q[[y_col]] <- roles$y_num
  }

  d_g <- d
  if (roles$treatment_type == "binary") {
    d_g[[a_col]] <- roles$a_factor
  } else {
    d_g[[a_col]] <- roles$a_num
  }

  d_C <- d
  d_C[[".delta"]] <- factor(roles$delta, levels = c("0", "1"))

  # Attach hardhat importance weights column if weights are present
  if (roles$weights_kind != "none") {
    w_obj <- if (roles$weights_kind == "frequency") {
      hardhat::frequency_weights(roles$weights)
    } else {
      hardhat::importance_weights(roles$weights)
    }
    d_Q[[".tmle_wt"]] <- w_obj
    d_g[[".tmle_wt"]] <- w_obj
    d_C[[".tmle_wt"]] <- w_obj
  }

  list(Q = d_Q, g = d_g, C = d_C)
}

#' Build or Adapt Role-Safe Causal Recipes for Q, g, and C
#' @keywords internal
split_causal_recipes <- function(base_recipe, roles, nuisance_data) {
  rec_Q <- adapt_recipe_for_role(base_recipe, role_type = "Q", roles = roles, template_data = nuisance_data$Q)
  rec_g <- adapt_recipe_for_role(base_recipe, role_type = "g", roles = roles, template_data = nuisance_data$g)
  rec_C <- if (roles$has_censoring) {
    adapt_recipe_for_role(base_recipe, role_type = "C", roles = roles, template_data = nuisance_data$C)
  } else {
    NULL
  }

  list(Q = rec_Q, g = rec_g, C = rec_C)
}

#' Adapt a User Recipe (or Build Default Recipe) for a Specific Causal Nuisance Role
#'
#' Guarantees zero leakage across causal roles:
#' - `Q`: outcome = `outcome_col`, predictors = `c(treatment_col, covariate_cols)`.
#' - `g`: outcome = `treatment_col`, predictors = `covariate_cols` (`treatment_col` and `outcome_col` strictly excluded from predictors).
#' - `C`: outcome = `.delta`, predictors = `c(treatment_col, covariate_cols)`.
#' - For LMTP stacked density ratio (`g_dr`): outcome = `.lambda`, predictors = `c(treatment_col, covariate_cols)`.
#'
#' @keywords internal
adapt_recipe_for_role <- function(base_recipe, role_type, roles, template_data) {
  y_col <- roles$outcome_col
  a_col <- roles$treatment_col
  w_cols <- roles$covariate_cols

  target_outcome <- switch(
    role_type,
    Q = y_col,
    g = a_col,
    C = ".delta",
    g_dr = ".lambda"
  )

  target_predictors <- switch(
    role_type,
    Q = c(a_col, w_cols),
    g = w_cols,
    C = c(a_col, w_cols),
    g_dr = c(a_col, w_cols)
  )

  has_tmle_wt <- ".tmle_wt" %in% names(template_data)

  if (is.null(base_recipe)) {
    # Construct clean formula recipe using only target outcome and target predictors
    # Include .tmle_wt if present (hardhat weights are auto-assigned role = "case_weights")
    terms_rhs <- if (has_tmle_wt) c(target_predictors, ".tmle_wt") else target_predictors
    fml <- stats::reformulate(termlabels = terms_rhs, response = target_outcome)
    rec <- recipes::recipe(fml, data = utils::head(template_data, 20L))

    # Auto-impute covariates if NAs exist in baseline predictors
    has_cov_na <- any(vapply(roles$data[w_cols], anyNA, logical(1L)))
    if (has_cov_na) {
      rec <- rec |>
        recipes::step_impute_median(recipes::all_numeric_predictors()) |>
        recipes::step_impute_mode(recipes::all_nominal_predictors())
    }
    rec <- rec |>
      recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = FALSE) |>
      recipes::step_zv(recipes::all_predictors())
    return(rec)
  }

  # Clone and retarget user-supplied recipe while preserving feature engineering steps
  rec <- base_recipe
  # Reset template slice with clean target types from template_data
  templ_slice <- utils::head(template_data, min(nrow(template_data), 20L))
  # Ensure all columns originally in rec$var_info exist in templ_slice
  missing_in_templ <- setdiff(rec$var_info$variable, names(templ_slice))
  if (length(missing_in_templ) > 0L && length(setdiff(missing_in_templ, c(".delta", ".lambda"))) > 0L) {
    # Keep existing template columns if any auxiliary columns were absent
    for (mc in missing_in_templ) {
      if (mc %in% names(rec$template)) {
        templ_slice[[mc]] <- utils::head(rec$template[[mc]], nrow(templ_slice))
      }
    }
  }
  rec$template <- templ_slice

  # Rebuild var_info and term_info cleanly for the target causal role
  all_vars <- unique(c(target_outcome, target_predictors, if (has_tmle_wt) ".tmle_wt", rec$var_info$variable))
  # Only retain variables present in templ_slice
  all_vars <- intersect(all_vars, names(templ_slice))

  get_var_type <- function(v) {
    x <- templ_slice[[v]]
    if (hardhat::is_importance_weights(x) || hardhat::is_frequency_weights(x)) {
      list("case_weights")
    } else if (is.numeric(x) || is.integer(x)) {
      list(c("double", "numeric"))
    } else if (is.factor(x) || is.character(x)) {
      list(c("nominal", "unordered", "factor"))
    } else if (is.logical(x)) {
      list(c("logical", "binary"))
    } else {
      list("other")
    }
  }

  assign_role <- function(v) {
    if (identical(v, target_outcome)) {
      return("outcome")
    }
    if (identical(v, ".tmle_wt")) {
      return("case_weights")
    }
    if (v %in% target_predictors) {
      return("predictor")
    }
    if (!is.null(roles$id_col) && identical(v, roles$id_col)) {
      return("id")
    }
    "ignored"
  }

  vi <- tibble::tibble(
    variable = all_vars,
    type = unlist(lapply(all_vars, get_var_type), recursive = FALSE),
    role = vapply(all_vars, assign_role, character(1L)),
    source = "original"
  )

  rec$var_info <- vi
  rec$term_info <- vi
  rec$ptype <- templ_slice[0L, all_vars, drop = FALSE]
  rec$orig_lvls <- lapply(templ_slice[all_vars], function(x) {
    if (is.factor(x)) list(values = levels(x), ordered = is.ordered(x)) else list(values = NA, ordered = NA)
  })
  names(rec$orig_lvls) <- all_vars

  # Ensure ignored/id/case_weights roles do not fail bake checks when omitted
  extra_roles <- setdiff(unique(vi$role), c("outcome", "predictor"))
  for (er in extra_roles) {
    rec <- recipes::update_role_requirements(rec, role = er, bake = FALSE)
  }

  # Strip id/ignored columns cleanly at the end of preprocessing
  if (any(vi$role %in% c("id", "ignored"))) {
    rec <- recipes::step_rm(rec, recipes::has_role(c("id", "ignored")))
  }

  rec
}

#' Build Smart Default Resamples for Cross-Fitting
#' @keywords internal
build_default_resamples <- function(roles, v_folds = 10L, strata = TRUE) {
  d <- roles$data
  n <- nrow(d)
  v_eff <- min(as.integer(v_folds), n)

  if (!is.null(roles$id_col)) {
    n_clusters <- length(unique(roles$id_vec))
    v_grp <- min(v_eff, n_clusters)
    if (v_grp < 2L) {
      cli::cli_abort("At least 2 clusters are required for grouped cross-validation.")
    }
    return(rsample::group_vfold_cv(d, group = !!rlang::sym(roles$id_col), v = v_grp))
  }

  if (!isTRUE(strata)) {
    return(rsample::vfold_cv(d, v = v_eff))
  }

  if (roles$treatment_type == "binary") {
    # Double-stratify on (A, Y) when binary outcome and all strata have >= v_eff rows
    if (roles$outcome_type == "binary") {
      strata_vec <- ifelse(
        roles$delta == 1L,
        paste(roles$a_num, roles$y_num, sep = "_"),
        paste(roles$a_num, "cens", sep = "_")
      )
      min_cell <- min(table(strata_vec))
      if (min_cell >= v_eff) {
        d_tmp <- d
        d_tmp[[".tmle_strata"]] <- strata_vec
        res <- rsample::vfold_cv(d_tmp, v = v_eff, strata = ".tmle_strata")
        res$splits <- lapply(res$splits, function(sp) {
          sp$data[[".tmle_strata"]] <- NULL
          sp
        })
        return(res)
      }
    }
    # Otherwise stratify on treatment A if both arms have >= v_eff rows
    if (min(table(roles$a_num)) >= v_eff) {
      return(rsample::vfold_cv(d, v = v_eff, strata = !!rlang::sym(roles$treatment_col)))
    }
  }

  # Fallback unstratified or quantile-binned continuous treatment
  rsample::vfold_cv(d, v = v_eff)
}
