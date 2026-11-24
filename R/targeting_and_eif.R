#' Compute Adaptive or Explicit Propensity Score Truncation Bounds
#' @noRd
compute_adaptive_g_bounds <- function(n_eff, g_bounds = "adaptive") {
  if (identical(g_bounds, "adaptive")) {
    n_safe <- max(as.numeric(n_eff), 10)
    b <- 5 / (sqrt(n_safe) * log(n_safe))
    b <- min(0.15, max(0.01, b))
    return(c(b, 1 - b))
  }
  if (is.numeric(g_bounds) && length(g_bounds) == 1L) {
    return(c(g_bounds, 1 - g_bounds))
  }
  as.numeric(g_bounds)
}

#' Truncate Propensity Scores and Record Diagnostics
#' @noRd
truncate_propensity <- function(g_hat, bounds) {
  g_trunc <- pmin(pmax(g_hat, bounds[1L]), bounds[2L])
  is_trunc <- (g_hat < bounds[1L]) | (g_hat > bounds[2L])
  list(
    g_trunc = g_trunc,
    is_truncated = is_trunc,
    pct_truncated = mean(is_trunc),
    g_min = min(g_hat),
    g_max = max(g_hat),
    bounds = bounds
  )
}

#' Scale Outcome to (0, 1) and Invert Scaling
#' @noRd
scale_outcome_01 <- function(y, bounds) {
  y_min <- bounds[1L]
  y_max <- bounds[2L]
  span <- y_max - y_min
  if (span <= 0) {
    span <- 1
  }
  (y - y_min) / span
}

unscale_outcome_01 <- function(y_scaled, bounds) {
  y_min <- bounds[1L]
  y_max <- bounds[2L]
  span <- y_max - y_min
  y_min + y_scaled * span
}

clamp_prob <- function(p, eps = 1e-5) {
  pmin(pmax(p, eps), 1 - eps)
}

#' Fit Weighted One-Parameter Logistic Fluctuation Submodel
#' @noRd
fit_logistic_fluctuation_1d <- function(y01, q01_init, clever_h, weights) {
  q_clamped <- clamp_prob(q01_init, eps = 1e-5)
  offset_eta <- stats::qlogis(q_clamped)
  valid <- is.finite(y01) & is.finite(clever_h) & is.finite(offset_eta) & (weights > 0)
  if (!any(valid) || all(clever_h[valid] == 0)) {
    return(0)
  }
  fit_glm <- suppressWarnings(
    stats::glm(
      y_sub ~ -1 + h_sub + offset(off_sub),
      data = data.frame(
        y_sub = y01[valid],
        h_sub = clever_h[valid],
        off_sub = offset_eta[valid]
      ),
      family = stats::quasibinomial(),
      weights = weights[valid],
      control = stats::glm.control(maxit = 100L, epsilon = 1e-10)
    )
  )
  eps <- as.numeric(stats::coef(fit_glm)[1L])
  if (!is.finite(eps)) 0 else eps
}

#' Fit Weighted Intercept-Only Logistic Fluctuation Submodel (Submodel with Weights)
#' @noRd
fit_weighted_logistic_intercept <- function(y01, q01_init, fluc_weights) {
  q_clamped <- clamp_prob(q01_init, eps = 1e-5)
  offset_eta <- stats::qlogis(q_clamped)
  valid <- is.finite(y01) & is.finite(offset_eta) & is.finite(fluc_weights) & (fluc_weights > 0)
  if (!any(valid)) {
    return(0)
  }
  fit_glm <- suppressWarnings(
    stats::glm(
      y_sub ~ 1 + offset(off_sub),
      data = data.frame(
        y_sub = y01[valid],
        off_sub = offset_eta[valid]
      ),
      family = stats::quasibinomial(),
      weights = fluc_weights[valid],
      control = stats::glm.control(maxit = 100L, epsilon = 1e-10)
    )
  )
  eps <- as.numeric(stats::coef(fit_glm)[1L])
  if (!is.finite(eps)) 0 else eps
}

#' Targeted Fluctuation Step for Binary Point-Treatment or Policy Estimands
#' @noRd
run_tmle_targeting <- function(y_num,
                               a_num,
                               delta,
                               q0_init,
                               q1_init,
                               qa_init,
                               g_hat,
                               c1_hat,
                               c0_hat,
                               weights,
                               estimand,
                               outcome_type,
                               y_bounds,
                               g_bounds_arg,
                               fluctuation = "logistic",
                               cluster_id = NULL,
                               conf.level = 0.95,
                               simultaneous = FALSE,
                               e_value = FALSE) {
  n <- length(y_num)
  n_eff <- if (!is.null(cluster_id)) length(unique(cluster_id)) else n
  bounds <- compute_adaptive_g_bounds(n_eff = n_eff, g_bounds = g_bounds_arg)

  # Normalize weights so mean(w) = 1
  w <- weights / mean(weights)
  obs_mask <- (delta == 1L)

  if (estimand$kind == "policy") {
    return(run_policy_targeting(
      y_num = y_num,
      delta = delta,
      qa_init = qa_init,
      q_shift_init = q1_init,
      p_lambda1 = g_hat,
      c_obs_hat = c1_hat,
      weights = w,
      policy = estimand$policies[[1L]],
      outcome_type = outcome_type,
      y_bounds = y_bounds,
      g_bounds = bounds,
      fluctuation = fluctuation,
      cluster_id = cluster_id,
      conf.level = conf.level
    ))
  }

  # Binary point-treatment targeting
  if ("or" %in% estimand$names && outcome_type != "binary") {
    cli::cli_abort(c(
      "Odds Ratio ({.val or}) requires a binary outcome in {{0, 1}}.",
      "i" = "For continuous outcomes, use {.val ate}, {.val att}, {.val atc}, {.val ey1}, or {.val ey0}."
    ))
  }

  g_info <- truncate_propensity(g_hat, bounds)
  g_t <- g_info$g_trunc

  # Truncate censoring probabilities conservatively away from 0
  c1_t <- pmax(c1_hat, 0.025)
  c0_t <- pmax(c0_hat, 0.025)

  # Scale outcomes and initial predictions into (0, 1)
  if (outcome_type == "continuous") {
    y_01 <- scale_outcome_01(y_num, y_bounds)
    q0_01 <- clamp_prob(scale_outcome_01(q0_init, y_bounds), 1e-4)
    q1_01 <- clamp_prob(scale_outcome_01(q1_init, y_bounds), 1e-4)
  } else {
    y_01 <- y_num
    q0_01 <- clamp_prob(q0_init, 1e-4)
    q1_01 <- clamp_prob(q1_init, 1e-4)
  }
  qa_01 <- ifelse(a_num == 1L, q1_01, q0_01)

  # Clever covariates for joint (EY1, EY0) targeting
  h1_obs <- ifelse(a_num == 1L & obs_mask, 1 / (g_t * c1_t), 0)
  h0_obs <- ifelse(a_num == 0L & obs_mask, 1 / ((1 - g_t) * c0_t), 0)
  clever_h <- ifelse(a_num == 1L, 1 / (g_t * c1_t), -1 / ((1 - g_t) * c0_t))

  if (identical(fluctuation, "none")) {
    q1_star_01 <- q1_01
    q0_star_01 <- q0_01
    eps_vec <- c(eps1 = 0, eps0 = 0)
  } else if (identical(fluctuation, "logistic")) {
    # Joint 2-arm weighted logistic fluctuation:
    # Solves both sum w_i * H1_i * (Y_i - Q*_1(W_i)) = 0 and sum w_i * H0_i * (Y_i - Q*_0(W_i)) = 0!
    eps1 <- fit_weighted_logistic_intercept(
      y01 = y_01[a_num == 1L & obs_mask],
      q01_init = q1_01[a_num == 1L & obs_mask],
      fluc_weights = (w * h1_obs)[a_num == 1L & obs_mask]
    )
    eps0 <- fit_weighted_logistic_intercept(
      y01 = y_01[a_num == 0L & obs_mask],
      q01_init = q0_01[a_num == 0L & obs_mask],
      fluc_weights = (w * h0_obs)[a_num == 0L & obs_mask]
    )
    q1_star_01 <- stats::plogis(stats::qlogis(q1_01) + eps1)
    q0_star_01 <- stats::plogis(stats::qlogis(q0_01) + eps0)
    eps_vec <- c(eps1 = eps1, eps0 = eps0)
  } else if (identical(fluctuation, "logistic_covariate")) {
    # Classical unweighted Clever-Covariate logistic submodel (matching tmle::tmle target.gwt = FALSE)
    h1_cov <- ifelse(a_num == 1L, 1 / (g_t * c1_t), 0)
    h0_cov <- ifelse(a_num == 0L, 1 / ((1 - g_t) * c0_t), 0)
    qa_clamped <- clamp_prob(qa_01, eps = 1e-5)
    off_cov <- stats::qlogis(qa_clamped)
    valid_cov <- is.finite(y_01) & obs_mask & is.finite(off_cov) & (w > 0)
    if (!any(valid_cov)) {
      eps1 <- 0
      eps0 <- 0
    } else {
      fit_cov <- suppressWarnings(
        stats::glm(
          y_sub ~ -1 + h1_sub + h0_sub + offset(off_sub),
          data = data.frame(
            y_sub = y_01[valid_cov],
            h1_sub = h1_cov[valid_cov],
            h0_sub = h0_cov[valid_cov],
            off_sub = off_cov[valid_cov]
          ),
          family = stats::quasibinomial(),
          weights = w[valid_cov],
          control = stats::glm.control(maxit = 100L, epsilon = 1e-10)
        )
      )
      coefs_cov <- stats::coef(fit_cov)
      eps1 <- if (is.finite(coefs_cov["h1_sub"])) as.numeric(coefs_cov["h1_sub"]) else 0
      eps0 <- if (is.finite(coefs_cov["h0_sub"])) as.numeric(coefs_cov["h0_sub"]) else 0
    }
    q1_star_01 <- stats::plogis(stats::qlogis(clamp_prob(q1_01, 1e-5)) + eps1 / (g_t * c1_t))
    q0_star_01 <- stats::plogis(stats::qlogis(clamp_prob(q0_01, 1e-5)) + eps0 / ((1 - g_t) * c0_t))
    eps_vec <- c(eps1 = eps1, eps0 = eps0)
  } else {
    # Weighted linear fluctuation
    res1 <- y_01[a_num == 1L & obs_mask] - q1_01[a_num == 1L & obs_mask]
    w1 <- (w * h1_obs)[a_num == 1L & obs_mask]
    eps1 <- stats::weighted.mean(res1, w = w1)
    res0 <- y_01[a_num == 0L & obs_mask] - q0_01[a_num == 0L & obs_mask]
    w0 <- (w * h0_obs)[a_num == 0L & obs_mask]
    eps0 <- stats::weighted.mean(res0, w = w0)
    q1_star_01 <- q1_01 + eps1
    q0_star_01 <- q0_01 + eps0
    eps_vec <- c(eps1 = eps1, eps0 = eps0)
  }

  # Rescale back to original outcome scale
  if (outcome_type == "continuous") {
    q1_star <- unscale_outcome_01(q1_star_01, y_bounds)
    q0_star <- unscale_outcome_01(q0_star_01, y_bounds)
  } else {
    q1_star <- q1_star_01
    q0_star <- q0_star_01
  }
  qa_star <- ifelse(a_num == 1L, q1_star, q0_star)

  # Separate ATT/ATC targeting when requested
  att_res <- NULL
  if ("att" %in% estimand$names) {
    att_res <- target_att_or_atc(
      y_01 = y_01, a_num = a_num, obs_mask = obs_mask,
      q1_01 = q1_01, q0_01 = q0_01, g_t = g_t, c1_t = c1_t, c0_t = c0_t,
      w = w, y_bounds = y_bounds, outcome_type = outcome_type,
      fluctuation = fluctuation, target = "att"
    )
  }
  atc_res <- NULL
  if ("atc" %in% estimand$names) {
    atc_res <- target_att_or_atc(
      y_01 = y_01, a_num = a_num, obs_mask = obs_mask,
      q1_01 = q1_01, q0_01 = q0_01, g_t = g_t, c1_t = c1_t, c0_t = c0_t,
      w = w, y_bounds = y_bounds, outcome_type = outcome_type,
      fluctuation = fluctuation, target = "atc"
    )
  }

  # Compute Joint EIF matrix for EY1 and EY0
  psi1 <- stats::weighted.mean(q1_star, w = w)
  psi0 <- stats::weighted.mean(q0_star, w = w)

  y_safe <- ifelse(obs_mask, y_num, 0)
  d1 <- w * (ifelse(a_num == 1L & obs_mask, (y_safe - q1_star) / (g_t * c1_t), 0) + q1_star - psi1)
  d0 <- w * (ifelse(a_num == 0L & obs_mask, (y_safe - q0_star) / ((1 - g_t) * c0_t), 0) + q0_star - psi0)
  d_ate <- d1 - d0

  # Plug-in initial ATE before targeting
  plugin_ate <- stats::weighted.mean(q1_init - q0_init, w = w)

  # Build inference table across all requested estimands
  est_list <- lapply(estimand$names, function(est) {
    compute_single_estimand_inference(
      est_name = est,
      psi1 = psi1, psi0 = psi0,
      d1 = d1, d0 = d0,
      att_res = att_res, atc_res = atc_res,
      y_safe = y_safe, a_num = a_num, obs_mask = obs_mask,
      w = w, cluster_id = cluster_id,
      conf.level = conf.level,
      e_value = e_value,
      outcome_type = outcome_type
    )
  })
  estimates_tbl <- dplyr::bind_rows(est_list)

  # Simultaneous uniform confidence bands via Rademacher wild multiplier bootstrap
  if (isTRUE(simultaneous) && nrow(estimates_tbl) > 1L) {
    estimates_tbl <- add_simultaneous_bands(
      estimates_tbl = estimates_tbl,
      eif_list = lapply(est_list, function(r) attr(r, "eif_scale")),
      cluster_id = cluster_id,
      conf.level = conf.level
    )
  }

  # Primary EIF for augment() and glance()
  primary_est <- estimand$names[1L]
  primary_eif <- switch(
    primary_est,
    ey1 = d1,
    ey0 = d0,
    att = if (!is.null(att_res)) att_res$eif else d_ate,
    atc = if (!is.null(atc_res)) atc_res$eif else d_ate,
    rr = d1 / psi1 - d0 / psi0,
    or = d1 / (psi1 * (1 - psi1)) - d0 / (psi0 * (1 - psi0)),
    d_ate
  )

  # Effective sample sizes (Kish's ESS on IPW weights among observed)
  ipw_1 <- w[a_num == 1L & obs_mask] / (g_t[a_num == 1L & obs_mask] * c1_t[a_num == 1L & obs_mask])
  ipw_0 <- w[a_num == 0L & obs_mask] / ((1 - g_t[a_num == 0L & obs_mask]) * c0_t[a_num == 0L & obs_mask])
  ess_treated <- if (length(ipw_1) > 0L) (sum(ipw_1)^2) / sum(ipw_1^2) else 0
  ess_control <- if (length(ipw_0) > 0L) (sum(ipw_0)^2) / sum(ipw_0^2) else 0

  eif_sq <- primary_eif^2
  eif_leverage <- if (sum(eif_sq) > 0) eif_sq / sum(eif_sq) else rep(0, n)

  # Representative scalar epsilon, targeting outputs, and plug-in shift for primary estimand
  if (primary_est == "att" && !is.null(att_res)) {
    q1_out <- att_res$q1_star
    q0_out <- att_res$q0_star
    qa_out <- ifelse(a_num == 1L, q1_out, q0_out)
    clever_out <- att_res$clever_h
    eps_arms_out <- c(eps1 = att_res$eps1, eps0 = att_res$eps0)
    eps_scalar <- mean(abs(eps_arms_out)) * sign(eps_arms_out[1L] - eps_arms_out[2L] + 1e-15)
    primary_shift <- att_res$estimate - att_res$plugin_estimate
  } else if (primary_est == "atc" && !is.null(atc_res)) {
    q1_out <- atc_res$q1_star
    q0_out <- atc_res$q0_star
    qa_out <- ifelse(a_num == 1L, q1_out, q0_out)
    clever_out <- atc_res$clever_h
    eps_arms_out <- c(eps1 = atc_res$eps1, eps0 = atc_res$eps0)
    eps_scalar <- mean(abs(eps_arms_out)) * sign(eps_arms_out[1L] - eps_arms_out[2L] + 1e-15)
    primary_shift <- atc_res$estimate - atc_res$plugin_estimate
  } else {
    q1_out <- q1_star
    q0_out <- q0_star
    qa_out <- qa_star
    clever_out <- clever_h
    eps_arms_out <- eps_vec
    eps_scalar <- if (primary_est == "ey1") {
      eps_vec["eps1"]
    } else if (primary_est == "ey0") {
      eps_vec["eps0"]
    } else {
      mean(abs(eps_vec)) * sign(eps_vec[1L] - eps_vec[2L] + 1e-15)
    }
    primary_shift <- switch(
      primary_est,
      ey1 = psi1 - stats::weighted.mean(q1_init, w = w),
      ey0 = psi0 - stats::weighted.mean(q0_init, w = w),
      rr = (psi1 / psi0) - (stats::weighted.mean(q1_init, w = w) / stats::weighted.mean(q0_init, w = w)),
      or = {
        p1_0 <- clamp_prob(stats::weighted.mean(q1_init, w = w), 1e-6)
        p0_0 <- clamp_prob(stats::weighted.mean(q0_init, w = w), 1e-6)
        p1_s <- clamp_prob(psi1, 1e-6)
        p0_s <- clamp_prob(psi0, 1e-6)
        (p1_s / (1 - p1_s)) / (p0_s / (1 - p0_s)) - (p1_0 / (1 - p1_0)) / (p0_0 / (1 - p0_0))
      },
      (psi1 - psi0) - plugin_ate
    )
  }

  list(
    estimates = estimates_tbl,
    q0_init = q0_init,
    q1_init = q1_init,
    qa_init = ifelse(a_num == 1L, q1_init, q0_init),
    q0_star = q0_out,
    q1_star = q1_out,
    qa_star = qa_out,
    g_hat = g_hat,
    g_trunc = g_t,
    c1_hat = c1_hat,
    c0_hat = c0_hat,
    clever_h = clever_out,
    eif = primary_eif,
    eif_1 = d1,
    eif_0 = d0,
    eif_ate = d_ate,
    eif_leverage = eif_leverage,
    epsilon = as.numeric(eps_scalar),
    epsilon_arms = eps_arms_out,
    mean_eif = mean(primary_eif),
    plugin_vs_tmle_shift = as.numeric(primary_shift),
    g_info = g_info,
    ess_treated = ess_treated,
    ess_control = ess_control,
    conf_level = conf.level
  )
}

#' Target ATT or ATC Parameter via TMLE Fluctuation
#' @noRd
target_att_or_atc <- function(y_01, a_num, obs_mask, q1_01, q0_01, g_t, c1_t, c0_t,
                              w, y_bounds, outcome_type, fluctuation, target = c("att", "atc")) {
  target <- match.arg(target)
  p_target <- if (target == "att") mean(w * a_num) else mean(w * (1 - a_num))
  p_target <- max(p_target, 1e-8)

  if (target == "att") {
    # Clever covariate for ATT:
    # H(1, W) = 1 / (p_A * C(1, W)), H(0, W) = -g(W) / ((1 - g(W)) * p_A * C(0, W))
    h1 <- ifelse(a_num == 1L & obs_mask, 1 / (p_target * c1_t), 0)
    h0 <- ifelse(a_num == 0L & obs_mask, -g_t / ((1 - g_t) * p_target * c0_t), 0)
    h_obs <- h1 + h0
    # Fluctuate arm 1 and arm 0 separately on log-likelihood to solve ATT EIF = 0
    if (fluctuation == "logistic") {
      eps1 <- fit_weighted_logistic_intercept(
        y_01[a_num == 1L & obs_mask], q1_01[a_num == 1L & obs_mask],
        (w / (p_target * c1_t))[a_num == 1L & obs_mask]
      )
      eps0 <- fit_weighted_logistic_intercept(
        y_01[a_num == 0L & obs_mask], q0_01[a_num == 0L & obs_mask],
        (w * g_t / ((1 - g_t) * p_target * c0_t))[a_num == 0L & obs_mask]
      )
      q1_s01 <- stats::plogis(stats::qlogis(q1_01) + eps1)
      q0_s01 <- stats::plogis(stats::qlogis(q0_01) + eps0)
    } else if (fluctuation == "logistic_covariate") {
      h1_vec <- 1 / (p_target * c1_t)
      h0_vec <- -g_t / ((1 - g_t) * p_target * c0_t)
      eps1 <- fit_logistic_fluctuation_1d(
        y01 = y_01[a_num == 1L & obs_mask],
        q01_init = q1_01[a_num == 1L & obs_mask],
        clever_h = h1_vec[a_num == 1L & obs_mask],
        weights = w[a_num == 1L & obs_mask]
      )
      eps0 <- fit_logistic_fluctuation_1d(
        y01 = y_01[a_num == 0L & obs_mask],
        q01_init = q0_01[a_num == 0L & obs_mask],
        clever_h = h0_vec[a_num == 0L & obs_mask],
        weights = w[a_num == 0L & obs_mask]
      )
      q1_s01 <- stats::plogis(stats::qlogis(clamp_prob(q1_01, 1e-5)) + eps1 * h1_vec)
      q0_s01 <- stats::plogis(stats::qlogis(clamp_prob(q0_01, 1e-5)) + eps0 * h0_vec)
    } else if (fluctuation == "linear") {
      w1 <- (w / (p_target * c1_t))[a_num == 1L & obs_mask]
      eps1 <- stats::weighted.mean((y_01 - q1_01)[a_num == 1L & obs_mask], w = w1)
      w0 <- (w * g_t / ((1 - g_t) * p_target * c0_t))[a_num == 0L & obs_mask]
      eps0 <- stats::weighted.mean((y_01 - q0_01)[a_num == 0L & obs_mask], w = w0)
      q1_s01 <- q1_01 + eps1
      q0_s01 <- q0_01 + eps0
    } else {
      eps1 <- 0
      eps0 <- 0
      q1_s01 <- q1_01
      q0_s01 <- q0_01
    }

    q1_0_orig <- if (outcome_type == "continuous") unscale_outcome_01(q1_01, y_bounds) else q1_01
    q0_0_orig <- if (outcome_type == "continuous") unscale_outcome_01(q0_01, y_bounds) else q0_01
    plugin_est <- stats::weighted.mean((q1_0_orig - q0_0_orig)[a_num == 1L], w = w[a_num == 1L])

    q1_s <- if (outcome_type == "continuous") unscale_outcome_01(q1_s01, y_bounds) else q1_s01
    q0_s <- if (outcome_type == "continuous") unscale_outcome_01(q0_s01, y_bounds) else q0_s01
    y_orig <- if (outcome_type == "continuous") unscale_outcome_01(y_01, y_bounds) else y_01

    psi <- stats::weighted.mean((q1_s - q0_s)[a_num == 1L], w = w[a_num == 1L])
    eif <- w * (
      h_obs * (ifelse(obs_mask, y_orig, 0) - ifelse(a_num == 1L, q1_s, q0_s)) +
        (a_num / p_target) * (q1_s - q0_s - psi)
    )
    return(list(
      estimate = psi, eif = eif,
      eps1 = eps1, eps0 = eps0,
      q1_star = q1_s, q0_star = q0_s,
      clever_h = h_obs, plugin_estimate = plugin_est
    ))
  } else {
    h1 <- ifelse(a_num == 1L & obs_mask, (1 - g_t) / (g_t * p_target * c1_t), 0)
    h0 <- ifelse(a_num == 0L & obs_mask, -1 / (p_target * c0_t), 0)
    h_obs <- h1 + h0
    if (fluctuation == "logistic") {
      eps1 <- fit_weighted_logistic_intercept(
        y_01[a_num == 1L & obs_mask], q1_01[a_num == 1L & obs_mask],
        (w * (1 - g_t) / (g_t * p_target * c1_t))[a_num == 1L & obs_mask]
      )
      eps0 <- fit_weighted_logistic_intercept(
        y_01[a_num == 0L & obs_mask], q0_01[a_num == 0L & obs_mask],
        (w / (p_target * c0_t))[a_num == 0L & obs_mask]
      )
      q1_s01 <- stats::plogis(stats::qlogis(q1_01) + eps1)
      q0_s01 <- stats::plogis(stats::qlogis(q0_01) + eps0)
    } else if (fluctuation == "logistic_covariate") {
      h1_vec <- (1 - g_t) / (g_t * p_target * c1_t)
      h0_vec <- -1 / (p_target * c0_t)
      eps1 <- fit_logistic_fluctuation_1d(
        y01 = y_01[a_num == 1L & obs_mask],
        q01_init = q1_01[a_num == 1L & obs_mask],
        clever_h = h1_vec[a_num == 1L & obs_mask],
        weights = w[a_num == 1L & obs_mask]
      )
      eps0 <- fit_logistic_fluctuation_1d(
        y01 = y_01[a_num == 0L & obs_mask],
        q01_init = q0_01[a_num == 0L & obs_mask],
        clever_h = h0_vec[a_num == 0L & obs_mask],
        weights = w[a_num == 0L & obs_mask]
      )
      q1_s01 <- stats::plogis(stats::qlogis(clamp_prob(q1_01, 1e-5)) + eps1 * h1_vec)
      q0_s01 <- stats::plogis(stats::qlogis(clamp_prob(q0_01, 1e-5)) + eps0 * h0_vec)
    } else if (fluctuation == "linear") {
      w1 <- (w * (1 - g_t) / (g_t * p_target * c1_t))[a_num == 1L & obs_mask]
      eps1 <- stats::weighted.mean((y_01 - q1_01)[a_num == 1L & obs_mask], w = w1)
      w0 <- (w / (p_target * c0_t))[a_num == 0L & obs_mask]
      eps0 <- stats::weighted.mean((y_01 - q0_01)[a_num == 0L & obs_mask], w = w0)
      q1_s01 <- q1_01 + eps1
      q0_s01 <- q0_01 + eps0
    } else {
      eps1 <- 0
      eps0 <- 0
      q1_s01 <- q1_01
      q0_s01 <- q0_01
    }

    q1_0_orig <- if (outcome_type == "continuous") unscale_outcome_01(q1_01, y_bounds) else q1_01
    q0_0_orig <- if (outcome_type == "continuous") unscale_outcome_01(q0_01, y_bounds) else q0_01
    plugin_est <- stats::weighted.mean((q1_0_orig - q0_0_orig)[a_num == 0L], w = w[a_num == 0L])

    q1_s <- if (outcome_type == "continuous") unscale_outcome_01(q1_s01, y_bounds) else q1_s01
    q0_s <- if (outcome_type == "continuous") unscale_outcome_01(q0_s01, y_bounds) else q0_s01
    y_orig <- if (outcome_type == "continuous") unscale_outcome_01(y_01, y_bounds) else y_01

    psi <- stats::weighted.mean((q1_s - q0_s)[a_num == 0L], w = w[a_num == 0L])
    eif <- w * (
      h_obs * (ifelse(obs_mask, y_orig, 0) - ifelse(a_num == 1L, q1_s, q0_s)) +
        ((1 - a_num) / p_target) * (q1_s - q0_s - psi)
    )
    return(list(
      estimate = psi, eif = eif,
      eps1 = eps1, eps0 = eps0,
      q1_star = q1_s, q0_star = q0_s,
      clever_h = h_obs, plugin_estimate = plugin_est
    ))
  }
}

#' Targeted Fluctuation Step for Modified Treatment Policies (LMTP)
#' @noRd
run_policy_targeting <- function(y_num, delta, qa_init, q_shift_init, p_lambda1,
                                 c_obs_hat, weights, policy, outcome_type,
                                 y_bounds, g_bounds, fluctuation, cluster_id, conf.level) {
  n <- length(y_num)
  obs_mask <- (delta == 1L)
  c_t <- pmax(c_obs_hat, 0.025)

  # Density ratio r(A, W) = P(Lambda=1 | A, W) / P(Lambda=0 | A, W)
  p_clamped <- clamp_prob(p_lambda1, 1e-4)
  r_raw <- p_clamped / (1 - p_clamped)
  # Convert g_bounds to ratio cap: ratio cap corresponds to 1 / g_bounds[1]
  max_r <- 1 / g_bounds[1L]
  is_trunc <- r_raw > max_r
  r_trunc <- pmin(r_raw, max_r)
  clever_h <- ifelse(obs_mask, r_trunc / c_t, 0)

  if (outcome_type == "continuous") {
    y_01 <- scale_outcome_01(y_num, y_bounds)
    qa_01 <- clamp_prob(scale_outcome_01(qa_init, y_bounds), 1e-4)
    qs_01 <- clamp_prob(scale_outcome_01(q_shift_init, y_bounds), 1e-4)
  } else {
    y_01 <- y_num
    qa_01 <- clamp_prob(qa_init, 1e-4)
    qs_01 <- clamp_prob(q_shift_init, 1e-4)
  }

  if (fluctuation == "logistic") {
    eps <- fit_weighted_logistic_intercept(
      y01 = y_01[obs_mask],
      q01_init = qa_01[obs_mask],
      fluc_weights = (weights * clever_h)[obs_mask]
    )
    qa_star_01 <- stats::plogis(stats::qlogis(qa_01) + eps)
    qs_star_01 <- stats::plogis(stats::qlogis(qs_01) + eps)
  } else if (fluctuation == "logistic_covariate") {
    h_full <- r_trunc / c_t
    eps <- fit_logistic_fluctuation_1d(
      y01 = y_01[obs_mask],
      q01_init = qa_01[obs_mask],
      clever_h = h_full[obs_mask],
      weights = weights[obs_mask]
    )
    qa_star_01 <- stats::plogis(stats::qlogis(clamp_prob(qa_01, 1e-5)) + eps * h_full)
    qs_star_01 <- stats::plogis(stats::qlogis(clamp_prob(qs_01, 1e-5)) + eps * h_full)
  } else if (fluctuation == "linear") {
    eps <- stats::weighted.mean((y_01 - qa_01)[obs_mask], w = (weights * clever_h)[obs_mask])
    qa_star_01 <- qa_01 + eps
    qs_star_01 <- qs_01 + eps
  } else {
    eps <- 0
    qa_star_01 <- qa_01
    qs_star_01 <- qs_01
  }

  qa_star <- if (outcome_type == "continuous") unscale_outcome_01(qa_star_01, y_bounds) else qa_star_01
  qs_star <- if (outcome_type == "continuous") unscale_outcome_01(qs_star_01, y_bounds) else qs_star_01

  psi_policy <- stats::weighted.mean(qs_star, w = weights)
  y_safe <- ifelse(obs_mask, y_num, 0)
  eif_policy <- weights * (clever_h * (y_safe - qa_star) + qs_star - psi_policy)

  is_obs_contrast <- identical(policy$contrast, "observed")
  if (is_obs_contrast) {
    if (all(obs_mask)) {
      psi_obs <- stats::weighted.mean(y_num, w = weights)
      eif_obs <- weights * (y_num - psi_obs)
    } else {
      c_inv_obs <- ifelse(obs_mask, 1 / c_t, 0)
      aug_obs <- c_inv_obs * (y_safe - qa_star) + qa_star
      psi_obs <- stats::weighted.mean(aug_obs, w = weights)
      eif_obs <- weights * (aug_obs - psi_obs)
    }
    est_val <- psi_policy - psi_obs
    eif_out <- eif_policy - eif_obs
    plugin_est <- stats::weighted.mean(q_shift_init - qa_init, w = weights)
  } else {
    est_val <- psi_policy
    eif_out <- eif_policy
    plugin_est <- stats::weighted.mean(q_shift_init, w = weights)
  }

  se <- compute_eif_se(eif_out, cluster_id = cluster_id)
  z_crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  stat <- est_val / se
  pval <- 2 * stats::pnorm(-abs(stat))

  est_label <- if (!is.null(policy$base_label)) {
    if (is_obs_contrast) paste0(policy$base_label, " - observed") else policy$base_label
  } else {
    policy$label
  }

  estimates_tbl <- tibble::tibble(
    estimand = est_label,
    term = est_label,
    estimate = est_val,
    std.error = se,
    statistic = stat,
    p.value = pval,
    conf.low = est_val - z_crit * se,
    conf.high = est_val + z_crit * se
  )

  eif_sq <- eif_out^2
  eif_leverage <- if (sum(eif_sq) > 0) eif_sq / sum(eif_sq) else rep(0, n)

  ipw_p <- (weights * clever_h)[obs_mask]
  ess_p <- if (length(ipw_p) > 0L) (sum(ipw_p)^2) / sum(ipw_p^2) else 0

  list(
    estimates = estimates_tbl,
    q0_init = qa_init,
    q1_init = q_shift_init,
    qa_init = qa_init,
    q0_star = qa_star,
    q1_star = qs_star,
    qa_star = qa_star,
    g_hat = p_lambda1,
    g_trunc = r_trunc,
    c1_hat = c_obs_hat,
    c0_hat = c_obs_hat,
    clever_h = clever_h,
    eif = eif_out,
    eif_1 = eif_out,
    eif_0 = rep(0, n),
    eif_ate = eif_out,
    eif_leverage = eif_leverage,
    epsilon = eps,
    epsilon_arms = c(eps = eps),
    mean_eif = mean(eif_out),
    plugin_vs_tmle_shift = est_val - plugin_est,
    g_info = list(
      g_trunc = r_trunc,
      is_truncated = is_trunc,
      pct_truncated = mean(is_trunc),
      g_min = min(r_raw),
      g_max = max(r_raw),
      bounds = c(0, max_r)
    ),
    ess_treated = ess_p,
    ess_control = n,
    conf_level = conf.level
  )
}

#' Compute Cluster-Robust Influence Function Standard Error
#' @noRd
compute_eif_se <- function(eif_vec, cluster_id = NULL) {
  n <- length(eif_vec)
  if (!is.null(cluster_id)) {
    cl_fac <- droplevels(factor(cluster_id))
    cl_sums <- as.numeric(tapply(eif_vec, cl_fac, sum))
    cl_sums <- cl_sums[is.finite(cl_sums)]
    J <- length(cl_sums)
    if (J < 2L) {
      cli::cli_warn("Subgroup contains only {J} cluster; falling back to observation-level EIF standard error.")
      return(sqrt(max(stats::var(eif_vec) / n, 1e-16)))
    }
    # Variance of (1/n) * sum_c D_c = (J / (J - 1)) * (1 / n^2) * sum_c (D_c - mean(D_c))^2
    var_est <- (J / (J - 1L)) * sum((cl_sums - mean(cl_sums))^2) / (n^2)
    return(sqrt(max(var_est, 1e-16)))
  }
  sqrt(max(stats::var(eif_vec) / n, 1e-16))
}

#' Compute Single-Estimand Inference (Linear or Log-Scale Delta Method)
#' @noRd
compute_single_estimand_inference <- function(est_name,
                                              psi1,
                                              psi0,
                                              d1,
                                              d0,
                                              att_res,
                                              atc_res,
                                              y_safe,
                                              a_num,
                                              obs_mask,
                                              w,
                                              cluster_id,
                                              conf.level,
                                              e_value,
                                              outcome_type) {
  z_crit <- stats::qnorm(1 - (1 - conf.level) / 2)

  if (est_name == "ate") {
    est <- psi1 - psi0
    eif_s <- d1 - d0
    se <- compute_eif_se(eif_s, cluster_id)
    stat <- est / se
    pval <- 2 * stats::pnorm(-abs(stat))
    ci <- c(est - z_crit * se, est + z_crit * se)
    term <- "ATE"
  } else if (est_name == "ey1") {
    est <- psi1
    eif_s <- d1
    se <- compute_eif_se(eif_s, cluster_id)
    stat <- est / se
    pval <- 2 * stats::pnorm(-abs(stat))
    ci <- c(est - z_crit * se, est + z_crit * se)
    term <- "E[Y(1)]"
  } else if (est_name == "ey0") {
    est <- psi0
    eif_s <- d0
    se <- compute_eif_se(eif_s, cluster_id)
    stat <- est / se
    pval <- 2 * stats::pnorm(-abs(stat))
    ci <- c(est - z_crit * se, est + z_crit * se)
    term <- "E[Y(0)]"
  } else if (est_name == "att") {
    est <- att_res$estimate
    eif_s <- att_res$eif
    se <- compute_eif_se(eif_s, cluster_id)
    stat <- est / se
    pval <- 2 * stats::pnorm(-abs(stat))
    ci <- c(est - z_crit * se, est + z_crit * se)
    term <- "ATT"
  } else if (est_name == "atc") {
    est <- atc_res$estimate
    eif_s <- atc_res$eif
    se <- compute_eif_se(eif_s, cluster_id)
    stat <- est / se
    pval <- 2 * stats::pnorm(-abs(stat))
    ci <- c(est - z_crit * se, est + z_crit * se)
    term <- "ATC"
  } else if (est_name == "rr") {
    if (psi1 <= 0 || psi0 <= 0) {
      cli::cli_abort(c(
        "Risk/Rate Ratio ({.val rr}) requires strictly positive expected potential outcomes E[Y(1)] and E[Y(0)].",
        "x" = "Observed targeted means: E[Y(1)] = {round(psi1, 4)}, E[Y(0)] = {round(psi0, 4)}.",
        "i" = "A ratio such as E[Y(1)] / E[Y(0)] is undefined or flips sign when outcomes cross zero; use {.code estimand = \"ate\"} for zero-centered/signed outcomes."
      ))
    }
    p1 <- psi1
    p0 <- psi0
    log_rr <- log(p1 / p0)
    eif_s <- d1 / p1 - d0 / p0
    se_log <- compute_eif_se(eif_s, cluster_id)
    est <- exp(log_rr)
    se <- se_log # SE on log scale for ratio estimands (or Delta SE)
    stat <- log_rr / se_log
    pval <- 2 * stats::pnorm(-abs(stat))
    ci <- exp(c(log_rr - z_crit * se_log, log_rr + z_crit * se_log))
    term <- "RR"
  } else if (est_name == "or") {
    p1 <- clamp_prob(psi1, 1e-6)
    p0 <- clamp_prob(psi0, 1e-6)
    log_or <- log((p1 / (1 - p1)) / (p0 / (1 - p0)))
    eif_s <- d1 / (p1 * (1 - p1)) - d0 / (p0 * (1 - p0))
    se_log <- compute_eif_se(eif_s, cluster_id)
    est <- exp(log_or)
    se <- se_log
    stat <- log_or / se_log
    pval <- 2 * stats::pnorm(-abs(stat))
    ci <- exp(c(log_or - z_crit * se_log, log_or + z_crit * se_log))
    term <- "OR"
  }

  out <- tibble::tibble(
    estimand = est_name,
    term = term,
    estimate = est,
    std.error = se,
    statistic = stat,
    p.value = pval,
    conf.low = ci[1L],
    conf.high = ci[2L]
  )

  if (isTRUE(e_value)) {
    ev <- compute_vanderweele_evalue(est_name, est, se, ci[1L], ci[2L], y_safe[obs_mask])
    out$e_value <- ev$e_value
    out$e_value_ci <- ev$e_value_ci
  }

  attr(out, "eif_scale") <- eif_s
  attr(out, "is_ratio") <- est_name %in% c("rr", "or")
  out
}

#' Compute VanderWeele & Ding (2017) E-values for Unmeasured Confounding Sensitivity
#' @noRd
compute_vanderweele_evalue <- function(est_name, est, se, ci_low, ci_high, y_obs) {
  calc_ev <- function(rr) {
    if (!is.finite(rr) || rr <= 0) return(NA_real_)
    rr_adj <- if (rr < 1) 1 / rr else rr
    rr_adj + sqrt(rr_adj * max(rr_adj - 1, 0))
  }

  if (est_name %in% c("rr", "or")) {
    ev_pt <- calc_ev(est)
    ev_ci <- if (ci_low <= 1 && ci_high >= 1) 1 else if (ci_low > 1) calc_ev(ci_low) else calc_ev(ci_high)
    return(list(e_value = ev_pt, e_value_ci = ev_ci))
  }

  # Approximate RR conversion for continuous/risk-difference estimands (VanderWeele & Ding 2017)
  sd_y <- stats::sd(y_obs, na.rm = TRUE)
  if (!is.finite(sd_y) || sd_y <= 0) {
    return(list(e_value = NA_real_, e_value_ci = NA_real_))
  }
  d_pt <- est / sd_y
  rr_pt <- exp(0.91 * d_pt)
  ev_pt <- calc_ev(rr_pt)

  if (ci_low <= 0 && ci_high >= 0) {
    ev_ci <- 1
  } else {
    d_bound <- if (ci_low > 0) ci_low / sd_y else ci_high / sd_y
    ev_ci <- calc_ev(exp(0.91 * d_bound))
  }
  list(e_value = ev_pt, e_value_ci = ev_ci)
}

#' Add Simultaneous Uniform Confidence Bands via Fast EIF Multiplier Bootstrap
#' @noRd
add_simultaneous_bands <- function(estimates_tbl, eif_list, cluster_id, conf.level, n_boot = 2500L) {
  k <- nrow(estimates_tbl)
  n <- length(eif_list[[1L]])
  E <- do.call(cbind, eif_list)

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
  set.seed(12345L)

  cl_fac <- if (!is.null(cluster_id)) droplevels(factor(cluster_id)) else NULL
  if (!is.null(cl_fac) && nlevels(cl_fac) >= 2L) {
    # Aggregate EIFs by active clusters for cluster-valid multiplier bootstrap
    E_cl <- apply(E, 2L, function(col) as.numeric(tapply(col, cl_fac, sum)))
    J <- nrow(E_cl)
    E_centered <- sweep(E_cl, 2L, colMeans(E_cl), "-")
    ses <- sqrt(pmax((J / (J - 1L)) * colSums(E_centered^2) / (n^2), 1e-16))
    # Wild Rademacher multipliers (+1 / -1)
    xi <- matrix(sample(c(-1, 1), size = J * n_boot, replace = TRUE), nrow = J, ncol = n_boot)
    T_boot <- abs(crossprod(E_centered, xi) / n) / ses * sqrt(J / (J - 1L))
  } else {
    E_centered <- sweep(E, 2L, colMeans(E), "-")
    ses <- sqrt(pmax(apply(E_centered, 2L, stats::var) / n, 1e-16))
    xi <- matrix(sample(c(-1, 1), size = n * n_boot, replace = TRUE), nrow = n, ncol = n_boot)
    T_boot <- abs(crossprod(E_centered, xi) / n) / ses
  }

  max_t <- apply(T_boot, 2L, max)
  crit_sim <- as.numeric(stats::quantile(max_t, probs = conf.level, names = FALSE))

  for (i in seq_len(k)) {
    est_name <- estimates_tbl$estimand[i]
    se_i <- estimates_tbl$std.error[i]
    if (est_name %in% c("rr", "or")) {
      log_est <- log(estimates_tbl$estimate[i])
      estimates_tbl$conf.low[i] <- exp(log_est - crit_sim * se_i)
      estimates_tbl$conf.high[i] <- exp(log_est + crit_sim * se_i)
    } else {
      est_i <- estimates_tbl$estimate[i]
      estimates_tbl$conf.low[i] <- est_i - crit_sim * se_i
      estimates_tbl$conf.high[i] <- est_i + crit_sim * se_i
    }
  }
  estimates_tbl$crit.value <- crit_sim
  estimates_tbl
}
