test_that("Weighted and covariate logistic fluctuation submodels both solve EIF and agree on binary/continuous parity", {
  df_bin <- sim_causal_cohort(n = 350, outcome_type = "binary", seed = 5101)
  spec_bin <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, lab_1, lab_2, lab_3),
    estimand = c("ate", "rr", "or", "att", "atc")
  ) |>
    set_outcome_learners(parsnip::logistic_reg() |> parsnip::set_engine("glm")) |>
    set_treatment_learners(parsnip::logistic_reg() |> parsnip::set_engine("glm")) |>
    set_targeting(cv_tmle = FALSE, g_bounds = c(0.025, 0.975), fluctuation = "logistic") |>
    set_resamples(v_folds = 4L)

  fit_gwt <- fit(spec_bin, df_bin)
  fit_cov <- retarget(fit_gwt, fluctuation = "logistic_covariate")

  expect_lt(abs(glance(fit_gwt)$mean_eif), 1e-6)
  expect_lt(abs(glance(fit_cov)$mean_eif), 1e-6)
  expect_lt(abs(mean(fit_cov$targeting$eif_1)), 1e-6)
  expect_lt(abs(mean(fit_cov$targeting$eif_0)), 1e-6)
  # Both fluctuation submodels target the exact same semiparametric functional
  expect_equal(tidy(fit_gwt)$estimate, tidy(fit_cov)$estimate, tolerance = 5e-3)

  # Independent textbook base-R reference TMLE (mirrors tmle::tmle target.gwt = TRUE / FALSE)
  q_glm <- stats::glm(y ~ a + age + sex + lab_1 + lab_2 + lab_3, data = df_bin, family = stats::binomial())
  g_glm <- stats::glm(a ~ age + sex + lab_1 + lab_2 + lab_3, data = df_bin, family = stats::binomial())
  g_ref <- pmin(pmax(stats::predict(g_glm, type = "response"), 0.025), 0.975)
  d1 <- df_bin
  d1$a <- 1L
  d0 <- df_bin
  d0$a <- 0L
  q1_0 <- stats::predict(q_glm, newdata = d1, type = "response")
  q0_0 <- stats::predict(q_glm, newdata = d0, type = "response")
  qa_0 <- ifelse(df_bin$a == 1, q1_0, q0_0)

  # Weighted-intercept fluctuation (target.gwt = TRUE)
  e1 <- stats::coef(suppressWarnings(stats::glm(
    df_bin$y ~ 1 + offset(stats::qlogis(q1_0)),
    weights = as.numeric(df_bin$a == 1) / g_ref, family = stats::quasibinomial()
  )))[1L]
  e0 <- stats::coef(suppressWarnings(stats::glm(
    df_bin$y ~ 1 + offset(stats::qlogis(q0_0)),
    weights = as.numeric(df_bin$a == 0) / (1 - g_ref), family = stats::quasibinomial()
  )))[1L]
  ref_ate_gwt <- mean(stats::plogis(stats::qlogis(q1_0) + e1) - stats::plogis(stats::qlogis(q0_0) + e0))

  # Two-parameter clever-covariate fluctuation (target.gwt = FALSE)
  h1 <- as.numeric(df_bin$a == 1) / g_ref
  h0 <- as.numeric(df_bin$a == 0) / (1 - g_ref)
  ec <- stats::coef(suppressWarnings(stats::glm(
    df_bin$y ~ -1 + h1 + h0 + offset(stats::qlogis(qa_0)),
    family = stats::quasibinomial()
  )))
  ref_ate_cov <- mean(
    stats::plogis(stats::qlogis(q1_0) + ec[1L] / g_ref) -
      stats::plogis(stats::qlogis(q0_0) + ec[2L] / (1 - g_ref))
  )

  expect_equal(tidy(fit_gwt)$estimate[tidy(fit_gwt)$estimand == "ate"], unname(ref_ate_gwt), tolerance = 1e-6)
  expect_equal(tidy(fit_cov)$estimate[tidy(fit_cov)$estimand == "ate"], unname(ref_ate_cov), tolerance = 1e-6)

  # Continuous LMTP covariate fluctuation also solves the EIF score equation cleanly
  df_pol <- sim_causal_cohort(n = 250, outcome_type = "continuous", treatment_type = "continuous", seed = 5102)
  fit_pol_gwt <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = shift(1.0)) |>
    set_resamples(v_folds = 3L) |>
    fit(df_pol)
  fit_pol_cov <- retarget(fit_pol_gwt, fluctuation = "logistic_covariate")
  expect_lt(abs(glance(fit_pol_cov)$mean_eif), 1e-6)
  expect_equal(tidy(fit_pol_gwt)$estimate, tidy(fit_pol_cov)$estimate, tolerance = 0.05)
})

test_that("Sorted character covariates, fold-isolated rare factors, and .pred_res fallback work without warnings", {
  set.seed(5202)
  n <- 240L
  site <- rep(c("East", "North", "South", "West"), each = 60L)
  stage <- factor(c(rep("Rare", 3L), rep(c("I", "II", "III"), length.out = n - 3L)))
  x1 <- stats::rnorm(n)
  a <- stats::rbinom(n, 1L, stats::plogis(0.4 * x1 + ifelse(site == "West", 0.5, -0.3)))
  y <- 1.5 + 1.2 * a + 0.8 * x1 + ifelse(site == "North", 0.6, 0) + stats::rnorm(n)
  df_edge <- tibble::tibble(y = y, a = a, site = site, stage = stage, x1 = x1)

  expect_no_warning({
    fit_edge <- tmle_spec(y, a, c(site, stage, x1)) |>
      set_resamples(v_folds = 5L) |>
      fit(df_edge)
  })
  expect_s3_class(fit_edge, "tmle_fit")
  expect_false(anyNA(collect_metrics(fit_edge, summarize = FALSE)$.estimate))
  expect_lt(abs(glance(fit_edge)$mean_eif), 1e-6)
})

test_that("add_ranger_always_split preserves existing engine args and stratify = TRUE T-learner works cleanly", {
  ms <- parsnip::rand_forest(trees = 250L, min_n = 5L) |>
    parsnip::set_engine("ranger", num.threads = 1L, seed = 1L) |>
    parsnip::set_mode("regression")
  ms_upd <- tidytmle:::add_ranger_always_split(ms, "a")
  expect_true("always.split.variables" %in% names(ms_upd$eng_args))
  expect_true("num.threads" %in% names(ms_upd$eng_args))
  expect_true("seed" %in% names(ms_upd$eng_args))
  expect_equal(rlang::eval_tidy(ms_upd$eng_args$always.split.variables), "a")
  expect_equal(rlang::eval_tidy(ms_upd$eng_args$num.threads), 1L)

  # Arm-stratified Q fitting on continuous heterogeneous cohort & binary cohort
  df_het <- sim_causal_cohort(n = 300, outcome_type = "continuous", seed = 5301)
  truth_het <- attr(df_het, "truth")
  fit_strat <- tmle_spec(y, a, c(age, sex, lab_1, lab_2, lab_3), estimand = c("ate", "att")) |>
    set_outcome_learners("default", stratify = TRUE) |>
    set_resamples(v_folds = 4L) |>
    fit(df_het)
  expect_s3_class(fit_strat, "tmle_fit")
  expect_lt(abs(glance(fit_strat)$mean_eif), 1e-6)
  td_s <- tidy(fit_strat)
  expect_lt(abs(td_s$estimate[td_s$estimand == "ate"] - truth_het$true_ate), 3.5 * td_s$std.error[1L])
  expect_s3_class(check_tmle(fit_strat), "tmle_check")
  expect_s3_class(autoplot(fit_strat, type = "learners"), "ggplot")

  df_bin <- sim_causal_cohort(n = 280, outcome_type = "binary", seed = 5302)
  fit_strat_bin <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = c("ate", "rr")) |>
    set_outcome_learners("default", stratify = TRUE) |>
    set_resamples(v_folds = 3L) |>
    fit(df_bin)
  expect_lt(abs(glance(fit_strat_bin)$mean_eif), 1e-6)
})

test_that("Zero-mean scale_tx recovery and status-quo contrast EIF work cleanly", {
  set.seed(5401)
  n <- 400L
  w1 <- stats::rnorm(n)
  w2 <- stats::rnorm(n)
  a <- 0.8 * w1 - 0.8 * w2 + stats::rnorm(n)
  y_fn <- function(aa, x1, x2) 1.5 + 1.0 * aa + 0.6 * aa^2 + 0.7 * aa * x1 + x1 - 0.8 * x2
  y <- y_fn(a, w1, w2) + stats::rnorm(n)
  df_sc <- tibble::tibble(y = y, a = a, w1 = w1, w2 = w2)
  true_policy <- mean(y_fn(0.75 * a, w1, w2))
  true_contrast <- mean(y_fn(0.75 * a, w1, w2) - y_fn(a, w1, w2))

  fit_pol <- tmle_spec(y, a, c(w1, w2), estimand = scale_tx(0.75)) |>
    set_resamples(v_folds = 4L) |>
    fit(df_sc)
  td_pol <- tidy(fit_pol)
  expect_equal(td_pol$estimand, "scale_tx(0.75)")
  expect_lt(abs(td_pol$estimate - true_policy), 3.0 * td_pol$std.error)

  td_con <- tidy(fit_pol, contrast = "observed")
  expect_equal(td_con$estimand, "scale_tx(0.75) - observed")
  expect_lt(abs(td_con$estimate - true_contrast), 3.0 * td_con$std.error)
  expect_equal(td_con$estimate, td_pol$estimate - mean(y), tolerance = 1e-6)

  # Direct constructor shift(..., contrast = "observed") works end-to-end across S3 verbs
  fit_con_direct <- tmle_spec(y, a, c(w1, w2), estimand = shift(1.0, contrast = "observed")) |>
    set_resamples(v_folds = 3L) |>
    fit(df_sc)
  td_dir <- tidy(fit_con_direct)
  expect_equal(td_dir$estimand, "shift(1) - observed")
  expect_lt(abs(glance(fit_con_direct)$mean_eif), 1e-6)
  expect_s3_class(autoplot(fit_con_direct, type = "effect"), "ggplot")
  expect_output(summary(fit_con_direct), "relative to observed exposure")
  rep_con <- report_tmle(fit_con_direct, outcome_label = "score", treatment_label = "dose")
  expect_match(rep_con$results, "relative to observed dose")
})
