test_that("check_tmle() evaluates causal trust pillars and handles positivity stress & LMTP", {
  df <- sim_causal_cohort(n = 300, outcome_type = "binary", seed = 2026)
  fit <- tmle_spec(y, a, c(age, sex, starts_with("lab_")), estimand = c("ate", "rr")) |>
    set_resamples(v_folds = 3L) |>
    fit(df)

  chk <- check_tmle(fit)
  expect_s3_class(chk, "tmle_check")
  expect_equal(nrow(chk), 5L)
  expect_true(all(c("check", "status", "metric", "summary", "why_it_matters", "next_step") %in% names(chk)))
  expect_output(print(chk), "TMLE Causal Health Audit")

  # Positivity stress triggers warn or fail with actionable next_step
  df_stress <- sim_causal_cohort(n = 300, outcome_type = "binary", positivity_stress = TRUE, seed = 2027)
  fit_stress <- tmle_spec(y, a, c(age, sex, starts_with("lab_"))) |>
    set_resamples(v_folds = 3L) |>
    set_targeting(g_bounds = c(0.15, 0.85)) |>
    fit(df_stress)
  chk_stress <- check_tmle(fit_stress)
  expect_true(any(chk_stress$status %in% c("warn", "fail")))

  # Continuous LMTP returns 4 checks (omitting binary balance check)
  df_cont <- sim_causal_cohort(n = 200, outcome_type = "continuous", treatment_type = "continuous", seed = 2028)
  fit_lmtp <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = shift(1.0)) |>
    set_resamples(v_folds = 3L) |>
    fit(df_cont)
  chk_lmtp <- check_tmle(fit_lmtp)
  expect_equal(nrow(chk_lmtp), 4L)
})

test_that("report_tmle() generates publication-ready Methods, Results, and Assumptions prose", {
  df <- sim_causal_cohort(
    n = 250,
    outcome_type = "binary",
    censoring_rate = 0.12,
    n_clusters = 10L,
    seed = 2029
  )
  df_no_delta <- dplyr::select(df, -delta)

  suppressMessages({
    fit <- tmle_spec(
      outcome = y,
      treatment = a,
      covariates = c(age, sex, starts_with("lab_")),
      id = clinic_id,
      estimand = c("ate", "rr")
    ) |>
      set_resamples(v_folds = 3L) |>
      fit(df_no_delta)
  })

  rep_paper <- report_tmle(
    fit,
    outcome_label = "recovery",
    treatment_label = "therapy",
    units = "percentage points",
    style = "paper",
    simultaneous = TRUE,
    e_value = TRUE
  )
  expect_s3_class(rep_paper, "tmle_report")
  expect_match(rep_paper$methods, "10 clusters")
  expect_match(rep_paper$methods, "non-missingness indicator")
  expect_match(rep_paper$methods, "VanderWeele-Ding")
  expect_match(rep_paper$results, "therapy shifted recovery")
  expect_match(rep_paper$assumptions, "No unmeasured confounding")
  expect_output(print(rep_paper), "Publication-Ready TMLE Report")

  rep_plain <- report_tmle(fit, style = "plain")
  expect_match(rep_plain$methods, "predict, then target")

  rep_md <- report_tmle(fit, style = "markdown")
  expect_match(as.character(rep_md), "## Methods")
  expect_output(print(rep_md), "## Results")
})

test_that("Tufte plots (effect default, non-redundant balance, CATE/Qini formatting) work without defects", {
  df <- sim_causal_cohort(n = 300, outcome_type = "binary", seed = 2030)
  fit <- tmle_spec(y, a, c(age, sex, starts_with("lab_")), estimand = c("ate", "rr", "or")) |>
    set_resamples(v_folds = 3L) |>
    fit(df)

  # Default autoplot(fit) is the forest contrast plot
  p_eff <- autoplot(fit)
  expect_s3_class(p_eff, "ggplot")
  expect_equal(p_eff$labels$title, "Estimated Causal Effect & Uncertainty")

  # Binary factor sex should appear once (sexM), not duplicated as sexF + sexM
  bal_tbl <- tidytmle:::compute_covariate_balance(fit)
  expect_false("sexF" %in% bal_tbl$covariate)
  expect_true("sexM" %in% bal_tbl$covariate)
  p_bal <- autoplot(fit, type = "balance")
  expect_s3_class(p_bal, "ggplot")

  # p-value formatting must never emit double operator 'p = <'
  expect_identical(tidytmle:::format_p_human(1e-6), "p < 0.0001")
  expect_identical(tidytmle:::format_p_human(0.0434), "p = 0.0434")

  df_het <- sim_causal_cohort(n = 250, outcome_type = "continuous", seed = 2031)
  fit_c <- tmle_spec(y, a, c(age, sex, starts_with("lab_"))) |>
    set_resamples(v_folds = 3L) |>
    fit(df_het)
  cate <- fit_cate(fit_c, v_folds = 3L)
  p_cate <- autoplot(cate, type = "cate")
  p_qini <- autoplot(cate, type = "qini")
  expect_s3_class(p_cate, "ggplot")
  expect_s3_class(p_qini, "ggplot")
  expect_false(grepl("p = <", p_cate$labels$subtitle, fixed = TRUE))
})

test_that("Regression tests: NA balance alignment, 1-level subgroup factors, chained retarget, cluster subsets, fit_cate NAs, conf.level, and multi-policy guard", {
  # 1. Exact SMD preservation with missing covariates
  df_na <- sim_causal_cohort(n = 200, outcome_type = "binary", seed = 10)
  df_na$lab_1[1:20] <- NA_real_
  df_na$sex[21:35] <- NA
  fit_na <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
    set_resamples(v_folds = 3L) |>
    fit(df_na)

  bal_na <- tidytmle:::compute_covariate_balance(fit_na)
  i1 <- df_na$a == 1L
  i0 <- df_na$a == 0L
  true_age_smd <- abs(mean(df_na$age[i1]) - mean(df_na$age[i0])) /
    sqrt(0.5 * (stats::var(df_na$age[i1]) + stats::var(df_na$age[i0])))
  expect_equal(bal_na$smd_raw[bal_na$covariate == "age"], true_age_smd, tolerance = 1e-6)

  # 2. fit_cate handles missing covariates silently and yields zero NA predictions
  expect_no_warning({
    cate_na <- fit_cate(fit_na, v_folds = 3L)
    preds_na <- predict(cate_na, new_data = df_na[1:25, ])
    expect_equal(sum(is.na(preds_na$.pred)), 0L)
  })

  # 3. 1-level character covariate inside subgroup & chained retarget() subgroup preservation
  df_sub <- sim_causal_cohort(n = 240, outcome_type = "binary", seed = 77)
  df_sub$region <- as.character(df_sub$sex)
  df_sub$wt <- hardhat::importance_weights(runif(nrow(df_sub), 0.7, 1.4))
  fit_sub_base <- tmle_spec(y, a, c(age, region, lab_1, lab_2), weights = wt, estimand = c("ate", "rr")) |>
    set_resamples(v_folds = 3L) |>
    fit(df_sub)

  fit_f <- retarget(fit_sub_base, subgroup = region == "F")
  expect_s3_class(check_tmle(fit_f), "tmle_check")
  expect_s3_class(autoplot(fit_f, type = "balance"), "ggplot")

  fit_f_strict <- retarget(fit_f, g_bounds = c(0.05, 0.95), conf.level = 0.90)
  expect_equal(glance(fit_f_strict)$n, glance(fit_f)$n)
  expect_identical(fit_f_strict$subgroup_expr, fit_f$subgroup_expr)
  expect_output(print(fit_f_strict), "90% CI")
  expect_equal(tidy(fit_f_strict), tidy(fit_f_strict, conf.level = 0.90))

  rep_90 <- report_tmle(fit_f_strict)
  expect_match(rep_90$methods, "region == \"F\"", fixed = TRUE)
  expect_match(rep_90$methods, "importance case weights")
  expect_match(rep_90$markdown, "90% CI")
  expect_match(rep_90$results, "probability scale")

  fit_reset <- retarget(fit_f_strict, subgroup = TRUE)
  expect_equal(glance(fit_reset)$n, nrow(df_sub))
  expect_null(fit_reset$subgroup_expr)

  # 4. Clustered subgroup missing factor levels & 1-cluster fallback warning
  df_cl <- sim_causal_cohort(n = 240, outcome_type = "continuous", n_clusters = 12L, seed = 91)
  fit_cl <- tmle_spec(y, a, c(age, sex, lab_1), id = clinic_id) |>
    set_resamples(v_folds = 3L) |>
    fit(df_cl)

  keep_clinics <- levels(df_cl$clinic_id)[1:6]
  fit_cl_6 <- retarget(fit_cl, subgroup = clinic_id %in% keep_clinics)
  td_6 <- tidy(fit_cl_6)
  expect_true(is.finite(td_6$std.error) && is.finite(td_6$conf.low))

  expect_warning(
    fit_cl_1 <- retarget(fit_cl, subgroup = clinic_id == "clinic_01"),
    "only 1 cluster"
  )
  expect_true(is.finite(tidy(fit_cl_1)$std.error))

  # 5. Multi-policy list aborts clearly instead of silently dropping policies
  expect_error(
    tmle_spec(y, a, c(age, sex), estimand = list(shift(1.0), scale_tx(0.9))),
    "one continuous Modified Treatment Policy at a time"
  )
})

test_that("Scientific fidelity guards: RCT Pillar 5, ATT/ATC balance weights, common-OR/NA E-values, honest report prose, and ATE-aligned fit_cate HC1", {
  # 1. RCT-like unconfounded cohort passes Model Learning Quality when Q is strong and g AUC ~ 0.50
  set.seed(8101)
  n <- 300L
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  a <- stats::rbinom(n, 1L, 0.5)
  y <- 1.0 + 1.5 * x1 + 0.8 * x2 + 0.7 * a + stats::rnorm(n)
  df_rct <- tibble::tibble(age = 50 + 10 * x1, lab_1 = x2, a = a, y = y)

  fit_rct <- tmle_spec(y, a, c(age, lab_1), estimand = c("ate", "ey1", "ey0")) |>
    set_treatment_learners(parsnip::logistic_reg() |> parsnip::set_engine("glm")) |>
    set_resamples(v_folds = 3L) |>
    fit(df_rct)
  chk_rct <- check_tmle(fit_rct)
  expect_identical(chk_rct$status[chk_rct$check == "Model Learning Quality"], "pass")

  # Marginal means (ey1, ey0) have NA E-values while ATE contrast has finite E-value
  td_ev <- tidy(fit_rct, e_value = TRUE)
  expect_true(is.finite(td_ev$e_value[td_ev$estimand == "ate"]))
  expect_true(all(is.na(td_ev$e_value[td_ev$estimand %in% c("ey1", "ey0")])))

  # Uncontrasted shift() also yields NA from E-value helper
  ev_shift <- tidytmle:::compute_vanderweele_evalue("shift(1)", 4.2, 0.2, 3.8, 4.6, y)
  expect_true(is.na(ev_shift$e_value) && is.na(ev_shift$e_value_ci))

  # Common binary outcome (> 15% prevalence) uses sqrt(OR) E-value formula
  df_bin_common <- sim_causal_cohort(n = 300L, outcome_type = "binary", seed = 8102)
  fit_bc <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = c("ate", "or", "ey1")) |>
    set_resamples(v_folds = 3L) |>
    fit(df_bin_common)
  td_bc <- tidy(fit_bc, e_value = TRUE)
  or_est <- td_bc$estimate[td_bc$estimand == "or"]
  rr_eq <- sqrt(or_est)
  expect_equal(td_bc$e_value[td_bc$estimand == "or"], rr_eq + sqrt(rr_eq * (rr_eq - 1)), tolerance = 1e-6)
  expect_true(is.na(td_bc$e_value[td_bc$estimand == "ey1"]))
  expect_output(explain_tmle(fit_bc), "percentage points")

  # 2. Estimand-aware ATT and ATC balance weights match exact odds weighting
  fit_att <- retarget(fit_bc, estimand = "att")
  fit_atc <- retarget(fit_bc, estimand = "atc")
  bal_ate <- tidytmle:::compute_covariate_balance(fit_bc)
  bal_att <- tidytmle:::compute_covariate_balance(fit_att)
  bal_atc <- tidytmle:::compute_covariate_balance(fit_atc)
  expect_identical(attr(bal_ate, "weight_scheme"), "IPW")
  expect_identical(attr(bal_att, "weight_scheme"), "ATT odds")
  expect_identical(attr(bal_atc, "weight_scheme"), "ATC odds")
  expect_false(isTRUE(all.equal(bal_ate$smd_wt, bal_att$smd_wt)))
  expect_match(check_tmle(fit_att)$summary[3L], "ATT odds weighting")
  expect_match(autoplot(fit_att, type = "balance")$labels$title, "ATT Odds Weighting")

  # 3. report_tmle prose reflects actual fluctuation submodel, high clipping, SUTVA, and MAR
  fit_lin <- retarget(fit_bc, fluctuation = "linear")
  rep_lin <- report_tmle(fit_lin, style = "paper")
  expect_match(rep_lin$methods, "weighted linear fluctuation")
  expect_match(rep_lin$assumptions, "without unmodeled cross-unit spillover")

  fit_none <- retarget(fit_bc, fluctuation = "none")
  rep_none <- report_tmle(fit_none, style = "paper")
  expect_match(rep_none$methods, "Targeting fluctuation was disabled|fluctuation = 'none'")
  expect_identical(check_tmle(fit_none)$status[1L], "fail")

  fit_clip <- retarget(fit_bc, g_bounds = c(0.35, 0.65))
  rep_clip <- report_tmle(fit_clip)
  expect_false(grepl("Only [2-9][0-9]\\.[0-9]% of observations", rep_clip$assumptions))
  expect_match(rep_clip$assumptions, "Substantial positivity stress")

  df_cens <- sim_causal_cohort(n = 220L, outcome_type = "continuous", censoring_rate = 0.2, seed = 8103)
  fit_cens <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), censoring = delta) |>
    set_resamples(v_folds = 3L) |>
    fit(df_cens)
  rep_cens <- report_tmle(fit_cens)
  expect_match(rep_cens$assumptions, "Missing At Random \\(MAR\\)")

  # 4. fit_cate() on an ATT fit retargets internally to ATE DR pseudo-outcome & produces HC1 BLP SEs
  fit_hetero <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = "ate") |>
    set_resamples(v_folds = 3L) |>
    fit(sim_causal_cohort(n = 240L, outcome_type = "continuous", seed = 8104))
  fit_hetero_att <- retarget(fit_hetero, estimand = "att")
  set.seed(11); c1 <- fit_cate(fit_hetero, v_folds = 3L)
  set.seed(11); c2 <- fit_cate(fit_hetero_att, v_folds = 3L)
  expect_equal(c1$phi, c2$phi, tolerance = 1e-8)
  expect_true(all(is.finite(tidy(c2)$std.error)))
  expect_s3_class(c2, "tmle_cate")
})
