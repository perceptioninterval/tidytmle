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
