# tidytmle 0.1.0

* Initial release of **tidytmle**, a cross-fitted Targeted Maximum Likelihood
  Estimation (CV-TMLE) orchestrator for the `tidymodels` ecosystem.
* Added progressive-disclosure specification API (`tmle_spec()`,
  `set_preprocessing()`, `set_outcome_learners()`, `set_treatment_learners()`,
  `set_censoring_learners()`, `set_resamples()`, `set_targeting()`) accepting
  `parsnip` specs, `workflows`, and `workflowsets::workflow_set` objects,
  including optional arm-stratified outcome ensembles (`stratify = TRUE`) and
  default treatment-interaction / natural-spline / treatment-splitting workflows.
* Added role-aware `recipes` splitting (`split_causal_recipes()`) blocking
  treatment leakage into `g(W)`, automatic informative-missingness modeling
  (`C(A, W)`), smart stratified/grouped `rsample` cross-fitting, and native
  `hardhat` case weights support.
* Added bounded arm-specific weighted-intercept (`fluctuation = "logistic"`)
  and two-parameter clever-covariate (`fluctuation = "logistic_covariate"`)
  submodels, adaptive positivity truncation (`Petersen et al., 2012; Gruber
  et al., 2022`), joint EIF inference across `"ate"`, `"rr"`, `"or"`, `"ey1"`,
  `"ey0"`, `"att"`, and `"atc"`, cluster-robust EIF SEs, simultaneous uniform
  confidence bands, and VanderWeele-Ding E-values (including common-outcome
  square-root OR approximation).
* Added instantaneous `retarget()` without ML refitting, continuous Modified
  Treatment Policies (`shift()`, `scale_tx()`) with status-quo contrasts
  (`contrast = "observed"`), DR-Learner CATE estimation (`fit_cate()`,
  `qini_auc()`) with HC1-robust BLP heterogeneity inference, 5-pillar causal
  trust auditing (`check_tmle()`), copy-paste-ready paper/blog reporting
  (`report_tmle()`), plain-English `summary()` / `explain_tmle()`, and
  estimand-aware Tufte-styled `autoplot()` diagnostics (`"effect"`, `"overlap"`,
  `"targeting"`, `"balance"`, `"sensitivity"`, `"learners"`, `"eif"`,
  `"dashboard"`).
* Added cross-package validation & benchmarking vignette
  (`vignette("benchmarks-and-validation", package = "tidytmle")`) against
  `tmle`, `AIPW`, and `lmtp`.
