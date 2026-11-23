# tidytmle 0.1.0

* Initial release of **tidytmle**, a cross-fitted Targeted Maximum Likelihood
  Estimation (TMLE) orchestrator for the `tidymodels` ecosystem.
* Added progressive-disclosure specification API (`tmle_spec()`,
  `set_preprocessing()`, `set_outcome_learners()`, `set_treatment_learners()`,
  `set_censoring_learners()`, `set_resamples()`, `set_targeting()`) accepting
  `parsnip` specs, `workflows`, and `workflowsets::workflow_set` objects.
* Added role-aware `recipes` splitting (`split_causal_recipes()`), automatic
  informative-missingness modeling (`C(A, W)`), smart stratified/grouped
  `rsample` cross-fitting, and native `hardhat` case weights support.
* Added bounded canonical logistic fluctuation, adaptive positivity truncation,
  joint EIF inference across `"ate"`, `"rr"`, `"or"`, `"ey1"`, `"ey0"`, `"att"`,
  and `"atc"`, cluster-robust EIF SEs, simultaneous uniform confidence bands,
  and VanderWeele-Ding E-values.
* Added instantaneous `retarget()` without ML refitting, `shift()` /
  `scale_tx()` Modified Treatment Policies, DR-Learner CATE estimation
  (`fit_cate()`, `qini_auc()`), 5-pillar causal trust auditing (`check_tmle()`),
  copy-paste-ready paper/blog reporting (`report_tmle()`), plain-English
  `summary()` / `explain_tmle()`, and Tufte-styled `autoplot()` diagnostics
  (`"effect"`, `"overlap"`, `"targeting"`, `"balance"`, `"sensitivity"`,
  `"learners"`, `"eif"`, `"dashboard"`).
