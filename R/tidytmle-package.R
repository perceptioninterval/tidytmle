#' @keywords internal
"_PACKAGE"

#' @importFrom rlang enquo quo_is_null quo_is_missing quo_text eval_tidy as_label sym caller_env :=
#' @importFrom generics fit tidy glance augment
#' @importFrom ggplot2 autoplot
#' @importFrom tune collect_metrics
NULL

#' @export
generics::fit

#' @export
generics::tidy

#' @export
generics::glance

#' @export
generics::augment

#' @export
ggplot2::autoplot

#' @export
tune::collect_metrics

utils::globalVariables(c(
  ".estimate", ".estimator", ".metric", ".pred", "arm", "conf.high",
  "conf.low", "covariate", "density", "eif", "estimand", "estimate", "fold",
  "g", "group", "h_mag", "is_outlier", "learner", "leverage", "metric_type",
  "nuisance", "obs_id", "pct_truncated", "smd", "term", "threshold", "val",
  "value", "weight", "x_init", "y_star"
))
