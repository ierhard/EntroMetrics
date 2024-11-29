#' Wald confidence interval (based on asymptotic normality)
#'
#' @description
#' Confidence interval based on the Wald method (normality assumption) for the entropy of a distribution and standard error estimates of the form \eqn{\hat S_n = \sqrt{\frac 1n \left[\big( \sum_{i=1}^s \hat p_i(\log_2 \hat p_i)^2\big) - \hat H^2 \right]}} due to Basharin's proof of asymptotic normality and \eqn{O(1/n^2)} approximation of the variance of the ML entropy estimator.
#'
#' @param bin_counts Vector of observed bin counts
#' @param conf_level Confidence level
#' @param z_score_dist Choose distribution from which to derive z-values: `"t"` (\eqn{t_{n-1}} distribution, where \eqn{n} is sample size; more conservative, especially recommended when \eqn{n \leq 30}) and `"normal"` (standard normal distribution)
#' @param pt_est Optional: Point estimate of entropy or a function (of bin counts) that returns a point estimate of entropy in bits. (For validity, the point estimator should share asymptotic variance with ML estimator for entropy, e.g. Millwer-Madow estimator.) By default ML estimator of entropy is used.
#' @param prob_est Optional: Vector of bin probabilites or a function (of bin counts) that returns a vector of bin probabilities. By default ML estimator of probabilities is used.
#'
#' @return
#' @export
#'
#' @examples
entropy_ci_Wald <- function(

  bin_counts,
  conf_level = 0.95,
  z_score_dist = c("t", "normal"),
  pt_est = NULL,
  prob_est = NULL,
  unit = c("log2", "ln", "normalize")

){

  z_score_dist <- match.arg(z_score_dist)
  unit <- match.arg(unit)

  # Sample size
  n <- sum(bin_counts)

  alpha <- 1 - conf_level

  # Wald CI (based on asymptotic normality)

  # Relevant z-score
  z <- ifelse(z_score_dist == "t", stats::qt(1-alpha/2, df = n-1), stats::qnorm(1-alpha/2))

  # Standard error estimate (based on Basharin [1959; doi: 10.1137/1104033])

  p_hat <- if (is.null(prob_est)) {
    entropy::freqs(y = bin_counts, method = "ML")
  } else if (is.function(prob_est)) {
    prob_est(bin_counts)
  } else if (!is.null(prob_est) & !is.function(prob_est)) {
    prob_est
  } else {
    stop("Invalid `prob_est` argument.")
  }

  nu_hat <- sum(p_hat[p_hat > 0] * (log2(p_hat[p_hat > 0]))^2)

  H_hat <- if (is.null(pt_est)) {
    entropy_estimator_ML(bin_counts)$pt_est
  } else if (is.function(pt_est)) {
    pt_est(bin_counts)
  } else if (!is.null(pt_est) & !is.function(pt_est)) {
    pt_est
  } else {
    stop("Invalid `pt_est` argument.")
  }

  std_error <- sqrt((nu_hat - trim(H_hat)^2) / n) %>% max(0)

  # confidence interval

  ci.lower <- H_hat - z * std_error

  ci.upper <- H_hat + z * std_error

  # Return results

  results_list <- list(ci = c(ci.lower = ci.lower, ci.upper = ci.upper),
                       pt_est = H_hat,
                       prob_est = p_hat,
                       std_error = std_error)

  if(unit != "log2"){

    results_list[setdiff(names(results_list), "prob_est")] <- purrr::map(

      results_list[setdiff(names(results_list), "prob_est")],

      \(x) convert_entropy_units(x,
                                 new_unit = unit,
                                 alphabet_size = length(bin_counts))

      )

  }

  return(results_list)

}
