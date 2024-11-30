#' Wald confidence interval for Shannon entropy
#'
#' @description
#' This function computes a confidence interval for Shannon entropy using the Wald method, which relies on the assumption of asymptotic normality. The standard error estimate is based on Basharin's asymptotic normality proof and an \eqn{O(1/n^2)} approximation of the variance of the maximum likelihood (ML) entropy estimator.
#'
#' @param bin_counts A numeric vector of observed bin counts.
#' @param conf_level A numeric value indicating the confidence level for the interval (default is 0.95).
#' @param z_score_dist A character string specifying the distribution to derive z-values from: either "t" for the \eqn{t_{n-1}} distribution (recommended for sample sizes \eqn{n \leq 30}) or "normal" for the standard normal distribution. Defaults to "t".
#' @param pt_est Optional: A numeric value for the point estimate of entropy or a function that computes the entropy (in bits) given bin counts. If not provided, the ML estimator is used by default. The point estimator should share the asymptotic variance with the ML estimator (e.g., Miller-Madow estimator).
#' @param prob_est Optional: A numeric vector of bin probabilities or a function that returns a vector of probabilities given the bin counts. If not provided, the ML estimator for probabilities is used by default.
#' @param unit A character string specifying the unit for entropy: "log2" for bits, "ln" for nats, or "normalize" for normalized entropy/evenness. Defaults to "log2".
#' @param pt_est_output Optional: A list output from an entropy point estimation function (e.g., `entropy_estimator_ML`). If provided, this list can update `pt_est` and `prob_est` if they are included.
#'
#' @return A list with the following components:
#' \item{ci}{A numeric vector with the lower and upper bounds of the confidence interval for Shannon entropy.}
#' \item{pt_est}{The point estimate of Shannon entropy.}
#' \item{prob_est}{The estimated bin probabilities.}
#' \item{std_error}{The standard error of the entropy estimate.}
#'
#' @export
#'
#' @examples
#' # Compare t vs. normal distribution for z-values
#' entropy_ci_Wald(c(4, 1, 5), z_score_dist = "t")$ci
#' entropy_ci_Wald(c(4, 1, 5), z_score_dist = "normal")$ci
#'
#' # Using output from entropy_estimator_MM()
#' entropy_ci_Wald(1:10)
#' entropy_ci_Wald(1:10, pt_est_output = entropy_estimator_MM(1:10))


entropy_ci_Wald <- function(

  bin_counts,
  conf_level = 0.95,
  z_score_dist = c("t", "normal"),
  pt_est = NULL,
  prob_est = NULL,
  unit = c("log2", "ln", "normalize"),
  pt_est_output = NULL

){
  # Match argument values
  z_score_dist <- match.arg(z_score_dist)
  unit <- match.arg(unit)

  # Update pt_est and prob_est if pt_est_output is provided

  update_vars_in_current_scope(potential_update_list = pt_est_output,
                               possible_to_update_vec = c("pt_est", "prob_est"))

  # Sample size
  n <- sum(bin_counts)

  # Alpha level for CI
  alpha <- 1 - conf_level

  # Calculate relevant z-score
  z <- ifelse(z_score_dist == "t",
              stats::qt(1-alpha/2, df = n-1),
              stats::qnorm(1-alpha/2))

  # Standard error estimate (based on Basharin [1959; doi: 10.1137/1104033])

  p_hat <- if (is.null(prob_est)) {
    bin_counts/sum(bin_counts)
  } else if (is.function(prob_est)) {
    prob_est(bin_counts)
  } else if (!is.null(prob_est) & !is.function(prob_est)) {
    prob_est
  } else {
    stop("Invalid `prob_est` argument.")
  }

  nu_hat <- sum(p_hat[p_hat > 0] * (log2(p_hat[p_hat > 0]))^2)

  # Point estimate
  H_hat <- if (is.null(pt_est)) {
    entropy_estimator_ML(bin_counts)$pt_est
  } else if (is.function(pt_est)) {
    pt_est(bin_counts)
  } else if (!is.null(pt_est) & !is.function(pt_est)) {
    pt_est
  } else {
    stop("Invalid `pt_est` argument.")
  }

  # Calculate standard error
  std_error <- sqrt((nu_hat - trim(H_hat)^2) / n) %>% max(0)

  # Compute confidence interval bounds
  ci_lower <- H_hat - z * std_error
  ci_upper <- H_hat + z * std_error

  # Prepare results list
  results_list <- list(ci = c(ci_lower = ci_lower, ci_upper = ci_upper),
                       pt_est = H_hat,
                       prob_est = p_hat,
                       std_error = std_error)

  # Convert units if necessary
  if(unit != "log2"){

    results_list[setdiff(names(results_list), "prob_est")] <- purrr::map(

      results_list[setdiff(names(results_list), "prob_est")],

      \(x) convert_entropy_units(x,
                                 current_unit = "log2",
                                 new_unit = unit,
                                 alphabet_size = length(bin_counts))

      )

  }

  # Return results
  return(results_list)

}
