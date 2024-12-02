#' Credible interval for Shannon entropy using the Dirichlet posterior
#'
#' @description
#' This function estimates a credible interval (the Bayesian counterpart of a confidence interval) for Shannon entropy based on the Dirichlet posterior distribution. The posterior is either derived from observed bin counts and a Dirichlet prior, or from a pre-calculated Dirichlet posterior distribution passed via `pt_est_outputs`. Monte Carlo simulations are used to estimate the posterior distribution for entropy.
#'
#' @param bin_counts A numeric vector of observed bin counts. Unnecessary if `Dir_alpha_posterior` is provided.
#' @param Dir_alpha_prior A numeric value or vector representing the parameters of the Dirichlet prior distribution. Defaults to `1`, representing a uniform prior. Unnecessary if `Dir_alpha_posterior` is provided.
#' @param Dir_alpha_posterior Optional: A numeric vector specifying the parameters of the Dirichlet posterior distribution. If not provided, it is computed from the `bin_counts` and `Dir_alpha_prior`.
#' @param N_MC A numeric value indicating the number of Monte Carlo simulations to estimate the posterior entropy distribution. Defaults to `5000`.
#' @param unit A character string specifying the unit for entropy: `"log2"` for bits, `"ln"` for nats, or `"normalize"` for normalized entropy/evenness. Defaults to `"log2"`.
#' @param pt_est_output Optional: A list output from an entropy point estimation function (e.g., `entropy_estimator_Dirichlet`). If provided and contains `Dir_alpha_posterior`, the argument `Dir_alpha_posterior` will be updated.
#'
#' @return A list with the following components:
#' \item{ci}{A numeric vector containing the lower and upper bounds of the 95% highest posterior density interval for Shannon entropy (i.e. the shortest interval containing 95% of the probability mass of the posterior entropy distribution).}
#' \item{credible_mean}{The mean of the sampled entropy posterior distribution.}
#' \item{credible_median}{The median of the sampled entropy posterior distribution.}
#' \item{credible_cdf}{The estimated posterior entropy distribution as an empirical cumulative distribution function (ECDF) object.}
#'
#' @export
#'
#' @examples
#' # Example usage with uniform prior
#' entropy_ci_Dirichlet(bin_counts = c(4, 1, 5))
#'
#' # Example with a non-uniform Dirichlet prior
#' entropy_ci_Dirichlet(bin_counts = c(4, 1, 5), Dir_alpha_prior = c(1, 2, 1))
#'
#' # Example with output from an entropy estimator
#' pt_est <- entropy_estimator_Dirichlet(c(4, 1, 5))
#' entropy_ci_Dirichlet(pt_est_output = pt_est)

entropy_ci_Dirichlet <- function(bin_counts,
                                 Dir_alpha_prior = 1,
                                 Dir_alpha_posterior = NULL,
                                 N_MC = 5*10^3,
                                 unit = c("log2", "ln", "normalize"),
                                 pt_est_output = NULL
                                 ){

  unit <- match.arg(unit)

  # Extract Dir_alpha_posterior from pt_est_output if possible
  update_vars_in_current_scope(potential_update_list = pt_est_output,
                               possible_to_update_vec = c("Dir_alpha_posterior"))

  # Calculate Dir_alpha_posterior not given in Dir_alpha_posterior or pt_est_output
  if(is.null(Dir_alpha_posterior)){

    Dir_alpha_posterior <- entropy_estimator_Dirichlet(bin_counts = bin_counts,
                                                       Dir_alpha = Dir_alpha_prior)$Dir_alpha_posterior

  }

  # Estimate distribution of H( Dirichlet(Dir_alpha_posterior) ) via Monte Carlo simulation

  p_samples <- gtools::rdirichlet(N_MC, Dir_alpha_posterior) %>% t() %>% as.data.frame()

  H_samples <- purrr::map_dbl(p_samples,
                              \(p){
                                -sum( p[p > 0] * log2(p[p > 0]) )
                              })

  # Convert units if necessary
  if(unit != "log2"){

    alphabet_size <- length(Dir_alpha_posterior)

    H_samples <- purrr::map_dbl(

      as.list(H_samples),

      \(x) convert_entropy_units(x,
                                 current_unit = "log2",
                                 new_unit = unit,
                                 alphabet_size = alphabet_size)

    )

  }

  # Calcuate the 95% HPD interval

  hpdi <- coda::HPDinterval(as.mcmc(H_samples), prob = 0.95)

  ci <- c(ci_lower = hpdi[1, "lower"],
          ci_upper = hpdi[1, "upper"])


  # Save credible mean, median, and CDF of entropy
  credible_mean <- mean(H_samples)
  credible_median <- median(H_samples)
  credible_cdf <- stats::ecdf(H_samples)

  # Return results
  results_list <- list(ci = ci,
                       credible_mean = credible_mean,
                       credible_median = credible_median,
                       credible_cdf = credible_cdf)

  return(results_list)

}
