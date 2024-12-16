#' Bayesian entropy estimator with Dirichlet prior
#'
#' @description
#' Bayesian entropy estimator based on the prior \eqn{\mathbf p \sim \operatorname{Dirichlet}(\boldsymbol \alpha)}
#'
#'
#' @param bin_counts Vector of observed bin counts
#' @param Dir_alpha Either (a) a numeric vector for the \eqn{\boldsymbol \alpha} parameter of the Dirichlet distriubiton or (b) a positive real \eqn{\alpha}, which is interpreted as the parameter \eqn{(\alpha , \dots, \alpha )}
#' @param est_type Type of entropy estimator to use: either `"mean_entropy_of_p"` for \eqn{\mathbb E[H(p)]} or `"entropy_of_mean_p"` for \eqn{H(\mathbb E[p])}
#' @param monte_carlo_samples Number of Monte Carlo samples to use for the `"mean_entropy_of_p"` estimator (default: \eqn{1000})
#' @param unit Desired unit (`"log2"` for bits, `"ln"` for nats, and `"normalize"` for normalized entropy/evenness)
#'
#' @return A list containing three elements: entropy point estimate in the specified unit, an estimate of the probabilities (mean of \eqn{\mathbf p \sim \operatorname{Dirichlet}(\boldsymbol \alpha_{\text{posterior}})}), and the parameter \eqn{\boldsymbol \alpha_{\text{posterior}}} from the posterior Dirichlet distribution
#' @export
#'
#' @examples
#' set.seed(123)
#' p <- c(.7, .15, .05, .05)
#' bin_counts <- rmultinom(1, 15, p)[,1]
#' print(bin_counts)
#' entropy_estimator_Dirichlet(bin_counts)
#' entropy_estimator_Dirichlet(bin_counts, Dir_alpha = c(4, 2, 1, 1), unit = "ln")
#' entropy_estimator_Dirichlet(bin_counts, unit = "normalize")
entropy_estimator_Dirichlet <- function(bin_counts,
                                        Dir_alpha = 1,
                                        est_type = c("mean_entropy_of_p",
                                                     "entropy_of_mean_p"),
                                        monte_carlo_samples = 10^3, # if estimator_type = "mean_entropy_of_p"
                                        unit = c("log2", "ln", "normalize")){

  unit <- match.arg(unit)

  est_type <- match.arg(est_type)

  # Validate Dir_alpha
  if (any(Dir_alpha <= 0)) {
    stop("Error: Dir_alpha (or all entries thereof) must be positive")
  }

  if (length(Dir_alpha) != 1 && length(Dir_alpha) != length(bin_counts)) {
    stop("Error: Dir_alpha must either be a scalar or a vector of the same length as bin_counts")
  }

  # If Dir_alpha scalar, expand to match bin_counts length
  if (length(Dir_alpha) == 1) {
    Dir_alpha <- rep(Dir_alpha, length(bin_counts))
  }

  # Ensure that bin counts are valid (throws error if not)
  valid_bin_counts(bin_counts)

  # Posterior Dir_alpha
  Dir_alpha_posterior <- Dir_alpha + bin_counts

  # Calculate the Bayesian point estimate of probabilities from posterior Dir(alpha)

  # Same as entropy::freqs.Dirichlet(y = bin_counts, a = Dir_alpha)
  p_hat <- Dir_alpha_posterior / sum(Dir_alpha_posterior)

  # Calculate entropy estimator in bits
  if (estimator_type == "mean_entropy_of_p"){
    MC_samples <- gtools::rdirichlet(n = monte_carlo_samples,
                                     alpha = Dir_alpha_posterior))
    MC_samples <- as.data.frame(t(MC_samples))

    MC_samples <- purrr::map_dbl(MC_samples,
                                 \(p) - sum(p[p > 0] * log2(p[p > 0])))

    H_hat <- mean(MC_samples)

    rm(MC_samples)

  } else if (estimator_type == "entropy_of_mean_p") {
    # Same as entropy::entropy.Dirichlet(y = bin_counts, a = Dir_alpha, unit = "log2")
    H_hat <- - sum(p_hat[p_hat > 0] * log2(p_hat[p_hat > 0]))
  }

  # Convert to desired units
  if(unit != "log2"){
    H_hat <- convert_entropy_units(entropy = H_hat,
                                   current_unit = "log2",
                                   new_unit = unit,
                                   alphabet_size = length(bin_counts)) # Does nothing unless unit = "normalize"
  }

  # Results
  results_list <- list(pt_est = H_hat,
                       prob_est = p_hat,
                       Dir_alpha_posterior = Dir_alpha_posterior)

  return(results_list)
}
