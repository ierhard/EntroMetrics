#' Bayesian entropy estimator with Dirichlet prior
#'
#' @description
#' Bayesian entropy estimator based on the prior \eqn{\mathbf p \sim \operatorname{Dirichlet}(\boldsymbol \alpha)}
#'
#'
#' @param bin_counts Vector of observed bin counts
#' @param Dir_alpha Either (a) a numeric vector for the \eqn{\boldsymbol \alpha} parameter of the Dirichlet distriubiton or (b) a positive real \eqn{\alpha}, which is interpreted as the parameter \eqn{(\alpha , \dots, \alpha )}
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
                                        unit = c("log2", "ln", "normalize")){

  unit <- match.arg(unit)

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

  # Calculate entropy estimator in bits
  H_hat <- entropy::entropy.Dirichlet(y = bin_counts,
                                      a = Dir_alpha,
                                      unit = "log2")

  # Calculate the Bayesian point estimate of probabilities from posterior Dir(alpha)
  p_hat <- entropy::freqs.Dirichlet(y = bin_counts,
                                    a = Dir_alpha)

  # Convert to desired units
  if(unit != "log2"){
    H_hat <- convert_entropy_units(entropy = H_hat,
                                   current_unit = "log2",
                                   new_unit = unit,
                                   alphabet_size = length(bin_counts)) # Does nothing unless unit = "normalize"
  }

  # Posterior Dir_alpha

  Dir_alpha_posterior <- Dir_alpha + bin_counts

  # Results
  results_list <- list(pt_est = H_hat,
                       prob_est = p_hat,
                       Dir_alpha_posterior = Dir_alpha_posterior)

  return(results_list)
}
