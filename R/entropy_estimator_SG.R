#' Bayesian entropy estimator with the Schürmann-Grassberger Dirichlet prior
#'
#' @description
#' Bayesian entropy estimator based on the Schürmann-Grassberger prior \eqn{\mathbf p \sim \operatorname{Dirichlet}\big((1/K, \dots, 1/K)\big)}, where \eqn{k \hat= } alphabet size
#'
#'
#' @param bin_counts Vector of observed bin counts
#' @param est_type Type of entropy estimator to use: either `"mean_entropy_of_p"` for \eqn{\mathbb E[H(p)]} or `"entropy_of_mean_p"` for \eqn{H(\mathbb E[p])},
#' where the expectation is taken over the posterior Dirichlet distribution of \eqn{p}
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
#' entropy_estimator_SG(bin_counts)
#' entropy_estimator_SG(bin_counts, unit = "ln")
#' entropy_estimator_SG(bin_counts, unit = "normalize")
entropy_estimator_SG <- function(bin_counts,
                                 est_type = c("mean_entropy_of_p",
                                              "entropy_of_mean_p"),
                                 unit = c("log2", "ln", "normalize")){

  return(

    entropy_estimator_Dirichlet(bin_counts = bin_counts,
                                Dir_alpha = 1/length(bin_counts),
                                est_type = est_type,
                                unit = unit)

  )

}
