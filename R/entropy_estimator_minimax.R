#' Bayesian entropy estimator with the minimax Dirichlet prior
#'
#' @description
#' Bayesian entropy estimator based on the minimax prior \eqn{\mathbf p \sim \operatorname{Dirichlet}\big((\sqrt{n}/K, \dots, \sqrt{n}/K)\big)}, where \eqn{K} is the alphabet size and \eqn{n} is the sample size
#'
#'
#' @param bin_counts Vector of observed bin counts
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
#' entropy_estimator_minimax(bin_counts)
#' entropy_estimator_minimax(bin_counts, unit = "ln")
#' entropy_estimator_minimax(bin_counts, unit = "normalize")
entropy_estimator_minimax <- function(bin_counts,
                                 unit = c("log2", "ln", "normalize")){

  return(

    entropy_estimator_Dirichlet(bin_counts = bin_counts,
                                Dir_alpha = sqrt(sum(bin_counts))/length(bin_counts),
                                unit = unit)

  )

}
