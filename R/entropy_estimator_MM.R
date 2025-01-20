#' Miller-Madow entropy estimator
#'
#' @param bin_counts Vector of observed bin counts
#' @param unit Desired unit (`"log2"` for bits, `"ln"` for nats, and `"normalize"` for normalized entropy/evenness)
#'
#' @return List containing one element: the entropy point estimate (`"pt_est"`) in the specified unit
#' @export
#'
#' @examples
#' set.seed(123)
#' p <- c(.7, .15, .05, .05)
#' bin_counts <- rmultinom(1, 15, p)[,1]
#' print(bin_counts)
#' entropy_estimator_MM(bin_counts)
#' entropy_estimator_MM(bin_counts, unit = "ln")
#' entropy_estimator_MM(bin_counts, unit = "normalize")
entropy_estimator_MM <- function(bin_counts,
                                 unit = c("log2", "ln", "normalize")){

  unit <- match.arg(unit)

  # Ensure that bin counts are valid (throws error if not)
  valid_bin_counts(bin_counts)

  # Calculate MM estimator in bits

  # Approximate bias up to O(1/N^2), i.e. Miller-Madow bias correction
  K <- length(bin_counts) # Alphabet size, i.e. number of bins
  N <- sum(bin_counts) # Sample size
  approx_bias <- - log2(exp(1)) * (K-1)/(2*N)

  # ML entropy estimate
  H_ML <- entropy_estimator_ML(bin_counts, unit = "log2")$pt_est

  # Miller-Madow entropy estimate
  H_hat <- H_ML - approx_bias

  # Convert to desired units
  if(unit != "log2"){
    H_hat <- convert_entropy_units(entropy = H_hat,
                                   current_unit = "log2",
                                   new_unit = unit,
                                   alphabet_size = length(bin_counts)) # Does nothing unless unit = "normalize"
  }

  # Results
  results_list <- list(pt_est = H_hat)

  return(results_list)
}
