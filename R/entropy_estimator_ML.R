#' Maximum likelihood entropy estimator (plugin estimator)
#'
#' @param bin_counts Vector of observed bin counts
#' @param unit Desired unit ("log2" for bits, "ln" for nats, and "normalize" for normalized entropy/evenness)
#'
#' @return List containing two elements: "pt_est" in the specified unit; "prob_est_ML" the ML point estimates of bin probabilities
#' @export
#'
#' @examples
#' set.seed(123)
#' p <- c(.7, .15, .05, .05)
#' bin_counts <- rmultinom(1, 15, p)[,1]
#' print(bin_counts)
#' entropy_estimator_ML(bin_counts)
#' entropy_estimator_ML(bin_counts, unit = "ln")
#' entropy_estimator_ML(bin_counts, unit = "normalize")
entropy_estimator_ML <- function(bin_counts,
                                 unit = c("log2", "ln", "normalize")){

  unit <- match.arg(unit)

  # Ensure that bin counts are valid (throws error if not)
  valid_bin_counts(bin_counts)

  # Calculate ML estimator in bits
  # ML estimate of probabilities
  p_hat <- bin_counts/sum(bin_counts)

  # ML estimate of entropy
  H_hat <- - sum( p_hat[p_hat > 0] * log2(p_hat[p_hat > 0]) )

  # Convert to desired units
  if(unit != "log2"){
    H_hat <- convert_entropy_units(entropy = H_hat,
                                   current_unit = "log2",
                                   new_unit = unit,
                                   alphabet_size = length(bin_counts)) # Does nothing unless unit = "normalize"
  }

  # Results
  results_list <- list(pt_est = H_hat,
                       prob_est = p_hat)

  return(results_list)
}
