#' @name entropy_ci_bootstrap_bca
#' @title Bootstrap bias-corrected and accelerated (BCa) confidence interval for entropy estimators
#'
#' @description
#' Computes the bias-corrected and accelerated (BCa) bootstrap confidence interval for an entropy estimator based on bin counts.
#'
#' @param bin_counts A numeric vector of bin counts representing the observed frequencies of discrete outcomes.
#' @param pt_est_fct A function to compute the point estimate of entropy. This function must take `bin_counts` as its sole argument and return the entropy estimate in bits as a scalar.
#' @param B An integer specifying the number of bootstrap samples to generate (default is \eqn{1000}).
#' @param conf_level A numeric value specifying the confidence level for the interval (default is \eqn{0.95}).
#' @param unit A character string specifying the unit of entropy. Options are `"log2"` (default), `"ln"`, or `"normalize"`. The result will be converted to the specified unit if needed.
#'
#' @return A list with the following components <- <- <- <- <- <- <- <- <- <- <- :
#' \item{ci}{A numeric vector of length 2 containing the lower and upper bounds of the confidence interval, expressed in the specified entropy unit.}
#' \item{pt_est}{The point estimate of entropy computed from the input `bin_counts`, expressed in the specified entropy unit.}
#'
#' @examples
#' # Define a simple entropy estimation function (e.g., Shannon entropy in bits)
#' shannon_entropy_fct <- function(bin_counts) {
#'   proportions <- bin_counts / sum(bin_counts)
#'   -sum(proportions[proportions > 0] * log2(proportions[proportions > 0]))
#' }
#'
#' # Example bin counts representing observed frequencies of discrete outcomes
#' bin_counts <- c(10, 20, 30, 40)
#'
#' # Compute the BCa bootstrap confidence interval for Shannon entropy
#' result <- entropy_ci_bootstrap_bca(
#'   bin_counts = bin_counts,
#'   pt_est_fct = shannon_entropy_fct,
#'   B = 1000,
#'   conf_level = 0.95,
#'   unit = "log2"
#' )
#'
#' # View the results
#' print(result)
#'
#' @export
#'
entropy_ci_bootstrap_bca <- function(bin_counts,
                                     pt_est_fct, # Point estimation function (must take bin_counts as sole argument and return an entropy estimate in bits as a scalar)
                                     B = 10^3,
                                     conf_level = 0.95,
                                     unit = c("log2", "ln", "normalize")
){

  # Match argument values
  unit <- match.arg(unit)

  # Convert bin_counts to data vector x_vec for bootstrap_t function
  x_vec <- rep(1:length(bin_counts), bin_counts)

  # Call general bootstrap BCa function
  results_list <- bootstrap_bca(
    x_vec = x_vec,
    pt_est_fct = \(x) pt_est_fct(as.numeric(table(x))),
    conf_level = 0.95,
    B = B
  )

  # Convert units if necessary
  if(unit != "log2"){

    alphabet_size <- length(bin_counts)

    results_list$ci <- purrr::map_dbl(

      results_list$ci,

      \(x) convert_entropy_units(x,
                                 current_unit = "log2",
                                 new_unit = unit,
                                 alphabet_size = alphabet_size)

    )

    results_list$pt_est <- convert_entropy_units(results_list$pt_est,
                                                 current_unit = "log2",
                                                 new_unit = unit,
                                                 alphabet_size = alphabet_size)
  }

  # Return results

  return(results_list)
}
