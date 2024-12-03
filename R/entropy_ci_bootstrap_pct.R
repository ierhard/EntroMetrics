#' Confidence interval for Shannon entropy using the percentage bootstrap method
#'
#' @description
#' This function estimates a confidence interval for Shannon entropy using the percentage bootstrap method. Bootstrap samples are generated based on the observed bin counts, and entropy is estimated for each bootstrap sample. The interval is determined using the empirical cumulative distribution function (ECDF) of the bootstrap entropy estimates.
#'
#' @param bin_counts A numeric vector of observed bin counts.
#' @param pt_est_fct A function that computes the point estimate of entropy. The function must take `bin_counts` as its sole argument and return a scalar entropy estimate (in bits by default).
#' @param B A numeric value specifying the number of bootstrap resamples. Defaults to `1000`.
#' @param conf_level A numeric value representing the desired confidence level for the interval. Defaults to `0.95` (95% confidence level).
#' @param unit A character string specifying the unit for entropy: `"log2"` for bits, `"ln"` for nats, or `"normalize"` for normalized entropy/evenness. Defaults to `"log2"`.
#' @param pre_calc_bootstraps Optional: A matrix of bootstrap bin counts, such as that produced by `rmultinom()`. Each column should represent the bin counts for one bootstrap sample. If provided, these pre-calculated bootstrap samples will be used instead of generating new ones.
#'
#' @return A list with the following components:
#' \item{ci}{A numeric vector containing the lower and upper bounds of the confidence interval for Shannon entropy.}
#' \item{bootstrap_cdf}{The ECDF object representing the bootstrap distribution of entropy estimates.}
#'
#' @export
#'
#' @examples
#' # Example usage with observed bin counts
#' entropy_ci_bootstrap_pct(bin_counts = c(4, 1, 5),
#'                          pt_est_fct = \(bin_counts) entropy_estimator_ML(bin_counts)$pt_est)
#'
#' # Example with custom confidence level and number of bootstrap samples
#' entropy_ci_bootstrap_pct(bin_counts = c(4, 1, 5),
#'                          pt_est_fct = \(bin_counts) entropy_estimator_ML(bin_counts)$pt_est,
#'                          B = 2000,
#'                          conf_level = 0.90)
#'
#' # Example with pre-calculated bootstrap samples
#' pre_boot <- rmultinom(1000, sum(c(4, 1, 5)), c(4, 1, 5)/sum(c(4, 1, 5)))
#' entropy_ci_bootstrap_pct(pt_est_fct = \(bin_counts) entropy_estimator_ML(bin_counts)$pt_est,
#'                          pre_calc_bootstraps = pre_boot)

entropy_ci_bootstrap_pct <- function(bin_counts,
                                     pt_est_fct, # Point estimation function (must take bin_counts as sole argument and return an entropy estimate in bits as a scalar)
                                     B = 10^3,
                                     conf_level = 0.95,
                                     unit = c("log2", "ln", "normalize"),
                                     pre_calc_bootstraps = NULL # as output of rmultinom() function, i.e. matrix with each column representing the bin counts for a bootstrap subsample
                                     ){

  # Match argument values
  unit <- match.arg(unit)

  # Bootstrap bin counts, if not provided via pre_calc_bootstraps argument

  if(is.null(pre_calc_bootstraps)){
    bootstrap_bin_counts <- rmultinom(B,
                                      sum(bin_counts),
                                      bin_counts/sum(bin_counts))
  } else {
    bootstrap_bin_counts <- pre_calc_bootstraps
  }

  # Calculate bootstrap point estimates
  bootstrap_H_hats <- purrr::map_dbl(as.data.frame(bootstrap_bin_counts),
                                     pt_est_fct)

  # Convert to desired units

  # Convert units if necessary
  if(unit != "log2"){

    alphabet_size <- dim(bootstrap_bin_counts)[2]

    bootstrap_H_hats <- purrr::map_dbl(

      as.list(bootstrap_H_hats),

      \(x) convert_entropy_units(x,
                                 current_unit = "log2",
                                 new_unit = unit,
                                 alphabet_size = alphabet_size)

    )

  }

  # Calculate confidence intervals and bootstrap CDF
  alpha <- 1 - conf_level

  ci <- stats::quantile(bootstrap_H_hats,
                        c(alpha/2, 1-alpha/2),
                        type = 1) # type = 1 for pseudoinverse quantiles, i.e. standard quantile definition from probability theory

  names(ci) <- c("ci_lower", "ci_upper")

  bootstrap_cdf <- stats::ecdf(bootstrap_H_hats)

  # Return results
  results_list <- list(ci = ci,
                       bootstrap_cdf = bootstrap_cdf)

  return(results_list)
}


