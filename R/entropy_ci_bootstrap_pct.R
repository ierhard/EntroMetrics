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


