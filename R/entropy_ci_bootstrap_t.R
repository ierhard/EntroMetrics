entropy_ci_bootstrap_t <- function(bin_counts,
                                   pt_est_fct, # Point estimation function (must take bin_counts as sole argument and return an entropy estimate in bits as a scalar)
                                   B = 10^3,
                                   B_subboot = 10^3, # Number of bootstrap subsamples drawn from each bootstrap sample to estimate standard errors
                                   variance_stabilization = TRUE, # If TRUE, use the variance-stabilized bootstrap-t method
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

  # Calculate the original point estimate
  H_hat <- pt_est_fct(bin_counts)

  # Calculate bootstrap point estimates
  bootstrap_H_hats <- purrr::map_dbl(as.data.frame(bootstrap_bin_counts),
                                     pt_est_fct)

  # Calculate bootstrap standard errors
  bootstrap_ses <- purrr::map_dbl(as.data.frame(bootstrap_bin_counts),
                                  \(bin_counts_boot) {
                                    sqrt(var(
                                      purrr::map_dbl(
                                        as.data.frame(
                                          rmultinom(
                                            B_subboot,
                                            sum(bin_counts_boot),
                                            bin_counts_boot / sum(bin_counts_boot)
                                          )
                                        ),
                                        pt_est_fct
                                      )
                                    ))
                                  })

  # Compute the bootstrap-t statistics
  t_stats <- (bootstrap_H_hats - H_hat) / bootstrap_ses

  # Calculate confidence intervals using the bootstrap-t method
  alpha <- 1 - conf_level

  t_quantiles <- stats::quantile(t_stats, c(alpha/2, 1-alpha/2), type = 1) %>% as.numeric()

  se_H_hat <- sd(bootstrap_H_hats)

  ci <- c(ci_lower = H_hat - t_quantiles[2] * se_H_hat,
          ci_upper = H_hat - t_quantiles[1] * se_H_hat)

  # Convert units if necessary
  if(unit != "log2"){

    alphabet_size <- dim(bootstrap_bin_counts)[2]

    ci <- purrr::map_dbl(

      as.list(bootstrap_H_hats),

      \(x) convert_entropy_units(x,
                                 current_unit = "log2",
                                 new_unit = unit,
                                 alphabet_size = alphabet_size)

    )

    H_hat <- convert_entropy_units(H_hat,
                                   current_unit = "log2",
                                   new_unit = unit,
                                   alphabet_size = alphabet_size)
  }

  # Return results
  results_list <- list(ci = ci,
                       pt_est = H_hat)

  return(results_list)
}
