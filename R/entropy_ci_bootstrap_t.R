entropy_ci_bootstrap_t <- function(bin_counts,
                                   pt_est_fct, # Point estimation function (must take bin_counts as sole argument and return an entropy estimate in bits as a scalar)
                                   B = 10^3,
                                   B_sub = 10^2, # Number of nested bootstrap samples used for std. error estimates, if se_fct not provided
                                   se_fct = NULL, # Optional: analytic formula for the standard error of estimator
                                   variance_stabilization = FALSE, # If TRUE, use the variance-stabilized bootstrap-t method
                                   conf_level = 0.95,
                                   unit = c("log2", "ln", "normalize")#,
                                   # pre_calc_bootstraps = NULL # as output of rmultinom() function, i.e. matrix with each column representing the bin counts for a bootstrap subsample
){

  # # FOR TESTING ----------------------------------------------------------------
  # bin_counts <- 5*1:5
  # pt_est_fct <- \(bin_counts) {
  #   entropy_estimator_MM(bin_counts)$pt_est
  #   }
  # B = 10^3
  # B_sub = 10^2 # Number of nested bootstrap samples used for std. error estimates, if se_fct not provided
  # se_fct = NULL # Optional: analytic formula for the standard error of estimator
  # variance_stabilization = FALSE # If TRUE, use the variance-stabilized bootstrap-t method
  # conf_level = 0.95
  # unit = "normalize"
  #
  # entropy_ci_bootstrap_t(bin_counts = 5*1:5,
  #                        pt_est_fct = \(bin_counts) {
  #                          entropy_estimator_MM(bin_counts)$pt_est}
  #                        )
  # #-----------------------------------------------------------------------------

  # Match argument values
  unit <- match.arg(unit)

  # Convert bin_counts to data vector x_vec for bootstrap_t function
  x_vec <- rep(1:length(bin_counts), bin_counts)

  # Call general bootstrap-t function
  results_list <- bootstrap_t(
    x_vec = x_vec,
    pt_est_fct = \(x) pt_est_fct(as.numeric(table(x))),
    conf_level = 0.95,
    B = B,
    B_sub = B_sub, # Number of nested bootstrap samples used for std. error estimates, if se_fct not provided
    se_fct = se_fct, # optional: analytic formula for the standard error of estimator
    variance_stabilization = variance_stabilization # If TRUE, use the variance-stabilized bootstrap-t method
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
