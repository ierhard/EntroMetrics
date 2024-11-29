entropy_ci <- function(bin_counts,
                       method = c("bootstrap_pct",
                                  "bootstrap_t",
                                  "bootstrap_bca",
                                  "Wald",
                                  "multiple"), # Available CI methods; "multiple" allows multiple methods.
                       method_args = NULL, # Arguments specific to the selected CI method.
                       multiple_methods = NULL, # Specify multiple methods and their arguments. Format: list(m1 = list(method = "bootstrap_pct", method_args = list(...)), m2 = ...)
                       pre_calc_bootstraps = NULL, # Optional: precomputed bootstrap samples to avoid redundant computations.
                       unit = c("log2", "ln", "normalize"),
                       conf_level = 0.95, # Confidence level for CI methods.
                       ... ){ # Shared args (will try to match), e.g. B = 10^3 or pre_calc_bootstraps for bootstrap methods or
                              # conf_level = 0.95 or pt_est_fct = entropy_estimator_ML

  # Ensure bin_counts is a vector of non-negative integers (e.g., counts in histogram bins).
  if(!all(bin_counts == floor(bin_counts) & bin_counts >= 0)){
    stop("Error: bin_counts must be a vector of non-negative whole numbers.")
  }

  # Match method argument to one of the allowed options; default is the first ("bootstrap_pct").
  method <- match.arg(method)

  # If multiple_methods is not provided, create a default list with one method.
  if(is.null(multiple_methods)){
    multiple_methods <- list(list(method = method, method_args = method_args))
  }

  # Helper function to calculate confidence intervals for a single method.
  # Assumes input of the form list(method = "bootstrap_pct", method_args = list(...)).
  calculate_fct <- function(multiple_methods_elt){

    # Extract method name and arguments from the list element.
    method <- multiple_methods_elt$method
    method_args <- multiple_methods_elt$method_args

    # Combine method_args with bin_counts for computation.
    all_args <- method_args
    all_args$bin_counts <- bin_counts

    # If using a bootstrap method and pre-calculated bootstraps are provided, include them in the arguments.
    if(stringr::str_starts(method, "bootstrap") & !is.na(pre_calc_bootstraps)){
      all_args$pre_calc_bootstraps <- pre_calc_bootstraps
    }

    # Use switch to call the appropriate function based on the method name.
    result <- switch(method,
                     "bootstrap_pct"  = do.call(entropy_ci_bootstrap_pct, all_args),
                     "bootstrap_t"    = do.call(entropy_ci_bootstrap_t,   all_args),
                     "bootstrap_bca"  = do.call(entropy_ci_bootstrap_bca, all_args),
                     "Wald"           = do.call(entropy_ci_Wald,          all_args)
    )

    return(result)

  }

  # Apply the calculate_fct function to each method in multiple_methods using purrr::map.
  result <- purrr::map(multiple_methods, calculate_fct)

  # If only one method is specified, simplify the output from a list of lists to a single list.
  if(length(result) == 1){
    result <- result[[1]]
  }

  return(result)

}
