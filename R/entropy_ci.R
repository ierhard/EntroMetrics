entropy_ci <- function(bin_counts,
                       method = c("Wald",
                                  "bootstrap_pct",
                                  "bootstrap_t",
                                  "bootstrap_bca",
                                  "multiple"), # Available CI methods; "multiple" allows multiple methods.
                       method_args = NULL, # Arguments specific to the selected CI method.
                       multiple_methods = NULL, # Specify multiple methods and their arguments. Format: list(m1 = list(method = "bootstrap_pct", method_args = list(...)), m2 = ...)
                       unit = c("log2", "ln", "normalize"),
                       conf_level = 0.95, # Confidence level for CI methods.
                       pt_est_fcts = NULL,
                       pt_est_outputs = NULL,
                       # pre_calc_bootstraps = NULL, # Optional: precomputed bootstrap samples to avoid redundant computations.
                       ...){ # Shared args (will try to match), e.g. B = 10^3 or pre_calc_bootstraps for bootstrap methods or
                             # conf_level = 0.95 or pt_est_fct = entropy_estimator_ML

  # Since entropy_pt_est() delists output when only one method is specified,
  # we "relist" this here to avoid errors.
  if("pt_est" %in% names(pt_est_outputs)){
    pt_est_outputs <- list(pt_est_outputs)
  }

  # Match arguments
  method <- match.arg(method)
  unit <- match.arg(unit)

  # Ensure that bin counts are valid (throws error if not)
  valid_bin_counts(bin_counts)

  # If no multiple_methods list provided, create a default list with the selected method
  if(is.null(multiple_methods)){
    multiple_methods <- list(list(method = method, method_args = method_args))
  }

  # If share arguments ... non-empty,
  # then assign these to the method_args of each element in multiple_methods,
  # (unless already defined in method_args, in which case the method_args value takes precedence).

  multiple_methods <- update_multiple_methods_with_shared_args(
    multiple_methods = multiple_methods,
    shared_args_list = list(...)
  )

  # If bootstrap methods are present, precalculate bootstrap subsamples for shared use
  # across multiple methods.

  method_names <- purrr::map_chr(multiple_methods,
                                 function(multiple_methods_elt){
                                   return(multiple_methods_elt[[1]]) # Extract method name
                                 }
                                 )

  if(any(stringr::str_starts(method_names, "bootstrap"))){

    # Extract B from method_args if present
    B_max <- purrr::map_dbl(multiple_methods,
                            \(multiple_methods_elt){
                              return(multiple_methods_elt$method_args$B)
                              }
                            ) %>% max()

    bootstrap_bin_counts <- rmultinom(B_max,
                                      sum(bin_counts),
                                      bin_counts/sum(bin_counts))

  }

  # Update multiple_methods with point estimate functions and bootstrap subsamples
  # as necessary.

  multiple_methods <- purrr::pmap(

    list(multiple_methods, pt_est_fcts, pt_est_outputs),

    function(multiple_methods_elt, pt_est_fct_elt, pt_est_output_elt){

      # Add point estimate output to each method_args list
      multiple_methods_elt$method_args$pt_est_output <- pt_est_output_elt

      # Add point estimate function where applicable
      elt_arg_names <- names(formals(get(

        paste0("entropy_ci_", multiple_methods_elt$method)

        )))

      if("pt_est_fct" %in% elt_arg_names){
        multiple_methods_elt$method_args$pt_est_fct <- pt_est_fct_elt
      }

      # Add precalculated bootstrap subsamples (for bin_counts) where applicable
      if(stringr::str_starts(multiple_methods_elt$method, "bootstrap")){
        multiple_methods_elt$method_args$pre_calc_bootstraps <- bootstrap_bin_counts
      }

      return(multiple_methods_elt)

    }

  )

  # Helper function to calculate confidence intervals for a single method.
  # Assumes input of the form list(method = "bootstrap_pct", method_args = list(...)).
  calculate_fct <- function(multiple_methods_elt){

    # Extract method name and arguments from the list element.
    method <- multiple_methods_elt$method
    method_args <- multiple_methods_elt$method_args

    # Combine method_args with bin_counts for computation.
    all_args <- method_args
    all_args$bin_counts <- bin_counts
    all_args$unit <- unit

    # Call the appropriate function based on the method name
    result <- do.call(

      get(paste0("entropy_ci_", method)),

      all_args

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
