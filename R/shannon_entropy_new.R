shannon_entropy_new <- function(bin_counts,
                                pt_method = c("ML",
                                              "MM",
                                              "Dirichlet",
                                                "Jeffrey",
                                                "Laplace",
                                                "SG",
                                                "minimax",
                                              "multiple"),
                                pt_method_args = NULL,
                                ci_method = c("automatic",
                                              "Wald",
                                              "bootstrap_pct",
                                              "bootstrap_t",
                                              "bootstrap_bca",
                                              "multiple"),
                                ci_method_args = NULL,
                                multiple_methods = NULL,
                                unit = c("log2", "ln", "normalize"),
                                conf_level = 0.95,
                                ...){ # Shared args (will try to match)

  # Checks: Ensure desired format ----------------------------------------------

  # Ensure that bin counts are valid (throws error if not)
  valid_bin_counts(bin_counts)

  # Ensure pt_method and ci_method are valid by matching them with provided options
  pt_method <- match.arg(pt_method)
  ci_method <- match.arg(ci_method)

  # Indicator: Calculating multiple methods?
  multiple <- "multiple" %in% c(pt_method, ci_method)

  # Handling multiple methods -------------------------------------------------
  if(multiple == TRUE){

    # If only pt_method == "multiple", append ci_method and args to multiple_methods list
    if(ci_method != "multiple"){

      # # For testing - Example of how multiple_methods might be defined
      # multiple_methods <- list(m1 = list(pt_method = "ML", pt_method_args = NULL),
      #                          m2 = list(pt_method = "MM", pt_method_args = NULL))

      # Append ci_method and ci_method_args to each element of multiple_methods
      multiple_methods <- map(
        multiple_methods,
        \(method_list_elt){

          # Add ci_method and ci_method_args to each method in the list
          method_list_elt$ci_method <- ci_method
          method_list_elt$ci_method_args <- ci_method_args

          return(method_list_elt)
        }
      )
    }
    # If only ci_method == "multiple", append pt_method and args to multiple_methods list
    if((pt_method != "multiple") & (ci_method == "multiple")){

      # # For testing - Example of how multiple_methods might be defined
      # pt_method <- "ML"
      # pt_method_args <- list(NULL)
      # multiple_methods <- list(m1 = list(ci_method = "bootstrap_pct",
      #                                    ci_method_args = list(B = 10^3)),
      #                          m2 = list(ci_method = "bootstrap_bca",
      #                                    ci_method_args = list(B = 10^3)))

      # Append pt_method and pt_method_args to each element of multiple_methods
      multiple_methods <- map(
        multiple_methods,
        \(method_list_elt){

          # Add pt_method and pt_method_args to each method in the list
          method_list_elt$pt_method <- pt_method
          method_list_elt$pt_method_args <- pt_method_args

          return(method_list_elt)
        }
      )
    }
  }

  # If no multiple_methods list provided, create one with single method --------
  if(is.null(multiple_methods)){
    multiple_methods <- list(list(pt_method = pt_method,
                                  pt_method_args = pt_method_args,
                                  ci_method = ci_method,
                                  ci_method_args = ci_method_args))
  }

  # Point estimation -----------------------------------------------------------

  # Compute point estimates for all methods in multiple_methods list
  pt_est_outputs <- entropy_pt_est(bin_counts,
                                   method = "multiple",
                                   multiple_methods = multiple_methods[1:2],
                                   ...)

  # CI calculations ------------------------------------------------------------

  # Save list of point estimator functions (one for each element of multiple_methods)
  # with sole argument bin_counts (useful for resampling methods)
  pt_est_fcts <- purrr::map(

    multiple_methods,

    function(method_list_elt){
      fct <- function(bin_counts){
        # Return point estimate for given bin_counts using specified method
        entropy_pt_est(bin_counts,
                       method = method_list_elt$pt_method,
                       method_args = method_list_elt$pt_method_args,
                       ...)$pt_est
      }
      return(fct)
    }

  )

  # Calculate confidence intervals

  # Use entropy_ci() function to calculate confidence intervals
  ci_outputs <- entropy_ci(bin_counts,
                           method = "multiple",
                           multiple_methods = multiple_methods,
                           conf_level = conf_level,
                           unit = unit,
                           pt_est_fcts = pt_est_fcts,
                           pt_est_outputs = pt_est_outputs,
                           ...)

  # Merge point estimates with confidence interval results --------------------

  if(length(multiple_methods) == 1){

    result_list <- merge_lists_with_names(pt_est_outputs$pt_est,
                                          ci_outputs$ci)

  } else {

    result_list <- purrr::map2(pt_est_outputs, ci_outputs,
                               merge_lists_with_names)

  }

  # Return final result -------------------------------------------------------
  return(result_list)
}
