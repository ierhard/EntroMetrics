shannon_entropy_new <- function(bin_counts,
                                pt_method = c("ML",
                                              "MM",
                                              "SHR",
                                              "multiple"),
                                pt_method_args = list(NULL),
                                ci_method = c("bootstrap_pct",
                                              "bootstrap_t",
                                              "bootstrap_bca",
                                              "Wald",
                                              "multiple"),
                                ci_method_args = list(NULL),
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
  pt_estimates_list <- entropy_pt_est(bin_counts,
                                      method = "multiple",
                                      multiple_methods = multiple_methods[1:2])

  # CI calculations ------------------------------------------------------------

  # Save list of point estimator functions (one for each element of multiple_methods)
  # with sole argument bin_counts (useful for bootstrapping)
  pt_estimator_fcts <- purrr::map(

    multiple_methods,

    function(method_list_elt){
      fct <- function(bin_counts){
        # Return point estimate for given bin_counts using specified method
        entropy_pt_est(bin_counts,
                       method = method_list_elt$pt_method,
                       method_args = method_list_elt$pt_method_args)$entropy_pt_est
      }
      return(fct)
    }

  )

  # Generate bootstraps if applicable -----------------------------------------

  # Check if any of the methods require bootstrapping
  is_bootstrap <- purrr::map_chr(multiple_methods,
                                 \(method_list_elt){
                                   method_list_elt$ci_method %>%
                                     stringr::str_starts("bootstrap") %>%
                                     return()
                                 })

  # If any bootstrap methods are required, generate bootstrap samples
  if(any(is_bootstrap)){
    # Find maximum number of bootstraps needed among all methods
    B_max <- purrr::map_int(multiple_methods[is_bootstrap],
                            \(method_list_elt){
                              B_temp <- method_list_elt$ci_method_args$B
                              B_temp <- if_else(is.null(B_temp),
                                                10^3,
                                                B_temp)
                              return(B_temp)
                            }) %>% max()

    # Generate multinomial bootstrap samples
    bootstrap_subsamples <- rmultinom(n = B_max,
                                      size = sum(bin_counts),
                                      prob = bin_counts / sum(bin_counts)) %>%
      as.data.frame()

    # Name each bootstrap sample column
    names(bootstrap_subsamples) <- paste0("b", 1:B_max)
  }

  # Append additional arguments for CI calculation ----------------------------

  multiple_methods <- purrr::map(
    multiple_methods,
    function(multiple_methods_elt){

      ci_method <- multiple_methods_elt$ci_method
      ci_fct_name <- paste0("entropy_ci_", ci_method)
      ci_fct_args <- names(formals(get(ci_fct_name)))

      # Add pt_est_fct to ci_method_args if needed and not provided
      if(("pt_est_fct" %in% ci_fct_args) &
         is.null(multiple_methods_elt$ci_method_args$pt_est_fct)){

        multiple_methods_elt$ci_method_args$pt_est_fct <- function(bin_counts){
          entropy_pt_est(bin_counts,
                         method = multiple_methods_elt$pt_method,
                         method_args = multiple_methods_elt$pt_method_args)$entropy_pt_est
        }
      }
      # Add pt_est to ci_method_args if needed and not provided
      if(("pt_est" %in% ci_fct_args) &
         is.null(multiple_methods_elt$ci_method_args$pt_est)){

        multiple_methods_elt$ci_method_args$pt_est <- entropy_pt_est(
          bin_counts,
          method = multiple_methods_elt$pt_method,
          method_args = multiple_methods_elt$pt_method_args)$entropy_pt_est
      }

      return(multiple_methods_elt)
    }
  )

  # Calculate confidence intervals --------------------------------------------

  # Use entropy_ci() function to calculate confidence intervals
  ci_list <- entropy_ci(bin_counts,
                        method = "multiple",
                        multiple_methods = multiple_methods,
                        pre_calc_bootstraps = bootstrap_subsamples)

  # Merge point estimates with confidence interval results --------------------

  if(length(multiple_methods) == 1){
    result_list <- c(pt_estimates_list, ci_list)
  } else {
    result_list <- purrr::map2(pt_estimates_list, ci_list,
                               function(L1, L2){
                                 return(c(L1, L2))
                               })
  }

  # Return final result -------------------------------------------------------
  return(result_list)
}
