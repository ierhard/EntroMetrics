#' Shannon Entropy Point and Confidence Interval Estimation
#'
#' This function calculates the Shannon entropy of a set of binned data using
#' specified point estimation and confidence interval (CI) methods.
#'
#' @param bin_counts A numeric vector of non-negative values representing the
#'   counts or frequencies in each bin. Zero counts are allowed.
#' @param pt_method Character string specifying the point estimation method for
#'   entropy. Options are:
#'   - `"ML"`: Maximum Likelihood (plugin estimator)
#'   - `"MM"`: Miller-Madow
#'   - `"Dirichlet"`: Dirichlet prior
#'   - `"Jeffrey"`: Jeffrey prior
#'   - `"Laplace"`: Laplace prior
#'   - `"SG"`: Schürmann-Grassberger
#'   - `"minimax"`: Minimax prior
#'   - `"multiple"`: Combine multiple methods.
#' @param pt_method_args Optional list of arguments specific to the chosen
#'   `pt_method`. Default is `NULL`.
#' @param ci_method Character string specifying the confidence interval (CI)
#'   estimation method. Options are:
#'   - `"automatic"`: Automatically selects the best CI method.
#'   - `"Wald"`: Wald interval.
#'   - `"Wald_alt"`: Alternative method based on asymptotic normality Basharin's variance approximation.
#'   - `"bootstrap_pct"`: Percentile bootstrap.
#'   - `"bootstrap_t"`: Bootstrap-t interval.
#'   - `"bootstrap_bca"`: Bias-corrected and accelerated bootstrap.
#'   - `"multiple"`: Combine multiple CI methods.
#' @param ci_method_args Optional list of arguments specific to the chosen
#'   `ci_method`. Default is `NULL`.
#' @param multiple_methods A list of lists defining combinations of point and
#'   CI methods to calculate. Each element of the list should specify a
#'   `pt_method`, `pt_method_args`, `ci_method`, and `ci_method_args`. Default
#'   is `NULL`.
#' @param unit Character string specifying the unit of entropy. Options are:
#'   - `"log2"`: Bits (base 2 logarithm).
#'   - `"ln"`: Nats (natural logarithm).
#'   - `"normalize"`: Normalized entropy.
#' @param conf_level Numeric value between 0 and 1 specifying the confidence
#'   level for interval estimation. Default is `0.95`.
#' @param ... Shared arguments. These will be passed to all entropy point estimator
#' functions and CI method functions with matching argument names, but will be overriden by
#' method-specific arguments from `pt_method_args`, `ci_method_args`, and `multiple_methods`.
#'
#' @details This function computes Shannon entropy and its associated
#' confidence intervals. It supports various estimation methods, allowing
#' users to compute results for multiple methods simultaneously using the
#' `multiple_methods` argument. If no `multiple_methods` list is provided,
#' the function creates a single-method list based on the specified arguments.
#'
#' @return A list containing:
#' \item{`pt_est`}{point estimate of entropy}
#' \item{`ci`}{bounds of the confidence interval}
#' If multiple methods are used, the result is a list of lists, each matching the structure of the single-method output.
#'
#' @examples
#' # Example bin counts
#' counts <- c(10, 20, 30, 40)
#'
#' # Basic usage with default settings
#' shannon_entropy(counts)
#'
#' # Specify point and CI methods
#' shannon_entropy(counts, pt_method = "ML", ci_method = "bootstrap_bca")
#'
#' # Using multiple methods
#' methods <- list(
#'   list(pt_method = "ML", ci_method = "Wald"),
#'   list(pt_method = "MM", ci_method = "bootstrap_pct")
#' )
#' shannon_entropy(counts, pt_method = "multiple", ci_method = "multiple",
#'                 multiple_methods = methods)
#'
#' # Using multiple methods with the same point method
#' bin_counts <- 10*1:5
#' multiple_methods <- list(m1 = list(ci_method = "Wald",
#'                                    ci_method_args = list(z_score_dist = "normal")),
#'                          m2 = list(ci_method = "Wald_alt",
#'                                    ci_method_args = list(z_score_dist = "normal")))
#' shannon_entropy(bin_counts,
#'                 pt_method = "MM",
#'                 ci_method = "multiple",
#'                 multiple_methods = multiple_methods)
#'
#' # Using multiple methods with the same CI method
#' bin_counts <- 10*1:5
#' multiple_methods <- list(m1 = list(pt_method = "ML",
#'                                    pt_method_args = NULL),
#'                          m2 = list(pt_method = "MM",
#'                                    pt_method_args = NULL))
#' shannon_entropy(bin_counts,
#'                 pt_method = "multiple",
#'                 ci_method = "Wald_alt",
#'                 multiple_methods = multiple_methods)
#'
#' @seealso [entropy_pt_est()], [entropy_ci()]
#'
#' @export
shannon_entropy <- function(bin_counts,
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
                                          "Wald_alt",
                                          "bootstrap_pct",
                                          "bootstrap_t",
                                          "bootstrap_bca",
                                          "multiple"),
                            ci_method_args = NULL,
                            multiple_methods = NULL,
                            unit = c("log2", "ln", "normalize"),
                            conf_level = 0.95,
                            ...){ # Shared args (will try to match)

  # # FOR TESTING ----------------------------------------------------------------
  # # Single method
  # bin_counts <- c(10, 20, 30, 40)
  # pt_method <- "ML"
  # pt_method_args <- NULL
  # ci_method <- "Wald"
  # ci_method_args <- NULL
  # multiple_methods <- NULL
  # unit <- "ln"
  # conf_level <- 0.95

  # # Multiple Methods
  # bin_counts <- c(10, 20, 30, 40)
  # pt_method <- "multiple"
  # pt_method_args <- NULL
  # ci_method <- "multiple"
  # ci_method_args <- NULL
  # multiple_methods <- list(
  #   m1 = list(pt_method = "ML",
  #             pt_method_args = NULL,
  #             ci_method = "Wald",
  #             ci_method_args = NULL),
  #   m2 = list(pt_method = "MM",
  #             pt_method_args = NULL,
  #             ci_method = "Wald",
  #             ci_method_args = NULL)
  # )
  # unit <- "ln"
  # conf_level <- 0.95

  # # Multiple Methods with same point method and different CI methods
  # bin_counts <- 10*1:5
  # pt_method <- "MM"
  # pt_method_args <- NULL
  # ci_method <- "multiple"
  # ci_method_args <- NULL
  # unit <- "log2"
  # conf_level <- 0.95
  # multiple_methods <- list(m1 = list(ci_method = "Wald", ci_method_args = NULL),
  #                          m2 = list(ci_method = "Wald_alt", ci_method_args = NULL))


  # #-----------------------------------------------------------------------------


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
      multiple_methods <- purrr::map(
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

      # Append pt_method and pt_method_args to each element of multiple_methods
      multiple_methods <- purrr::map(
        multiple_methods,
        \(method_list_elt){

          # Add pt_method and pt_method_args to each method in the list
          method_list_elt$pt_method <- pt_method
          method_list_elt$pt_method_args <- pt_method_args

          method_list_elt <- method_list_elt[c("pt_method", "pt_method_args", "ci_method", "ci_method_args")]

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
                                   multiple_methods = multiple_methods,
                                   unit = "log2",
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
                       unit = "log2",
                       ...)$pt_est
      }
      return(fct)
    }

  )

  # Calculate confidence intervals

  # Extract ci_method and ci_method_args from each elt in multiple_methods list
  ci_multiple_methods <- purrr::map(
    multiple_methods,
    function(method_list_elt){
      method_list_elt <- method_list_elt[3:4] # Extract ci_method and ci_method_args
      names(method_list_elt) <- c("method", "method_args")
      return(method_list_elt)
    }
  )

  # Use entropy_ci() function to calculate confidence intervals
  ci_outputs <- entropy_ci(bin_counts,
                           method = "multiple",
                           multiple_methods = ci_multiple_methods,
                           conf_level = conf_level,
                           unit = unit,
                           pt_est_fcts = pt_est_fcts,
                           pt_est_outputs = pt_est_outputs,
                           ...)

  # Merge point estimates with confidence interval results --------------------

  if(length(multiple_methods) == 1){

    result_list <- merge_lists_with_names(ci_outputs,
                                          pt_est_outputs)

    result_list <- result_list[c("pt_est", "ci", setdiff(names(result_list), c("pt_est", "ci")))]

  } else {

    result_list <- purrr::map2(ci_outputs, pt_est_outputs,
                               merge_lists_with_names)

    result_list <- purrr::map(result_list,
                              \(x) x[c("pt_est", "ci", setdiff(names(x), c("pt_est", "ci")))]
                              )

  }

  # Return final result -------------------------------------------------------
  return(result_list)
}
