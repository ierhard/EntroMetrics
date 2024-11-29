#' Shannon entropy point estimation
#'
#' @param bin_counts Vector of observed bin counts
#' @param method Entropy estimator
#' @param method_args List of arguments required for method
#' @param multiple_methods List with one element for each desired estimation. This element is itself a list of two elements: (i) the method (character) and (ii) a list of method arguments
#' @param unit Desired unit (`"log2"` for bits, `"ln"` for nats, and `"normalize"` for normalized entropy/evenness)
#' @param ... Shared arguments. These will be passed to all entropy estimator functions with matching argument names, but will be overriden by method-specific arguments
#'
#' @return List of results of the entropy estimate(s). If `method = "multiple"`, the function returns a nested list where each primary element is the output from a user-specified `method`/`method_args` pair.
#' @export
#'
#' @examples
#' entropy_pt_est(bin_counts = 1:4, method = "MM")
#'
#' mult_methods_list <- list(
#'   m1 = list(method = "ML",
#'             method_args = NULL),
#'   m2 = list(method = "MM",
#'             method_args = NULL),
#'   m3 = list(method = "Dirichlet",
#'             method_args = list(Dir_alpha = 10))
#' )
#' entropy_pt_est(bin_counts = 1:20, method = "multiple", multiple_methods = mult_methods_list)
#'
#' mult_methods_list <- list(
#'   m1 = list(method = "ML",
#'             method_args = NULL),
#'   m2 = list(method = "MM",
#'             method_args = NULL),
#'   m3 = list(method = "Dirichlet",
#'             method_args = NULL)
#' )
#' entropy_pt_est(bin_counts = 1:20,
#'                method = "multiple",
#'                multiple_methods = mult_methods_list)
#'
#' entropy_pt_est(bin_counts = 1:20,
#'                method = "multiple",
#'                multiple_methods = mult_methods_list,
#'                Dir_alpha = 10) # Should be passed to Dirichlet estimator

entropy_pt_est <- function(bin_counts,
                           method = c("ML",
                                      "MM",
                                      "Dirichlet",
                                        "Jeffrey",
                                        "Laplace",
                                        "SG",
                                        "minimax",
                                      "multiple"),
                           method_args = NULL,
                           multiple_methods = NULL,
                           unit = c("log2", "ln", "normalize"),
                           ...){

  method <- match.arg(method)
  unit <- match.arg(unit)

  # Ensure that bin counts are valid (throws error if not)
  valid_bin_counts(bin_counts)

  # If method = "own",

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

  # Calculate the point estimate for each element of multiple_methods --
  # i.e. for each method, method_args pair
  result <- purrr::map(

    multiple_methods,

    function(multiple_methods_elt){

      # Extract method and arguments
      method <- multiple_methods_elt[[1]] # method name
      method_args <- multiple_methods_elt[[2]] # list of arguments

      # Add bin_counts to method arguments
      all_args <- method_args
      all_args$bin_counts <- bin_counts
      all_args$unit <- unit

      # Calculate point estimate using the appropriate estimator function
      result <- switch(method,
                       "ML" = do.call(entropy_estimator_ML, all_args),  # Maximum Likelihood estimator
                       "MM" = do.call(entropy_estimator_MM, all_args),  # Method of Moments estimator
                       "Dirichlet" = do.call(entropy_estimator_Dirichlet, all_args),
                         "Jeffrey" = do.call(entropy_estimator_Jeffrey, all_args),
                         "Laplace" = do.call(entropy_estimator_Laplace, all_args),
                         "SG" = do.call(entropy_estimator_SG, all_args),
                         "minimax" = do.call(entropy_estimator_minimax, all_args)
      )

      return(result)
    }

    )

  # If only one method specified, return result directly instead of a list
  if(length(result) == 1){
    result <- result[[1]]
  }

  return(result)
}
