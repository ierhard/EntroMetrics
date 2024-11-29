#' Checks validity of `bin_counts` parameter
#'
#' @description
#' Checks if `bin_counts` parameter is a vector of non-negative integers, not all zero.
#'
#'
#' @param bin_counts Parameter bin_counts to check
#'
#' @return `NULL` if valid, stops function with an error message otherwise
#'
#' @examples
#' valid_bin_counts(-1:10)
#' valid_bin_counts(c(3, 3.1, 4))
#' valid_bin_counts(1:5)
valid_bin_counts <- function(bin_counts){

  if(!all(bin_counts == floor(bin_counts) & bin_counts >= 0)){
    stop("Error: bin_counts must be a vector of non-negative whole numbers.")
  }
  if(all(bin_counts == 0)){
    stop("Error: bin_counts must contain at least one positive entry.")
  }

}
