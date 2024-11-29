#' Merge Two Named Lists Without Duplicates
#'
#' @description
#' This function merges two named lists while ensuring that elements with the
#' same name appear only once in the resulting list. Only the first occurrence
#' of each name is kept.
#'
#' @param l1 A named list. The first input list to be merged.
#' @param l2 A named list. The second input list to be merged.
#'
#' @return A named list containing all unique elements from both input lists. If
#'   duplicate names exist, only the first occurrence is retained.
#' @examples
#' list1 <- list(a = list(1, 2, 3), b = list(4, 5))
#' list2 <- list(b = list(5, 6), c = list(7, 8))
#' merge_lists_with_names(list1, list2)
#' # Returns:
#' # $a
#' # [[1]]
#' # [1] 1
#' #
#' # [[2]]
#' # [1] 2
#' #
#' # [[3]]
#' # [1] 3
#' #
#' # $b
#' # [[1]]
#' # [1] 4
#' #
#' # [[2]]
#' # [1] 5
#' #
#' # $c
#' # [[1]]
#' # [1] 7
#' #
#' # [[2]]
#' # [1] 8
#' @export
merge_lists_with_names <- function(l1, l2) {
  # Combine the two lists
  combined <- c(l1, l2)

  # Deduplicate by names
  merged <- combined[!duplicated(names(combined))]

  # Return the merged list
  return(merged)
}
