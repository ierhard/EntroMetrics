#' Update variables in the current scope based on a provided list
#'
#' @description
#' This function updates variables in the current scope based on a provided list of potential updates. If the list contains values corresponding to the names in the vector of possible updates, those values are assigned to the respective variables in the current environment.
#'
#' @param potential_update_list A list containing potential updates for variables. Each element of the list should be named, corresponding to the variable intended to be updated.
#' @param possible_to_update_vec A character vector of variable names that are eligible for updating. This specifies which variables can be updated from `potential_update_list`.
#'
#' @return None. The function updates variables directly in the parent environment.
#' @export
#'
#' @examples
#' # Example usage:
#' pt_est <- 0.1
#' prob_est <- c(0.2, 0.8)
#' var3 <- 0.5
#' potential_updates <- list(pt_est = 0.5, prob_est = c(0.1, 0.9), other_var = 0.3)
#' possible_updates <- c("pt_est", "prob_est")
#' update_vars_in_current_scope(potential_updates, possible_updates)
#' print(pt_est)
#' print(prob_est)
#' print(var3)
update_vars_in_current_scope <- function(potential_update_list, possible_to_update_vec){

  update_exists_names <- generics::intersect(names(potential_update_list), possible_to_update_vec)

  if(length(update_exists_names) > 0){

    for(name in update_exists_names){

      assign(name, potential_update_list[[name]], envir = parent.frame())

    }

  }

}
