# If ... non-empty,
# then assign shared arguments from ... to the method_args of each element in multiple_methods,
# (unless already defined in method_args, in which case the method_args value takes precedence).

#' Update multiple methods with shared arguments
#'
#' @param multiple_methods List of lists, each containing method as a character and the method arguments as a list (although names can differ, e.g. `"pt_method"` and `"pt_method_args"`)
#' @param shared_args_list List of shared arguments to be assigned to each element in multiple_methods (if not locally defined)
#'
#' @return Updated multiple_methods list
#'
#' @examples
#' multiple_methods <- list(m1 = list(pt_method = "ML", pt_method_args = NULL),
#'                          m2 = list(pt_method = "Dirichlet", pt_method_args = NULL))
#' print(multiple_methods)
#' shared_args_list <- list(Dir_alpha = 10)
#' multiple_methods <- update_multiple_methods_with_shared_args(multiple_methods, shared_args_list)
#' print(multiple_methods)
#'
#' multiple_methods <- list(m1 = list(pt_method = "ML", pt_method_args = NULL),
#'                          m2 = list(pt_method = "Dirichlet", pt_method_args = list(Dir_alpha = 5)))
#' print(multiple_methods)
#' shared_args_list <- list(Dir_alpha = 10)
#' multiple_methods <- update_multiple_methods_with_shared_args(multiple_methods, shared_args_list)
#' print(multiple_methods)
update_multiple_methods_with_shared_args <- function(multiple_methods,
                                                     shared_args_list){

  # If shared_args_list non-empty,
  # then assign shared arguments from ... to the method_args of each element in multiple_methods,
  # (unless already defined in method_args, in which case the method_args value takes precedence).

  if(shared_args_list %>% length() > 0){

    multiple_methods <- purrr::map(

      multiple_methods,

      function(multiple_methods_elt){

        # Store original names of the list element (e.g. "pt_method", "pt_method_args", "ci_method", "ci_method_args")
        orig_names <- names(multiple_methods_elt)

        # Rename list element to "method" and "method_args"
        names(multiple_methods_elt) <- c("method", "method_args")

        # Extract names of arguments for the method
        method_arg_names <- get(paste0("entropy_estimator_", multiple_methods_elt$method)) %>%
          formals() %>% names()

        # Find overlapping names
        shared_arg_names <- generics::intersect(names(shared_args_list), method_arg_names)

        # If relevant shared argument not locally defined, add it to method_args
        for(shared_arg_name in shared_arg_names){

          if(!(shared_arg_name %in% names(multiple_methods_elt$method_args))){
            multiple_methods_elt$method_args[[shared_arg_name]] <- shared_args_list[[shared_arg_name]]
          }

        }

        # Restore original names
        names(multiple_methods_elt) <- orig_names

        return(multiple_methods_elt)
      }
    )
  }

  return(multiple_methods)

}
