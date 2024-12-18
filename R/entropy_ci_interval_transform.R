entropy_ci_interval_transform <- function(bin_counts,
                                          conf_level = 0.95,
                                          multinom_prop_ci = c("Sison-Glaz", "Goodman"),
                                          unit = c("log2", "ln", "normalize")
){
#
#   # TESTING --------------------------------------------------------------------
#
#   bin_counts <- 5*1:5
#   conf_level = 0.95
#   multinom_prop_ci = "Sison-Glaz"
#   unit = "log2"
#
#   # ----------------------------------------------------------------------------

  # Match argument values
  multinom_prop_ci <- match.arg(multinom_prop_ci)
  unit <- match.arg(unit)

  # Ensure that bin_counts is provided
  if(is.null(bin_counts)){
    stop("Please provide the argument bin_counts.")
  }

  # Ensure that bin_counts is valid
  valid_bin_counts(bin_counts)

  # Calculate simultaneous CI for multinomial proportions (i.e. probability vector parameter)

  switch(

    multinom_prop_ci,

    "Sison-Glaz" = {

      ci_multinomial_props <- MultinomialCI::multinomialCI(x = bin_counts,
                                                            alpha = 1 - conf_level)

    },

    "Goodman" = {

      stop("Goodman method not yet implemented.")

    }

  )
  # Calculate CI for multinomial proportions

  ci <- simult_multinom_ci_to_diversity_ci(
    simultaneous_prob_ci_lower = ci_multinomial_props[,1],
    simultaneous_prob_ci_upper = ci_multinomial_props[,2],
    diversity_measure = c("entropy")
  )$ci

  # Convert units if necessary
  if(unit != "ln"){

    alphabet_size <- length(bin_counts)

    ci <- convert_entropy_units(ci,
                                current_unit = "ln",
                                new_unit = unit,
                                alphabet_size = alphabet_size)

  }


  # Return results
  results_list <- list(ci = ci)

  return(results_list)
}


