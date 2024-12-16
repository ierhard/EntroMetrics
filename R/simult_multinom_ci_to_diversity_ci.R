#===============================================================================
# Confidence intervals for diversity measures
#===============================================================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function definition
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Goal: Calculate a confidence interval for a diversity measure f()
#       of the proportion parameters of a multinomial distribution.
#       Namely, one of the following: entropy, efficiency, Gini-Simpson index,
#                                     or normalized Gini-Simpson index.
#
# Input: Simultaneous CI for proportions of multinomial distribution, say C,
#        at desired level (e.g. 95%).
#
# Output: Calculate ci_lower = min(f(C)) and ci_upper = max(f(C)).
#         This is a conservative CI of the same level as the simultaneous CI.

simult_multinom_ci_to_diversity_ci <- function(
    simultaneous_prob_ci_lower,
    simultaneous_prob_ci_upper,
    diversity_measure = c("entropy",
                          "efficiency", # i.e. normalized entropy
                          "Gini-Simpson",
                          "Gini-Simpson_normalized")
  ){

  # library(tidyverse)

  # Number of categories in multinomial distribution
  n_p <- length(simultaneous_prob_ci_lower)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Maximization problem

  # Note: This is concave maximization problem, which is equivalent to a convex
  #       minimization problem. Therefore, we can use disciplined convex programming
  #       via {CVXR}.

  # Load package for disciplined convex optimization in R
  # library(CVXR)

  # Define variable to maximize over
  p <- Variable(n_p)

  # Define constraints
  const <- list(sum(p) == 1,
                p >= lower_vec,
                p <= upper_vec)

  # Define objective
  if (diversity_measure == "entropy") {
    obj <- CVXR::Maximize(sum(entr(p)))
  } else if (diversity_measure == "efficiency") {
    obj <- CVXR::Maximize(sum(entr(p))/log(n_p))
  } else if (diversity_measure == "Gini-Simpson") {
    obj <- CVXR::Maximize(1 - sum_squares(p))
  } else if (diversity_measure == "Gini-Simpson_normalized") {
    obj <- CVXR::Maximize((1 - sum_squares(p)) / (1 - 1/n_p))
  }

  # Define problem
  prob_max <- CVXR::Problem(obj, const)


  sol_max <- CVXR::solve(prob_max)

  # Get maximal point
  # sol_max$getValue(p)

  # Store maximum as the upper bound of the confidence interval.
  ci_upper <- sol_max$value

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Minimization problem

  # Note: This is no longer a convex minimization (or concave maximization) problem,
  #       so we can't use disciplined convex programming via {CVXR}.
  #       However, since the feasible set is compact and convex, we can use
  #       Bauer's maximum principle to find the minimal point of the concave function:
  #       - Identify the extreme points of the set (finite!).
  #       - Find the minimum and minimizing point(s) among the the extreme points.

  # Load {rcdd} package. We use this to identify extreme points.
  # library(rcdd)

  # Define the matrix of coefficients and the right-hand side vector that define
  # the polyhedron that we are maximizing over via the equation:
  # Ax <= b.

  A <- rbind(

    matrix(rep(1, n_p), ncol = n_p),

    matrix(rep(-1, n_p), ncol = n_p),

    diag(n_p),

    -diag(n_p)
  )

  b <- rbind(
    matrix(1,1,1),
    matrix(-1,1,1),
    matrix(upper_vec, nrow = n_p),
    matrix(-lower_vec, nrow = n_p)
  )

  # H-representation
  # This defines a digital formalization of the set {x : Ax = b}
  hrep <- rcdd::makeH(A, b)

  # Compute the V-representation (extreme points and rays)
  # This converts the H-representation of the feasible set defined above into
  # a representation of this set defined via the extreme points of the set.
  vrep <- rcdd::scdd(hrep)

  # Extract the extreme points from the V-representation
  if (nrow(vrep$output) == 1) {
    extreme_points <- vrep$output[,-(1:2)] %>% data.frame()
  } else if (nrow(vrep$output) > 1) {
    extreme_points <- t(vrep$output[,-(1:2)]) %>% data.frame()
  }
  names(extreme_points) <- paste0("v",1:ncol(extreme_points))

  # Print the extreme points
  # print(extreme_points)

  # Find the minimal point(s)

  # Define objective function
  if (diversity_measure == "entropy") {
    obj_fct <- function(p) {
      p <- p[p > 0]
      return( -sum(p*log(p)) )
    }
  } else if (diversity_measure == "efficiency") {
    obj_fct <- function(p) {
      # n_p <- length(p) # Commented out, since already defined above
      p <- p[p > 0]
      return( -sum(p*log(p)) / log(n_p) )
    }
  } else if (diversity_measure == "Gini-Simpson") {
    obj_fct <- function(p) {
      return( (1 - sum(p^2)) )
    }
  } else if (diversity_measure == "Gini-Simpson_normalized") {
    obj_fct <- function(p) {
      return( (1 - sum(p^2)) / (1 - 1/length(p)) )
    }
  }

  # Calculate the function values of the extreme points
  extreme_point_values <- map_dbl(extreme_points, obj_fct)

  # Store minimum as the lower bound of the confidence interval.
  ci_lower <- min(extreme_point_values)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Output

  results_list <- list(ci = c(ci_lower = ci_lower, ci_upper = ci_upper),
                       sol_max = sol_max,
                       extreme_points = extreme_points,
                       extreme_point_values = extreme_point_values
  )

  return(results_list)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # Example: lower and upper bounds in a simultaneous confidence interval
# lower_vec <- c(.05, .2, .1, .45)
# upper_vec <- c(.10, .3, .3, .50)
#
# simult_multinom_ci_to_diversity_ci(lower_vec,
#                                    upper_vec,
#                                    "efficiency")$ci
#
#
# # No further restrictions beyond p[i] in [0,1] and sum(p) == 1
# lower_vec <- rep(0, n_p)
# upper_vec <- rep(1, n_p)
#
# simult_multinom_ci_to_diversity_ci(lower_vec,
#                                    upper_vec,
#                                    "efficiency")$ci
#
# # Another test (with values from actual survey)
#
# lower_vec <- c(0.61111111, 0, 0, 0, 0.03703704)
#
# upper_vec <- c(0.8317557, 0.1836075, 0.2021261, 0.1650890, 0.2576816)
#
# simult_multinom_ci_to_diversity_ci(lower_vec,
#                                    upper_vec,
#                                    "efficiency")$ci
#







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# End
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Code below can be ignored (used in developing the code above)

# # Load package for disciplined convex optimization in R
# library(CVXR)
#
# # Number of question variables
# n_p <- 4

# # Solve maximum entropy problem
# p <- Variable(n_p)
#
# obj <- Maximize(sum(entr(p))/log(n_p))
#
# const <- list(sum(p) == 1,
#               p >= lower_vec,
#               p <= upper_vec)
#
# prob_maxEnt <- Problem(obj, const)
#
# sol_maxEnt <- solve(prob_maxEnt)
#
# sol_maxEnt$getValue(p)
# sol_maxEnt$value
#
# # Solve maximum Gini-Simpson index problem
#
# prob_GS <- prob_maxEnt
#
# objective(prob_GS) <- Maximize((1 - sum_squares(p))/(1-1/n_p))
#
# sol_GS <- solve(prob_GS)
#
# sol_GS$getValue(p)
# sol_GS$value
#
# # Minimization problems
#
# # Note: These are no longer convex minimization / concave maximization problems,
# #       so we can't use {CVXR}.
# #       However, since the feasible set is compact and convex, we can use
# #       Bauer's maximum principle to find the optimal minimum point of the concave function:
# #       - Identify the extreme points of the set (finite!)
# #       - Find the minimal point(s) among the the extreme points
#
# # Load {rcdd} package. We use this to identify extreme points.
# library(rcdd)
#
# # Define the matrix of coefficients and the right-hand side vector
# # For a polyhedron Ax <= b, 'A' is the matrix of coefficients, and 'b' is the right-hand side
# # Include your constraints here, considering the format
#
# A <- rbind(
#
#   matrix(rep(1, n_p), ncol = n_p),
#
#   matrix(rep(-1, n_p), ncol = n_p),
#
#   diag(n_p),
#
#   -diag(n_p)
# )
#
# b <- rbind(
#   matrix(1,1,1),
#   matrix(-1,1,1),
#   matrix(upper_vec, nrow = n_p),
#   matrix(-lower_vec, nrow = n_p)
# )
#
# # H-representation
# hrep <- makeH(A, b)
#
# # Compute the V-representation (extreme points and rays)
# vrep <- scdd(hrep)
#
# # Extract the extreme points
# if (nrow(vrep$output) == 1) {
#   extreme_points <- vrep$output[,-(1:2)] %>% data.frame()
# } else if (nrow(vrep$output) > 1) {
#   extreme_points <- t(vrep$output[,-(1:2)]) %>% data.frame()
# }
# names(extreme_points) <- paste0("v",1:ncol(extreme_points))
#
# # Print the extreme points
# print(extreme_points)
#
# # Find the minimal point
#
# entropy.normed <- function(p) {
#   n_p <- length(p)
#   p <- p[p > 0]
#   -sum(p*log(p))/log(n_p)
# }
#
# GiniSimpson.normed <- function(p) {(1 - sum(p^2)) / (1 - 1/length(p))}
#
# my_fct <- entropy.normed
# # my_fct <- GiniSimpson.normed
#
# fct.val_extreme.pts <- map_dbl(extreme_points, my_fct)
# fct.val_extreme.pts
#
# min_value <- min(fct.val_extreme.pts)
# min_value








