# Bootstrap-t method for confidence intervals of statistics

# cf. section 12.5 in *Introduction to the Bootstrap* by Efron & Tibshirani (1993)
bootstrap_t_basic <- function(x_vec,
                              pt_est_fct,
                              conf_level = 0.95,
                              B = 10^3,
                              B_sub = 10^2, # Number of nested bootstrap samples used for std. error estimates, if se_fct not provided
                              se_fct = NULL, # optional: analytic formula for the standard error of estimator
                              verbose = FALSE
                              ){

  # pt_est_fct <- function(x) pt_est_fct(x, ...)

  theta_hat <- pt_est_fct(x_vec)

  # Get bootstraps (each column of x_bootstraps_df is a bootstrap sample)
  x_bootstraps_df <- replicate(B,
                               sample(x_vec, replace = TRUE)) %>% as.data.frame()

  # Calculate bootstrap estimates
  theta_boot_hats <- purrr::map_dbl(x_bootstraps_df, pt_est_fct)

  # Calculate standard error estimate

  se_hat <- ifelse(is.null(se_fct),
                   sd(theta_boot_hats), # via bootstrap standard error, if no analytic formula provided
                   se_fct(x_vec)) # via analytic formula otherwise

  # Calculate standard errors estimates of the statistic for each bootstrap sample
  se_boots <- purrr::map_dbl(

    x_bootstraps_df,

    \(x_boot) {

      ifelse(is.null(se_fct),
             # If no analytic formula provided, use bootstrap standard error,
             # cf. Alorithm 6.1 in *Introduction to the Bootstrap*, Efron & Tibshirani (1993)
             sd(
               purrr::map_dbl(
                 replicate(B_sub, sample(x_boot, replace = TRUE)) %>% as.data.frame(),
                 pt_est_fct
               )
             ),
             # Otherwise, use the analytic formula
             se_fct(x_boot)
      )

    }

  )

  # Calculate bootstrap t-statistics
  t_stats <- (theta_boot_hats - theta_hat) / se_boots

  # Calculate confidence intervals
  alpha <- 1 - conf_level
  q_lower <- quantile(t_stats, alpha/2, type = 1) %>% as.numeric()
  q_upper <- quantile(t_stats, 1 - alpha/2, type = 1) %>% as.numeric()
  ci <- c(ci_lower = theta_hat - q_upper * se_hat,
          ci_upper = theta_hat - q_lower * se_hat)

  # Return results
  if (verbose == TRUE) {
    list(ci = ci,
         pt_est = theta_hat,
         se_est = se_hat,
         pt_boots = theta_boot_hats,
         se_boots = se_boots,
         t_stats = t_stats)
  } else {
    results_list <- list(ci = ci,
                         pt_est = theta_hat)
  }

  return(results_list)

}

# # Example
#
# set.seed(123)
# probs <- gtools::rdirichlet(1, alpha = rep(.5, 26))[1,]
# n_sample <- 50
#
# entropy_in_bits <- function(x){
#   bin_counts <- as.numeric(table(x))
#   p <- bin_counts/sum(bin_counts)
#   H <- -sum(p[p > 0] * log2(p[p > 0]))
#   return(H)
# }
#
# true_val <- -sum(probs[probs > 0] * log2(probs[probs > 0]))
#
# in_ci_vec <- logical(20)
#
# counter <- 0
#
# for(i in 1:length(in_ci_vec)){
#
#   x_vec <- sample(letters, n_sample, replace = TRUE, prob = probs)
#
#   ci <- bootstrap_t_basic(x_vec = x_vec,
#                           pt_est_fct = entropy_in_bits)$ci
#
#   in_ci_vec[i] <- ci[1] <= true_val & true_val <= ci[2]
#
#   counter <- counter + 1
#
#   print(
#     counter/length(in_ci_vec)
#   )
#
#   print(
#     mean(in_ci_vec[1:counter])
#   )
# }
#
# mean(in_ci_vec)

################################################################################
# Bootstrap-t method for confidence intervals of statistics with variance stabilization

bootstrap_t <- function(x_vec,
                        pt_est_fct,
                        conf_level = 0.95,
                        B = 10^3,
                        B_sub = 10^2, # Number of nested bootstrap samples used for std. error estimates, if se_fct not provided
                        se_fct = NULL, # optional: analytic formula for the standard error of estimator
                        variance_stabilization = TRUE # If TRUE, use the variance-stabilized bootstrap-t method
                        ){

  # # FOR TESTING-----------------------------------------------------------------
  # set.seed(123)
  # probs <- gtools::rdirichlet(1, alpha = rep(.5, 26))[1,]
  # n_sample <- 50
  # x_vec <- sample(letters, n_sample, replace = TRUE, prob = probs)
  #
  # entropy_in_bits <- function(x){
  #   bin_counts <- as.numeric(table(x))
  #   p <- bin_counts/sum(bin_counts)
  #   H <- -sum(p[p > 0] * log2(p[p > 0]))
  #   return(H)
  # }
  #
  # pt_est_fct <- entropy_in_bits
  #
  # conf_level <- .95
  # B <- 10^3
  # B_sub <- 25
  # se_fct <- NULL
  # #-----------------------------------------------------------------------------

  if (variance_stabilization == FALSE) {
    # Save results
    results_list <- bootstrap_t_basic(x_vec = x_vec,
                                      pt_est_fct = pt_est_fct,
                                      conf_level = conf_level,
                                      B = B,
                                      B_sub = B_sub,
                                      se_fct = se_fct,
                                      verbose = FALSE)

  } else {

    boot_t_basic_output <- bootstrap_t_basic(x_vec = x_vec,
                                             pt_est_fct = pt_est_fct,
                                             conf_level = conf_level,
                                             B = B,
                                             B_sub = B_sub,
                                             se_fct = se_fct,
                                             verbose = TRUE)

    x <- boot_t_basic_output[["pt_boots"]]
    y <- boot_t_basic_output[["se_boots"]]

    # # Plot
    # plot(x, y)

    # Fit curve using linear regression on cubic spline basis functions
    # ns stands for natural splines, which enforces linearity at the boundaries
    spline_model <- lm(y ~ splines::ns(x, df = 10))

    # Create a function from the fitted model
    se_spline_fct <- function(new_x) {
      predict(spline_model, newdata = data.frame(x = new_x))
    }

    # # Plot the original data and the fitted curve
    # plot(x, y,
    #      main = "Cubic Spline Fit",
    #      col = "blue", pch = 16,
    #      xlim = c(0, max(x) + 1))
    # curve(se_spline_fct, add = TRUE, col = "red", lwd = 2)
    #
    # plot(x, y,
    #      main = "Cubic Spline Fit",
    #      col = "blue", pch = 16)
    # curve(se_spline_fct, add = TRUE, col = "red", lwd = 2)

    fct_to_integrate <- function(x) {
      1/se_spline_fct(x)
    }

    # # Plot g'
    # plot(x, 1/y,
    #      main = "Inverse of cubic spline fit of standard error of estimator",
    #      col = "blue", pch = 16,
    #      xlim = c(0, max(x) + 1))
    # curve(fct_to_integrate, add = TRUE, col = "red", lwd = 2)

    # Integrate the function from 0 to t
    g <- function(t){
      integrate(fct_to_integrate, lower = 0, upper = t)$value
    }

    # # Plot the curve for the function g
    # grid_pts <- seq(min(x), max(x), length.out = 100)
    #
    # plot(grid_pts,
    #      purrr::map_dbl(grid_pts, g),
    #      type = "l",
    #      main = "Transformation g",
    #      col = "blue", pch = 16)

    # Transformed statistic function g(theta_hat)
    g_pt_est_fct <- function(x) {
      g(pt_est_fct(x))
    }

    # Bootstrap-t output with transformed statistic
    output_transformed <- bootstrap_t_basic(x_vec = x_vec,
                                            pt_est_fct = g_pt_est_fct,
                                            conf_level = conf_level,
                                            B = B,
                                            B_sub = B_sub,
                                            se_fct = se_fct)

    # Define the inverse of g
    inverse_g <- function(y,
                          lower = 0,
                          upper = log2(30)) {
      # Wrapper for uniroot to find x such that g(x) = y
      root <- uniroot(function(x) g(x) - y, lower = lower, upper = upper)

      # Return the root (inverse of g at y)
      root$root
    }

    # Transform the confidence interval back to the original scale
    ci <- c(ci_lower = inverse_g(output_transformed$ci[1]),
            ci_upper = inverse_g(output_transformed$ci[2]))

    # Save results
    results_list <- list(ci = ci,
                         pt_est = boot_t_basic_output$pt_est,
                         ci_basic = boot_t_basic_output$ci)

  }
  # Return results
  return(results_list)
}

# # Example
#
# set.seed(123)
# probs <- gtools::rdirichlet(1, alpha = rep(.5, 26))[1,]
# n_sample <- 50
#
# entropy_in_bits <- function(x){
#   bin_counts <- as.numeric(table(x))
#   p <- bin_counts/sum(bin_counts)
#   H <- -sum(p[p > 0] * log2(p[p > 0]))
#   return(H)
# }
#
# true_val <- -sum(probs[probs > 0] * log2(probs[probs > 0]))
#
# N_test <- 10
#
# coverage_vec_basic <- logical(N_test)
# coverage_vec_stable <- logical(N_test)
#
# counter <- 0
#
# for(i in 1:length(in_ci_vec)){
#
#   x_vec <- sample(letters, n_sample, replace = TRUE, prob = probs)
#
#   bootstrap_t_output <- bootstrap_t(x_vec = x_vec,
#                                     pt_est_fct = entropy_in_bits)
#
#   ci_stable <- bootstrap_t_output$ci
#   ci_basic <- bootstrap_t_output$ci_basic
#
#   coverage_vec_basic[i] <- ci_basic[1] <= true_val & true_val <= ci_basic[2]
#   coverage_vec_stable[i] <- ci_stable[1] <= true_val & true_val <= ci_stable[2]
#
#   counter <- counter + 1
#
#   print(
#     counter/length(in_ci_vec)
#   )
#
#   print(
#     mean(in_ci_vec[1:counter])
#   )
# }
#
# mean(coverage_vec_basic)
# mean(coverage_vec_stable)

















