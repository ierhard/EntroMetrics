# install.packages("extraDistr", "randtoolbox", "stats")

library(randtoolbox)   # For QMC sampling
library(stats)         # For Kolmogorov-Smirnov test

# Shannon entropy functional
discrete_entropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
}

# Monte Carlo method
mc_entropy <- function(alpha, num_samples) {
  entropy_values <- numeric(num_samples)
  for (i in 1:num_samples) {
    p <- gtools::rdirichlet(1, alpha)[1,]
    entropy_values[i] <- discrete_entropy(p)
  }
  return(entropy_values)
}

# Quasi-Monte Carlo method
qmc_entropy <- function(alpha, num_samples) {
  K <- length(alpha)

  # Generate Sobol sequence
  qmc_points <- sobol(num_samples, dim = K)

  # Transform Sobol points to Dirichlet samples
  entropy_values <- numeric(num_samples)
  for (i in 1:num_samples) {
    u <- qmc_points[i, ]
    g <- qgamma(u, shape = alpha, rate = 1)
    p <- g / sum(g)
    entropy_values[i] <- discrete_entropy(p)
  }
  return(entropy_values)
}

library(coda)  # For HPDI computation

# Function to compute 95% HPDI
compute_hpd_interval <- function(samples, prob = 0.95) {
  hpd <- HPDinterval(as.mcmc(samples), prob = prob)
  return(c(hpd[1, "lower"], hpd[1, "upper"]))
}

# # Speed of convergence test with HPDI comparison
# speed_of_convergence_hpd <- function(alpha, num_samples, step = 100) {
#   # Monte Carlo
#   mc_entropies <- mc_entropy(alpha, num_samples)
#
#   # Quasi-Monte Carlo
#   qmc_entropies <- qmc_entropy(alpha, num_samples)
#
#   # Compute L2 norm between HPDIs at each step
#   l2_norms <- numeric(num_samples / step)
#   for (i in seq(step, num_samples, by = step)) {
#     mc_subset <- mc_entropies[1:i]
#     qmc_subset <- qmc_entropies[1:i]
#
#     # Compute HPDIs
#     mc_hpd <- compute_hpd_interval(mc_subset)
#     qmc_hpd <- compute_hpd_interval(qmc_subset)
#
#     # Compute L2 norm of HPDI difference
#     l2_norms[i / step] <- sqrt(sum((mc_hpd - qmc_hpd)^2))
#   }
#
#   return(l2_norms)
# }
#
# # Plot L2 norms
# alpha <- rep(1, 100)  # Symmetric Dirichlet
# num_samples <- 10000  # Total number of samples
# step <- 100          # Step size for HPDI calculation
#
# l2_norms <- speed_of_convergence_hpd(alpha, num_samples, step)
#
# # Plot results
# hpd_steps <- seq(step, num_samples, by = step)
# plot(x = hpd_steps, y = l2_norms, type = "b", col = "purple", pch = 16,
#      xlab = "Number of Samples", ylab = "L2 Norm of HPDI Difference",
#      main = "L2 Norm Between MC and QMC HPDIs")
# abline(h = 0, col = "gray", lty = 2)


plot_ci_bounds_with_timing <- function(alpha, num_samples, step = 100) {
  # Measure timing
  timing_results <- compare_computation_speed(alpha, num_samples)

  # Monte Carlo
  mc_entropies <- mc_entropy(alpha, num_samples)

  # Quasi-Monte Carlo
  qmc_entropies <- qmc_entropy(alpha, num_samples)

  # Store CI bounds
  mc_lower <- numeric(num_samples / step)
  mc_upper <- numeric(num_samples / step)
  qmc_lower <- numeric(num_samples / step)
  qmc_upper <- numeric(num_samples / step)

  # Compute CI bounds at each step
  for (i in seq(step, num_samples, by = step)) {
    mc_subset <- mc_entropies[1:i]
    qmc_subset <- qmc_entropies[1:i]

    mc_hpd <- compute_hpd_interval(mc_subset)
    qmc_hpd <- compute_hpd_interval(qmc_subset)

    mc_lower[i / step] <- mc_hpd[1]
    mc_upper[i / step] <- mc_hpd[2]
    qmc_lower[i / step] <- qmc_hpd[1]
    qmc_upper[i / step] <- qmc_hpd[2]
  }

  # Plot the CI bounds
  ci_steps <- seq(step, num_samples, by = step)
  plot(ci_steps, mc_lower, type = "l", col = "blue", lwd = 2,
       xlab = "Number of Samples", ylab = "CI Bounds",
       main = "MC vs QMC Confidence Interval Bounds",
       ylim = range(c(mc_lower, mc_upper, qmc_lower, qmc_upper)))
  lines(ci_steps, mc_upper, type = "l", col = "blue", lwd = 2, lty = 2)
  lines(ci_steps, qmc_lower, type = "l", col = "red", lwd = 2)
  lines(ci_steps, qmc_upper, type = "l", col = "red", lwd = 2, lty = 2)

  # Add legend
  legend("right", legend = c("MC Lower", "MC Upper", "QMC Lower", "QMC Upper"),
         col = c("blue", "blue", "red", "red"), lty = c(1, 2, 1, 2), lwd = 2)
  abline(h = mean(c(qmc_upper[length(qmc_upper)], mc_upper[length(mc_upper)])), col = "gray", lty = 2)
  abline(h = mean(c(qmc_lower[length(qmc_lower)], mc_lower[length(mc_lower)])), col = "gray", lty = 1)

  # Add timing as text
  text(x = num_samples * 0.5, y = max(mc_upper, qmc_upper),
       labels = sprintf("MC Time: %.2f sec\nQMC Time: %.2f sec",
                        timing_results$mc_time, timing_results$qmc_time),
       col = "black", cex = 0.8, pos = 4)
}

# Example usage
plot_ci_bounds_with_timing(alpha, num_samples, step)



