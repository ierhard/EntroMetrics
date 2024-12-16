#' @name bootstrap_bca
#'
#' @title Bootstrap bias-corrected and accelerated (BCa) confidence interval
#'
#' @description
#' Computes the bias-corrected and accelerated (BCa) bootstrap confidence interval for a given point estimate.
#'
#' @param x_vec A numeric vector containing the data. Missing values (`NA`) will be removed.
#' @param pt_est_fct A function to compute the point estimate of interest. This function should take `x_vec` as input and return a single numeric value.
#' @param conf_level A numeric value specifying the confidence level for the interval (default is \eqn{0.95}).
#' @param B An integer specifying the number of bootstrap samples to generate (default is 1000).
#'
#' @return A list with the following components:
#' \item{ci}{A numeric vector of length 2 containing the lower and upper bounds of the confidence interval.}
#' \item{pt_est}{The point estimate computed from the input data.}
#'
#' @examples
#' # Define a function for the sample median
#' median_func <- function(x) median(x, na.rm = TRUE)
#'
#' # Generate sample data from an exponential distribution
#' set.seed(123)
#' data <- rexp(100, rate = 1) # Exponential with rate = 1
#'
#' # Compute BCa confidence interval for the median
#' bootstrap_bca(data, pt_est_fct = median_func, conf_level = 0.95, B = 1000)
#'
#' # Note: True median of the Exponential(rate = lambda) distribution is log(2)/lambda
#' #       Thus, the true median here is log(2) = 0.6931472
#'
#' @export
bootstrap_bca <- function(x_vec,
                          pt_est_fct,
                          conf_level = 0.95,
                          B = 10^3){

  # Compute alpha
  alpha <- 1 - conf_level

  # Remove NA's
  x_vec <- x_vec[!is.na(x_vec)]

  # Sample size
  n <- length(x_vec)

  # Calculate point estimate
  theta_hat <- pt_est_fct(x_vec)

  # Generate bootstrap estimates
  theta_boot <- rep(NA, B)

  for(i in 1:B){

    # Select a bootstrap sample
    boot_sample <- sample(x_vec, length(x_vec), replace=TRUE)

    # Save bootstrap estimate
    theta_boot[i] <- pt_est_fct(boot_sample)

  }

  # Desired quantiles
  u <- c(alpha/2, 1-alpha/2)

  # Compute constants
  z0 <- qnorm(mean(theta_boot <= theta_hat))

  zu <- qnorm(u)

  # Acceleration constant

  I <- rep(NA, n)

  for (i in 1:n) {

    # Remove ith data point
    xnew <- x_vec[-i]

    # Estimate theta
    theta_jack <- pt_est_fct(xnew)

    I[i] <- (n-1)*(theta_hat - theta_jack)

  }

  # Estimate acceleration constant a
  a_hat <- (sum(I^3) / sum(I^2)^1.5)/6

  # Adjusted quantiles
  u_adjusted <- pnorm(z0 + (z0+zu)/(1-a_hat*(z0+zu)))

  # Accelerated Bootstrap CI
  ci <- quantile(theta_boot, u_adjusted, na.rm = TRUE) %>% as.numeric()

  names(ci) <- c("ci_lower", "ci_upper")

  # Return results
  results_list <- list(ci = ci,
                       pt_est = theta_hat)

  return(results_list)

}
