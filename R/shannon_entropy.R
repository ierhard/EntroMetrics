#' Shannon entropy estimation
#'
#' @param bin_counts Vector of bin counts
#' @param total_bins_known Logical vector indicating if number of bins known
#' @param point_estimator.method Point estimator
#' @param confidence_interval.method Confidence interval estimator
#' @param confidence_level Confidence level (default = 0.95)
#' @param B Number of bootstrap re-samplings (if applicable)
#' @param normalize Logical vector for normalization (TRUE -> normalized entropy)
#'
#' @return Point estimate and confidence interval
#'
#' @export
#'
#' @examples
#' shannon_entropy(c(100,9,0), TRUE, "ML", "Wald", normalize = TRUE)
shannon_entropy <- function(bin_counts,
                            total_bins_known = TRUE,
                            point_estimator.method = c(
                              "ML", # Max. likelihood
                              "MM", # Miller-Madow
                              "HS"  # Hausser-Strimmer
                            ),
                            confidence_interval.method = c(
                              "Wald",                     # Normality assumption
                              "Asym.Norm_quadratic",      # Normality assumption -> solve quadratic (Erhard)
                              "pct_boot",                 # percent bootstrap
                              "multinom_ci_based_Wallis", # Sean Wallis's (pairs with ML-estimator)
                              "multinom_ci_based_Erhard"  # Erhard (pairs with ML-estimator)
                            ),
                            confidence_level = 0.95,
                            B = 10^3,
                            normalize = TRUE){ # Number of bootstrap re-samplings (if CI method is "pct_boot")

  # Point estimate -------------------------------------------------------------

  # Number of categories
  K <- length(bin_counts)

  # Maximum entropy (in bits) for discrete distribution over K categories
  H_max <- log2(K)

  H_estimator <- function(bin_counts){

    # Entropy estimate (in bits)
    if (point_estimator.method == "HS"){
      H_hat <- entropy::entropy(y = bin_counts, method = "shrink", unit = "log2", verbose = FALSE)
    } else {
      H_hat <- entropy::entropy(y = bin_counts, method = point_estimator.method, unit = "log2")
    }

    # Ensure estimate lies between 0 and H_max
    H_hat <- max(min(H_hat, H_max), 0)

    return(H_hat)

  }

  H_hat <- H_estimator(bin_counts)

  # Confidence interval --------------------------------------------------------

  # Sample size
  n <- sum(bin_counts)

  alpha <- 1 - confidence_level

  # Wald CI (based on asymptotic normality)
  if(confidence_interval.method == "Wald"){

    # Relevant z-score
    z <- stats::qnorm(1-alpha/2)

    # Standard error estimate

    p_hat <- entropy::freqs(y = bin_counts, method = "ML")

    nu_hat <- sum(p_hat[p_hat > 0] * log2(p_hat[p_hat > 0])^2)

    H.empirical <- entropy::entropy.empirical(bin_counts, unit = "log2")

    std_error <- sqrt((nu_hat - H.empirical^2) / n) %>% max(0)

    # confidence interval

    ci.lower <- (H_hat - z*std_error) %>% min(H_max) %>% max(0)

    ci.upper <- (H_hat + z*std_error) %>% min(H_max) %>% max(0)

  }

  # Wald CI (based on asymptotic normality)
  if(confidence_interval.method == "Asym.Norm_quadratic"){

    # Relevant z-score
    z <- stats::qnorm(1-alpha/2)

    # Standard error estimate

    p_hat <- entropy::freqs(y = bin_counts, method = "ML")

    nu_hat <- sum(p_hat[p_hat > 0] * log2(p_hat[p_hat > 0])^2)

    H.empirical <- entropy::entropy.empirical(bin_counts, unit = "log2")

    a <- (1 + z^2 / n)

    b <- -2*H.empirical

    c <- H.empirical^2 - z^2*nu_hat

    roots <- c( (- b - sqrt( b^2 - 4*a*c ))/(2*a),
                (- b + sqrt( b^2 - 4*a*c ))/(2*a) )

    # confidence interval

    ci.lower <- min(roots)

    ci.upper <- max(roots)

  }

  # Percentile bootstrap CI
  if(confidence_interval.method == "pct_boot"){

    # Data for resampling
    data_vec <- rep(1:K, bin_counts)

    # Perform bootstraps
    bootstrap_distribution <- numeric(B)

    bootstrap_distribution[1] <- H_hat

    for(b in 2:B){
      data_boot <- sample(data_vec, replace = TRUE)

      bin_counts_boot <- purrr::map_int(1:K, \(i) sum(data_boot == i))

      bootstrap_distribution[b] <- H_estimator(bin_counts_boot)
    }

    # Bootstrap distribution
    bootstrap_distribution <- sort(bootstrap_distribution)

    # Confidence interval

    lower_index <- (1 + (B-1)*alpha/2) %>% floor()

    upper_index <- (1 + (B-1)*(1-alpha/2)) %>% ceiling()

    ci.lower <- bootstrap_distribution[lower_index]

    ci.upper <- bootstrap_distribution[upper_index]

  }

  # Normalize (if necessary) ---------------------------------------------------

  if(normalize == FALSE){
    result <- data.frame(entropy = H_hat,
                         ci.lower = ci.lower,
                         ci.upper = ci.upper)
  } else {
    result <- data.frame(normalized.entropy = H_hat / H_max,
                         ci.lower = ci.lower / H_max,
                         ci.upper = ci.upper / H_max)
  }


  # Return result: Point estimate and CI
  return(result)

}
