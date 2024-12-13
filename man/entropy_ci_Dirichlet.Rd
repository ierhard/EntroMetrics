% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropy_ci_Dirichlet.R
\name{entropy_ci_Dirichlet}
\alias{entropy_ci_Dirichlet}
\title{Credible interval for Shannon entropy using the Dirichlet posterior}
\usage{
entropy_ci_Dirichlet(
  bin_counts,
  Dir_alpha_prior = 1,
  Dir_alpha_posterior = NULL,
  N_MC = 5 * 10^3,
  unit = c("log2", "ln", "normalize"),
  pt_est_output = NULL
)
}
\arguments{
\item{bin_counts}{A numeric vector of observed bin counts. Unnecessary if \code{Dir_alpha_posterior} is provided.}

\item{Dir_alpha_prior}{A numeric value or vector representing the parameters of the Dirichlet prior distribution. Defaults to \code{1}, representing a uniform prior. Unnecessary if \code{Dir_alpha_posterior} is provided.}

\item{Dir_alpha_posterior}{Optional: A numeric vector specifying the parameters of the Dirichlet posterior distribution. If not provided, it is computed from the \code{bin_counts} and \code{Dir_alpha_prior}.}

\item{N_MC}{A numeric value indicating the number of Monte Carlo simulations to estimate the posterior entropy distribution. Defaults to \code{5000}.}

\item{unit}{A character string specifying the unit for entropy: \code{"log2"} for bits, \code{"ln"} for nats, or \code{"normalize"} for normalized entropy/evenness. Defaults to \code{"log2"}.}

\item{pt_est_output}{Optional: A list output from an entropy point estimation function (e.g., \code{entropy_estimator_Dirichlet}). If provided and contains \code{Dir_alpha_posterior}, the argument \code{Dir_alpha_posterior} will be updated.}
}
\value{
A list with the following components:
\item{ci}{A numeric vector containing the lower and upper bounds of the 95\% highest posterior density interval for Shannon entropy (i.e. the shortest interval containing 95\% of the probability mass of the posterior entropy distribution).}
\item{credible_mean}{The mean of the sampled entropy posterior distribution.}
\item{credible_median}{The median of the sampled entropy posterior distribution.}
\item{credible_cdf}{The estimated posterior entropy distribution as an empirical cumulative distribution function (ECDF) object.}
}
\description{
This function estimates a credible interval (the Bayesian counterpart of a confidence interval) for Shannon entropy based on the Dirichlet posterior distribution. The posterior is either derived from observed bin counts and a Dirichlet prior, or from a pre-calculated Dirichlet posterior distribution passed via \code{pt_est_outputs}. Monte Carlo simulations are used to estimate the posterior distribution for entropy.
}
\examples{
# Example usage with uniform prior
entropy_ci_Dirichlet(bin_counts = c(4, 1, 5))

# Example with a non-uniform Dirichlet prior
entropy_ci_Dirichlet(bin_counts = c(4, 1, 5), Dir_alpha_prior = c(1, 2, 1))

# Example with output from an entropy estimator
pt_est <- entropy_estimator_Dirichlet(c(4, 1, 5))
entropy_ci_Dirichlet(pt_est_output = pt_est)
}
