#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bootstrap BCa
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bootstrap.bca <- function(estimator_fct, data_vec, B = 10^3, alpha = 0.05) {
  
  # Remove NA's
  data_vec <- data_vec[!is.na(data_vec)]
  
  # Sample size
  n <- length(data_vec)
  
  # Calculate point estimate
  theta_hat <- estimator_fct(data_vec)

  # Generate bootstrap estimates
  
  theta_boot <- rep(NA, B)
  
  for(i in 1:B){
    
    # Select a bootstrap sample
    boot_sample <- sample(data_vec, length(data_vec), replace=TRUE)
    
    # Save bootstrap estimate
    theta_boot[i] <- estimator_fct(boot_sample)
    
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
    xnew <- data_vec[-i]
    
    # Estimate theta
    theta_jack <- estimator_fct(xnew)
    
    I[i] <- (n-1)*(theta_hat - theta_jack)
    
  }
  
  # Estimate a
  a_hat <- (sum(I^3) / sum(I^2)^1.5)/6
  
  # Adjusted quantiles 
  u_adjusted <- pnorm(z0 + (z0+zu)/(1-a_hat*(z0+zu))) 
  
  # Accelerated Bootstrap CI
  
  ci <- quantile(theta_boot, u_adjusted, na.rm = TRUE) %>% as.numeric()
  
  return(ci)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bootstrap BCa confidence intervals for many confidence levels
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


bootstrap.bca.2arg.df <- function(estimator_fct, data_vec1, data_vec2,
                                  B = 10^3, alpha_vec = (1:(10^3-1))/10^3) {
  
  # # For testing
  # estimator_fct <- spearman_cor
  # data_vec1 <- cluster_data$TM1_SQ001_ordered
  # data_vec2 <- cluster_data$TM1_SQ003_ordered
  # B = 10^3
  # alpha_vec = (1:(10^3-1))/10^3
  
  # Remove NA's
  keep_ind <- (!is.na(data_vec1)) & (!is.na(data_vec2))
  data_vec1 <- data_vec1[keep_ind]
  data_vec2 <- data_vec2[keep_ind]
  
  # Sample size
  n <- length(data_vec1)
  
  # Calculate point estimate
  theta_hat <- estimator_fct(data_vec1, data_vec2)
  
  # Generate bootstrap estimates
  
  theta_boot <- rep(NA, B)
  
  for(i in 1:B){
    
    # Select a bootstrap sample
    boot_ind <- sample(1:n, replace=TRUE)
    
    boot_v1 <- data_vec1[boot_ind]
    boot_v2 <- data_vec2[boot_ind]
    
    # Save bootstrap estimate
    theta_boot[i] <- estimator_fct(boot_v1, boot_v2)
    
  }
  
  # List of desired quantiles
  u_list <- map(alpha_vec, function(alpha) c(alpha/2, 1 - alpha/2))
  
  # Compute constants
  z0 <- qnorm(mean(theta_boot <= theta_hat))
  
  zu_list <- map(u_list, qnorm)
  
  # Acceleration constant
  I <- rep(NA, n)
  
  for (i in 1:n) {
    
    # Remove ith data point
    data_vec1_new <- data_vec1[-i]
    data_vec2_new <- data_vec2[-i]
    
    # Estimate theta
    theta_jack <- estimator_fct(data_vec1_new, data_vec2_new)
    
    I[i] <- (n-1)*(theta_hat - theta_jack)
    
  }
  
  # Estimate a
  a_hat <- (sum(I^3) / sum(I^2)^1.5)/6
  
  # Calculate all CIs
  ci_df <- data.frame(alpha = numeric(),
                      ci.lower = numeric(),
                      ci.upper = numeric())
  
  for (i in 1:length(alpha_vec)) {
    
    zu <- zu_list[[i]]
    
    # Adjusted quantiles
    
    u_adjusted <- pnorm(z0 + (z0 + zu) / (1 - a_hat*(z0 + zu))) 
    
    # Accelerated Bootstrap CI
    
    ci_temp <- quantile(theta_boot, u_adjusted, na.rm = TRUE) %>% as.numeric()
    
    row_df <- data.frame(alpha = alpha_vec[i],
                         ci.lower = ci_temp[1],
                         ci.upper = ci_temp[2])
    
    ci_df <- bind_rows(ci_df, row_df)
    
  }
  
  return(ci_df)
  
}



