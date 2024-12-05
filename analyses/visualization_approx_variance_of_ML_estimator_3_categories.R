# Plot approximate variance of ML entropy estimator for 3-category distribution


# Define z(p) function
z <- function(p1, p2) {
  p3 <- 1 - p1 - p2
  if (p3 < 0 || p1 < 0 || p2 < 0 || p3 > 1) return(NA)
  sum_term <- p1 * (log(p1)^2) + p2 * (log(p2)^2) + p3 * (log(p3)^2)
  entropy_term <- -(p1 * log(p1) + p2 * log(p2) + p3 * log(p3))
  return(sum_term - entropy_term)
}

# Generate grid for p1 and p2
p1_vals <- seq(0, 1, length.out = 100)
p2_vals <- seq(0, 1, length.out = 100)
z_vals <- outer(p1_vals, p2_vals, Vectorize(function(p1, p2) z(p1, p2)))

# Mask invalid values where p1 + p2 > 1
z_vals[outer(p1_vals, p2_vals, "+") > 1] <- NA

# Plot the 3D surface
library(rgl)
persp3d(p1_vals, p2_vals, z_vals, xlab = "p1", ylab = "p2", zlab = "z(p)", col = "lightblue")

################################################################################

# Plot approximate variance of ML entropy estimator vs entropy
# (for 3-category distribution)

N_fineness <- 200

p1 <- seq(0, 1, length.out = N_fineness)
p2 <- seq(0, 1, length.out = N_fineness)

data_df <- data.frame(p1 = NULL, p2 = NULL, var_approx = NULL, entropy = NULL)

counter <- 0
for(i in 1:length(p1)){
  for(j in 1:length(p2)){
    p <- c(p1[i], p2[j], 1 - p1[i] - p2[j])
    if (p[1] + p[2] <= 1){
      var_approx <- sum(p[p > 0] * log(p[p > 0])^2) - sum(p[p > 0] * log(p[p > 0]))
      entropy <- -sum(p[p > 0] * log(p[p > 0]))
      data_df <- rbind(data_df,
                       data.frame(p1 = p[1],
                                  p2 = p[2],
                                  var_approx = var_approx,
                                  entropy = entropy))
    }
    counter <- counter + 1
    print(paste0("Progress: ", counter/(N_fineness^2) * 100, "%"))
  }
}

plot(data_df$entropy, data_df$var_approx, xlab = "Entropy", ylab = "Approximate Variance")

