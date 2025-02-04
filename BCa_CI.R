
################## The bias-corrected and accelerated (BCa) bootstrap interval ####################
BC.CI <- function(theta, conf.level = 0.95) {
  # Remove NA values
  theta <- theta[!is.na(theta)]
  alpha <- (1 - conf.level)/2
  low <- alpha
  high <- 1 - alpha
  
  # Compute bias correction factor z0
  # z0 is the quantile of the standard normal corresponding to the proportion of
  # bootstrap estimates less than the observed mean.
  obs <- mean(theta) # using mean as the point estimate (you could use median or original estimate)
  z.inv <- length(theta[theta < obs]) / length(theta)
  z0 <- qnorm(z.inv)
  
  # Compute acceleration factor a using the skewness of the distribution
  U <- (length(theta)-1)*(obs - theta)
  top <- sum(U^3)
  under <- (1/6)*(sum(U^2))^(3/2)
  a <- top/under
  
  # Adjusted quantiles for lower bound
  alpha1 <- pnorm(z0 + (z0 + qnorm(low)) / (1 - a*(z0 + qnorm(low))))
  # Adjusted quantiles for upper bound
  alpha2 <- pnorm(z0 + (z0 + qnorm(high)) / (1 - a*(z0 + qnorm(high))))
  
  # Compute BCa confidence interval
  bca.lower <- quantile(theta, alpha1, na.rm = TRUE)
  bca.upper <- quantile(theta, alpha2, na.rm = TRUE)
  return(c(bca.lower, bca.upper))
}

