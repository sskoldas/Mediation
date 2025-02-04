# mediation_analysis_for_alpha_diversity Function 
#
# Inputs --
#   X               : Exposure variable matrix.
#   Y               : Response variable.
#   alpha_diversity : Alpha Diversity variable matrix.
#   grp             : Vector indicating group membership for each observation.
#   COV             : Covariate variable matrix.
#   n_boot          : Number of bootstrap samples. Relevant only if method = "bootstrap".
#   n_cores         : Number of cores for parallel processing. Relevant only if method = "bootstrap".
# Outputs --
#   Direct effect, Total effect, indirect_effect : Data frame containing estimates with confidence interval and p-values.



# Community-Level Mediation
mediation_analysis_for_alpha_diversity <- function( 
    X, 
    Y, 
    COV,
    grp,
    alpha_diversity,
    n_boot = 1000,
    n_cores = parallel::detectCores()-1) {
  
  X <- matrix(X, ncol = 1)
  COV <- as.matrix(COV)
  if (!is.factor(grp)) grp <- as.factor(grp)
  alpha_diversity <- as.matrix(alpha_diversity)
  
  if (any(table(grp) < 2)) {
    stop("All groups must contain at least two observations.")
  }
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  bootstrap_results <- foreach(
    b=seq_len(n_boot),
    .packages = c('nlme'),
    .errorhandling = 'pass') %dopar% {
      
      Clusters <- levels(grp)
      G <- length(Clusters)
      resample_successful <- FALSE
      attempt <- 1
      max_attempt <- 10
      while(!resample_successful && attempt <= max_attempt){
        sampled_clusters <- sample(Clusters, size = G, replace = TRUE)
        X_boot_list <- list()
        Y_boot_list <- list()
        alphaD_boot_list <- list()
        COV_boot_list <- list()
        grp_boot <- c()
        
        for(i in seq_along(sampled_clusters)){
          cluster <- sampled_clusters[i]
          indices <- which(grp == cluster)
          X_boot_list[[i]] <- X[indices, , drop =FALSE]
          Y_boot_list[[i]] <- Y[indices]
          alphaD_boot_list[[i]] <- alpha_diversity[indices, , drop=FALSE]
          COV_boot_list[[i]] <- COV[indices, , drop=FALSE]
          grp_boot <- c(grp_boot, rep(paste0("boot_grp_",i), length(indices)))
        }
        X_boot <- do.call(rbind, X_boot_list)
        Y_boot <- unlist(Y_boot_list)
        alphaD_boot <- do.call(rbind, alphaD_boot_list)
        COV_boot <- do.call(rbind, COV_boot_list)
        grp_boot <- factor(grp_boot)
        
        grp_counts <- table(grp_boot)
        if (any(grp_counts < 2)) {
          # Resample again
          attempt <- attempt + 1
          next
        } else {
          resample_successful <- TRUE
        }
      }
      if(!resample_successful){
        return(NULL)
      }
      
      alphaD_boot <- scale(alphaD_boot)
      COV_boot <- scale(COV_boot)
      Y_boot <- Y_boot - mean(Y_boot)
      
      # Outcome Model
      outcome_model <- nlme::lme(Y_boot ~ X_boot + alphaD_boot + COV_boot,
                                 random = (~1 | grp_boot), method = "REML",
                                 control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim", niterEM = 100))
      outcome_summary <- summary(outcome_model)
      
      # Mediator Model
      mediator_model <- nlme::lme(alphaD_boot ~ X_boot + COV_boot, 
                                  random = (~1 | grp_boot), method = "REML",
                                  control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim", niterEM = 100))
      mediator_summary <- summary(mediator_model)
      
      # Indirect Effect
      alpha_est <- outcome_summary$tTable[3, 1]
      beta_est <- mediator_summary$tTable[2, 1]
      indirect_est_boot <- alpha_est * beta_est
      
      # Direct effect
      direct_est_boot <- outcome_summary$tTable[2, 1]
      
      # Total Effect
      total_est_boot <- indirect_est_boot + direct_est_boot
      
      list(indirect_effect = indirect_est_boot,
           direct_effect = direct_est_boot,
           total_effect = total_est_boot)
    }
  stopCluster(cl)
  
  successes <- bootstrap_results[!sapply(bootstrap_results, inherits, "error") & !sapply(bootstrap_results, is.null)]
  errors <- bootstrap_results[sapply(bootstrap_results, inherits, "error")]
  
  print(length(successes))
  print(length(errors))
  if (length(successes) < n_boot * 0.75) {
    warning("Less than 50% of bootstrap samples were successful. Results may be unreliable.")
  }
  
  if (length(successes) == 0) {
    stop("All bootstrap samples failed. Cannot compute estimates.")
  }
  
  n_successes <- length(successes)
  
  bootstrap_direct_effects <- numeric(n_successes)
  bootstrap_total_effects <- numeric(n_successes)
  bootstrap_indirect_effects <- numeric(n_successes)
  
  for (i in seq_along(successes)) {
    res <- successes[[i]]
    bootstrap_direct_effects[i] <- res$direct_effect
    bootstrap_total_effects[i] <- res$total_effect
    bootstrap_indirect_effects[i] <- res$indirect_effect
  }
  
  direct_effect_est <- mean(bootstrap_direct_effects, na.rm = TRUE)
  direct_CI_bca <- BC.CI(bootstrap_direct_effects, conf.level = 0.95)
  P_direct <- mean(bootstrap_direct_effects <= 0, na.rm = TRUE)
  P_direct <- 2 * min(P_direct, 1 - P_direct)
  
  total_effect_est <- mean(bootstrap_total_effects, na.rm = TRUE)
  total_CI_bca <- BC.CI(bootstrap_total_effects, conf.level = 0.95)
  P_total <- mean(bootstrap_total_effects <= 0, na.rm = TRUE)
  P_total <- 2 * min(P_total, 1 - P_total)
  
  indirect_effect_est <- mean(bootstrap_indirect_effects, na.rm = TRUE)
  indirect_CI_bca <- BC.CI(bootstrap_indirect_effects, conf.level = 0.95) 
  P_indirect <- mean(bootstrap_indirect_effects <= 0, na.rm = TRUE)
  P_indirect <- 2 * min(P_indirect, 1 - P_indirect)
  
  results <- data.frame(
    Effects = c("Total Effect", "Direct Effect", "Indirect Effect"),
    Estimate = c(total_effect_est, direct_effect_est, indirect_effect_est),
    CI_Lower = c(total_CI_bca[1], direct_CI_bca[1], indirect_CI_bca[1]),
    CI_Upper = c(total_CI_bca[2], direct_CI_bca[2], indirect_CI_bca[2]),
    p_value = c(P_total, P_direct, P_indirect),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  results[, sapply(results, is.numeric)] <- lapply(
    results[, sapply(results, is.numeric)],
    round,
    digits = 3
  )
  return(results)
}
