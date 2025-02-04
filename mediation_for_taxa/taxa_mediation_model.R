# mediation_with_multilevel Function 
#
# Inputs --
#   X        : Design matrix for fixed effects.
#   Y        : Response variable.
#   OTU      : Mediator variables matrix.
#   grp      : Vector indicating group membership for each observation.
#   COV      : Covariate variables matrix.
#   n_boot   : Number of bootstrap samples. Relevant only if method = "bootstrap".
#   n_cores  : Number of cores for parallel processing. Relevant only if method = "bootstrap".
# Outputs --
#   Direct effect, Total effect, indirect_effects : Data frame containing estimates with confidence interval and p-values.

mediation_with_multilevel <- function(X, 
                                      Y, 
                                      OTU, 
                                      grp,
                                      COV, 
                                      n_boot = 1000,
                                      n_cores = parallel::detectCores() - 1)
{
  X <- matrix(X, ncol = 1)
  print(table(X))
  M.clr <- as.matrix(OTU)
  COV <- as.matrix(COV)
  M_ID_name <- colnames(M.clr)
  
  if(is.null(M_ID_name)) M_ID_name <- paste0("otu", seq_len(ncol(M.clr)))
  
  if (any(apply(COV, 2, var) == 0)) {
    stop("Some covariates have zero variance.")
  }
  
  if (!is.factor(grp)) grp <- as.factor(grp)
  
  if (any(table(grp) < 2)) {
    stop("All groups must contain at least two observations.")
  }
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Run bootstrap in parallel
  bootstrap_results <- foreach(b = seq_len(n_boot), 
                               .packages = c('vegan', 'nlme', 'glmnet', 'scalreg'), 
                               .export = c('outcome_model_with_dlasso', 'optimal_a', 'Varcomp.est', 'BC.CI'),
                               .errorhandling = 'pass') %dopar% {
                                 
                                 # Resample clusters with replacement
                                 Clusters <- levels(grp)
                                 G <- length(Clusters)
                                 resample_successful <- FALSE
                                 attempt <- 1
                                 max_attempt <- 10
                                 
                                 while(!resample_successful && attempt <= max_attempt){
                                   sampled_clusters <- sample(Clusters, size = G, replace = TRUE)
                                   
                                   # Initialize lists to store data
                                   X_boot_list <- list()
                                   Y_boot_list <- list()
                                   M_clr_boot_list <- list()
                                   COV_boot_list <- list()
                                   grp_boot <- c()
                                   
                                   for (i in seq_along(sampled_clusters)) {
                                     cluster <- sampled_clusters[i]
                                     indices <- which(grp == cluster)
                                     X_boot_list[[i]] <- X[indices, , drop = FALSE]
                                     Y_boot_list[[i]] <- Y[indices]
                                     M_clr_boot_list[[i]] <- M.clr[indices, , drop = FALSE]
                                     COV_boot_list[[i]] <- COV[indices, , drop = FALSE]
                                     # Assign new unique group labels
                                     grp_boot <- c(grp_boot, rep(paste0("boot_grp_", i), length(indices)))
                                   }
                                   # Combine data
                                   X_boot <- do.call(rbind, X_boot_list)
                                   Y_boot <- unlist(Y_boot_list)
                                   M.clr_boot <- do.call(rbind, M_clr_boot_list)
                                   COV_boot <- do.call(rbind, COV_boot_list)
                                   grp_boot <- factor(grp_boot)
                                   
                                   # Check for singleton groups
                                   grp_counts <- table(grp_boot)
                                   if (any(grp_counts < 2)) {
                                     # Resample again
                                     attempt <- attempt + 1
                                     next
                                   } else {
                                     resample_successful <- TRUE
                                   }
                                   if (attempt == max_attempt && !resample_successful){
                                     warning(paste("Bootstrap sample", b, "failed after", max_attempt, "attempts."))
                                   }
                                 }
                                 
                                 if(!resample_successful){
                                   return(NULL)
                                 }
                                 
                                 # Apply scaling and centering
                                 MT_boot <- scale(M.clr_boot)
                                 colnames(MT_boot) <- colnames(M.clr_boot)
                                 COV.T_boot <- scale(COV_boot)
                                 colnames(COV.T_boot) <- colnames(COV_boot)
                                 MX_boot <- cbind(MT_boot, X_boot)
                                 MXCOV_boot <- cbind(MX_boot, COV.T_boot)
                                 Y.C_boot <- Y_boot - mean(Y_boot)
                                 z <- matrix(1, nrow = length(Y_boot), ncol = 1)
                                 exposure_col_boot <- ncol(MT_boot) + 1
                                 inf.coord_boot <- c(1:ncol(MT_boot), exposure_col_boot)
                                 # Define the basis matrices for random intercept
                                 #G.list <- list(diag(1, 1))
                                 
                                 # Find the optimum constant 'a'
                                 opt.a_boot <- optimal_a(X = MXCOV_boot, y = Y.C_boot, z = z,  grp = grp_boot, a.seq = seq(1, 100, 1), train_frac=0.80)
                                 print(paste("optimal a:", opt.a_boot))
                                 if(is.null(opt.a_boot) || is.na(opt.a_boot)){
                                   stop("optimal_a returned NULL or NA.")
                                 }
                                 
                                 # Fit the outcome model including all mediators with De-biased Lasso estimates
                                 fit.dlasso_boot  <- outcome_model_with_dlasso(X = MXCOV_boot, y = Y.C_boot, z = z, grp = grp_boot, a = opt.a_boot, inf.coord = inf.coord_boot)
                                 
                                 if(is.null(fit.dlasso_boot) || is.null(fit.dlasso_boot$beta.hat)){
                                   stop("outcome_model_with_dlasso failed or returned NULL.")
                                 }
                                 
                                 # Estimates variance components 
                                 #r.hat= Y.C_boot-MXCOV_boot%*%fit.dlasso_boot$beta.hat
                                 #Varcomp.est_boot <- Varcomp.est(r.hat = r.hat, z = z, G.list = G.list, alpha = opt.a_boot, grp = grp_boot)
                                 #print(paste("variance component:", Varcomp.est_boot$eta.hat))
                                 #if(is.null(Varcomp.est_boot)){
                                 #   stop("Varcomp.est failed.")
                                 #}
                                 
                                 beta_est_boot <- fit.dlasso_boot$beta.db
                                 beta_EST_boot <- beta_est_boot[1:ncol(MT_boot)]
                                 
                                 M_ID_name_boot <- colnames(MT_boot)
                                 n_mediators_boot <- length(M_ID_name_boot)
                                 
                                 # Calculate the direct effect (exposure variable)
                                 direct_effect_est_boot <- beta_est_boot[exposure_col_boot]
                                 
                                 # Fit mediator models for all mediators
                                 mediator_data_boot <- MT_boot
                                 alpha_est_boot <- numeric(n_mediators_boot)
                                 
                                 for(i in seq_len(n_mediators_boot)){
                                   mediator_values_boot <- mediator_data_boot[, i]
                                   lme_mediator_boot <- nlme::lme(
                                     mediator_values_boot ~ X_boot + COV.T_boot,
                                     random = (~1 | grp_boot), method = "REML",
                                     control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim", niterEM = 200))
                                   lme_mediator_summary_boot <- summary(lme_mediator_boot)
                                   alpha_est_boot[i] <- lme_mediator_summary_boot$tTable[2,1]
                                 }
                                 # Calculate indirect effects for all mediators
                                 indirect_effects_boot <- alpha_est_boot * beta_EST_boot
                                 
                                 # Calculate the total effect (Direct effect + Indirect effects)
                                 total_effect_est_boot <- direct_effect_est_boot + sum(indirect_effects_boot)
                                 
                                 # Return the estimates
                                 list(direct_effect = direct_effect_est_boot,
                                      total_effect = total_effect_est_boot,
                                      indirect_effects = indirect_effects_boot)
                               }
  
  # Stop the cluster
  stopCluster(cl)
  
  # Separate successful results and errors
  successes <- bootstrap_results[!sapply(bootstrap_results, inherits, "error") & !sapply(bootstrap_results, is.null)]
  errors <- bootstrap_results[sapply(bootstrap_results, inherits, "error")]
  
  print(length(successes))
  print(length(errors))
  
  # Check if we have enough bootstrap samples
  if (length(successes) < n_boot * 0.75) {
    warning("Less than 75% of bootstrap samples were successful. Results may be unreliable.")
  }
  
  if (length(successes) == 0) {
    stop("All bootstrap samples failed. Cannot compute estimates.")
  }
  
  # Proceed with processing successes
  n_successes <- length(successes)
  n_mediators <- ncol(M.clr)
  mediator_names <- M_ID_name
  
  bootstrap_direct_effects <- numeric(n_successes)
  bootstrap_total_effects <- numeric(n_successes)
  bootstrap_indirect_effects_matrix <- matrix(NA, nrow = n_successes, ncol = n_mediators)
  colnames(bootstrap_indirect_effects_matrix) <- mediator_names
  
  for (i in seq_along(successes)) {
    res <- successes[[i]]
    bootstrap_direct_effects[i] <- res$direct_effect
    bootstrap_total_effects[i] <- res$total_effect
    bootstrap_indirect_effects_matrix[i, ] <- res$indirect_effects
  }
  
  # Calculate estimates, confidence intervals, and p-values
  direct_CI_bca <- BC.CI(bootstrap_direct_effects, conf.level = 0.95)
  direct_effect_est <- mean(bootstrap_direct_effects, na.rm = TRUE)
  P_direct <- mean(bootstrap_direct_effects <= 0, na.rm = TRUE)
  P_direct <- 2 * min(P_direct, 1 - P_direct)
  
  total_effect_est <- mean(bootstrap_total_effects, na.rm = TRUE)
  total_CI_bca <- BC.CI(bootstrap_total_effects, conf.level = 0.95)
  P_total <- mean(bootstrap_total_effects <= 0, na.rm = TRUE)
  P_total <- 2 * min(P_total, 1 - P_total)
  
  indirect_effects_est <- colMeans(bootstrap_indirect_effects_matrix, na.rm = TRUE)
  indirect_CI_bca_list <- lapply(seq_len(ncol(bootstrap_indirect_effects_matrix)), function(i) {
    BC.CI(bootstrap_indirect_effects_matrix[, i], conf.level = 0.95)
  })
  indirect_CI_lower <- sapply(indirect_CI_bca_list, `[`, 1)
  indirect_CI_upper <- sapply(indirect_CI_bca_list, `[`, 2)
  P_indirect <- apply(bootstrap_indirect_effects_matrix, 2, function(x) {
    x <- x[!is.na(x)]
    p_value <- mean(x <= 0)
    2 * min(p_value, 1 - p_value)
  })
  
  # Prepare the results
  indirect_effects_df <- data.frame(
    Mediators = mediator_names,
    Estimate = indirect_effects_est,
    CI_Lower = indirect_CI_lower,
    CI_Upper = indirect_CI_upper,
    p_value = P_indirect,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  # Adjust p-values for multiple comparisons
  indirect_effects_df$p_adj <- p.adjust(indirect_effects_df$p_value, method = "BH")
  
  # Direct effect
  direct_effect_df <- data.frame(
    Estimate = direct_effect_est,
    CI_Lower = direct_CI_bca[1],
    CI_Upper = direct_CI_bca[2],
    p_value = P_direct,
    row.names = NULL
  )
  
  # Total effect
  total_effect_df <- data.frame(
    Estimate = total_effect_est,
    CI_Lower = total_CI_bca[1],
    CI_Upper = total_CI_bca[2],
    p_value = P_total,
    row.names = NULL
  )
  
  # Prepare the final results
  results <- list(
    indirect_effects = indirect_effects_df,
    direct_effect = direct_effect_df,
    total_effect = total_effect_df
  )
  return(results)
}

