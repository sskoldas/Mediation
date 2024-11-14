# mediation_with_multilevel Function 
#
# Inputs --
#   X        : Design matrix for fixed effects.
#   Y        : Response variable.
#   OTU      : Mediator variables matrix.
#   grp      : Vector indicating group membership for each observation.
#   COV      : Covariate variables matrix.
#   method   : Method for mediation analysis. Options are "normal" (default) or "bootstrap".
#   n_boot   : Number of bootstrap samples. Relevant only if method = "bootstrap".
#   seed     : Random seed for reproducibility.
#   n_cores  : Number of cores for parallel processing. Relevant only if method = "bootstrap".
# Outputs --
#   Direct effect, Total effect, indirect_effects : Data frame containing estimates with confidence interval and p-values.

mediation_with_multilevel <- function(X, 
                                      Y, 
                                      OTU, 
                                      grp,
                                      COV, 
                                      method = "bootstrap",
                                      n_boot = 1000,
                                      seed = 1243,
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
  
  if (method == "normal") {
    
    set.seed(seed)
    d <- ncol(M.clr)
    Y.C <- Y - mean(Y)
    MT <- scale(M.clr)
    COV.T<- scale(COV)
    colnames(COV.T) <- colnames(COV)
    colnames(MT) <- colnames(M.clr)
    MX <- cbind(MT, X)
    MXCOV <- cbind(MX, COV.T)
    
    exposure_col <- ncol(MT) + 1
    inf.coord <- c(1:ncol(MT), exposure_col)
    
    # Find the optimum constant 'a'
    opt.a <- optimal_a(X = MXCOV, y = Y.C, z = X,  grp = grp, a.seq = seq(1, 100, 1), train_frac=0.95)
    
    # Fit the outcome model including all mediators with De-biased Lasso estimates
    fit.dlasso  <- outcome_model_with_dlasso(X = MXCOV, y = Y.C, z = X,  grp = grp, a = opt.a, inf.coord = inf.coord)
    
    beta_est <- fit.dlasso$beta.db
    beta_se  <- fit.dlasso$beta.db.sd
    beta_EST <- beta_est[1:ncol(MT)]
    beta_SE  <- beta_se[1:ncol(MT)]
    
    M_ID_name_all <- M_ID_name
    n_mediators <- length(M_ID_name_all)
    
    # Calculate the direct effect (exposure variable)
    direct_effect_est <- beta_est[exposure_col]
    direct_effect_se <- beta_se[exposure_col]
    direct_CI_lower <- direct_effect_est - 1.96*direct_effect_se
    direct_CI_upper <- direct_effect_est + 1.96*direct_effect_se
    P_direct <- 2*(1-pnorm(abs(direct_effect_est/direct_effect_se)))
    
    # Fit mediator models for all mediators
    alpha_est <- numeric(n_mediators)
    alpha_se <- numeric(n_mediators)
    
    mediator.data <- MT[, M_ID_name_all, drop = FALSE]
    
    for(i in seq_len(n_mediators)){
      mediator_values <- mediator.data[, i]
      lme_mediator <- nlme::lme(
        mediator_values ~ X + COV.T,
        random = (~1 + X | grp), method = "REML",
        control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim", niterEM = 200))
      lme_mediator_summary <- summary(lme_mediator)
      alpha_est[i] <- lme_mediator_summary$tTable[2,1]
      alpha_se[i] <- lme_mediator_summary$tTable[2,2]
    }
    
    # Calculate indirect effects for all mediators
    indirect_effects_est <- alpha_est * beta_EST
    indirect_se <- sqrt((beta_EST^2) * (alpha_se^2) + (alpha_est^2) * (beta_SE^2))
    z_indirect <- indirect_effects_est / indirect_se
    P_indirect <- 2 * (1 - pnorm(abs(z_indirect)))
    indirect_CI_lower <- indirect_effects_est - 1.96 * indirect_se
    indirect_CI_upper <- indirect_effects_est + 1.96 * indirect_se
    P_indirect_adj <- p.adjust(P_indirect, method = "BH")
    
    # Calculate the total effect by fitting Y ~ X (without mediators)
    lme_total_effect <- nlme::lme(Y.C ~ X + COV.T, 
                                  random = (~1 + X | grp), method = "REML",
                                  control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim", niterEM = 200))
    lme_total_summary <- summary(lme_total_effect)
    total_effect_est <- lme_total_summary$tTable[2,1]
    total_effect_se <- lme_total_summary$tTable[2,2]
    total_CI_lower <- total_effect_est - 1.96 * total_effect_se
    total_CI_upper <- total_effect_est + 1.96 * total_effect_se
    P_total <- 2 * (1 - pnorm(abs(total_effect_est / total_effect_se)))
    
    # Output results including all mediators
    results <- list(
      indirect_effects = data.frame(
        Mediators = M_ID_name_all,
        Estimate = indirect_effects_est,
        StdErr = indirect_se,
        CI_Lower = indirect_CI_lower,
        CI_Upper = indirect_CI_upper,
        p_value = P_indirect,
        p_adj = P_indirect_adj,
        row.names = NULL
      ),
      direct_effect = data.frame(
        Estimate = direct_effect_est,
        StdErr = direct_effect_se,
        CI_Lower = direct_CI_lower,
        CI_Upper = direct_CI_upper,
        p_value = P_direct,
        row.names = NULL
      ),
      total_effect = data.frame(
        Estimate = total_effect_est,
        StdErr = total_effect_se,
        CI_Lower = total_CI_lower,
        CI_Upper = total_CI_upper,
        p_value = P_total,
        row.names = NULL
      )
    )
  } else if (method == "bootstrap") {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Run bootstrap in parallel
    bootstrap_results <- foreach(b = seq_len(n_boot), 
                                 .packages = c('vegan', 'nlme', 'glmnet', 'scalreg'), 
                                 .export = c('outcome_model_with_dlasso', 'optimal_a'),
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
                                   
                                   exposure_col_boot <- ncol(MT_boot) + 1
                                   inf.coord_boot <- c(1:ncol(MT_boot), exposure_col_boot)
                                   
                                   # Find the optimum constant 'a'
                                   opt.a_boot <- optimal_a(X = MXCOV_boot, y = Y.C_boot, z = X_boot,  grp = grp_boot, a.seq = seq(1, 100, 1), train_frac=0.95)
                                   
                                   # Fit the outcome model including all mediators with De-biased Lasso estimates
                                   fit.dlasso_boot  <- outcome_model_with_dlasso(X = MXCOV_boot, y = Y.C_boot, z = X_boot, grp = grp_boot, a = opt.a_boot, inf.coord = inf.coord_boot)
                                   
                                   beta_est_boot <- fit.dlasso_boot$beta.db
                                   beta_EST_boot <- beta_est_boot[1:ncol(MT_boot)]
                                   
                                   M_ID_name_boot <- colnames(MT_boot)
                                   n_mediators_boot <- length(M_ID_name_boot)
                                   
                                   # Calculate the direct effect (exposure variable)
                                   direct_effect_est_boot <- beta_est_boot[exposure_col_boot]
                                   
                                   # Fit the total effect model Y ~ X (without mediators)
                                   lme_total_effect_boot <- nlme::lme(Y.C_boot ~ X_boot + COV.T_boot, 
                                                                      random = (~1 + X_boot | grp_boot), method = "REML",
                                                                      control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim", niterEM = 200))
                                   lme_total_summary_boot <- summary(lme_total_effect_boot)
                                   total_effect_est_boot <- lme_total_summary_boot$tTable[2,1]
                                   
                                   # Fit mediator models for all mediators
                                   mediator_data_boot <- MT_boot
                                   alpha_est_boot <- numeric(n_mediators_boot)
                                   
                                   for(i in seq_len(n_mediators_boot)){
                                     mediator_values_boot <- mediator_data_boot[, i]
                                     lme_mediator_boot <- nlme::lme(
                                       mediator_values_boot ~ X_boot + COV.T_boot,
                                       random = (~1 + X_boot | grp_boot), method = "REML",
                                       control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim", niterEM = 200))
                                     lme_mediator_summary_boot <- summary(lme_mediator_boot)
                                     alpha_est_boot[i] <- lme_mediator_summary_boot$tTable[2,1]
                                   }
                                   # Calculate indirect effects for all mediators
                                   indirect_effects_boot <- alpha_est_boot * beta_EST_boot
                                   
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
    if (length(successes) < n_boot * 0.5) {
      warning("Less than 50% of bootstrap samples were successful. Results may be unreliable.")
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
    direct_effect_est <- mean(bootstrap_direct_effects, na.rm = TRUE)
    direct_CI_lower <- quantile(bootstrap_direct_effects, probs = 0.025, na.rm = TRUE)
    direct_CI_upper <- quantile(bootstrap_direct_effects, probs = 0.975, na.rm = TRUE)
    P_direct <- mean(bootstrap_direct_effects <= 0, na.rm = TRUE)
    P_direct <- 2 * min(P_direct, 1 - P_direct)
    
    total_effect_est <- mean(bootstrap_total_effects, na.rm = TRUE)
    total_CI_lower <- quantile(bootstrap_total_effects, probs = 0.025, na.rm = TRUE)
    total_CI_upper <- quantile(bootstrap_total_effects, probs = 0.975, na.rm = TRUE)
    P_total <- mean(bootstrap_total_effects <= 0, na.rm = TRUE)
    P_total <- 2 * min(P_total, 1 - P_total)
    
    indirect_effects_est <- colMeans(bootstrap_indirect_effects_matrix, na.rm = TRUE)
    indirect_CI_lower <- apply(bootstrap_indirect_effects_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
    indirect_CI_upper <- apply(bootstrap_indirect_effects_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
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
      CI_Lower = direct_CI_lower,
      CI_Upper = direct_CI_upper,
      p_value = P_direct,
      row.names = NULL
    )
    
    # Total effect
    total_effect_df <- data.frame(
      Estimate = total_effect_est,
      CI_Lower = total_CI_lower,
      CI_Upper = total_CI_upper,
      p_value = P_total,
      row.names = NULL
    )
    
    # Prepare the final results
    results <- list(
      indirect_effects = indirect_effects_df,
      direct_effect = direct_effect_df,
      total_effect = total_effect_df
    )
  } else {
    stop("Invalid method specified. Choose 'normal' or 'bootstrap'.")
  }
  return(results)
}
