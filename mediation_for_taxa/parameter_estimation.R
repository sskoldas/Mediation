# optimal_a Function
#
# Inputs --
#   X        : Design matrix for fixed effects predictors.
#   y        : Response variable vector.
#   z        : Design matrix for random effects covariates.
#   grp      : Group membership factor vector.
#   a.seq    : (Optional) Sequence of tuning parameter 'a' values to evaluate.
#   train_frac: (Optional) Fraction of unique groups to include in the training set.
# Outputs --
#   best.a   : Optimal tuning parameter 'a' that minimizes the prediction error.

optimal_a <- function(X, y, z, grp, 
                      a.seq = seq(1, 200, 1),
                      train_frac=0.9) {
  best.err <- sum(y^2)
  best.a <- a.seq[1]
  
  data <- data.frame(Y = y, GRP = grp, X, Z = z, stringsAsFactors = FALSE)
  data$GRP <- as.factor(data$GRP)
  
  unique_grps <- unique(data$GRP)
  shuffle_grps <- sample(unique_grps)
  
  n_tr <- ceiling(train_frac * length(shuffle_grps))
  train_grp <- shuffle_grps[1:n_tr]
  test_grp <- shuffle_grps[(n_tr + 1):length(shuffle_grps)]
  
  train_data <- dplyr::filter(data, GRP %in% train_grp)
  test_data <- dplyr::filter(data, GRP %in% test_grp)
  
  n_X <- ncol(X)
  n_z <- ncol(z)
  
  X_train <- as.matrix(train_data[, 3:(2+n_X)])
  Y_train <- train_data$Y
  Z_train <- as.matrix(train_data[, (3+n_X):(2+n_X+n_z)])
  
  X_test <- as.matrix(test_data[, 3:(2+n_X)])
  Y_test <- test_data$Y
  
  for (i in seq_along(a.seq)) {
    a_val <- a.seq[i]
    est.re <- outcome_model_with_dlasso(X_train, Y_train, Z_train, train_data$GRP, a = a_val, inf.coord=NULL)
    predictions <- X_test %*% est.re$beta.hat
    pred.err <- sum((Y_test - predictions)^2)
    if (pred.err < best.err) {
      best.err <- pred.err
      best.a <- a_val
    }
  }
  
  return(best.a)
}

