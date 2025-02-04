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





##################variance components estimation ####################
# Inputs-- r.hat: empirical residuals computed from y-X%*%beta.hat; z: design for random effects;
#  G.list: a list of basis functions for Psi_eta; grp: a vector of length N indicating group membership; alpha : tuning parameter in $\Sig_a$
# Outputs-- eta.hat: estimator of eta; sig.e.sq.hat: estimator of sig_e^2;
#  is.spd: is Psi_{eta.hat} positive definite?


Varcomp.est<- function(r.hat,z, G.list, alpha, grp){
  q=ncol(z)
  n = length(unique(grp))
  ugrp = unique(grp)
  z <- as.matrix(z)
  ###estimate sig^2.e## residual variance
  denom<-0
  num<-0
  for (lev in ugrp){
    cur.mem<-which(grp==lev)
    if(length(cur.mem) <= q) next
    fit.res <- lm(r.hat[cur.mem]~z[cur.mem,]-1)$res
    denom = denom + length(cur.mem)-q
    num <- num + sum(fit.res^2)
  }
  sig.e.sq.hat<- num/denom
  sig.e.sq.hat = max(sig.e.sq.hat, 10^(-4))
  ######estimate eta ## variance components
  G=G.list
  d<-length(G)
  B.a = matrix(0,nrow=d,ncol=d)
  omega.hat<-rep(0,d)
  for (lev in ugrp){
    cur.mem=which(grp==lev)
    zi = as.matrix(z[cur.mem, , drop = FALSE])
    ri<-r.hat[cur.mem]
    mi=length(cur.mem)
    cat("Group i =", lev, "; size =", mi, "\n")
    if (mi == 0) cat("WARNING: empty group found!\n")
    Sig.ai.inv = solve(alpha*(zi%*%t(zi))+ diag(1,mi)) #(Sig_a^i)^{-1}
    for (j in 1:d){
      for (k in 1:d){
        B.a[j,k] = B.a[j,k]+sum(diag(as.matrix(Sig.ai.inv%*%zi%*%G[[j]]%*%t(zi)%*%Sig.ai.inv%*%zi%*%G[[k]]%*%t(zi))))
      }
      
      omega.hat[j] = omega.hat[j] + t(ri)%*%Sig.ai.inv%*%zi%*%G[[j]]%*%t(zi)%*%Sig.ai.inv%*%ri-
        sig.e.sq.hat*sum(diag(Sig.ai.inv%*%zi%*%G[[j]]%*%t(zi)%*%Sig.ai.inv))
    }
  }
  eta.hat = solve(B.a)%*%omega.hat
  psi.hat<- sum(sapply(1:d, function(j) eta.hat[j]*G[[j]]))
  if(length(G)==1){
    eta.hat<-max(eta.hat,0)
  }
  is.spd<- (min(psi.hat)>=0)
  return(list(eta.hat = eta.hat, sig.e.sq.hat = sig.e.sq.hat, is.spd=is.spd))
}



