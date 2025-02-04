# outcome_model_with_dlasso function
#
# Inputs-- 
#  X: design for fixed effects; y: response; z: design for random effects;
#  grp: a vector of length N indicating group membership; a : tuning parameter
#  inf.coord: the coordinates of beta for inference.
# Outputs-- 
#  beta.hat: the Lasso estimator based on psedo-likelihood; 
#  beta.db:debiased Lasso estimators for the fixed effects in inf.coord;
#  beta.db.sd: standard deviation of the debiased Lasso estimators for the fixed effects in inf.coord.

outcome_model_with_dlasso <- function(X, y, z, grp, a=10, inf.coord=NULL){
  
  X <- as.matrix(X)
  N=length(y)
  cluster_ids = unique(grp)
  n=length(cluster_ids)
  q=ncol(z)
  p=ncol(X)
  
  #mixed effect model fitting
  X.a<-X
  y.a<-y
  for (i in 1:n){
    cur.grp = cluster_ids[i]
    cur.mem=which(grp == cur.grp)
    mi=length(cur.mem)
    
    zi = as.matrix(z[cur.mem,])
    sigmai = a*zi%*%t(zi)+diag(rep(1,mi))
    # Compute the inverse square root of the covariance matrix
    sigmai.svd <- svd(sigmai)
    Sig.a.inv.half <- sigmai.svd$u %*% diag(1/sqrt(sigmai.svd$d)) %*% t(sigmai.svd$u)
    X.a[cur.mem,] =  Sig.a.inv.half %*% as.matrix(X[cur.mem,]) 
    y.a[cur.mem] =  Sig.a.inv.half %*% as.matrix(y[cur.mem])
    
  }
  cv.init<-cv.glmnet(X.a, y.a, lambda=seq(5, 0.1, -0.1)*sqrt(2*log(p)/N))
  beta.hat<-coef(cv.init, s=cv.init$lambda.min)[-1]
  
  #Debiased Lasso for Mixed-Effect Model with Cluster-Robust Standard Errors
  if(is.null(inf.coord)){
    return( list(beta.hat = beta.hat))
  }
  res<-y.a-X.a %*% beta.hat
  lam.seq<-seq(5, 0.1, -0.1)*sqrt(2*log(p)/N)
  beta.db.sd.mlm<-rep(NA,length=length(inf.coord))
  beta.db.mlm<-rep(NA,length=length(inf.coord))
  
  for(j in 1:length(inf.coord)){
    col.j<-inf.coord[j]
    # Fit Lasso regression of X[, col.j] on other predictors
    cv.x<-cv.glmnet(X.a[,-col.j], X.a[,col.j], lambda=lam.seq)
    gam.j<- coef(cv.x, s=cv.x$lambda.min)[-1]
    wj.mlm <- X.a[,col.j]-  X.a[,-col.j]%*%gam.j
    denom <- sum(wj.mlm * X.a[, col.j])
    beta.db.mlm[j] = beta.hat[col.j] + sum( wj.mlm * res)/(denom)
    
    num = 0
    for( i in 1:n){
      cur.grp = cluster_ids[i]
      cur.mem <- which(grp == cur.grp)
      num <- num + sum(wj.mlm[cur.mem]*res[cur.mem])^2
    }
    beta.db.sd.mlm[j] <- sqrt(num) / sum(wj.mlm * X.a[,col.j])
    
  }
  return(list(beta.hat = beta.hat, beta.db=beta.db.mlm, beta.db.sd=beta.db.sd.mlm))
}







