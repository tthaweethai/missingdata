# Inverse probability weighting for missing data
# Linear analysis model
# Logistic selection model

# Tanayott (Tony) Thaweethai
# Tuesday December 24, 2019
  
get.ipw.se <- function(dset,analysis.model.formula,selection.model.formula){
  
  selection.model.fit <- glm(selection.model.formula, data = dset, family = "binomial")
  dset$pi.hat <- predict(selection.model.fit, dset, type = "response")
  analysis.model.fit <- glm(analysis.model.formula, data = dset, weights = 1/pi.hat)
  
  beta.hat <- matrix(coef(analysis.model.fit),ncol=1)
  alpha.hat <- matrix(coef(selection.model.fit),ncol=1)
  
  N <- nrow(dset)
  
  R.var <- all.vars(selection.model.formula)[1]
  R.vec <- dset[,R.var]
  R.indices <- which(R.vec==1)
  
  Z.mat <- model.matrix(selection.model.fit)
  q <- ncol(Z.mat)
  
  X.mat.R <- model.matrix(analysis.model.formula,data=dset,na.action='na.pass')
  p <- ncol(X.mat.R)
  X.mat <- do.call("rbind",lapply(1:N,function(i){
    if (i %in% R.indices){
      return(X.mat.R[which(rownames(X.mat.R) == i),])
    } else {
      return(matrix(NA,nrow=1,ncol=p))
    }
  }))
  
  Y.var <- all.vars(analysis.model.formula)[1]
  Y.vec <- dset[,Y.var]
  pi.hat <- dset[,"pi.hat"]
  
  S.theta <- do.call("cbind",lapply(1:N,function(i){
    if (R.vec[i] == 0){
      S.theta.i <- matrix(0,nrow=p,ncol=1)
    } else {
      X.i <- matrix(X.mat[i,],nrow=1)
      Y.i <- Y.vec[i]
      pi.hat.i <- pi.hat[i]
      S.theta.i <- 1/pi.hat.i*t(X.i)%*%(Y.i-X.i%*%beta.hat)
    }
    return(S.theta.i)
  }))
  
  S.alpha <- do.call("cbind",lapply(1:N,function(i){
    Z.i <- matrix(Z.mat[i,],nrow=1)
    R.i <- R.vec[i]
    return(t(Z.i)%*%(R.i-plogis(Z.i%*%alpha.hat)))
  }))
  
  I.alpha.hat <- 1/N * Reduce("+",lapply(1:N,function(i){
    Z.i <- matrix(Z.mat[i,],nrow=1)
    R.i <- R.vec[i]
    return(t(Z.i)%*%((exp(Z.i%*%alpha.hat)/((1+exp(Z.i%*%alpha.hat))^2)))%*%Z.i)
  }))
  
  delta.hat <- 1/N * Reduce("+",lapply(1:N,function(i){
    if (R.vec[i] == 1){
      Z.i <- matrix(Z.mat[i,],nrow=1)
      X.i <- matrix(X.mat[i,],nrow=1)
      Y.i <- Y.vec[i]
      pi.hat.i <- pi.hat[i]
      return(as.numeric(1/exp(Z.i%*%alpha.hat))*t(X.i)%*%(Y.i-X.i%*%beta.hat)%*%Z.i)
    } else {
      return(matrix(0,nrow=p,ncol=q))
    }
  }))
  
  tau.hat <- 1/N * Reduce("+",lapply(1:N,function(i){
    if (R.vec[i] == 1){
      X.i <- matrix(X.mat[i,],nrow=1)
      pi.hat.i <- pi.hat[i]
      return(1/pi.hat.i*t(X.i)%*%X.i)
    } else {
      return(matrix(0,nrow=p,ncol=p))
    }
  }))
  
  I.alpha.hat.inv <- solve(I.alpha.hat)
  Omega.hat <- 1/N * Reduce("+",lapply(1:N,function(i){
    nu.i <- S.theta[,i] - delta.hat%*%I.alpha.hat.inv%*%S.alpha[,i]
    return(nu.i%*%t(nu.i))
  }))
  
  Sigma.hat <- solve(tau.hat) %*% Omega.hat %*% t(solve(tau.hat))
  se <- sqrt(diag(Sigma.hat)/N)

  return(c(beta.hat,se))
  
}

# Test dataset
# N <- 10^3
# beta <- c(1.25, 0.50, 1.25, -0.75)
# d <- data.frame(matrix(nrow=N,ncol=0))
# d$X <- rbinom(N,1,0.5)
# d$Y <- rnorm(N,beta[1]+beta[2]*d$X)
# d$R <- rbinom(N,1,plogis(beta[3]+beta[4]*d$Y))
# d.obs <- d
# d.obs$X <- ifelse(d.obs$R, d.obs$X, NA)
# get.ipw.se(dset=d.obs,
#            analysis.model.formula=Y~X,
#            selection.model.formula=R~Y)

