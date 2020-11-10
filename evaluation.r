est <- function(H, X, Zs) {
  Sigmainvs <- list()
  for(j in 1:ncol(X)) {
    x <- X[,j]
    n <- length(x)
    D <- diag(x)
    R <- t(H)%*%H
    Rinv <- solve(R)
    Rx <- t(H)%*%D%*%H
    Sigma <- R-Rx%*%Rinv%*%Rx
    Sigmainvs[[j]] <- solve(Sigma)
  }
  zs <- matrix(0, 0, ncol(X))
  for(i in 1:nrow(Zs)) {
      z <- Zs[i, ]
      zs <- rbind(zs, sapply(Sigmainvs, function(U) t(z)%*%U%*%z))
      
  }
  return(apply(zs, 1, mean))
}



est_lp <- function(H, X, Zs) {
 p <- ncol(H)
 ss <- ncol(X)
 max.obj <- rep(999, ss)
 max.obj.approx <- rep(999, ss)
 if(p <= 50) {
  for(j in 1:ss) {
    x <- X[,j]
    n <- length(x)
    D <- diag(x)
    R <- t(H)%*%H
    Rinv <- solve(R)
    Rx <- t(H)%*%D%*%H
    Sigma <- R-Rx%*%Rinv%*%Rx
    Sigmainv <- try(solve(Sigma), TRUE) 
    if(is.matrix(Sigmainv)) {
     zast <- lp_sdp(Sigmainv)
     max.obj[j] <- t(zast)%*% Sigmainv %*% zast
    } 
    Sigma.approx <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv
    zast <- lp_sdp(Sigma.approx)
    max.obj.approx[j] <- t(zast)%*% Sigma.approx %*% zast
  }
 re <- data.frame(obj=max.obj, obj.approx=max.obj.approx)
 }
 if(p > 50) {
    Sigmainvs <- list()
    Sigma.approx <- list()
    for(j in 1:ncol(X)) {
      x <- X[,j]
      n <- length(x)
      D <- diag(x)
      R <- t(H)%*%H
      Rinv <- solve(R)
      Rx <- t(H)%*%D%*%H
      Sigma <- R-Rx%*%Rinv%*%Rx
      Sigmainvs[[j]] <- solve(Sigma)
      Sigma.approx[[j]] <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv
   }
   zs <- matrix(0, 0, ncol(X))
   for(i in 1:nrow(Zs)) {
     z <- Zs[i, ]
     zs <- rbind(zs, sapply(Sigmainvs, function(U) t(z)%*%U%*%z))
   }
   max.obj <- apply(zs, 2, max)
   zs <- matrix(0, 0, ncol(X))
   for(i in 1:nrow(Zs)) {
     z <- Zs[i, ]
     zs <- rbind(zs, sapply(Sigma.approx, function(U) t(z)%*%U%*%z))
   }
   max.obj.approx <- apply(zs, 2, max)
   re <- data.frame(obj=max.obj, obj.approx=max.obj.approx)
 }
 return(re)
}