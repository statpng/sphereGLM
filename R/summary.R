#' @export summary.sphereGLM
summary.sphereGLM <- function(fit){
  
  X <- fit$X
  Y <- fit$Y
  
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)
  
  beta2 <- fit$beta2
  
  F_nj.list <- FisherInfoMatrix(X, beta2)
  
  Wj.list <- sapply(1:(p+1), function(j){
    t(beta2[j,]) %*% F_nj.list[[j]] %*% beta2[j,]
  })
  
  
  pv.list <- sapply(Wj.list, function(x){
    1 - pchisq(x, q)
  })
  
  
  
  
  list(beta=beta2,
       norm=apply(beta2, 1, norm, "2"), 
       Wald=Wj.list,
       p.value=pv.list)
  
}





FisherInfoMatrix <- function(X, beta2, OffSet=NULL){
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(beta2)
  
  if(is.null(OffSet)){
    OffSet <- matrix(0, n, q)
  }
  
  bhat <- as.vector(beta2)
  
  X1 <- cbind(1, X)
  Xt.list <- apply(X1, 1, function(Xi) kronecker(diag(1,q), Xi), simplify=FALSE) # [(p+1)q x q] x n
  Xt <- do.call("cbind", Xt.list) # (p+1)q x qn
  b2.list <- lapply(1:n, function(i) b2.vMF( OffSet[i,] + crossprod(Xt.list[[i]], bhat) ))
  
  
  Fnj <- NULL
  for( j in 1:(p+1) ){
    out <- 0
    for( i in 1:nrow(X1) ){
      # Xj.list <- lapply(idx.list, function(jj){
      #   Xt.list[[i]][jj,]
      # })
      Xj <- do.call("rbind", lapply(1:q, function(jj){
        idx <- j + (jj-1)*(p+1)
        Xt.list[[i]][idx,]
      }))
      out <- out + Xj %*% b2.list[[i]] %*% t(Xj)
    }
    
    Fnj[[j]] <- out
  }
  
  Fnj
}

