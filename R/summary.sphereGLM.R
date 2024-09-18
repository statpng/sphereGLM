#' @method summary sphereGLM
#' @export 
summary.sphereGLM <- function(fit, scale=FALSE, orthogonal=NULL){
  
  if( class(fit) == "try-error" ){
    return( list(beta=NA,
                 norm=NA, 
                 fisher=NA,
                 wald=NA,
                 pvalue=NA) )
  }
  
  if(is.null(orthogonal)){
    orthogonal <- fit$params$orthogonal
  }
  
  X <- fit$X
  Y <- fit$Y
  
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)
  
  if( scale ){
    beta_new <- fit$beta2
  } else {
    beta_new <- fit$beta
  }
  
  
  F_nj.list <- FisherInfoMatrix(X, beta_new)
  
  Wald.j.list <- sapply(1:(p+1), function(j){
    t(beta_new[j,]) %*% F_nj.list[[j]] %*% beta_new[j,]
  })
  
  
  if(orthogonal){
    df <- q-1
  } else {
    df <- q
  }
  
  
  pv.list <- sapply(Wald.j.list, function(x){
    1 - pchisq(x, df)
  })
  
  
  
  
  list(beta=beta_new,
       norm=apply(beta_new, 1, norm, "2"), 
       fisher=F_nj.list,
       orthogonal=fit$params$orthogonal,
       wald=Wald.j.list,
       df=df,
       pvalue=pv.list)
  
  
  
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

