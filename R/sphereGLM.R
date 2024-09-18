#' @importFrom manifold frechetMean
#' @importFrom manifold createM
#' @importFrom magic adiag
#' @importFrom dplyr `%>%`
#' @useDynLib sphereGLM
#' 
#' @export sphereGLM
sphereGLM <- function(X, Y, MU=NULL, orthogonal=FALSE, standardize=TRUE, Offset=NULL, maxit=100, eps=1e-6, lambda=5e-4, use.nlm=FALSE){
  
  
  
  if(FALSE){
    
    library(dplyr)
    set.seed(1)
    simdata <- sphereGLM::sim.sphereGLM(n=100, p=1, q=3, mu=c(0,100,0), snr=50, s=10, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1])
    
    system.time(fit0 <- sphereGLM(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=1e-10, orthogonal=F))
    system.time(fit1 <- sphereGLM.R(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=1e-6, orthogonal=F))
    
    fit0
    
    fit0 %>% class
    fit0 %>% summary
    fit0 %>% plot
    
    fit1 %>% class
    fit1 %>% summary
    fit1 %>% plot
    
    
    
    #
    
    
    fit0 %>% summary.sphereGLM()
    fit1 %>% summary.sphereGLM()
    
    fit0$beta
    fit1$beta
    
    plot.sphereGLM(fit0)
    plot.sphereGLM(fit1)
    
    
    
    res.benchmark <- microbenchmark::microbenchmark(
      fit0=sphereGLM(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=5e-4, orthogonal=T),
      fit1=sphereGLM.R(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=5e-4, orthogonal=T), times=10
    )
    res.benchmark
    
    
    
    
    
    set.seed(1)
    simdata <- sphereGLM::sim.sphereGLM(n=1000, p=1, q=3, mu=c(0,100,0), snr=50, s=10, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1])
    
    res.benchmark.1000 <- microbenchmark::microbenchmark(
      fit0=sphereGLM_cpp(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=1e-6, orthogonal=F), times=10
    )
    
    res.benchmark.1000
    #
    
    
    #
    
  }
  
  
  
  FrechetMean_eqiv <- function(X, ...){
    manifold::frechetMean(manifold::createM("Sphere"), t(X), ...)
  }
  
  diag.matrix <- function(X, lambda){
    n=dim(X)[1]; p=dim(X)[2]
    if(n != p) stop("The matrix is not square; it has unequal dimensions.")
    
    diag(lambda, n, p)
  }
  
  
  
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  
  
  X00 <- X
  X0 <- scale(X, center=TRUE, scale=FALSE)
  
  if(standardize){
    sdx.inv <- apply(X0,2,sd) %>% {diag(1/., length(.), length(.))}
    X <- X0 %*% sdx.inv
  } else {
    X <- X0
  }
  
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  
  
  {
    
    
    
    if(is.null(Offset)){
      Offset <- matrix(0,n,q)
    }
    
    
    if(is.null(MU)){
      MU <- FrechetMean_eqiv(Y)
      
      muhat <- colSums(Y) %>% {./norm(.,"2")} # mle of mean direction in vMF
      rbar <- colSums(Y) %>% {norm(.,"2")/nrow(Y)}
      
      kappa.type <- 4
      
      if(is.null(kappa.type)){
        kappa0 <- 1
      } else {
        kappa0 <- list(
          (q-1) / 2*(1-rbar),
          q*rbar*( 1+q/(q+2)*rbar^2 + q^2*(q+8)/((q+2)^2*(q+4))*rbar^4 ),
          (rbar*q) / (1-rbar^2),
          (rbar*q - rbar^3) / (1-rbar^2)
        )[[kappa.type]]
      }
      
      
      # kappa0 <- (q-1) / 2*(1-rbar) # approximation
      # kappa0 <- q*rbar*( 1+q/(q+2)*rbar^2 + q^2*(q+8)/((q+2)^2*(q+4))*rbar^4 ) # approximation
      # kappa0 <- (rbar*q) / (1-rbar^2) # approximation
      # kappa0 <- (rbar*q - rbar^3) / (1-rbar^2) # approximation
      
      MU <- muhat * kappa0
    }
    
    
    
    
    
    {
      B0 <- lm(Y ~ -1 + X, offset=tcrossprod(rep(1,n),MU)+Offset) %>% coef
      beta <- as.vector(t(rbind(MU, B0)))
    }
    
    
    
    
    
    X1 <- cbind(1, X)
    
    Xt.list <- apply(X1, 1, function(Xi) kronecker(Xi, diag(1,q)), simplify=FALSE) # [pq x q]
    Xt <- do.call("cbind", Xt.list) # pq x qn
    yt <- as.vector(t(Y))
    
    
    l1 <- 1
    beta_old <- beta
    beta_new <- beta + 1
    beta.list <- Fn.list <- loglik.list <- crit.list <- NULL
    crit <- 1
    
    
    fit0 <- sphereGLM_iteration(X=X, Y=Y, Offset=Offset, beta=beta, Xt=Xt, Xt_list=Xt.list, eps=eps, maxit=maxit, lambda=lambda, orthogonal=orthogonal)
    
    beta.list <- fit0$beta_list
    beta_new <- fit0$beta
    
    
    loglik.list <- fit0$loglik_list
    crit.list <- fit0$crit_list
    Fn.list <- fit0$Fn_list
  
    
    beta.list2 <- lapply(beta.list, function(bb) matrix(bb, p+1, q, byrow=TRUE) )
    beta_new <- matrix(beta_new, p+1, q, byrow=TRUE)
    
    if(standardize){
      beta.list2 <- lapply(beta.list2, function(b) magic::adiag(1, sdx.inv) %*% b)
      beta_new <- magic::adiag(1, sdx.inv) %*% beta_new
    }
    
    
    
    
    result <- list(X=X0, Y=Y, mu=beta_new[1,],
                   beta=beta_new,
                   beta2=beta_new,
                   beta.list=beta.list2,
                   offset=Offset,
                   loglik.list=loglik.list[-1], 
                   crit.list=crit.list[-1])  
    
  }
  
  
  
  
  params <- list( MU=MU, Offset=Offset, maxit=maxit, eps=eps, standardize=standardize, orthogonal=orthogonal, use.nlm=use.nlm )
  
  
  result$beta2[-1,] <- t( apply(result$beta2[-1,,drop=F], 1, function(x) x/norm(result$beta2[1,], "2")) )
  result$beta2[1,] <- result$beta2[1,]/norm(result$beta2[1,], "2")
  
  result$params <- params
  
  class(result) <- "sphereGLM"
  
  return(result)
  
}






