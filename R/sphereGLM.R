#' @importFrom manifold frechetMean
#' @importFrom manifold createM
#' @importFrom magic adiag
#' @importFrom dplyr `%>%`
#' @useDynLib sphereGLM
#' 
#' @export sphereGLM
sphereGLM <- function(X, Y, MU=NULL, orthogonal=FALSE, penalty.factor=rep(0,ncol(X)), standardize=TRUE, Offset=NULL, maxit=100, eps=1e-4, lambda=1e-4, use.nlm=FALSE, gamma=9999){
  
  if(FALSE){
    MU=NULL; orthogonal=FALSE; standardize=TRUE; Offset=NULL; maxit=100; eps=1e-4; lambda=1e-4; use.nlm=FALSE; gamma=0
    penalty.factor=rep(0,ncol(X))
  }
  
  
  
  if(FALSE){
    devtools::load_all()
    
    simdata <- sim.sphereGLM(n=100, p=2, q=3, mu=c(0,0,20), s=10, s0=0, type="vMF", seed.UDV=1, qr.U=FALSE, qr.V=FALSE)
    
    res <- microbenchmark::microbenchmark(
      fit.nlm = with(simdata, sphereGLM.R(X=X, Y=Y, orthogonal=FALSE, lambda=1e-10, use.nlm = TRUE)),
      fit.R = with(simdata, sphereGLM.R(X=X, Y=Y, orthogonal=FALSE, lambda=1e-10, use.nlm = FALSE)),
      fit.nonortho = try(with(simdata, sphereGLM(X=X, Y=Y, orthogonal=FALSE, eps=1e-8))),
      fit.ortho = with(simdata, sphereGLM(X=X, Y=Y, orthogonal=TRUE, eps=1e-8)),
      times=10
    )

    res
    
    # Unit: milliseconds
    # expr       min        lq      mean    median        uq       max neval
    # fit.nlm 806.07033 811.06639 826.83488 815.26673 817.06096 940.69174    10
    # fit.R 388.73970 400.26500 421.86182 412.03739 416.28784 541.79122    10
    # fit.nonortho  29.77621  30.86558  33.97734  33.82582  37.11218  38.96271    10
    # fit.ortho  21.77768  26.24324  27.73491  27.33673  29.47982  31.81202    10
    # 
  }
  
  
  if(FALSE){
    
    library(dplyr)
    set.seed(1)
    simdata <- sphereGLM::sim.sphereGLM(n=100, p=1, q=3, mu=c(0,100,0), snr=50, s=50, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1])
    
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
  

  if(is.null(Offset)){
    Offset <- matrix(0,n,q)
  }
  
  
  
  {
    
    if(is.null(MU)){
      MU <- FrechetMean(Y)
      
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
    
    
    fit0 <- sphereGLM_iteration(X=X, Y=Y, Offset=Offset, beta=beta, Xt=Xt, Xt_list=Xt.list, eps=eps, maxit=maxit, lambda=lambda, orthogonal=orthogonal, gamma=gamma, zero_beta = which(penalty.factor==1))
    
    
    
    beta.list <- lapply( fit0$beta_list, function(b) ifelse(abs(b) < 1e-10, 0, b) )
    beta_new <- ifelse(abs(fit0$beta) < 1e-10, 0, fit0$beta)
    
    
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
  
  
  
  
  params <- list( MU=MU, Offset=Offset, maxit=maxit, eps=eps, standardize=standardize, orthogonal=orthogonal, lambda=lambda, gamma=gamma, use.nlm=use.nlm )
  
  
  result$beta2[-1,] <- t( apply(result$beta2[-1,,drop=F], 1, function(x) x/norm(result$beta2[1,], "2")) )
  result$beta2[1,] <- result$beta2[1,]/norm(result$beta2[1,], "2")
  
  result$params <- params
  
  class(result) <- "sphereGLM"
  
  return(result)
  
}













#' @export sphereGLM.R
sphereGLM.R <- function(X, Y, MU=NULL, orthogonal=FALSE, standardize=TRUE, Offset=NULL, maxit=100, eps=1e-6, lambda=5e-4, use.nlm=FALSE, gamma=9999){
  
  # Figure: How good is the initial value?
  if(FALSE){
    
    set.seed(1)
    simdata <- sim.sphereGLM(n=150, p=2, q=3, mu=c(0,0,1), s=5, s0 = 0.01)
    
    # true theta
    plot.sphere(simdata$Theta, col="black")
    
    # initial estimate for theta
    plot.sphere(simdata %>% {.$X %*% solve(crossprod(.$X), crossprod(.$X, .$Y))}, add=TRUE, col="red")
    
    # sphereGLM estimate for theta
    fit <- sphereGLM(X=simdata$X, Y=simdata$Y)
    plot.sphere(cbind(1, fit$X) %*% fit$beta2, add=TRUE, col="blue")
    
  }
  
  
  
  # For Plotting
  if(FALSE){
    set.seed(1)
    simdata <- sim.sphereGLM(n=150, mu=c(0,0,1), p=2, s=5, s0 = 0.01)
    X <- simdata$X;  Y <- simdata$Y
    fit <- sphereGLM(X[,1:2], Y)
    plot.sphereGLM(fit, plot.mu = T)
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
  
  
  
  if(use.nlm){
    
    
    x <- model.matrix(~., data.frame(X))
    y <- Y
    
    
    
    regvmf <- function(be, y, x) {
      
      be <- matrix(be, ncol = q)
      theta <- x %*% be
      ki <- sqrt(Rfast::rowsums(theta^2))
      
      # Cq <- function(NORM, q){
      #   (  (q/2-1) * log(NORM) - ( (q/2) * log(2*pi) + log( besselI(NORM, q/2-1, expon.scaled=TRUE) ) + NORM )  )
      # }
      # - sum( Cq(ki, q) ) - sum(theta * y)
      
      - sum( apply(theta, 1, Cq, log=TRUE) ) - sum(theta * y)
      # -sum(  log(ki)) + sum(log(sinh(ki))) - sum(theta * y)
    }
    
    
    ini <- solve(crossprod(x), crossprod(x, y))
    # val1 <- nlm(regvmf, ini, y = y, x = x, iterlim = 1000)
    # val2 <- nlm(regvmf, val1$estimate, y = y, x = x, iterlim = 1000)
    
    
    
    
    suppressWarnings({
      val1 <- nlm(regvmf, ini, y = y, x = x, iterlim = maxit)
      val2 <- nlm(regvmf, val1$estimate, y = y, x = x, iterlim = maxit)
      # while(val1$minimum - val2$minimum > eps) {
      #   val1 <- val2
      #   val2 <- nlm(regvmf, val1$estimate, y = y, x = x, 
      #               iterlim = 1000)
      # }
      da <- optim(val2$estimate, regvmf, y = y, x = x, 
                  control = list(maxit = maxit), 
                  hessian = TRUE)
    })
    
    
    
    result <- NULL
    result$X <- X0
    result$Y <- Y
    result$mu <- matrix(da$par, ncol = q)[1,]
    result$beta <- matrix(da$par, ncol = q)
    if(standardize){
      result$beta <- magic::adiag(1, sdx.inv) %*% result$beta
    }
    result$beta2 <- result$beta
    result$loglik <- -da$value
    
    
    
  } else {
    
    
    
    if(is.null(Offset)){
      Offset <- matrix(0,n,q)
    }
    
    
    if(is.null(MU)){
      MU <- FrechetMean(Y)
      
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
    
    Xt.list <- apply(X1, 1, function(Xi) kronecker(Xi, diag(1,q)), simplify=FALSE) # [pq x q] x 100
    Xt <- do.call("cbind", Xt.list) # pq x qn
    yt <- as.vector(t(Y))
    
    
    l1 <- 1
    beta_old <- beta
    beta_new <- beta + 1
    beta.list <- Fn.list <- loglik.list <- crit.list <- NULL
    crit <- 1
    while( crit > eps & l1 <= maxit ){
      
      beta.list[[l1]] <- beta
      
      beta_old <- beta
      
      (b1.list <- lapply(1:n, function(i) b1.vMF( Offset[i,] + crossprod(Xt.list[[i]], beta) )))
      eta <- do.call("rbind", b1.list)
      
      (b2.list <- lapply(1:n, function(i) b2.vMF( Offset[i,] + crossprod(Xt.list[[i]], beta) )))
      W <- do.call("adiag", b2.list) # magic::adiag
      # W.inv <- do.call("adiag", lapply(b2.list, function(x) solve(x)) ) # magic::adiag
      
      # (Z <- crossprod(Xt, beta) + chol2inv(chol(W)) %*% (yt - eta))
      Z <- crossprod(Xt, beta) + solve(W + diag.matrix(W, lambda)) %*% (yt - eta)
      
      # (beta <- chol2inv(chol(Xt %*% W %*% t(Xt))) %*% Xt %*% W %*% Z)
      
      mu <- beta[1:q]
      XWX <- Xt %*% W %*% t(Xt)
      
      if(orthogonal){
        
        beta <- solve( XWX +  gamma * Matrix::bdiag(append(list(matrix(0,q,q)), replicate(p, tcrossprod(mu), simplify=FALSE))) ) %*% Xt %*% W %*% Z
        
      } else {
        
        beta <- solve(XWX + Matrix::bdiag( 0, diag.matrix(XWX[-1,-1], lambda) ) ) %*% Xt %*% W %*% Z
        
      }
      
      beta_new <- beta
      
      (loglik <- Reduce("+", lapply(1:n, function(i){
        theta <- Offset[i,] + crossprod(Xt.list[[i]], beta)
        crossprod(theta, Y[i,]) + Cq(theta, log=TRUE)
      })))
      
      crit <- norm(beta_old-beta_new, "2")
      
      loglik.list[l1] <- loglik
      crit.list[l1] <- crit
      Fn.list[[l1]] <- XWX
      
      
      l1 <- l1 + 1
    }
    
    # print( matrix(beta,p,q) )
    # print( l1 )
    
    
    # plot(loglik.list, type="l")
    
    
    
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







diag.matrix <- function(X, lambda){
  n=dim(X)[1]; p=dim(X)[2]
  if(n != p) stop("The matrix is not square; it has unequal dimensions.")
  
  diag(lambda, n, p)
}





Cq <- function(theta, logarithm=TRUE){
  
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  
  # (NORM)^(q/2-1) / ( (2*pi)^(q/2) * besselI(NORM, q/2-1) )
  
  if(logarithm){
    
    (q/2-1) * log(NORM) - ( (q/2) * log(2*pi) + log( besselI(NORM, q/2-1, expon.scaled=TRUE) ) + NORM )
    
  } else {
    
    exp(  (q/2-1) * log(NORM) - ( (q/2) * log(2*pi) + log( besselI(NORM, q/2-1, expon.scaled=TRUE) ) + NORM )  )
    
  }
  
}




Bq <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  # besselI(NORM, nu=q/2, expon.scaled=FALSE) == besselI(NORM, nu=q/2, expon.scaled=TRUE)/exp(-NORM)
  # besselI(NORM, nu=q/2-1, expon.scaled=FALSE) == besselI(NORM, nu=q/2-1, expon.scaled=TRUE)/exp(-NORM)
  
  besselI(NORM, nu=q/2, expon.scaled=TRUE) / besselI(NORM, nu=q/2-1, expon.scaled=TRUE)
}



Hq <- function(theta){
  
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  I0 <- besselI(NORM, nu=q/2-1, expon.scaled=TRUE)
  I1 <- besselI(NORM, nu=q/2, expon.scaled=TRUE)
  
  1 - (I1/I0)^2 - (q-1)/NORM * I1/I0
  # numerator <- I0^2 - I1^2 - (q-1)/(NORM) * I0 * I1
  # denumerator <- I0^2
  # 
  # numerator / denumerator
}




subgrad <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  if(NORM == 0){
    v <- runif(q,-1,1)
    v / (norm(v, "2")*1.1)
  } else {
    theta / NORM
  }
}



#' @export b1.vMF
b1.vMF <- function(theta){
  Bq(theta) * subgrad(theta)
}



#' @export b2.vMF
b2.vMF <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  s <- subgrad(theta)
  
  # if( NORM < 1e-10 ){
  #   tcrossprod(s) * Hq(theta) + Bq(theta)
  # } else {
  tcrossprod( theta/NORM ) * Hq(theta) + Bq(theta) / NORM * (diag(1,q) - tcrossprod(theta/NORM) )
  # }
  
}


