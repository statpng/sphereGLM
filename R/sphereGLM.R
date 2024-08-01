#' @export sphereGLM
sphereGLM <- function(X, Y, MU=NULL, Offset=NULL, lambda1=1e-12, lambda2=1e-12, beta0=NULL, maxit=100, eps=1e-8, standardize=TRUE, kappa.type=4, use.nlm=TRUE){
  
  
  # Figure: Initial value가 얼마나 괜찮은것인지?
  if(FALSE){
    
    # true theta
    png.sphere(simdata$Theta, col="black")
    
    # initial estimate for theta
    png.sphere(x %*% solve(crossprod(x), crossprod(x, y)), add=TRUE, col="red")
    
    # sphereGLM estimate for theta
    fit <- sphereGLM(X=simdata$X, Y=simdata$Y)
    png.sphere(cbind(1, fit$X) %*% fit$beta2, add=TRUE, col="blue")
    
  }
  
  
  
  if(FALSE){
    simdata <- sim.sphereGLM()
    X=simdata$X; Y=simdata$Y; MU=NULL; Offset=NULL; lambda1=1e-12; lambda2=1e-12; 
    beta0=NULL; maxit=100; eps=1e-8; standardize=TRUE; kappa.type=4; use.nlm=TRUE
  }
  
  
  # For Plotting
  if(FALSE){
    set.seed(1)
    simdata <- sim.sphere.v2(n=150, mu=c(0,1,0), r=2, s=5, s0 = 0.01)
    X <- simdata$U;  Y <- simdata$X
    fit <- glm.vmf(X[,1:2], Y)
    glm.vmf.plot(fit, plot.mu = T)
  }
  
  
  
  if(FALSE){
    library(png.Directional)
    
    Y <- rvmf(150, rnorm(3), 5)
    X <- iris[, 3:4]
    
    vmf.reg2(Y, X)$beta
    glm.vmf(as.matrix(X), Y, use.nlm=TRUE)$beta
    glm.vmf(as.matrix(X), Y, use.nlm=FALSE)$beta
    
    
    microbenchmark::microbenchmark(
      vmf.reg2(Y, X)$beta,
      glm.vmf(as.matrix(X), Y)$beta,
      glm.vmf(as.matrix(X), Y, use.nlm=FALSE)$beta,
      times=10
    )
    
    
    
    
  }
  
  
  if(FALSE){
    ii <- pvec.list[[d]]
    X=U.joint[,-1]
    Y=X_joint[,ii]
    Offset=Theta.inds[,ii]
    
    beta0=NULL; maxit=1000; eps=1e-8; standardize=TRUE
  }
  
  
  
  if(use.nlm){
    
    x <- model.matrix(~., data.frame(X))
    y <- Y
    n <- dim(Y)[1]
    q <- ncol(Y)
    p <- ncol(X)
    
    
    
    regvmf <- function(be, y, x) {
      be <- matrix(be, ncol = q)
      mu <- x %*% be
      ki <- sqrt(Rfast::rowsums(mu^2))
      
      - sum( apply(mu, 1, Cq, log=TRUE) ) - sum(mu * y)
      # -sum(  log(ki)) + sum(log(sinh(ki))) - sum(mu * y)
    }
    
    
    # regvmf <- function(be, y, x) {
    #   be <- matrix(be, ncol = q)
    #   mu <- x %*% be
    #   ki <- sqrt(Rfast::rowsums(mu^2))
    # 
    #   - sum( apply(mu, 1, Cq, log=TRUE) ) - sum(mu * y)
    #   # -sum(  log(ki)) + sum(log(sinh(ki))) - sum(mu * y)
    # }
    
    
    ini <- solve(crossprod(x), crossprod(x, y))
    val1 <- nlm(regvmf, ini, y = y, x = x, iterlim = 1000)
    # val2 <- nlm(regvmf, val1$estimate, y = y, x = x, iterlim = 1000)
    
    result <- NULL
    result$X <- X
    result$Y <- Y
    result$beta <- matrix(val1$estimate, ncol = q)
    result$beta2 <- result$beta
    result$beta2[-1,] <- t( apply(result$beta[-1,,drop=F], 1, function(x) x/norm(result$beta[1,], "2")) )
    result$beta2[1,] <- result$beta[1,]/norm(result$beta[1,], "2")
    result$loglik <- -val1$minimum
    
    
  } else {
    
    
    # start_time <- proc.time()[3]
    # tt <- 0
    # 
    # 
    # part_time <- proc.time()[3]
    # tt <- tt+1
    # cat("Part ", tt, " time:", part_time - start_time, "\n")
    # start_time <- part_time
    
    
    
    
    
    
    library(magic)
    
    n <- nrow(Y)
    q <- ncol(Y)
    p <- ncol(X)
    
    if(is.null(Offset)){
      Offset <- matrix(0,n,q)
    }
    
    X0 <- scale(X, center=TRUE, scale=FALSE)
    
    if(standardize){
      sdx.inv <- apply(X0,2,sd) %>% {diag(1/., length(.), length(.))}
      X <- X0 %*% sdx.inv
    } else {
      X <- X0
    }
    
    
    
    
    if(is.null(MU)){
      MU <- FrechetMean(Y)
      
      muhat <- colSums(Y) %>% {./norm(.,"2")} # mle of mean direction in vMF
      rbar <- colSums(Y) %>% {norm(.,"2")/nrow(Y)}
      
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
      beta <- as.vector(rbind(MU, B0))
    }
    
    
    
    
    
    X1 <- cbind(1, X)
    
    Xt.list <- apply(X1, 1, function(Xi) kronecker(diag(1,q), Xi), simplify=FALSE) # [pq x q] x 100
    Xt <- do.call("cbind", Xt.list) # pq x qn
    yt <- as.vector(t(Y))
    
    
    l1 <- 1
    beta_old <- beta
    beta_new <- beta + 1
    beta.list <- loglik.list <- crit.list <- NULL
    crit <- 1
    while( crit > eps & l1 <= maxit ){
      
      beta.list[[l1]] <- beta
      
      beta_old <- beta
      
      (b1.list <- lapply(1:n, function(i) b1.vMF( Offset[i,] + crossprod(Xt.list[[i]], beta) )))
      eta <- do.call("rbind", b1.list)
      
      (b2.list <- lapply(1:n, function(i) b2.vMF( Offset[i,] + crossprod(Xt.list[[i]], beta) )))
      W <- do.call("adiag", b2.list) # magic::adiag
      # W.inv <- do.call("adiag", lapply(b2.list, function(x) solve(x)) ) # magic::adiag
      
      (Z <- crossprod(Xt, beta) + chol2inv(chol(W)) %*% (yt - eta))
      # Z <- crossprod(Xt, beta) + solve(W + diag(lambda1, n*q, n*q)) %*% (yt - eta)
      
      # beta <- solve(Xt %*% W %*% t(Xt) + diag(c(0,rep(lambda2,p)),q*(p+1),q*(p+1))) %*% Xt %*% W %*% Z
      (beta <- chol2inv(chol(Xt %*% W %*% t(Xt))) %*% Xt %*% W %*% Z)
      
      # beta <- solve(Xt %*% W %*% t(Xt) + diag(lambda2,q*(p+1),q*(p+1))) %*% Xt %*% W %*% Z
      
      # beta[1:3] <- beta[1:3] %>% {./norm(.,"2")}
      
      beta_new <- beta
      
      (loglik <- Reduce("+", lapply(1:n, function(i){
        theta <- Offset[i,] + crossprod(Xt.list[[i]], beta)
        crossprod(theta, Y[i,]) + Cq(theta, log=TRUE)
      })))
      
      crit <- norm(beta_old-beta_new, "2")
      
      loglik.list[l1] <- loglik
      crit.list[l1] <- crit
      
      l1 <- l1 + 1
    }
    
    # print( matrix(beta,p,q) )
    # print( l1 )
    
    
    # plot(loglik.list, type="l")
    
    
    
    
    
    beta.list <- lapply(beta.list, function(bb) matrix(bb, p+1, q, byrow=FALSE) )
    beta_new <- matrix(beta_new, p+1, q, byrow=FALSE)
    
    if(standardize){
      beta.list <- lapply(beta.list, function(b) magic::adiag(1, sdx.inv) %*% b)
      beta_new <- magic::adiag(1, sdx.inv) %*% beta_new
    }
    
    
    
    result <- list(X=X, Y=Y,
                   mu=beta_new[1,], 
                   beta=beta_new,
                   beta.list=beta.list,
                   offset=Offset,
                   lambda1=lambda1,
                   lambda2=lambda2,
                   # beta=matrix(beta,p,q,byrow=TRUE),
                   # beta.list=beta.list %>% lapply(function(x) matrix(x,p,q,byrow=TRUE)),
                   # beta=matrix(beta,p+1,q,byrow=TRUE),
                   # beta.list=beta.list %>% lapply(function(x) matrix(x,p+1,q,byrow=TRUE)),
                   loglik.list=loglik.list[-1], crit.list=crit.list[-1])  
    
  }
  
  result$beta2 <- result$beta
  result$beta2[-1,] <- t( apply(result$beta[-1,,drop=F], 1, function(x) x/norm(result$beta[1,], "2")) )
  result$beta2[1,] <- result$beta[1,]/norm(result$beta[1,], "2")
  
  
  return(result)
  
}








#' @export Cq
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



#' @export Bq
Bq <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  # besselI(NORM, nu=q/2, expon.scaled=FALSE) == besselI(NORM, nu=q/2, expon.scaled=TRUE)/exp(-NORM)
  # besselI(NORM, nu=q/2-1, expon.scaled=FALSE) == besselI(NORM, nu=q/2-1, expon.scaled=TRUE)/exp(-NORM)
  
  besselI(NORM, nu=q/2, expon.scaled=TRUE) / besselI(NORM, nu=q/2-1, expon.scaled=TRUE)
}


#' @export Hq
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



#' @export subgrad
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
  tcrossprod(s) * Hq(theta) + Bq(theta) / NORM * (diag(1,q) - tcrossprod(theta/NORM) )
  # }
  
}


