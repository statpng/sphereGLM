#' @importFrom manifold frechetMean createM
#' @importFrom magic adiag
#' @importFrom dplyr `%>%`
#' 
#' @export deprecated.sphereGLM
deprecated.sphereGLM <- function(X, Y, MU = NULL, Offset = NULL, beta0 = NULL, 
                         maxit = 100, eps = 1e-6, standardize = TRUE, lambda = 1e-4, 
                         use.nlm = TRUE) {
  if(FALSE){
    X <- simdata1$X
    Y <- simdata1$Y
    MU = NULL; Offset = NULL; beta0 = NULL;
    maxit = 100; eps = 1e-6; standardize = TRUE; lambda = 1e-4;
    use.nlm = TRUE
  }
  
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  X1 <- cbind(1, X)
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  if (standardize) {
    X_scaled <- scale(X, center = TRUE, scale = TRUE)
    sdx_inv <- 1 / attr(X_scaled, "scaled:scale")
    X <- X_scaled
  } else {
    X <- scale(X, center = TRUE, scale = FALSE)
  }
  
  
  
  if (use.nlm) {
    
    regvmf <- function(be, y, x) {
      theta <- x %*% matrix(be, ncol = q)
      ki <- sqrt(rowSums(theta^2))
      -sum(apply(theta, 1, Cq, log = TRUE)) - sum(theta * y)
    }
    
    ini <- solve(crossprod(X1), crossprod(X1, Y))
    
    suppressWarnings({
      val1 <- nlm(regvmf, ini, y = Y, x = X1, iterlim = maxit)
      val2 <- nlm(regvmf, val1$estimate, y = Y, x = X1, iterlim = maxit)
      da <- optim(val2$estimate, regvmf, y = Y, x = X1, 
                  control = list(maxit = maxit), 
                  hessian = TRUE)
    })
    
    result <- list(
      X = X,
      Y = Y,
      mu = matrix(da$par, ncol = q)[1,],
      beta = matrix(da$par, ncol = q),
      loglik = -da$value
    )
    
    if (standardize) {
      result$beta <- adiag(1, diag(sdx_inv)) %*% result$beta
    }
    
    result$beta2 <- result$beta
    
  } else {
    
    if (is.null(Offset)) Offset <- matrix(0, n, q)
    
    if (is.null(MU)) {
      MU <- manifold::frechetMean(manifold::createM("Sphere"), t(Y))
      muhat <- colSums(Y) %>% {./norm(.,"2")} # mle of mean direction in vMF
      rbar <- colSums(Y) %>% {norm(.,"2")/nrow(Y)}
      kappa0 <- (rbar * q - rbar^3) / (1 - rbar^2)
      MU <- muhat * kappa0
    }
    
    B0 <- coef(lm(Y ~ X - 1, offset = rep(1, n) %o% MU + Offset))
    beta <- c(MU, as.vector(B0))
    
    
    Xt_list <- lapply(seq_len(n), function(i) kronecker(diag(q), X1[i,]))
    Xt <- do.call(cbind, Xt_list)
    yt <- as.vector(t(Y))
    
    beta_old <- beta + 1
    beta_new <- beta
    beta_list <- list()
    loglik_list <- numeric(maxit)
    crit_list <- numeric(maxit)
    
    for (l1 in seq_len(maxit)) {
      beta_old <- beta_new
      beta_list[[l1]] <- beta_old
      
      theta_list <- lapply(seq_len(n), function(i) Offset[i,] + crossprod(Xt_list[[i]], beta_old))
      eta <- do.call(rbind, lapply(theta_list, b1.vMF))
      W <- Matrix::bdiag(lapply(theta_list, b2.vMF))
      
      Z <- crossprod(Xt, beta_old) + solve(W + diag(lambda, nrow(W))) %*% (yt - as.vector(eta))
      
      XWX <- Xt %*% W %*% t(Xt)
      beta_new <- solve(XWX + Matrix::bdiag(0, diag(lambda, nrow(XWX)-1))) %*% (Xt %*% W %*% Z)
      
      loglik <- sum(mapply(function(theta, y) sum(theta * y) + Cq(theta, log = TRUE), theta_list, apply(Y, 1, c, simplify=FALSE)))
      
      crit <- sqrt(sum((beta_old - beta_new)^2))
      
      loglik_list[l1] <- loglik
      crit_list[l1] <- crit
      
      if (crit <= eps) break
    }
    
    beta_matrix <- matrix(beta_new, p + 1, q, byrow = FALSE)
    if (standardize) {
      beta_matrix <- adiag(1, diag(sdx_inv)) %*% beta_matrix
    }
    
    result <- list(
      X = X,
      Y = Y,
      mu = beta_matrix[1,],
      beta = beta_matrix,
      beta2 = beta_matrix,
      beta_list = lapply(beta_list, function(b) matrix(b, p + 1, q, byrow = FALSE)),
      offset = Offset,
      loglik_list = loglik_list[seq_len(l1)],
      crit_list = crit_list[seq_len(l1)]
    )
  }
  
  result$beta2[-1,] <- t(apply(result$beta2[-1,, drop = FALSE], 1, function(x) x / sqrt(sum(result$beta2[1,]^2))))
  result$beta2[1,] <- result$beta2[1,] / sqrt(sum(result$beta2[1,]^2))
  
  result$params <- list(MU = MU, Offset = Offset, beta0 = beta0, maxit = maxit, eps = eps, standardize = standardize, use.nlm = use.nlm)
  
  return(result)
}