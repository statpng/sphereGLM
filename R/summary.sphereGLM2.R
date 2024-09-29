#' @export summary.sphereGLM2
summary.sphereGLM2 <- function(fit, scale = FALSE, orthogonal = NULL) {
  if (inherits(fit, "try-error")) {
    return(list(
      beta = NA, norm = NA, Fisher = NA, orthogonal = NA,
      wald = NA, df = NA,
      statistic = list(wald = NA, LRT = NA, score = NA),
      pvalue = list(wald = NA, LRT = NA, score = NA)
    ))
  }
  
  if (is.null(orthogonal)){
    orthogonal <- fit$params$orthogonal 
  }
  
  X <- fit$X
  Y <- fit$Y
  Offset <- fit$offset
  
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)
  
  df <- ifelse(orthogonal, q - 1, q)
  gamma <- ifelse(orthogonal, fit$params$gamma, 0)
  
  if( scale ){
    beta_new <- fit$beta2
  } else {
    beta_new <- fit$beta
  }
  
  # Fisher Information Matrix
  fit.Fisher <- FisherMatrix(X, beta_new, orthogonal, gamma)
  F_nj.list <- fit.Fisher$Fnj
  
  
  
  # Wald test ----
  Wald.global <- as.numeric(crossprod(as.vector(t(beta_new[-1,])), fit.Fisher$Fn %*% as.vector(t(beta_new[-1,]))))
  pv.wald.global <- pchisq(Wald.global, df = p * df, lower.tail = FALSE)
  
  Wald.j.list <- sapply(seq_len(p), function(j) {
    crossprod(beta_new[j + 1,], F_nj.list[[j]] %*% beta_new[j + 1,])
  })
  pv.wald.ind <- pchisq(Wald.j.list, df, lower.tail = FALSE)
  
  pv.wald <- c(global = pv.wald.global, ind = pv.wald.ind)
  stat.wald <- c(global = Wald.global, ind = Wald.j.list)
  
  
  # LR test ----
  beta0 <- rbind(fit$mu, matrix(0, p, q))
  
  TestStat.LRT.global <- -2 * (loglik.vmf(X, Y, beta0, gamma) - loglik.vmf(X, Y, beta_new, gamma))
  
  if (p == 1) {
    TestStat.LRT.ind <- TestStat.LRT.global
  } else {
    betahat0 <- lapply(seq_len(p), function(j) {
      
      X2 <- X
      X2[,j] <- 0
      fit0 <- sphereGLM(X = X2, Y = Y, orthogonal = orthogonal, gamma = gamma)
      
      fit0 <- sphereGLM(X = X[, -j, drop = FALSE], Y = Y, orthogonal = orthogonal, gamma = gamma)
      beta0 <- matrix(0, p + 1, q)
      beta0[1,] <- fit0$beta[1,]
      beta0[-c(1, j + 1),] <- fit0$beta[-1,]
      beta0
    })
    
    TestStat.LRT.ind <- sapply(betahat0, function(b0) {
      -2 * (loglik.vmf(X, Y, b0, gamma) - loglik.vmf(X, Y, beta_new, gamma))
    })
  }
  
  pv.LRT.global <- pchisq(TestStat.LRT.global, p * df, lower.tail = FALSE)
  pv.LRT.ind <- pchisq(TestStat.LRT.ind, df, lower.tail = FALSE)
  
  pv.LRT <- c(global = pv.LRT.global, ind = pv.LRT.ind)
  stat.LRT <- c(global = TestStat.LRT.global, ind = TestStat.LRT.ind)
  
  
  # Rao's Score Test ----
  beta0 <- rbind(fit$mu, matrix(1e-12, p, q))
  fit.Fisher <- FisherMatrix(X, beta0, orthogonal, gamma)
  
  b1.list <- lapply(seq_len(n), function(i) {
    b1.vMF(Offset[i,] + t(beta0) %*% c(1, X[i,]))
  })
  
  ScoreVector <- Reduce("+", lapply(seq_len(n), function(i) {
    kronecker(X[i,], diag(1, q)) %*% (Y[i,] - b1.list[[i]])
  })) + 2 * gamma * Matrix::bdiag(replicate(p, tcrossprod(beta0[1,]), simplify = FALSE)) %*% as.vector(t(beta0[-1,]))
  
  TestStat.Score.global <- crossprod(ScoreVector, solve(fit.Fisher$Fn) %*% ScoreVector)
  
  
  if (p == 1) {
    TestStat.Score.ind <- TestStat.Score.global
  } else {
    TestStat.Score.ind <- sapply(seq_len(p), function(j) {
      b0 <- betahat0[[j]]
      Xnew <- X
      pnew <- ncol(Xnew)
      
      fit.Fisher <- FisherMatrix(Xnew, b0, orthogonal, gamma)
      
      b1.list <- lapply(seq_len(n), function(i) {
        b1.vMF(Offset[i,] + t(b0) %*% c(1, Xnew[i,]))
      })
      
      ScoreVector <- Reduce("+", lapply(seq_len(n), function(i) {
        kronecker(Xnew[i,], diag(1, q)) %*% (Y[i,] - b1.list[[i]])
      })) + 2 * gamma * Matrix::bdiag(replicate(pnew, tcrossprod(b0[1,]), simplify = FALSE)) %*% as.vector(t(b0[-1,]))
      
      crossprod(ScoreVector, solve(fit.Fisher$Fn, ScoreVector))
    })
  }
  
  stat.score <- c(global = TestStat.Score.global, ind = TestStat.Score.ind)
  pv.score <- pchisq(stat.score, df, lower.tail = FALSE)
  
  list(
    beta = beta_new,
    norm = apply(beta_new, 1, norm, "2"),
    Fisher = F_nj.list,
    orthogonal = fit$params$orthogonal,
    wald = Wald.j.list,
    df = df,
    statistic = list(wald = stat.wald, LRT = stat.LRT, score = stat.score),
    pvalue = list(wald = pv.wald, LRT = pv.LRT, score = pv.score)
  )
}