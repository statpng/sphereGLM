#' @method summary sphereGLM
#' @export 
summary.sphereGLM <- function(fit, scale=FALSE, orthogonal=NULL, type="all"){
  
  if(FALSE){
    devtools::load_all()
    
    simdata.list <- lapply(c(50, 100, 200, 500, 1000), function(n){
      simdata <- sim.sphereGLM(n=n, p=2, q=3, mu=c(0,0,10), s=0, s0=0, type="vMF", seed.UDV=1, seed.E=2, qr.U=FALSE, qr.V=TRUE)
    })
    
    fit.list <- lapply(1:length(simdata.list), function(i){
      simdata <- simdata.list[[i]]
      fit <- sphereGLM(X=simdata$X, Y=simdata$Y, orthogonal=simdata$params$qr.V, gamma=100)
      fit
    })
    
    summary.list <- fit.list %>% lapply(function(x) summary(x))
    
    lapply(summary.list, function(x) do.call("rbind", x$statistic))
    
    
    cov1.list <- summary.list %>% lapply(function(x) x$cov[[1]] %>% {eigen(.)$values[1]})
    cov2.list <- summary.list %>% lapply(function(x) x$cov[[2]] %>% {eigen(.)$values[1]})
    stat.list <- summary.list %>% lapply(function(x) x$stat)
    Fisher.list <- summary.list %>% lapply(function(x) x$Fisher)
    
    
    
    
    
    simdata.list2 <- lapply(c(100, 500, 1000, 2000, 5000, 10000), function(n){
      simdata <- sim.sphereGLM(n=n, p=2, q=3, mu=c(0,0,20), s=5, s0=0, type="vMF", seed.UDV=1, qr.U=FALSE, qr.V=FALSE)
    })
    
    eigvalue.list <- lapply(simdata.list2, function(simdata){
      orthogonal <- TRUE
      FisherMatrix(X=simdata$X, beta_new=rbind(simdata$mu,simdata$B), orthogonal=orthogonal, gamma=ifelse(orthogonal, 9999, 0))$Fn %>% eigen %>% .$values %>% min
    })
    
    eigvalue.list
    
    
    
    
    #
    #
    #
    
    summary.sphereGLM(fit)$pvalue
    
    summary.sphereGLM(fit)$stat
    
    
    #
    #
  }
  
  
  
  # Asymptotic Equivalence between three types of test ----
  if(FALSE){
    
    
    # TRUE beta_j is not orthogonal to mu
    set.seed(2)
    simdata1 <- sim.sphereGLM(n=500, p=2, q=10, mu=c(rep(0,9),10), s=c(0,1), s0=0, type="vMF", seed.UDV=1, qr.U=FALSE, qr.V=FALSE)
    
    fit1 <- sphereGLM(X=simdata1$X, Y=simdata1$Y, orthogonal=FALSE)
    fit.summary1 <- summary( fit1 )
    fit.summary1$stat
    
    fit2 <- sphereGLM(X=simdata1$X, Y=simdata1$Y, orthogonal=TRUE, gamma=99)
    fit.summary2 <- summary( fit2 )
    fit.summary2$stat
    
    
    
    # TRUE beta_j is orthogonal to mu
    set.seed(2)
    simdata2 <- sim.sphereGLM(n=500, p=2, q=3, mu=c(0,0,20), s=c(0,1), s0=0, type="vMF", seed.UDV=1, qr.U=FALSE, qr.V=TRUE)
    
    fit3 <- sphereGLM(X=simdata2$X, Y=simdata2$Y, orthogonal=FALSE)
    fit.summary3 <- summary( fit3 )
    fit.summary3$stat
    
    fit4 <- sphereGLM(X=simdata2$X, Y=simdata2$Y, orthogonal=TRUE)
    fit.summary4 <- summary( fit4 )
    fit.summary4$stat
    
    
    # bad result
    set.seed(2)
    simdata2 <- sim.sphereGLM(n=500, p=2, q=10, mu=c(rep(0,9),10), s=c(0,1), s0=0, type="vMF", seed.UDV=1, qr.U=FALSE, qr.V=TRUE) 
    
    
    fit4 <- sphereGLM(X=simdata2$X, Y=simdata2$Y, orthogonal=TRUE)
    fit.summary4 <- summary( fit4 )
    fit.summary4$stat
    
    
  }
  
  
  
  
  if( type == "all" ){
    
    type.wald <- type.LRT <- type.score <- TRUE
    
  } else if ( type %in%  c("wald", "LRT", "score") ){
    
    type.wald <- ifelse(type == "wald", TRUE, FALSE)
    type.LRT <- ifelse(type == "LRT", TRUE, FALSE)
    type.score <- ifelse(type == "score", TRUE, FALSE)
    
  }
  
  if( class(fit) == "try-error" ){
    
    return( list(beta=NA,
                 norm=NA, 
                 Fisher=NA,
                 orthogonal=NA,
                 wald=NA,
                 df=NA,
                 statistic=list(wald=NA,
                                LRT=NA,
                                score=NA),
                 pvalue=list(wald=NA,
                             LRT=NA,
                             score=NA
                 )) )
    
    
  }
  
  
  if(is.null(orthogonal)){
    orthogonal <- fit$params$orthogonal
  }
  
  
  
  
  {
    
    X <- fit$X
    Y <- fit$Y
    Offset <- fit$offset
    
    p <- ncol(X)
    q <- ncol(Y)
    n <- nrow(Y)
    
    df <- ifelse(orthogonal, q-1, q)
    gamma <- ifelse(orthogonal, fit$params$gamma, 0)
    
    beta_new <- fit$beta
    # if( scale ){
    #   beta_new <- fit$beta2
    # }
    
  }
  
  
  
  
  
  # gamma
  if(orthogonal){
    bbeta <- as.vector(t(beta_new[-1,]))
    bmu <- beta_new[1,]
    ScoreVector <- score2(X, Y, beta_new, gamma=0, Offset=Offset)
    # ==
    # Xt.list <- apply(cbind(1,X), 1, function(Xi) kronecker(Xi, diag(1,q)), simplify=FALSE)
    # b1.list <- lapply(1:n, function(i) b1.vMF( Offset[i,] + t(Xt.list[[i]]) %*% as.vector(t(beta_new)) ))
    # t(Xt) %*% (as.vector(t(Y)) - do.call("rbind", b1.list))
    
    gamma <- (kronecker(diag(p), tcrossprod(bmu)) %*% bbeta) %>% {
      solve(t(.) %*% (.) + diag(1e-16, ncol(.), ncol(.)) ) %*% t(.)
    } %*% ScoreVector %>% {c(./2)}
    
    # gamma <- 0
    #
    
  }
  
  
  
  betahat0 <- lapply(1:p, function(j){
    # fit0 <- vmf.reg(y=Y, x=X[,j], tol=1e-4, orthogonal=orthogonal, gamma=gamma)
    
    penalty.factor <- rep(0,p)
    penalty.factor[j] <- 1
    fit0 <- sphereGLM(X=X, Y=Y, orthogonal=orthogonal, gamma=gamma, penalty.factor=penalty.factor)
    fit0$beta
    # fit0 <- sphereGLM(X=X[,-j], Y=Y, orthogonal=orthogonal, gamma=gamma)
    # beta0 <- matrix(0, p+1, q)
    # beta0[1,] <- fit0$beta[1,]
    # beta0[-c(1,j+1),] <- fit0$beta[-1,]
    # beta0
  })
  
  
  
  
  
  # Wald test ----
  
  fit.Fisher <- FisherMatrix(X=X, beta_new=beta_new, orthogonal=orthogonal, gamma=gamma)
  
  cov <- lapply(1:p, function(j) solve(fit.Fisher$Fnj[[j]]))
  
  
  if(type.wald){
    
    Fn <- fit.Fisher$Fn
    F_nj.list <- fit.Fisher$Fnj
    
    Wald.global <- {
      bhat <- as.vector(t(beta_new[-1,]))
      t(bhat) %*% fit.Fisher$Fn[-(1:q),-(1:q)] %*% bhat
    } %>% as.numeric
    
    pv.wald.global <- 1-pchisq(Wald.global, df = p*df)
    
    Wald.j.list <- sapply(2:(p+1), function(j){
      t(beta_new[j,]) %*% F_nj.list[[j]] %*% beta_new[j,]
    })
    
    pv.wald.ind <- sapply(Wald.j.list, function(x){
      1 - pchisq(x, df)
    })
    
    pv.wald <- c(global=pv.wald.global, ind=pv.wald.ind)
    
    stat.wald <- c( global=Wald.global, ind=Wald.j.list )
  }
  
  
  
  
  
  # LR test ----
  if(type.LRT){
    beta0 <- rbind(fit$mu, matrix(0, p, q))
    
    TestStat.LRT.global <- -2 * ( loglik.vmf( X, Y, beta0, gamma=gamma) - loglik.vmf( X, Y, beta_new, gamma=gamma) )
    
    if( p == 1 ){
      
      TestStat.LRT.ind <- TestStat.LRT.global
      
    } else {
      
      TestStat.LRT.ind <- sapply(betahat0, function(b0){
        -2 * ( loglik.vmf( X, Y, b0, gamma=gamma) - loglik.vmf( X, Y, beta_new, gamma=gamma) )
      })
      
    }
    
    
    pv.LRT.global <- 1 - pchisq(TestStat.LRT.global, p*df)
    pv.LRT.ind <- sapply(1:p, function(j) 1 - pchisq(TestStat.LRT.ind[j], df) )
    
    
    pv.LRT <- c(global=pv.LRT.global, ind=pv.LRT.ind)
    stat.LRT <- c( global=TestStat.LRT.global, ind=TestStat.LRT.ind )
    
  }
  
  
  
  # Rao's Score Test ----
  if(type.score){
    
    beta0 <- rbind(fit$mu, matrix(1e-12, p, q))
    fit.Fisher <- FisherMatrix(X=X, beta_new=beta0, orthogonal=orthogonal, gamma=gamma)
    ScoreVector <- score(X, Y, beta0, gamma=gamma, Offset=Offset)
    
    TestStat.Score.global <- c( t(ScoreVector) %*% solve(fit.Fisher$Fn) %*% ScoreVector )
    TestStat.Score.global <- c( t(ScoreVector[-(1:q)]) %*% solve( fit.Fisher$Fn[-(1:q),-(1:q)]) %*% ScoreVector[-(1:q)] )
    TestStat.Score.global
    
    
    if(p == 1){
      
      TestStat.Score.ind <- TestStat.Score.global
      
    } else {
      
      TestStat.Score.ind <- sapply(1:p, function(j){
        beta0 <- betahat0[[j]]
        fit.Fisher <- FisherMatrix(X=X, beta_new=beta0, orthogonal=orthogonal, gamma=gamma)
        ScoreVector <- score(X, Y, beta0, gamma=gamma, Offset=Offset)
        
        # ScoreStat.ind <- c( t(ScoreVector) %*% solve(fit.Fisher$Fn) %*% ScoreVector )
        
        ScoreStat.ind <- c( t(ScoreVector[1:q+q*j]) %*% solve( fit.Fisher$Fnj[[j+1]]) %*% ScoreVector[1:q+q*j] )
        
        ScoreStat.ind
        
      })
      
      
    }
    
    stat.score <- c( global=TestStat.Score.global, ind=TestStat.Score.ind )
    
    pv.score.global <- 1 - pchisq(TestStat.Score.global, p*df)
    pv.score.ind <- sapply(1:p, function(j) 1 - pchisq(TestStat.Score.ind[j], df) )
    
    pv.score <- c( global=pv.score.global, ind=pv.score.ind )
    
  }
  
  
  
  
  list(beta=beta_new,
       norm=apply(beta_new, 1, norm, "2"), 
       Fisher=list(Fn=Fn, Fnj=F_nj.list),
       orthogonal=fit$params$orthogonal,
       wald=Wald.j.list,
       df=df,
       cov=cov,
       statistic=list(wald=stat.wald,
                      LRT=stat.LRT,
                      score=stat.score),
       pvalue=list(wald=pv.wald,
                   LRT=pv.LRT,
                   score=pv.score
       ))
  
  
  
}




#' #' @export summary.sphereGLM2
#' summary.sphereGLM2 <- function(fit, scale=FALSE, orthogonal=NULL){
#'   
#'   if( class(fit) == "try-error" ){
#'     return( list(beta=NA,
#'                  norm=NA, 
#'                  fisher=NA,
#'                  wald=NA,
#'                  pvalue=NA) )
#'   }
#'   
#'   if(is.null(orthogonal)){
#'     orthogonal <- fit$params$orthogonal
#'   }
#'   
#'   X <- fit$X
#'   Y <- fit$Y
#'   
#'   p <- ncol(X)
#'   q <- ncol(Y)
#'   n <- nrow(Y)
#'   
#'   if( scale ){
#'     beta_new <- fit$beta2
#'   } else {
#'     beta_new <- fit$beta
#'   }
#'   
#'   
#'   F_nj.list <- FisherInfoMatrix(X, beta_new)
#'   
#'   Wald.j.list <- sapply(1:(p+1), function(j){
#'     t(beta_new[j,]) %*% F_nj.list[[j]] %*% beta_new[j,]
#'   })
#'   
#'   
#'   if(orthogonal){
#'     df <- q-1
#'   } else {
#'     df <- q
#'   }
#'   
#'   
#'   pv.list <- sapply(Wald.j.list, function(x){
#'     1 - pchisq(x, df)
#'   })
#'   
#'   
#'   
#'   
#'   list(beta=beta_new,
#'        norm=apply(beta_new, 1, norm, "2"), 
#'        fisher=F_nj.list,
#'        orthogonal=fit$params$orthogonal,
#'        wald=Wald.j.list,
#'        df=df,
#'        pvalue=pv.list)
#'   
#'   
#'   
#' }




# #' @export FisherInfoMatrix
# FisherInfoMatrix <- function(X, beta_new, OffSet=NULL){
#   
#   n <- nrow(X)
#   p <- ncol(X)
#   q <- ncol(beta_new)
#   
#   if(is.null(OffSet)){
#     OffSet <- matrix(0, n, q)
#   }
#   
#   bhat <- as.vector(t(beta_new))
#   
#   X1 <- cbind(1, X)
#   Xt.list <- apply(X1, 1, function(Xi) kronecker(Xi, diag(1,q)), simplify=FALSE) # [(p+1)q x q] x n
#   Xt <- do.call("cbind", Xt.list) # (p+1)q x qn
#   b2.list <- lapply(1:n, function(i) b2.vMF( OffSet[i,] + crossprod(Xt.list[[i]], bhat) ))
#   
#   
#   Fnj <- NULL
#   for( jj in 1:(p+1) ){
#     out <- 0
#     for( i in 1:nrow(X1) ){
#       
#       idx <- 1:q + (jj-1)*q
#       Xi.jj <- kronecker(X1[i,], diag(1,q,q))[idx,,drop=F]
#       out <- out + Xi.jj %*% b2.list[[i]] %*% t(Xi.jj)
#       
#     }
#     
#     Fnj[[jj]] <- out
#   }
#   
#   Fnj
# }




#' @export loglik.vmf
loglik.vmf <- function(X, Y, beta, gamma=0){
  
  q <- ncol(Y)
  X1 <- cbind(1, X)
  
  Cq0 <- function(NORM, q){
    (  (q/2-1) * log(NORM) - ( (q/2) * log(2*pi) + log( besselI(NORM, q/2-1, expon.scaled=TRUE) ) + NORM )  )
  }
  
  B <- matrix(beta, ncol = q)
  Theta <- X1 %*% B
  NORM <- sqrt(Rfast::rowsums(Theta^2))
  
  sum( Cq0(NORM, q) ) + sum(Theta * Y) + gamma * sum( apply(beta[-1,,drop=FALSE], 1, function(x) crossprod(x, beta[1,])^2 ) )
}




#' @export score
score <- function(X, Y, beta_new, gamma, Offset){
  n=nrow(X); p=ncol(X); q=ncol(Y)
  bmu=beta_new[1,]
  bB=beta_new[-1,,drop=F]
  bbeta=as.vector(t(beta_new[-1,]))
  
  b1.list <- lapply(1:n, function(i){
    b1.vMF( Offset[i,] + bmu + t(bB) %*% X[i,] )
  } )
  
  
  ScoreVector <- matsum(1:n, function(i){
    # \partial\bmu
    s1 <- -2*gamma * matsum(1:p, function(j) tcrossprod(bB[j,])) %*% bmu
    # \partial\bbeta^*
    s2 <- -2*gamma * kronecker(diag(p), tcrossprod(bmu)) %*% bbeta
    
    kronecker(c(1,X[i,]), diag(1,q,q)) %*% ( Y[i,] - b1.list[[i]] ) + c(s1, s2)
    
  })
  
  ScoreVector
}



#' @export score2
score2 <- function(X, Y, beta_new, gamma, Offset){
  n=nrow(X); p=ncol(X); q=ncol(Y)
  bmu=beta_new[1,]
  bB=beta_new[-1,,drop=F]
  bbeta=as.vector(t(beta_new[-1,]))
  
  b1.list <- lapply(1:n, function(i){
    b1.vMF( Offset[i,] + bmu + t(bB) %*% X[i,] )
  } )
  
  
  ScoreVector <- matsum(1:n, function(i){
    # \partial\bmu
    s1 <- -2*gamma * matsum(1:p, function(j) tcrossprod(bB[j,])) %*% bmu
    # \partial\bbeta^*
    s2 <- -2*gamma * kronecker(diag(p), tcrossprod(bmu)) %*% bbeta
    
    kronecker(c(X[i,]), diag(1,q,q)) %*% ( Y[i,] - b1.list[[i]] ) + s2 # c(s1, s2)
    
  })
  
  ScoreVector
}






#' @export matsum
matsum <- function(seq, FUN){
  Reduce("+", lapply(seq, FUN))
}



eigvalue.b2 <- function(theta){
  # theta <- (cbind(1,fit$X) %*% fit$beta)[i,]
  
  NORM <- norm(theta, "2")
  B <- Bq_cpp(theta)
  q <- length(theta)
  
  B / NORM + 1 - B^2 - q/NORM*B
  
  # b2.list[[i]] %>% {eigen(.)$values}
}



#' @export FisherMatrix
FisherMatrix <- function(X, beta_new, orthogonal=NULL, gamma=0, OffSet=NULL){
  
  if(FALSE){
    X=Xnew; beta_new=b0; orthogonal=orthogonal; gamma=gamma
  }
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(beta_new)
  
  if(is.null(OffSet)){
    OffSet <- matrix(0, n, q)
  }
  
  bhat <- as.vector( t(beta_new) )
  bbeta <- as.vector( t(beta_new[-1,,drop=F]) )
  bB <- beta_new[-1,,drop=F]
  bmu <- beta_new[1,]
  
  # X1 <- cbind(1, X)
  Xt.list <- apply(X, 1, function(Xi) kronecker(Xi, diag(1,q)), simplify=FALSE) # [(p+1)q x q] x n
  Xt <- t( do.call("cbind", Xt.list) ) # (p+1)q x qn
  b2.list <- lapply(1:n, function(i) b2.vMF( OffSet[i,] + bmu + crossprod(Xt.list[[i]], bbeta) ))
  
  
  F11 <- -matsum(1:n, function(i) b2.list[[i]] ) - gamma*matsum(1:p, function(j) tcrossprod(bB[j,]))
  
  F22 <- -matsum(1:n, function(i) Xt.list[[i]] %*% b2.list[[i]] %*% t(Xt.list[[i]])) - gamma*kronecker(diag(p), tcrossprod(bmu))
  
  F21 <- -matsum(1:n, function(i) Xt.list[[i]] %*% b2.list[[i]] ) - gamma * do.call("rbind", lapply(1:p, function(j) c(crossprod(bmu, bB[j,]))*diag(q) + tcrossprod(bmu, bB[j,])))
  
  Fn <- matrix(NA, (p+1)*q, (p+1)*q)
  
  Fn[1:q,1:q] <- -F11
  Fn[(q+1):((p+1)*q),(q+1):((p+1)*q)] <- -F22
  Fn[1:q,(q+1):((p+1)*q)] <- -F21
  Fn[(q+1):((p+1)*q),1:q] <- -t(F21)
  
  idx.list <- lapply(1:(p+1), function(j) (j-1)*q + 1:q)
  Fnj <- lapply(idx.list, function(idx) Fn[idx,idx])
  
  list(b2.list=b2.list, Fn=Fn, Fnj=Fnj)
}







# FisherMatrix2 <- function(X, beta_new, OffSet=NULL){
#   
#   n <- nrow(X)
#   p <- ncol(X)
#   q <- ncol(beta_new)
#   
#   if(is.null(OffSet)){
#     OffSet <- matrix(0, n, q)
#   }
#   
#   bhat <- as.vector(beta_new)
#   
#   X1 <- cbind(1, X)
#   Xt.list <- apply(X1, 1, function(Xi) kronecker(diag(1,q), Xi), simplify=FALSE) # [(p+1)q x q] x n
#   Xt <- do.call("cbind", Xt.list) # (p+1)q x qn
#   b2.list <- lapply(1:n, function(i) b2.vMF( OffSet[i,] + crossprod(Xt.list[[i]], bhat) ))
#   
#   
#   Fnj <- NULL
#   for( j in 1:(p+1) ){
#     out <- 0
#     for( i in 1:nrow(X1) ){
#       # Xj.list <- lapply(idx.list, function(jj){
#       #   Xt.list[[i]][jj,]
#       # })
#       Xj <- do.call("rbind", lapply(1:q, function(jj){
#         idx <- j + (jj-1)*(p+1)
#         Xt.list[[i]][idx,]
#       }))
#       out <- out + Xj %*% b2.list[[i]] %*% t(Xj)
#     }
#     
#     Fnj[[j]] <- out
#   }
#   
#   list(b2.list=b2.list, Fn=Fn, Fnj=Fnj)
# }























# Rank of FIM -------------------------------------------------------------
if(FALSE){
  
  set.seed(1)
  simdata <- sim.sphereGLM(n=200, p=3, q=3, mu=c(rep(0,2),10), s=5, s0=0, type="vMF", seed.UDV=1, qr.U=FALSE, qr.V=FALSE)
  
  fit.list <- lapply(c(FALSE, TRUE), function(ortho) sphereGLM(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=1e-10, gamma=1e+4, orthogonal=ortho))
  
  fit.list[[1]] %>% {summary(.)$Fisher[[2]]} %>% { eigen(.)$values } %>% {max(.)/min(.)}
  fit.list[[2]] %>% {summary(.)$Fisher[[2]]} %>% { eigen(.)$values } %>% {max(.)/min(.)}
  
  
  fit.list[[1]]$beta[-1,] %>% tcrossprod() %>% {svd(.)$d}
  fit.list[[2]]$beta[-1,] %>% tcrossprod() %>% {svd(.)$d}
  
  #
  #
  
  summary(fit)$Fisher
  
  FisheroMatrix(fit$X, fit$beta)
  
  #
  #
  #
  
  
  plot(fit)
  summary(fit)$Fisher %>% lapply(function(.) eigen(.)$values)
  
  fit$beta 
  
  
}


# check the minimum eigenvalue of FIM ----
if(FALSE){
  
  set.seed(1)
  simdata <- sim.sphereGLM(n=100, p=1, q=3, mu=c(0,0,100), s=10, s0=0, type="vMF", seed.UDV=1, qr.U=FALSE, qr.V=F)
  
  fit <- sphereGLM(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=1e-10, orthogonal=T)
  
  
  bhat <- as.vector(t(fit$beta))
  X <- fit$X
  X1 <- cbind(1, X)
  Xt.list <- apply(X1, 1, function(Xi) kronecker(Xi, diag(1,q)), simplify=FALSE) # [(p+1)q x q] x n
  Xt <- t(do.call("cbind", Xt.list))
  b2.list <- lapply(1:n, function(i) b2.vMF( OffSet[i,] + crossprod(Xt.list[[i]], bhat) ))
  
  summary(fit)$Fisher
  (t(Xt) %*% Matrix::bdiag(b2.list) %*% Xt)
  
  
  summary(fit)$Fisher[[1]] %>% {eigen(.)$values} %>% {
    print(.)
    max(.)/min(.)
  }
  
  
  b2.list[[i]] %>% {eigen(.)$values}
  
  (cbind(1,fit$X) %*% fit$beta)[i,] %>% {
    NORM <- norm(.,"2")
    B <- Bq_cpp(.)
    q <- length(.)
    
    B / NORM + 1 - B^2 - q/NORM*B
  }
  
  
  #
  
  
  
  
  
}
