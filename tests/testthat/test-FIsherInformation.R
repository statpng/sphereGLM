# detach("package:sphereGLM",unload=TRUE)
library(sphereGLM)
library(dplyr)
library(png.utils)
# detach("package:png.utils", unload=TRUE)


{
  n.seq <- c(50,200)
  p.seq <- c(2)
  q.seq <- c(3)
  m.seq <- c(10, 20)
  s.seq <- c(0, 2, 4)
  ortho.seq <- c(TRUE, FALSE)
  seed.seq <- c(1:50)
  nrep <- 500
  
  params <- list(n.seq=n.seq,
                 p.seq=p.seq,
                 q.seq=q.seq,
                 m.seq=m.seq,
                 s.seq=s.seq,
                 ortho.seq=ortho.seq,
                 seed.seq=seed.seq)
  
  grid <- expand.grid(params)
}





print(paste0(ii, " / ", max(AA)))

# Simulation Parameters
{
  n=grid[ii,]$n.seq
  p=grid[ii,]$p.seq
  q=grid[ii,]$q.seq
  s=grid[ii,]$s.seq
  ortho=grid[ii,]$ortho.seq
  mu=grid[ii,]$m.seq
  seed=grid[ii,]$seed.seq
  
  mu=c(0,0,mu)
  s0=0
  type="vMF"
  snr=NULL
}






simdata.list <- fit.list <- NULL
for( jj in 1:nrep ){
  if( jj %% 50 == 0 ) print(jj)
  simdata <- sim.sphereGLM(n=n, p=p, q=q, mu=mu, s=s, s0=s0, type="vMF",
                           seed.UDV=seed, seed.E=100*seed+jj, qr.U=FALSE, qr.V=ortho)
  fit.nonortho <- try(with(simdata, sphereGLM(X=X, Y=Y, orthogonal=FALSE, lambda=1e-12)))
  fit.ortho <- with(simdata, sphereGLM(X=X, Y=Y, orthogonal=TRUE, gamma=9999))
  
  # fit.nonortho %>% summary
  # fit.ortho %>% summary
  # fit.nonortho %>% plot
  # fit.ortho %>% plot
  
  # simdata.list[[jj]] <- simdata
  fit.list[[jj]] <- list(nonortho = fit.nonortho, ortho = fit.ortho)
}










if(FALSE){
  
  
  matrix.scale <- function(X){
    sqrt.mat <- t(chol(X))
    solve(sqrt.mat)
  }
  
  simdata <- sim.sphereGLM(n=100, p=1, q=3, mu=c(0,0,20), s=10, s0=0, type="vMF",
                           seed.UDV=sample(100,1), seed.E=123, qr.U=FALSE, qr.V=TRUE)
  fit.nonortho <- try(with(simdata, sphereGLM(X=X, Y=Y, orthogonal=FALSE, lambda=1e-12)))
  fit.ortho <- with(simdata, sphereGLM(X=X, Y=Y, orthogonal=TRUE, gamma=9999))
  
  
  
  FisherInfoMatrix( simdata$X, rbind(simdata$mu, t(simdata$B)) )
  
  
  inv.sqrt.mat <- matrix.scale( Matrix::bdiag( FisherInfoMatrix(simdata$X, rbind(simdata$mu, t(simdata$B))) ) )
  
  Matrix::bdiag( summary(fit.ortho, orthogonal=TRUE)$fisher ) %>% { inv.sqrt.mat %*% . %*% t(inv.sqrt.mat) }
  Matrix::bdiag( summary(fit.nonortho, orthogonal=FALSE)$fisher ) %>% { inv.sqrt.mat %*% . %*% t(inv.sqrt.mat) }
  
  
}




if(FALSE){
  
  
  Matrix::bdiag( FisherInfoMatrix(simdata$X, rbind(simdata$mu, t(simdata$B))) ) %>% {
    # N <- nrow(fit.ortho$Y)
    c(
      sum(abs(Matrix::bdiag( summary(fit.ortho, orthogonal=TRUE)$fisher ) - .)),
      sum(abs(Matrix::bdiag( summary(fit.nonortho, orthogonal=F)$fisher ) - .))
    )
  }
  
  
  
  summary(fit.ortho, orthogonal=TRUE)$fisher %>% lapply(function(x) eigen(x)$values)
  summary(fit.nonortho, orthogonal=FALSE)$fisher %>% lapply(function(x) eigen(x)$values)
  
  #
  #
  
  summary(fit.list[[jj]]$ortho, orthogonal=TRUE)$fisher %>% lapply(function(x) x %>% eigen %>% .$values %>% {max(.)/min(.)})
  summary(fit.list[[jj]]$nonortho, orthogonal=FALSE)$fisher %>% lapply(function(x) x %>% eigen %>% .$values %>% {max(.)/min(.)})
  
  #
  #
  
  mean(sapply(fit.list, function(FIT) summary(FIT$ortho, scale=FALSE)$pvalue[2])<0.05)
  mean(sapply(fit.list, function(FIT) summary(FIT$nonortho, scale=FALSE)$pvalue[2])<0.05, na.rm=TRUE)
  
  
  
}