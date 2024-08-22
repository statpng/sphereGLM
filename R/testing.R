if(FALSE){
  
  library(sphereGLM)
  detach("package:sphereGLM", unload=TRUE)
  
  
  {
    set.seed(1)
    simdata1 <- sphereGLM::sim.sphereGLM(n=50, p=1, q=3, mu=c(0,1000,0), snr=50, s=50, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1])
    fit1 <- with(simdata1, sphereGLM(X=X, Y=Y, use.nlm=F))
    fit <- fit1
    
    
    summary.sphereGLM(fit1)
    
    plot.sphereGLM(fit1)
    
  }
  
  
  {
    Fit.list <- NULL
    for( ii in 1:200 ){
      if( ii %% 10 == 0 ) print(ii)
      simdata <- sim.sphereGLM(n=50, p=1, q=3, mu=c(0,50,0), snr=50, s=100, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1], seed.UDV=1, seed.E=ii)
      
      fit <- with(simdata, sphereGLM(X=X, Y=Y, use.nlm=F))
      Fit.list[[ii]] <- fit
    }
    
    
    sapply(Fit.list, function(x) summary.sphereGLM(x)$Wald[2] ) %>% hist
    
    
    #
    #
    
    
    
    Fit.list2 <- NULL
    for( ii in 1:50 ){
      if( ii %% 10 == 0 ) print(ii)
      simdata <- sim.sphereGLM(n=300, p=1, q=3, mu=c(0,50,0), snr=50, s=100, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1], seed.UDV=1, seed.E=ii)
      
      fit <- with(simdata, sphereGLM(X=X, Y=Y, use.nlm=F))
      Fit.list2[[ii]] <- fit
    }
    
    
    sapply(Fit.list2, function(x) summary.sphereGLM(x)$Wald[2] ) %>% hist
    
    #
    #
    #
    #
    
    plot.sphereGLM(fit2)
    
    fit1$testing$beta %>% apply(1,norm,"2")
    fit2$testing$beta %>% apply(1,norm,"2")
    
  }
  #
  
  fit1$beta
  fit1$testing$Fnj[[2]] %>% solve
  fit1$testing
  fit2$testing
  
  
  
  
  plot.sphereGLM(fit1)
  #
  #
  #
  #
  #
  
  
  fit2 <- with(simdata, sphereGLM(X=X, Y=Y, use.nlm = FALSE))
  plot.sphereGLM(fit1, plot.mu=TRUE, col="blue")
  plot.sphereGLM(fit2, plot.mu=TRUE, add = T, col="red")
  
  
}