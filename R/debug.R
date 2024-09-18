



function(){
  
  library(png.Directional)
  
  Directional::vmf.reg()
  
  y <- rvmf(150, rnorm(3), 5)
  a1 <- iag.reg(y, iris[, 4])
  a2 <- iag.reg(y, iris[, 4:5])
  
  b1 <- vmf.reg(y, iris[, 4])
  b2 <- vmf.reg(y, iris[, 3:4])
  
  simdata1 <- sim.vMFglm(n=100, p=1, q=3, sd=50, mu=c(0,0,100), orthogonal=T)
  vmf.reg(simdata1$Y, simdata1$X)$beta
  vmf.reg2(simdata1$Y, simdata1$X)$beta
  iag.reg(simdata1$Y, simdata1$X)$beta
  sipc.reg(simdata1$Y, simdata1$X)$beta
  esag.reg(simdata1$Y, simdata1$X)$beta
  sespc.reg(simdata1$Y, simdata1$X, lati=30, longi=30)$beta
  
  
  
  
  
  #
  
  #
  
  
  
  fit <- simdata1 %>% { glm_vmf_FixedMean_Offset3(X=.$X, Y=.$Y, eps=1e-6) }
  fit$beta
  simdata1$B
  
  fit <- glm_vmf_FixedMean_Offset3(X=iris[,3:4,drop=F], Y=y, eps=1e-16)
  b2$beta
  fit$beta
  
  
  # 1: p=1, sd=50, mu=(0,0,100)
  function(){
    devtools::load_all()
    library(dplyr)
    
    simdata1 <- sim.vMFglm(n=100, p=1, q=3, sd=50, mu=c(0,0,100), orthogonal=T)
    # simdata1$Y %>% png.sphere()
    fit1 <- simdata1 %>% { glm_vmf_FixedMean_Offset2(X=.$X, Y=.$Y, lambda=1e-4, eps=1e-16) }
    list(simdata1, fit1) %>% {
      rbind(.[[1]]$B, .[[2]]$beta)
    }
    
  }
  
  
  # 2: p=1, sd=50, mu=(0,100,100)
  function(){
    
    simdata1 <- sim.vMFglm(n=100, p=1, q=3, sd=50, mu=c(0,100,100), orthogonal=T)
    # simdata1$Y %>% png.sphere()
    fit1 <- simdata1 %>% { glm_vmf_FixedMean_Offset2(X=.$X, Y=.$Y, lambda=1e-4, eps=1e-16) }
    list(simdata1, fit1) %>% {
      rbind(.[[1]]$B, .[[2]]$beta)
    }
    
  }
  
  
  # 3: p=1, sd=15, mu=(0,100,100)
  function(){
    
    simdata3 <- sim.vMFglm(n=100, p=1, q=3, sd=15, mu=c(0,100,100), orthogonal=T)
    # simdata3$Y %>% png.sphere()
    fit3 <- simdata3 %>% { glm_vmf_FixedMean_Offset2(X=.$X, Y=.$Y, lambda=1e-4, eps=1e-16) }
    list(simdata3, fit3) %>% {
      rbind(.[[1]]$B, .[[2]]$beta)
    }
    
  }
  
  
  # 4: p=1, sd=5, mu=(0,100,100)
  function(){
    
    simdata4 <- sim.vMFglm(n=100, p=1, q=3, sd=5, mu=c(0,100,100), orthogonal=T)
    # simdata4$Y %>% png.sphere()
    fit4 <- simdata4 %>% { glm_vmf_FixedMean_Offset2(X=.$X, Y=.$Y, lambda=1e-4, eps=1e-16) }
    list(simdata4, fit4) %>% {
      rbind(.[[1]]$B, .[[2]]$beta)
    }
    
    list(as.vector(simdata4$B), fit4$beta[2,]) %>% 
      lapply(function(x) x/norm(x,"2"))
    
  }
  
  
  
  # 5: p=1, sd=1, mu=(0,100,100)
  sim5 <- function(){
    # library(png.utils)
    # detach("package:png.utils", unload=TRUE)
    library(tidyverse)
    
    sim_func <- function(n, sd, tau){
      sim.vMFglm(n=n, p=1, q=3, sd=sd, mu=c(0,0,tau), orthogonal=T)
    }
    
    fit_func <- function(simdata){
      glm_vmf_FixedMean_Offset2(X=simdata$X, Y=simdata$Y, lambda=1e-6, eps=1e-12)
    }
    
    result_func <- function(simdata, fit){
      err1 <- mean( (simdata$mu - fit$beta[1,])^2 )
      err2 <- mean( (simdata$B - fit$beta[2,])^2 )
      c(err1, err2)
    }
    
    out.df <- png.utils::png.sim.run(
      GRID=expand.grid(n=c(50, 100), sd=c(5, 10, 50), tau=c(10, 50, 100), iter=1:10), 
      sim_func=sim_func, fit_func=fit_func, result_func=result_func, n.results=2, columns=c("mu", "beta"))
    
    out.df %>% 
      ggplot() + 
      geom_boxplot(aes(n, value, group=n)) + 
      facet_grid(type~sd+tau, scales="free", labeller=png.labeller()) +
      ylab("Mean Squared Error")
    
    filename <- png.utils::png.sim.filename(GRID, title="glm_vmf", measure="MSE")
    ggsave(filename=filename, width=10, height=5)
    
    out.df %>% reshape2::dcast(sd+tau~type+paste0("n=",n), mean, margins="n", value.var="value")
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  library(dplyr)
  devtools::load_all()
  
  simdata1 <- sim.vMFglm(n=100, p=1, q=3, sd=50, mu=c(0,0,100), orthogonal=T)
  # simdata1$Y %>% png.sphere()
  fit <- glm_vmf_FixedMean_Offset2(X=simdata1$X, Y=simdata1$Y, lambda=1e-4, eps=1e-16)
  rbind(simdata1$B, fit$beta)
  
  
  simdata2 <- sim.vMFglm(n=50, p=1, q=3, sd=50, mu=c(0,0,100), orthogonal=T)
  # simdata2$Y %>% png.sphere()
  fit <- glm_vmf_FixedMean_Offset2(X=simdata2$X, Y=simdata2$Y, lambda=1e-6, eps=1e-16)
  rbind(simdata2$B, fit$beta)
  
  
  simdata3 <- sim.vMFglm(n=50, p=1, q=3, sd=50, mu=c(0,0,500), orthogonal=T)
  # simdata2$Y %>% png.sphere()
  fit <- glm_vmf_FixedMean_Offset2(X=simdata3$X, Y=simdata3$Y, lambda=0, eps=1e-16)
  rbind(simdata3$B, fit$beta)
  
  
  
  simdata4 <- sim.vMFglm(n=50, p=1, q=3, sd=1000, mu=c(0,0,500), orthogonal=T)
  # simdata4$Y %>% png.sphere()
  fit <- glm_vmf_FixedMean_Offset2(X=simdata4$X, Y=simdata4$Y, lambda=0, eps=1e-16)
  rbind(simdata4$B, fit$beta)
  
  
  simdata5 <- sim.vMFglm(n=50, p=2, q=4, sd=500, mu=c(0,0,0,100), orthogonal=T)
  # simdata5$Y %>% png.sphere()
  fit <- glm_vmf_FixedMean_Offset2(X=simdata5$X, Y=simdata5$Y, lambda=0, eps=1e-16)
  rbind(simdata5$B, fit$beta)
  
  
  
  simdata6 <- sim.vMFglm(n=50, p=1, q=3, sd=5, mu=c(0,0,50), orthogonal=T)
  # simdata6$Y %>% png.sphere()
  fit <- glm_vmf_FixedMean_Offset2(X=simdata5$X, Y=simdata5$Y, lambda=0, eps=1e-16)
  rbind(simdata5$B, fit$beta)
}







function(){
  
  devtools::load_all()
  
  
  set.seed(1)
  out.beta <- array(NA, c(3, 4, 5, 100))
  out <- array(NA, c(1, 3, 5, 100))
  for( j in 1:dim(out)[3] ){
    print(j)
    
    for( i in 1:dim(out)[4] ){
      
      tau <- c(5,10,20,50,100)[j]
      simdata <- sim.vMFglm(n=50, p=1, q=3, sd=100, mu=c(0,0,100), orthogonal=TRUE)  # png.sphere(simdata$Y)
      
      fit.GeodRegr <- with(simdata, GeodRegr::geo_reg("sphere", X, t(Y), estimator="l2", max_iter=1e+2))
      fit.vmf <- with(simdata, glm_vmf_FixedMean_Offset(X=X, Y=Y, MU=c(0,0,tau), lambda=1e-2, maxit=100, eps=1e-8, standardize=TRUE))
      
      fit.vmf <- glm_vmf_FixedMean_Offset2(X=simdata$X, Y=simdata$Y, lambda=1e-12, maxit=100, eps=1e-16, standardize=TRUE)
      
      
      fit.vmf$beta
      simdata$mu
      simdata$B
      
      # install.packages("Directional")
      # library(Directional)
      # y <- rvmf(150, rnorm(3), 5)
      # b1 <- vmf.reg(y, iris[, 4])
      # b2 <- glm_vmf_FixedMean(iris[,4,drop=F], y)
      # b1$beta
      # rbind(b2$mu, b2$beta)
      # GeodRegr::geo_reg(manifold="sphere", x=as.matrix(iris[,4,drop=F]), y=t(y), estimator="l2") %>% {
      #   t(cbind(.$p, .$V))
      # }
      
      
      
      result <- list(simdata$B, 
                     t(fit.GeodRegr$V), 
                     fit.vmf$beta) %>% 
        lapply(t) %>% lapply(function(x) x/norm(x,"F"))
      
      result
      
      
      out[,1,j,i] <- png.utils::png.angle(result[[1]], result[[2]])$max
      out[,2,j,i] <- png.utils::png.angle(result[[1]], result[[3]])$max
      # out[,3,j,i] <- png.utils::png.angle(result[[1]], result[[4]])$max
      
      out.beta[,1,j,i] <- simdata$B
      out.beta[,2,j,i] <- t(fit.GeodRegr$V)
      out.beta[,3,j,i] <- fit1$beta
      # out.beta[,4,j,i] <- fit2$beta[-1,]
    }
  }
  
  if(FALSE){
    save(out, out.beta, file="out.glm_vmf1.RData")
  }
  
  load("out.glm_vmf1.RData")
  
  
  apply(out,2:3,mean)
  #          [,1]     [,2]     [,3]     [,4]     [,5]
  # [1,] 10.68148 5.796712 3.926280 4.011802 4.258179
  # [2,] 12.87231 7.023528 4.696356 3.490202 2.928996
  
  
  
  for( idx.sd in 1:dim(out.beta)[3] ){
    pdf(file=paste0("glm_vmf1-ScatterPlot", idx.sd, ".pdf"), width=5, height=5)
    
    par(mfrow=c(2,2), omi=c(.2,.4,.2,.2), mai=c(.8,.6,.2,.2))
    cbind.data.frame(true=out.beta[1,1,idx.sd,], est=out.beta[1,2,idx.sd,]) %>% {plot(., pch=18); cor(.) %>% {.[2]} %>% round(4) %>% mtext(., side=3, line=-1, adj=.01)}
    mtext("GeodRegr", 3, line=.5)
    mtext("Y1", 2, line=4.5, outer=F)
    cbind.data.frame(true=out.beta[1,1,idx.sd,], est=out.beta[1,3,idx.sd,]) %>% {plot(., pch=18); cor(.) %>% {.[2]} %>% round(4) %>% mtext(., side=3, line=-1, adj=.01)}
    mtext("IRLS with vMF", 3, line=.5)
    cbind.data.frame(true=out.beta[1,1,idx.sd,], est=out.beta[1,2,idx.sd,]) %>% {plot(., pch=18); cor(.) %>% {.[2]} %>% round(4) %>% mtext(., side=3, line=-1, adj=.01)}
    mtext("Y2", 2, line=4.5, outer=F)
    cbind.data.frame(true=out.beta[1,1,idx.sd,], est=out.beta[1,3,idx.sd,]) %>% {plot(., pch=18); cor(.) %>% {.[2]} %>% round(4) %>% mtext(., side=3, line=-1, adj=.01)}
    
    dev.off()
  }
  
  
  
  
  
  
}
