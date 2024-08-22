if(FALSE){
  
  
  
  
  
  
  
  
  
  load("/Users/png/Library/CloudStorage/Dropbox/1. JSK/3. GLM using vMF/RealData/Sungkyu_Data_2024/LP_DS_Rep_Data.Rdata")
  
  library(dplyr)
  
  spokesLengths_G1<-LP_DS_Rep_Data$spokesLengths_G1
  spokesLengths_G2<-LP_DS_Rep_Data$spokesLengths_G2
  
  connectionsLengths_G1<-LP_DS_Rep_Data$connectionsLengths_G1
  connectionsLengths_G2<-LP_DS_Rep_Data$connectionsLengths_G2
  
  spokesDirectionsBasedOnFrames_G1<-LP_DS_Rep_Data$spokesDirectionsBasedOnFrames_G1
  spokesDirectionsBasedOnFrames_G2<-LP_DS_Rep_Data$spokesDirectionsBasedOnFrames_G2
  
  connectionsDirectionsBasedOnParentFrames_G1<-LP_DS_Rep_Data$connectionsDirectionsBasedOnParentFrames_G1
  connectionsDirectionsBasedOnParentFrames_G2<-LP_DS_Rep_Data$connectionsDirectionsBasedOnParentFrames_G2
  
  framesBasedOnParentsUnitQuaternion_G1<-LP_DS_Rep_Data$framesBasedOnParentsUnitQuaternion_G1
  framesBasedOnParentsUnitQuaternion_G2<-LP_DS_Rep_Data$framesBasedOnParentsUnitQuaternion_G2
  
  
  # Quarternion mapping by treating negative sign
  {
    
    idx <- 5
    framesBasedOnParentsUnitQuaternion_G1[idx,,] %>% t %>% pairs()
    framesBasedOnParentsUnitQuaternion_G1[idx,2,] <- -abs(framesBasedOnParentsUnitQuaternion_G1[idx,2,])
    framesBasedOnParentsUnitQuaternion_G1[idx,4,] <- -abs(framesBasedOnParentsUnitQuaternion_G1[idx,4,])
    
    
    idx <- 8
    framesBasedOnParentsUnitQuaternion_G1[idx,,] %>% t %>% pairs()
    framesBasedOnParentsUnitQuaternion_G1[idx,2,] <- -abs(framesBasedOnParentsUnitQuaternion_G1[idx,2,])
    framesBasedOnParentsUnitQuaternion_G1[idx,4,] <- -abs(framesBasedOnParentsUnitQuaternion_G1[idx,4,])
    
    idx <- 11
    framesBasedOnParentsUnitQuaternion_G1[idx,,] %>% t %>% pairs()
    framesBasedOnParentsUnitQuaternion_G1[idx,2,] <- abs(framesBasedOnParentsUnitQuaternion_G1[idx,2,])
    framesBasedOnParentsUnitQuaternion_G1[idx,4,] <- abs(framesBasedOnParentsUnitQuaternion_G1[idx,4,])
    
    idx <- 40
    framesBasedOnParentsUnitQuaternion_G1[idx,,] %>% t %>% pairs()
    framesBasedOnParentsUnitQuaternion_G1[idx,2,] <- -abs(framesBasedOnParentsUnitQuaternion_G1[idx,2,])
    framesBasedOnParentsUnitQuaternion_G1[idx,4,] <- -abs(framesBasedOnParentsUnitQuaternion_G1[idx,4,])
    
    idx <- 42
    framesBasedOnParentsUnitQuaternion_G1[idx,,] %>% t %>% pairs()
    framesBasedOnParentsUnitQuaternion_G1[idx,2,] <- abs(framesBasedOnParentsUnitQuaternion_G1[idx,2,])
    framesBasedOnParentsUnitQuaternion_G1[idx,4,] <- abs(framesBasedOnParentsUnitQuaternion_G1[idx,4,])
    
  }
  
  
  
  {
    
    idx <- 5
    framesBasedOnParentsUnitQuaternion_G2[idx,,] %>% t %>% pairs()
    framesBasedOnParentsUnitQuaternion_G2[idx,2,] <- -abs(framesBasedOnParentsUnitQuaternion_G2[idx,2,])
    framesBasedOnParentsUnitQuaternion_G2[idx,4,] <- -abs(framesBasedOnParentsUnitQuaternion_G2[idx,4,])
    
    
    idx <- 8
    framesBasedOnParentsUnitQuaternion_G2[idx,,] %>% t %>% pairs()
    framesBasedOnParentsUnitQuaternion_G2[idx,2,] <- -abs(framesBasedOnParentsUnitQuaternion_G2[idx,2,])
    framesBasedOnParentsUnitQuaternion_G2[idx,4,] <- -abs(framesBasedOnParentsUnitQuaternion_G2[idx,4,])
    
    idx <- 11
    framesBasedOnParentsUnitQuaternion_G2[idx,,] %>% t %>% pairs()
    framesBasedOnParentsUnitQuaternion_G2[idx,2,] <- abs(framesBasedOnParentsUnitQuaternion_G2[idx,2,])
    framesBasedOnParentsUnitQuaternion_G2[idx,4,] <- abs(framesBasedOnParentsUnitQuaternion_G2[idx,4,])
    
    idx <- 40
    framesBasedOnParentsUnitQuaternion_G2[idx,,] %>% t %>% pairs()
    framesBasedOnParentsUnitQuaternion_G2[idx,2,] <- -abs(framesBasedOnParentsUnitQuaternion_G2[idx,2,])
    framesBasedOnParentsUnitQuaternion_G2[idx,4,] <- -abs(framesBasedOnParentsUnitQuaternion_G2[idx,4,])
    
    idx <- 42
    framesBasedOnParentsUnitQuaternion_G2[idx,,] %>% t %>% pairs()
    framesBasedOnParentsUnitQuaternion_G2[idx,2,] <- abs(framesBasedOnParentsUnitQuaternion_G2[idx,2,])
    framesBasedOnParentsUnitQuaternion_G2[idx,4,] <- abs(framesBasedOnParentsUnitQuaternion_G2[idx,4,])
    
  }
  
  
  # Note that the 16th observation has a single value.
  {
    framesBasedOnParentsUnitQuaternion_G2[16,,] %>% t %>% pairs(xlim=c(-1,1), ylim=c(-1,1))
  }
  
  
  
  
  
  
  
  
  
  
  library(dplyr)
  library(sphereGLM)
  detach("package:sphereGLM", unload=TRUE)
  
  
  j <- 2
  
  X <- connectionsLengths_G1[j,]
  Y <- framesBasedOnParentsUnitQuaternion_G1[j,,] %>% t
  
  
  fit1 <- sphereGLM(X=X, Y=Y, use.nlm=TRUE, maxit=10)
  fit2 <- sphereGLM(X=X, Y=Y, use.nlm=FALSE)
  system.time( fit3 <- png.Directional::vmf.reg2(Y, X) )
  
  fit1$beta
  fit2$beta
  fit3$beta
  
  
  
  #
  #
  #
  # $beta2
  # [,1]       [,2]       [,3]       [,4]
  # [1,] 0.5035829 -0.5031596 -0.4945819 -0.4986215
  # [2,] 0.3765303 -0.3773492 -0.3633077 -0.3703704
  # 
  # $loglik
  # [1] 1143.905
  
  fit$loglik
  
  X
  
  
  
  
  {
    j <- 5
    
    X <- connectionsLengths_G1[j,]
    Y <- framesBasedOnParentsUnitQuaternion_G1[j,,] %>% t %>% {.[,1:3]} %>% apply(1,function(x) x/norm(x,"2")) %>% t
    
    system.time( fit <- sphereGLM(X=X, Y=Y, use.nlm=TRUE, maxit=10) )
    fit$beta
    plot.sphereGLM(fit, plot.mu=TRUE)
    
    
  }
  
  
  
  {
    j <- 8
    
    X <- connectionsLengths_G1[j,]
    Y <- framesBasedOnParentsUnitQuaternion_G1[j,,] %>% t %>% {.[,1:3]} %>% apply(1,function(x) x/norm(x,"2")) %>% t
    
    system.time( fit <- sphereGLM(X=X, Y=Y, use.nlm=TRUE, maxit=10) )
    fit$beta2
    plot.sphereGLM(fit, plot.mu=TRUE)
  }
  
  
  
  {
    j <- 11
    
    X <- connectionsLengths_G1[j,]
    Y <- framesBasedOnParentsUnitQuaternion_G1[j,,] %>% t %>% {.[,1:3]} %>% apply(1,function(x) x/norm(x,"2")) %>% t
    
    system.time( fit <- sphereGLM(X=X, Y=Y, use.nlm=TRUE, maxit=10) )
    fit$beta2
    plot.sphereGLM(fit, plot.mu=TRUE)
  }
  
  
  {
    j <- 40
    
    X <- connectionsLengths_G1[j,]
    Y <- framesBasedOnParentsUnitQuaternion_G1[j,,] %>% t %>% {.[,1:3]} %>% apply(1,function(x) x/norm(x,"2")) %>% t
    
    system.time( fit <- sphereGLM(X=X, Y=Y, use.nlm=TRUE, maxit=10) )
    fit$beta
    plot.sphereGLM(fit, plot.mu=TRUE)
  }
  
  
  
  {
    j <- 42
    
    X <- connectionsLengths_G1[j,]
    Y <- framesBasedOnParentsUnitQuaternion_G1[j,,] %>% t %>% {.[,1:3]} %>% apply(1,function(x) x/norm(x,"2")) %>% t
    
    system.time( fit <- sphereGLM(X=X, Y=Y, use.nlm=TRUE, maxit=10) )
    fit$beta
    plot.sphereGLM(fit, plot.mu=TRUE)
  }
  
  
  
  
  
  
  
  
  # 3d sphere ---------------------------------------------------------------
  
  
  {
    sphereGLM( X=connectionsLengths_G1[j,], Y=t(framesBasedOnParentsUnitQuaternion_G1[j,,]))
    
    X=connectionsLengths_G1[j,]; Y=t(framesBasedOnParentsUnitQuaternion_G1[j,,])
  }
  
  
  
  LP_DS_Rep_Data$connectionsDirectionsBasedOnParentFrames_G1 %>% str
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # GroupTesting ------------------------------------------------------------
  
  rm(list=ls())
  
  library(sphereGLM)
  
  RData <- R.utils::loadToEnv("/Users/png/Library/CloudStorage/Dropbox/1. JSK/3. GLM using vMF/RealData/FinalData.RData")
  
  group <- RData$group
  idx <- RData$idx$Aligned
  
  
  
  {
    i <- 3
    
    fit1 <- RData %>% {
      list(Y=.$framesBasedOnParentsUnitQuaternion[group==1,,idx[i]] %>% .[,1:3] %>% apply(1, function(x) x/norm(x, "2")) %>% t,
           X=.$connectionsLengths[group==1, idx[i]])
    } %>% with( sphereGLM(X=X, Y=Y) )
    
    fit2 <- RData %>% {
      list(Y=.$framesBasedOnParentsUnitQuaternion[group==2,,idx[i]] %>% .[,1:3] %>% apply(1, function(x) x/norm(x, "2")) %>% t,
           X=.$connectionsLengths[group==2, idx[i]])
    } %>% with( sphereGLM(X=X, Y=Y) )
    
    fit1 %>% plot.sphereGLM(plot.mu=TRUE, col="red", cex=0.008)
    fit2 %>% plot.sphereGLM(plot.mu=TRUE, col="blue", cex=0.008, add=TRUE)
    
  }
  
  
  
  {
    i <- 4
    
    fit1 <- RData %>% {
      list(Y=.$framesBasedOnParentsUnitQuaternion[group==1,,idx[i]] %>% .[,1:3] %>% apply(1, function(x) x/norm(x, "2")) %>% t,
           X=.$connectionsLengths[group==1, idx[i]])
    } %>% with( sphereGLM(X=X, Y=Y) )
    
    fit2 <- RData %>% {
      list(Y=.$framesBasedOnParentsUnitQuaternion[group==2,,idx[i]] %>% .[,1:3] %>% apply(1, function(x) x/norm(x, "2")) %>% t,
           X=.$connectionsLengths[group==2, idx[i]])
    } %>% with( sphereGLM(X=X, Y=Y) )
    
    fit1 %>% plot.sphereGLM(plot.mu=TRUE, col="red", cex=0.008)
    fit2 %>% plot.sphereGLM(plot.mu=TRUE, col="blue", cex=0.008, add=TRUE)
    
  }
  
  
  
  {
    i <- 5
    
    fit1 <- RData %>% {
      list(Y=.$framesBasedOnParentsUnitQuaternion[group==1,,idx[i]] %>% .[,1:3] %>% apply(1, function(x) x/norm(x, "2")) %>% t,
           X=.$connectionsLengths[group==1, idx[i]])
    } %>% with( sphereGLM(X=X, Y=Y) )
    
    fit2 <- RData %>% {
      list(Y=.$framesBasedOnParentsUnitQuaternion[group==2,,idx[i]] %>% .[,1:3] %>% apply(1, function(x) x/norm(x, "2")) %>% t,
           X=.$connectionsLengths[group==2, idx[i]])
    } %>% with( sphereGLM(X=X, Y=Y) )
    
    fit1 %>% plot.sphereGLM(plot.mu=TRUE, col="red", cex=0.008)
    fit2 %>% plot.sphereGLM(plot.mu=TRUE, col="blue", cex=0.008, add=TRUE)
    
  }
  
  
  
  
}