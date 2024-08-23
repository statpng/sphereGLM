if(FALSE){
  library(sphereGLM)
  
  
  # 1. Simulation data ----
  set.seed(1)
  simdata <- sphereGLM::sim.sphereGLM(n=50, p=1, q=3, mu=c(0,100,0), snr=50, s=100, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1])
  
  
  # the estimated coefficient vectors (beta_j) are orthogonal (or non-orthogonal) to the mean vector (mu)
  fit1 <- with(simdata, sphereGLM(X=X, Y=Y, orthogonal=TRUE))
  fit2 <- with(simdata, sphereGLM(X=X, Y=Y, orthogonal=FALSE))
  
  plot.sphereGLM(fit1)
  plot.sphereGLM(fit2)
  
  # Color mapping:
  # - Larger values: Blue
  # - Smaller values: Red
  # - Intermediate values: Purple
  
  
  
  
  # 2. Real data ----
  
  RData <- R.utils::loadToEnv("/Users/png/Library/CloudStorage/Dropbox/1. JSK/3. GLM using vMF/RealData/FinalData.RData")
  
  
  print( RData$idx$Aligned )
  
  idx <- 1 # 1~5  based on length(RData$idx$Aligned)
  # >> the 5th spoke seems to look better than the others
  
  
  group <- RData$group
  id <- RData$idx$Aligned[idx] # spoke
  
  X <- RData$connectionsLengths[group==1,id]
  Y <- RData$framesBasedOnParentsUnitQuaternion[group==1,,id]
  table( apply(Y,1,norm,"2") )
  # >> All have the unit norm
  
  print( apply(Y, 2, sd) ) 
  # >> The fourth dimension has the smallest variation >> visualize the data without it.
  
  Ynew <- t( apply(Y[,1:3], 1, function(x) x/norm(x, "2")) )
  png.sphere(Ynew)
  
  
  fit1 <- sphereGLM(X=X, Y=Ynew, orthogonal=TRUE) # this takes about 5.5 seconds
  fit2 <- sphereGLM(X=X, Y=Ynew, orthogonal=FALSE) # this takes about 2.6 seconds
  
  plot.sphereGLM(fit1)
  plot.sphereGLM(fit2)
  
  
  
  
}
