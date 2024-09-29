if(FALSE){
  
  
  
  
  # TwoSampleTest ----
  
  library(sphereGLM)
  library(dplyr)
  
  RData <- R.utils::loadToEnv("/Users/png/Library/CloudStorage/Dropbox/1. JSK/3. GLM using vMF/RealData/FinalData.RData")
  
  group <- RData$group
  idx <- RData$idx$Aligned
  
  
  
  
  Y <- RData$framesBasedOnParentsUnitQuaternion[,,1] %>% .[,1:3] %>% apply(1, function(x) x/norm(x, "2")) %>% t
  X <- group
  
  fit <- sphereGLM(X=X, Y=Y, orthogonal=TRUE, lambda=1)
  
  
  
  
  
}