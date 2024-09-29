if(FALSE){
  
  library(rgl)
  df <- data.frame(x=rnorm(100),
                   y=rnorm(100),
                   z=rnorm(100))
  plot3d(df)
  grid3d('z')
  
  df2 <- t(apply(df,1,function(x) x/norm(x,"2"))) %>% as.data.frame()
  df2 %>% plot.sphere(opacity = FALSE)
  
  
  cross_product(mu, mu)
  
  show2d({
    par(mar=c(0,0,0,0))
    plot(x = df2$x, y = df2$y, 
         col = "black")
  })
  
  
}