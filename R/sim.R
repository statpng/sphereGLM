#' @importFrom Directional rvmf
#' @export sim.sphereGLM
sim.sphereGLM <- function(n=50, p=3, r=1, mu=c(0,0,1), snr=NULL, s=2, s0=0, type=c("vMF", "Proj", "ExpMap")){
  # typeA = "sparse", typeB = "all", n = 50, p = 10, q = 10, 
  #  d = 3, rvec = NULL, nuA = 0.2, nuB = 0.5, d0 = 3, es = "1", 
  #  es.B = 1, snr = 1, simplify = TRUE, sigma = NULL, rho_X = 0.5, 
  #  rho_E = 0
  
  if(FALSE){
    with( sim.sphereGLM(r=1, mu=c(0,0,50), s=200), png.sphere(Y))
    with( sim.sphereGLM(r=1, mu=c(0,0,1), snr=200, s=5, type="Proj"), png.sphere(Y))
  }
  
  if(FALSE){
    n=50; p=3; r=1; mu=c(0,0,1); snr=10; s=2; s0=0; use.ExpMap=FALSE
  }
  
  
  if(length(type)>1){
    type <- "Proj"
  }
  
  V.tmp <- do.call("cbind", lapply(1:r, function(x) rnorm(p)))
  V <- qr.Q(qr( cbind(mu, V.tmp) ))[,-1]
  D <- sapply(1:r, function(k) (s/k^2 + s0))
  
  U.tmp <- do.call("cbind", lapply(1:r, function(x) rnorm(n)))
  U <- qr.Q(qr( U.tmp )) %*% diag(D,r,r)
  
  # U <- matrix(0,n,r)
  # for( k in 1:r ){
  #   U[,k] <- rnorm(n, mean=0, sd=D[k])
  # }
  
  
  Theta0 <- tcrossprod(U, V)
  Theta <- tcrossprod(rep(1,n), mu) + Theta0
  
  
  
  E <- NULL
  if( type %in% c("ExpMap", "Proj") ){
    if(is.null(snr)){
      snr <- 10
    }
    
    E0 <- matrix(rnorm(n*(p-1),0,1),n,p-1)
    sigma <- sqrt(sum(as.numeric(Theta0)^2)/sum(as.numeric(E0)^2)/snr)
    E0 <- E0 * sigma
    
    E.basis <- qr.Q(qr( cbind(mu, do.call("cbind", lapply(1:(p-1), function(x) rnorm(p))) ) ))[,-1]
    E <- E0 %*% t(E.basis)
    
    if(FALSE){
      rgl.sphgrid(radaxis=F, radlab=F)
      rgl::points3d(E, col="red")
      rgl::points3d(tcrossprod(rep(1,n),mu)+Theta0+E)
    }
    
  }
  
  
  
  
  if( type == "ExpMap" ){
    
    Y <- Expmu(mu, Theta0 + E )
    
  } else if( type == "Proj" ){
    
    Y <- (Theta + E) %>% apply(1, function(theta) theta/norm(theta,"2")) %>% t
    
  } else if( type == "vMF" ){
    
    Y <- Theta %>% apply(1, function(theta) rvmf(1, theta, k=norm(theta,"2")) ) %>% t
    
    # Yproj <- (Theta + E) %>% apply(1, function(x) x/norm(x,"2")) %>% t
    # Yvmf <- Theta %>% apply(1, function(x) rvmf(1, x, k=norm(x,"2")) ) %>% t
    # {
    #   png.sphere(Yvmf, opacity = FALSE, add=FALSE, cex=1)
    #   png.sphere(Yproj, opacity = FALSE, add=TRUE, col="blue")
    #   rgl::points3d(Theta+E, add=TRUE)
    # }
    
  }
  
  
  params <- c(n=n, p=p, r=r, mu=mu, snr=snr, s=s, s0=s0, type=type)
  result <- list(mu=mu, X=U, B=V, Y=Y, Theta=Theta, E=E, params=params)
  
  result
}




#' @export sim.sphere.runif
sim.sphere.runif <- function(n){
  aa <- seq(0, 360, length.out = n)
  bb <- seq(0, 180, length.out = n)
  grid <- expand.grid(theta = aa, phi = bb)
  
  x <- sin(grid$theta * pi / 180) * cos(grid$phi * pi / 180)
  y <- sin(grid$theta * pi / 180) * sin(grid$phi * pi / 180)
  z <- cos(grid$theta * pi / 180)
  
  list(x=x,y=y,z=z)
}



#' @export sim.sphere.circle
sim.sphere.circle <- function(n){
  
  ndata = n
  theta = seq(from=-0.99,to=0.99,length.out=ndata)*pi
  tmpx  = cos(theta) + rnorm(ndata,sd=0.1)
  tmpy  = sin(theta) + rnorm(ndata,sd=0.1)
  
  data <- cbind.data.frame(x=tmpx, y=tmpy, z=1)
  data <- t( apply(data, 1, function(x) x/norm(x,"2")) )
  
  # data = riemfactory(data, name="sphere")
  
  return(data)
}






#' @export sim.sphere.nonlinear
sim.sphere.nonlinear <- function(n){
  
  ndata = n
  theta = seq(from=-0.99,to=0.99,length.out=ndata)*pi
  tmpx  = cos(theta) + rnorm(ndata,sd=0.1)
  tmpy  = sin(theta) + rnorm(ndata,sd=0.1)
  
  data <- cbind.data.frame(x=tmpx, y=tmpy, z=1)
  data <- t( apply(data, 1, function(x) x/norm(x,"2")) )
  
  # data = riemfactory(data, name="sphere")
  
  return(data)
}
















