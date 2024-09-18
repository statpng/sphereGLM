#' @importFrom Directional rvmf
#' @export sim.sphereGLM
sim.sphereGLM <- function(n=50, p=1, q=3, mu=c(0,0,1), snr=NULL, s=2, s0=0, type=c("vMF", "Proj", "ExpMap"), seed.UDV=1, seed.E=NULL, qr.U=FALSE, qr.V=FALSE){
  # typeA = "sparse", typeB = "all", n = 50, p = 10, q = 10, 
  #  d = 3, rvec = NULL, nuA = 0.2, nuB = 0.5, d0 = 3, es = "1", 
  #  es.B = 1, snr = 1, simplify = TRUE, sigma = NULL, rho_X = 0.5, 
  #  rho_E = 0
  
  if(FALSE){
    with( sim.sphereGLM(p=1, mu=c(0,0,50), s=200), plot.sphere(Y))
    with( sim.sphereGLM(p=1, mu=c(0,0,1), snr=200, s=5, type="Proj"), plot.sphere(Y))
  }
  
  if(FALSE){
    n=50; p=1; q=3; mu=c(0,0,1); snr=10; s=2; s0=0; use.ExpMap=FALSE
  }
  
  
  
  if(is.null(seed.E)){
    seed.E <- sample( setdiff(1:1000, seed.UDV), 1 )
  }
  
  
  
  if(length(type)>1){
    type <- "Proj"
  }
  
  set.seed(seed.UDV)
  D <- sapply(1:p, function(k) (s/k^2 + s0))
  
  V <- do.call("cbind", lapply(1:p, function(x) rnorm(q)))
  if(qr.V){
    # original_lengths <- sqrt(colSums(V^2))
    V <- qr.Q(qr( cbind(mu, V) ))[,-1] #%*% diag(original_lengths, ncol(V), ncol(V))
  } else {
    V <- apply(V, 2, function(v) v/norm(v, "2"))
  }
  V <- V %*% diag(D,p,p)
  
  U <- do.call("cbind", lapply(1:p, function(x) rnorm(n)))
  if(qr.U) U <- qr.Q(qr( U ))
  
  # U <- matrix(0,n,p)
  # for( k in 1:p ){
  #   U[,k] <- rnorm(n, mean=0, sd=D[k])
  # }
  
  
  Theta0 <- tcrossprod(U, V)
  Theta <- tcrossprod(rep(1,n), mu) + Theta0
  
  
  
  
  set.seed(seed.E)
  
  E <- NULL
  if( type %in% c("ExpMap", "Proj") ){
    if(is.null(snr)){
      snr <- 10
    }
    
    E0 <- matrix(rnorm(n*(q-1),0,1),n,q-1)
    sigma <- sqrt(sum(as.numeric(Theta0)^2)/sum(as.numeric(E0)^2)/snr)
    E0 <- E0 * sigma
    
    E.basis <- qr.Q(qr( cbind(mu, do.call("cbind", lapply(1:(q-1), function(x) rnorm(q))) ) ))[,-1]
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
    #   plot.sphere(Yvmf, opacity = FALSE, add=FALSE, cex=1)
    #   plot.sphere(Yproj, opacity = FALSE, add=TRUE, col="blue")
    #   rgl::points3d(Theta+E, add=TRUE)
    # }
    
  }
  
  
  params <- c(n=n, q=q, p=p, mu=mu, snr=snr, s=s, s0=s0, type=type)
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
















#' @export orthogonalize_vectors
orthogonalize_vectors <- function(mu, vectors) {
  # Ensure mu is a column vector
  mu <- as.matrix(mu)
  if (ncol(mu) != 1) mu <- t(mu)
  
  # Combine mu and vectors into a matrix
  X <- cbind(mu, vectors)
  
  # Perform QR decomposition
  qr_result <- qr(X)
  Q <- qr.Q(qr_result)
  
  # Extract the orthogonalized vectors (excluding the first column which corresponds to mu)
  orthogonalized <- Q[, -1, drop = FALSE]
  
  # Rescale the orthogonalized vectors to preserve their original lengths
  original_lengths <- sqrt(colSums(vectors^2))
  rescaled <- orthogonalized %*% diag(original_lengths)
  
  return(rescaled)
}


