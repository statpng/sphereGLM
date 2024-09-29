#' @import rgl
#' @method plot sphereGLM
#' @export
plot.sphereGLM <- function(fit, plot.mu=TRUE, plot.conf=TRUE, plot.tangent=TRUE, use.PNS=FALSE, conf.level=0.95, scale=TRUE, cex=0.01, opacity=FALSE, col=NULL, add=FALSE, main=NULL, arrow.type="rotation", ...){
  
  # if( !arrow.type %in% c("lines", "rotation") ) stop( 'arrow.type %in% c("lines", "rotation")' )
  
  if(FALSE){
    opacity=FALSE
    plot.mu = T
    plot.mu=TRUE; cex=0.01; opacity=FALSE; col=NULL; add=FALSE; main=NULL; ...=NULL
    plot.conf=TRUE
    scale=TRUE
  }
  
  
  if(FALSE){
    
    simdata.list <- lapply(c(100, 500, 1000), function(n){
      simdata <- sim.sphereGLM(n=n, p=2, q=3, mu=c(0,0,20), s=0, s0=0, type="vMF", seed.UDV=1, qr.U=FALSE, qr.V=FALSE)
    })
    
    fit.list <- lapply(1:length(simdata.list), function(i){
      simdata <- simdata.list[[i]]
      fit <- sphereGLM(X=simdata$X, Y=simdata$Y, orthogonal=simdata$params$qr.V)
      fit
    })
    
    
    fit.list[[1]] %>% plot
    fit.list[[2]] %>% plot
    fit.list[[3]] %>% plot
    
    
  }
  
  
  if(FALSE){
    
    devtools::load_all()
    
    set.seed(1)
    simdata <- sim.sphereGLM(n=100, p=2, q=3, mu=c(0,0,20), s=c(10,5), s0=0, type="vMF", seed.UDV=1, qr.U=FALSE, qr.V=TRUE)
    
    fit <- with(simdata, sphereGLM(X=X, Y=Y, orthogonal=F))
    plot(fit)
    fit <- with(simdata, sphereGLM(X=X, Y=Y, orthogonal=T))
    plot(fit, scale=T)
    
    fit.summary <- summary( fit )
    fit.summary$stat
    
  }
  
  
  
  if(FALSE){
    
    set.seed(1)
    fit1 <- sim.sphereGLM(n=150, p=1, q=3, mu=c(0,0,10), s=5, s0=0, type="vMF", qr.U=FALSE, qr.V=TRUE) %>% with(sphereGLM(X, Y))
    set.seed(1)
    fit2 <- sim.sphereGLM(n=150, p=1, q=3, mu=c(0,0,10), s=2, s0=0, type="vMF", qr.U=FALSE, qr.V=TRUE) %>% with(sphereGLM(X, Y))
    
    plot(fit1, plot.mu = T, col="blue", opacity=F, cex=0.01, main="abc")
    plot(fit2, plot.mu = T, col="red", add=TRUE, opacity=F)
    
  }
  
  if(FALSE){
    
    set.seed(123)
    fit1 <- sim.sphereGLM(n=150, p=1, q=3, mu=c(0,0,10), s=5, s0=0, type="vMF", qr.U=FALSE, qr.V=TRUE) %>% with(sphereGLM(X, Y, orthogonal=TRUE))
    fit2 <- sim.sphereGLM(n=150, p=1, q=3, mu=c(0,0,10), s=2, s0=0, type="vMF", qr.U=FALSE, qr.V=TRUE) %>% with(sphereGLM(X, Y, orthogonal=TRUE))
    
    plot(fit1, plot.mu = T, col="blue", opacity=F, cex=0.01, main="abc")
    plot(fit2, plot.mu = T, col="red", add=TRUE, opacity=F)
    
  }
  
  
  
  
  if(FALSE){
    
    set.seed(1)
    fit1 <- sim.sphereGLM(n=150, mu=c(0,0,1), r=2, s=5, s0 = 0.01) %>% 
      with(sphereGLM(X[,1], Y))
    fit2 <- sim.sphereGLM(n=150, mu=c(0,1,0), r=2, s=5, s0 = 0.01) %>% 
      with(sphereGLM(X[,1], Y))
    
    plot.sphereGLM(fit1, plot.mu = T, col="red", opacity=F, cex=0.01, main="abc")
    plot.sphereGLM(fit2, plot.mu = T, col="blue", add=TRUE, opacity=F)
    #
    
  }
  
  
  
  
  if(FALSE){
    
    library(png.Directional)
    library(sphereGLM)
    library(dplyr)
    
    
    # Figure 1
    {
      set.seed(2)
      simdata <- sim.sphereGLM(n=150, mu=c(1/sqrt(3),1/sqrt(3),1/sqrt(3)) %>% {./norm(.,"2")}, snr=100, r=2, s=3, s0 = 0.0001, type="Proj")
      X <- simdata$X;  Y <- simdata$Y
      fit <- sphereGLM(X, Y)
      
      plot.sphereGLM(fit, plot.mu=TRUE)
      
      fit
    }
    
    
    # Figure 2
    {
      set.seed(2)
      simdata <- sim.sphereGLM(n=150, mu=c(1/sqrt(3),1/sqrt(3),1/sqrt(3)) %>% {./norm(.,"2")}, snr=100, r=1, s=10, s0 = 0.0001, type="Proj")
      X <- simdata$X;  Y <- simdata$Y
      fit <- sphereGLM(X, Y)
      
      plot.sphereGLM(fit, plot.mu=TRUE)
      
      fit
    }
    
    
    # Figure 3
    {
      set.seed(1)
      simdata1 <- sim.sphereGLM(n=150, mu=c(1/sqrt(3),1/sqrt(3),1/sqrt(3))*50, r=1, s=100, s0 = 0.0001, type="Proj")
      
      set.seed(1)
      simdata1 <- sim.sphereGLM(n=150, mu=c(1/sqrt(3),1/sqrt(3),1/sqrt(3))*50, r=1, s=100, s0 = 0.0001, type="vMF")
      set.seed(1)
      simdata2 <- sim.sphereGLM(n=150, mu=c(10/sqrt(3),10/sqrt(3),10/sqrt(3))*150, r=1, s=500, s0 = 0.0001, type="vMF")
      
      fit1 <- glm.vmf(simdata1$U, simdata1$X)
      fit2 <- glm.vmf(simdata2$U, simdata2$X)
      
      glm.vmf.plot(fit1, plot.mu=TRUE)
      glm.vmf.plot(fit2, plot.mu=FALSE)
      
      fit1
      fit2
      
    }
    
    
  }
  

  
  
  # start ----
  
  
  # use.PNS ----
  # if(use.PNS){
  #   fit.pns <- pns(t(sphere4d), output=FALSE)
  # }
  
  {
    
    X <- fit$X
    Y <- fit$Y
    mu <- fit$beta[1,,drop=T]
    beta <- fit$beta[-1,,drop=FALSE]
    
    n <- nrow(Y)
    
    norm.mu <- ifelse( abs(norm(mu, "2")-1) < 0.1, norm(mu, "2"), norm(mu, "2") * 0.8 )
    
  }
  
  
  # For better visualization
  if(scale) mu <- mu / norm.mu
  if(scale) beta <- apply(beta, 1, function(x) x / norm.mu ) %>% t()
  
  
  
  
  # rgl::planes3d(0,0,1,1, alpha=0.2)
  MU <- tcrossprod(rep(1,nrow(beta)), mu)
  
  
  if(is.null(col)){
    col1 <- "red"
    col2 <- "blue"
    col3 <- "purple"
    
    value_to_color <- function(value, min_val, max_val) {
      scaled <- (value - min_val) / (max_val - min_val)
      r <- 1 - scaled
      g <- 0
      b <- scaled
      rgb(r, g, b)
    }
    # blue: high; red: low
    
    col1 <- value_to_color(X, min(X), max(X))
    
  } else {
    
    
    {
      a1 <- col2rgb(col)
      a2 <- rgb2hsv(a1)
      pastel_col1 <- hsv(a2[1,], a2[2,]*0.9, a2[3,])
      pastel_col2 <- hsv(a2[1,], a2[2,]*0.7, a2[3,])
      pastel_col3 <- hsv(a2[1,], a2[2,]*0.5, a2[3,])
      hue <- a2["h",]*360
      n <- 0.4
      pastel_col4 <- hcl(hue, 35, 85)
    }
    
    col1 <- pastel_col1
    col2 <- pastel_col2
    col3 <- pastel_col3
  }
  
  
  
  arrow.num <- ifelse(arrow.type=="lines", 100, 3)
  
  
  
  plot.sphere(Y, col=col1, opacity=opacity, add=add, cex=cex, main=main, ...)
  
  
  rgl::spheres3d(expand.grid(replicate(3, c(-1.5,1.5), simplify=F)), alpha=0.01, col="white", radius=0.001, lit=FALSE)
  
  
  if(nrow(beta)>1){
    normal_vector <- project_to_plane(c(0,0,0), cross_product(beta[1,], beta[2,]), mu)
    rgl::planes3d(a=normal_vector, d=-sum(normal_vector * mu), alpha=0.1, add=TRUE)
    # rgl::abclines3d(mu[1],mu[2],mu[3], beta, alpha=0.5, col="blue", lwd=2)
  } else {
    
    if( fit$params$orthogonal & plot.tangent ){
      rgl::planes3d(a=mu, d=-sum(mu * mu), alpha=0.1, add=TRUE)
      
      # rgl::abclines3d(mu[1],mu[2],mu[3], beta, alpha=0.5, col="blue", lwd=2)
    }
    
  }
  
  
  
  
  
  
  for( k in 1:nrow(MU) ){
    # rgl::arrow3d(MU[k,]-beta[k,], MU[k,]+beta[k,], type=c("lines", "rotation")[1], col=col2, s=0.05, barblen=0.03, width=0.5, add=TRUE, n=100)
    rgl::arrow3d(MU[k,]-beta[k,], MU[k,]+beta[k,], type=arrow.type, col=col2, s=0.05, barblen=0.04*sqrt(norm(beta[k,], "2")), width=0.1, add=TRUE, n=arrow.num)
  }
  
  
  
  # plot.mu ----
  if(plot.mu){
    
    switch(3,
           `1`=rgl::arrow3d(c(0,0,0), mu, type=c("lines", "rotation")[1], col=col3, s=0.1, barblen=0.02, width=0.5, add=TRUE),
           `2`=rgl::arrow3d(c(0,0,0), mu, type=c("lines", "rotation")[2], col=col3, s=0.1, barblen=0.02, width=0.1, add=TRUE),
           `3`=rgl::arrow3d(c(0,0,0), mu, type=arrow.type, col=col3, s=0.02, barblen=0.02, width=0.1, add=TRUE, n=arrow.num)
    )
    
  }
  
  
  
  
  
  # plot.conf ----
  # plot.conf & !scale
  if(plot.conf){
    
    fit.summary <- summary(fit)
    df <- fit.summary$df
    p <- ncol(fit$X)
    
    
    for(j in 1:p){
      
      
      if(scale){
        S <- solve(fit.summary$Fisher$Fnj[[j+1]]) / (norm.mu^2)
      } else {
        S <- solve(fit.summary$Fisher$Fnj[[j+1]])
      }
      
      
      fit.eig <- eigen(S)
      chisq.val <- qchisq(conf.level, df = df)
      
      
      theta <- seq(0, 2*pi, length.out = 100)
      
      for( idx in list(c(1,2), c(1,3), c(2,3)) ){
        a <- sqrt(chisq.val * fit.eig$values[idx[1]])
        b <- sqrt(chisq.val * fit.eig$values[idx[2]])
        
        ellipse_data <- cbind(a * cos(theta), b * sin(theta))
        
        rotated_data <- t(fit.eig$vectors[,idx] %*% t(ellipse_data))
        rotated_data2 <- tcrossprod(rep(1,100), MU[j,]+beta[j,]) + rotated_data
        
        lines3d(rotated_data2[,1], rotated_data2[,2], rotated_data2[,3], lwd=2, 
                col=col2, barblen=0.03, width=0.5, add=TRUE, n=200)
      }
      
      
      
    }
    
  }
  
  
  
}






# Function to project a point onto a plane
project_to_plane <- function(point, plane_normal, plane_point) {
  # Extract the components of the point
  x <- point[1]
  y <- point[2]
  z <- point[3]
  
  # Extract the components of the plane normal vector
  a <- plane_normal[1]
  b <- plane_normal[2]
  c <- plane_normal[3]
  
  # Extract the components of a point on the plane
  x0 <- plane_point[1]
  y0 <- plane_point[2]
  z0 <- plane_point[3]
  
  # Calculate the projection
  t <- (a * (x0 - x) + b * (y0 - y) + c * (z0 - z)) / (a^2 + b^2 + c^2)
  
  # Calculate the projected point
  x_proj <- x + t * a
  y_proj <- y + t * b
  z_proj <- z + t * c
  
  return(c(x_proj, y_proj, z_proj))
}



# 외적을 계산하는 함수 정의
cross_product <- function(a, b) {
  c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}








#' @export plot.sphere
plot.sphere <- function(df, X=NULL, col="red", cex=0.01, opacity=FALSE, add=FALSE, main=NULL, main.cex=2, main.position="topleft"){
  
  library(rgl)
  
  
  if(FALSE){
    df=Y; col=col1; opacity=opacity; add=add
    cex=0.01; axis.arrange=FALSE
    main <- "This is the main title"
    main.cex <- 2
  }
  
  if(FALSE){
    df = sim.sphereGLM(n=100,s=1)$Y
    
    X=NULL; col="red"; cex=0.01; opacity=FALSE; add=FALSE; main=NULL; main.cex=2; main.position="topleft"
  }
  
  
  if(FALSE){
    devtools::load_all()
    
    set.seed(1)
    
    df <- sim.sphereGLM(n=100, mu = c(10,0,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.sphereGLM(n=100, mu = c(0,-10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.sphereGLM(n=100, mu = c(10,10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.sphereGLM(n=100, mu = c(-10,10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.sphereGLM(n=100, mu = c(-10,-10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.sphereGLM(n=100, mu = c(-10,10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    
    
    df <- sim.sphereGLM(n=100, mu = c(-10,10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.sphereGLM(n=100, mu = c(10,-10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.sphereGLM(n=100, mu = c(-10,-10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.sphereGLM(n=100, mu = c(10,10,10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    df <- sim.sphereGLM(n=100, mu = c(-10,10,-10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.sphereGLM(n=100, mu = c(10,-10,-10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.sphereGLM(n=100, mu = c(-10,-10,-10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    df <- sim.sphereGLM(n=100, mu = c(10,10,-10))$Y
    plot.sphere(df, col="red", opacity=FALSE, main="This is the main title")
    
    
    
    
    plot.sphere(sim.sphereGLM(n=100)$Y, col="blue", opacity=TRUE, add=TRUE)
    
    set.seed(1)
    plot.sphere(sim.sphereGLM(n=100, mu=c(1,0,0), s=10)$Y, opacity=T)
    set.seed(1)
    plot.sphere(sim.sphereGLM(n=100, mu = c(0,1,0), s=10)$Y, col="blue", add=TRUE)
  }
  
  if(FALSE){
    df <- sim.sphereGLM(n=100, s=10)$X
    col="red"; cex=0.01; opacity=TRUE; add=FALSE
  }
  
  
  
  
  
  if( ncol(df) != 3 ) stop("Check the number of columns")
  
  
  if(!is.null(X)){
    
    col1 <- "red"
    col2 <- "blue"
    col3 <- "purple"
    
    value_to_color <- function(value, min_val, max_val) {
      scaled <- (value - min_val) / (max_val - min_val)
      r <- 1 - scaled
      g <- 0
      b <- scaled
      rgb(r, g, b)
    }
    # blue: high; red: low
    
    col <- value_to_color(X, min(X), max(X))
    
  }
  
  
  
  
  # if(axis.arrange){
  #   x <- df[,1]
  #   y <- df[,3]
  #   z <- df[,2]
  # } else {
  x <- df[,1]
  y <- df[,2]
  z <- df[,3]
  # }
  
  
  if(!add){
    
    open3d()
    par3d(windowRect = c(20, 30, 800, 800))
    
    # getr3dDefaults()
    # rgl.par3d.names
    # par3d()$viewport
    
    
    # par3d(userMatrix = rotationMatrix(45*pi/180, 1, 0, 0),
    #       # userProjection = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=T),
    #       # windowRect=c(0,45,800,800)
    #       windowRect=c(0,45,500,500)
    # )
    
    
    # if(!is.null(main) | !is.null(submain)){
    #   
    #   bgplot3d({
    #     plot.new()
    #     if(!is.null(main)) title(main = main, line = 3)
    #     if(!is.null(submain)) mtext(side = 1, submain, line = 4)
    #     # use here any other way you fancy to write your title
    #   }, magnify = 0.5)
    #   
    # }
    
    if(opacity){
      spheres3d(0,0,0,lit=FALSE,color="white")
      spheres3d(0,0,0,radius=1.0,lit=FALSE,color="black",front="lines")
    } else {
      plot.sphere.grid(add=TRUE)
    }
    
    
  }
  
  
  if(!is.null(main)) legend3d(main.position, legend = main, adj=0.1, cex=main.cex, bty="n", border = FALSE)
  
  # legend3d("topright", legend = paste('Type', c('A', 'B', 'C')), pch = 16, col = rainbow(3), cex=1, inset=c(0.02))
  
  
  spheres3d(x,y,z,col=col,radius=cex)
  
  rgl::text3d(1.1,0,0, texts="(1,0,0)")
  rgl::text3d(0,1.1,0, texts="(0,1,0)")
  rgl::text3d(0,0,1.1, texts="(0,0,1)")
  
  rgl::text3d(-1.1,0,0, texts="(-1,0,0)")
  rgl::text3d(0,-1.1,0, texts="(0,-1,0)")
  rgl::text3d(0,0,-1.1, texts="(0,0,-1)")
  
  
  
  
  # viewpoint ----
  # Rotate viewpoint to the direction with the data points
  {
    angles <- calculate_view_angles(df)
    view3d(theta = angles['theta'], phi = angles['phi'], zoom = 0.75)
  }
  # {
  #   xbar <- colMeans(df) %>% {./norm(.,"2")}
  # 
  #   theta1 <- acos( (c(xbar[1],0,xbar[3]) %>% {./norm(.,"2")}) %*% c(0,0,1) ) * 180 / pi
  #   theta2 <- acos( -(c(xbar[1],0,xbar[3]) %>% {./norm(.,"2")}) %*% c(0,0,1) ) * 180 / pi - 180
  # 
  #   v1 <- (theta1 *pi/180) %>% {
  #     matrix(
  #     c(cos(.), 0, sin(.),
  #       0, 1, 0,
  #       -sin(.), 0, cos(.)),
  #     nrow = 3, byrow = TRUE
  #     )
  #   } %>% { . %*% c(0,0,1) }
  #   
  #   v2 <- (theta2 *pi/180) %>% {
  #     matrix(
  #       c(cos(.), 0, sin(.),
  #         0, 1, 0,
  #         -sin(.), 0, cos(.)),
  #       nrow = 3, byrow = TRUE
  #     )
  #   } %>% { . %*% c(0,0,1) }
  #   
  #   theta <- ifelse( sum(xbar*v1) > sum(xbar*v2), theta1, theta2 )
  #   
  # 
  #   rotation_matrix <- (theta *pi/180) %>% {
  #     matrix(
  #       c(cos(.), 0, sin(.),
  #         0, 1, 0,
  #         -sin(.), 0, cos(.)),
  #       nrow = 3, byrow = TRUE
  #     )
  #   }
  #   
  #   tmp <- rotation_matrix %*% cbind(c(0,0,1), c(0,0.1,-1) %>% {./norm(.,"2")})
  #   xbar2 <- tmp %*% solve(crossprod(tmp)) %*% t(tmp) %*% xbar
  # 
  #   phi1 <- acos( t(xbar) %*% (rotation_matrix %*% c(0,0,1)) ) * 180 / pi
  #   phi2 <- acos( -t(xbar) %*% (rotation_matrix %*% c(0,0,1)) ) * 180 / pi - 180
  #   
  #   
  #   v1 <- (phi1 *pi/180) %>% {
  #     matrix(
  #       c(1, 0, 0,
  #         0, cos(.), sin(.),
  #         0, -sin(.), cos(.)),
  #       nrow = 3, byrow = TRUE
  #     )
  #   } %>% { . %*% tmp[,1] }
  #   
  #   v2 <- (phi2 *pi/180) %>% {
  #     matrix(
  #       c(1, 0, 0,
  #         0, cos(.), sin(.),
  #         0, -sin(.), cos(.)),
  #       nrow = 3, byrow = TRUE
  #     )
  #   } %>% { . %*% tmp[,1] }
  #   
  #   phi <- ifelse( sum(xbar*v1) > sum(xbar*v2), phi1, phi2 )
  #   
  #   # phi <- atan( t(xbar) %*% (rotation_matrix %*% c(0,0,1)) ) * 180 / pi
  #   
  #   view3d(theta = 0, phi = 0, zoom = 0.75)
  #   view3d(theta = theta, phi = 0, zoom = 0.75)
  #   view3d(theta = theta, phi = phi, zoom = 0.75)
  # 
  # }
  # 
  
  
}







#' @export calculate_view_angles
calculate_view_angles <- function(df) {
  
  rotation_matrix_y <- function(angle) {
    angle_rad <- angle * pi / 180
    matrix(c(cos(angle_rad), 0, sin(angle_rad),
             0, 1, 0,
             -sin(angle_rad), 0, cos(angle_rad)), 
           nrow = 3, byrow = TRUE)
  }
  
  rotation_matrix_x <- function(angle) {
    angle_rad <- angle * pi / 180
    matrix(c(1, 0, 0,
             0, cos(angle_rad), sin(angle_rad),
             0, -sin(angle_rad), cos(angle_rad)), 
           nrow = 3, byrow = TRUE)
  }
  
  # 정규화된 평균 벡터 계산
  xbar <- colMeans(df)
  xbar <- xbar / sqrt(sum(xbar^2))
  
  # theta 계산
  v_xz <- c(xbar[1], 0, xbar[3])
  v_xz <- v_xz / sqrt(sum(v_xz^2))
  theta1 <- acos(sum(v_xz * c(0,0,1))) * 180 / pi
  theta2 <- acos(-sum(v_xz * c(0,0,1))) * 180 / pi - 180
  
  # 최적의 theta 선택
  v1 <- rotation_matrix_y(theta1) %*% c(0,0,1)
  v2 <- rotation_matrix_y(theta2) %*% c(0,0,1)
  theta <- ifelse(sum(xbar * v1) > sum(xbar * v2), theta1, theta2)
  
  # phi 계산을 위한 중간 단계
  R_theta <- rotation_matrix_y(theta)
  TT <- R_theta %*% c(0,0,1)
  
  phi1 <- acos(sum(xbar * (R_theta %*% c(0,0,1)))) * 180 / pi
  phi2 <- acos(-sum(xbar * (R_theta %*% c(0,0,1)))) * 180 / pi - 180
  
  # 최적의 phi 선택
  v1 <- rotation_matrix_x(phi1) %*% TT
  v2 <- rotation_matrix_x(phi2) %*% TT
  phi <- ifelse(sum(xbar * v1) > sum(xbar * v2), phi1, phi2)
  
  if(abs(theta) > 90){
    phi <- -phi
  }
  
  return(c(theta = theta, phi = phi))
}








#' @export plot.sphere.grid
plot.sphere.grid <- function (radius = 1, col.long = "red", col.lat = "blue", deggap = 15, longtype = "H", add = FALSE, radaxis = TRUE, radlab = "Radius"){
  if(FALSE){
    radius = 1; col.long = "red"; col.lat = "blue"; deggap = 15;
    longtype = "H"; add = FALSE; radaxis = TRUE; radlab = "Radius"
  }
  
  
  # if(FALSE){
  #   library(sphereplot)
  #   sphereplot::pointsphere()
  #   
  #   rgl.sphgrid()
  #   rgl.sphpoints(pointsphere(100,c(0,90),c(0,45),c(0.25,0.8)),deg=T)
  #   
  #   rgl.sphgrid(radaxis=F, radlab=F)
  #   rgl.sphpoints(40,50,0.5,deg=TRUE,col='red',cex=2)
  # }
  
  
  
  if (add == F) {
    open3d()
    
    par3d(userMatrix = rotationMatrix(-15*pi/180, 1, 0, 0),
          # userProjection = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=T),
          # windowRect=c(0,45,800,800)
          windowRect=c(0,45,800,800)
    )
    
  }
  for (lat in seq(-90, 90, by = deggap)) {
    if (lat == 0) {
      col.grid = "grey50"
    }
    else {
      col.grid = "grey"
    }
    plot3d(png.sph2coord(long = seq(0, 360, len = 100), lat = lat,
                         radius = radius, deg = T), col = col.grid, add = T,
           type = "l")
  }
  for (long in seq(0, 360 - deggap, by = deggap)) {
    if (long == 0) {
      col.grid = "grey50"
    }
    else {
      col.grid = "grey"
    }
    plot3d(png.sph2coord(long = long, lat = seq(-90, 90, len = 100),
                         radius = radius, deg = T), col = col.grid, add = T,
           type = "l")
  }
  if (longtype == "H") {
    scale = 15
  }
  if (longtype == "D") {
    scale = 1
  }
  # rgl.sphtext(long = 0, lat = seq(-90, 90, by = deggap), radius = radius,
  #             text = seq(-90, 90, by = deggap), deg = TRUE, col = col.lat)
  # rgl.sphtext(long = seq(0, 360 - deggap, by = deggap), lat = 0,
  #             radius = radius, text = seq(0, 360 - deggap, by = deggap)/scale,
  #             deg = TRUE, col = col.long)
  # if (radaxis) {
  #   radpretty = pretty(c(0, radius))
  #   radpretty = radpretty[radpretty <= radius]
  #   lines3d(c(0, 0), c(0, max(radpretty)), c(0, 0), col = "grey50")
  #   for (i in 1:length(radpretty)) {
  #     lines3d(c(0, 0), c(radpretty[i], radpretty[i]),
  #             c(0, 0, radius/50), col = "grey50")
  #     text3d(0, radpretty[i], radius/15, radpretty[i],
  #            col = "darkgreen")
  #   }
  #   text3d(0, radius/2, -radius/25, radlab)
  # }
  
}



png.coord2sph <- function(x,y,z){
  R=1
  lat = asin(z / R)
  long = atan2(y, x)
  cbind.data.frame(long=long, lat=lat)
}





png.sph2coord <- function (long, lat, radius = 1, deg = TRUE){
  # if(FALSE){
  #   x = R * cos(lat) * cos(long)
  #   y = R * cos(lat) * sin(long)
  #   z = R * sin(lat)
  #   cbind.data.frame(x=x,y=y,z=z)
  # }
  
  
  if (is.matrix(long) || is.data.frame(long)) {
    if (ncol(long) == 1) {
      long = long[, 1]
    }
    else if (ncol(long) == 2) {
      lat = long[, 2]
      long = long[, 1]
    }
    else if (ncol(long) == 3) {
      radius = long[, 3]
      lat = long[, 2]
      long = long[, 1]
    }
  }
  if (missing(long) | missing(lat)) {
    stop("Missing full spherical 3D input data.")
  }
  if (deg) {
    long = long * pi/180
    lat = lat * pi/180
  }
  return = cbind(x = radius*cos(lat)*cos(long),
                 y = radius*cos(lat)*sin(long),
                 z = radius*sin(lat))
}








