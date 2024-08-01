#' @export png.sphere.grid
png.sphere.grid <- function (radius = 1, col.long = "red", col.lat = "blue", deggap = 15, longtype = "H", add = FALSE, radaxis = TRUE, radlab = "Radius"){
  if(FALSE){
    radius = 1; col.long = "red"; col.lat = "blue"; deggap = 15;
    longtype = "H"; add = FALSE; radaxis = TRUE; radlab = "Radius"
  }
  
  
  if(FALSE){
    library(sphereplot)
    sphereplot::pointsphere()
    
    rgl.sphgrid()
    rgl.sphpoints(pointsphere(100,c(0,90),c(0,45),c(0.25,0.8)),deg=T)
    
    rgl.sphgrid(radaxis=F, radlab=F)
    rgl.sphpoints(40,50,0.5,deg=TRUE,col='red',cex=2)
  }
  
  
  
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


#' @export png.sphere
png.sphere <- function(df, col="red", cex=0.01, opacity=TRUE, add=FALSE, axis.arrange=FALSE){
  
  
  library(rgl)
  
  if(FALSE){
    set.seed(1)
    png.sphere(sim.sphereGLM(n=2)$Y, opacity=F)
    set.seed(1)
    png.sphere(sim.sphereGLM(n=100, s=10)$Y, opacity=T)
  }
  
  if(FALSE){
    df <- sim.sphereGLM(n=100, s=10)$X
    col="red"; cex=0.01; opacity=TRUE; add=FALSE
  }
  
  
  
  
  if( ncol(df) != 3 ) stop("Check the number of columns")
  
  x <- df[,1]
  y <- df[,2]
  z <- df[,3]
  
  
  if(opacity){
    if(add){
      
      # par3d(userMatrix = rotationMatrix(75*pi/180, 1, 0, 0)
      #       # userProjection = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=T),
      #       # windowRect=c(0,45,800,800)
      #       # windowRect=c(0,45,800,800)
      # )
      
      if(axis.arrange){
        spheres3d(x,z,y,col=col,radius=cex)
      } else {
        spheres3d(x,y,z,col=col,radius=cex)
      }
      
      
    } else {
      
      
      
      open3d()
      # getr3dDefaults()
      # rgl.par3d.names
      # par3d()$viewport
      par3d(userMatrix = rotationMatrix(45*pi/180, 1, 0, 0),
            # userProjection = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=T),
            # windowRect=c(0,45,800,800)
            windowRect=c(0,45,800,800)
      )
      spheres3d(0,0,0,lit=FALSE,color="white")
      spheres3d(0,0,0,radius=1.0,lit=FALSE,color="black",front="lines")
      
      if(axis.arrange){
        spheres3d(y,z,x,col=col,radius=cex)
      } else {
        spheres3d(x,y,z,col=col,radius=cex)
      }
      
      
      
      
      
      # axis.arrange = FALSE
      # {
      #   
      #   open3d()
      #   # getr3dDefaults()
      #   # rgl.par3d.names
      #   # par3d()$viewport
      #   par3d(userMatrix = rotationMatrix(45*pi/180, 1, 0, 0),
      #         # userProjection = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,byrow=T),
      #         # windowRect=c(0,45,800,800)
      #         windowRect=c(0,45,800,800)
      #   )
      #   spheres3d(0,0,0,lit=FALSE,color="white")
      #   spheres3d(0,0,0,radius=1.0,lit=FALSE,color="black",front="lines")
      #   spheres3d(x,z,y,col=col,radius=cex)
      #   
      # }
      
    }
  } else {
    
    
    if(axis.arrange){
      df <- df[,c(2,3,1)]
      df[,3] <- df[,3]
      x <- df[,1]
      y <- df[,2]
      z <- df[,3]
      
      
      png.sphere.grid(add=F)
      points3d(x,y,z,col=col,radius=cex)
    } else {
      png.sphere.grid(add=F)
      spheres3d(x,y,z,col=col,radius=cex)
    }
    
  }
  
  rgl::text3d(1.1,0,0, texts="(1,0,0)")
  rgl::text3d(0,1.1,0, texts="(0,1,0)")
  rgl::text3d(0,0,1.1, texts="(0,0,1)")
  
  rgl::text3d(-1.1,0,0, texts="(-1,0,0)")
  rgl::text3d(0,-1.1,0, texts="(0,-1,0)")
  rgl::text3d(0,0,-1.1, texts="(0,0,-1)")
  
}




#' @export png.coord2sph
png.coord2sph <- function(x,y,z){
  R=1
  lat = asin(z / R)
  long = atan2(y, x)
  cbind.data.frame(long=long, lat=lat)
}





#' @export png.sph2coord
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











#' @export plot.sphereGLM
plot.sphereGLM <- function(fit, opacity=FALSE, plot.mu=FALSE){
  
  Y <- fit$Y
  
  if(FALSE){
    
    library(png.Directional)
    library(sphereGLM)
    library(dplyr)
    
    
    # Figure 1
    {
      set.seed(2)
      simdata <- sim.sphereGLM(n=150, mu=c(1/sqrt(3),1/sqrt(3),1/sqrt(3)) %>% {./norm(.,"2")}, snr=100, r=2, s=3, s0 = 0.0001, type="Proj")
      X <- simdata$U;  Y <- simdata$X
      fit <- sphereGLM(X, Y)
      
      plot.sphereGLM(fit, plot.mu=TRUE)
      
      fit
    }
    
    
    # Figure 2
    {
      set.seed(2)
      simdata <- sim.sphereGLM(n=150, mu=c(1/sqrt(3),1/sqrt(3),1/sqrt(3)) %>% {./norm(.,"2")}, snr=100, r=1, s=10, s0 = 0.0001, type="Proj")
      X <- simdata$U;  Y <- simdata$X
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
  
  

  mu <- fit$beta[1,,drop=T]
  beta <- fit$beta[-1,,drop=FALSE]
  beta <- apply(beta, 1, function(x) x / norm(mu, "2")) %>% t()
  mu <- mu / norm(mu, "2")
  
  # For better visualization
  beta <- apply(beta, 1, function(x) x / norm(x, "2") ) %>% t()
  
  
  # rgl::planes3d(0,0,1,1, alpha=0.2)
  MU <- tcrossprod(rep(1,nrow(beta)), mu)
  
  
  
  if(plot.mu){
    
    png.sphere(Y, opacity=opacity)
    rgl::arrow3d(c(0,0,0), mu, type=c("lines", "rotation")[1], col="purple", s=0.1,barblen=0.05, width=0.5, add=TRUE)
    
    if(nrow(beta)>1){
      normal_vector <- project_to_plane(c(0,0,0), cross_product(beta[1,], beta[2,]), mu)
      rgl::planes3d(a=normal_vector, d=-1, alpha=0.1, add=TRUE)
      # rgl::abclines3d(mu[1],mu[2],mu[3], beta, alpha=0.5, col="blue", lwd=2)
    }
    
    for( k in 1:nrow(MU) ){
      rgl::arrow3d(MU[k,]-beta[k,], MU[k,]+beta[k,], type=c("lines", "rotation")[1], col="blue", s=0.1, barblen=0.05, width=0.5, add=TRUE)
    }
    
  } else {
    
    png.sphere(Y, opacity=opacity)
    
    if(nrow(beta)>1){
      normal_vector <- project_to_plane(c(0,0,0), cross_product(beta[1,], beta[2,]), mu)
      rgl::planes3d(a=normal_vector, d=-1, alpha=0.1, add=TRUE)
      # rgl::abclines3d(mu[1],mu[2],mu[3], beta, alpha=0.5, col="blue", lwd=2)
    }
    
    for( k in 1:nrow(MU) ){
      rgl::arrow3d(MU[k,]-beta[k,], MU[k,]+beta[k,], type=c("lines", "rotation")[1], col="blue", s=0.1, barblen=0.05, width=0.5, add=TRUE)
    }
    
  }
  
  
  
  
}
