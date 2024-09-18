
sim01 <- function(){
  # Simulation for controlling type I error rate
  
  II=1
  
  
  
  n.thread <- 8
  n.save <- 1 # 50
  
  # detach("package:sphereGLM",unload=TRUE)
  library(sphereGLM)
  library(dplyr)
  
  
  {
    n.seq <- c(20,30,50,100)
    p.seq <- c(1,2)
    q.seq <- c(3)
    m.seq <- c(20, 50)
    s.seq <- c(0, 5, 10, 20, 50)
    seed.seq <- c(1:100)
    
    params <- list(n.seq=n.seq,
                   p.seq=p.seq,
                   q.seq=q.seq,
                   m.seq=m.seq,
                   s.seq=s.seq,
                   seed.seq=seed.seq)
    
    grid <- expand.grid(params)
  }
  
  
  GRID <- grid[,-which(colnames(grid) == "seed.seq")] %>% 
    {.[!duplicated(.),]}
  
  
  
  
  # GRID %>% head(20)
  # idx <- 5
  # simdata <- sphereGLM::sim.sphereGLM(n=GRID[idx,]$n.seq, p=GRID[idx,]$p.seq, q=GRID[idx,]$q.seq, mu=c(0,0,GRID[idx,]$mu.seq), snr=NULL, s=GRID[idx,]$s.seq, s0=0, type="vMF", seed.UDV=1, seed.E=idx)
  # plot.sphere(simdata$Y)
  # 
  # plot.sphereGLM(sphereGLM(X=simdata$X, Y=simdata$Y))
  
  
  
  
  
  AA.list <- png.utils::png.sim.split(nrow(grid), n.thread)
  AA <- AA.list[[II]]
  
  out.list <- NULL
  for( ii in AA ){
    
    print(paste0(ii, " / ", max(AA)))
    
    # Simulation Parameters
    {
      n=grid[ii,]$n.seq
      p=grid[ii,]$p.seq
      q=grid[ii,]$q.seq
      s=grid[ii,]$s.seq
      mu=grid[ii,]$m.seq
      seed=grid[ii,]$seed.seq
      
      mu=c(0,0,mu)
      s0=0
      type="vMF"
      snr=NULL
    }
    
    
    
    devtools::load_all()
    
    
    
    result <- lapply(c(20, 50, 100), function(n){
      simdata.list <- fit.list <- NULL
    for( jj in 1:100 ){
      if( jj %% 10 == 0 ) print(jj)
      simdata <- sphereGLM::sim.sphereGLM(n=n, p=p, q=q, mu=mu, s=s, s0=s0, type="vMF", seed.UDV=seed, seed.E=100*seed+jj)
      fit <- with(simdata, sphereGLM(X=X, Y=Y, use.nlm=F))
      simdata.list[[jj]] <- simdata
      fit.list[[jj]] <- fit
    }
    
    pv.list <- sapply(fit.list, function(FIT) summary.sphereGLM(FIT, scale=FALSE)$p.value[2])
    mean(pv.list < 0.05)
    hist(pv.list)
    return(pv.list)
    })
    
    
    result2 <- lapply(c(20, 50, 100), function(n){
      simdata.list <- fit.list <- NULL
      for( jj in 1:100 ){
        if( jj %% 10 == 0 ) print(jj)
        simdata <- sim.sphereGLM(n=n, p=p, q=q, mu=mu, s=s, s0=s0, type="vMF", seed.UDV=seed, seed.E=100*seed+jj, qr.U=TRUE)
        fit <- with(simdata, sphereGLM(X=X, Y=Y, use.nlm=F))
        simdata.list[[jj]] <- simdata
        fit.list[[jj]] <- fit
      }
      
      pv.list <- sapply(fit.list, function(FIT) summary.sphereGLM(FIT, scale=FALSE)$p.value[2])
      mean(pv.list < 0.05)
      hist(pv.list)
      return(pv.list)
    })
    
    
    lapply( result, function(pv.list) mean(pv.list < 0.05) )
    lapply( result2, function(pv.list) mean(pv.list < 0.05) )
    
    
    #
    
    # fit.list0 <- fit.list;  simdata.list0 <- simdata.list
    # fit.list1 <- fit.list;  simdata.list1 <- simdata.list
    # fit.list10 <- fit.list;  simdata.list10 <- simdata.list
    # fit.list20 <- fit.list;  simdata.list20 <- simdata.list
    # fit.list50 <- fit.list;  simdata.list50 <- simdata.list
    # 
    # mean(sapply(fit.list, function(FIT) summary.sphereGLM(FIT, scale=FALSE)$p.value[2]) < 0.05)
    # 
    # sapply(fit.list, function(FIT) summary.sphereGLM(FIT, scale=FALSE)$p.value[2]) %>% hist
    # 
    # sapply(fit.list, function(FIT) summary.sphereGLM(FIT, scale=FALSE)$p.value[2]) %>% hist
    
    #
    #
    #
    
    out.list[[ii]] <- list( simdata.list=simdata.list,
                            fit.list=fit.list,
                            params=params,
                            grid=grid,
                            AA=AA)
    
    
    check.start = ( ( ii-AA[1]+1 ) / length(AA)*100 ) == 1
    check.10xNperc = ( ( ii-AA[1]+1 ) / length(AA)*100 ) %% 10 == 0
    if( ii == 1 | check.10xNperc ){
      save(out.list, file=paste0("sim01-1_II=",II,"_ii=",ii,".RData"))
    }
    
  }
  
  
  
  
  
}