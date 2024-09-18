# tests/testthat/test-sphereGLM.R

library(testthat)
library(sphereGLM)

test_that("sphereGLM creates an object of class sphereGLM", {
  set.seed(1)
  simdata <- sphereGLM::sim.sphereGLM(n=100, p=1, q=3, mu=c(0,100,0), snr=50, s=10, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1])
  
  start <- proc.time()
  fit <- sphereGLM(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=1e-6, orthogonal=F)
  end <- proc.time()
  
  expect_s3_class(fit, "sphereGLM")
})


test_that("plot.sphereGLM executes without errors", {
  
  set.seed(1)
  simdata <- sphereGLM::sim.sphereGLM(n=100, p=1, q=3, mu=c(0,100,0), snr=50, s=10, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1])
  
  start <- proc.time()
  fit <- sphereGLM(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=1e-6, orthogonal=F)
  end <- proc.time()
  
  expect_silent(plot(fit))
})

test_that("summary.sphereGLM returns expected output", {
  
  
  set.seed(1)
  simdata <- sphereGLM::sim.sphereGLM(n=100, p=1, q=3, mu=c(0,100,0), snr=50, s=10, s0=0.0, type=c("vMF", "Proj", "ExpMap")[1])
  
  start <- proc.time()
  fit <- sphereGLM(X=simdata$X, Y=simdata$Y, eps=1e-6, maxit=100, lambda=1e-6, orthogonal=F)
  end <- proc.time()
  
  expect_output(summary(fit), "Summary of sphereGLM Object")
})
