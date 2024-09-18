#' @method print sphereGLM
#' @export 
print.sphereGLM <- function(fit, all=FALSE){
  
  
  if(all){
    
    print(fit)
    
  } else{
    
    print(
      list(beta = fit$beta,
         beta.scaled = fit$beta2,
         loglik = fit$loglik.list,
         params = c(maxit = fit$params$maxit,
                    eps = fit$params$eps,
                    orthogonal = fit$params$orthogonal,
                    use.nlm = fit$params$use.nlm)
      )
    )
  }
  
  
}

