#' @export eval.CovMat
eval.CovMat <- function(mat){
  out <- NULL
  
  out$fnorm <- norm(mat)
  out$gvar <- det(mat + diag(1e-16,nrow(mat),nrow(mat)))
  out$tvar <- sum(diag(mat))
  out$cnum <- eigen(mat)$values %>% {max(.)/min(.)}
  out$eig1 <- eigen(mat)$values %>% max()
  out$eig1 <- eigen(mat)$values %>% max()
  unlist(out)
}
