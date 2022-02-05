#' Sampling from MVN using Cholesky decoposition
#'
#' @param N, number of samples
#' @param mu mean vector
#' @param sigma covariance matrix
#' @return N realizations from MVN(mu,sigma)
#' @examples
#' n<-10
#' N<-100
#' mu<- mu<-matrix(rnorm(n),nrow=n,ncol=1)
#' sigma<- matrix(rWishart(1,n,diag(n)),nrow=n,ncol=n)
#' mvn_cholesky(N,mu,sigma)
#' @export

mvn_cholesky<-function(N,mu,sigma){
  n<-ncol(sigma)
  L<-t(chol(sigma))
  Z<-matrix(rnorm(N*n),nrow=n,ncol=N)
  X<-matrix(rep(mu,N),nrow=n,ncol=N)+L%*%Z
  return(X)
}
