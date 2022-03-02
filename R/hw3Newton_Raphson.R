#' Newton-Raphson/Fisher scoring/IRLS method
#'
#' @param X covariance matrix
#' @param Y response vector
#' @param epilson stopping criteria
#' @param n.iter max iterations
#'
#' @return beta parameters,corresponding asymptotic confidence intervals, and a vector of the log-likelihoods from Newton_Raphson
#' @examples
#' X<- matrix(x1,x2,..xn)
#' Y<-Y
#' Newton_Raphson(X,Y,1e-8,1000)
#' @export

Newton_Raphson <- function(X,y,epilson,n.iter) {
  l_beta <- c()
  beta <- rep(0,ncol(X))
  diff <- Inf
  for (i in 1:n.iter){
    theta<- X%*%beta
    p<-1/(1+exp(-theta))
    l_beta<-c(l_beta,t(y) %*% X %*% beta -sum(log(1+exp(theta))))
    delta_l_beta <- t(X)%*%(y-1/(1+exp(-theta)))
    W<-diag(as.vector(p*(1-p)))
    d2l_beta <- t(X)%*%W%*%X
    beta<-beta+solve(d2l_beta)%*%(t(X)%*%(y-p))
    diff<- sqrt(sum((solve(d2l_beta)%*%(t(X)%*%(y-p)))^2))
    if (diff <= epilson){
      break
    }
  }
  theta<- X%*%beta
  p<-1/(1+exp(-theta))
  W <- diag(as.vector(p*(1-p)))
  lowerCI<-beta-1.96*sqrt(diag(solve(t(X)%*%W%*%X)))
  upperCI<-beta+1.96*sqrt(diag(solve(t(X)%*%W%*%X)))
  return(list(beta,l_beta,lowerCI,upperCI))
}
