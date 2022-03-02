#' Gradient descent/steepest ascent method
#'
#' @param X covariance matrix
#' @param Y response vector
#' @param alpha steep size
#' @param epilson stopping criteria
#' @param n.iter max iterartions
#'
#' @return beta parameters,corresponding asymptotic confidence intervals, and a vector of the log-likelihoods from Newton_Raphson
#' @examples
#' X<- matrix(x1,x2,..xn)
#' Y<-Y
#' gradient_descent(X,Y,1e-5,1e-8,1000)
#' @export

gradient_descent <- function(X,y,alpha,epilson,n.iter) {
  l_beta <- c()
  M<--diag(ncol(X))
  beta <- rep(0,ncol(X))
  for(i in 1:n.iter){
    theta<- X%*%beta
    p<-1/(1+exp(-theta))
    l_beta<-c(l_beta,t(y) %*% X %*% beta -sum(log(1+exp(theta))))
    delta_l_beta <- t(X)%*%(y-1/(1+exp(-theta)))
    beta<- beta-alpha*solve(M)%*%delta_l_beta
    diff<-sqrt(sum((delta_l_beta)^2))
    if(diff <= epilson){
      break
    }
  }
  theta<-X%*%beta
  p <- 1/(1+exp(-theta))
  W <- diag(as.vector(p*(1-p)))
  lowerCI<- beta[1,1]-1.96*sqrt(diag(solve(t(X)%*%W%*%X)))
  upperCI<- beta[1,1]+1.96*sqrt(diag(solve(t(X)%*%W%*%X)))
  return(list(beta,l_beta,lowerCI,upperCI))
}

