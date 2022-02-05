#' QLS estimation using SVD decomposition
#'
#' @param X covariance matrix
#' @param Y response vector
#'
#' @return beta parameters from QLS estimation using QR decomposition
#' @examples
#' X<- matrix(x1,x2,..xn)
#' Y<-Y
#' QLS_SVD(X,Y)
#' @export
QLS_SVD<-function(X,Y){
  X_SVD<-svd(X)
  D<-diag(X_SVD$d)
  invD<-solve(D)
  U<-X_SVD$u
  V<-X_SVD$v
  beta<-V%*%invD%*%t(U)%*%Y
  #b2<-matrix(beta[,1],nrow=1)
  #colnames(b2)<- c(n)
  n<-paste0('x', 1:(ncol(X)),sep = "")
  return(setNames(beta[,1],c(n)))
}
