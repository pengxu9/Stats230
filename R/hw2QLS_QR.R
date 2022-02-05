#' QLS estimation using QR decomposition
#'
#' @param X covariance matrix
#' @param Y response vector
#'
#' @return beta parameters from QLS estimation using QR decomposition
#' @examples
#' X<- matrix(x1,x2,..xn)
#' Y<-Y
#' QLS_QR(X,Y)
#' @export

QLS_QR<-function(X,Y){
  QR<-qr(X)
  Q<-qr.Q(QR)
  R<-qr.R(QR)
  beta<-solve.qr(QR,Y)
  #beta<-solve(R)%*%t(Q)%*%Y
  return(beta)
}
