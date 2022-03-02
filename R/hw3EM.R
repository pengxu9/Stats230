#' ABO Blood inference by EM Algorithm
#'
#' @param n different blood type counts
#' @param p nitial values of the allele frequency probabilities
#' @param epilson stopping criteria
#'
#' @return beta parameters,corresponding asymptotic confidence intervals, and a vector of the log-likelihoods from Newton_Raphson
#' @examples
#' n<-c(2,4,50,1)
#' p<-c(1/3,1/3,1/3)
#' EM(n,p,epilson=1e-8)
#' @export


EM<-function(n,p,epilson){
  p<-p
  pA<-p[1]
  pB<-p[2]
  pO<-p[3]
  nA<-n[1]
  nAB<-n[2]
  nB<-n[3]
  nO<-n[4]
  diff <-Inf
  m<-sum(n)
  while (diff > epilson){
    p0<-p
    mAA<-nA*pA^2/(pA^2+2*pA*pO)
    mAO<-nA*2*pA*pO/(pA^2+2*pA*pO)
    mBB<-nB*pB^2/(pB^2+2*pB*pO)
    mBO<-nB*2*pB*pO/(pB^2+2*pB*pO)
    mAB<-nAB
    mOO<-nO
    pA<-(2*mAA+mAO+mAB)/(2*m)
    pB<-(2*mBB+mBO+mAB)/(2*m)
    pO<-(2*n[4]+mAO+mBO)/(2*m)
    diff <- sqrt(sum(p-p0)^2)
    pk <- (pA^2+2*pA*pO)^(nA)*(pB^2+2*pB*pO)^(nB)*(pO^2)^(nO)
  }
  return(list(pA=pA,pB=pB,pO=pO,pk=pk))
}

