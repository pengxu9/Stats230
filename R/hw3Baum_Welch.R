#' Baum Welch HMM
#'
#' @param V Random initial distribution
#' @param P Random Transition Probability
#' @param E Random Emission Probability
#' @param S hidden states
#' @param obs simulated observations
#'
#' @return estimated True initial distribution,Transition Probability and Random Emission Probability
#' @examples
#' S<-c("1","2")
#' P=matrix(c(0.75,0.25,0.10,0.9),nrow=2)
#' E=matrix(c(1/6,1/6,1/6,1/6,1/6,1/6,1/6,1/6,1/6,1/6,1/6,1/6,1/6,1/6),nrow=2)
#' V=matrix(c(1/4,3/4))
#' obs=c(1,2,3,4,5,6,7)
#'
#' Baum_Welch(V,P,E,S,obs)
#' @export



Baum_Welch = function(V0,P0,E0,S,obs,n.iter,epilson){
  diff<- c()
  for (iter in 1:n.iter){
    T=length(obs)
    M=nrow(P0)
    K=ncol(E0)
    a<-forward_HMM(V0,P0,E0,S,obs)
    b<-backward_HMM(V0,P0,E0,S,obs)
    gamma<- a*b/rowSums(a*b,1)
    V<-gamma[1,]
    P<-P0*0
    g<-matrix(0, nrow=M,ncol=M)
    for (t in 2:T){
      for (i in 1:M){
        for (j in 1:M){
          g[i,j]<-b[t,j]*E[j,obs[t]]*P0[i,j]*a[t-1,i]
        }
      }
      P<-P+g/sum(g)
    }
    P<-P/rowSums(P)
    for (i in 1:M) {
      for (j in obs) {
        E[i,j]<-sum(gamma[,i]*(obs == j)) /sum(gamma[, i])
      }
    }
    d = sqrt(sum((P0-P)^2)) + sqrt(sum((E0-E)^2))
    diff=c(diff,d)
    if(d < epilson)
    {
      break
    }
  }
  return(list(V,P,E))
}

