#' Forward HMM simulation
#'
#' @param V initial distribution
#' @param P Transition Probability
#' @param E Emission Probability
#' @param S hidden states
#' @param obs simulated obersvations
#'
#' @return forwards probabilities
#' @examples
#' S<-c("1","2")
#' P=matrix(c(0.98,0.05,0.02,0.95),nrow=2)
#' E=matrix(c(1/6,1/10,1/6,1/10,1/6,1/2,1/6,1/10,1/6,1/10,1/6,1/10,1/6,1/10),nrow=2)
#' V=matrix(c(1/2,1/2))
#' obs=c(1,2,3,4,5,6,7)
#'
#' forward_HMM(V,P,E,S,obs)
#' @export


forward_HMM <- function(V,P,E,S,obs) {
  names(V)<-S
  rownames(P)<-S
  colnames(P)<-S
  rownames(E)<-S
  colnames(E)<-seq_len(dim(E)[2])
  V<-V
  P<-P
  E<-E
  N <-length(obs)
  M<-length(V)
  a <-matrix(0, nrow=N,ncol=M)
  colnames(a)<-S
  a[1, ] <- V[S]*E[S,obs[1]]
  for (t in 1:(N - 1)) {
    for (i in 1:M) {
      a[t+1,i]<- E[i, obs[t + 1]]*sum(a[t, 1:M]*P[1:M,i])
    }
  }
  return(a)
}
