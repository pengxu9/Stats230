#' Forward HMM simulation
#'
#' @param V initial distribution
#' @param P Transition Probability
#' @param E Emission Probability
#' @param S hidden states
#' @param obs simulated obersvations
#'
#' @return backward probabilities
#' @examples
#' S<-c("1","2")
#' P=matrix(c(0.98,0.05,0.02,0.95),nrow=2)
#' E=matrix(c(1/6,1/10,1/6,1/10,1/6,1/2,1/6,1/10,1/6,1/10,1/6,1/10,1/6,1/10),nrow=2)
#' V=matrix(c(1/2,1/2))
#' obs=c(1,2,3,4,5,6,7)
#'
#' backward_HMM(V,P,E,S,obs)
#' @export



backward_HMM <- function(V,P,E,S,obs) {
  rownames(P)<-S
  colnames(P)<-S
  rownames(E)<-S
  colnames(E)<-seq_len(dim(E)[2])
  P<-P
  E<-E
  N<-length(obs)
  M<-nrow(P)
  b <- matrix(1,nrow=N,ncol=M)
  colnames(b)<- S
  for (t in (N - 1):1) {
    for (i in 1:M) {
      b[t,i] <- sum(P[i,] *E[1:M,obs[t+1]] *b[t+1,])
    }
  }
  return(b)
}
