#' HMM simulation
#'
#' @param V initial distribution
#' @param P Transition Probability
#' @param E Emission Probability
#' @param S hidden states
#' @param N number of simulations
#'
#' @return N of simulated X's and Y's
#' @examples
#' S<-c("1","2")
#' P=matrix(c(0.98,0.05,0.02,0.95),nrow=2)
#' E=matrix(c(1/6,1/10,1/6,1/10,1/6,1/2,1/6,1/10,1/6,1/10,1/6,1/10,1/6,1/10),nrow=2)
#' V=matrix(c(1/2,1/2))
#'
#' sim_HMM(V,P,E,S,N=100)
#' @export


sim_HMM=function(V,P,E,S,N){
  rownames(P)<-S
  colnames(P)<-S
  rownames(E)<-S
  colnames(E)<-seq_len(dim(E)[2])
  states   = c()
  obs = c()
  states   = c(states, sample(S,1,prob=V))
  for(i in 2:N)
  {
    s0  = sample(S, 1, prob=P[states[i-1],])
    states = c(states, s0)
  }
  for(i in 1:N)
  {
    obs0      = sample(seq_len(dim(E)[2]), 1, prob=E[states[i],])
    obs = c(obs, obs0)
  }
  return(list(states=states,obs=obs))
}
