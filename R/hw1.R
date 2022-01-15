#' Matrices and Vector Mulitiplication
#'
#' @param A a square matrix
#' @param B a square matrix
#' @param x a vector
#' @param matrices, if matrices is true, matrices A and B will be multiplied first, else, Matrix B and vector x will be multiplied first
#' @return product A*B*x
#' @examples
#' A<-matrix(c(1,5,3,8), ncol=2, nrow=2)
#' B<- matrix(c(1,7,4,8), ncol=2, nrow=2)
#' x<-matrix(c(1,5), ncol=2, nrow=1)
#' matrices_vector_mult (A,B,x,matrices = T), matrices_vector_mult (A,B,x,matrices = F)

matrice_vector_mult <- function(A,B,x, matrices=T){
  if (matrices){
    c <- A %*% B
    return(c %*% x)
    }
  else {
    c<-B %*% x
    return (a %*% c)
  }

}

