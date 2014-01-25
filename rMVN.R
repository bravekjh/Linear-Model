rMVN <- function(A.chol, b){
  n <- length(b)
  devs <- rnorm(n)
  return(backsolve(A.chol, backsolve(A.chol, b, transpose = TRUE) + devs))
}