#' Positive Semidefinite Correlation Matrix
#' @description
#' Generate a positive semidefinite correlation coefficients matrix
#'
#' @param d the dimension of matrix.
#'
#' @return a correlation coefficients matrix.
#'
#' @export
#'
#' @examples
#' X = gen_CCM(4)
#' print(X)
gen_CCM = function(d){
  A = matrix(rnorm(d^2), d, d)
  A = A %*% t(A)
  sigma1 = sqrt(diag(A))
  A = A/(sigma1%*%t(sigma1))
  return(A)
}
