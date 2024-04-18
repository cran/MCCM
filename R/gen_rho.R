#' Generate Specific Binormal Distribution
#' @description
#' Generate random number of binormal distribution with 0 mean unit variance and correlation coefficient rho.
#' @param n sample size.
#' @param rho correlation coefficient.
#'
#' @return Binormal random number with length n(in a 2*n matrix).
#'
#' @seealso
#' [gen_polyseries]
#' [gen_polychoric]
#' @export
#'
#' @examples
#' gen_rho(100,0.5)
gen_rho = function(n,rho){
  if(abs(rho)>=1){
    stop('rho must be in (-1,1)')
  }
  x = rnorm(2*n)
  x = matrix(x,2,n)
  R = matrix(c(1,rho,rho,1),2,2)
  R = chol(R)
  return(t(R)%*%x)
}
