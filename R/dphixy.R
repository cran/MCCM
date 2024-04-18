#' Scaled Bivariate Normal Density
#'
#' @importFrom mvtnorm dmvnorm
#' @description Bivariate normal density with mean 0 variance 1.
#'
#' @param x,y points value.
#' @param rho correlation coefficient.
#'
#' @return the density value.
#'
#' @export
#'
#' @examples
#' library(mvtnorm)
#' dmvnorm(c(1,-1),sigma = matrix(c(1,0.5,0.5,1),2,2))
#' dphixy(1,-1,0.5)
dphixy = function(x,y,rho){
  return(
    1/(2*pi*sqrt(1-rho^2))*exp(-(x^2+y^2-2*x*y*rho)/(2*(1-rho^2)))
  )
}
