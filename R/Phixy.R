#' Scaled Bivariate Normal Approximation
#'
#' @importFrom mvtnorm pmvnorm
#' @description
#' Standard bivariate normal distribution approximated with Legendre polynomials.
#'
#' @param x,y P(X<=x,Y<=y).
#' @param rho correlation coefficient.
#' @param korder order of Legendre approximation.
#' @param app bool value TRUE for approximation, FALSE for integral.
#'
#' @return P(X<=x,Y<=y).
#'
#' @export
#'
#' @examples
#' library(mvtnorm)
#' pmvnorm(upper = c(1,-1),sigma = matrix(c(1,0.5,0.5,1),2,2))
#' Phixy(1,-1,0.5,2,app=TRUE)
#' Phixy(1,-1,0.5,app=TRUE)
Phixy = function(x,y,rho,korder=3,app=TRUE){
  if(app){
    if(korder==2){
      return(
        0.5*rho*(dphixy(x,y,(3-sqrt(3))/6*rho)+dphixy(x,y,(3+sqrt(3))/6*rho))+pnorm(x)*pnorm(y)
      )
    }else if(korder==3){
      return(
        1/18*rho*(5*dphixy(x,y,(1-sqrt(3/5))/2*rho)+5*dphixy(x,y,(1+sqrt(3/5))/2*rho)+8*dphixy(x,y,1/2*rho))+pnorm(x)*pnorm(y)
      )
    }
  }else{
    return(
      pmvnorm(upper = c(x,y),sigma = matrix(c(1,rho,rho,1),2,2))
    )
  }
}
