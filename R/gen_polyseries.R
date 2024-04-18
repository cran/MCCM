#' Generate Polyseries Sample
#' @description
#' Generate polyseries sample with hidden distribution: binormal with correlation coefficient rho.
#'
#' @param n sample size.
#' @param rho correlation coefficient.
#' @param a the cutoff points array.
#'
#' @return Polyseries sample with size n(in a 2*n matrix).
#' @seealso
#' [gen_rho]
#' [gen_polychoric]
#' @export
#'
#' @examples
#' gen_polyseries(100,0.5,-1:1)
gen_polyseries = function(n,rho,a){
  x = gen_rho(n,rho)
  y = x[2,]
  for (i in 1:length(x[2,])) {
    y[i] = sum(y[i]>a)+1
  }
  x[2,] = y
  return(x)
}
