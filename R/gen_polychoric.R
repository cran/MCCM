#' Generate Polychoric Sample
#' @description
#' Generate polychoric sample with hidden distribution: binormal with correlation coefficient rho.
#'
#' @param n sample size.
#' @param rho correlation coefficient.
#' @param a the cutoff points array.
#' @param b the cutoff points array.
#' @return Polychoric sample with size n(in a 2*n matrix).
#' @seealso
#' [gen_polyseries]
#' [gen_rho]
#' @export
#'
#' @examples
#' gen_polychoric(100,0.5,-1:1,1:2)
gen_polychoric = function(n,rho,a,b){
  x = gen_rho(n,rho)
  y1 = x[1,]
  y2 = x[2,]
  for (i in 1:length(x[2,])) {
    y1[i] = sum(y1[i]>a)+1
    y2[i] = sum(y2[i]>b)+1
  }
  x[1,] = y1
  x[2,] = y2
  return(x)
}
