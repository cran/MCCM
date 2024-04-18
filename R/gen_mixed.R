#' Continuous and Ordinal Simulated Data
#'
#' @importFrom mvtnorm rmvnorm
#'
#' @description
#' Generate multi-normal sample and segment it into ordinal.
#'
#' @param n the sample size.
#' @param R the correlation coefficient matrix.
#' @param indc vector to indicate whether variables are continuous or categorical.
#' @param thresholds list contains thresholds for ordinal variables
#'
#' @return \item{latent}{the original normal data.}
#' @return \item{observed}{the observed ordinal data.}
#'
#' @export
#'
#' @examples
#' library(mvtnorm)
#' set.seed(1997)
#' R1 = gen_CCM(6)
#' n = 1000
#' indc = 4:6
#' thresholds = list(
#'   c(),
#'   c(),
#'   c(),
#'   c(0),
#'   c(-0.5,0),
#'   c(0,0.5)
#' )
#' data1 = gen_mixed(n,R1,indc,thresholds)$observed
#' data1 = data.frame(data1)
#' table(data1$X4,data1$X5)
#' table(data1$X5,data1$X6)
gen_mixed = function(n,R,indc,thresholds){
  data1 = rmvnorm(n,sigma = R)
  data2 = data1
  for (i in indc) {
    a = c(-Inf,thresholds[[i]],Inf)
    data2[,i] = as.numeric(cut(data1[,i],a))
  }
  return(list('latent'=data1,
              'observed'=data2))
}
