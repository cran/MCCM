#' Thresholds Estimation
#'
#' @importFrom mvtnorm rmvnorm
#' @description
#' Function to calculate thresholds from ordinal variables.
#'
#' @param X a ordinal series.
#'
#' @return the estimated value for thresholds.
#'
#' @export
#'
#' @examples
#' library(mvtnorm)
#' set.seed(1997)
#' R1 = gen_CCM(4)
#' n = 1000
#' indc = 3:4
#' thresholds = list(c(),c(),c(-1),c(1))
#' data1 = gen_mixed(n,R1,indc,thresholds=thresholds)$observed
#' est_thre(data1[,3])
#' est_thre(data1[,4])
est_thre = function(X){
  xt = table(X)
  P = xt/sum(xt)
  s = length(P)
  Pi = c(0,cumsum(xt)/sum(xt))
  ai = qnorm(Pi)
  ai[1] = -1e8
  ai[length(ai)] = 1e8
  return(ai)
}
