#' Mean Bias
#' @description
#' Calculate the MB of an array of estimates relative to the true value.
#'
#' @param rhohat an array of estimators of rho.
#' @param rho the true value of rho.
#'
#' @return the mean bias of rhohat array.
#' @seealso
#' [mrb]
#' [rmse]
#'
#' @export
#'
#' @examples
#' rho = 0.5
#' rhohat = 0.5 + rnorm(10)
#' mb(rhohat,rho)
mb = function(rhohat,rho){
  return(mean(rhohat-rho))
}
