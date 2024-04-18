#' Mean Relative Bias
#' @description
#' Calculate the MRB of an array of estimates relative to the true value.
#'
#' @param rhohat an array of estimators of rho.
#' @param rho the true value of rho.
#'
#' @return the mean relative bias of rhohat array.
#' @seealso
#' [mb]
#' [rmse]
#'
#' @export
#'
#' @examples
#' rho = 0.5
#' rhohat = 0.5 + rnorm(10)
#' mrb(rhohat,rho)
mrb = function(rhohat,rho){
  return(mean((rhohat-rho)/rho))
}
