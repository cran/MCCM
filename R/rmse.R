#' Root Mean Squared Error
#' @description
#' Calculate the RMSE of an array of estimates relative to the true value.
#'
#' @param rhohat an array of estimators of rho.
#' @param rho the true value of rho.
#'
#' @return the root mean squared error of rhohat array.
#' @seealso
#' [mb]
#' [mrb]
#'
#' @export
#'
#' @examples
#' rho = 0.5
#' rhohat = 0.5 + rnorm(10)
#' rmse(rhohat,rho)
rmse = function(rhohat,rho){
  return(sqrt(mean((rhohat-rho)^2)))
}
