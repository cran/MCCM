#' Polyserial Correlation
#' @description
#' Estimate the polyserial correlation coefficient.
#'
#' @param X a matrix(2*N) or dataframe contains two polyserial variable(Continuous variable first).
#' @param maxn the maximum iterations times.
#' @param e the maximum tolerance of convergence.
#'
#' @return \item{rho}{estimated value of polyserial correlation coefficient.}
#' @return \item{std}{standard deviation of rho.}
#' @return \item{iter}{times of iteration convergence.}
#' @return \item{Ex,Ey}{the support point of regression model.}
#'
#' @seealso
#' [esti_polychoric]
#'
#' @references
#' Zhang, P., Liu, B., & Pan, J. (2024). Iteratively Reweighted Least Squares Method for Estimating Polyserial and Polychoric Correlation Coefficients. Journal of Computational and Graphical Statistics, 33(1), 316–328. https://doi.org/10.1080/10618600.2023.2257251
#'
#' @examples
#' X = gen_polyseries(1000,0.5,-1:1)
#' result = esti_polyserial(X)
#' result
#' @export
esti_polyserial = function(X,maxn=100,e=1e-8){
  X = as.matrix(X)
  d = dim(X)
  if(min(d)<2){stop('You must have two variable.')}
  if(d[1]!=2){X = t(X)}
  y = X[1,]
  x = X[2,]
  xt = table(x)
  P = xt/sum(xt)
  s = length(P)
  Pi = c(1e-8,cumsum(xt)/sum(xt))
  Pi[s+1] = Pi[s+1]-1e-8
  xl = as.numeric(names(xt))
  ai = qnorm(Pi)


  Ex = (dnorm(qnorm(Pi[1:(length(Pi)-1)]))-dnorm(qnorm(Pi[2:length(Pi)])))/P
  Ey = rep(0,length(Ex))
  for (i in 1:length(xl)) {
    Ey[i] = mean(y[x==xl[i]])
  }

  rho = cor(Ex,Ey)
  iter =0;dif = 1;n=maxn;ep=e;rho0=rho
  while((iter<n)&(dif>ep)){
    sigma = (1 + rho^2*(ai[1:s]*dnorm(ai[1:s])-ai[2:(s+1)]*dnorm(ai[2:(s+1)]))/P-rho^2*(dnorm(ai[1:s])-dnorm(ai[2:(s+1)]))^2/P^2)/xt
    if(sum(is.na(sigma)>0)){
      warning("NA appears during the calculation of covariance matrix, and the result may be unstable.")
      break
    }
    #if(is.singular.matrix(diag(sigma))){sigma = runif(length(xt))}
    #rep(1,length(xt))/length(xt)}
    rho = sum(Ex*Ey/sigma)/sum(Ex^2/sigma)
    if(rho>=1){rho=1-1e-4}
    else if(rho<=-1){rho = -1+1e-4}
    dif = abs(rho - rho0)
    rho0 = rho
    iter = iter +1
  }
  varrho = 1/sum(Ex^2/sigma)
  if(iter>=maxn){
    warning('The number of iterations reaches the upper limit. Your dataset may be fragile.')
  }
  return(list(rho = rho,
              std = sqrt(varrho),
              iter = iter,
              Ex = Ex,
              Ey = Ey))
}
