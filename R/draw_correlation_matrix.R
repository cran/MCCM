#' Draw the Correlation Matrix
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom MASS ginv
#' @importFrom lavaan sem
#' @importFrom polycor polychor
#' @importFrom polycor polyserial
#' @import graphics
#'
#' @description
#' Estimate the MCCM from dataframe and draw it with scatter plot of matrices (SPLOM). With bivariate scatter plots below the diagonal, histograms on the diagonal, and the polychoric correlation coefficients with standard errors above the diagonal. Correlation ellipses are drawn in the same graph. The red lines below the diagonal are the LOESS smoothed lines, fitting a smooth curve between two variables.
#'
#' @param data1 a dataframe containing continuous or ordinal variable.
#' @param order_indx a vector to indicate the ordinal variables.
#' @param pair_est bool value, TRUE for pairwise estimation, FALSE for simultaneous estimation.
#' @param MLE bool value, TRUE for maximum likelihood estimation, FALSE for IRLS (pairwise) or IGMM (simultaneous) estimation.
#' @inheritParams est_mixedGMM
#'
#' @return the SPLOM plot.
#'
#' @seealso
#' [MCCM_est],
#' [summary_MCCM_est]
#'
#' @export
#'
#' @examples
#'
#' library(mvtnorm)
#' library(MASS)
#' library(polycor)
#' library(lavaan)
#' set.seed(1997)
#' n = 10000
#' rho12=0.3
#' rho13=0.4
#' rho14=0.5
#' rho23=0.6
#' rho24=0.7
#' rho34=0.8
#'
#' R = matrix(c(1,rho12,rho13,rho14,rho12,1,rho23,rho24,rho13,rho23,1,rho34,
#' rho14,rho24,rho34,1),4,4)
#' indc = c(3,4)
#' thresholds = list(c(),c(),0,0)
#' data1 = gen_mixed(n=n,R=R,indc=indc,thresholds=thresholds)
#' data2 = data.frame(data1$observed)
#'
#' # pairwise MLE estimation
#' draw_correlation_matrix(data2,indc,TRUE,TRUE)
#' # pairwise IRLS estimation
#' draw_correlation_matrix(data2,indc,TRUE,FALSE)
#' # simultaneous MLE estimation
#' draw_correlation_matrix(data2,indc,FALSE,TRUE)
#' # simultaneous IGMM estimation
#' draw_correlation_matrix(data2,indc,FALSE,FALSE)
draw_correlation_matrix = function(data1,order_indx,pair_est=FALSE,MLE=FALSE,R0=NULL,app=TRUE,korder=2,max_iter=1000,max_tol=1e-8,show_log = FALSE){
  diag_hist_density = function(x,ordinal,nameofx){
    usr = par(mar = c(0.5,0.5,0.5,0.5))
    hist(x, main=NULL,xlab = "", ylab = "", col = "lightblue", border = "black",prob=T,axes = FALSE)
    text(x=mean( par("usr")[1:2]),y=1.7*mean( par("usr")[3:4]),label=nameofx,font = 2,cex=1.2)
    lines(density(x), col = "black")
    if(ordinal){
      rug(x)
    }
  }

  lower_points_ellipse <- function(x, y, rhoxy, pch = 20, col.smooth = "red",
                                   span = 2/3, iter = 3, ...) {
    xm <- mean(x, na.rm = TRUE)
    ym <- mean(y, na.rm = TRUE)
    xs <- sd(x, na.rm = TRUE)
    ys <- sd(y, na.rm = TRUE)
    lines(stats::lowess(x, y, f = span,iter = iter), col = col.smooth)
    draw_ellipse(xm, ym, xs, ys, rhoxy, col.smooth = col.smooth, ...)
  }

  draw_ellipse <- function(x = 0, y = 0, xs = 1, ys = 1,
                           r = 0, col.smooth, add = TRUE, segments = 51, ...) {
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))

    if (abs(r) > 0)
      theta <- sign(r)/sqrt(2)
    else theta = 1/sqrt(2)
    shape <- diag(c(sqrt(1 + r), sqrt(1 - r))) %*% matrix(
      c(theta, theta, -theta, theta), ncol = 2, byrow = TRUE)
    ellipse <- unit.circle %*% shape
    ellipse[, 1] <- ellipse[, 1] * xs + x
    ellipse[, 2] <- ellipse[, 2] * ys + y

    points(x, y, pch = 19, col = col.smooth, cex = 1)
    lines(ellipse, ...)
  }
  n_variable = ncol(data1)

  estimation = MCCM_est(dataYX=data1,order_indx=order_indx,pair_est=pair_est,MLE=MLE,app=app,R0=R0,korder=korder,max_iter=max_iter,max_tol=max_tol,show_log=show_log)
  Rhatm = estimation$Rmatrix
  stdhatm = estimation$std_matrix

  op <- par(mfrow = c(n_variable, n_variable), mar = c(1, 1, 1, 1))
  on.exit(par(op))
  namesdata1 = names(data1)
  for (i in 1:n_variable) {
    for (j in 1:n_variable) {
      if(i==j){
        diag_hist_density(data1[,i],ordinal=(i %in% order_indx),namesdata1[i])
        usr <- par("usr")
        rect(usr[1], usr[3], usr[2], usr[4], col = NA, border = "black", lwd = 0)
      }else if(i>j){
        usr = par(mar = c(0.5,0.5,0.5,0.5))
        plot(data1[,i],data1[,j],pch=20,cex=0.5)
        lower_points_ellipse(data1[,i],data1[,j],Rhatm[i,j])
      }else{
        zstat = Rhatm[j,i]/stdhatm[j,i]
        pvalue = 2*(1-pnorm(abs(zstat)))
        if(pvalue<=0.001){
          starlabel = '***'
        }else if(pvalue<=0.01){
          starlabel = '**'
        }else if(pvalue<=0.05){
          starlabel = '*'
        }else if(pvalue<=0.1){
          starlabel = '.'
        }else{
          starlabel = ''
        }
        plot.new()
        usr = par("usr")
        text(x=mean( par("usr")[1:2]),y=mean( par("usr")[3:4]),label=paste(
          sprintf("%.4f",round(Rhatm[j,i],4)),'\n','(',
          sprintf("%.4f",round(stdhatm[j,i],4)),')','\n',starlabel),cex=1.2,font = 1)
        rect(usr[1], usr[3], usr[2], usr[4], col = NA, border = "black", lwd = 0)
      }
    }
  }
}
