#' General Function to Estimate Mixed Correlation Coefficient Matrix
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom MASS ginv
#' @importFrom lavaan sem
#' @importFrom lavaan lavInspect
#' @importFrom polycor polychor
#' @importFrom polycor polyserial
#'
#' @description
#' Estimate the correlation matrix for dataframes containing continuous and ordinal variable, in pairs or simultaneously, using MLE, IRLS, or IGMM.
#'
#' @param dataYX a dataframe or matrix containing both continuous and ordinal variables.
#' @param order_indx a vector to indicate the ordinal variables.
#' @param pair_est bool value, TRUE for pairwise estimation, FALSE for simultaneous estimation.
#' @param MLE bool value, TRUE for maximum likelihood estimation, FALSE for IRLS (pairwise) or IGMM (simultaneous) estimation.
#' @inheritParams est_mixedGMM
#'
#' @return \item{Rmatrix}{Estimated mixed correlation coefficient matrix.}
#' @return \item{std_matrix}{Estimated standard deviation for each mixed correlation coefficient.}
#' @return \item{COV}{The covariance matrix for MCCM (simultaneous estimation only).}
#'
#' @seealso
#' [esti_polyserial],
#' [esti_polychoric],
#' [est_mixedGMM],
#' [summary_MCCM_est],
#' [draw_correlation_matrix]
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
#' out_pair_MLE = MCCM_est(dataYX=data2,order_indx=indc,pair_est=TRUE,MLE=TRUE)
#' # pairwise IRLS estimation
#' out_pair_IRLS = MCCM_est(dataYX=data2,order_indx=indc,pair_est=TRUE,MLE=FALSE)
#' # simultaneous MLE estimation
#' out_sim_MLE = MCCM_est(dataYX=data2,order_indx=indc,pair_est=FALSE,MLE=TRUE)
#' # simultaneous IGMM estimation
#' out_sim_IGMM = MCCM_est(dataYX=data2,order_indx=indc,pair_est=FALSE,MLE=FALSE)
#'
#' summary_MCCM_est(out_pair_MLE)
#' summary_MCCM_est(out_pair_IRLS)
#' summary_MCCM_est(out_sim_MLE)
#' summary_MCCM_est(out_sim_IGMM)
MCCM_est = function(dataYX,order_indx,pair_est=FALSE,MLE=FALSE,R0=NULL,app=TRUE,korder=2,max_iter=1000,max_tol=1e-8,show_log = FALSE){
  n_variable = ncol(dataYX)
  n_x = length(order_indx)
  n_y = n_variable - n_x
  nc_yy = n_y*(n_y-1)/2
  nc_xx = n_x*(n_x-1)/2
  nc_yx = n_y*n_x
  R0yy_indx = 1:nc_yy
  R0yx_indx = (nc_yy+1):(nc_yy+nc_yx)
  R0xx_indx = (nc_yy+nc_yx+1):(nc_yy+nc_yx+nc_xx)
  names_data = names(dataYX)

  if(pair_est){
    Rmatrix = diag(n_variable)
    colnames(Rmatrix) = names_data
    rownames(Rmatrix) = names_data
    std_matrix = Rmatrix
    if(MLE){
      for (i in 1:(n_variable-1)) {
        for (j in (i+1):n_variable) {
          x1 = dataYX[,i]
          y1 = dataYX[,j]
          if((i %in% order_indx)&(j %in% order_indx)){
            pair_mle = polychor(x1,y1,ML=TRUE,std.err = TRUE)
            Rmatrix[j,i] = pair_mle$rho
            std_matrix[j,i] = sqrt(pair_mle$var[1,1])
          }else if((i %in% order_indx)&(!(j %in% order_indx))){
            pair_mle = polyserial(y1,x1,ML = TRUE,std.err = TRUE)
            Rmatrix[j,i] = pair_mle$rho
            std_matrix[j,i] = sqrt(pair_mle$var[1,1])
          }else if((!(i %in% order_indx))&(j %in% order_indx)){
            pair_mle = polyserial(x1,y1,ML = TRUE,std.err = TRUE)
            Rmatrix[j,i] = pair_mle$rho
            std_matrix[j,i] = sqrt(pair_mle$var[1,1])
          }else{
            pair_mle = cor.test(x1,y1)
            Rmatrix[j,i] = pair_mle$estimate
            std_matrix[j,i] = pair_mle$estimate/pair_mle$statistic
          }
        }
      }
      return(list(
        names_data = names_data,
        pair_est = pair_est,
        MLE = MLE,
        Rmatrix = Rmatrix,
        std_matrix = std_matrix
      ))
    }else{
      for (i in 1:(n_variable-1)) {
        for (j in (i+1):n_variable) {
          if((i %in% order_indx)&(j %in% order_indx)){
            X = t(as.matrix(dataYX[,c(i,j)]))
            pair_IRLS = esti_polychoric(X)
            Rmatrix[j,i] = pair_IRLS$rho
            std_matrix[j,i] = pair_IRLS$std
          }else if((i %in% order_indx)&(!(j %in% order_indx))){
            X = t(as.matrix(dataYX[,c(j,i)]))
            pair_IRLS = esti_polyserial(X)
            Rmatrix[j,i] = pair_IRLS$rho
            std_matrix[j,i] = pair_IRLS$std
          }else if((!(i %in% order_indx))&(j %in% order_indx)){
            X = t(as.matrix(dataYX[,c(i,j)]))
            pair_IRLS = esti_polyserial(X)
            Rmatrix[j,i] = pair_IRLS$rho
            std_matrix[j,i] = pair_IRLS$std
          }else{
            x1 = dataYX[,i]
            y1 = dataYX[,j]
            pair_mle = cor.test(x1,y1)
            Rmatrix[j,i] = pair_mle$estimate
            std_matrix[j,i] = pair_mle$estimate/pair_mle$statistic
          }
        }
      }
      return(list(
        names_data = names_data,
        pair_est = pair_est,
        MLE = MLE,
        Rmatrix = Rmatrix,
        std_matrix = std_matrix
      ))
    }
  }else{
    if(MLE){
      names_matrix = outer(names_data,names_data,function(x, y) paste(paste(x, y, sep = " ~~ "),'\n',sep=''))
      model1 = t(names_matrix)[lower.tri(names_matrix)]
      model1 = paste(model1,collapse = " ")
      m1 = sem(model1,data=dataYX,ordered = names_data[order_indx])
      est1 = lavInspect(m1,'est')
      Rmatrix = as.matrix(est1$theta)
      Rhat = Rmatrix[lower.tri(Rmatrix)]

      COV = (m1@vcov$vcov)[1:length(Rhat),1:length(Rhat)]
      names_matrix1 = t(outer(names_data,names_data,function(x, y) paste(x, y, sep = "-")))
      names_cov = names_matrix1[lower.tri(names_matrix1)]
      colnames(COV) = names_cov
      rownames(COV) = names_cov
      names(Rhat) = names_cov

      stdhat = sqrt(diag(COV))
      stdhatm = diag(n_variable)

      indm = matrix(FALSE,n_variable,n_variable)
      indm[1:n_y,1:n_y] = lower.tri(indm[1:n_y,1:n_y])
      stdhatm[indm] = stdhat[R0yy_indx]

      stdhatm[(n_y+1):(n_y+n_x),1:n_y] = stdhat[R0yx_indx]

      indm = matrix(FALSE,n_variable,n_variable)
      indm[(n_y+1):(n_y+n_x),(n_y+1):(n_y+n_x)] = lower.tri(indm[(n_y+1):(n_y+n_x),(n_y+1):(n_y+n_x)])
      stdhatm[indm] = stdhat[R0xx_indx]
      colnames(stdhatm) = names_data
      rownames(stdhatm) = names_data

      return(list(
        names_data = names_data,
        names_cov = names_cov,
        pair_est = pair_est,
        MLE = MLE,
        Rmatrix = Rmatrix,
        Rhat = Rhat,
        std_matrix = stdhatm,
        COV = COV
      ))
    }else{
      out1 = est_mixedGMM(dataYX=dataYX,order_indx=order_indx,R0=R0,app=app,korder=korder,max_iter=max_iter,max_tol=max_tol,show_log=show_log)
      out1[['pair_est']] = FALSE
      out1[['MLE']] = FALSE
      Rhat = out1$Rhat
      stdhat = sqrt(diag(out1$COV))
      Rhatm = diag(n_variable)
      stdhatm = Rhatm

      indm = matrix(FALSE,n_variable,n_variable)
      indm[1:n_y,1:n_y] = lower.tri(indm[1:n_y,1:n_y])
      Rhatm[indm] = cor(dataYX[,setdiff(1:n_variable,order_indx)])[lower.tri(indm[1:n_y,1:n_y])] + rnorm(sum(indm),0,0.01)
      stdhatm[indm] = stdhat[R0yy_indx]

      Rhatm[(n_y+1):(n_y+n_x),1:n_y] = Rhat[R0yx_indx]
      stdhatm[(n_y+1):(n_y+n_x),1:n_y] = stdhat[R0yx_indx]

      indm = matrix(FALSE,n_variable,n_variable)
      indm[(n_y+1):(n_y+n_x),(n_y+1):(n_y+n_x)] = lower.tri(indm[(n_y+1):(n_y+n_x),(n_y+1):(n_y+n_x)])
      Rhatm[indm] = Rhat[R0xx_indx]
      stdhatm[indm] = stdhat[R0xx_indx]

      colnames(Rhatm) = names_data
      rownames(Rhatm) = names_data

      colnames(stdhatm) = names_data
      rownames(stdhatm) = names_data

      out1[['Rmatrix']] = Rhatm
      out1[['std_matrix']] = stdhatm

      return(out1)
    }
  }
}
