#' Estimating Mixed Correlation Matrix by IGMM
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom MASS ginv
#' @importFrom lavaan sem
#' @description
#' An accelerated function to estimate a mixed correlation coefficient matrix, as well as its covariance matrix, for dataframes containing continuous and ordinal variable.
#'
#' @param dataYX a dataframe or matrix containing both continuous and ordinal variables.
#' @param order_indx a vector to indicate the ordinal variables.
#' @param R0 the initial value for correlation vector, default Pearson correlation matrix.
#' @param app bool value for approximation, TRUE for Legendre approximation, FALSE for common integral.
#' @param korder the order of Legendre approximation.
#' @param max_iter max iteration number for IGMM.
#' @param max_tol max tolerance for iteration algorithm.
#' @param show_log bool value, TRUE for showing calculation log.
#'
#' @return \item{Rhat}{The estimated correlation coefficients.}
#' @return \item{COV}{The estimated covariance matrix for Rhat}
#'
#' @references
#' arXiv:2404.06781
#'
#' @export
#'
#' @examples
#' library(mvtnorm)
#' library(MASS)
#' set.seed(1997)
#' n = 500
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
#' out1 = est_mixedGMM(dataYX = data2,order_indx = indc)
#' print(out1$Rhat)
#' print(out1$COV)
est_mixedGMM = function(dataYX,order_indx,R0=NULL,app=TRUE,korder=2,max_iter=1000,max_tol=1e-8,show_log = FALSE){
  # moments function
  mRar_fun = function(R0yy,R0yx,R0xx){
    mr_yy = R0yy
    mr_yx = rep(R0yx,rep(si,n_y))*rep(unlist(lapply(a[order_indx],function(x){-diff(dnorm(x))})),n_y)
    mr_xx = c()
    rxx_indx = 1
    for (i2 in 1:(n_x-1)) {
      for(j2 in (i2+1):n_x){
        rxx = R0xx[rxx_indx]
        thre1 = as.numeric(a[[i2+n_y]])
        thre2 = as.numeric(a[[j2+n_y]])
        if(app){
          MPhi1 = outer(thre2,thre1,Phixy,rho=rxx,app=app)
        }else{
          MPhi1 = matrix(0,ncol = length(thre1),nrow = length(thre2))
          for (k in 1:ncol(MPhi1)) {
            for (l in 1:nrow(MPhi1)) {
              MPhi1[l,k] = Phixy(thre1[k],thre2[l],rho=rxx,app=app)
            }
          }
        }
        dim1 = dim(MPhi1)[1]
        dim2 = dim(MPhi1)[2]
        mr_xx1 = MPhi1[2:dim1,2:dim2]-MPhi1[1:(dim1-1),2:dim2]-MPhi1[2:dim1,1:(dim2-1)]+MPhi1[1:(dim1-1),1:(dim2-1)]
        mr_xx = c(mr_xx,as.numeric(mr_xx1))
        rxx_indx = rxx_indx + 1
      }
    }
    mRar = c(mr_yy,mr_yx,mr_xx)
    return(mRar)
  }
  # gradient function, dR
  G22_fun = function(R0yy,R0yx,R0xx){
    G22yy = -diag(length(R0yy))

    g22yx1 = lapply(a[order_indx],function(x){-diff(dnorm(x))})
    g22indx = cumsum(lengths(g22yx1))
    ncg22yx = length(g22indx)
    g22yx = matrix(0,g22indx[ncg22yx],ncg22yx)
    g22indx = c(0,g22indx)
    for (i in 1:ncg22yx) {
      g22yx[(g22indx[i]+1):g22indx[i+1],i] = -g22yx1[[i]]
    }
    dimyx = dim(g22yx)
    G22yx = matrix(0,dimyx[1]*n_y,dimyx[2]*n_y)
    for (i in 1:n_y) {
      G22yx[((i-1)*dimyx[1]+1):(i*dimyx[1]),((i-1)*dimyx[2]+1):(i*dimyx[2])] = g22yx
    }

    g22xx = list()
    rxx_indx = 1
    for (i2 in 1:(n_x-1)) {
      for (j2 in (i2+1):(n_x)) {
        rxx = R0xx[rxx_indx]
        thre1 = a[[i2+n_y]]
        thre2 = a[[j2+n_y]]
        mphi1 = outer(thre2,thre1,dphixy,rho=rxx)
        dim1 = dim(mphi1)[1]
        dim2 = dim(mphi1)[2]
        g22xx1 = mphi1[2:dim1,2:dim2]-mphi1[1:(dim1-1),2:dim2]-mphi1[2:dim1,1:(dim2-1)]+mphi1[1:(dim1-1),1:(dim2-1)]
        g22xx[[rxx_indx]] = -as.numeric(g22xx1)
        rxx_indx = rxx_indx + 1
      }
    }
    G22indx = cumsum(lengths(g22xx))
    ncG22xx = length(G22indx)
    G22xx = matrix(0,G22indx[ncG22xx],ncG22xx)
    G22indx = c(0,G22indx)
    for (i in 1:ncG22xx) {
      G22xx[(G22indx[i]+1):G22indx[i+1],i] = g22xx[[i]]
    }

    dimyy = dim(G22yy)
    dimyx = dim(G22yx)
    dimxx = dim(G22xx)
    G22 = matrix(0,dimyy[1]+dimyx[1]+dimxx[1],dimyy[2]+dimyx[2]+dimxx[2])
    G22[1:dimyy[1],1:dimyy[2]] = G22yy
    G22[(dimyy[1]+1):(dimyy[1]+dimyx[1]),(dimyy[2]+1):(dimyy[2]+dimyx[2])] = G22yx
    G22[(dimyy[1]+dimyx[1]+1):(dimyy[1]+dimyx[1]+dimxx[1]),(dimyy[2]+dimyx[2]+1):(dimyy[2]+dimyx[2]+dimxx[2])] = G22xx

    return(G22)
  }
  # gradient function, da
  G21_fun = function(R0yy,R0yx,R0xx){
    G21yy = array(0,dim = c(length(R0yy),sum(si)-n_x))

    G21yx = c()
    G21yx_list = lapply(a[order_indx], function(x){
      dx=x*dnorm(x)
      dx = diag(dx)
      dx = dx[-c(1,length(x)),-c(1,length(x))]
      return(rbind(dx,0)+rbind(0,-dx))
    })
    for (i in 1:n_y) {
      R0yy1 = R0yx[(n_x*(i-1)+1):(n_x*i)]
      G21yx1 = Map("*",R0yy1,G21yx_list)
      G21indx = cumsum(sapply(G21yx1, function(mat) nrow(mat)))
      ncG21yx1 = cumsum(sapply(G21yx1, function(mat) ncol(mat)))
      g21yx = matrix(0,G21indx[length(G21indx)],ncG21yx1[length(ncG21yx1)])
      G21indx = c(0,G21indx)
      ncG21yx1 = c(0,ncG21yx1)
      for (j in 1:(length(ncG21yx1)-1)) {
        g21yx[(G21indx[j]+1):G21indx[j+1],(ncG21yx1[j]+1):ncG21yx1[j+1]] = G21yx1[[j]]
      }
      G21yx = rbind(G21yx,g21yx)
    }

    G21xx = c()
    rxx_indx = 1
    for (i2 in 1:(n_x-1)) {
      for (j2 in (i2+1):n_x) {
        g21xx1 = array(0,dim = c(si[i2]*si[j2],sum(si+1)))
        rxx = R0xx[rxx_indx]
        thre1 = a[[i2+n_y]]
        thre2 = a[[j2+n_y]]
        for (k in 1:si[i2]) {
          for (l in 1:si[j2]) {
            aik = thre1[k+1]
            ajl = thre2[l+1]
            aik1 = thre1[k]
            ajl1 = thre2[l]

            daik = dnorm(aik)*(pnorm((ajl-rxx*aik)/sqrt(1-rxx^2))-pnorm((ajl1-rxx*aik)/sqrt(1-rxx^2)))
            dajl = dnorm(ajl)*(pnorm((aik-rxx*ajl)/sqrt(1-rxx^2))-pnorm((aik1-rxx*ajl)/sqrt(1-rxx^2)))
            daik1 = dnorm(aik1)*(pnorm((ajl1-rxx*aik1)/sqrt(1-rxx^2))-pnorm((ajl-rxx*aik1)/sqrt(1-rxx^2)))
            dajl1 = dnorm(ajl1)*(pnorm((aik1-rxx*ajl1)/sqrt(1-rxx^2))-pnorm((aik-rxx*ajl1)/sqrt(1-rxx^2)))

            g21xx1[(k-1)*si[j2]+l,cumsum(c(0,si+1))[i2]+k] = daik1
            g21xx1[(k-1)*si[j2]+l,cumsum(c(0,si+1))[i2]+k+1] = daik
            g21xx1[(k-1)*si[j2]+l,cumsum(c(0,si+1))[j2]+l] = dajl1
            g21xx1[(k-1)*si[j2]+l,cumsum(c(0,si+1))[j2]+l+1] = dajl
          }
        }
        rxx_indx = rxx_indx+1
        G21xx = rbind(G21xx,g21xx1)
      }
    }
    G21xx = G21xx[,(dnorm(unlist(a))!=0)]

    return(-rbind(G21yy,G21yx,G21xx))
  }
  # cov matrix of h
  Sigma_fun = function(){
    epdXh = epdX[,-cumsum(si)]
    har = unlist(lapply(a[order_indx], function(x){y=diff(pnorm(x));return(y[-length(y)])}))
    res_ha = t(epdXh) - har
    return(res_ha%*%t(res_ha)/n_sample)
  }
  # total loss function
  Loss_function = function(R0,W0){
    R0yy = R0[R0yy_indx]
    R0yx = R0[R0yx_indx]
    R0xx = R0[R0xx_indx]

    mRar = mRar_fun(R0yy=R0yy,R0yx=R0yx,R0xx=R0xx)
    mRa = mRal - mRar
    loss = t(mRa)%*%W0%*%mRa/2
    return(loss)
  }
  # total gradient function
  Gradient_function = function(R0,W0){
    R0yy = R0[R0yy_indx]
    R0yx = R0[R0yx_indx]
    R0xx = R0[R0xx_indx]

    mRar = mRar_fun(R0yy=R0yy,R0yx=R0yx,R0xx=R0xx)
    mRa = mRal - mRar
    G22 = G22_fun(R0yy=R0yy,R0yx=R0yx,R0xx=R0xx)
    GR = t(G22)%*%W0%*%mRa

    return(GR)
  }

  names_data = names(dataYX)
  dataYX = as.matrix(dataYX)
  n_sample = nrow(dataYX)
  n_yx = ncol(dataYX)
  conti_indx =setdiff(1:n_yx,order_indx)
  n_x = length(order_indx)
  n_y = length(conti_indx)
  dataY = dataYX[,conti_indx]
  dataX = dataYX[,order_indx]

  # initialized R
  if(is.null(R0)){
    R0 = cor(dataYX)
  }
  names_matrix = t(outer(names_data,names_data,function(x, y) paste(x, y, sep = "-")))

  R0yym = R0[1:n_y,1:n_y]
  R0yy = R0yym[lower.tri(R0yym)]
  names(R0yy) = (names_matrix[1:n_y,1:n_y])[lower.tri(R0yym)]

  R0yxm = R0[(n_y+1):n_yx,1:n_y]
  R0yx = as.numeric(R0yxm)
  names(R0yx) = as.vector(names_matrix[(n_y+1):n_yx,1:n_y])

  R0xxm = R0[(n_y+1):n_yx,(n_y+1):n_yx]
  R0xx = R0xxm[lower.tri(R0xxm)]
  names(R0xx) = (names_matrix[(n_y+1):n_yx,(n_y+1):n_yx])[lower.tri(R0xxm)]

  R0 = c(R0yy,R0yx,R0xx)
  R0yy_indx = 1:length(R0yy)
  R0yx_indx = 1:length(R0yx)+length(R0yy)
  R0xx_indx = 1:length(R0xx)+length(R0yy)+length(R0yx)

  # estimate thresholds a as a list
  a = list()
  for (i in order_indx) {
    a[[i]] = est_thre(dataYX[,i])
  }
  si = unlist(lapply(a, length))[order_indx]-1
  si_list = 1:sum(si)
  si_list = split(si_list,rep(1:n_x,si))

  # initialized W
  mxx = si%*%t(si)
  n_equations = n_y*(n_y-1)/2 + n_y*sum(si) + sum(mxx[lower.tri(mxx)])
  W0 = diag(n_equations)

  # gyy left
  gl_yy = c()
  for (i1 in 1:(n_y-1)) {
    for (j1 in (i1+1):n_y) {
      gl_yy = rbind(gl_yy,dataY[,i1]*dataY[,j1])
    }
  }

  # expanded one-hot X
  epdX = c()
  for (i in 1:n_x) {
    epdX = cbind(epdX,model.matrix(~factor(dataX[,i])-1))
  }

  # gyx left
  gl_yx = c()
  for (i1 in 1:n_y) {
    for (i2 in 1:n_x) {
      gl_yx = rbind(gl_yx,t(dataY[,i1]*epdX[,si_list[[i2]]]))
    }
  }

  # gxx left
  gl_xx = c()
  for (i2 in 1:(n_x-1)) {
    for (j2 in (i2+1):(n_x)) {
      glxx1 = c()
      for (k in 1:si[i2]) {
        glxx1 = rbind(glxx1,t((epdX[,si_list[[i2]]])[,k]*epdX[,si_list[[j2]]]))
      }
      gl_xx = rbind(gl_xx,glxx1)
    }
  }

  # gRa left and mRa left, since they are unchanged, so we calculate them in advance
  gRal = rbind(gl_yy,gl_yx,gl_xx)
  mRal = rowMeans(gRal)

  # the iterative GMM (renew the W)
  n_iter = 0
  ebound=10

  while ((n_iter<=max_iter)&(ebound>=max_tol)) {
    Rhat = optim(R0,fn=Loss_function,gr=Gradient_function,W0=W0,method = 'BFGS')$par
    resgRa = gRal - mRar_fun(Rhat[R0yy_indx],Rhat[R0yx_indx],Rhat[R0xx_indx])
    Omegahat = resgRa%*%t(resgRa)/n_sample
    What = ginv(Omegahat)
    # ebound = sum((Rhat-R0)^2)+sum((Omegahat-ginv(W0))^2)
    if(app){
      ebound = max((Omegahat-ginv(W0))^2)/length(Rhat)^2
    }else{
      ebound = sum((Rhat-R0)^2)+sum((Omegahat-ginv(W0))^2)
    }
    R0 = Rhat
    W0 = What
    n_iter = n_iter+1
    if(show_log){
      print(n_iter)
      print(Rhat)
      print(ebound)
    }
  }

  ## modified the variance of 2-step estimators
  # calculate Ehh'
  Sigmahat = Sigma_fun()
  # gradient
  g22 = G22_fun(Rhat[R0yy_indx],Rhat[R0yx_indx],Rhat[R0xx_indx])
  g21 = G21_fun(Rhat[R0yy_indx],Rhat[R0yx_indx],Rhat[R0xx_indx])

  Lambda = ginv(t(g22)%*%What%*%g22)
  Gamma = t(g22)%*%What%*%g21
  COV = Lambda+Lambda%*%Gamma%*%Sigmahat%*%t(Gamma)%*%Lambda
  COV = COV/n_sample
  colnames(COV) = names(Rhat)
  rownames(COV) = names(Rhat)
  return(list(Rhat=Rhat,
              What = What,
              Sigmahat = Sigmahat,
              g22 = g22,
              g21 = g21,
              Lambda = Lambda,
              COV = COV,
              names_data = names_data
  ))

}
