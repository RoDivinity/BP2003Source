#load required functions
source('OLS.R')
source('mdiv.R')
source('mpower.R')
#function to calculate bandwidth based on AR(1) approximation for each
#vector of the matrix vhat, with equal weight 1
bandw = function(vhat) {
  #####
  # vhat = vector of matrix \hat{v}
  # *****
  # st = bandwidth
  ####
  nt = dim(vhat)[1]
  d = dim(vhat)[2]
  
  a2n = 0
  a2d = 0
  
  for (i in seq(1,d,1)){
    b = OLS(vhat[seq(2,nt,1),i,drop=FALSE] , vhat[seq(1,nt-1,1),i,drop=FALSE]) #obtain AR(1) coef
    res = vhat[seq(2,nt,1),i,drop=FALSE] - vhat[seq(1,nt-1,1),i,drop=FALSE] %*% b
    sig = t(res) %*% res
    sig = sig/(nt-1)
    a2n = a2n + 4 * b * b * sig * sig / (1-b)^8
    a2d = a2d + sig * sig / (1-b)^4
  }
  
  a2 = a2n/a2d
  st = 1.3221 * (a2 * nt) ^.2
  return (st)
}

