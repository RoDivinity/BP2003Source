#load required functions
source('bandw.R')
source('kern.R')

#function to compute long-run covariance matrix of vmat
jhatpr = function(vmat) {
  ####
  # vmat = variance matrix
  # ****
  # jhat = long run covariance matrix
  ###
  
  nt = dim(vmat)[1]
  d = dim(vmat)[2]
  #automatic bandwidth selection
  st = bandw(vmat)
  #lag 0 covariance  
  jhat = t(vmat) %*% vmat
  
  #forward sum
  for (j in seq(1,nt-1,1)){
    vl = vmat[(j+1):nt,,drop=FALSE]
    vr = vmat[1:(nt-j),,drop=FALSE]
    vmat_s = t(vl) %*% vr
    jhat = jhat + as.vector(kern(j/st)) * vmat_s
  }
  
  #backward sum
  for (j in seq(1,nt-1,1)) {
    vl = vmat[seq(1,nt-j,1),,drop=FALSE]
    vr = vmat[seq(j+1,nt,1),,drop=FALSE]
    vmat_s = t(vl) %*% vr
    jhat = jhat + as.vector(kern(j/st)) * vmat_s
  }
  
  #small sample correction
  jhat = jhat/(nt-d)
  
  return(jhat)
}
# #test code
# tMatrix = matrix(c(1:100),20,5)
# jhatpr(tMatrix)

