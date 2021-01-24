#function to compute a diagonal matrix of dimension i+1 with i-th entry
#is the estimate of the variance of residuals of segment i

psigmq = function (res,b,q,i,nt) {
  ####
  # res = residual
  # b = betas
  # q = # of regressors z (time-variant coeffs)
  # i = segment i
  # nt = total length
  # ****
  # sigmat = (i+1)x(i+1) diagonal matrix with i-th entry = estimate of variance
  ###
  
  sigmat = matrix(0L, nrow = i+1, ncol = i+1)
  resid_t = res[seq(1,b[1,1],1),1,drop=FALSE] #check this line if has problem
  sigmat[1,1] = t(resid_t) %*% resid_t / b[1,1]
  
  kk = 2
  while (kk <= i) {
    #check this loop for problem if errors
    bf = b[kk-1,1]
    bl = b[kk,1]
    resid_temp = res[seq(bf,bl,1),1]
    sigmat[kk,kk] = t(resid_temp) %*% resid_temp / (bl - bf)
    kk = kk+1
  }
  
  res_f = res[seq(b[i,1]+1,nt,1),1,drop = FALSE]
  sigmat[i+1,i+1] = t(res_f) %*% res_f / (nt - b[i,1])
  
  return(sigmat)
}