source('/Users/linhnguyen/Desktop/BP2003/OLS.R')
#recursive residuals from dataset from: start - end: last
ssr = function(start,y,z,h,last) {
  ###
  #start: starting index of the sample used
  #last: ending index of the last segment considered
  #y: dependent vars
  #z: matrix of time-variant regressors (dimension q)
  #h: minimal segment length
  #*****************
  #out: vector of sum of SSR of length last-start+1
  ###
  out = matrix(0L,nrow = last, ncol = 1)
  
  #initialize the recursion
  index_0 = seq(start,start+h-1,1)
  z_0 = z[index_0,,drop=FALSE]
  y_0 = y[index_0,,drop=FALSE]
  inv1 = solve(t(z_0) %*% z_0)
  delta1 = OLS(y_0,z_0) #calculate initial beta
  res = y_0 - z_0 %*% delta1 #calculate initial residuals
  out[start+h-1,1] = t(res) %*% res # store initial SSR
  
  #loop to construct the recursive residuals and update SSR
  r = start+h
  while (r <= last) {
    v = y[r,1] - z[r,]%*%delta1  #residual of next obs
    v = drop(v)
    invz = inv1 %*% matrix(t(z[r,]))
    f = 1 + z[r,,drop=FALSE] %*% invz
    f = drop(f)
    delta2 = delta1 + (invz * v) / f
    inv2 = inv1 - (invz %*% t(invz)) / f
    inv1 = inv2
    delta1 = delta2
    out_t = out[r-1,1] + v*v/f
    out[r,1] = out_t
    r=r+1
  }
  return(out)
}

# #### TEST CODE ####
# y = matrix(c(1:12),nrow=12)
# start = 2
# last = 10
# h = 2
# z = matrix(c(1:24),nrow = 12, ncol =2)
# print(y)
# print(z)
#out1 = ssr(1,y,z,h,bigT)
#out = ssr(1,y,z,h,bigT)
#print(out)
#print((out-out1))
