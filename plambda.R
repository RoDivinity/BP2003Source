#function to construct a diagonal matrix of dimension m+1 
# with i-th entry (T_{i} - T_{i-1} )/ T

plambda = function(b,m,bigT) {
  ####
  # b = beta
  # m = max num of break
  # bigT = sample period T
  # ***********
  # lambda = (m+1)x(m+1) diagonal matrix with i-th entry (T_{i} - T_{i-1} )/ T
  ####
  
  lambda = matrix(0L , nrow = m+1, ncol = m+1)
  lambda[1,1] = b[1,1]/bigT
  
  k = 2
  while (k<= m){
    lambda[k,k] = ( b[k,1] - b[k-1,1] ) / bigT
    k = k+1
  }
  lambda[m+1,m+1] = ( bigT - b[m,1] ) / bigT
  
  return (lambda)
}