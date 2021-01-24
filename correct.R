#load required functions
source('OLS.R')
source('jhatpr.R')
#function to compute robust standard errors
correct = function (reg,res,prewhit){
  ####
  # reg = regression model
  # res = residuals
  # prewhit = option of using prewhitening
  # *****
  # hac = robust standard error
  ####
  
  #initialize storage
  nt = dim(reg)[1]
  d = dim(reg)[2]
  b = matrix(0L, nrow = d, ncol = 1)
  bmat = matrix(0L, nrow = d, ncol = d)
  vstar = matrix(0L, nrow = nt-1, ncol = d)
  vmat = matrix(0L, nrow = nt, ncol = d)
  
  #Construct matrix z_t * u_t
  for (i in 1:d){
    vmat[,i] = reg[,i,drop=FALSE] * res #element-wise multiplication/ need to check
  }
  
  #Prewhitening to matrix vmat by filtering with a VAR(1). If prewhit = 0, skip
  if (prewhit == 1){
    for(i in 1:d){
      #check carefully if errors
      b = OLS(vmat[seq(2,nt,1),i,drop=FALSE],vmat[seq(1,nt-1,1),,drop=FALSE])
      bmat[i,] = t(b) 
      vstar[,i] = vmat[seq(2,nt,1),i,drop=FALSE] - vmat[seq(1,nt-1,1),,drop=FALSE] %*% b
    }
    #kernel on the residuals
    jh = jhatpr(vstar)
    
    #recolor
    inv = solve( diag(1,d) - bmat)
    hac = inv %*% jh %*% t(inv)
  }
  else {
    
    hac = jhatpr(vmat)
  }
    
  return (hac)
}

