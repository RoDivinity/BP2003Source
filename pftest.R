#load required functions
source('pvdel.R')
source('kron.R')
source('OLS.R')
source('diag_par.R')

#function to conduct supF test for i breaks
pftest = function(y,z,i,q,bigT,datevec,prewhit,robust,x,p,hetdat,hetvar){
  #######
  # y = dep vars
  # z = ind vars with time-variant coefs
  # q = num of z regressors
  # x = ind vars with time-invariant coefs
  # p = num of x regressors
  # i = num of max breaks
  # bigT = sample period T
  # datevec = estimated date
  # prewhit,robust,hetdat,hetvar = options for assumptions on error terms
  # ***********
  # ftest = supF test results
  ######
  
  #construct matrix R
  rsub = matrix(0L, nrow = i , ncol = i+1)
  j = 1
  while(j<=i){
    rsub[j,j] = -1
    rsub[j,j+1] = 1
    j=j+1
  }
  rmat = kron(rsub,diag(1,q))
  date = datevec[seq(1,i,1),i,drop=FALSE]
  zbar = diag_par(z,i,date)
  
  if (p==0){
    delta = OLS(y,zbar)
  }
  else {
    dbdel = OLS(y, cbind(zbar,x))
    delta = dbdel[seq(1,(i+1)*q) , 1,drop=FALSE]
  }
  
  vdel = pvdel(y,z,i,q,bigT,date,prewhit,robust,x,p,0,hetdat,hetvar)
  fstar = t(delta) %*% t(rmat) %*% solve(rmat %*% vdel %*% t(rmat)) %*%
    rmat %*% delta
  ftest = (bigT - (i+1)*q - p) %*% fstar / (bigT*i)
  
  return(ftest)
}

