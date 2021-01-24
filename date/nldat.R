#load required functions
source('/Users/linhnguyen/Desktop/BP2003/date/dating.R')
source('/Users/linhnguyen/Desktop/BP2003/diag_par.R')
source('/Users/linhnguyen/Desktop/BP2003/OLS.R')
#function to estimate partial structural change model
nldat = function(y,z,x,h,m,p,q,bigT,fixb,eps,maxi,betaini,printd){
  ###
  # y: dependent var
  # z: time-variant coeff independent vars
  # x: time-invariant coeff independent vars
  # h: minimal length of segment
  # m: MAX num of breaks
  # p: # of x regressors
  # q: # of z regressors
  # bigT: sample periods T
  # fixb: use given betas option
  # eps: trimming
  # maxi: maximum number of iterations
  # betaini: initial betas
  # printd: print results option
  #********
  # glb: minimum SSR
  # datevec: optimal dates 
  # bigvec: associated SSR
  ###
  
  #initialize storage
  glb = matrix(0L , nrow = m, ncol = 1)
  globnl = matrix(0L , nrow = m, ncol = 1)
  datevec = matrix(0L, nrow = m, ncol = m)
  datenl = matrix(0L, nrow = m, ncol = m)
  
  #initialize current max break
  mi = 1
  while(mi <= m){
    if (printd == 1){
      print(paste('Breaks of model',mi))
    }
    if (fixb == 0){
      qq = p+q
      zz = cbind(x,z)
      
      #initial estimate of the model 
      out = dating(y,zz,h,mi,qq,bigT)
      date = out$datevec[1:mi,mi,drop=FALSE]
      bigvec = out$bigvec
      
      #partition regressors with initial estimated date
      xbar = diag_par(x,mi,date)
      zbar = diag_par(z,mi,date)
      
      #calculate initial values of estimate
      teta = OLS(y,cbind(zbar,xbar))
      delta1 = teta[seq(1,q*(mi+1),1),1,drop=FALSE]
      beta1 = OLS(y - zbar %*% delta1, x)
      
      #calculate initial SSR of the model
      resid1 = y - x%*%beta1 - zbar%*%delta1
      ssr1 = t(resid1) %*% resid1
      
      if(printd==1) {
        print('The iterations are initialized with')
        prmatrix(delta1)
        prmatrix(beta1)
        print('With break date')
        print(date)}
    }
    else {
      beta1 = betaini
      ssr1 = -5
    }
    
    #start the iterations
   
    
    length = 999999999999
    i = 1
    
    while (length > eps) {
      
      out = dating(y-x%*%beta1,z,h,mi,q,bigT)
      #store the date vector for current max
      date = out$datevec[1:mi,mi,drop=FALSE]
      bigvec = out$bigvec
      zbar = diag_par(z,mi,date)
      
      #update estimates based on new partition
      teta1 = OLS(y,cbind(x,zbar))
      beta1 = teta1[seq(1,p,1),1,drop=FALSE]
      delta1 = teta1[seq(p+1,p+q*(mi+1),1),1,drop=FALSE]
     
     
      #check convergence condition 
      resid_n = y - cbind(x,zbar) %*% teta1
      ssrn = t(resid_n) %*% resid_n
      length = abs(ssrn - ssr1)
      
      if(printd==1){
        print(paste('Iteration',i))
      }
      
      #check upper bound of iterations and update
      if (i >= maxi){
        print('Number of iterations has reached MAX ITER')
      }
      else {
        i = i+1
        ssr1 = ssrn
        glb[mi,1]=ssrn
        datevec[1:mi,mi] = date
      }
      
    }
    #finished current max breaks & update
    mi = mi + 1
     
  }
  
  
  out = list('glb' = glb, 'datevec' = datevec, 'bigvec' = bigvec)
  return(out)
}



###Test code####
