#load required functions
source('/Users/linhnguyen/Desktop/BP2003/cvg.R')
source('/Users/linhnguyen/Desktop/BP2003/mpower.R')
source('/Users/linhnguyen/Desktop/BP2003/OLS.R')

# function to compute confidence intervals for the break dates based on
# the "shrinking shifts" asymptotic framework

interval = function(y,z,zbar,b,q,m,robust,prewhit,hetomega,hetq,x,p){
  ######
  #
  ##
  # bound = CIs for break dates
  #####
  
  cvf = matrix(0L,nrow = 4, ncol = 1)
  nt = dim(z)[1]
  bound = matrix(0L, nrow = m, ncol = 4)
  bf = matrix(0L, nrow = m+2, ncol = 1)
  
  #specify resid
  if(p==0){
    delta = OLS(y,zbar)
    res = y - zbar %*% delta
  }
  else{
    dbdel = OLS(y, cbind(zbar,x))
    res = y - cbind(zbar,x) %*% dbdel
    delta = dbdel[seq(1,(m+1)*q),1]
  }
  
  bf[1,1] = 0
  bf[seq(2,m+1),1] = b[seq(1,m),1,drop=FALSE]
  bf[m+2,1] = nt
  
  for(i in 1:m){
    delv = delta[seq(i*q+1,(i+1)*q),1,drop=FALSE] - 
      delta[seq((i-1)*q+1,i*q),1,drop=FALSE]
    
    if (robust == 0){
      if (hetq == 1){
        index = seq(bf[i,1]+1,bf[i+1,1])
        index1 = seq(bf[i+1,1]+1,bf[i+2,1])
        qmat = t(z[index,,drop=FALSE]) %*% z[index,,drop=FALSE] / (bf[i+1,1] - bf[i,1])
        qmat1 = t(z[index1,,drop=FALSE]) %*% z[index1,,drop=FALSE] / (bf[i+2,1] - bf[i+1,1])
      }
      else{
        qmat = t(z) %*% z / nt
        qmat1 = qmat
      }
      
      if(hetomega==1){
        index = seq(bf[i,1]+1,bf[i+1,1])
        index1 = seq(bf[i+1,1]+1,bf[i+2,1])
        phi1s = t(res[index,1,drop=FALSE]) %*% res[index,1,drop=FALSE] / (bf[i+1,1] - bf[i,1])
        phi2s = t(res[index1,1,drop=FALSE]) %*% res[index1,1,drop=FALSE] / (bf[i+2,1] - bf[i+1,1])
      }
      else{
        phi1s=t(res)%*%res/nt;
        phi2s=phi1s;
      }
      eta = t(delv) %*% qmat1 %*% delv /  (t(delv) %*% qmat %*% delv)
      cvf=cvg(eta,phi1s,phi2s)
      
      a = (t(delv) %*% qmat %*% delv) / phi1s
      
      bound[i,1] = b[i,1] - cvf[4,1]/a
      bound[i,2] = b[i,1] - cvf[1,1]/a
      bound[i,3] = b[i,1] - cvf[3,1]/a
      bound[i,4] = b[i,1] - cvf[2,1]/a
      
    }
    else{
      #robust SE
      if (hetq == 1){
        #check if the same with robust = 0 & hetq = 1
        index = seq(bf[i,1]+1,bf[i+1,1])
        index1 = seq(bf[i+1,1]+1,bf[i+2,1])
        qmat = t(z[index,,drop=FALSE]) %*% z[index,,drop=FALSE] / (bf[i+1,1] - bf[i,1])
        qmat1 = t(z[index1,,drop=FALSE]) %*% z[index1,,drop=FALSE] / (bf[i+2,1] - bf[i+1,1])
      }
      else{
        qmat=t(z)%*%z/nt;
        qmat1=qmat;
      }
      
      if(hetomega == 1){
        index = seq(bf[i,1]+1,bf[i+1,1])
        index1 = seq(bf[i+1,1]+1,bf[i+2,1])
        omega = correct(z[index,,drop=FALSE],res[index,1,drop=FALSE],prewhit)
        omega1 = correct(z[index1,,drop=FALSE],res[index1,1,drop=FALSE],prewhit)
      }
      else{
        omega=correct(z,res,prewhit);
        omega1=omega
      }
      #check matrix operations
      phi1s = t(delv) %*% omega %*% delv / (t(delv) %*% qmat %*% delv)
      phi2s = t(delv) %*% omega1 %*% delv / (t(delv) %*% qmat %*% delv)
      
      eta =  t(delv) %*% qmat1 %*% delv / (t(delv) %*% qmat %*% delv)
      
      cvf = cvg(eta,phi1s,phi2s)
      
      a = mpower(t(delv) %*% qmat %*% delv,2) / (t(delv) %*% omega %*% delv)
      
      bound[i,1] = b[i,1] - cvf[4,1]/a
      bound[i,2] = b[i,1] - cvf[1,1]/a
      bound[i,3] = b[i,1] - cvf[3,1]/a
      bound[i,4] = b[i,1] - cvf[2,1]/a
    }
    
    #round inside or outside loop?
    bound[,1] = round(bound[,1,drop=FALSE])
    bound[,2] = round(bound[,2,drop=FALSE])+1
    bound[,3] = round(bound[,3,drop=FALSE])
    bound[,4] = round(bound[,4,drop=FALSE])+1
  } 
  
  return (bound)
}