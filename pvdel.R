#load required functions
source('diag_par.R')
source('plambda.R')
source('psigmq.R')
source('kron.R')
source('correct.R')
#function to compute covariance matrix of estimate delta
pvdel = function(y,z,i,q,bigT,b,prewhit,robust,x,p,withb,hetdat,hetvar) {
  #####
  # y = dep vars
  # z = ind vars
  # i = max num breaks
  # q = # of z regressors
  # bigT = sample period T
  # b = date vector
  # p = # of x regressors
  # prewhit,robust,withb,hetdat,hetvar = options for assumptions
  #********
  # vdel = covariance matrix of delta
  ####
  ev = matrix(0L , nrow = i+1, ncol = 1)
  zbar = diag_par(z,i,b)
  
  if (p == 0){
    delv = OLS(y,zbar)
    res = y - zbar %*% delv
    reg = zbar
  }
  else{
      delv = OLS(y, cbind(zbar,x))
      res = y - cbind(zbar,x) %*% delv
      
      if (withb == 0) {reg = zbar - x %*% solve(t(x) %*% x) %*% t(x) %*% zbar}
      else {reg = cbind(x,zbar)}
  }
  
  
  if (robust == 0) {
    #testing with no serial correlation in errors 
    if(p==0) {
      if (hetdat==1 && hetvar == 0){
        sig = t(res) %*% res / bigT
        vdel = drop(sig) * solve(t(reg) %*% reg)
        
      }
      if (hetdat == 1 && hetvar == 1){
        sig = psigmq(res,b,q,i,bigT)
        vdel = kron(sig,diag(1,q)) %*% solve(t(reg) %*% reg)
      }
      if (hetdat == 0 && hetvar == 0){
        lambda = plambda(b,i,bigT)
        sig = t(res) %*% res / bigT
        vdel = drop(sig) * solve(kron(lambda,t(z) %*% z))
      }
      if (hetdat == 0 && hetvar == 1) {
        lambda = plambda(b,i,bigT)
        sig = psigmq(res,b,q,i,bigT)
        vdel = kron(sig,diag(1,q)) %*% solve(kron(lambda, t(z) %*% z))
      }
    }
    else {
      if (hetdat == 0) {
        print(paste('hetdat == 0 is not allowed','vdel is returned zeros'))
        vdel = matrix (0L, nrow = q*(i+1), ncol = q*(i+1))
        }
      if (hetdat == 1 && hetvar == 0) {
        sig = t(res) %*% res / bigT
        vdel = drop(sig) * solve(t(reg) %*% reg)
      }
      if (hetdat == 1 && hetvar == 1) {
        wbar = diag_par(reg,i,b)
        ww = t(wbar) %*% wbar
        sig = psigmq(res,b,q,i,bigT)
        gg = matrix (0L, nrow = (i+1)*q + p*withb, ncol = (i+1)*q + p*withb)
        ie = 1
        while(ie <= i + 1){
          index = seq((ie-1)*((i+1)*q+p*withb)+1,ie*((i+1)*q+p*withb),1)
          increment = sig[ie,ie] * ww[index,index]
          gg = gg + increment
          ie = ie + 1
        }
        vdel = solve(t(reg) %*% reg) %*% gg %*% solve(t(reg) %*% reg)
      }
    }
  }
  else {
   #testing with serial correlation in errors
    if(hetdat == 0) {
      print(paste('hetdat = 0 is not allowed','vdel is returned zeros'))
      vdel = matrix(0L,nrow = q*(i+1),ncol = q*(i+1))
    }
  
    if (p==0){
      if (hetvar == 1){
        hac = matrix(0L, nrow = q*(i+1), ncol = q*(i+1))
        vdel = matrix(0L, nrow = q*(i+1), ncol = q*(i+1))
        ind_temp = seq(1,b[1,1],1)
        temp = correct(z[ind_temp,,drop=FALSE], res[ind_temp,1,drop=FALSE] , prewhit)
        hac[1:q,1:q] = b[1,1] * temp
        if(dim(b)[1] > 1){
        for (j in 2:i) {
          ind_hac = seq((j-1)*q+1,j*q,1)
          ind_temp = seq(b[j-1,1]+1,b[j,1],1)
          temp = correct(z[ind_temp,,drop=FALSE], res[ind_temp,1,drop=FALSE], prewhit)
          hac[ind_hac,ind_hac] = (b[j,1] - b[j-1,1]) * temp
        }
        }
        ind_hac = seq(i*q+1,(i+1)*q,1)
        ind_temp = seq(b[i,1]+1,bigT,1)
        temp = correct(z[ind_temp,,drop=FALSE],res[ind_temp,1,drop=FALSE],prewhit)
        hac[ind_hac,ind_hac] = (bigT - b[i,1]) * temp
        vdel = solve(t(reg) %*% reg) %*% hac %*% solve(t(reg) %*% reg)
        
      }
      else {
        hac = correct(z,res,prewhit)
        lambda = plambda(b,i,bigT)
        vdel = bigT * solve(t(reg) %*% reg) %*% kron(lambda,hac) %*% solve(t(reg) %*% reg)
      }
    }
    else{
      hac = correct(reg,res,prewhit)
      vdel = bigT * solve(t(reg) %*% reg) %*% hac %*% solve(t(reg) %*% reg) 
    }
  }
  
  return(vdel)
}

