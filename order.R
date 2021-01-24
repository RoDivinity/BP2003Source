#function to calculate the order using BIC and the criterion of 
#Liu, Wu and Zidek

order = function(ssr0,globl,bigT,m,q) {
  ####
  # ssr0 = SSR under null hypothesis of no break
  # globl= SSR of current model
  # bigT = sample size T
  # m = maximum number of break
  # q = # of regressors z
  # ********
  # mBIC: number of breaks selected by BIC
  # mLWZ: number of breaks selected by LWZ
  ####
  
  delta0 = 0.1 #optimal parameters in LWZ paper
  c0 = 0.299
  glob= matrix(0L, nrow = m+1, ncol=1)
  glob[1,1] = ssr0
  glob[seq(2,m+1),1] = globl
  
  bic = matrix(0L,nrow = m+1, ncol = 1)
  lwz = matrix(0L,nrow = m+1, ncol = 1)
  
  for (i in seq(1,m+1)){
    bic [i,1] = log(glob[i,1]/bigT) + log(bigT)*(i-1)*(q+1)/bigT
    lwz[i,1] = log(glob[i,1]/(bigT-i*q-i+1)) + 
      ((i-1)*(q+1)*c0*(log(bigT))^(2+delta0))/bigT
    
    print(paste('With',i-1,'breaks:'))
    print(paste('BIC=',bic[i,1]))
    print(paste('LWZ=',lwz[i,1]))
  }
  
  mBIC = which.min(bic) - 1
  mLWZ = which.min(lwz) - 1
  out = list('mBIC' = mBIC, 'mLWZ' = mLWZ)
  return(out)
}