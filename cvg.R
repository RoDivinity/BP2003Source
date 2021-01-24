#load required functions
source('/Users/linhnguyen/Desktop/BP2003/funcg.R')

#function to calculate critical values for break date
cvg = function(eta,phi1s,phi2s){
  #cvec = critical values of break dates
  
  cvec = matrix(0L,nrow = 4, ncol = 1)
  a=phi1s/phi2s
  gam=((phi2s/phi1s)+1)*eta/2
  b=sqrt(phi1s*eta/phi2s)
  deld=sqrt(phi2s*eta/phi1s)+b/2
  alph=a*(1+a)/2
  bet=(1+2*a)/2
  sig = c(0.025,0.05,0.95,0.975)
  
  isig = 1
  while(isig <= 4){
    #initialize upper, lower bound and critical value
    upb = 2000
    lwb = -2000
    crit = 999999
    cct = 1
    
    while(abs(crit) >= 0.000001) {
      cct = cct + 1
      if (cct > 100){
        print('the procedure to get critical values for the break dates has reached the upper bound on the number of iterations. This may happens in the procedure cvg. The resulting confidence interval for this break date is incorect')
        break
      }
      else{
        xx = lwb + (upb-lwb)/2
        pval=funcg(xx,bet,alph,b,deld,gam)
        crit = pval - sig[isig]
        if (crit <= 0) {
          lwb = xx}
        else {
          upb = xx}
      }
    }
    cvec[isig,1] = xx
    isig = isig+ 1
  }
  
  return(cvec)
}