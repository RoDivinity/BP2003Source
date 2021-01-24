#load required functions
source('/Users/linhnguyen/Desktop/BP2003/diag_par.R')
source('/Users/linhnguyen/Desktop/BP2003/OLS.R')

#function to compute optimal one break partition in partial structural
#change model by searching over all possible breaks
onebp = function(y,z,x,h,start,last) {
  #### 
  # y = dep vars
  # z = ind vars with coeffs allowed to change
  # x = ind vars with coeffs constant
  # h = minimal segment length
  # start = initial date to search
  # last = last date to search
  # *************
  # ssrind = associated SSR of break date
  # bd = break date
  ####
  
  ssrind = 999999999999999
  i = matrix(h,ncol=1)
  
  while(i <= last - start + 1 - h){
    
    zb = diag_par(z[start:last,,drop=FALSE],1,i)
    y_reg = y[start:last,1]
    x_reg = cbind(x[start:last,],zb)
    bb = OLS(y_reg,x_reg) #check if beta estimate correct
    resid = y_reg - x_reg %*% bb
    ssrn = t(resid) %*% resid
    
    if (ssrn < ssrind){
      ssrind = ssrn
      bdat = i
    }
    
    i = i + 1
  }
  bd = bdat + start - 1
  out = list('ssrind' = ssrind,'bd' = bd)
  return(out)
}