#load required functions
source('/Users/linhnguyen/Desktop/BP2003/date/parti.R')
source('/Users/linhnguyen/Desktop/BP2003/onebp.R')

# function to compute the supF(l+1|l) test l=nseg-1. The l breaks under
# the null are taken from the global minimization (in dt).
spflp1 = function(bigvec,dt,nseg,y,z,h,q,prewhit,robust,x,p,hetdat,hetvar){
  ##########
  # bigvec = associated SSR of estimated break date under H0
  # dt = estimated date under H0
  # nseg = number of segment under H1
  # y = dep vars
  # z,x = ind vars
  # p,q = # of regressors
  # prewhit,robust,hetdat,hetvar = options on residuals/errors
  # *******
  # maxf: maximum value of F test
  # newd: additional date in l+1 
  ########
  ssr = matrix(0L,nrow = nseg, ncol = 1)
  ftestv = matrix(0L, nrow = nseg, ncol = 1)
  bigT = dim(z)[1]
  dv = matrix(0L,nrow = nseg+1, ncol = 1)
  dv[2:nseg,1] = dt
  dv[nseg+1,1] = bigT
  ds = matrix(0L, nrow = nseg, ncol = 1)
  
  i_n = 0
  for (is in 1:nseg){
    length = dv[is+1,1] - dv[is,1]
  
    if(length >= 2*h){
      if (p == 0){
        out = parti(dv[is,1]+1,dv[is,1]+h,dv[is+1,1]-h,dv[is+1,1],bigvec,bigT)
        ssr[is,1] = out$ssrmin
        ds[is,1] = out$dx
        y_test = y[seq(dv[is,1]+1,dv[is+1,1],1),1,drop=FALSE]
        z_test = z[seq(dv[is,1]+1,dv[is+1,1],1),,drop=FALSE]
        ftestv[is,1] = pftest(y_test,z_test,1,q,length,ds[is,1,drop=FALSE]-dv[is,1,drop=FALSE],
                              prewhit,robust,0,p,hetdat,hetvar)
      }
      else{
        out = onebp(y,z,x,h,dv[is,1]+1,dv[is+1,1])
        ssr[is,1] = out$ssrind
        ds[is,1] = out$bd
        y_test = y[seq(dv[is,1]+1,dv[is+1,1],1),1,drop=FALSE]
        z_test = z[seq(dv[is,1]+1,dv[is+1,1],1),,drop=FALSE]
        x_test = x[seq(dv[is,1]+1,dv[is+1,1],1),,drop=FALSE]
        ftestv[is,1] = pftest(y_test,z_test,1,q,length,ds[is,1,drop=FALSE]-dv[is,1,drop=FALSE],
                              prewhit,robust,x_test,p,hetdat,hetvar)
      }
    }
    else {
      i_n = i_n+1
      ftestv[is,1] = 0.0
    }
  }
  
  if (i_n == nseg) {
    print(paste('Given the location of the breaks from the global optimization with',
       nseg,'breaks there was no more place to insert an additional breaks that satisfy the minimal length requirement.'))
  }
  
  maxf = max(ftestv[1:nseg,1])
  newd = ds[which.max(ftestv[1:nseg,1]),1]
  out = list('maxf' = maxf, 'newd' = newd)
  return(out)
}