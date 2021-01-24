#load required functions
source('/Users/linhnguyen/Desktop/BP2003/date/parti.R')
source('/Users/linhnguyen/Desktop/BP2003/date/ssr.R')
#function to calculate break points that globally minimizes SSR
dating = function(y,z,h,m,q,bigT){
  ####
  # y = dependent var 
  # z = independent vars matrix
  # h = minimal length segment
  # m = maximum num of breaks
  # q = # of regressors z
  # bigT = sample period T
  # ************
  # glb = min global SSR
  # datevec = vector of optimal dates
  # bigvec =  associated SSRs
  ###
  
  #initialize arrays to store results
  datevec = matrix(0L, nrow = m, ncol = m)
  optdat = matrix(0L, nrow = bigT, ncol = m)
  optssr = matrix(0L, nrow = bigT, ncol = m)
  dvec = matrix(0L, nrow = bigT, ncol = 1)
  glb = matrix(0L,nrow = m, ncol = 1)
  bigvec = matrix(0L,nrow = bigT*(bigT+1)/2,ncol = 1)
 
  #calculate all possible SSR and store
  for (i in 1:(bigT-h+1)) {
    vecssr = ssr(i,y,z,h,bigT)
    bigvec[seq((i-1)*bigT+i-(i-1)*i/2, i*bigT - (i-1)*i/2 ,1),1] = vecssr[seq(i,bigT,1),1]
  }
  
  #base case: 1 break
  if (m == 1) {
    out = parti(1,h,bigT-h,bigT,bigvec,bigT)
    datevec[1,1] = out$dx
    glb[1,1] = out$ssrmin
  }
  #
  #more than 1 break
  else {
    #change the end point from smallest to full sample T, with case m = 1
    for (j1 in seq(2*h,bigT,1)){
      out = parti(1,h,j1-h,j1,bigvec,bigT)
      optssr[j1,1] = out$ssrmin
      optdat[j1,1] = out$dx
    }
    glb[1,1] = optssr[bigT,1]
    datevec[1,1] = optdat[bigT,1]
    
    #with case m >= 2
    for (ib in 2:m){
      if (ib == m) { 
        jlast = bigT
        for (jb in seq(ib*h,jlast-h,1)) {
          dvec[jb,1] = optssr[jb,ib-1] + bigvec[(jb+1)*bigT - jb*(jb+1)/2,1]
        }
        optssr[jlast,ib] = matrix(t(min(dvec[seq(ib*h,jlast-h,1)])))
        optdat[jlast,ib] = ib*h-1 + which.min(dvec[seq(ib*h,jlast-h,1)])
      }
      
      else {
        for (jlast in seq((ib+1)*h,bigT,1)){
          for (jb in seq(ib*h,jlast-h,1)){
            dvec[jb,1] = optssr[jb,ib-1] + bigvec[jb*bigT - jb*(jb-1)/2 + jlast -jb,1]  
          }
        optssr[jlast,ib] = min(dvec[seq(ib*h,jlast-h,1),1])
        optdat[jlast,ib] = ib*h-1 + which.min(dvec[seq(ib*h,jlast-h,1),1])
        }
      }
      
      datevec[ib,ib] = optdat[bigT,ib]
      
      for (i in seq(1,ib-1,1)){
        xx = ib-i
        datevec[xx,ib] = optdat[datevec[xx+1,ib],xx]
      }
      glb[ib,1] = optssr[bigT,ib]
    }
  }
  
  out = list('glb' = glb, 'datevec' = datevec, 'bigvec' = bigvec)  
  return (out)  
}

