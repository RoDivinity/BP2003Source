#load required functions
source('/Users/linhnguyen/Desktop/BP2003/rot90.R')
source('/Users/linhnguyen/Desktop/BP2003/date/ssr.R')
source('/Users/linhnguyen/Desktop/BP2003/partione.R')
source('/Users/linhnguyen/Desktop/BP2003/getcv2.R')
source('/Users/linhnguyen/Desktop/BP2003/pftest.R')
source('/Users/linhnguyen/Desktop/BP2003/onebp.R')
#function to apply sequential procedure to obtain number of breaks and break
#dates. Current version only allows pure structural changes. This will be
#generalized
sequa = function(m,signif,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1){
  
  
  dv = matrix(0L, nrow = m+2, ncol = 1)
  dv2 = matrix(0L, nrow = m+2, ncol = 1)
  ftestv = matrix(0L, nrow = m+1,ncol = 1)
  
  cv = getcv2(signif,eps1)
  dv[1,1] = 0
  
  if (p == 0){
    y_rev = rot90(rot90(y))
    z_rev = rot90(rot90(z))
    vssrev = ssr(1,y_rev,z_rev,h,bigT)
    vssr = ssr(1,y,z,h,bigT)
    out = partione(h,bigT-h,bigT,vssr,vssrev)
    datx = out$dx
    ssrmin = out$ssrmin
  }
  else{
    out = onebp(y,z,x,h,1,bigT)
    datx = out$bd
    ssrmin = out$ssrind
  }
  
  dv[2,1] = datx
  
  ftest=pftest(y,z,1,q,bigT,dv[2,1,drop=FALSE],prewhit,robust,x,p,hetdat,hetvar)
  
  if (ftest < cv[q,1]) {
    nbreak = 0
    dv[2,1] = 0
    #dv0 = dv[seq(2,nbreak+1,1),1]
    nseg = 1
  }
  else{
    print(paste('First break found at:',datx))
    nbreak = 1
    nseg = 2
    dv[nseg+1,1] = bigT
  }
  
  while(nseg <= m){
    ds = matrix(0L,nseg+1,1)
    ftestv = matrix(0L,nseg+1,1)
    
    i_s = 1
    
    while(i_s <= nseg){
      length = dv[i_s+1,1] - dv[i_s,1]
      
      if(length >= 2*h){
        if(p==0){
          y_temp = y[seq(dv[i_s,1]+1,dv[i_s+1,1]),1,drop=FALSE]
          z_temp = z[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          vssr = ssr(1,y_temp,z_temp,h,length)
          y_temp_rev = rot90(rot90(y_temp))
          z_temp_rev = rot90(rot90(z_temp))
          vssrev = ssr(1,y_temp_rev,z_temp_rev,h,length)
          out = partione(h,length-h,length,vssr,vssrev)
          ds[i_s,1] = out$dx
          ftestv[i_s,1] = pftest(y_temp,z_temp,1,q,length,ds[i_s,1,drop=FALSE],prewhit,
                                 robust,0,p,hetdat,hetvar)
        }
        else{
          y_temp = y[seq(dv[i_s,1]+1,dv[i_s+1,1]),1,drop=FALSE]
          z_temp = z[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          x_temp = x[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          out = onebp(y,z,x,h,dv[i_s,1]+1,dv[is+1,1])
          ds[i_s,1] = out$bd - dv[i_s,1]
          ftestv[i_s,1] = pftest(y_temp,z_temp,1,q,length,ds[i_s,1],
                                 prewhit,robust,x_temp,p,hetdat,hetvar)
        }
      }
      else{
        ftestv[i_s,1] = 0.0
      }
      i_s = i_s + 1
    }
    
    maxf = max(ftestv[seq(1,nseg,1),1])
    
    if (maxf < cv[q,nseg]){
      #print(nbreak)
      #dv0 = dv[seq(2,nbreak+1,1),1]
    }
    else {
      newseg = which.max(ftestv[seq(1,nseg),1])
      print(paste('Next break is found at:',ds[newseg,1]+dv[newseg,1]))
      dv[nseg+2,1] = ds[newseg,1] + dv[newseg,1]
      nbreak = nbreak + 1
      #check this sort 
      dv2 = sort(dv[seq(2,nseg+2,1),1])
      dv2 = matrix(dv2, ncol = 1)
      dv[1,1] = 0
      dv[seq(2,nseg+2,1),1] = dv2
      
    }
    nseg = nseg + 1
  }
  
  print('The sequential procedure has reached the upper limit')
  if (nbreak < 1) {dv0 = c()}
  else{
  dv0 = dv[seq(2,nbreak+1,1),1]}
  out = list('nbreak' = nbreak, 'dv0' = dv0)
}