#load required functions
source('/Users/linhnguyen/Desktop/BP2003/onebp.R')
source('/Users/linhnguyen/Desktop/BP2003/partione.R')
source('/Users/linhnguyen/Desktop/BP2003/date/ssr.R')
source('/Users/linhnguyen/Desktop/BP2003/rot90.R')

#function to construct modified sequential estimates using the
#repartition method of Bai (1995), Estimating Breaks one at a time
preparti = function(y,z,nbreak,dateseq,h,x,p) {
  bigT = dim(z)[1]
  q = dim(z)[2]
  
  #care if nbreak is matrix or scalar
  dv = matrix(0L, nrow = nbreak+2, ncol = 1)
  dv[1,1] = 0
  dv[seq(2,nbreak+1,1),1] = dateseq
  
  dv[nbreak+2,1] = bigT
  ds = matrix(0L,nrow = nbreak, ncol = 1)
  dr = matrix(0L,nrow = nbreak, ncol = 1)
  
  for (is in 1:nbreak){
    length = dv[is+2,1] - dv[is,1]
    if (p == 0){
      index = seq(dv[is,1]+1,dv[is+2,1],1)
      y_temp = y[index,1,drop=FALSE]
      z_temp = z[index,,drop=FALSE]
      vssr = ssr(1,y_temp,z_temp,h,length)
      y_temp_rev = rot90(rot90(y_temp))
      z_temp_rev = rot90(rot90(z_temp))
      vssrev = ssr(1,y_temp_rev,z_temp_rev,h,length)
      out = partione(h,length-h,length,vssr,vssrev)
      ds[is,1] = out$dx
      dr[is,1] = ds[is,1] + dv[is,1]
    }
    else{
      out = onebp(y,z,x,h,dv[is,1]+1,dv[is+2,1])
      ds[is,1] = out$bd
      dr[is,1] = ds[is,1]
    }
    
  }
  return(dr)
  
}