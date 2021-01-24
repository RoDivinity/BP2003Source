#function to obtain optimal one break partition for a segment
parti = function (start,b_start,b_end,last,bigvec,bigT) {
  ###
  # start: start date index of the segment
  # last: end date index of the segment
  # b_start: first possible break date
  # b_end: last possible breakdate
  # bigT : sample periods T
  #***************
  # ssrmin: return associated SSR of optimal break
  # dx:  optimal date (global minimizer)
  ###
  
  #initialize array to store betas and value of index for procedure
  dvec = matrix(0L , nrow = bigT, ncol = 1)
  ini = (start-1)*bigT - (start-2)*(start-1)/2 + 1
  j = b_start
  
  #start the loop on the segment to find optimal break
  while (j <= b_end){
    jj = j - start
    k = j*bigT - (j-1)*j/2+last-j
    dvec[j,1] = bigvec[ini+jj,1] + bigvec[k,1]
    j = j+1
  }
  #get min SSR
  ssrmin = min( dvec[seq(b_start,b_end,1),1] )
  #get the date with lowest SSR
  dx = (b_start - 1) + which.min(dvec[seq(b_start,b_end,1)]) 
  out = list('ssrmin' = ssrmin, 'dx' = dx)
  #print(paste('min SSR',ssrmin))
  return(out)
}


#### TEST CODE ######
#bigvec = matrix(c(1,0.5,0.2,2,1,6),ncol=1)
#out = parti(1,1,2,3,bigvec,3)
# print(out)
