#function to get critical values of supF(l+1|l) test
getcv2 = function(signif,eps1){
  ###
  # eps1: trimming level
  # signif: significant level
  # the critical values of the test is stored in /Asset/supF_next/cv.csv
  ###
  # cv: critical values
  ####
  
  if(eps1 == .05){
    out = read.csv('~/Desktop/BP2003/Asset/supF_next/cv_1.csv',header = FALSE)
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }
  
  if(eps1 == .10) {
    out = read.csv('~/Desktop/BP2003/Asset/supF_next/cv_2.csv',header = FALSE)
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }
  
  if(eps1 == .15) {
    out = read.csv('~/Desktop/BP2003/Asset/supF_next/cv_3.csv',header = FALSE)
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }
  
  if(eps1 == .20) {
    out = read.csv('~/Desktop/BP2003/Asset/supF_next/cv_4.csv',header = FALSE)
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }
  
  if(eps1 == .25) {
    out = read.csv('~/Desktop/BP2003/Asset/supF_next/cv_5.csv',header = FALSE)
    cv = out[seq((signif-1)*10+1,signif*10,1),]
  }
  cv = data.matrix(cv)
  colnames(cv) = NULL
  rownames(cv) = NULL
  return(cv)
}