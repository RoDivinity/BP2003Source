#load required functions
source('/Users/linhnguyen/Desktop/BP2003/interval.R')
source('/Users/linhnguyen/Desktop/BP2003/pvdel.R')
source('/Users/linhnguyen/Desktop/BP2003/OLS.R')
source('/Users/linhnguyen/Desktop/BP2003/diag_par.R')

# Function to estimate the model by OLS given the obtained break dates
# It also computes and reports confidence intervals for the break dates
# The method used depends on specification for robust

estim = function(m,q,z,y,b,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar){
  ####
  # m = number of break
  # q = # of regressors z
  # z = matrix of regressor z
  # p = # of regressors x
  # x = matrix of regressor x
  # b = break date
  # prewhit,hetomega,hetq,hetdat,hetvar = options for assumptions on res
  # ***********
  # CI = list of Confidence Intervals for each corresponding break
  ########
  
  #print string
  sep = '----------------------------------------------'
  
  if (m == 0){print('There are no breaks in this model and estimation is skipped')}
  else{
    bigT = dim(z)[1]
    d = (m+1)*q + p
    vdel = matrix(0L,nrow = d, ncol = d)
    
    #construct zbar matrix. Diagonal partition of Z 
    #at the estimated break date
    zbar = diag_par(z,m,b)
    
    #estimation and printing
    if (p == 0){
      reg = zbar
    }
    else{
      reg = cbind(x,zbar)
    }
    
    #estimation of β and δ in pure/partial model
    #lm ...
    
    vdel = pvdel(y,z,m,q,bigT,b,prewhit,robust,x,p,1,hetdat,hetvar)
    
    print(sep)
    print('Corrected standard errors for the coefficients')
    print(sep)
    
    for (i in 1:d){
      print(paste('Corrected SE for coefficient',i,'is',sqrt(vdel[i,i])))
    }
    
    if (robust == 0 && hetdat == 1 && hetvar == 0){
      print('In this case robust=0, hetdat=1 and hetvar=0, the "corrected" are the same as that of the printout except for a different small sample correction.')
    }
    
    #confidence interval for break date
    bound = interval(y,z,zbar,b,q,m,robust,prewhit,hetomega,hetq,x,p)
    print(sep)
    print('Confidence interval for the break dates')
    print(sep)
    
    for (i in 1:m){
      print(paste('The 95% C.I for the',i,'th break is:',bound[i,1],' ',bound[i,2]))
      print(paste('The 90% C.I for the',i,'th break is:',bound[i,3],' ',bound[i,4]))
    }
  }
  
}