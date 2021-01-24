#load required functions to run all procedures
source('/Users/linhnguyen/Desktop/BP2003/date/nldat.R')
source('/Users/linhnguyen/Desktop/BP2003/pftest.R')
source('/Users/linhnguyen/Desktop/BP2003/getcv1.R')
source('/Users/linhnguyen/Desktop/BP2003/getcv2.R')
source('/Users/linhnguyen/Desktop/BP2003/getdmax.R')
source('/Users/linhnguyen/Desktop/BP2003/spflp1.R')
source('/Users/linhnguyen/Desktop/BP2003/ssrnull.R')
source('/Users/linhnguyen/Desktop/BP2003/order.R')
source('/Users/linhnguyen/Desktop/BP2003/sequa.R')
source('/Users/linhnguyen/Desktop/BP2003/preparti.R')
source('/Users/linhnguyen/Desktop/BP2003/estim.R')
#main function to carry out configurations from BPcode setup
pbreak = function(bigT,y,z,q,m,h,eps1,robust,prewhit,hetomega,hetq,doglobal,dotest,dospflp1,doorder,dosequa,dorepart,estimbic,estimlwz,estimseq,estimrep,estimfix,fixb,x,p,eps,maxi,betaini,printd,hetdat,hetvar,fixn){

if(TRUE){
  #print initial configurations
  sep = '***********************************************'
  sep_s = '------------------------------------------------'
  print('The option chosen are: ')
  print(paste('h =',h))
  print(paste('eps1 =',eps1))
  print(paste('hetdat =',hetdat))
  print(paste('hetvar =',hetvar))
  print(paste('hetomega =',hetomega))
  print(paste('hetq =',hetq))
  print(paste('hetvar =',hetvar))
  print(paste('robust =',robust,'prewhit =',prewhit))
  print(paste('Maximum number of breaks is:',m))
  print(sep)
}
  
#significance level
siglev=matrix(c(10,5,2.5,1),4,1);

if(doglobal == 1){
  #procedure to obtain break dates and the associated SSR for all 
  #numbers of breaks between 1 to m
  print('Output from the global optimization')
  print(sep)
  
  if(p == 0) {
    print('This is a pure structural change model with the following specifications:')
    print(paste(q,'regressors z with allowed to change coefficients'))
    print(paste('maximum number of breaks:',m))
    out = dating(y,z,h,m,q,bigT)
  }
  else{
  print('This is a partial structural change model with the following specifications:')
    print(paste(q,'regressors z with allowed to change coefficients'))
    print(paste(p,'regressors x with fixed coefficients'))
    print(paste('maximum number of breaks:',m))
    print(paste('initial values of β option (1)=TRUE/(0)=FALSE',fixb))
    if(fixb == 1) {print('initial values β')
      print(betaini)}
    print(paste('convergence criterion:',eps))
    print(paste('print iteration option (1)=TRUE/(0)=FALSE',printd))
  out = nldat(y,z,x,h,m,p,q,bigT,fixb,eps,maxi,betaini,printd)
  }
}

#store all relevant information about model
glb = out$glb
datevec = out$datevec
bigvec = out$bigvec

#printing results
for (i in 1:m){
  print(paste('Model with',i,'breaks has SSR:',glb[i,1]))
  print('The dates of breaks are:')
  print(datevec[1:i,i])
}


if (dotest == 1) {
  print(sep)
  print('Output from testing procedure')
  
  #procedure for F test
  print(sep_s)
  print('a) supF tests against a fixed number of breaks')
  print(sep_s)
  
  ftest = matrix(0L, nrow = m, ncol = 1)
  wftest = matrix(0L, nrow = m, ncol = 1)
  
  for (i in 1:m){
    ftest[i,1] = pftest(y,z,i,q,bigT,datevec,prewhit,robust,x,p,hetdat,hetvar)
    print(paste('supF test for 0 versus',i,'breaks (scaled by q):',ftest[i,1]))
  }
  print(sep_s)
  
  for (c in 1:4){
    #critical values for supF test
    cv = getcv1(c,eps1)
    print(paste('The critical values at the',siglev[c,1],'% level are (for k = 1 to',m,'):'))
    print(cv[q,1:m,drop=FALSE])
  }
  
  #procedure for Dmax and UDmax test
  print(sep_s)
  print('b) Dmax test against an unknown number of breaks')
  print(sep_s)
  print(paste('The UDmax test is:',max(ftest)))
  
  for (c in 1:4) {
    #critical values for Dmax test
    cvm = getdmax(c,eps1)
    print(paste('The critical values at the',siglev[c,1],'% level is:',
                cvm[q,1]))
  }
  print(sep)
  
  for (c in 1:4){
    #computation of WDmax test
    cv = getcv1(c,eps1)
    for( i in 1:m){
      wftest[i,1] = cv[q,1] * ftest[i,1] / cv[q,1]
    }
    print(sep_s)
    print(paste('WDmax test at the',siglev[c,1],'% level is:',max(wftest)))
  }
  
}

if (dospflp1 == 1) {
  # the procedure to calculate supF(l+1|l) tests where
  # first l breaks are taken from the global minimization
  print(sep)
  print('supF(l+1|l) tests using global optimizers under the null')
  print(sep_s)
  
  for (i in seq(1,m-1,1)){
    out = spflp1(bigvec,datevec[1:i,i,drop=FALSE],i+1,y,z,h,q,prewhit,robust,x,p,hetdat,hetvar)
    supfl = out$maxf
    ndat = out$newd
    print(paste('The supF(',i+1,'|',i,') test is',supfl))
    print(paste('It corresponds to a new break at:',ndat))
  }
  
  for (c in 1:4){
    #critical values for supF(l+1|l) test
    cv = getcv2(c,eps1)
    print(paste('The critical values at the',siglev[c,1],'% level are (for k = 1 to',m,'):'))
    print(cv[q,1:m,drop=FALSE])
  }
}

if (doorder == 1) {
  #procedure to estimate number of breaks using information criteria
  #BIC and LWZ (Liu, Wu & Zidek). It returns the orders selected
  print(sep)
  print('Output from the application of Information Criteria')
  print(sep_s)
  
  if (p == 0){zz = z}
  else{zz = cbind(z,x)}
  
  ssr0 = ssrnull(y,zz)
  out = order(ssr0,glb,bigT,m,q)
  #store result
  mbic = out$mBIC
  mlwz = out$mLWZ
  print(paste('The number of breaks chosen by BIC is:',mbic))
  print(paste('The number of breaks chosen by LWZ is:',mlwz))
}

if (dosequa == 1) {
  #the procedure to sequentially estimate each break one at a time.
  #It stops when supF(l+1|l) test is not significant
  #It returns the number of breaks found and the break dates. 
  #Note that it can be used independently of the other procedures, 
  #i.e. global minimizers need not be obtained since 
  #it used a method to compute the breaks in O(T) operations.
  nbreak = matrix(0L, nrow = 4, ncol = 1)
  dateseq = matrix(0L,nrow = 4, ncol = m)
  
  for (j in 1:4){
    print(sep)
    print(paste('Output from the sequential procedure at significance level',
                siglev[j,1],'%'))
    print(sep_s)
    
    out = sequa(m,j,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)
    nbr = out$nbreak
    datese = out$dv0
    print(sep_s)
    print(paste('The sequential procedure estimated the number of breaks at:',nbr))
    if (nbr > 0) {
      print('The break dates are:')
      print(datese)
    }
    nbreak[j,1] =nbr
    
    if (nbr!=0){
      dateseq[j,seq(1,nbreak[j,1])] = t(datese)
    }
    
  }
  
}

if (dorepart == 1){
  # The following procedure constructs the so-called repartition
  # estimates of the breaks obtained by the sequential method (see Bai
  # (1995), Estimating Breaks one at a time, Econometric Theory, 13,
  # 315-352. It alows estimates that have the same asymptotic
  # distribution as those obtained by global minimization. Otherwise, the
  # output from the procedure "estim" below do not deliver asymptotically
  # correct confidence intervals for the break dates.
  
  reparv = matrix (0L,4,m)
  
  for (j in 1:4){
    print(sep)
    print(paste('Output from the repartition procedure for the',
                siglev[j,1],'% significance level'))
    if (nbreak[j,1] == 0){
      print(sep)
      print(('The sequential procedure found no break and the repartition procedure is skipped.'))
      print(sep)
    }
    else {
      repartda = preparti(y,z,nbreak[j,1,drop=FALSE],
                          t(dateseq[j,seq(1:nbreak[j,1]),drop=FALSE]),
                          h,x,p)
      reparv[j,seq(1:nbreak[j,1])] = repartda
      }
    }
    
  }
  
if(estimbic == 1){
  print(sep)
  print('Output from the estimation of the model selected by BIC')
  print(sep_s)
  estim(mbic,q,z,y,datevec[,mbic,drop=FALSE],robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
}

if(estimlwz == 1){
  print(sep)
  print('Output from the estimation of the model selected by LWZ')
  print(sep_s)
  estim(mlwz,q,z,y,datevec[,mlwz,drop=FALSE],robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
}

if(estimseq == 1){
  print(sep)
  ii = 0
  j = 1
  while(j<=4){
    if (ii == 0){
      if (nbreak[j,1] != 0){
        print(paste('Output from the estimation of the model selected by the sequential method at significance level',
                    siglev[j,1], '%'))
        print(sep_s)
        
        estim(nbreak[j,1],q,z,y,t(dateseq[j,seq(1,nbreak[j,1]),drop=FALSE]),robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
      }
    }
    
    j = j+1
    
    if(j <= 4){
      if(nbreak[j,1] == nbreak[j-1,1]){
        if (nbreak[j,1] == 0){
          print(paste('For the',siglev[j,1], '% level, the model is the same as for the',
                      siglev[j-1,1], '% level.'))
          print('The estimation is not repeated')
          print(sep_s)
          ii = 1
        }
        else{
          if (identical(dateseq[j,seq(1,nbreak[j,1]),drop = FALSE],dateseq[j-1,seq(1,nbreak[j-1,1]),drop = FALSE])){
            print(paste('For the',siglev[j,1], '% level, the model is the same as for the',
                        siglev[j-1,1], '% level.'))
            print('The estimation is not repeated')
            print(sep_s)
            ii = 1
          }
        }
      }
    }
    else{ii = 0}
  }
}

if(estimrep == 1){
  ii = 0
  print(sep)
  
  while (j <= 4){
    if (ii==0){
    if (nbreak[j,1] == 0){
      print(sep)
      print(paste('The sequential procedure at the significance level',
                  siglev[j,1],'% found no break and'))
      print('the repartition procedure was skipped')
      print(sep)
    }
    else{
      print(paste('Output from the estimation of the model selected by the repartition method from the sequential procedure at the significance level',
                  siglev[j,1],'%'))
      print(sep_s)
      
      estim(nbreak[j,1],q,z,y,t(reparv[j,seq(1,nbreak[j,1]),drop=FALSE]),
            robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
      }
    }
    
    j = j+1
    if (j <= 4) {
      if (nbreak[j,1] == nbreak[j-1,1]){
        if (nbreak[j,1] == 0){
          print(paste('For the',siglev[j,1], '% level, the model is the same as for the',
                      siglev[j-1,1], '% level.'))
          print('The estimation is not repeated')
          print(sep_s)
          ii = 1
        }
        else {
          if (identical(dateseq[j,seq(1,nbreak[j,1]),drop = FALSE],dateseq[j-1,seq(1,nbreak[j-1,1]),drop = FALSE])){
            print(paste('For the',siglev[j,1], '% level, the model is the same as for the',
                        siglev[j-1,1], '% level.'))
            print('The estimation is not repeated')
            print(sep_s)
            ii = 1
          }
        }
      }
      else {
        ii = 0
      }
    
  }
  }
}

if(estimfix == 1){
  print(sep)
  print(paste('Output from the estimation of the model with',fixn,'breaks'))
  print(sep_s)
  estim(fixn,q,z,y,datevec[,fixn,drop=FALSE],robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
}


}



