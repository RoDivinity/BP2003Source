#load required functions
source('/Users/linhnguyen/Desktop/BP2003/OLS.R')

#function to compute SSR under H0: no break
ssrnull = function(y,zz) {
  delta = OLS(y,zz)
  resid = y - zz %*% delta
  ssrn = t(resid) %*% resid
  return(ssrn)
}