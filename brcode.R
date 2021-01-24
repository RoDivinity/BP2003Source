###################  R version, January 2021 ################
# These routines can be used by and distributed for non-profit academic
# purposes without any royalty except that the users must cite:
# Bai, Jushan and Pierre Perron (1998): "Estimating and Testing Linear
# Models with Multiple Structural Changes," Econometrica, vol 66,@47-78
# and
#  Bai, Jushan and Pierre Perron (2003): "Computation and Analysis of
#  Multiple Structural Change Models," Journal of Applied Econometrics, 18, 1-22.
# For any other commercial use and/or comments, please contact Pierre
# Perron at perron at bu.edu.
# Even though we tried to make this program error-free we cannot be held
# responsible for any consequences that could result from remaining errors.
# Copyright, Pierre Perron (1999, 2004) , 
# Matlab Code
##########################################################
rm(list=ls())
cat("\f")

#### Read data #######
data = read.csv('~/Desktop/BP2003/real.csv',header = FALSE)

######### Initial data setup #############
# y: dependent var
# z: T x q matrix of q regressors with coefficients are allowed to change
# x: T x p matrix of z regressors with coefficients are allowed to change
# bigT: effective sample size T
##########################################
y = matrix(data[,1])
bigT = as.numeric(dim(y)[1]) #set effective sample size
z = matrix(1,bigT,1)
#z = matrix(data[,4])
#z = data.matrix(z)
x = c()
#x = matrix(data[,2])
if (is.null(x)){print('no x regressors')} else {x = data.matrix(x)}


######### Basic configuration ##########
# m:    maximum number of structural changes allowed
# eps1: value of trimming in percentage for the construction 
#     and critical values of the supF type tests (used in the 
#     supF test, the Dmax, the supF(l+1|l) and the sequential 
#     procedure). If these tests are used, h below should be set 
#     at int(eps1*bigt). But if the tests are not required, 
#     estimation can be done with an arbitrary h. There are five 
#     options: eps1= .05, .10, .15, .20, or .25. For each option, 
#     the maximal value of m above is: 10 for eps=.05, 8 for 
#     eps1=.10, 5 for eps1=.15, 3 for eps1=.20, and 2 for eps1=.25.
# h:    minimal length of a segment (h>=q). Note: if robust=1, h 
#     should be set at a larger value. 
########
m = 5
q = as.numeric(dim(z)[2])
p = as.numeric(dim(x)[2])
if (is.null(x)){p = 0}
eps1 = 0.15
h = round(eps1*bigT) #minimal segment length

######### Partial structural change configuration ######
# fixb:     option to use fixed initial values of β
# betaini:  initial fixed values of β (Only available when fixb = 1)
# maxi:     maximum number of iterations for the nonlinear procedure
#       to obtain global minimizers
# printd:   option to print output from the iterations
# eps:      criterion for convergence
##########
fixb = 0
betaini = 0
maxi = 20
printd = 0
eps = 0.0001

######### Error term/Residual assumptions #######
# robust: allow for heterogeneity and autocorrelation 
#       in the residuals, 0 otherwise. The method used is Andrews(1991) 
#       automatic bandwidth with AR(1) approximation and the quadratic 
#       kernel. Note: Do not set to 1 if lagged dependent variables are 
#       included as regressors.
# prewhit: apply AR(1) prewhitening prior to estimating 
#       the long run covariance matrix.
# hetdat: allow different moment matrices of the regressors across segments. 
#       If hetdat=0, the same moment matrices are assumed for each segment 
#       and estimated from the ful sample. It is recommended to set 
#       hetdat=1 if p>0. (option for the construction of the F tests)
# hetvar: allow for the variance of the residuals to be different 
#       across segments. If hetvar=0, the variance of the residuals is assumed 
#       constant across segments and constructed from the ful sample. This 
#       option is not available when robust=1. (option for the construction
#       of the F tests) 
# hetomega:If hetomega=0, the long run covariance matrix of zu is assumed
#       identical across segments (the variance of the errors u if robust=0)
#       (option for the construction of the confidence interval of break dates)
# hetq:   If hetq=0, the moment matrix of the data is assumed identical 
#       across segments. (option for construction of the confidence
#       interval of break dates)
#################
prewhit = 1
robust = 1
hetdat = 1
hetvar = 1
hetomega = 1
hetq = 1


######### Procedure Configuration #########
# 1 = TRUE, 0 = FALSE
# doglobal: procedure to obtain global minimizers
# dotest:   procedure to construct the supF, UDmax and WDmax tests 
#           (doglobal must be set to 1 to run this procedure.)*
# dospflp1: procedure to construct the supF(l+1|l) tests where under
#         the null the l breaks are obtained using global minimizers. 
#         (doglobal must be set to 1 to run this procedure.)*
# doorder:  procedure that selects the number of breaks using 
#         information criteria. 
#         (doglobal must be set to 1 to run this procedure.)*
# dosequa:  procedure to estimate the breaks sequentially and estimate 
#         the number of breaks using supF(l+1|l) test
# dorepart: procedure to modify the break dates obtained from the 
#         sequential method using the repartition method of Bai (1995),
#         Estimating breaks one at a time. Required for the confidence 
#         intervals obtained with estim below to be valid.**
# (**): The estim procedure is in the last configuration section 
#################
doglobal = 1
dotest = 1
dospflp1 = 1
doorder = 1
dosequa = 1
dorepart = 1

############ Estimation Procedure Configuration ##############
# estimbic: procedure to estimate the model with the number of breaks 
#           selected by BIC.
# estimlwz: procedure to estimate the model with the number of breaks  
#           selected by LWZ
# estimseq: procedure to estimate the model with the number of breaks
#           selected using the sequential procedure
# estimrep: procedure to estimate the model with the breaks selected
#           using the repartition method
# estimfix: procedure to estimate the model with a prespecified number
#           of breaks equal to fixn set below
# fixn:     prespecified number of breaks
###########################################
estimbic = 1
estimlwz = 1
estimseq = 1
estimrep = 1
estimfix = 1
fixn = 2

##### main backend function ########
setwd('/Users/linhnguyen/Desktop/BP2003/')
source('/Users/linhnguyen/Desktop/BP2003/pbreak.R')

pbreak(bigT,y,z,q,m,h,eps1,robust,prewhit,hetomega,hetq,doglobal,dotest,dospflp1,doorder,dosequa,dorepart,estimbic,estimlwz,estimseq,estimrep,estimfix,fixb,x,p,eps,maxi,betaini,printd,hetdat,hetvar,fixn)






