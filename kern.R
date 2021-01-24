#load required functions
source('mdiv.R')

#procedure to evaluate quadratic kernel at x
kern = function(x) {
  #####
  # x = value of some x
  ###
  # ker = quadratic kernel
  ####
  del=6*pi*x/5;
  ker=3*(sin(del)/del-cos(del))/(del*del);
  return (ker)
}
