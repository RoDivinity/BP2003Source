#function to calculate simple OLS regression in matrix form
OLS = function(y,x) {
  ###
  # x = matrix of independent vars
  # y = vector of dependent vars
  # **********
  # b = beta of OLS
  ###
  #b = lm(formula = y~x)
  b1 = solve((t(x) %*% x)) %*% t(x) %*% y
  b1 = matrix(b1)
  return (b1)
}


# ### Test code  ####


