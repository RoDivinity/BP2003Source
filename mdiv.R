#function to calculate "right matrix division" in Matlab
mdiv = function(A,B){
  ###
  # A = left matrix 
  # B = right matrix
  ###
  # out = A/B
  ###
  
  ## 2 cases: 
  # if B is a scalar, element-wise operator
  # if B is a matrix, use the following steps:
  #  # 1) x = A\B, then A\B is LS solution to A*x = B
  # => x = inv(A'*A) * A'B
  # 2) A/B = (B'\A') = inv((B')'*B') * (B')' * A'
  if(length(class(B)) == 1) {return (A/B)}
  else{
  A_t = t(A)
  B_t = t(B)
  out = solve(t(B_t) %*% B_t) %*% t(B_t) %*% A_t
  out = matrix(out)
  }
  return (out)  
}
