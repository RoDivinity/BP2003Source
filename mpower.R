#function to calculate power of a matrix
mpower = function(A,t){
  ####
  # A = input matrix
  # t = exponent
  # ****
  # out = A^t
  ####
  B = A
  if (t == 1) {return(B)}
  else{for (i in 2:t) {B = B %*% A}}
  return (B)
}
