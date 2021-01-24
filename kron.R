#function to compute Kronecker Tensor Product
kron = function(A,B){
  ###
  # A,B = m x n matrix A & p x q matrix B
  # ****
  # out = A (x) B
  ###
  rA = dim(A)[1]
  cA = dim(A)[2]
  rB = dim(B)[1]
  cB = dim(B)[2]
  
  out = matrix(0L, nrow = rA * rB, ncol = cA * cB)
  
  for (i in 1:rA){
    for (j in 1:cA) {
      r_index = seq((i-1)*rB+1,i*rB)
      c_index = seq((j-1)*cB+1,j*cB)
      out[r_index,c_index] = A[i,j] * B
    }
  }
  return(out)
}
