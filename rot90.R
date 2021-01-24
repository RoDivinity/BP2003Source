#function to rotate a matrix 90 degree counterclockwise
rot90 = function(A){
  rA = dim(A)[1]
  cA = dim(A)[2]
  B = matrix (0L, nrow = cA, ncol = rA)
  for (i in 1:rA){
    tempA = A[i,,drop=FALSE]
    tempA = tempA[,seq(dim(A)[2],1,-1)] #reverse the row elements
    tempB = t(tempA) #transpose the row 
    B[,i] = tempB
  }
  return(B)
}
