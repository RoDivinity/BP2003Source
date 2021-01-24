#function to generate m x m identity matrix
eye = function(m){
  ###
  # m = dimension of matrix
  # ***
  # eye = m x m identity matrix 
  ###
  
  eye = matrix (0L , nrow = m, ncol = m)
  for (i in 1:m){
    eye[i,i] = 1
  }
  
  return(eye)
}
