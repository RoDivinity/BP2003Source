#diagonal partition given a break date 
diag_par = function(input,m,date) {
  ###
  # input = matrix of independent vars z (no change) to diagonally partite
  # m = number of breaks in the series
  # date = vector of date for partition
  # ******************
  # output = partitioned dependent vars correspond to # breaks & break dates
  ###
  nt = dim(input)[1] 
  q1 = dim(input)[2]
  #create output matrix of zeros with size: nt x (break+1)*q1
  output = matrix(0L,nrow = nt, ncol = (m+1)*q1)
  #copy values of 1st segment of z to the output matrix
  output[c(1:date[1,1]),c(1:q1)] = input[c(1:date[1,1]),,drop=FALSE]
  i = 2
  while (i <= m){
    #copy values of i-th segment of z to output matrix corresponding to date vector
    r_i = seq(date[i-1,1]+1,date[i,1],1) #indices of rows to copy input values
    c_i = seq((i-1)*q1+1,i*q1,1)     #indices of cols to copy input values
    output[r_i,c_i]=input[r_i,,drop=FALSE]
    i = i+1
  }
  rl_i = seq(date[m,1]+1,nt,1)
  c_i = seq(m*q1+1,(m+1)*q1,1)
  output[rl_i,c_i] = input[rl_i,,drop=FALSE]
  return (output)
}


# ### TEST CODE ### IT'S WORKING
# tMatrix = matrix(c(1:15),nrow = 3, ncol = 5)
# tDate = 2
# tM = 1
# print(tMatrix)
# #print(tMatrix[c(2:5),c(4:6)])
# test = diag_par(tMatrix,tM,tDate)
# print(test)

