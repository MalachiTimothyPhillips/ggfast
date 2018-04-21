append_design_matrix <-
function(df){
  X<-matrix(NA_real_, nr=dim(df)[1], nc=dim(df)[2]+1, byrow=TRUE)
  X[,1] = 1.0
  for(i in 1 : dim(df)[2]){
    for(j in 1 : dim(df)[1]){
      X[j,i+1] = df[j,i]
    }
  }
  list(X=X)
}
