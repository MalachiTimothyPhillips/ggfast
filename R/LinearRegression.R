#'
#' Point estimate for \beta
#' @param X Design matrix
#' @param Y Response vector
#' @example
#' betah(matrix(c(1,2,3,4,5,6), nr=2, byrow=TRUE),
#'       matrix(c(1,2,3), nc=3))
#'
betah<-function(X,Y){
  solve(t(X)%*%X)%*%t(X)%*%Y
}

#'
#' Prepare the data matrix by appending a column of 1's top it
#' @param df Data frame, minus the response
#' @examples
#' append(matrix(c(1,2,3,4,5,6), nr=2, byrow=TRUE))
#'
append<-function(df){
  X<-matrix(NA_real_, nr=dim(df)[1], nc=dim(df)[2]+1, byrow=TRUE)
  X[,1] = 1.0
  for(i in 1 : dim(df)[2]){
    for(j in 1 : dim(df)[1]){
      X[j,i+1] = df[j,i]
    }
  }
  list(X=X)
}



