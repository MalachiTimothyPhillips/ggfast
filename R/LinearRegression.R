#'
#' Point estimate for beta, used in linear regression
#'
#' Notice that this is done by solving the linear system:
#' (X'X) beta = X'Y
#'
#' Rather than using a direct method, ie.
#'
#' beta = (X'X)^(-1)X'Y
#'
#' We shall instead use a package known as 'pcig' that solves
#' a pre-conditioned conjugate gradient formulation, and is more
#' numerically stable whenever the matrix approaches singularity.
#'
#' @param X Design matrix
#' @param Y Response vector
#' @examples
#' betah(matrix(c(1,2,3,4,5,6), nr=3, byrow=TRUE),
#'       matrix(c(1,2,3), nr=3))
#'
betah<-function(X,Y){
  if(!require("pcg")) install.packages("pcg")
  pcg::pcg(t(X)%*%X, t(X)%*%Y)
}

#'
#' Prepare the data matrix by appending a column of 1's top it
#' @param df Data frame, minus the response
#' @examples
#' append_design_matrix(matrix(c(1,2,3,4,5,6), nr=2, byrow=TRUE))
#'
#' df<-iris
#' df<-df[,2:dim(df)[2]] # remove first column, using it as the response vector
#' append_design_matrix(df)
#'
append_design_matrix<-function(df){

  # strip off non numeric data
  df<-df[,sapply(df, is.numeric)]
  X<-matrix(NA_real_, nr=dim(df)[1], nc=dim(df)[2]+1, byrow=TRUE)
  X[,1] = 1.0
  for(i in 1 : dim(df)[2]){
    for(j in 1 : dim(df)[1]){
      X[j,i+1] = df[j,i]
    }
  }
  return (X)
}

#'
#' Perform boot strapping
#' @param Y response vector
#' @param X design matrix
#' @param iter(=1) number of iterations to run bootstrap
#' @param alpha (=0.05) significance level for generating confidence intervals
#' @examples
#' df<-iris
#'
#' # scrape off all non-numeric data first
#' df<-df[,sapply(df, is.numeric)]
#' df
#' X<-df[,-1]
#' Y<-as.matrix(df[,1])
#' X
#' Y
#' dim(X)
#' dim(Y)
#' X<-append_design_matrix(X)
#' dim(X)
#' bh=bootreg(Y,X,1000) # works just fine! Whoo!
#' bh
bootreg=function(Y,X,iter=1, alpha = 0.05)
{
  n=dim(X)[1]
  # n is the number of rows
  nb=dim(X)[2]
  # nb is the number of columns and is equal to number of betas to be estimated
  beta.mat=matrix(NA_real_,nrow = nb,ncol=iter)
  index.mat=matrix(NA_real_,nrow=n,ncol=iter,byrow=FALSE)
  for(i in 1:iter)
  {
    index.mat[,i] = sample(1:n,n,replace=TRUE)
    #index.mat will be a matrix of random indices
    X.iter<-matrix(NA_real_, nrow=n, ncol = nb, byrow=TRUE)
    Y.iter<-matrix(NA_real_, nrow=n, ncol = 1, byrow = TRUE)
    for(j in 1 : n){
      X.iter[j,] = X[index.mat[j,i],]
      Y.iter[j,1] = Y[index.mat[j,i]]
    }
    beta.mat[,i]= betah(X.iter, Y.iter)
    #beta.mat is a matrix filled with bootstrap beta estimates
  }
  layout(matrix(1:nb, nr=nb,nc=1,byrow=TRUE))
  CIs = c();
  for(i in 1:nb)
  {
    j=i-1
    CI = quantile(beta.mat[i,], c(alpha, 1-alpha))
    CIs = c(CIs, CI)
    hist(beta.mat[i,],freq=FALSE,
         ylab="Density"  , main=substitute(beta[j]),
         xlab = substitute(beta[j]))
    mtext(paste("CI = (", round(CI[1],4), ",",round(CI[2],4),")" ), cex=0.6)
  }


  list(bh=beta.mat, Y=Y,X=X, Confidence_Intervals = CIs)
}

# Some other cools ideas: add in a predictor method

#'
#' Predict the next value using the same multilinear regression model
#'
#' @param X design matrix
#' @param model the result of running bootreg on the data set, ie. a great big collection
#' of different betas, upon which the column mean is used to get the most accurate possible
#' estimate for beta
#' @examples
#' # predict a singular value, given the description of the linear model and the design matrix of the new points
#' df<-iris
#'
#' # scrape off all non-numeric data first
#' df<-df[,sapply(df, is.numeric)]
#' df
#' X<-df[,-1]
#' Y<-as.matrix(df[,1])
#' X
#' Y
#' dim(X)
#' dim(Y)
#' X<-append_design_matrix(X)
#' dim(X)
#' bh=bootreg(Y,X,1000) # works just fine! Whoo!
#' data.matrix(rowMeans(bh$bh))
#' predictor(matrix(c(1, 2, 2, 3), nr=1), bh)
#'
predictor<-function(X, model){
  betas = model$bh

  # Use the point estimate for the mean of the betas
  beta_means<-data.matrix(rowMeans(betas))

  # Y = Xbeta
  X%*%beta_means

}

############################################################################################

"
NOTICE: This work is not originally my own, but rather is a culmination of
work that was inspired by a repository full of fantastic Bayesian inspired stuff
that is possible within R, for more details, see https://github.com/timothyfrasier/Bayesian-Stats-Course
"


