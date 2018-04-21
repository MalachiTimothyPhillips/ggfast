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
#' @param iter number of iterations to run bootstrap
#' @param alpha significance level for generating confidence intervals
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

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = floor( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}
# openGraphSaveGraph.R
# John K. Kruschke, 2013

openGraph = function( width=7 , height=7 , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    X11( width=width , height=height , type="cairo" , ... )
  } else { # Windows OS
    windows( width=width , height=height , ... )
  }
}

saveGraph = function( file="saveGraphOutput" , type="pdf" , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    if ( any( type == c("png","jpeg","jpg","tiff","bmp")) ) {
      sptype = type
      if ( type == "jpg" ) { sptype = "jpeg" }
      savePlot( file=paste(file,".",type,sep="") , type=sptype , ... )
    }
    if ( type == "pdf" ) {
      dev.copy2pdf(file=paste(file,".",type,sep="") , ... )
    }
    if ( type == "eps" ) {
      dev.copy2eps(file=paste(file,".",type,sep="") , ... )
    }
  } else { # Windows OS
    savePlot( file=file , type=type , ... )
  }
}

plotPost = function( paramSampleVec , credMass=0.95 , compVal=NULL ,
                     HDItextPlace=0.7 , ROPE=NULL , yaxt=NULL , ylab=NULL ,
                     xlab=NULL , cex.lab=NULL , cex=NULL , xlim=NULL , main=NULL ,
                     col=NULL , border=NULL , showMode=F , showCurve=F , breaks=NULL ,
                     ... ) {
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlab) ) xlab="Parameter"
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(xlim) ) xlim=range( c( compVal , paramSampleVec ) )
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="skyblue"
  if ( is.null(border) ) border="white"

  postSummary = matrix( NA , nrow=1 , ncol=11 ,
                        dimnames=list( c( xlab ) ,
                                       c("mean","median","mode",
                                         "hdiMass","hdiLow","hdiHigh",
                                         "compVal","pcGTcompVal",
                                         "ROPElow","ROPEhigh","pcInROPE")))
  postSummary[,"mean"] = mean(paramSampleVec)
  postSummary[,"median"] = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  postSummary[,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]

  #("HDIofMCMC.R")
  HDI = HDIofMCMC( paramSampleVec , credMass )
  postSummary[,"hdiMass"]=credMass
  postSummary[,"hdiLow"]=HDI[1]
  postSummary[,"hdiHigh"]=HDI[2]

  # Plot histogram.
  if ( is.null(breaks) ) {
    breaks = c( seq( from=min(paramSampleVec) , to=max(paramSampleVec) ,
                     by=(HDI[2]-HDI[1])/18 ) , max(paramSampleVec) )
  }
  if ( !showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                     freq=F , border=border , col=col ,
                     xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                     breaks=breaks , ... )
  }
  if ( showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , plot=F )
    densCurve = density( paramSampleVec , adjust=2 )
    plot( densCurve$x , densCurve$y , type="l" , lwd=5 , col=col , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
          main=main , cex=cex , cex.lab=cex.lab , ... )
  }
  cenTendHt = 0.9*max(histinfo$density)
  cvHt = 0.7*max(histinfo$density)
  ROPEtextHt = 0.55*max(histinfo$density)
  # Display mean or mode:
  if ( showMode==F ) {
    meanParam = mean( paramSampleVec )
    text( meanParam , cenTendHt ,
          bquote(mean==.(signif(meanParam,3))) , adj=c(.5,0) , cex=cex )
  } else {
    dres = density( paramSampleVec )
    modeParam = dres$x[which.max(dres$y)]
    text( modeParam , cenTendHt ,
          bquote(mode==.(signif(modeParam,3))) , adj=c(.5,0) , cex=cex )
  }
  # Display the comparison value.
  if ( !is.null( compVal ) ) {
    cvCol = "darkgreen"
    pcgtCompVal = round( 100 * sum( paramSampleVec > compVal )
                         / length( paramSampleVec )  , 1 )
    pcltCompVal = 100 - pcgtCompVal
    lines( c(compVal,compVal) , c(0.96*cvHt,0) ,
           lty="dotted" , lwd=1 , col=cvCol )
    text( compVal , cvHt ,
          bquote( .(pcltCompVal)*"% < " *
                    .(signif(compVal,3)) * " < "*.(pcgtCompVal)*"%" ) ,
          adj=c(pcltCompVal/100,0) , cex=0.8*cex , col=cvCol )
    postSummary[,"compVal"] = compVal
    postSummary[,"pcGTcompVal"] = ( sum( paramSampleVec > compVal )
                                    / length( paramSampleVec ) )
  }
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    ropeCol = "darkred"
    pcInROPE = ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                 / length( paramSampleVec ) )
    lines( c(ROPE[1],ROPE[1]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol)
    text( mean(ROPE) , ROPEtextHt ,
          bquote( .(round(100*pcInROPE))*"% in ROPE" ) ,
          adj=c(.5,0) , cex=1 , col=ropeCol )

    postSummary[,"ROPElow"]=ROPE[1]
    postSummary[,"ROPEhigh"]=ROPE[2]
    postSummary[,"pcInROPE"]=pcInROPE
  }
  # Display the HDI.
  lines( HDI , c(0,0) , lwd=4 )
  text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
        adj=c(.5,-1.7) , cex=cex )
  text( HDI[1] , 0 , bquote(.(signif(HDI[1],3))) ,
        adj=c(HDItextPlace,-0.5) , cex=cex )
  text( HDI[2] , 0 , bquote(.(signif(HDI[2],3))) ,
        adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
  par(xpd=F)
  #
  return( postSummary )
}
BMLRmcmc = function( dataMat , numSavedSteps=250000 , checkConvergence=FALSE ) {
  # Program written in the style of
  require(rjags)     # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
  # A Tutorial with R and JAGS. Academic Press / Elsevier.
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  # JAGS model specification begins here...
  model {
  # Likelihood:
  for( i in 1:N ) {
  y[i] ~ dnorm( y.hat[i] , tau )
  y.hat[i] <- b0 +  inprod( b[1:nPred] , x[i,1:nPred] )
  }
  # Prior (assumes standardized data):
  tau <- 1/pow(sigma,2)
  sigma ~ dunif( 0 , 10 )
  b0 ~ dnorm( 0 , 1.0E-2 )
  for ( j in 1:nPred ) {
  b[j] ~ dnorm( 0 , 1.0E-2 )
  }
  }
  # ... end JAGS model specification
  " # close quote for modelstring
  writeLines(modelstring,con="model.txt")

  #------------------------------------------------------------------------------
  # THE DATA.

  # dataMat is supplied as an argument to this function.
  # The program assumes that the first column is y, and the remaining columns
  # are x (predictors).

  # Now re-name variables from dataMat:
  N = NROW(dataMat)
  y = cbind(as.matrix(dataMat[,1]))
  predictedName = colnames(dataMat)[1]
  x = as.matrix(dataMat[,-1])
  predictorNames = colnames(dataMat)[-1]
  nPred = NCOL(x)

  # Prepare data for JAGS:
  # Re-center data at mean, to reduce autocorrelation in MCMC sampling.
  # Divide by SD to make prior specification generic.
  standardizeCols = function( dataMat ) {
    zDataMat = dataMat
    for ( colIdx in 1:NCOL( dataMat ) ) {
      mCol = mean( dataMat[,colIdx] )
      sdCol = sd( dataMat[,colIdx] )
      zDataMat[,colIdx] = ( dataMat[,colIdx] - mCol ) / sdCol
    }
    return( zDataMat )
  }
  zx = standardizeCols( x )
  zy = standardizeCols( y )

  # Get the data into JAGS:
  dataList = list(
    x = zx ,
    y = as.vector( zy ) , # JAGS does not treat 1-column mat as vector
    N = N ,
    nPred = nPred
  )

  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.

  # Start the chains at the least-squares fit:
  lmInfo = lm( dataList$y ~ dataList$x )
  initsList = list(
    b0 = lmInfo$coef[1] ,
    b = lmInfo$coef[-1] ,
    sigma = sqrt(mean(lmInfo$resid^2))
  )

  #------------------------------------------------------------------------------
  # RUN THE CHAINS

  parameters = c("b0" , "b" , "sigma" )
  adaptSteps = 1000         # Number of steps to "tune" the samplers.
  burnInSteps = 100        # Number of steps to "burn-in" the samplers.
  nChains = 3               # Number of chains to run.
  numSavedSteps=250000       # Total number of steps in chains to save.
  thinSteps=1               # Number of steps to "thin" (1=keep every step).
  nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
  # Create, initialize, and adapt the model:
  jagsModel = jags.model( "model.txt" , data=dataList , inits=initsList ,
                          n.chains=nChains , n.adapt=adaptSteps )
  # Burn-in:
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel , n.iter=burnInSteps )
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples( jagsModel , variable.names=parameters ,
                              n.iter=nPerChain , thin=thinSteps )
  # resulting codaSamples object has these indices:
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

  #------------------------------------------------------------------------------
  # EXAMINE THE RESULTS

  if ( checkConvergence ) {
    openGraph(width=7,height=7)
    autocorr.plot( codaSamples[[1]] , ask=FALSE )
    show( gelman.diag( codaSamples ) )
    effectiveChainLength = effectiveSize( codaSamples )
    show( effectiveChainLength )
  }

  # Convert coda-object codaSamples to matrix object for easier handling.
  # But note that this concatenates the different chains into one long chain.
  # Result is mcmcChain[ stepIdx , paramIdx ]
  mcmcChain = as.matrix( codaSamples )
  chainLength = NROW(mcmcChain)

  # Rename results, for convenience:
  zsigma = mcmcChain[,"sigma" ]
  zb0 = mcmcChain[,"b0"]
  # changed lines below to adapt for single regression
  zb = matrix(mcmcChain[,"b"])
  #   zb = matrix( 0 , nrow=chainLength , ncol=nPred )
  #   for ( j in 1:nPred ) {
  #     zb[,j] = mcmcChain[,paste("b[",j,"]",sep="")]
  #   }
  colnames(zb) = "zb1"
  # end of changed lines

  # Convert to original scale (see book Eqn 17.1):
  sigma = zsigma * sd(as.vector(y))
  # for b0:
  subtractTerm = rep(0,chainLength)
  for ( j in 1:nPred ) {
    subtractTerm = subtractTerm + zb[,j] * sd(as.vector(y)) * mean(x[,j]) / sd(x[,j])
  }
  b0 = zb0 * sd(as.vector(y)) + mean(y) - subtractTerm
  # for b:
  b = 0*zb
  for ( j in 1:nPred ) {
    b[,j] = zb[,j] * sd(as.vector(y)) / sd(x[,j])
  }
  colnames(b) = paste("b",1:nPred,sep="")

  # Compute R^2 for credible parameters:
  Rsq = rep(0,chainLength)
  YcorX = cor( y , x ) # correlation of y with each x predictor
  for ( stepIdx in 1:chainLength ) {
    Rsq[stepIdx] = sum( zb[stepIdx,] * YcorX )
  }
  ## Old version, perhaps useless:
  #Rsq = rep(0,chainLength)
  #for ( stepIdx in 1:chainLength ) {
  #  predY = dataList$x %*% cbind(zb[stepIdx,]) + zb0[stepIdx]
  #  Rsq[stepIdx] = cor(dataList$y,predY)^2
  #}

  # Combine results into one big matrix:
  mcmcChain = cbind( zsigma,zb0,zb, sigma,b0,b, Rsq )
  return( mcmcChain )
}


BMLRplot = function( mcmcChain , ROPEbeta=NULL , ROPEbetaDiff=NULL ) {

  summaryMatrix = matrix( NA , nrow=0 , ncol=11 ,
                          dimnames=list(ParameterName=NULL,PosteriorInfo=NULL) )

  # Compute number of predictors from NCOL of mcmcChain matrix.
  # Assumes mcmcChain has columns zsigma,zb0,zb,sigma,b0,b,Rsq
  nPred = (NCOL(mcmcChain)-1)/2 - 2

  # Display scatter plots of parameter values, pairwise:
  openGraph(width=7,height=7)
  chainLength = NROW(mcmcChain)
  nToPlot = 700
  plotIdx = seq(1,chainLength,by=ceiling(chainLength/nToPlot))
  pairs( mcmcChain[plotIdx,c("Rsq","zsigma","zb0",paste("zb",1:nPred,sep=""))] ,
         col="skyblue" )
  openGraph(width=7,height=7)
  pairs( mcmcChain[plotIdx,c("Rsq","sigma","b0",paste("b",1:nPred,sep=""))] ,
         col="skyblue" )

  #   # Show correlation matrix on console:
  #   cat("\nCorrlations of posterior sigma, b0, and b's:\n")
  #   show( round( cor( mcmcChain ) , 3) )
  #   cat("\nCovariances of posterior sigma, b0, and b's:\n")
  #   show( round( cov( mcmcChain ) , 5) )

  # Display the marginals of posterior:
  nPostCol = 3 # arbitrary number of columns for display
  nPostRow = ceiling((3+nPred)/nPostCol)
  openGraph(width=nPostCol*2.5,height=nPostRow*2.0)
  layout(matrix(1:(nPostCol*nPostRow),nrow=nPostRow,byrow=TRUE))
  par( mar=c(4,2,3,1) , mgp=c(2.5,0.7,0) )
  histInfo = plotPost( mcmcChain[,"Rsq"] , xlab=bquote(R^2) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  summaryMatrix = rbind( summaryMatrix , histInfo )
  histInfo = plotPost( mcmcChain[,"zsigma"] , xlab=bquote(sigma[z]) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  summaryMatrix = rbind( summaryMatrix , histInfo )
  histInfo = plotPost( mcmcChain[,"zb0"] , xlab=bquote(beta[0]) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  summaryMatrix = rbind( summaryMatrix , histInfo )
  for ( j in 1:nPred ) {
    histInfo = plotPost( mcmcChain[,paste("zb",j,sep="")] ,
                         xlab=bquote(beta[.(j)]) , main="" ,
                         compVal = 0.0 , ROPE=ROPEbeta ,
                         cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
    summaryMatrix = rbind( summaryMatrix , histInfo )
  }
  #
  nPostCol = 3 # arbitrary number of columns for display
  nPostRow = ceiling((3+nPred)/nPostCol)
  openGraph(width=nPostCol*2.5,height=nPostRow*2.0)
  layout(matrix(1:(nPostCol*nPostRow),nrow=nPostRow,byrow=TRUE))
  par( mar=c(4,2,3,1) , mgp=c(2.5,0.7,0) )
  histInfo = plotPost( mcmcChain[,"Rsq"] , xlab=bquote(R^2) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  #summaryMatrix = rbind( summaryMatrix , histInfo )
  histInfo = plotPost( mcmcChain[,"sigma"] , xlab=bquote(sigma[y]) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  summaryMatrix = rbind( summaryMatrix , histInfo )
  histInfo = plotPost( mcmcChain[,"b0"] , xlab=bquote(b[0]) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  summaryMatrix = rbind( summaryMatrix , histInfo )
  for ( j in 1:nPred ) {
    histInfo = plotPost( mcmcChain[,paste("b",j,sep="")] ,
                         xlab=bquote(b[.(j)]) , main="" ,
                         compVal = 0.0 ,
                         cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
    summaryMatrix = rbind( summaryMatrix , histInfo )
  }

  # Display differences of standardized regression coefficients:
  nPostCol = 3
  nPostCell = nPred*(nPred-1)/2
  nPostRow = ceiling(nPostCell/nPostCol)
  openGraph(width=nPostCol*2.5,height=nPostRow*2.0)
  layout(matrix(1:(nPostCol*nPostRow),nrow=nPostRow,byrow=TRUE))
  par( mar=c(4,2,3,1) , mgp=c(2.5,0.7,0) )
  for ( j in 1:(nPred-1) ) {
    for ( k in (j+1):nPred ) {
      histInfo = plotPost( mcmcChain[,paste("zb",j,sep="")]
                           - mcmcChain[,paste("zb",k,sep="")],
                           xlab=bquote(beta[.(j)] - beta[.(k)]) , main="" ,
                           compVal = 0.0 , ROPE=ROPEbetaDiff ,
                           cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
      summaryMatrix = rbind( summaryMatrix , histInfo )
    }
  }

  return( summaryMatrix )

} # end of function BMLRplot


n     <- 20                  # length of vector
rho   <- 0.55                 # desired correlation = cos(angle)
theta <- acos(rho)             # corresponding angle
x1    <- rnorm(n, 1, 1)        # fixed given data
x2    <- rnorm(n, 1, 1)      # new random data
X     <- cbind(x1, x2)         # matrix
Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)

Id   <- diag(n)                               # identity matrix
Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
cor(x1, x)
df <- data.frame(x1, x)
plot(df)
cor(x1, x)
summary(lm(df[,2] ~ df[,1]))


graphics.off()          # Clears all graphical displays.
# rm(list=ls(all=TRUE))   # Removes all variables from memory!

# Get the Bayesian functions into R's working memory:
#source("BMLR.R")

# Load the data into R:
# dataMat = read.csv( file="BMLRexampleData.csv" )
# Important: The matrix dataMat must have the criterion y in its first column,
# and the predictors in the subsequent columns!

# Run the Bayesian analysis and put the results in variable named mcmcChain:
# dataMat <- dataMat[,1:2]
dataMat <- df # use fake data to explore and compare lm output
mcmcChain = BMLRmcmc( dataMat )

# Plot the posterior distribution and put summary in variable named postInfo:
postInfo = BMLRplot( mcmcChain )
# Display the summary on the console:
# show(postInfo)

# check
summary(lm(dataMat[,1] ~ dataMat[,2]))


# # Another example, using a large N data set:
# dataMat = read.csv( file="BMLRexampleLargeNdata.csv" )
# mcmcChain = BMLRmcmc( dataMat , numSavedSteps=100000 )
# postInfo = BMLRplot( mcmcChain ,
#                      ROPEbeta=c(-0.05,0.05) , ROPEbetaDiff=c(-0.05,0.05) )
# show(postInfo)
