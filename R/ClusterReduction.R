#'
#' Reduce the amount of data through using kmeans clustering algorithm
#' Call the in R kmean cluster function in order to reduce a pool of data into a few, clustered points
#'
#' @param x Data in the form of a data frame or data matrix, for safety, please use data.matrix(...) as a wrapper
#' @param k Number of clusters to form
#' @param fileName Name of file you wish to save a .csv formatted file to
#' @param iter.max Number of iterations to use at a maximum (default = 10)
#' @param nstart How many random sets to be chosen
#' @param algorithm What algorithm to use (default="Lloyd").
#'                  choices include: "Hartigan-Wong", "LLoyd", "MacQueen", and "Forgy"
#' @param trace logical or integer number, currently only used in the default method ("Hartigan-Wong"):
#'  if positive (or true), tracing information on the progress of the algorithm is produced.
#'  Higher values may produce more tracing information.
#'
#'  @examples
#'  reduce_cluster(data.matrix(iris), 10, "myFile.csv")
#'  reduce_cluster(df, 15, "myManyIterationsFile.csv", iter.max=1000, algorithm="MacQueen")
#'
reduce_cluster<-function(x, k, fileName, iter.max = 10, nstart = 1, algorithm="Lloyd", trace = FALSE)
{
  # perform cluster analysis to find means of each cluster group
  cl<-kmeans(x,k,iter.max=iter.max,nstart=nstart,algorithm=algorithm,trace=trace)
  reduced_data<-cl$centers

  # write output data in the form of a .csv file
  write.csv(reduced_data, fileName)

}

#'
#' Visualize the amount of data reduction in the first component from using clustering
#' This verifies that the clustering analysis can be of some use for data reduction
#'
#' @param x Data in the form of a data frame or data matrix, for safety, please use data.matrix(...) as a wrapper
#' @param k Number of clusters to form
#'
#' @examples
#' visualize_reduction(data.matrix(iris), 19)
visualize_reduction<-function(x,k)
{
  cl<-kmeans(x,k)
  reduced_data<-cl$centers
  windows()
  plot(x)
  points(reduced_data,col="Blue", add=TRUE)
}

#'
#' Perform the types of anaylsis which do a reduction on the data set,
#' finding the optimum value of k
#'
#' NOTE: by default, this will perform the analysis on the scaled data
#'
#' @param x Data in the form of a data frame or data matrix, for safety, please use data.matrix(...) as a wrapper
#' @param kmax (=10) Maximum number of clusters to form
#' @param elbowPlotFileName (="elbowPlot.jpeg"), name of file for elbow plot output
#' @param silhouetteFileName (="silhouette.jpeg"), name of file for silhouette plot output
#' @param gapFileName (="gap.jpeg"), name of file for output of gap statistics plot
#' @param iter.max (=100) Maximum number of iterations
#' @param nstart How many random sets to be chosen
#' @param algorithm What algorithm to use (default="Lloyd").
#'                  choices include: "Hartigan-Wong", "LLoyd", "MacQueen", and "Forgy"
#' @param trace logical or integer number, currently only used in the default method ("Hartigan-Wong"):
#'  if positive (or true), tracing information on the progress of the algorithm is produced.
#'  Higher values may produce more tracing information.
#' @param debugOutput (=FALSE) Boolean value for whether or not to output silhouette plots for
#' each and every cluster value used. This generates a TON of output! YOU HAVE BEEN WARNED!
#' @param debugOutFileNameBase (="sil") Base file name for debug output, if enabled.
#'
#' @examples
#' cluster_k_analysis(data.matrix(iris), kmax=19)
#' cluster_k_analysis(data.matrix(iris), kmax=19, iter.max=200, debugOutput=T)
cluster_k_analysis<-function(x,
                             kmax=10,
                             elbowPlotFileName = "elbowPlot.jpeg",
                             silhouetteFileName = "silhouette.jpeg",
                             gapFileName = "gap.jpeg",
                             iter.max = 100,
                             nstart = 1,
                             algorithm="Lloyd",
                             trace = FALSE,
                             debugOutput=FALSE,
                             debugOutFileNameBase="sil")
{
  # normalize the data
  x<-scale(x)
  wss<-function(k){
    kmeans(x,k,iter.max,nstart,algorithm,trace)$tot.withinss
  }

  # Compute from k = 1 to k = n
  n = dim(x)[1]-1
  k.values <- 2 : kmax

  # extract wss for 2-n clusters
  wss_values<-purrr::map_dbl(k.values,wss)

  # Output the 'elbow plot' approach to determining the ideal
  # number of clusters to use in this type of analysis
  jpeg(elbowPlotFileName)

  # make a plot
  plot(k.values, wss_values,
       type="b", pch = 19, frame = FALSE,
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")

  dev.off()

  # function to compute average silhouette for k clusters
  avg_sil <- function(k) {
    km.res <- kmeans(x, centers = k, nstart = 25)
    ss <- cluster::silhouette(km.res$cluster, dist(x))
    if(debugOutput){
      name=paste(debugOutFileNameBase, k, ".jpeg", sep="")
      jpeg(name)
      plot(ss)
      dev.off()
    }
    mean(ss[,3])
  }

  avg_sil_values <- purrr::map(k.values, avg_sil)

  jpeg(silhouetteFileName)
  plot(k.values, avg_sil_values,
       type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters K",
       ylab = "Average Silhouettes")
  dev.off()

  gap_stat<-cluster::clusGap(x,FUN=kmeans, nstart=nstart, K.max=kmax, B=50)
  jpeg(gapFileName)
  plot(gap_stat, xlab="Number of clusters k")
  dev.off()

  list(GapStatistics = gap_stat)
}

