#'
#' After calling whatever means of cluster reductions, now do a PCA
#' on the already reduced data set.
#'
#' @param x Data set to be analyzed
#' @param k Number of components considerd
#' @param barPlotFileName (="barPlotComp.jpeg") file name used to store bar plot output
#' @examples
#' PCAPlot(data.matrix(iris), 2)
#'
PCAPlot<-function(x, k, barPlotFileName="barPlotComp.jpeg")
{
  unscaled<-x # keep a local copy for later
  x<-scale(x)
  pca<-determine_greatest_contributors(x,k)
  PC<-x%*% pca$EigenVectors
  Component<-colnames(x)
  Proportion<-pca$Proportions
  Errors<-1-pca$Proportions
  df<-data.frame(Component,Proportion,Errors)

  theme<-ggplot2::theme_set(cowplot::theme_cowplot()) + ggplot2::theme(legend.position="none")
  g<-ggpubr::ggbarplot(df, x="Component", y="Proportion", fill="Component",palette="jco",ggtheme=theme, sort.val="desc")
  p<-ggpubr::ggbarplot(df, x="Component", y="Errors", fill="Component", palette="jco",ggtheme=theme, sort.val="desc")
  graph<-cowplot::plot_grid(g,p)
  cowplot::save_plot(barPlotFileName, graph, base_aspect_ratio=3.0)

  # Now, in the case of k = 1,2, or 3, we can visualize the data and apply SLR on it
  data<-cbind(Component, Proportion)
  data<-data[order(-Proportion),]

  # Take the top k parts only
  names<-data[(1:k),1]
  #print(data)
  print(names)
  reduced_scaled_data_set<-subset(x, select=names) # Additionally have the option of outputting the data

  reduced_unscaled_data_set<-subset(x, select=names)
  df<-reshape2::melt(reduced_unscaled_data_set)

  if(k == 2){
    windows()
    plot(reduced_unscaled_data_set,
         main="Reduced Data Set Scatter Plot",
         xlab=colnames(reduced_unscaled_data_set)[1],
         ylab=colnames(reduced_unscaled_data_set)[2])

    # After performing this analysis, we can also do a cluster analysis
    # Typically, k = 5 clusters is an appropriate enough choice for now
    k = 5
    clus<-kmeans(reduced_unscaled_data_set,centers=k)
    windows()
    cluster::clusplot(reduced_unscaled_data_set, clus$cluster)
  }

  # Additional plot -- look at the clusplot for all of the data
  k = 5
  clus<-kmeans(x, centers=k)
  windows()
  cluster::clusplot(x, clus$cluster)

  # We may now extend this into three dimensions instead
  kdf<-kmeans(x,4)
  newdf<-data.frame(x, K=kdf$cluster)
  pcdf<-princomp(x,cor=T,score=T)
  summary(pcdf)
  rgl::plot3d(pcdf$scores, col=newdf$K)

}

#'
#' Determine greatest data factors
#' @param x Data set to be analyzes
#' @param k Number of components considered
#'
#' @examples
#' determine_greatest_contributors(data.matrix(iris), 2)
#'
determine_greatest_contributors<-function(x, k)
{
  x<-scale(x)
  S<-cov(x)
  eig_res<-eigen(S)

  lambda_sum = sum(eig_res$value)
  n = length(eig_res$value)
  prop<-matrix(NA, nr=n)
  for(i in 1 : n){
    prop[i]=eig_res$value[i]/lambda_sum
  }

  err<-1-sum(prop[1:k])

  list(EigenValues = eig_res$value, EigenvalueSum= lambda_sum, Proportions=prop, Error=err, EigenVectors=eig_res$vectors)
}

