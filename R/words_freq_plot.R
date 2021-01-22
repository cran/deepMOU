#' @title Graph of most frequent words of each cluster
#'
#' @description Graphical plot of the most frequent words of each cluster
#'
#' @param x               Document-term matrix describing the frequency of terms that occur in a collection of documents. Rows correspond to documents in the collection and columns correspond to terms.
#' @param clusters        Integer vector of length of the number of cases, which indicates a clustering. The clusters have to be numbered from 1 to the number of clusters. 
#' @param clust_label     Vector of length of the number of cluster containing the cluster names to be displayed (by default \code{"Cluster_1"}, \code{"Cluster_2"}, ...).
#' @param n_words         Number of words to display.
#' @param words_size      A numerical value giving the amount by which plotting words should be magnified with respect to the default setting.
#' @param axis_size       Magnification to be used for axis annotation with respect to the default setting.
#' @param set_colors      Choose palette for word colors.
#' @param main            A title for the plot. Default is \code{"Most frequent words for each cluster"}.
#' @param xlabel             A title for the x-axis. Default is empty.
#' @param ylabel             A title for the y-axis. Default is empty.
#' @return A graphical aid for visualizing the most frequent terms for each cluster.
#'
#' @details The number of most frequent words to be shown can be set by
#' \code{n_words} and also clusters names can be passed beforehand as a character
#' vector to \code{clust_label}
#'
#' @examples
#' # Load the CNAE2 dataset
#' data("CNAE2")
#'
#' # Perform parameter estimation and clustering
#' mou_CNAE2 = mou_EM(x = CNAE2, k = 2)
#'
#' # Usage of the function
#' words_freq_plot(mou_CNAE2$x, mou_CNAE2$clusters,n_words = 4, words_size = 2, main = "Example" )
#'
#' @export

words_freq_plot = function (x, clusters, clust_label = NULL, n_words = 5,
                            words_size = 2, axis_size = 1, set_colors = NA, 
                            main="Most frequent words for each cluster",xlabel="",ylabel="")  {
  
  vocabulary = colnames(x)
  k = max(clusters)
  
  if (is.null(clust_label))    clust_label = paste("Cluster_", sep = "", 1:k)
  xx = x/matrix(rowSums(x),nrow(x),ncol(x),byrow=F)
  omega = sapply(by(xx,clusters,colMeans,simplify = F),c)
  
  
  indicator_order = apply(omega, 2, function(omega_) order(omega_, decreasing = T)[1:n_words])
  words = matrix(vocabulary[as.numeric(indicator_order)], n_words, k)
  
  if ( any(!is.na(set_colors))  &  length(set_colors) != k ){
    warning("The vector of colors must have length equal to the number of clusters")
  }
  if ( any(is.na(set_colors) )  ) set_colors = c("dodgerblue","grey62","red","green4","coral", "deepskyblue","deeppink3","mediumaquamarine","lightsalmon2","blue4","palegreen","blueviolet","orange","mediumorchid3",  "khaki2")
  
  
  text_position = expand.grid(1:n_words, 1:k)
  text_position = text_position[, c(2, 1)]
  plot(text_position, axes = FALSE, type = "n", xlab =xlabel,
       ylab =ylabel, ylim = c(0.5, n_words + 0.5), xlim = c(0.5,
                                                         k + 0.5),
       main=main)
  axis(1, seq(1, k, 1), cex.axis = axis_size, labels = clust_label)
  text_position[1:n_words, 1] = text_position[1:n_words, 1] +
    0.05
  text(text_position[, 1], text_position[, 2], as.vector(words),
       cex = words_size, col = rep(set_colors[1:k], each = n_words))

}
