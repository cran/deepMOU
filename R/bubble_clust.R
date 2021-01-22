#' @title Bubble plot
#'
#' @description Bi-dimensional representation (via Multi-Dimensional Scaling) of the clusters, where each bubble corresponds to a cluster,
#' its size is proportional to the relative frequency of documents and color saturation reflects cluster cohesion.
#'
#' @param x                  Document-term matrix describing the frequency of terms that occur in a collection of documents. Rows correspond to documents in the collection and columns correspond to terms.
#' @param clusters           Integer vector of length of the number of cases, which indicates a clustering. The clusters have to be numbered from 1 to the number of clusters. 
#' @param bubble_size        Graphical parameter for bubbles size.
#' @param bubble_col         Choose palette with two colors (default \code{"red"} and \code{"white"}).  The first (darker) color will denote homogeneous clusters, while the latter (lighter) more heterogeneous ones.
#' @param cex_text           Size of texts inside bubbles.
#' @param main               A title for the plot.
#' @return A graphical aid to visualize and to describe the obtained clusters.
#'
#' @importFrom MASS isoMDS
#' @importFrom skmeans skmeans_xdist
#' @importFrom grDevices colorRamp palette rgb
#' @importFrom graphics par symbols
#' @importFrom stats cmdscale runif
#'
#' @examples
#' # Load the CNAE2 dataset
#' data("CNAE2")
#'
#' # Perform parameter estimation and clustering
#' mou_CNAE2 = mou_EM(x = CNAE2, k = 2)
#'
#' # Usage of the function
#' bubble_clust(mou_CNAE2$x,mou_CNAE2$clusters, bubble_size = 5 )
#'
#' @export

bubble_clust = function (x, clusters, bubble_size = 1, bubble_col = c("red", "white"), 
                        cex_text = 1, main=NULL) 
{
  oldpar = par(no.readonly = TRUE)   
  on.exit(par(oldpar)) 
  par(mar = c(0.5, 0.5, 2, 0.5))
  
  k = max(clusters)
  p = ncol(x)
  clusters = factor(clusters, levels = c(1:k))
  
  cos_dissimilarity_x = skmeans_xdist(x)
  average_dist = rep(0, k)
  
  for (i in 1:k) {
    D_i = cos_dissimilarity_x[clusters == i, clusters == i]
    if (is.matrix(D_i)) D_i = D_i[lower.tri(D_i)]
    average_dist[i] = mean(D_i)^2
  }
  
  average_dist = ifelse(is.na(average_dist), 0, average_dist)
  color_fill = average_dist
  radius_bubbles = sqrt(table(clusters)/pi)
  ramp = colorRamp(bubble_col)
  
  clust_means = matrix(0, k, p)
  for (i in 1:k) clust_means[i, ] = colMeans(x[clusters == i, ])
  
  D_means = skmeans_xdist(clust_means)
  D_means = ifelse(is.na(D_means), 0, D_means)
  
  if (k == 1) {
    fit = matrix(0, 1, 2)
    op = palette(rgb(ramp(c(color_fill, color_fill + 0.2)), maxColorValue = 255))
    fit = list(points = fit)
  }
  
  if (k == 2) {
    fit = try(cmdscale(D_means, k = 1), silent = TRUE)
    if (is.character(fit))  fit = cbind(seq(-0.3, 0.3, length = k), runif(k, -0.3, 0.3))
    if (!is.character(fit)) fit = cbind(fit, seq(-0.3, 0.3, length = k))
    if (dim(fit)[2] == 1)   fit = cbind(fit, runif(k, -0.3, 0.3))
    fit = list(points = fit)
  }
  
  if (k > 2) {
    fit = try(isoMDS(D_means, k = 2,trace=FALSE), silent = TRUE)
    if (is.character(fit)) fit = cbind(seq(-0.3, 0.3, length = k), runif(k, -0.3, 0.3))
  }
  
  if (k > 1)  op = palette(rgb(ramp(color_fill), maxColorValue = 255))
  x_margins = c(min(fit$points[, 1]) - 0.03 * bubble_size, 
                max(fit$points[, 1]) + 0.03 * bubble_size)
  y_margins = c(min(fit$points[, 2]) - 0.03 * bubble_size, 
                max(fit$points[, 2]) + 0.03 * bubble_size)
  symbols(fit$points[, 1], fit$points[, 2], circles = radius_bubbles, 
          inches = 0.1 * bubble_size, fg = "white", bg = 1:k, xlab = "", 
          ylab = "", ylim = y_margins, xlim = x_margins, xaxt = "n", 
          yaxt = "n",main=main)
  palette(op)
  vocabulary = colnames(x)
  words_to_plot = matrix("", k, 2)
  for (i in 1:k) {
    freq_words = colSums(x[clusters == i, , drop = F])
    most_important = order(freq_words, decreasing = T)[1:2]
    selected_words = vocabulary[most_important]
    words_to_plot[i, ] = paste(selected_words)
    if (radius_bubbles[i] > 0) {
      text(fit$points[i, 1], fit$points[i, 2], words_to_plot[i, 
                                                             1], cex = cex_text, adj = 0.1, pos = 3, offset = 0)
      text(fit$points[i, 1], fit$points[i, 2], words_to_plot[i, 
                                                             2], cex = cex_text, adj = 0.1, pos = 1, offset = 0.3)
    }
  }
  text(fit$points[, 1], fit$points[, 2], 1:k, cex = 0.8, adj = 1, 
       pos = 3, col = "darkred", offset = 0.7)
}
