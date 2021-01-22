#' @title Plotting method for "shallow" and deep mixtures of Unigrams and mixtures of Dirichlet-Multinomials
#'
#' @description Bi-dimensional representation (via Multi-Dimensional Scaling) of the clusters, where each bubble corresponds to a cluster,
#' its size is proportional to the relative frequency of documents and color saturation reflects cluster cohesion.
#'
#' @param x                Output from \code{mou_EM}, \code{dir_mult_GD} or \code{deep_mou_gibbs}.
#' @param bubble_size        Graphical parameter for bubbles size.
#' @param bubble_col         Choose palette with two colors (default \code{"red"} and \code{"white"}).  The first (darker) color will denote homogeneous clusters, while the latter (lighter) more heterogeneous ones.
#' @param cex_text           Size of texts inside bubbles.
#' @param main               A main title for the plot.
#' @param y                Parameter not used
#' @param ...              Parameter not used
#' @return A graphical aid to visualize and to describe the obtained clusters.

#' @details The default graphical representation of \code{mou_EM}, \code{dir_mult_GD} and \code{deep_mou_gibbs} is the bubble plot. 
#' Namely, a bi-dimensional representation (via Multi-Dimensional Scaling) of the clusters, each bubble corresponds to a cluster,
#' its size is proportional to the relative frequency of documents and color saturation reflects cluster cohesion.


#' @examples
#' # Load the CNAE2 dataset
#' data("CNAE2")
#'
#' # Perform parameter estimation and clustering
#' mou_CNAE2 = mou_EM(x = CNAE2, k = 2)
#'
#' # Usage of the function
#' plot(mou_CNAE2, bubble_size = 5 )
#'
#' @export

plot.deepMOU<-function(x, y, bubble_size = 1, bubble_col = c("red", "white"), 
                       cex_text = 1, main=NULL,...){
  bubble_clust(x$x, x$clusters, bubble_size = bubble_size, bubble_col = bubble_col, cex_text = cex_text,main=main )
}