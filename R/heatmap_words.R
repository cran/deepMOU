#' @title Heatmap of word frequencies by cluster
#'
#' @description Displays the heatmap of the cluster frequency distributions of the most frequent terms sorted by the most informative ones.
#'
#' @param x                  Document-term matrix describing the frequency of terms that occur in a collection of documents. Rows correspond to documents in the collection and columns correspond to terms.
#' @param clusters           Integer vector of length of the number of cases, which indicates a clustering. The clusters have to be numbered from 1 to the number of clusters. 
#' @param n_words            Number of words to include in the heatmap (default is 50).
#' @param legend_position    Position of the legend (\code{"none"}, \code{"left"}, \code{"right"}, \code{"bottom"}, \code{"top"}, or two-element numeric vector as in \link[ggplot2]{theme}). Default is \code{"bottom"}.
#' @param font_size          Text size in pts (default is 12).
#' @param low_color          Base color for terms with no occurrence in a cluster (default is \code{"grey92"})
#' @param top_color          Base color for terms concentrated in a single cluster (default is \code{"red"})
#' @param main               A title for the plot. Default is \code{"Row frequencies of terms distribution"}.
#' @param xlabel             A title for the x-axis. Default is \code{NULL}.
#' @param ylabel             A title for the y-axis. Default is \code{NULL}.
#' @param legend_title       A title for the legend. Default is \code{"Entropy"}.
#' 
#' @return A graphical aid to describe the clusters according to the most informative words.
#'
#' @importFrom entropy entropy
#' @import ggplot2
#'
#'
#' @details Takes as input the bag-of-words matrix and returns a heatmap displaying the row
#' frequency distribution of terms according to the clusters. Words are sorted by entropy.
#'
#'
#' @examples
#' # Load the CNAE2 dataset
#' data("CNAE2")
#'
#' # Get document labels by clustering using mou_EM
#' mou_CNAE2 = mou_EM(x = CNAE2, k = 2)
#'
#' # Usage of the function
#' heatmap_words(x = mou_CNAE2$x, clusters = mou_CNAE2$clusters)
#'
#' @export

heatmap_words = function (x, clusters, n_words = 50, legend_position = "bottom", font_size=12,
                          low_color = "grey92", top_color = "red", main="Row frequencies of terms distribution",xlabel=NULL,ylabel=NULL,legend_title="Entropy") {
  numobs = dim(x)[1]
  p = dim(x)[2]
  if (is.null(colnames(x)))
    colnames(x) = paste("W", sep = "", 1:p)
  clust_label = levels(factor(clusters))
  k = length(clust_label)
  clusters_words_count = matrix(NA, p, k, dimnames = list(colnames(x),
                                                          clust_label))
  for (i in 1:k) {
    docs_cluster = x[clusters == clust_label[i], ]
    clusters_words_count[, i] = colSums(docs_cluster)
  }
  freq_treshold = sort(rowSums(clusters_words_count), decreasing = T)[n_words]

  row_selection = rowSums(clusters_words_count) >= freq_treshold
  clusters_words_count = clusters_words_count[row_selection,
                                              ]
  clusters_words_freq = clusters_words_count/rowSums(clusters_words_count)
  words_entropy = apply(clusters_words_freq, 1, entropy)
  clusters_words_freq = clusters_words_freq[order(words_entropy, decreasing = T), ]


  words_grid = expand.grid(Terms = rownames(clusters_words_freq),
                           Clusters = clust_label)
  words_grid$Freq = c(clusters_words_freq)
  ggplot(words_grid, aes(Clusters, Terms, fill = Freq)) + labs(subtitle =main,fill=legend_title) + xlab(xlabel)+ylab(ylabel)+
    geom_tile() + scale_fill_gradient(low = low_color, high = top_color) +
    theme(legend.position = legend_position, text = element_text(size=font_size))
}
