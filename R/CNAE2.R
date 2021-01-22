#' @title CNAE dataset on classes 4 and 9
#'
#' @description The data set \code{CNAE2} is a subset of the original CNAE-9 data, that 
#'  comprises 1080 documents categorized into 9 topics of free text business 
#'  descriptions of Brazilian companies. 
#'  
#' Specifically, \code{CNAE2} contains only the documents belonging to topics "4" and "9".
#'  The data set is already pre-processed and provides the bag-of-words representation of
#'  the documents; the columns with null counts are removed leading to a matrix with 240 documents
#'  on a vocabulary with cardinality 357. This data set is highly sparse
#'  (98% of the matrix is filled with zeros).
#'  
#'  Class labels are stored in \link[deepMOU]{cl_CNAE}
#'
#' @docType data
#'
#' @usage data(CNAE2)
#'
#' @format A matrix for the bag-of-words representation of the CNAE2 dataset.
#'
#' @keywords datasets
#'
#'
#' @source \href{https://archive.ics.uci.edu/ml/datasets/CNAE-9}{Original CNAE9 dataset}
#'
#' @examples
#' x = data(CNAE2)
#' print(head(x))
"CNAE2"

