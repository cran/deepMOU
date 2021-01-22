#' @title Abstracts dataset
#'
#' @description Dataset composed by titles and abstracts of articles published in 2015 in five statistical journals.
#'
#' @details  These are the titles and abstracts of all the articles published in 2015 by the following journals:
#'   - Journal of American Statistical Association (JASA)
#'   - Journal of the Royal Statistical Society - Series B
#'   - Annals of Statistics
#'   - Biometrika
#'   - Statistical Science
#'   
#' The dataset comprises 379 articles with a vocabulary of 606 words already
#' pre-processed (stemmed, lemmatized, stopwords removal etc.); terms with entropy
#' less than 0.3 were discarded (rule-of-thumb threshold).
#'
#' @docType data
#'
#' @usage data(Abstracts)
#'
#' @format A matrix for the bag-of-words representation of the Abstracts dataset; terms are sorted in alphabetical order.
#'
#' @keywords datasets
#'
#' @examples
#' x = data(Abstracts)
#' print(head(x))
"Abstracts"
