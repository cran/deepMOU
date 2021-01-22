#' @title Deep Mixture of Unigrams
#'
#' @description Performs parameter estimation by means of Gibbs sampling and cluster allocation
#' for the Deep Mixture of Unigrams.
#'
#' @param x            Document-term matrix describing the frequency of terms that occur in a collection of documents. Rows correspond to documents in the collection and columns correspond to terms.
#' @param k            Number of clusters/groups at the top layer.
#' @param g            Number of clusters at the bottom layer.
#' @param n_it         Number of Gibbs steps.
#' @param seed_choice  Set seed for reproducible results.
#' @param burn_in      Number of initial Gibbs samples to be discarded and not included in the computation of final estimates.
#' @return A list containing the following elements:
#' \item{x}{The data matrix.}
#' \item{clusters}{the clustering labels.}
#' \item{k}{the number of clusters at the top layer.}
#' \item{g}{the number of clusters at the bottom layer.}
#' \item{numobs}{the sample size.}
#' \item{p}{the vocabulary size.}
#' \item{z1}{the allocation variables at the top layer.}
#' \item{z2}{the allocation variables at the bottom layer.}
#' \item{Alpha}{the estimates of Alpha parameters.}
#' \item{Beta}{the estimates of the Beta parameters.}
#' \item{pi_hat}{estimated probabilities of belonging to the \code{k} clusters at the top layer conditional to the \code{g} clusters at the bottom layer.}
#' \item{pi_hat_2}{estimated probabilities of belonging to the \code{g} clusters at the bottom layer.}
#' 
#'
#' @importFrom Rfast Lgamma
#' @importFrom stats runif rmultinom
#' @importFrom extraDistr rdirichlet rdirmnom
#'
#' @details Starting from the data matrix \code{x}, the Deep Mixture of Unigrams is fitted
#' and \code{k} clusters are obtained.
#' The algorithm for the estimation of the parameters is the Gibbs sampling.
#' In particular, the function assigns initial values to all the parameters to be estimated. Then \code{n_it} samples for the parameters are obtained using
#' conditional distributions on all the other parameters. The final estimates are obtained by averaging the samples given that initial \code{burn_in} samples are
#' discarded. Clustering is eventually performed by maximizing the posterior distribution of the latent variables.
#' For further details see the references.
#'
#' @references
#' Viroli Cinzia, Anderlucci Laura (2020). "Deep mixtures of Unigrams for uncovering topics in textual data". \emph{Statistics and Computing}, forthcoming, 1-18.
#'
#' @examples
#' # Load the CNAE2 dataset
#' data("CNAE2")
#'
#' # Perform parameter estimation and clustering, very few iterations used for this example
#' deep_CNAE2 = deep_mou_gibbs(x = CNAE2, k = 2, g = 2, n_it = 5, burn_in = 2)
#'
#' # Shows cluster labels to documents
#' deep_CNAE2$clusters
#' @export
deep_mou_gibbs = function (x, k, g, n_it = 500, seed_choice = 1, burn_in = 200)
{
  set.seed(seed_choice)
  progress_bar = txtProgressBar(min = 0, max = n_it, style = 3)
  setTxtProgressBar(progress_bar, 0)
  if (is.null(colnames(x)))
    colnames(x) = paste("W", sep = "", 1:dim(x)[2])
  max_x = max(colSums(x))
  numobs = nrow(x)
  p = ncol(x)
  lik = NULL
  Alpha_dep = Beta_dep = pi_hat_1_dep = pi_hat_2_dep = 0
  delta = 1
  z1 = matrix(0, numobs, k)
  z2 = matrix(0, numobs, g)
  gamma_prev = array(0, c(numobs, k, g))
  Beta = matrix(runif(k * p, 0, 20), p, k)
  Alpha = matrix(runif(g * p, -0.1, 0.1), p, g)
  if (g == 1)
    Alpha = matrix(0, p, g)
  pi_hat_2 = rep(1/g, g)
  pi_hat_1 = matrix(rep(1/k, k), k, g)
  edu_guess = mou_EM(x, k)
  z1[cbind(seq_along(edu_guess$cl), edu_guess$cl)] = 1
  for (i in 1:k) Beta[, i] = colSums(x[edu_guess$cl == i, ])
  Beta = ifelse(Beta == 0, 0.001, Beta)
  Beta = ifelse(Beta > max_x, max_x, Beta)
  for (iter in 1:n_it) {
    setTxtProgressBar(progress_bar, iter)
    for (i in 1:k) for (j in 1:g) {
      beta_alf_prod = Beta[, i] * (Alpha[, j] + 1)
      mat_beta_alf = matrix(beta_alf_prod, numobs, p, byrow = TRUE)
      den = rowSums(Lgamma(mat_beta_alf) - Lgamma(x + mat_beta_alf))
      num = Lgamma(rep(sum(beta_alf_prod), numobs)) - Lgamma(rowSums(x +
                                                                       mat_beta_alf))
      gamma_prev[, i, j] = num - den
    }
    gammaa = gamma_prev - max(gamma_prev) + 700
    gammaa = exp(gammaa)
    if (g > 1) {
      prob = matrix(0, numobs, g)
      for (j in 1:g) {
        prob[, j] = rowSums((matrix(pi_hat_1[, j, drop = FALSE] *
                                      pi_hat_2[j], numobs, k, byrow = T) * matrix(gammaa[,
                                                                                           , j], numobs, k)) * z1)
      }
      prob = prob/rowSums(prob)
      prob = ifelse(is.na(prob), 1/g, prob)
      z2 = t(apply(prob, 1, function(x) rmultinom(n = 1,
                                                  size = 1, prob = x)))
      if (colSums(z2)[j] == 0) {
        index = sample(1:numobs, 1)
        z2[index, j] = 1
      }
      n_j = colSums(z2)
      pi_hat_2 = rdirichlet(1, n_j + delta)
    }
    for (j in 1:g) {
      n_i = colSums(z1 * matrix(z2[, j], numobs, k))
      pi_hat_1[, j] = rdirichlet(1, n_i + delta)
    }
    if (g > 1) {
      Beta_n = t(Beta[, apply(z1, 1, which.max)])
      for (j in 1:g) {
        Beta_ni = Beta_n[z2[, j] == 1, , drop = FALSE]
        xx = x[z2[, j] == 1, , drop = FALSE]
        xx = ifelse(xx == 0, 0.001, xx)
        n_j = nrow(xx)
        size = (Beta_ni %*% (1 + Alpha[, j]))
        prop = rdirmnom(length(size), round(size), xx)
        prop = ifelse(is.na(prop), 0.001, prop)
        prop = ifelse(prop == 0, 0.001, prop)
        prop = prop/Beta_ni - 1
        prop = colMeans(prop)
        prop = ifelse(prop > 1, 1, prop)
        Alpha[, j] = prop
      }
    }
    if (g > 1)
      Alpha_n = t(Alpha[, apply(z2, 1, which.max)])
    else Alpha_n = matrix(0, numobs, p)
    for (i in 1:k) {
      Alpha_ni = Alpha_n[z1[, i] == 1, , drop = FALSE]
      xx = x[z1[, i] == 1, , drop = FALSE]
      xx = ifelse(xx == 0, 0.001, xx)
      size = Beta[, i] %*% t(1 + Alpha_ni)
      prop = rdirmnom(length(size), round(size), xx)
      prop = ifelse(is.na(prop), 0.001, prop)
      prop = ifelse(prop == 0, 0.001, prop)
      prop = prop/(1 + Alpha_ni)
      prop = colMeans(prop)
      prop = ifelse(is.na(prop), 0.001, prop)
      prop = ifelse(prop > max_x, max_x, prop)
      Beta[, i] = prop
    }
    if (k > 1) {
      prob = matrix(0, numobs, k)
      for (i in 1:k) {
        prob[, i] = rowSums((matrix(pi_hat_1[i, ],
                                    numobs, g, byrow = TRUE) * gammaa[, i, ]) *
                              z2)
      }
      prob = prob/rowSums(prob)
      prob = ifelse(is.na(prob), 1/k, prob)
      z1 = t(apply(prob, 1, function(x) rmultinom(1, 1,
                                                  x)))
    }
    if (iter > burn_in) {
      Beta_dep = (Beta_dep * (iter - burn_in - 1) + Beta)/(iter -
                                                             burn_in)
      Alpha_dep = (Alpha_dep * (iter - burn_in - 1) + Alpha)/(iter -
                                                                burn_in)
      pi_hat_1_dep = (pi_hat_1_dep * (iter - burn_in -
                                            1) + pi_hat_1)/(iter - burn_in)
      pi_hat_2_dep = (pi_hat_2_dep * (iter - burn_in -
                                            1) + pi_hat_2)/(iter - burn_in)
    }
  }
  Beta = Beta_dep
  Alpha = Alpha_dep
  pi_hat_1 = pi_hat_1_dep
  pi_hat_2 = pi_hat_2_dep
  for (i in 1:k) for (j in 1:g) {
    den = rowSums(lgamma(matrix(Beta[, i] * (Alpha[, j] +
                                               1), numobs, p, byrow = TRUE)) - lgamma(x + matrix(Beta[,
                                                                                                      i] * (Alpha[, j] + 1), numobs, p, byrow = TRUE)))
    num = lgamma(rep(sum(Beta[, i] * (Alpha[, j] + 1)), numobs)) -
      lgamma(rowSums(x + matrix(Beta[, i] * (Alpha[, j] +
                                               1), numobs, p, byrow = TRUE)))
    gamma_prev[, i, j] = num - den
  }
  gammaa = gamma_prev - max(gamma_prev) + 700
  gammaa = exp(gammaa)
  if (g > 1) {
    prob = matrix(0, numobs, g)
    for (j in 1:g) prob[, j] = rowSums((matrix(pi_hat_1[,
                                                          j] * pi_hat_2[j], numobs, k, byrow = TRUE) * gammaa[,
                                                                                                                , j]) * z1)
    prob = prob/rowSums(prob)
    prob = ifelse(is.na(prob), 1/g, prob)
    z2 = t(apply(prob, 1, function(x) rmultinom(1, 1, x)))
    if (colSums(z2)[j] == 0) {
      index = sample(1:numobs, 1)
      z2[index, j] = 1
    }
    n_j = colSums(z2)
  }
  prob = matrix(0, numobs, k)
  for (i in 1:k) prob[, i] = rowSums((matrix(pi_hat_1[i,
                                                        ], numobs, g, byrow = TRUE) * gammaa[, i, ]) * z2)
  prob = prob/rowSums(prob)
  prob = ifelse(is.na(prob), 1/k, prob)
  z1 = t(apply(prob, 1, function(x) rmultinom(1, 1, x)))
  clusters = apply(z1, 1, which.max)
  output = list(x = x, k = k, g = g, numobs = numobs, p = p,
                clusters = clusters, z1 = z1, z2 = z2, Alpha = Alpha_dep,
                Beta = Beta_dep, pi_hat = pi_hat_1_dep, pi_hat_2 = pi_hat_2_dep)
  close(progress_bar)
  class(output) ="deepMOU"
  return(output)
}
