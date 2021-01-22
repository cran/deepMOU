#' @title Mixture of Unigrams by Expectation-Maximization algorithm
#'
#' @description Performs parameter estimation by means of the Expectation-Maximization (EM) algorithm and cluster allocation
#' for the Mixture of Unigrams.
#'
#' @param x            Document-term matrix describing the frequency of terms that occur in a collection of documents. Rows correspond to documents in the collection and columns correspond to terms.
#' @param k            Number of clusters/groups.
#' @param n_it         Number of iterations for the Expectation-Maximization algorithm.
#' @param eps          Tolerance level for the convergence of the algorithm. Default is \code{1e-07}.
#' @param seed_choice  Set seed for reproducible results
#' @param init         Vector containing the initial document allocations for the initialization of the algorithm. If \proglang{NULL} (default) initialization is carried out via spherical k-means (\link[skmeans]{skmeans}).
#' @return A list containing the following elements:
#' \item{x}{The data matrix.}
#' \item{clusters}{the clustering labels.}
#' \item{k}{the number of clusters.}
#' \item{numobs}{the sample size.}
#' \item{p}{the vocabulary size.}
#' \item{likelihood}{vector containing the likelihood values at each iteration.}
#' \item{pi_hat}{estimated probabilities of belonging to the \code{k} clusters.}
#' \item{omega}{matrix containing the estimates of the omega parameters for each cluster.}
#' \item{f_z_x}{matrix containing the posterior probabilities of belonging to each cluster.}
#' \item{AIC}{Akaike Information Criterion (AIC) value of the fitted model.}
#' \item{BIC}{Bayesian Information Criterion (BIC) value of the fitted model.}
#'
#' @importFrom skmeans skmeans
#'
#' @details Starting from the data given by \code{x} the Mixture of Unigrams is fitted
#' and \code{k} clusters are obtained.
#' The algorithm for the parameter estimation is the Expectation-Maximization (EM).
#' In particular, the function assigns initial values to weights of the Multinomial distribution for each cluster
#' and inital weights for the components of the mixture. The estimates are obtained with maximum \code{n_it} steps of the
#' EM algorithm or until the tolerance level \code{eps} is reached; by using the posterior distribution
#' of the latent variable z, the documents are allocated to the cluster which maximizes the
#' posterior distribution.
#' For further details see the references.
#'
#' @references
#' Nigam, K., McCallum, A., Thrun, S., Mitchell, T.: Text classification from labeled and unlabeled documents using EM. Machine learning 39, 103-134 (2000).
#'
#' @examples
#' # Load the CNAE2 dataset
#' data("CNAE2")
#'
#' # Perform parameter estimation and clustering
#' mou_CNAE2 = mou_EM(x = CNAE2, k = 2)
#'
#' # Shows cluster labels to documents
#' mou_CNAE2$clusters
#'
#' @export
mou_EM = function (x, k, n_it = 100, eps = 1e-07, seed_choice = 1, init = NULL)
{
  set.seed(seed_choice)
  if (is.null(colnames(x)))
    colnames(x) = paste("W", sep = "", 1:dim(x)[2])
  numobs = nrow(x)
  p = ncol(x)
  ratio = 1000
  lik = -10^10
  if (is.null(init)) {
    out = skmeans(x, k, method = "pclust", control = list(maxiter = 20, nruns = 2, maxchains = 2))
    omega = abs(out$prototypes)
    omega = ifelse(omega == 0, 9.99988867182683e-321, omega)
  } else {
    for (i in 1:k) omega[, i] = colSums(x[init == i, ])
    omega = ifelse(omega == 0, 9.99988867182683e-321, omega)
  }
  pi_hat = rep(1/k, k)
  likelihood = NULL
  f_x_z = matrix(0, numobs, k)
  f_z_x = matrix(0, numobs, k)
  h = 0
  while ((h < n_it) & (ratio > eps)) {
    h = h + 1
    for (i in 1:k) f_x_z[, i] = rowSums(x * matrix(log(omega[ i, ]), numobs, p, byrow = TRUE))
    for (i in 1:k) f_x_z[, i] = exp(f_x_z[, i])
    for (i in 1:k) f_z_x[, i] = f_x_z[, i] * pi_hat[i]/(pi_hat %*% t(f_x_z))
    f_z_x = ifelse(is.na(f_z_x), mean(f_z_x, na.rm = T), f_z_x)
    for (i in 1:k) omega[i, ] = colSums(x * matrix(f_z_x[, i], numobs, p))
    omega = omega/matrix(rowSums(omega), k, p)
    omega = ifelse(omega == 0, 9.99988867182683e-321, omega)
    for (i in 1:k) pi_hat[i] = sum(f_z_x[, i])/numobs
    temp = pi_hat %*% t(f_x_z)
    temp = ifelse(temp == 0, 9.99988867182683e-321, temp)
    temp = sum(log(temp))
    likelihood = c(likelihood, temp)
    ratio = (temp - lik)/abs(lik)
    if (h < 3)  ratio = 2 * eps
    lik = temp
    if (is.na(lik))  ratio = eps/2
  }

  likelihood=likelihood[-1]
  lik_final = tail(likelihood, 1)
  param_aic =  p * k + k - 1
  AIC =  -2 * lik_final + param_aic * 2
  BIC =  -2 * lik_final + param_aic * log(numobs)

  clusters = apply(f_z_x, 1, which.max)
  output = list(x = x, k = k, numobs = numobs, p = p, likelihood = likelihood,
                clusters = clusters, pi_hat = pi_hat, omega = t(omega),
                f_z_x = f_z_x, AIC = AIC, BIC = BIC)
  class(output) = "deepMOU"
  invisible(output)
}
