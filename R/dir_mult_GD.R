#' @title Dirichlet-Multinomial mixture model by Gradient Descend algorithm
#'
#' @description Performs parameter estimation by means of a Gradient Descend algorithm and cluster allocation
#' for the Dirichlet-Multinomial mixture model.
#'

#' @param x            Document-term matrix describing the frequency of terms that occur in a collection of documents. Rows correspond to documents in the collection and columns correspond to terms.
#' @param k            Number of clusters/groups.
#' @param n_it         Number of Gradient Descend steps.
#' @param eps          Tolerance level for the convergence of the algorithm. Default is \code{1e-05}.
#' @param seed_choice  Set seed for reproducible results.
#' @param KK           Maximum number of iterations allowed for the \link[stats]{nlminb} function (see below).
#' @param min_iter     Minimum number of Gradient Descend steps.
#' @param init         Vector containing the initial document allocations for the initialization of the algorithm. If \proglang{NULL} (default) initialization is carried out via spherical k-means (\link[skmeans]{skmeans}).
#' @return A list containing the following elements:
#' \item{x}{The data matrix.}
#' \item{clusters}{the clustering labels.}
#' \item{k}{the number of clusters.}
#' \item{numobs}{the sample size.}
#' \item{p}{the vocabulary size.}
#' \item{likelihood}{vector containing the likelihood values at each iteration.}
#' \item{pi_hat}{estimated probabilities of belonging to the \code{k} clusters.}
#' \item{Theta}{matrix containing the estimates of the Theta parameters for each cluster.}
#' \item{f_z_x}{matrix containing the posterior probabilities of belonging to each cluster.}
#' \item{AIC}{Akaike Information Criterion (AIC) value of the fitted model.}
#' \item{BIC}{Bayesian Information Criterion (BIC) value of the fitted model.}
#'
#'
#' @details Starting from the data given by \code{x} the Dirichlet-Multinomial mixture model is fitted
#' and \code{k} clusters are obtained.
#' The algorithm for the parameter estimation is the Gradiend Descend.
#' In particular, the function assigns initial values to weights of the Dirichlet-Multinomial distribution for each cluster
#' and inital weights for the elements of the mixture. The estimates are obtained with maximum \code{n_it} steps of the
#' Descent Algorithm algorithm or until a tolerance level \code{eps} is reached; by using the posterior distribution
#' of the latent variable z, the documents are allocated to the cluster which maximizes the
#' posterior distribution.
#' For further details see the references.
#'
#' @references
#' Anderlucci Laura, Viroli Cinzia (2020). "Mixtures of Dirichlet-Multinomial distributions for supervised and unsupervised classification of short text data". \emph{Advances in Data Analysis and Classification}. \doi{10.1007/s11634-020-00399-3}
#'
#' @examples
#' # Load the CNAE2 dataset
#' data("CNAE2")
#'
#' # Perform parameter estimation and clustering, very
#' # few iterations are used for this example
#' dir_CNAE2 = dir_mult_GD(x = CNAE2, k = 2, n_it = 2)
#'
#' # Shows cluster labels to documents
#' dir_CNAE2$clusters
#'
#' @export

dir_mult_GD = function (x, k, n_it = 100, eps = 1e-05, seed_choice = 1, KK = 20, min_iter = 2, init = NULL)  {

    set.seed(seed_choice)
    progress_bar = txtProgressBar(min = 0, max = n_it + 1, style = 3)
    setTxtProgressBar(progress_bar, 0)
    if (is.null(colnames(x)))
      colnames(x) = paste("W", sep = "", 1:dim(x)[2])
    numobs = nrow(x)
    p = ncol(x)
    if (is.null(init)) {
      Theta = matrix(0, p, k)
      rownames(Theta) = colnames(x)
      pi_hat = rep(1/k, k)
      out = skmeans(x, k, method = "pclust", control = list(maxiter = 20,
                                                            nruns = 2, maxchains = 2))
      for (i in 1:k) Theta[, i] = colSums(x[out$cl == i, ])
      Theta = ifelse(Theta == 0, exp(-500), Theta)
    } else {
      Theta = matrix(0, p, k)
      rownames(Theta) = colnames(x)
      pi_hat = table(init)
      pi_hat = pi_hat/sum(pi_hat)
      for (i in 1:k) Theta[, i] = colSums(x[init == i, ])
      Theta = ifelse(Theta == 0, exp(-500), Theta)
    }

    f_z_x = f_x_z = matrix(0, numobs, k)
    likelihood = NULL
    ratio = 2 * eps
    lik = -10^10
    h = 0
    while ((h < n_it) & (ratio > eps)) {
      h = h + 1
      setTxtProgressBar(progress_bar, h)
      for (i in 1:k) f_x_z[, i] = calc_log_f_x_z(x, Theta[,
                                                          i])
      f_x_z = exp(f_x_z)
      f_x_z = ifelse(f_x_z == 0, exp(-740), f_x_z)
      f_y = f_x_z %*% pi_hat
      f_z_x = f_x_z * matrix(pi_hat, numobs, k, byrow = T)/matrix(f_y,
                                                                    numobs, k)
      pi_hat = colMeans(f_z_x)
      for (i in 1:k) {
        out_Theta = try(stats::nlminb(Theta[, i], Theta_f, gradient = Theta_f_grad,
                               lower = 0.001, upper = 500, x = x, f_z_x = f_z_x[,
                                                                                i], control = list(iter.max = KK)))
        if (!is.character(out_Theta))
          Theta[, i] = out_Theta$par
      }
      temp = sum(log(f_y))
      likelihood = c(likelihood, temp)
      ratio = (temp - lik)/abs(lik)
      if (h < min_iter)
        ratio = 2 * eps
      lik = temp
      if (is.na(lik))
        ratio = eps/2
    }
    setTxtProgressBar(progress_bar, n_it + 1)
    clusters = apply(f_z_x, 1, which.max)

    lik_final = tail(likelihood, 1)
    param_aic =  p * k + k - 1
    AIC =  -2 * lik_final + param_aic * 2
    BIC =  -2 * lik_final + param_aic * log(numobs)

    output = list(x = x, k = k, numobs = numobs, p = p, likelihood = likelihood,
                  clusters = clusters, pi_hat = pi_hat, Theta = Theta,
                  f_z_x = f_z_x, AIC = AIC, BIC = BIC)
    class(output) = "deepMOU"
    return(output)
}



# Auxiliary functions for main algorithm

# f(x|z) calc
calc_log_f_x_z = function( x, Theta) {
  numobs = nrow(x)
  p = ncol(x)
  Theta= ifelse( Theta==0, 0.001, Theta)
  Theta_val = rep( sum(Theta), numobs)
  log_f_x_z =  Rfast::Lgamma(Theta_val) -
    Rfast::Lgamma( rowSums(x) + Theta_val) +
    rowSums( Rfast::Lgamma( x + matrix( Theta, numobs, p, byrow=T) ) ) -
    rowSums( Rfast::Lgamma( matrix( Theta, numobs, p, byrow = T) ) )
  return(log_f_x_z)
}

Theta_f = function( Theta, x, f_z_x){
  numobs = nrow(x)
  p = ncol(x)
  log_f_x_z = calc_log_f_x_z( x, Theta)
  return( -sum( log_f_x_z * f_z_x))
}


Theta_f_grad = function( Theta, x, f_z_x) {
  numobs = nrow(x)
  p = ncol(x)
  temp1 = array( NA, c( numobs, p))
  cost = matrix( f_z_x, numobs, p)
  temp1 = calc2grad_b( x, Theta)
  return( -colSums(temp1 * cost) )
}


calc2grad_b = function( x, Theta){
  numobs = nrow(x)
  p = ncol(x)
  alpha_star_v = Theta
  alpha_star = rep( sum(Theta), p)

  A = t( digamma( alpha_star ) ) ### 1 x p
  B = digamma( matrix( rowSums(x), numobs, p) + t( matrix( alpha_star, p, numobs) ) ) # n x p
  C = t( digamma( alpha_star_v) ) ### 1 x p
  D = digamma( x + matrix( alpha_star_v, numobs, p, byrow=TRUE)) ### n x p

  out = matrix( A, numobs, p, byrow = TRUE) - B - matrix( C, numobs, p, byrow=TRUE) + D
  return(out)
}
