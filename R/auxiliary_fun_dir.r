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



misc <-  function(classification, truth)
  {

    if (length(classification)>1) {

      q <- function(map, len, x) {
        x <- as.character(x)
        map <- lapply(map, as.character)
        y <- sapply(map, function(x) x[1])
        best <- y != x
        if (all(len) == 1)
          return(best)
        errmin <- sum(as.numeric(best))
        z <- sapply(map, function(x) x[length(x)])
        mask <- len != 1
        counter <- rep(0, length(len))
        k <- sum(as.numeric(mask))
        j <- 0
        while (y != z) {
          i <- k - j
          m <- mask[i]
          counter[m] <- (counter[m]%%len[m]) + 1
          y[x == names(map)[m]] <- map[[m]][counter[m]]
          temp <- y != x
          err <- sum(as.numeric(temp))
          if (err < errmin) {
            errmin <- err
            best <- temp
          }
          j <- (j + 1)%%k
        }
        best
      }
      if (any(isNA <- is.na(classification))) {
        classification <- as.character(classification)
        nachar <- paste(unique(classification[!isNA]), collapse = "")
        classification[isNA] <- nachar
      }
      MAP <- mapClass(classification, truth)
      len <- sapply(MAP[[1]], length)
      if (all(len) == 1) {
        CtoT <- unlist(MAP[[1]])
        I <- match(as.character(classification), names(CtoT),
                   nomatch = 0)
        one <- CtoT[I] != truth
      }
      else {
        one <- q(MAP[[1]], len, truth)
      }
      len <- sapply(MAP[[2]], length)
      if (all(len) == 1) {
        TtoC <- unlist(MAP[[2]])
        I <- match(as.character(truth), names(TtoC), nomatch = 0)
        two <- TtoC[I] != classification
      }
      else {
        two <- q(MAP[[2]], len, classification)
      }
      err <- if (sum(as.numeric(one)) > sum(as.numeric(two)))
        as.vector(one)
      else as.vector(two)
      bad <- seq(along = classification)[err]
      errorRate = length(bad)/length(truth)

    } else errorRate<-as.numeric(classification!=truth)

    return(errorRate)
  }


mapClass <-  function (a, b)
  {
    l <- length(a)
    x <- y <- rep(NA, l)
    if (l != length(b)) {
      warning("unequal lengths")
      return(x)
    }
    aChar <- as.character(a)
    bChar <- as.character(b)
    Tab <- table(a, b)
    Ua <- dimnames(Tab)[[1]]
    Ub <- dimnames(Tab)[[2]]
    aTOb <- rep(list(Ub), length(Ua))
    names(aTOb) <- Ua
    bTOa <- rep(list(Ua), length(Ub))
    names(bTOa) <- Ub
    k <- nrow(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 1, max)
    for (i in 1:k) {
      I <- match(Max[i], Tab[i, ], nomatch = 0)
      aTOb[[i]] <- Ub[I]
    }
    if (is.numeric(b))
      aTOb <- lapply(aTOb, as.numeric)
    k <- ncol(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 2, max)
    for (j in (1:k)) {
      J <- match(Max[j], Tab[, j])
      bTOa[[j]] <- Ua[J]
    }
    if (is.numeric(a))
      bTOa <- lapply(bTOa, as.numeric)
    list(aTOb = aTOb, bTOa = bTOa)
  }

#===========ADJ. RAND INDEX (Adelchi) ================
# compute "adjusted Rand index" to measure the degree
# of agreement between two classifications (02-11-2004)
#
adj.rand.index <- function (table){
  f2 <- function(n) n*(n-1)/2
  sum.f2 <- function(v) sum(f2(v))
  marg.1 <- apply(table,1,sum)
  marg.2 <- apply(table,2,sum)
  n<- sum(table)
  prod <- sum.f2(marg.1)*sum.f2(marg.2)/f2(n)
  num <- (sum.f2(as.vector(table))- prod)
  den <- 0.5*(sum.f2(marg.1)+sum.f2(marg.2))-prod
  num/den
}
