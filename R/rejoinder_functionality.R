#' Shrink towards a parametric covariance matrix
#' The one difference is that four parameters are printed instead of the usual
#' three used previously (non-zero mean for U). Therefore, downstream functions 
#' must be slightly updated.
#'
#' @param covariance parametric covariance matrix
#' @param m number of neighbors, must be specified
#' @param NNarray nearest neighbor matrix as usual
#' @param alpha alpha in the IG prior (related to the rejoinder c_d = sqrt(alpha - 2))
#' @param a_beta c_u in the rejoinder, the variance multiplier for U
#'
#' @return a list of the 4 prior elements: vectors a and b for the IG(a,b) prior on 
#' D elements, then matrices mean and var for means of U and the elements of the
#' diagonal prior covariance of U.
#' @export
#'
#' @examples
prior_covshrink <- function(covariance, m, NNarray, alpha = 3, a_beta = 0.5) {
  npt <- nrow(covariance)
  #initialize priors
  beta_mean <- matrix(NA, nr = npt, nc = m)
  a <- rep(alpha, npt)
  b <- rep(0, npt)
  b[1] <- (alpha - 1) * covariance[1,1]
  #loop over points and calculate priors using equations in Section 2.9 (basic Vecchia)
  for(i in 2:npt){
    gind <- na.omit(NNarray[i, 1:m])
    nn <- length(gind)
    temp <- -solve(covariance[gind, gind]) %*% covariance[gind, i]
    beta_mean[i, 1:nn] <- temp
    d_mean <- covariance[i,i] + t(temp) %*% covariance[gind, i]
    b[i] <- d_mean * (a[i] - 1)
  }
  #calculate beta variances
  beta_variances <- a_beta * beta_mean^2 / (b / (a - 1))
  return(list(a, b, beta_mean, beta_variances))
}

#' Updated posterior calculation with non-zero mean of beta
#'
#' @param datum data to use for getting posteriors
#' @param priors list of 4 (output from prior_covshrink)
#' @param NNarray nearest neighbor matrix as usual
#'
#' @return a list of the posteriors (same as priors but updated, except for
#' the covariance of the U elements isn't diagonal, making it a cube)
#' @export
#'
#' @examples
get_posts_shrink <- function(datum, priors, NNarray) {
  # check if priors is of correct length
  if(length(priors) != 4) {
    stop("Priors should be a list of length 4!")
  }
  #get b, g priors from the list
  b <- priors[[2]]
  bmean <- priors[[3]]
  g <- priors[[4]]
  # get n, N, m
  n <- ncol(datum)
  N <- nrow(datum)
  m <- min(ncol(g), ncol(NNarray))
  
  # make sure datum and NNarray are both matrices
  if(!(is.matrix(datum) && is.matrix(NNarray))) {
    stop("The data and NNarray must both be matrices")
  }
  # make sure NNarray is of correct size
  if(ncol(datum) != nrow(NNarray)){
    stop(paste("The number of locations (", ncol(datum), ") must equal the number of 
               rows of the neighbor matrix but given (", nrow(NNarray), ")", sep=""))
  }
  if(ncol(NNarray) < 2) {
    stop("At least 2 neighbors are required (2 or more columns in NNarray)")
  }
  # check if NNarray is an integer matrix
  if(! is.integer(NNarray)){
    warning("NNarray should consist of only integers (and/or NAs)!")
  }
  # check if a, b, g are correct sizes
  if(length(priors[[1]]) != n || length(b) != n || nrow(g) != n || ncol(g) < 2){
    stop("Please use priors created by thetas_to_priors for convenience. 
          The current priors are of the wrong size.")
  }
  
  # Create vectors/matrices/arrays to hold the posteriors
  a_post <- rep(0, n)
  b_post <- rep(0, n)
  muhat_post <- matrix(NA, nrow = n, ncol = m)
  G_post <- array(NA, dim = c(m, m, n))
  # Get posterior of a
  a_post <- priors[[1]] + N/2
  # Get first element of posterior of b
  b_post[1] <- b[1] + t(datum[, 1] %*% datum[, 1])/2
  for (i in 2:n) {
    # set up the regression as column i ~ nearest neighbors
    gind <- na.omit(NNarray[i, 1:m])
    nn <- length(gind)
    xi <- -datum[, gind]
    yi <- datum[, i]
    # Get inverse of posterior of G (coefficient variance)
    Ginv <- t(xi) %*% xi + diag(g[i, 1:nn]^(-1), nrow = nn)
    # take Cholesky
    Ginv_chol <- chol(Ginv)
    # G <- ginv(Ginv); muhat <- G%*%t(xi)%*%yi
    # Try to solve for muhat directly which occasionally has
    # numerical issues. If so, use generalized inverse of Ginv
    
    #the added non-zero prior mean of U
    prior_infl <- diag(g[i, 1:nn]^(-1), nrow = nn) %*% na.omit(bmean[i,])
    
    muhat <- tryCatch({
      solve(Ginv_chol, solve(t(Ginv_chol), prior_infl + crossprod(xi, yi)))
    }, error = function(e) {
      ginv(Ginv) %*% (prior_infl + crossprod(xi, yi))
    })
    # G_post[1:nn,1:nn,i] <- G
    # Fill in full posteriors with the elements from i'th regression
    muhat_post[i, 1:nn] <- muhat
    G_post[1:nn, 1:nn, i] <- Ginv_chol
    muhat_post[i, 1:nn] <- muhat
    b_post[i] <- b[i] + (t(yi) %*% yi - t(muhat) %*% Ginv %*% muhat + 
                           sum(na.omit(bmean[i,]) * prior_infl))/2
  }
  #Return list of posteriors
  return(list(a_post, b_post, muhat_post, G_post))
}

#' Prediction at new point
#'
#' @param covariance parametric covariance matrix (last element is the point to predict at)
#' @param m number of neighbors
#' @param NNarray nearest neighbor matrix as usual or vector of NN for the point to predict
#' @param dat data
#' @param alpha alpha in the IG prior (related to the rejoinder c_d = sqrt(alpha - 2))
#' @param a_beta c_u in the rejoinder, the variance multiplier for U
#'
#' @return a names list with the 3 parameters of the multivariate t distribution
#' @export
#'
#' @examples
covshrink_predictive <- function(covariance, m, NNarray, dat, alpha = 3, a_beta = 0.5) {
  #number of points
  npt <- nrow(covariance)
  #if its a matrix, use the last points NN or use the NN provided
  if(is.matrix(NNarray)){
    gind <- na.omit(NNarray[npt, 1:m])
  }else{
    gind <- na.omit(NNarray[1:m])
  }
  #number of NN
  nn <- length(gind)
  #temporary partial calculation of S^{-1} * cov, also exactly the mean of the new point
  temp <- -solve(covariance[gind, gind]) %*% covariance[gind, npt]
  bmean <- temp
  #update conditional variance
  d_mean <- covariance[npt,npt] + t(temp) %*% covariance[gind, npt]
  b <- c(d_mean * (alpha - 1))
  beta_variances <- a_beta * bmean^2 / (b / (alpha - 1))
  #as it is conjugate, only need to calculate the t distribution parameters
  return(list(prior_center = sum(-bmean * dat[gind]),df = 2*alpha, prior_prec = (b/alpha)*(1 + sum(dat[gind]^2*beta_variances))))
  # post_center = sum(-muhat * dat[gind]), post_df = (2*apost), post_prec = (bpost / apost) * (1 + c(t(dat[gind]) %*%solve(Ginv, dat[gind])))))
}


