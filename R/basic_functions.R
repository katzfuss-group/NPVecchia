#' Transforms hyperparameters to priors
#'
#' This function creates the priors on the variance and coefficients of the regressions from the three
#' hyperparameters described in the mathematical description. The coefficients are assumed to have zero
#' mean apriori, so only three prior elements are returned (the variance shape/scale, and the coefficnet variances).
#'
#' @param thetas 3 real numbers representing the three hyperparameters; as it is on the log scale,
#' it should generally be between -6 and 4 to avoid numerical issues (overflow or underflow)
#' @param n the number of locations
#' @param thresh the threshold for determining the number of neighbors based on the third
#' hyperparameter, defaults to 1e-3 
#'
#' @return List of priors, where 
#' 
#' the first element is a vector of length n containing the shape parameters of the IG prior on the variances, 
#' 
#' the second element similarly is of length n containing the corresponding scale parameters, and 
#' 
#' the last element is an n * m matrix, where each row contains the prior
#' variances for the regression coefficients (i.e. the diagonal of the prior covariance matrix)
#' @export
#'
#' @examples
#' 
#' thetas_ex <- c(1, 1, 1)
#' priors <- thetas_to_priors(thetas_ex, 100)
#' 
#' #with smaller threshold (leading to larger number of neighbors)
#' priors2 <- thetas_to_priors(thetas_ex, 500, thresh = 1e-6)
#' 
thetas_to_priors <- function(thetas, n, thresh = 1e-3) {
  if(length(thetas) != 3){
    stop("The number of hyperparameters (elements in the first argument) must be exactly 3!")
  }
  if(any(thetas > 4) || any(thetas < -6)){
    warning("A theta being too large/small will probably cause numerical issues.")
  }
  # Inverse-gamma scale prior parameter vector (prior on variances)
  b <- 5 * exp(thetas[[1]]) * (1 - exp(-exp(thetas[[2]])/sqrt(0:(n - 1))))
  # Inverse-gamma shape prior parameter vector (prior on variances)
  a <- rep(6, n)
  # temporary vector for determining number of neighbors
  # (based on the threshold provided as an input)
  tempor <- exp(-exp(thetas[[3]]) * (1:500))
  m <- which(tempor < thresh)[1] - 1
  # Force at least 2 neighbors, cannot be more than 500
  if (is.na(m) | m < 2) {
    m <- 2
  }
  # Create the prior on the coefficient variances
  g <- matrix(tempor[1:m], ncol = m, nrow = n, byrow = T)
  # Divide by mean of the IG prior for simplicity of derivations
  g <- g/(b/(a - 1))
  return(list(a, b, g))
}


#' Creates the simplist frequentist version of our method
#' 
#' For each regression, the MLE is used (i.e. the coefficients and residual standard error obtained from lm). 
#' This will only work when N > m (or lm will have issues), so call the function accordingly.
#'
#' @param datum an N * n matrix, where each row corresponds to N replications for that location
#' @param NNarray an n * m integer matrix giving the m nearest neighbors previous in the ordering (or 
#' outputting NAs if not available [i.e. there are not m previous points]) that are ordered 
#' from closest to furthest away. Required: m < N, m >= 2
#'
#' @return Sparse triangular matrix that is the Cholesky of the precision matrix \eqn{\Omega} 
#' such that \deqn{\Omega = U U'}
#' @export
#'
#' @examples
#' 
#' #create fake data and fake neighbor matrix
#' datum <- matrix(rnorm(1e4), nrow = 10)
#' NNarray <- matrix(NA, nrow = 1e3, ncol = 100)
#' #can only use previous points in ordering (this is actually 
#' #impossible in low dimensional space like this is designed for)
#' for(i in 1:100){
#'   NNarray[i:1e3, i] <- i
#' }
#' #Return sparse Cholesky of the precision matrix
#' uhat <- get_mle(datum, NNarray)
get_mle <- function(datum, NNarray) {
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
  # get n
  n <- ncol(datum)
  # get first regression variance
  d <- 1/sd(datum[, 1])
  # create sparse matrix with first entry d
  uhat <- sparseMatrix(i = 1, j = 1, x = d, dims = c(n, n), triangular = TRUE)
  for (i in 2:n) {
    # get indices of nearest neighbors (to regress column i on)
    gind <- na.omit(NNarray[i, ])
    # regress i on neighbors
    temp <- lm(datum[, i] ~ -1 + datum[, gind])
    # set the diagonal element to the regression SE
    d <- 1/(summary.lm(temp)$sigma)
    uhat[i, i] <- d
    # set the NN elements to the scaled coefficients
    uhat[gind, i] <- -coef(temp) * d
  }
  return(uhat)
}

#' Gets posteriors of Bayesian methodology
#' 
#' Given the priors calculated using thetas_to_priors, this function calculates the posterior 
#' distributions of the regression errors and coefficients.
#'
#' @param datum an N * n matrix of the data (N replications of n locations/variables)
#' @param priors a list of length 3 containing the priors for the shape of the IG prior,
#' the scale of the IG prior, and the prior variances of the coefficients (i.e. the output from
#' thetas_to_priors)
#' @param NNarray an n * m2 integer matrix giving the m nearest neighbors previous in the ordering (or 
#' outputting NAs if not available [i.e. there are not m previous points]) that are ordered 
#' from closest to furthest away. It is OK to have m2 > m, as it will be reduced to match the size
#' of the matrix g, but never have m2 < 2.
#'
#' @return List of posterior arguments, where
#' 
#' the first element is a vector of length n containing the posterior of the IG shape parameters,
#' 
#' the second element is a vector of length n containing the posteriors of the IG scale parameters,
#' 
#' the third element is an n * m matrix, where each row contains the posterior mean of the regression
#' coefficients (and if there are less than m, NAs fill the rest of the row),
#' 
#' the last element is an m * m \* n array, where each of the n slices contains the posterior variances
#' of the regression coefficients (and if there are less than m, NAs fill the rest of the slice)
#' 
#' @export
#'
#' @examples
#' 
#' #create fake data and fake neighbor matrix
#' datum <- matrix(rnorm(1e4), nrow = 10)
#' NNarray <- matrix(NA, nrow = 1e3, ncol = 100)
#' #can only use previous points in ordering (this is actually 
#' #impossible in low dimensional space like this is designed for)
#' for(i in 1:100){
#'   NNarray[i:1e3, i] <- i
#' }
#' priors <- thetas_to_priors(c(1, 1, 1), 1e3)
#' 
#' posteriors <- get_posts(datum, priors, NNarray)
#' 
get_posts <- function(datum, priors, NNarray) {
  # check if priors is of correct length
  if(length(priors) != 3) {
    stop("Priors should be a list of length 3!")
  }
  #get b, g priors from the list
  b <- priors[[2]]
  g <- priors[[3]]
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
    muhat <- tryCatch({
      solve(Ginv_chol, solve(t(Ginv_chol), crossprod(xi, yi)))
    }, error = function(e) {
      ginv(Ginv) %*% crossprod(xi, yi)
    })
    # G_post[1:nn,1:nn,i] <- G
    # Fill in full posteriors with the elements from i'th regression
    muhat_post[i, 1:nn] <- muhat
    G_post[1:nn, 1:nn, i] <- Ginv_chol
    muhat_post[i, 1:nn] <- muhat
    b_post[i] <- b[i] + (t(yi) %*% yi - t(muhat) %*% Ginv %*% muhat)/2
  }
  #Return list of posteriors
  return(list(a_post, b_post, muhat_post, G_post))
}

#' Creates posterior mean sparse matrix from posteriors
#' 
#' This function uses the posterior arguments to create a posterior mean estimate of the Cholesky of
#' the precision matrix. This function creates the Bayesian version of get_mle or its variants. See
#' the mathematical description for details of the methodology.
#'
#' @param posts a List of the posteriors from get_posts (or get_posts_c); alternatively it can be
#' custom values as long as the sizes match the output from get_posts.
#' @param NNarray an n * m2 integer matrix giving the m nearest neighbors previous in the ordering (or 
#' outputting NAs if not available [i.e. there are not m previous points]) that are ordered 
#' from closest to furthest away. It is OK to have m2 large, as it will be reduced to match the size
#' of the posterior means (i.e. number of columns in the third element of the posteriors), but
#' never have m2 < 2.
#'
#' @return Sparse triangular matrix that is the Cholesky of the precision matrix \eqn{\Omega} 
#' such that \deqn{\Omega = U U'}
#' @export
#'
#' @examples
#' 
#' #create fake data and fake neighbor matrix
#' datum <- matrix(rnorm(1e4), nrow = 10)
#' NNarray <- matrix(NA, nrow = 1e3, ncol = 100)
#' #can only use previous points in ordering (this is actually 
#' #impossible in low dimensional space like this is designed for)
#' for(i in 1:100){
#'   NNarray[i:1e3, i] <- i
#' }
#' priors <- thetas_to_priors(c(1, 1, 1), 1e3)
#' 
#' posteriors <- get_posts(datum, priors, NNarray)
#' 
#' uhat <- samp_posts(posteriors, NNarray)
#' 
samp_posts <- function(posts, NNarray) {
  # get n, m
  n <- nrow(NNarray)
  m <- ncol(posts[[3]])
  
  # make sure all elements in posts have a place in uhat
  if(ncol(NNarray) < m) {
    stop("The posteriors have more neighbors than the neighbor matrix accounts for!")
  }
  # check if NNarray is an integer matrix; warn if not
  if(! is.integer(NNarray)){
    warning("NNarray should consist of only integers (and/or NAs)!")
  }
  # make sure posts has 4 elements
  if(length(posts) != 4) {
    stop("There should be 4 elements in the list of posteriors! For simplicity,
         it is recommended to use get_posts to generate these posteriors.")
  }
  # make sure posts elements are of the correct sizes
  if(length(posts[[1]]) != n || length(posts[[2]]) != n || !all(dim(posts[[3]]) == c(n, m)) ||
  !all(dim(posts[[4]]) == c(m, m, n))){
    stop("Please use get_posts to generate the posteriors. Some of 
          the current posteriors have incorrect dimensions.")
  }
  
  # Get posterior mean of d and create the sparse matrix
  d <- (1/sqrt(posts[[2]][1])) * exp(lgamma((2 * posts[[1]][1] + 1)/2) - lgamma(posts[[1]][1]))
  uhat <- sparseMatrix(i = 1, j = 1, x = d, dims = c(n, n), triangular = TRUE)
  for (i in 2:n) {
    # get nearest neighbors and number of them (nn)
    gind <- na.omit(NNarray[i, 1:m])
    nn <- length(gind)
    # Fill in appropriate elements of sparse matrix with posterior means
    uhat[i, i] <- (1/sqrt(posts[[2]][i])) * exp(lgamma((2 * posts[[1]][i] + 1)/2) - lgamma(posts[[1]][i]))
    uhat[gind, i] <- posts[[3]][i, 1:nn] * uhat[i, i]
  }
  return(uhat)
}

#' Calculates the integrated log likelihood
#' 
#' Given the three hyperparameters, this calculates the negative of the integrated log-likelihood.
#' See the mathematical description for further explanation as needed. It is advised 
#' to use the C++ version minus_loglikeli_c to improve speed considerably. This 
#' function is often used with an optimizer or MCMC method to find the optimal 
#' hyperparameters.
#'
#' @param thetas 3 real numbers representing the three hyperparameters; as it is on the log scale,
#' it should generally be between -6 and 4 to avoid numerical issues (overflow or underflow)
#' @param datum an N * n matrix of the data (N replications of n locations/variables)
#' @param NNarray an n * m2 integer matrix giving the m nearest neighbors previous in the ordering (or 
#' outputting NAs if not available [i.e. there are not m previous points]) that are ordered 
#' from closest to furthest away. It is OK to have m2 large, as it will be reduced to match the size
#' of the posterior means (i.e. number of columns in the third element of the posteriors), but
#' never have m2 < 2.
#' @param threshh threshold for number of neighbors (for thetas_to_priors); defaults
#' to 1e-3
#'
#' @return a numeric value (the negative log likelihood)
#' @export
#'
#' @examples
#' 
#' #create fake data and fake neighbor matrix
#' datum <- matrix(rnorm(1e4), nrow = 10)
#' NNarray <- matrix(NA, nrow = 1e3, ncol = 100)
#' #can only use previous points in ordering (this is actually 
#' #impossible in low dimensional space like this is designed for)
#' for(i in 1:100){
#'   NNarray[i:1e3, i] <- i
#' }
#' 
#' #calculates log likelihood
#' minus_loglikeli(c(1, 1, 1), datum, NNarray)
#' 
#' 
minus_loglikeli <- function(thetas, datum, NNarray, threshh = 1e-3) {
  # get n, N
  n <- nrow(NNarray)
  N <- nrow(datum)
  # get alpha posterior
  a_post <- 6 + N/2
  # get priors
  pr <- thetas_to_priors(thetas, n, thresh = threshh)
  # get m
  # make sure m is not greater than max number of neighbors set
  # (from NNarray creation)
  m <- min(ncol(pr[[3]]), ncol(NNarray))
  # get needed priors from the list of priors
  b <- pr[[2]]
  g <- pr[[3]]
  # get the first element for the log-likelihood
  loglikelihood <- 6 * log(b[1]) - a_post * log(b[1] + crossprod(datum[, 1])/2)
  for (i in 2:n) {
    # get nearest neighbors and how many there are as nn
    gind <- na.omit(NNarray[i, 1:m])
    nn <- length(gind)
    # set up i'th regression with basic notation yi ~ xi
    xi <- -datum[, gind]
    yi <- datum[, i]
    # Get inverse of G posterior
    Ginv <- crossprod(xi) + diag(g[i,1:nn]^(-1), nrow = nn)
    # Take the Cholesky of Ginv
    Ginv_chol <- chol(Ginv)
    # Try to get muhat directly, but if Ginv_chol is not invertible, use 
    # the generalized inverse
    muhat <- tryCatch({
      solve(Ginv_chol, solve(t(Ginv_chol), crossprod(xi, yi)))
    }, error = function(e) {
      ginv(Ginv) %*% crossprod(xi, yi)
    })
    # get the posterior of b (IG scale parameter)
    b_post <- b[i] + (crossprod(yi) - t(muhat) %*% Ginv %*% muhat)/2
    # calculate the determinant term of the integrated likelihood
    log_det <- -0.5 * (2 * sum(log(na.omit(diag(Ginv_chol)))) + (sum(log(g[i, 1:nn]))))
    # calculate the term based on the ratio of IG parameters
    log_ig <- 6 * log(b[i]) - a_post * log(b_post)
    # Add these values to the log integrated likelihood
    loglikelihood <- loglikelihood + log_det + log_ig
  }
  # Return the negative of the log integrated likelihood
  return(c(-loglikelihood))
}
