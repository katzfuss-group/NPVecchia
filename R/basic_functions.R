#' Transforms hyperparameters to priors
#'
#' This function creates the priors on the variance and coefficients of the regressions from the three
#' hyperparameters described in the mathematical description. The coefficients are assumed to have zero
#' mean apriori, so only three prior elements are returned (the variance shape/scale, and the coefficnet variances).
#'
#' @param thetas 3 real numbers representing the log of the three hyperparameters
#' @param n the number of locations
#' @param thresh the threshold for determining the number of neighbors based on the third
#' hyperparameter, defaults to 1e-3 
#'
#' @return List of priors, where 
#' 
#' the first element is a vector of length n2 containing the shape parameters of the IG prior on the variances, 
#' 
#' the second element similarly is of length n2 containing the corresponding scale parameters, and 
#' 
#' the last element is an n2 * m matrix, where each row contains the prior
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
  b <- 5 * exp(thetas[[1]]) * (1 - exp(-exp(thetas[[2]])/sqrt(0:(n - 1))))
  a <- rep(6, n)
  tempor <- exp(-exp(thetas[[3]]) * (1:500))
  m <- which(tempor < thresh)[1] - 1
  if (is.na(m) | m < 2) {
    m <- 2
  }
  g <- matrix(tempor[1:m], ncol = m, nrow = n, byrow = T)
  g <- g/(b/(a - 1))
  return(list(a, b, g))
}

create_data <- function(covar_true, N) {
  
  # Lazy creation of data
  datum <- mvrnorm(N, rep(0, nrow(covar_true)), covar_true)
  return(datum)
}

#' Creates the simplist frequentist version of our method
#' 
#' For each regression, the MLE is used (i.e. the coefficients and residual standard error obtained from lm). 
#' This will only work when N > m (or lm will have issues), so call the function accordingly.
#'
#' @param dat an N * n matrix, where each row corresponds to N replications for that location
#' @param NNarray an n * m matrix giving the m nearest neighbors previous in the ordering (or 
#' outputting NAs if not available [i.e. there are not m previous points]) that are ordered 
#' for nearest to furthest. Required: m < N
#'
#' @return Sparse triangular matrix that is the Cholesky of the precision matrix \eqn{\Omega} 
#' such that \deqn{\Omega = U U'}
#' @export
#'
#' @examples
get_mle <- function(dat, NNarray) {
  n2 <- ncol(dat)
  d <- 1/sqrt(sd(dat[, 1]))
  uhat <- sparseMatrix(i = 1, j = 1, x = d, dims = c(n2, n2), triangular = TRUE)
  for (i in 2:n2) {
    gind <- na.omit(NNarray[i, ])
    temp <- lm(dat[, i] ~ -1 + dat[, gind])
    d <- 1/(summary.lm(temp)$sigma)
    uhat[i, i] <- d
    uhat[gind, i] <- -coef(temp) * d
  }
  return(uhat)
}

#' Gets posteriors of Bayesian methodology
#' 
#' Given the priors calculated using thetas_to_priors, this function calculates the posterior 
#' distributions of the regression errors and coefficients.
#'
#' @param x an N * n matrix of the data (N replications of n locations/variables)
#' @param a vector of length n of IG scale priors (thetas_to_priors()[[1]])
#' @param b vector of length n of IG shape priors (thetas_to_priors()[[2]])
#' @param g matrix of dimension n * m of prior coefficient variances (thetas_to_priors()[[3]])
#' @param NNarray an n * m2 matrix giving the m nearest neighbors previous in the ordering (or 
#' outputting NAs if not available [i.e. there are not m previous points]) that are ordered 
#' for nearest to furthest. It is OK to have m2 > m, as it will be reduced to match the size
#' of the matrix g.
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
#' data <- matrix(rnorm(1e4), nrow = 10)
#' NNarray <- matrix(NA, nrow = 1e3, ncol = 100)
#' #can only use previous points in ordering (this is actually impossible in low dimensional space like this is designed for)
#' for(i in 1:100){
#'   NNarray[i:1e3,i] <- i
#' }
#' priors <- thetas_to_priors(c(1, 1, 1), 1e3)
#' 
#' posteriors <- get_posts(data, priors[[1]], priors[[2]], priors[[3]], NNarray)
#' 
get_posts <- function(x, a, b, g, NNarray) {
  n2 <- ncol(x)
  N <- nrow(x)
  m <- ncol(g)
  a_post <- rep(0, n2)
  b_post <- rep(0, n2)
  muhat_post <- matrix(NA, nr = n2, nc = m)
  G_post <- array(NA, dim = c(m, m, n2))
  a_post <- a + N/2
  b_post[1] <- b[1] + t(x[, 1] %*% x[, 1])/2
  for (i in 2:n2) {
    gind <- na.omit(NNarray[i, 1:m])
    nn <- length(gind)
    xi <- -x[, gind]
    yi <- x[, i]
    Ginv <- t(xi) %*% xi + diag(g[i, 1:nn]^(-1), nrow = nn)
    Ginv_chol <- chol(Ginv)
    # G <- ginv(Ginv) muhat <- G%*%t(xi)%*%yi
    muhat <- solve(Ginv_chol, solve(t(Ginv_chol), t(xi) %*% yi))
    # G_post[1:nn,1:nn,i] <- G
    muhat_post[i, 1:nn] <- muhat
    G_post[1:nn, 1:nn, i] <- Ginv_chol
    muhat_post[i, 1:nn] <- muhat
    b_post[i] <- b[i] + (t(yi) %*% yi - t(muhat) %*% Ginv %*% muhat)/2
  }
  return(list(a_post, b_post, muhat_post, G_post))
}

samp_posts <- function(posts, NNarray) {
  n2 <- nrow(NNarray)
  m <- ncol(NNarray)
  d <- (1/sqrt(posts[[2]][1])) * exp(lgamma((2 * posts[[1]][1] + 1)/2) - lgamma(posts[[1]][1]))
  uhat <- sparseMatrix(i = 1, j = 1, x = d, dims = c(n2, n2), triangular = TRUE)
  for (i in 2:n2) {
    gind <- na.omit(NNarray[i, ])
    nn <- length(gind)
    # d <- rinvgamma(mcl, posts[[1]][i], posts[[2]][i])
    uhat[i, i] <- (1/sqrt(posts[[2]][i])) * exp(lgamma((2 * posts[[1]][i] + 1)/2) - lgamma(posts[[1]][i]))
    # uhat[i,i] <- mean(1/sqrt(d))
    uhat[gind, i] <- posts[[3]][i, 1:nn] * uhat[i, i]
  }
  return(uhat)
}

minus_loglikeli <- function(x, datum, NNarray, N) {
  n2 <- nrow(NNarray)
  ap <- 6 + N/2
  pr <- thetas_to_priors(x, n2)
  m <- min(ncol(pr[[3]]), ncol(NNarray))
  if (m < 2) {
    m <- 2
  }
  b <- pr[[2]]
  g <- pr[[3]]
  sums <- 6 * log(b[1]) - ap * log(b[1] + crossprod(datum[, 1])/2)
  for (i in 2:n2) {
    gind <- na.omit(NNarray[i, 1:m])
    nn <- length(gind)
    # browser()
    xi <- -datum[, gind]
    yi <- datum[, i]
    Ginv <- tryCatch({
      crossprod(xi) + diag(g[i, 1:nn]^(-1), nrow = nn)
    }, error = function(e) {
      show(dim(crossprod(xi)))
      show(m)
      show(nn)
      show(length(na.omit(NNarray[i, 1:m])))
      show(i)
      break
    })
    # Ginv <- crossprod(xi) + diag(g[i,1:nn]^(-1),nrow=nn)
    Ginv_chol <- chol(Ginv)
    # G <- ginv(Ginv) muhat <- G%*%t(xi)%*%yi muhat <- solve(Ginv_chol,
    # solve(t(Ginv_chol),initt[[i]][[2]]))
    muhat <- tryCatch({
      solve(Ginv_chol, solve(t(Ginv_chol), crossprod(xi, yi)))
    }, error = function(e) {
      ginv(Ginv) %*% crossprod(xi, yi)
    })
    # G_post[1:nn,1:nn,i] <- G
    b_post <- b[i] + (crossprod(yi) - t(muhat) %*% Ginv %*% muhat)/2
    ldet <- -0.5 * (2 * sum(log(na.omit(diag(Ginv_chol)))) + (sum(log(g[i, 1:nn]))))
    lb <- 6 * log(b[i]) - ap * log(b_post)
    sums <- sums + ldet + lb
  }
  return(c(-sums))
}