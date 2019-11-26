thetas_to_priors <- function(thetas, n2, thresh = 1e-3) {
  b <- 5 * exp(thetas[[1]]) * (1 - exp(-exp(thetas[[2]])/sqrt(0:(n2 - 1))))
  a <- rep(6, n2)
  tempor <- exp(-exp(thetas[[3]]) * (1:500))
  m <- which(tempor < thresh)[1] - 1
  if (is.na(m) | m < 2) {
    m <- 2
  }
  g <- matrix(tempor[1:m], nc = m, nr = n2, byrow = T)
  g <- g/(b/(a - 1))
  return(list(a, b, g))
}

create_data <- function(covar_true, N) {
  
  # Lazy creation of data
  datum <- mvrnorm(N, rep(0, nrow(covar_true)), covar_true)
  return(datum)
}

get_mle <- function(dat, NNarray) {
  n2 <- ncol(dat)
  d <- 1/sqrt(sd(dat[, 1]))
  uhat <- sparseMatrix(i = 1, j = 1, x = d, dims = c(n2, n2), triangular = TRUE)
  for (i in 2:n2) {
    gind <- na.omit(NNarray[i, ])
    temp <- lm(dat[, i] ~ -1 + dat[, gind])
    d <- 1/sd(resid(temp))
    uhat[i, i] <- d
    uhat[gind, i] <- -coef(temp) * d
  }
  return(uhat)
}

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