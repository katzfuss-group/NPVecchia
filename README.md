
<!-- README.md is generated from README.Rmd. Please edit that file -->
NPVecchia
=========

<!-- badges: start -->
<!-- badges: end -->
The goal of NPVecchia is to provide users with scalable GP covariance estimation without assuming its form. It provides some of the barebones functionality to accompany an upcoming paper with Dr. Katzfuss.

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("the-kidd17/NPVecchia")
```

Further Details
===============

This package was created specifically for spatial data, as the locations are known and in a low dimensional (usually 2 or 3) Euclidean space. If one already has data, it is still recommended to read the next section to understand how to order the locations and get the nearest neighbors.

Many of the functions have C++ counterparts by adding "\_c" to their names. However, there is no error checking, so errors that occur tend to be difficult to interpret and fix. That is the big use for the R versions: to understand the functions and usage better (with error-checking) before worrying about speed.

Data Creation
-------------

As a first step to experimenting with the package, data must be simulated. These examples will be on a 2-D grid, though the method generalizes to not gridded data and probably to somewhat higher dimensions. We use a unit square as the domain for simplicity; by changing the range (also known as scale) parameter in spatial covariance functions, the data can similarly be rescaled to a unit square.

It is recommended to use the provided function "orderMaxMinFast" to order the data in a max-min ordering. Basically, this ordering starts at a point in the center, then adds points one at a time that maximize the minimum distance from all previous points in the ordering. While not required, the decay of the both the regression errors and coefficients (i.e. the prior) is implicitly based on this ordering.

``` r
#seed for reproducibility
set.seed(1128)
#number of points
n <- 30^2
#Get the grid in 1-D
grd <- seq(0, 1, length.out = sqrt(n))
#expand to 2-D grid
grd2d <- expand.grid(grd, grd)
#Alternative for random (non-gridded) data
#grd2d <- matrix(runif(2 * n), ncol = 2)

#Order the grid points
order_mmd <- orderMaxMinFast(grd2d, n^2)
grd2d_ordered <- grd2d[order_mmd, ]
```

Now, with the locations, a covariance function is needed. As is common in spatial, this example will create a Matern covariance matrix. The package doesn't rely on functions for using other covariance matrices, but two useful packages for other covariance functions are 'geoR' and 'RandomFields'.

Also needed are the nearest neighbors for each point. This only needs to be done once, so it is useful to allow the number of neighbors to be large (it will often be decreased automatically by the method). This outputs a matrix where each row corresponds to that points' nearest neighbors (ordered from closest to furthest; NAs fill other spots). Our methodology does not include each point as a neighbor of itself, so the first column (that gives each point) is removed.

**Note**: if using an Anisotropic and/or non-stationary covariance, the ordering and/or neighborhood selection will not be optimal without accounting for it. Currently, these just use the spatial locations and are sub-optimal, though they are still better than random or coordinate-based orderings. Newer versions may include improved functionality for this.

``` r
#get nearest neighbors
NNarray <- GpGp::find_ordered_nn(grd2d_ordered, 40)
#Remove points as own neighbors
NNarray <- NNarray[, -1]
#get distances between points
distances <- fields::rdist(grd2d_ordered)
# #If one wants anisotropic distances instead
# mahala_dist <- function(a, b, Aniso_matrix) {
#   sqrt(c(t(a - b) %*% solve(Aniso_matrix, a - b)))
# }
# mahala_dist_vect <- Vectorize(mahala_dist)
# grd2d_ord_list <- as.list(data.frame(t(grd2d_ordered)))
# dists <- outer(grd2d_ord_list, grd2d_ord_list, mahala_dist_vect, 
#                Aniso_matrix = matrix(c(0.5, 0, 0, 1), ncol = 2))
#Specify parameters of the covariance matrix
marginal_variance <- 3
#for unit square, range < 0.5 is recommended to avoid numerical issues
range <- 0.25
#smoothness should similarly not be too large
smoothness <- 1
#create the true covariance matrix
covar_true <- marginal_variance * fields::Matern(distances, range = range, 
                                         nu = smoothness)
#Create data (can use whatever way one desires)
#Number of replications
N <- 50
datum <- MASS::mvrnorm(N, rep(0, nrow(covar_true)), covar_true)
```

Methodology
-----------

To keep the explanation simple, we estimate the Cholesky of the precision matrix using a series of simple linear regressions. The twist is that we only regress on the nearest neighbors spatially to drastically improve the computational cost. The frequentist version is "get\_mle", which performs the regressions with a known sparsity pattern to find the sparse matrix $\\hat{U}$ such that

$$
\\Sigma^{-1} = \\hat{U} \\hat{U}'
$$

``` r
uhat_mle <- get_mle(datum, NNarray)
#Or to decreate the number of neighbors
uhat_mle <- get_mle(datum, NNarray[, 1:9])
```

The Bayesian version is more involved, but adds a prior to further regularize the non-zero elements. See the [Math vignette](documents/math.pdf) for the mathematical details of the method. The priors on all of the regression rely on 3 hyperparameters that must be optimized for the best results. To avoid numerical instabilities, it is recommended to force the hyperparameters to be in the range $\[-6,4\]. One way to find good hyperparameters is to optimize the integrated log-likelihood as shown below.

**Note:** Another alternative is to use MCMC to get a distribution on the hyperparameters. This is very sensitive to inputs, is often not computationally feasible, and is not included. If one is set on MCMC, "adaptMCMC" is an easy to use package that includes an adaptive Metropolis-Hastings algorithm. It improves mixing and adaptively accounts for our correlated hyperparameters.

``` r
#define initialize thetas (starting values)
thetas_init <- c(1, -1, 0)

thetas_best <- optim(thetas_init, minus_loglikeli_c, datum = datum,
                     NNarray = NNarray, method = "L-BFGS-B",
                     lower = -6, upper = 4)
show(thetas_best)
#> $par
#> [1]  2.138334 -2.042955 -0.292671
#> 
#> $value
#> [1] 28264.86
#> 
#> $counts
#> function gradient 
#>       18       18 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
```

With these optimal hyperparameters, it is straightforward to get $\\hat{U}$.

``` r
#get priors
priors_final <- thetas_to_priors(thetas_best$par, n)
show(paste("The method found", ncol(priors_final[[3]]), "neighbors to be sufficient."))
#> [1] "The method found 9 neighbors to be sufficient."
#get posteriors
posteriors <- get_posts(datum, priors_final, NNarray)
#get uhat
uhat <- samp_posts(posteriors, NNarray)
```

Error metrics
-------------

To see the difference, we can evaluate Stein's loss (the exclusive Kullback-Leibler divergence). It is defined as

$$
KL(\\hat{\\Sigma} || \\Sigma) = tr(\\hat{\\Sigma} \\Sigma^{-1}) - log |\\hat{\\Sigma}\\Sigma^{-1}| - n
$$

When the number of neighbors is high, the "MLE" version of $\\hat{U}$ is not invertible, so some care is needed in choosing *m* &lt; *N* for the frequentist method. Alternatively, the regressions can be replaced with some sort of lasso to ensure smaller *m*.

``` r
covar_true_inv <- solve(covar_true)
get_kl <- function(uhat, covar_inv) {
  # get sigma hat
  sigma_hat <- crossprod(solve(uhat, tol= 1e-35))
  # get sigma hat * sigma inv 
  cov2 <- sigma_hat %*% covar_inv
  # get its determinant
  cov_det <- determinant(cov2)
  return((sum(diag(cov2)) - as.numeric(cov_det$modulus) - n))
}
show(paste("The KL-divergence for our method is", get_kl(uhat, covar_true_inv),
           " while it is", get_kl(uhat_mle, covar_true_inv), "for the MLE."))
#> [1] "The KL-divergence for our method is 194.505371978698  while it is 286.685939999608 for the MLE."
```

Another alternative to avoid having to invert the Cholesky is to use the frobenius norm (or other metric, such as singular value differences) to compare the estamates to the true precision matrix.

``` r
precision_frobenius <- function(uhat, covar_inv) {
  # get estimated precision
  precision_hat <- Matrix::tcrossprod(uhat)
  sqrt(sum((precision_hat - covar_inv)^2))
}
show(paste("The Frobenius norm for our method is", precision_frobenius(uhat, covar_true_inv),
           " while it is", precision_frobenius(uhat_mle, covar_true_inv), "for the MLE."))
#> [1] "The Frobenius norm for our method is 534.510013613234  while it is 774.209517038395 for the MLE."
```
