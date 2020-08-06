
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NPVecchia

<!-- badges: start -->

<!-- badges: end -->

The goal of NPVecchia is to provide users with scalable Gaussian Process
(GP) covariance estimation without assuming its form. We assume that the
GP is de-trended (has zero mean) and that the covariance between points
decreases with some kind of distance. It provides the functions to
accompany an upcoming paper with Dr. Katzfuss.

**TL;DR** This method reorders the data and computes a sparse estimate Û
of the inverse of the Cholesky of the (reordered) covariance matrix such
that

![equation](https://latex.codecogs.com/png.latex?%5CSigma%5E%7B-1%7D%20%3D%20%5Chat%7BU%7D%20%5Chat%7BU%7D%27).

``` r
ans <- run_npvecchia(data, locations)
#output contains ans$u, ans$order
```

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("the-kidd17/NPVecchia")
```

## Ordering

It is recommended to use the provided function “orderMaxMinFaster” to
order the data in a max-min ordering based on Euclidean distance.
Basically, this ordering starts at a point in the center, then adds
points one at a time that maximize the minimum distance from all
previous points in the ordering. While not required, the decay of the
both the regression errors and coefficients (i.e. the prior) was
specifically developed based on this ordering. There is also the option
to order based on a user defined distance matrix (order\_maximin\_dist)
or based on the sample correlation matrix (order\_mm\_tapered).

``` r
#number of locations, dimensions, samples
n_locs <- 50
d <- 2
num_reps <- 5
#random locations and data
locs <- matrix(runif(d*n_locs), nc = d)
dataa <- matrix(rnorm(num_reps*n_locs), nc = n_locs)

#Compute maximin ordering
order1 <- orderMaxMinFaster(locs)
#Correlation ordering
order2 <- order_mm_tapered(locs, dataa)

#Reorder data by correlation ordering
dataa <- dataa[, order2]
```

## Neighbors

The best previous neighbors (ordered from nearest to furthest) are also
needed. As our method trims down the number of neighbors, it is
advisable to initially compute the full neighbor matrix (or at least a
large number of neighbors). find\_nn computes the neighbors based on
correlation (analogous to order\_mm\_tapered), and find\_nn\_dist is
similarly analogous to order\_maximin\_dist.

``` r
#Correlation neighbors
nearest_neighbors <- find_nn(locs, dataa, m = n_locs)
#Euclidean neighbors
nearest_neighbors2 <- find_nn_dist(fields::rdist(dataa), n_locs)
```

## Details

The method relies on three hyperparameters related to the decays of the
regression coefficients and the marginal variances. One can use Bayesian
methods, optimization, etc. to determine hyperparameters.

### Frequentist Approach

The easiest is to use ?optim to obtain good values of the
hyperparameters and then get the maximum a posteriori (MAP) value for Û.

``` r
#initial thetas
init_theta = c(1,-1,0)
#Optimization, where one can change limits to match ones desired decay
thetas_f <- optim(init_theta, minus_loglikeli_c, datum=dataa,
                  NNarray = nearest_neighbors, method="L-BFGS-B",
                  lower=-6, upper=4)

#Get MAP
uhat <- get_map(thetas_f$par, dataa, nearest_neighbors)
```

### Bayesian Approach

While no Metropolis sampler for the hyperparameters is provided,
adaptMCMC::MCMC was used for our experiments and seemed to work.

``` r
#A scale matrix similar to this worked well for our application
scale_mat <- matrix(c(0.05,-0.04,0,-0.04,0.05,0,0,0,0.01),nc=3)

#number of samples
nruns <- 100

#Run the sampler
samp <- adaptMCMC::MCMC(minus_loglikeli_c, datum = dataa, 
                        NNarray = nearest_neighbors,
                        init=init_theta, negativ=FALSE,
                        scale=scale_mat, adapt=TRUE, 
                        acc.rate=0.234, n= nruns)
#>   generate 100 samples

#Sample z ~ N(0, (UU')^{-1}) for each theta
sampled <- apply(samp$samples,1,function(thetas_temp){
  #Get priors from thetas
  ptt <- thetas_to_priors(thetas_temp, n_locs)
  #get number of neighbors
  m <- min(ncol(ptt[[3]]),ncol(nearest_neighbors))
  #Get posteriors
  posts <- get_posts_c(dataa, ptt, nearest_neighbors[, 1:m])
  #Get a sample 
  sample1 <- samp_posts(posts, nearest_neighbors[, 1:m], bayesian = T, uhat = F)
  #return
  sample1
})
```
