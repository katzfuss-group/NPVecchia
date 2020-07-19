#' Maximin ordering based on a distance
#'
#' This function orders the points in a maximin ordering based on distance.
#'
#' @param d distance matrix (e.g. distances)
#'
#' @return Maximin ordering of the points corresponding to the matrix
#' @export
#'
#' @examples
#' n <- 100
#' fake_dists <- matrix(runif(n^2), nc = n)
#' #make it a valid distance matrix by making d(a,a) = 0
#' diag(fake_dists) <- 0
#' order_mm <- order_maximin_dist(fake_dists)
#' #reorder data
order_maximin_dist <- function(d){
  ## number of locations to order
  n = nrow(d)
  
  ## initialize ordering
  ord = numeric(length = n)
  #get first location (basically the most central point)
  ord[1] = which.min(rowSums(d))
  
  ## compute maxmin ordering
  #Get current column
  mint <- d[,ord[1]]
  for(i in 2:n){
    #Get the maximum minimum distance
    a <- which.max(mint)
    #Add it as the next point in the ordering
    ord[i] <- a
    #Update the minimums. This works because the current points distance from itself
    #is zero, so it gets ignored when taking the maximum in the future (as all distances
    #must be >= 0 to be a valid distance).
    mint <- pmin(mint, d[,a])
  }
  #returns the ordering
  return(ord)
}


#' Maximin ordering using a tapered covariance matrix
#'
#' This function is a wrapper for \code{\link{order_maximin_dist}} to take in data and
#' locations to calculate the ordering based on a tapered sample covariance matrix. It
#' provides an argument to control the amount of tapering. 
#' 
#' @param locs matrix of locations of points (to match input argument of fields::rdist)
#' @param datum Data where column i corresponds to observations at the i'th point of locs
#' @param tapering_range Percentage of the maximum distance for Exponential tapering, which
#' defaults to 0.4 * the maximum distance.
#'
#' @return Maximin ordering of the points
#' @export
#'
#' @examples
#' n <- 100
#' d <- 2
#' N <- 50
#' locationss <- matrix(runif(n * d), nc = d)
#' dataa <- matrix(rnorm(n * N), nr = N)
#' order_mm <- order_mm_tapered(locationss, dataa)
#' #reorder locations/data
#' locationss <- locationss[order_mm, ]
#' dataa <- dataa[, order_mm]
order_mm_tapered <- function(locs, datum, tapering_range = 0.4){
  #Get distances between locations
  ds <- rdist(locs)
  exp_const <- Exponential(ds, range = (tapering_range * max(ds)))
  cov_matrix <- cov(datum) * exp_const
  #Covariance matrix to a distance matrix
  d = 1 - cov2cor(cov_matrix)
  return(order_maximin_dist(d))
}

#' Faster maximin ordering by Euclidean distance
#' 
#' Assumes Euclidean distance and finds maximin ordering. It is equivalent to \code{\link{order_maximin_dist}}
#' but much faster.
#'
#' @param locs matrix of n locations in arbitrary dimension (to match input argument of fields::rdist)
#' @param low_mem flag that changes the memory requirements of the function. Defaults to false, which allows
#' the function to calculate the full distance matrix between all locations. If set to true, it will only calculate
#' one column of the distance matrix at a time, which slows down the function but also decreases memory requirements.
#'
#' @return Maximin ordering of the points
#' @export
#'
#' @examples
#' #' n <- 100
#' d <- 3
#' locationss <- matrix(runif(n * d), nc = d)
#' ordering <- orderMaxMinFaster(locationss)
orderMaxMinFaster <- function(locs, low_mem = FALSE){
  ## number of locations to order and number of dimensions
  n <- nrow(locs)
  dimens <- ncol(locs)
  
  ## initialize ordering
  ord = numeric(length = n)
  #get first location (point closest to center)
  mp <- matrix(colMeans(locs),1,dimens)
  distmp <- rdist(locs,mp)
  ord[1] = which.min(distmp)
  
  #get distances if low_mem = FALSE
  if(!low_mem){
    d = rdist(locs)
    mint <- d[,ord[1]]
    ## compute maxmin ordering
    #Get current column
    for(i in 2:n){
      #Get the maximum minimum distance
      a <- which.max(mint)
      #Add it as the next point in the ordering
      ord[i] <- a
      #Update the minimums. This works because the current points distance from itself
      #is zero, so it gets ignored when taking the maximum in the future (as all distances
      #must be >= 0 to be a valid distance).
      mint <- pmin(mint, d[,a])
    }
  }else{
    #Get current minimums
    mint <- rdist(locs, t(locs[ord[1],]))
    for(i in 2:n){
      #Get the maximum minimum distance
      a <- which.max(mint)
      #Add it as the next point in the ordering
      ord[i] <- a
      #Get distances from next point in ordering
      d_temp = rdist(locs, t(locs[a,]))
      #Update the minimums. This works because the current points distance from itself
      #is zero, so it gets ignored when taking the maximum in the future (as all distances
      #must be >= 0 to be a valid distance).
      mint = pmin(mint, d_temp)
    }
  }
  #returns the ordering
  return(ord)
}
