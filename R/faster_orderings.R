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
order_maximin_dist <- function(d){
  ## number of locations to order
  n = nrow(d)
  
  ## initialize ordering
  ord = numeric(length = n)
  #get first location (basically the most central point)
  ord[1] = which.min(rowSums(d))
  
  ## compute maxmin ordering
  #Get current column
  min <- d[,ord[1]]
  for(i in 2:n){
    #Get the maximum minimum distance
    a <- which.max(min)
    #Add it as the next point in the ordering
    ord[i] <- a
    #Update the minimums. This works because the current points distance from itself
    #is zero, so it gets ignored when taking the maximum in the future (as all distances
    #must be >= 0 to be a valid distance).
    min <- pmin(min, d[,a])
  }
  #returns the ordering
  return(ord)
}


#' Faster maximin ordering by locations
#' 
#' This function is a wrapper for \code{\link{order_maximin_dist}} to take in locations as
#' input and then calculate the whole distance matrix for ordering. The alternative is to use 
#' \code{\link{orderMaxMinFast}}, but that is much slower assuming the whole distance 
#' matrix fits in memory. However, the first point in the ordering is calculated differently.
#' This uses the smallest rowSum (so somewhat central to the point cloud), while
#' \code{\link{orderMaxMinFast}} uses the point closest to the average of locations.
#'
#' @param locs matrix of locations of points (to match input argument of rdist, see ?rdist)
#'
#' @return Maximin ordering of the locations
#' @export
#'
#' @examples
#' num_points <- 1000
#' locations <- matrix(runif(2 * num_points), nc = 2)
#' order_locs <- order_mm_locations(locations)
#' 
order_mm_locations <- function(locs){
  #Get distances
  d <- rdist(locs)
  #order using the ordering by distance function
  return(order_maximin_dist(d))
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
#' defaults to 1/2 the maximum distance.
#'
#' @return Maximin ordering of the points
#' @export
#'
#' @examples
order_mm_tapered <- function(locs, datum, tapering_range = 0.5){
  #Get distances between locations
  ds <- rdist(locs)
  exp_const <- Exponential(ds, range = (tapering_range * max(ds)))
  cov_matrix <- cov(datum) * exp_const
  #Covariance matrix to a distance matrix
  d = 1 - cov2cor(cov_matrix)
  return(order_maximin_dist(d))
}