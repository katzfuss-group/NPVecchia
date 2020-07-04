# These functions were not written by me! They are from another
# package in development that is not yet installable. Also, they 
# can be found in the BRISC package, but these are internal and 
# not callable from that package.

#' Max-min ordering function
#' 
#' This is the O(n^2) algorithm for AMMD ordering. Start with a point
#' in the middle, then propose a random set of the remaining points
#' (of size 'numpropose') and choose the one that has maximum minimum
#' distance to the already selected points. set 'numpropose' to n
#' to get the exact maximum minimum distance ordering
#' 
#' @param locs matrix of locations
#' @param numpropose the number of locations to order
#' @return orderinds a vector to order the locations
#' @importFrom matrixStats colMins
#' @importFrom matrixStats rowMins
orderMaxMinFast <- function( locs, numpropose ){
    
    n <- nrow(locs)
    d <- ncol(locs)
    remaininginds <- 1:n
    orderinds <- rep(0L,n)
    # pick a center point
    mp <- matrix(colMeans(locs),1,d)
    distmp <- rdist(locs,mp)
    ordermp <- order(distmp)
    orderinds[1] = ordermp[1]
    remaininginds <- remaininginds[remaininginds!=orderinds[1]]
    for( j in 2:(n-1) ){
        randinds <- sample(remaininginds,min(numpropose,length(remaininginds)))
        distarray <-  rdist(locs[orderinds[1:j-1],,drop=FALSE],locs[randinds,,drop=FALSE])
        bestind <- which(colMins(distarray) ==  max( colMins( distarray ) )) 
        orderinds[j] <- randinds[bestind[1]]    
        remaininginds <- remaininginds[remaininginds!=orderinds[j]]
    }    
    orderinds[n] <- remaininginds
    orderinds
}

#' Order by distance to some point loc0
#' @param locs locations to order by distance to point loc0
#'
#' @param loc0 reference point for ordering locations
#'
#' @keywords internal
orderDist2Point <- function( locs, loc0 ){
    d <- ncol(locs)
    loc0 <- matrix(c(loc0),1,d)
    distvec <- rdist(locs,loc0)
    orderinds <- order(distvec)
}

#' Order by distance to the center
#' 
#' @param locs locations to order
#'
#' @keywords internal
orderMiddleOut <- function( locs ){
    d <- ncol(locs)
    loc0 <- matrix(colMeans(locs),1,d)
    orderDist2Point(locs,loc0)
}

#' Order by one a single coordinate
#' 
#' @param locs locations to order
#'
#' @param coordinate index of the coordinate to order by
#'
#' @keywords internal
orderByCoordinate <- function( locs, coordinate ){
    # coordinate can be a single coordinate in {1,2,..,d}
    # in this case, the points are ordered increasingly
    # in that coordinate. if coordinate is a vector 
    # of coordinates, the points are ordered increasingly 
    # in the sum of the coordinates dicated by the vector 
    order(rowSums(locs[,coordinate,drop=FALSE]))
}

#' Find Nearest neighbors
#'
#' @param cov.matrix an (ordered) covariance matrix to use for finding NN based on
#' this distance (1 - correlation)
#' @param m number of (ordered) nearest neighbors to calculate
#'
#' @return a matrix of nearest neighbors of dimension n x m
#' @export
#'
#' @examples
find_nn <- function(cov.matrix,m)
{
    
    ## convert covariance matrix to correlation-based distance
    d=1-cov2cor(cov.matrix)
    n=nrow(cov.matrix)
    
    ## find ordered NN
    #initialize
    NN=matrix(NA,n,m)
    for(i in 2:n){
        # if((i %% 1000)==0) print(i)
        #get number of neighbors, m if that many previous points
        k=min(i-1,m)
        NN[i,1:k]=order(d[i,1:(i-1)])[1:k]
    }
    
    return(NN)
}
