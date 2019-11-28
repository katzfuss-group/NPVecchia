# These functions were not written by me! They are from another
# package in development that is not yet installable. Also, they 
# can be found in the BRISC package, but these are internal and 
# not callable from that package either.

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

# ordered by distance to some point loc0
#' @keywords internal
orderDist2Point <- function( locs, loc0 ){
    d <- ncol(locs)
    loc0 <- matrix(c(loc0),1,d)
    distvec <- rdist(locs,loc0)
    orderinds <- order(distvec)
}

# ordered by distance to the center
#' @keywords internal
orderMiddleOut <- function( locs ){
    d <- ncol(locs)
    loc0 <- matrix(colMeans(locs),1,d)
    orderDist2Point(locs,loc0)
}

#' @keywords internal
orderByCoordinate <- function( locs, coordinate ){
    # coordinate can be a single coordinate in {1,2,..,d}
    # in this case, the points are ordered increasingly
    # in that coordinate. if coordinate is a vector 
    # of coordinates, the points are ordered increasingly 
    # in the sum of the coordinates dicated by the vector 
    order(rowSums(locs[,coordinate,drop=FALSE]))
}

