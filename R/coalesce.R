#' Coalescent simulation
#'
#' @param int n number of samples
#' @param double theta mutation rate
#' 
#' @return tree coalescent tree
#'
coalescent <- function(n=10, theta=1.0) {
    l <- .Call("tree", n, theta)
    return (list(tmrca=l[[1]], tbl=l[[2]], nmut=l[[3]]))
}
