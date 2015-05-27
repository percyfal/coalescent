#' Coalescent simulation
#'
#' @param int n number of samples
#' @param double theta per generation mutation rate
#' @param int L fragment size
#'
#' @return vector containing TMRCA, TBL, S_n
#'
coalescent <- function(n=10, theta=1.0, L=0) {
    l <- .Call("tree", n, theta, L)
    return (c(tmrca=l[[1]], tbl=l[[2]], nmut=l[[3]]))
}
