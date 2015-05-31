#' Coalescent simulation
#'
#' @param n number of samples
#' @param theta per generation mutation rate
#' @param newick logical; return newick string representation of tree
#'
#' @return vector containing TMRCA, TBL, S_n, and possibly newick
#'     representation of tree
#'
coalescent <- function(n=10, theta=1.0, newick=FALSE) {
    l <- .Call("tree", n, theta, newick)
    return (as.list(c(tmrca=l[[1]], tbl=l[[2]], nmut=l[[3]], newick=l[[4]])))
}
