#' Coalescent simulation
#'
#' @param n number of samples
#' @param theta per generation mutation rate
#' @param newick logical; return newick string representation of tree
#' @param Tb time to "bottleneck" event (can also model expansion; see
#'     parameter f)
#' @param f fraction of population size at Tb compared to current size
#'     N_0; in other words, N_b = f * N_0. When f>1 a bottleneck is
#'     modelled, when f<1 an expansion
#'
#' @return vector containing TMRCA, TBL, S_n, and possibly newick
#'     representation of tree
#'
coalescent <- function(n=10, theta=1.0, newick=FALSE, Tb=1.0, f=1.0) {
    l <- .Call("tree", n, theta, newick, Tb, f)
    return (list(data=c(tmrca=l[[1]],
                        tbl=l[[2]],
                        nmut=l[[3]]),
                 newick=as.character(l[[4]]))
            )
}
