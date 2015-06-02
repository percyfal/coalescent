#' watterson estimator
#'
#' @param int Sn number of segregating sites
#' @param int n number of samples
#'
#' @return theta_w
#'
watterson <- function(Sn, n) {
    Sn / sum(1/seq(1, n-1))
}
#' watterson estimator for several fragments
#'
#' @param int x number of fragments
#' @param int n sample size
#' @param int L fragment size
#' @param double theta per site mutation rate
#'
#' @return theta_w for each fragment
#'
estimate_watterson <- function(x=10, n=100, L=200, theta=0.01, ...) {
    d <- as.data.frame(
        do.call("rbind",
                lapply(seq(1, x),
                       function(y) {
                           c(fragment=y, coalescent(n=n, theta = theta*L, ...)$data)
                       }
                       ))
    )
    theta_w = tapply(d$nmut, d$fragment, watterson, n)
    cbind(d, theta=theta_w)
}
