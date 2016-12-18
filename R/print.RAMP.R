#' Result summary of a fitted RAMP object.
#' 
#' Similar to the usual print methods, this function summarize results 
#' from a fitted \code{'RAMP'} object.
#' @export
#' @param x Fitted \code{'RAMP'} model object.
#' @param  digits  The number of significant digits for the coefficient estimates.
#' @param \dots Not used. Other arguments to predict. 
#' @return No value is returned.
#' @seealso \code{\link{RAMP}}
print.RAMP = function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("Important main effects:", x$mainInd, "\n")
    cat("Coefficient estimates for main effects:", signif(x$beta.m, digits), "\n")
    cat("Important interaction effects:", x$interInd, "\n")
    if (length(x$interInd) > 1) {
        cat("Coefficient estimates for interaction effects:", signif(x$beta.i, digits), 
            "\n")
    }
    cat("Intercept estimate:", signif(x$a0, digits), "\n")
}
