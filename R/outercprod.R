outercprod <- function(x, y) {
    x = as.matrix(x)
    y = as.matrix(y)
    do.call(cbind, lapply(1:ncol(x), function(i) x[, i] * y))
}
