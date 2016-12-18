cd.general <- function(X, y, a0, beta, epsilon, max.iter, lambda, family, bInd, pf, 
    pentype, gamma) {
    # pf: penalty factor all all components, 0 means no penalty
    X = cbind(X, 1)
    p = ncol(X)
    n = nrow(X)
    beta = c(beta, a0)
    pf = c(pf, 0)
    para.in = c(epsilon, max.iter, lambda, gamma)
    # cat('beta[1:10]', beta[1:min(length(beta),10)],'\n')
    if (family == 'gaussian') {
        out = .Fortran("cd_general_lin", X = as.double(X), y = as.double(y), p = as.integer(p), 
            n = as.integer(n), beta = as.double(beta), pf = as.double(pf), 
            pentype = as.integer(pentype), paraIn = as.double(para.in))
    } else if (family == 'binomial') {
        out = .Fortran("cd_general_bin", X = as.double(X), y = as.double(y), p = as.integer(p), 
            n = as.integer(n), beta = as.double(beta), pf = as.double(pf), 
            pentype = as.integer(pentype), paraIn = as.double(para.in))
    }
    return(list(a0 = out$beta[p], beta = out$beta[-p]))
}
