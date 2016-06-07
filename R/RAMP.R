#' Regularization Algorithm under Marginality Principle (RAMP) for high dimensional generalized quadratic regression.
#'
#' Regularization Algorithm under Marginality Principle (RAMP) for high dimensional generalized quadratic regression.
#' @useDynLib RAMP
#' @importFrom stats coef
#' @importFrom stats glm
#' @importFrom stats var
#' @export
#' @param X input matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y  response variable, of dimension nobs x 1. non-negative counts for
#' \code{family='poisson'}, binary for \code{family='binomial'}.
#' @param family response type. Default is 'gaussian'.
#' @param penalty Choose from \code{LASSO}, \code{SCAD} and \code{MCP}. Default is 'LASSO'. 
#' @param gamma concavity parameter. If missing, the code will use 3.7 for 'SCAD' and 2.7 for 'MCP'. 
#' @param inter whether to select interaction effects. Default is TRUE.
#' @param hier whether to enforce strong or weak heredity. Default is 'Strong'.
#' @param eps the precision used to test the convergence. Default is 1e-15.
#' @param tune tuning parameter selection method.
#' 'AIC', 'BIC', 'EBIC' and 'GIC' are available options. Default is EBIC.
#' @param lam.list a user supplied \eqn{\lambda} sequence.
#' typical usage is to have the program compute its own
#' \code{lambda} sequence based on \code{lambda.min.ratio} and \code{n.lambda}.
#'  supplying a value of \eqn{\lambda} overrides this.
#'  @param lambda.min.ratio optional input. smallest value for \code{lambda}, as a fraction of \code{max.lam}, the (data derived) entry value. the default depends on the sample size \code{n} relative to the number of variables \code{p}. if \code{n} > \code{p}, the default is 0.0001. otherwise, the default is 0.01.
#' @param max.iter maximum number of iteration in the computation. Default is 100.
#' @param max.num optional input. maximum number of nonzero coefficients.
#' @param n.lambda the number of \code{lambda} values. Default is 100.
#' @param ebic.gamma  the gamma parameter value in the EBIC criteria. Default is  1.
#' @param refit   whether to perform a MLE refit on the selected model. Default is TRUE.
#' @param trace  whether to trace the fitting process. Default is FALSE.
#' @return An object with S3 class RAMP.
#' \item{a0}{intercept vector of length(\code{lambda}).}
#' \item{mainInd}{index for the selected main effects.}
#' \item{interInd}{index for the selected interaction effects}
#' \item{beta.m}{coefficients for the selected main effects.}
#' \item{beta.i}{coefficients for the selected interaction effects.}
#' @seealso \code{\link{predict.RAMP}},\code{\link{print.RAMP}}
#' @examples
#' set.seed(0)
#' n = 500
#' p = 10 #Can be changed to a much larger number say 50000
#' x = matrix(rnorm(n*p),n,p)
#' eta = 1 * x[,1] + 2 * x[,3]  + 3*x[,6]  + 4*x[,1]*x[,3] + 5*x[,1]*x[,6]
#' y =  eta + rnorm(n)
#' xtest = matrix(rnorm(n*p),n,p) 
#' eta.test = 1 * xtest[,1] + 2 * xtest[,3]  + 3*xtest[,6] + 
#' 4*xtest[,1]*xtest[,3] + 5*xtest[,1]*xtest[,6]
#' ytest =  eta.test + rnorm(n)

#' fit1 = RAMP(x, y)
#' fit1    ###examine the results
#' ypred = predict(fit1, xtest)
#' mean((ypred-ytest)^2)
#' 
#' #fit1.scad = RAMP(x, y, penalty = 'SCAD')
#' #fit1.scad    ###examine the results
#' 
#' #fit1.mcp = RAMP(x, y, penalty = 'MCP')
#' #fit1.mcp    ###examine the results
#' 
#' ##Now, try a binary response
#' #y = rbinom(n, 1, 1/(1+exp(-eta)))
#' #fit2 = RAMP(x, y, family='binomial')  ###for binary response
#' 
#' ## Weak heredity
#' eta = 1 * x[,1] + 3*x[,6]  + 4*x[,1]*x[,3] + 5*x[,1]*x[,6]
#' y =  eta + rnorm(n)
#' eta.test = 1 * xtest[,1] +  3*xtest[,6] + 4*xtest[,1]*xtest[,3] + 
#' 5*xtest[,1]*xtest[,6]
#' ytest =  eta.test + rnorm(n)

#' 
#' fit3 = RAMP(x, y, hier = 'Strong')
#' fit3    ###examine the results
#' ypred3 = predict(fit3, xtest)
#' mean((ypred3-ytest)^2)
#' fit4 = RAMP(x, y, hier = 'Weak')
#' fit4   
#' ypred4 = predict(fit4, xtest)
#' mean((ypred4-ytest)^2)
RAMP <- function(X, y, family = "gaussian", penalty = "LASSO", gamma = 3.7, inter = TRUE, 
    hier = "Strong", eps = 1e-15, tune = "EBIC", lam.list, lambda.min.ratio, max.iter = 100, 
    max.num, n.lambda = 100, ebic.gamma = 1, refit = TRUE, trace = FALSE) {
    ## hier = 'Strong', or 'Weak', strong or weak heredity.
    if (penalty == "SCAD" & missing(gamma)) {
        gamma = 3.7
    }
    if (penalty == "MCP" & missing(gamma)) {
        gamma = 2.7
    }
    if (penalty == "LASSO" & !missing(gamma)) {
        stop("gamma is not needed")
    }
    pentype = switch(penalty, LASSO = 1, MCP = 2, SCAD = 3)
    
    n = dim(X)[1]
    p = dim(X)[2]
    oldX = X
    if (missing(max.num)) 
        max.num = p + 1 else max.num = min(max.num, p + 1)
    
    if (missing(lambda.min.ratio)) 
        lambda.min.ratio = ifelse(n < p, 0.01, 1e-04)
    
    ################ prepare variables
    lambda = NULL
    beta.mat = NULL  ##storing the beta coefficents along the path
    beta = rep(0, p)  ##the current beta estimates including the interation terms
    
    index = NULL
    
    
    ################ standardize design matrix
    X = scale(X)
    xm0 = attr(X, "scaled:center")
    xsd0 = attr(X, "scaled:scale")
    
    
    ############# lambda list
    if (family == "binomial") {
        max.lam = max(abs((1/n) * (t(X) %*% (y - mean(y)))))
        a0list = rep(log(mean(y)/(1 - mean(y))), n.lambda)
    }
    if (family == "poisson") {
        max.lam = max(abs((1/n) * (t(X) %*% (y - mean(y)))))
        a0list = rep(log(mean(y)), n.lambda)
        
    }
    
    if (family == "gaussian") {
        max.lam = max(abs((1/n) * (t(X) %*% (y - mean(y)))))
        a0list = rep(mean(y), n.lambda)
        
    }
    
    if (missing(lam.list)) {
        min.lam = max.lam * lambda.min.ratio
        lam.list = exp(seq(from = log(max.lam), to = log(min.lam), length.out = n.lambda))
    } else {
        lam.list = lam.list[lam.list <= max.lam]
        n.lambda = length(lam.list)
    }
    
    
    
    
    a0 = a0list[1]
    # loglik.list = rep(0, n.lambda) cri.list = rep(0,n.lambda) AIC.list =
    # rep(0,n.lambda) BIC.list = rep(0,n.lambda) EBIC.list = rep(0,n.lambda)
    # GIC.list = rep(0,n.lambda) df.list = rep(0,n.lambda) df.m.list =
    # rep(0,n.lambda) df.i.list = rep(0,n.lambda)
    loglik.list = NULL
    cri.list = NULL
    AIC.list = NULL
    BIC.list = NULL
    EBIC.list = NULL
    GIC.list = NULL
    df.list = NULL
    df.m.list = NULL
    df.i.list = NULL
    ind.list.inter = as.list(NULL)
    ind.list.main = as.list(NULL)
    ind.list.inter.xlab = as.list(NULL)
    beta.list.inter = as.list(NULL)
    beta.list.inter[[1]] = NULL
    ind.list.inter[[1]] = NULL
    ind.list.main[[1]] = NULL
    ind.list.inter.xlab[[1]] = NULL
    
    ############## main part
    
    
    ## expand X matrix dynamically each step by incorporating the new candidates of
    ## interactions, which are the interactions of the current active variables minus
    ## the current interaction effects, the interactions will play the same role as
    ## the main effect, except we keep track of the list through the lambda sequence.
    ## Keep track of other indices as well.colnames(X)
    
    colnames(X) = paste("X", 1:p, sep = "")
    
    nonPen = rep(0, p)
    for (k in 1:n.lambda) {
        ind.list.inter.xlab[[k]] = "None"
        # if(k==10)browser() print(k)
        old.beta = beta  ##back up the old beta
        nonPen = rep(0, p)
        
        break.ind = FALSE
        
        
        lam = lam.list[k]
        
        
        index = which(abs(beta[1:p]) > eps)  ###selected main effects from last step
        if (length(index) > 0) {
            ### find the candidate interation effects following strong heredity strong
            ### heredity
            aa = outer(index, index, f <- function(x, y) {
                paste("X", pmin(x, y), "X", pmax(x, y), sep = "")
            })
            aa = as.vector(aa)
            newinter = outercprod(X[, index], X[, index])
            bb = unique(as.vector(aa))
            loc = match(bb, aa)
            newinter = newinter[, loc]
            # newinter = as.matrix(newinter) colnames(newinter) = bb ##newinter contain the
            # candidate interaction matrix with the colnames the labels.
            
            if (hier == "Weak") {
                ### append strong heredity case to form the weak ones.
                rindex = setdiff(1:p, index)
                if (length(rindex) > 0) {
                  ## if there are candidate interactions under weak heredity
                  naa = outer(index, rindex, f <- function(x, y) {
                    paste("X", pmin(x, y), "X", pmax(x, y), sep = "")
                  })
                  naa = as.vector(naa)
                  nnewinter = outercprod(X[, rindex], X[, index])
                  # colnames(nnewinter) = naa
                  newinter = cbind(newinter, nnewinter)
                  bb = c(bb, naa)
                }
                
            }
            
            
            curInter = colnames(X)[-(1:p)]
            candInter = setdiff(bb, curInter)
            curloc = match(candInter, bb)
            newinter = as.matrix(newinter)
            newinter = newinter[, curloc]
            ncurInter = length(curInter)
            ncandInter = length(candInter)
            
            
            # cat('k=',k,'candidate interaction', candInter, '\n', sep=' ')
            if (ncurInter > 0 & hier == "Strong") {
                ## active interaction terms, setting the penalty for the parents to 0
                for (indInter in 1:ncurInter) {
                  pair = as.numeric(strsplit(curInter[indInter], "X")[[1]][2:3])
                  nonPen[pair[1]] = 1
                  nonPen[pair[2]] = 1
                }
            }
            # nonPen[index] = 1 Xinter = NULL
            if (ncandInter > 0 && inter) {
                ## candidate interaction terms
                
                
                # for (indInter in 1:ncandInter) { pair =
                # as.numeric(strsplit(candInter[indInter], 'X')[[1]][2:3]) Xinter =
                # cbind(Xinter, X[, pair[1]] * X[, pair[2]])
                
                
                # }
                Xnew = cbind(X, newinter)
                # allxm = c(xm,xinterm) allxvar = c(xvar,xintervar)
                
                colnames(Xnew) = c(colnames(X), candInter)
                X = Xnew
                
                beta = c(beta, rep(0, ncandInter))  #expand the beta coefficent vector to include the candiate interaction terms.
            }
        }
        X[, 1:p] = oldX
        X = scale(X)
        xm = attr(X, "scaled:center")
        xvar = attr(X, "scaled:scale")
        # xm[1:p] = m.m xvar[1:p] = m.var
        for (ite in 1:max.iter) {
            
            if (ite == 1) {
                cd.temp1 = cd.general(X = X, y = y, a0 = a0, beta = beta, epsilon = eps, 
                  max.iter = 1, lambda = lam, family = family, bInd = break.ind, 
                  nonPen = nonPen, pentype = pentype, gamma = gamma)
                a0 = cd.temp1$a0
                beta = cd.temp1$beta
                ind1 = which(abs(beta) > eps)
            }
            
            
            # #CD
            if (ite > 1) 
                {
                  cd.temp2 = cd.general(X = X[, ind1], y = y, a0 = a0, beta = beta[ind1], 
                    epsilon = eps, max.iter = max.iter, lambda = lam, family = family, 
                    bInd = break.ind, nonPen = nonPen[ind1[ind1 <= p]], pentype = pentype, 
                    gamma = gamma)
                  a0list[k] = cd.temp2$a0
                  beta2 = beta
                  beta2[ind1] = cd.temp2$beta
                  
                  
                  # ########redetect active set check.beta=new.beta
                  cd.temp3 = cd.general(X = X, y = y, a0 = a0, beta = beta2, epsilon = eps, 
                    max.iter = 1, lambda = lam, family = family, bInd = break.ind, 
                    nonPen = nonPen, pentype = pentype, gamma = gamma)
                  a0list[k] = cd.temp3$a0
                  beta = cd.temp3$beta
                  ind3 = which(abs(beta) > eps)
                  if (setequal(ind1, ind3)) {
                    break
                  }
                  ind1 = ind3
                  
                }  ##END iter>1
            
            
            
            
            # end check-cd-recheck
        }  ##END iter
        
        # cat('k=',k,'iter num',ite,'\n',sep=' ')
        
        # }
        
        # end case if sparse enough }
        
        ind.list.main[[k]] = which(abs(beta[1:p]) > eps)
        ind.list.inter[[k]] = which(abs(beta[-(1:p)]) > eps)  #record the interaction index pair location list
        beta.ols = beta
        size.main = length(ind.list.main[[k]])
        size.inter = length(ind.list.inter[[k]])
        
        index = which(abs(beta) > eps)
        beta[-index] = 0
        if (size.inter > 0 & hier == "Strong") {
            ### if interaction effects are detected, enforce the strong heredity in this step
            tmpindname = colnames(X)[ind.list.inter[[k]] + p]
            cur.main.ind = intertomain(tmpindname)
            if (trace == T) 
                cat("k=", k, "Enforced main effects", setdiff(cur.main.ind, ind.list.main[[k]]), 
                  "\n")
            
            ind.list.main[[k]] = union(ind.list.main[[k]], cur.main.ind)
            index = union(index, cur.main.ind)
        }
        size.main = length(ind.list.main[[k]])
        
        index = sort(index)
        
        df.list[k] = size.main + size.inter
        
        refit.beta = beta
        refit.a0 = a0
        loglik.list[k] = loglik(X, y, beta, family)
        if (df.list[k] > 0 && refit == TRUE && length(index) < n/2) {
            ## update the estimate using the result from an OLS refit.
            ols.fit = ols.refit(y, oldX, index[index <= p], colnames(X)[ind.list.inter[[k]] + 
                p], family = family)
            
            refit.beta[index] = ols.fit$beta[-1]
            refit.a0 = ols.fit$beta[1]
            loglik.list[k] = ols.fit$loglik
            
            if (family == "binomial" && min(abs(ols.fit$link)) > 10) 
                break
            
        }
        
        if (length(beta) > p) {
            tmp = which(abs(beta[-(1:p)]) > eps)
            beta = beta[c(1:p, p + tmp)]  ##update beta and X with only the main effect + important interaction effect
            refit.beta = refit.beta[c(1:p, p + tmp)]
            
            X = X[, c(1:p, p + tmp)]
            # ind.list.inter.xlab[[k]] = NULL
            if (length(tmp) > 0) {
                ## if interaction effects are selected.
                ind.list.inter.xlab[[k]] = colnames(X)[-(1:p)]
                beta.list.inter[[k]] = refit.beta[-(1:p)]
            } else {
                ind.list.inter.xlab[[k]] = "None"
            }
        }
        
        a0list[k] = refit.a0
        beta.mat = cbind(beta.mat, refit.beta[1:p])
        # break criteria
        
        
        
        AIC.list[k] = loglik.list[k] + 2 * length(index)
        BIC.list[k] = loglik.list[k] + log(n) * length(index)
        if (inter == FALSE) {
            p.eff = p
        } else if (hier == "Strong") {
            p.eff = p + size.main^2
        } else if (hier == "Weak") {
            p.eff = p + size.main * p
        }
        EBIC.list[k] = loglik.list[k] + log(n) * length(index) + 2 * ebic.gamma * 
            log(choose(p.eff, size.main + size.inter))
        GIC.list[k] = loglik.list[k] + log(log(n)) * log(p.eff) * length(index)
        if (tune == "AIC") 
            cri.list[k] = AIC.list[k]
        if (tune == "BIC") 
            cri.list[k] = BIC.list[k]
        if (tune == "EBIC") 
            cri.list[k] = EBIC.list[k]
        if (tune == "GIC") 
            cri.list[k] = GIC.list[k]
        
        df.m.list[k] = length(ind.list.main[[k]])
        df.i.list[k] = length(beta) - p
        
        # if (sum(beta != 0) >= min(n - 1,100))
        if (sum(beta != 0) >= n - 1) 
            break
        if (break.ind == TRUE) {
            print("Warning: Algorithm failed to converge for all values of lambda")
            break
        }
        
        
        if (trace == T) {
            cat("k=", k, "current main effects", ind.list.main[[k]], "\n", sep = " ")
            cat("k=", k, "current interaction", ind.list.inter.xlab[[k]], "\n", 
                sep = " ")
        }
        # end the outer loop: decrease in lambda
    }
    cri.loc = which.min(cri.list)
    AIC.loc = which.min(AIC.list)
    BIC.loc = which.min(BIC.list)
    EBIC.loc = which.min(EBIC.list)
    GIC.loc = which.min(GIC.list)
    all.locs = c(AIC.loc, BIC.loc, EBIC.loc, GIC.loc)
    
    lambda = lam.list[1:ncol(beta.mat)]
    # print(cri.loc) browser()
    if (length(ind.list.inter) == 0) {
        interInd = NULL
    } else {
        interInd = ind.list.inter.xlab[[cri.loc]]
    }
    if (length(beta.list.inter) == 0) {
        beta.i = NULL
    } else {
        beta.i = beta.list.inter[[cri.loc]]
    }
    val = list(a0 = a0list[cri.loc], a0.list = a0list, beta.m.mat = beta.mat, beta.i.mat = beta.list.inter, 
        beta.m = beta.mat[ind.list.main[[cri.loc]], cri.loc], beta.i = beta.i, df = df.list, 
        df.m = df.m.list, df.i = df.i.list, lambda = lambda[1:k], mainInd.list = ind.list.main, 
        mainInd = ind.list.main[[cri.loc]], cri.list = cri.list, loglik.list = loglik.list, 
        cri.loc = cri.loc, all.locs = all.locs, interInd.list = ind.list.inter.xlab, 
        interInd = interInd, family = family, X = oldX, y = y)
    class(val) = "RAMP"
    return(val)
} 
