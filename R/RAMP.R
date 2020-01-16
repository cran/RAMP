#' Regularization Algorithm under Marginality Principle (RAMP) for high
#' dimensional generalized quadratic regression.
#'
#' Regularization Algorithm under Marginality Principle (RAMP) for high
#' dimensional generalized quadratic regression.
#' @useDynLib RAMP, .registration = TRUE
#' @importFrom stats coef
#' @importFrom stats glm
#' @importFrom stats var
#' @export
#' @param X input matrix, of dimension nobs x nvars; each row is an observation
#' vector.
#' @param y  response variable, of dimension nobs x 1. continous for
#' \code{family='gaussian'}, binary for \code{family='binomial'}.
#' @param family response type. Default is 'gaussian'. The other choice is 'binomial'
#' for logistic regression.
#' @param penalty Choose from \code{LASSO}, \code{SCAD} and \code{MCP}. Default
#' is 'LASSO'.
#' @param gamma concavity parameter. If NULL, the code will use 3.7 for
#' 'SCAD' and 2.7 for 'MCP'.
#' @param inter whether to select interaction effects. Default is TRUE.
#' @param hier whether to enforce strong or weak heredity. Default is 'Strong'.
#' @param eps the precision used to test the convergence. Default is 1e-15.
#' @param tune tuning parameter selection method.
#' 'AIC', 'BIC', 'EBIC' and 'GIC' are available options. Default is EBIC.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param inter.penalty.factor the penalty factor for interactions effects. Default is 1. 
#' larger value discourage interaction effects. 
#' @param lam.list a user supplied \eqn{\lambda} sequence.
#' typical usage is to have the program compute its own
#' \code{lambda} sequence based on \code{lambda.min.ratio} and \code{n.lambda}.
#'  supplying a value of \eqn{\lambda} overrides this.
#' @param lambda.min.ratio optional input. smallest value for \code{lambda}, as
#' a fraction of \code{max.lam}, the (data derived) entry value. the default
#' depends on the sample size \code{n} relative to the number of variables
#' \code{p}. if \code{n} > \code{p}, the default is 0.0001. otherwise, the
#' default is 0.01.
#' @param max.iter maximum number of iteration in the computation. Default is
#' 100.
#' @param max.num optional input. maximum number of nonzero coefficients.
#' @param n.lambda the number of \code{lambda} values. Default is 100.
#' @param ebic.gamma  the gamma parameter value in the EBIC criteria. Default is
#'   1.
#' @param refit   whether to perform a MLE refit on the selected model. Default
#' is TRUE.
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
RAMP <- function(X, y, family = "gaussian", penalty = "LASSO", gamma = NULL, inter = TRUE, 
    hier = "Strong", eps = 1e-15, tune = "EBIC", penalty.factor = rep(1,ncol(X)), inter.penalty.factor = 1, lam.list, lambda.min.ratio, max.iter = 100, 
    max.num, n.lambda = 100, ebic.gamma = 1, refit = TRUE, trace = FALSE) {
    ## hier = 'Strong', or 'Weak', strong or weak heredity.
    if (penalty == "SCAD" & is.null(gamma)) {
        gamma = 3.7
    }
    if (penalty == "MCP" & is.null(gamma)) {
        gamma = 2.7
    }
    if(is.null(gamma)){
        gamma = 0
    }
    pentype = switch(penalty, LASSO = 1, MCP = 2, SCAD = 3)
    
    n = dim(X)[1]
    p = dim(X)[2]
    if (missing(max.num)) 
        max.num = p + 1 else max.num = min(max.num, p + 1)
    
    if (missing(lambda.min.ratio)) 
        lambda.min.ratio = ifelse(n < p, 0.01, 1e-04)
    
    ################ prepare variables
    lambda = NULL
    beta.mat = NULL  ##storing the beta coefficents along the path
    beta = rep(0, p)  ##the current beta estimates excluding the interation terms
 
    index = NULL
    
    ################ save the original X for generating interactions
    X0 = X
    ################ standardize design matrix
    X = scale(X0) 
    xm = attr(X, "scaled:center")
    xsd = attr(X, "scaled:scale")
    
    
    ############# lambda list
    if (family == "gaussian") {
    max.lam = max(abs((1/n) * (t(X) %*% (y - mean(y)))))
    a0list = rep(mean(y), n.lambda)
    } else if (family == "binomial") {
        max.lam = max(abs((1/n) * (t(X) %*% (y - mean(y)))))
        a0list = rep(log(mean(y)/(1 - mean(y))), n.lambda)
    } else if (family == "poisson") {
        max.lam = max(abs((1/n) * (t(X) %*% (y - mean(y)))))
        a0list = rep(log(mean(y)), n.lambda)
    }  
    
    if (missing(lam.list)) {
        min.lam = max.lam * lambda.min.ratio
        lam.list = exp(seq(from = log(max.lam), to = log(min.lam), length.out = n.lambda))
    } else {
        lam.list = lam.list[lam.list <= max.lam]
        n.lambda = length(lam.list)
    }
    
    
    ########initialie
    
    a0 = a0list[1]
    loglik.list = cri.list = AIC.list = BIC.list = EBIC.list = GIC.list = 
      df.list = df.m.list = df.i.list = vector('numeric', n.lambda)
    ind.list.inter = ind.list.main = ind.list.inter.xlab = beta.list.inter = 
      vector('list', n.lambda)

    
    ############## main part
    
    
    ## expand X matrix dynamically each step by incorporating the new candidates of
    ## interactions, which are the interactions of the current active variables minus
    ## the current interaction effects, the interactions will play the same role as
    ## the main effect, except we keep track of the list through the lambda sequence.
    ## Keep track of other indices as well.colnames(X)
    
    colnames(X) = paste("X", 1:p, sep = "")
    
    nonPen = rep(0, p)
    for (k in 1:n.lambda) {
        nonPen = rep(0, p) ### heredity force
        break.ind = FALSE
        lam = lam.list[k]
         index = which(abs(beta[1:p]) > eps)  ###selected main effects from last step
        if (length(index) > 0) {
            ### find the candidate interation effects following strong heredity strong heredity
            aa = outer(index, index, f <- function(x, y) {
                paste("X", pmin(x, y), "X", pmax(x, y), sep = "")
            })
            aa = as.vector(aa)
            newinter = outercprod(X0[, index], X0[, index])
            bb = unique(as.vector(aa))
            loc = match(bb, aa)
            newinter = newinter[, loc]
            if (hier == "Weak") {
                ### append strong heredity case to form the weak ones.
                rindex = setdiff(1:p, index)
                if (length(rindex) > 0) {
                  ## if there are candidate interactions under weak heredity
                  naa = outer(index, rindex, f <- function(x, y) {
                    paste("X", pmin(x, y), "X", pmax(x, y), sep = "")
                  })
                  naa = as.vector(naa)
                  nnewinter = outercprod(X0[, rindex], X0[, index])
                  # colnames(nnewinter) = naa
                  newinter = cbind(newinter, nnewinter)
                  bb = c(bb, naa)
                }
                
            }
            
            
            curInter = colnames(X)[-(1:p)]
            candInter = setdiff(bb, curInter)
            curloc = match(candInter, bb)
            if ( length(index) + length(curloc) > 5e5 ) 
            {
              k = k - 1
              break
            }
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
                xnewname = c(colnames(X), candInter)
                tmp = scale(newinter)
                X  = cbind(X,tmp)
                colnames(X) = xnewname
                xim1 = attr(tmp, "scaled:center")
                xisd1 = attr(tmp, "scaled:scale")
                xm = c(xm, xim1)
                xsd = c(xsd, xisd1)
                
                beta = c(beta, rep(0, ncandInter))  #expand the beta coefficent vector
                # to include the candiate interaction terms.
            }
        }
      
        pf = c(penalty.factor,rep(inter.penalty.factor,ncol(X) - p))
        nonpenind = which(nonPen != 0)
        pf[nonpenind] = 0
          for (ite in 1:max.iter) {
            if (ite == 1) {
                cd.temp1 = cd.general(X = X, y = y, a0 = a0, beta = beta, epsilon = eps, 
                  max.iter = 1, lambda = lam, family = family, bInd = break.ind, 
                  pf = pf, pentype = pentype, gamma = gamma)
                a0 = cd.temp1$a0
                beta = cd.temp1$beta
                ind1 = which(abs(beta) > eps)
            }
            
            
            # #CD
            if (ite > 1) {
                  cd.temp2 = cd.general(X = X[, ind1], y = y, a0 = a0, beta = beta[ind1], 
                    epsilon = eps, max.iter = max.iter, lambda = lam, family = family, 
                    bInd = break.ind, pf = pf[ind1], pentype = pentype, 
                    gamma = gamma)
                  a0 = cd.temp2$a0
                  beta2 = beta
                  beta2[ind1] = cd.temp2$beta
                  
                  
                  # ########redetect active set check.beta=new.beta
                  cd.temp3 = cd.general(X = X, y = y, a0 = a0, beta = beta2, epsilon = eps, 
                    max.iter = 1, lambda = lam, family = family, bInd = break.ind, 
                    pf = pf, pentype = pentype, gamma = gamma)
                  a0 = cd.temp3$a0
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
        ind.list.inter[[k]] = which(abs(beta[-(1:p)]) > eps)  #record the interaction
        # index pair location list
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
        
        beta.n = beta
        a0.n = a0
        if (refit == TRUE & length(index) > 0 ) {
          lmfit = glm(y ~ X[,index], family = family)
          beta.lmfit = coef(lmfit)
          a0.n = beta.lmfit[1]
          beta.n[index] = beta.lmfit[-1]
        }
    
        df.list[k] = size.main + size.inter
        
        loglik.list[k] = loglik(cbind(1,X), y, c(a0.n,beta.n), family)
        
         
        if (length(beta) > p) {
            tmp = which(abs(beta[-(1:p)]) > eps)
            beta = beta[c(1:p, p + tmp)] 
            beta.n = beta.n[c(1:p, p + tmp)]
            X = X[, c(1:p, p + tmp)]
            xm = xm[c(1:p, p + tmp)]
            xsd = xsd[c(1:p, p + tmp)]
            # ind.list.inter.xlab[[k]] = NULL
            if (length(tmp) > 0) {
                ## if interaction effects are selected.
                ind.list.inter.xlab[[k]] = colnames(X)[-(1:p)]
                beta.list.inter[[k]] = beta.n[-(1:p)]/xsd[-(1:p)]
            } 
        }
        if (family == 'binomial' & max(abs(beta)) > 10) {
          k = k - 1
          break
        }
        
        
        beta.s = beta.n/xsd
        a0.s = a0.n - sum(beta.s*xm)
        
        a0list[k] = a0.s
        beta.mat = cbind(beta.mat, beta.s[1:p])
    
        
        
        df.m.list[k] = length(ind.list.main[[k]])
        df.i.list[k] = length(beta) - p
        
        if (df.list[k] >= n - 1) 
            break
        if (break.ind == TRUE) {
            print("Warning: Algorithm failed to converge for all values of lambda")
            break
        }
        
        
        if (trace == T) {
            cat("k=", k, "current main effects", ind.list.main[[k]], "\n", sep = " ")
            cat("k=", k, "current interaction", ind.list.inter.xlab[[k]], "\n", sep = " ")
        }
        # end the outer loop: decrease in lambda
    }
    
    #df.m.max = max(df.m.list)
    # if (inter == FALSE) {
    #   p.eff = p
    # } else if (hier == "Strong") {
    #   p.eff = p + df.m.max *(df.m.max + 1)/2
    # } else if (hier == "Weak") {
    #   p.eff = p + df.m.max *(df.m.max + 1)/2 + df.m.max*(p - df.m.max) 
    # }
    if (inter == FALSE) {
      p.eff = p
    } else if (hier == "Strong") {
      p.eff = p + df.m.list *(df.m.list + 1)/2
    } else if (hier == "Weak") {
      p.eff = p + df.m.list *(df.m.list + 1)/2 + df.m.list*(p - df.m.list)
    }
    
    
      AIC.list = loglik.list + 2 * df.list
      BIC.list = loglik.list + log(n) * df.list
      EBIC.list = loglik.list + log(n) * df.list + 2 * ebic.gamma * 
        log(choose(p.eff, df.list))
      GIC.list = loglik.list + log(log(n)) * log(p.eff) * df.list
 
      cri.list = switch(tune, AIC = AIC.list, BIC = BIC.list, EBIC = EBIC.list,
                        GIC = GIC.list)
    region  = which(df.list[1:k] < sqrt(n))
  
   #browser()
    cri.loc = which.min(cri.list[region])
    AIC.loc = which.min(AIC.list[region])
    BIC.loc = which.min(BIC.list[region])
    EBIC.loc = which.min(EBIC.list[region])
    GIC.loc = which.min(GIC.list[region])
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
        beta.m = beta.mat[ind.list.main[[cri.loc]], cri.loc], beta.i = beta.i, df = df.list[1:k], 
        df.m = df.m.list[1:k], df.i = df.i.list[1:k], lambda = lambda[1:k], mainInd.list = ind.list.main[1:k], 
        mainInd = ind.list.main[[cri.loc]], cri.list = cri.list[1:k], loglik.list = loglik.list[1:k], 
        cri.loc = cri.loc, all.locs = all.locs, interInd.list = ind.list.inter.xlab[1:k], 
        interInd = interInd, family = family, X = X0, y = y)
    class(val) = "RAMP"
    return(val)
}
