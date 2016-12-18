#' Model prediction based on a fitted RAMP object.
#' 
#' Similar to the usual predict methods, 
#' this function returns predictions from a fitted \code{'RAMP'} object.
#' @export
#' @param   object Fitted \code{'RAMP'} model object.
#' @param newdata Matrix of new values for \code{x} at which predictions 
#' are to be made, without the intercept term.
#' @param type  Type of prediction required. 
#' Type \code{'response'} gives the fitted values for \code{'gaussian'}, 
#' fitted probabilities for \code{'binomial'}, 
#' fitted mean for \code{'poisson'}, 
#' and the fitted relative risk for \code{'cox'}. 
#' Type \code{'link'} returns the linear predictors for 
#' \code{'binomial'}, \code{'poisson'} and \code{'cox'} models;
#'  for \code{'gaussian'} models it is equivalent to type 
#'  \code{'response'}. Type \code{'class'} applies only to \code{'binomial'} models,
#'   and produces the class label corresponding to the maximum probability (0-1 labels).
#' @param  allpath  allpath = T will output all the predictions on the solution path. allpath = FALSE will only output the one the criterion selected in the \code{'RAMP'} object.
#' @param \dots Not used. Other arguments to predict.
#' @return The object returned depends on type.
#' @seealso \code{\link{RAMP}},\code{\link{print.RAMP}}
predict.RAMP <- function(object, newdata = NULL, type = c("link", "response", "class"), 
    allpath = FALSE, ...) {
    # object: a fitted object of class inheriting from 'RAMP'.  X: optionally, a data
    # matrix with which to predict. If omitted, the data in object will be used.
    # type: the type of prediction required. The default is on the scale of the
    # linear predictors; the #alternative 'response' is on the scale of the response
    # variable. Thus for a default binomial model the default predictions are of
    # log-odds (probabilities on logit scale) and type = 'response' gives the
    # predicted probabilities.  type='class' gives the class lable for binomail
    # distribution.
    X = newdata
    oldX = object$X
    if (is.null(X)) 
        X = oldX
    n = nrow(X)
    p = ncol(X)
    mainind.list = object$mainInd.list
    interind.list = object$interInd.list
    a0.list = object$a0.list
    k = length(mainind.list)  ###total number of models visited
    
    eta = matrix(0, n, k)
    model.list = as.list(NULL)
    count = 0
    for (i in 1:k) {
        # cat('i= ', i, ' out of ', k,'\n')
        mainind = mainind.list[[i]]
        Xmain = X[, mainind]
        count = count + 1
        Xinter = NULL
        Xi = Xmain
        coef = object$beta.m.mat[mainind, i]
        if (i <= length(interind.list)) {
            interind = interind.list[[i]]
            if (length(interind) > 0 ) {
                ## candidate interaction terms
                for (indInter in 1:length(interind)) {
                  pair = as.numeric(strsplit(interind[indInter], "X")[[1]][2:3])
                  Xinter = cbind(Xinter, X[, pair[1]] * X[, pair[2]])
                }
                
                
                if (!is.na(pair[1]) && !is.na(pair[2])) {
                  Xi = cbind(Xmain, Xinter)
                  coef = c(object$beta.m.mat[mainind, i], object$beta.i.mat[[i]])
                }
            }
        }
        eta[, i] = a0.list[i]
        if (length(mainind) > 0) {
            eta[, i] = a0.list[i] + as.matrix(Xi) %*% coef
        }
    }
    if (allpath == FALSE) {
        eta = eta[, object$cri.loc]
    }
    if (match.arg(type) == "response") 
        return(exp(eta)/(1 + exp(eta)))
    if (match.arg(type) == "link") 
        return(eta)
    if (match.arg(type) == "class") 
        return(eta > 0)
    
}
