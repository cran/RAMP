\name{predict.RAMP}
\alias{predict.RAMP}
\title{Model prediction based on a fitted RAMP object.
}
\description{
Similar to the usual predict methods, this function returns predictions from a fitted \code{"RAMP"} object.
}
\usage{
\method{predict}{RAMP}(object, newdata = NULL, type = c("link", "response", "class"), 
allpath = FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
 Fitted \code{"RAMP"} model object.
}
  \item{newdata}{
Matrix of new values for \code{x} at which predictions are to be made, without the intercept term.
}
  \item{type}{
Type of prediction required. Type \code{"response"} gives the fitted values for \code{"gaussian"}, fitted probabilities for
    \code{"binomial"}, fitted mean for 
    \code{"poisson"}, and the fitted relative risk for \code{"cox"}. Type \code{"link"} returns the linear predictors for 
    \code{"binomial"}, \code{"poisson"} and \code{"cox"} models; for \code{"gaussian"} models it is equivalent to type 
    \code{"response"}. Type \code{"class"} applies only to \code{"binomial"} models, and produces the class label corresponding 
    to the maximum probability (0-1 labels).
}
\item{allpath}{
	allpath = T will output all the predictions on the solution path. allpath = FALSE will only output the one the criterion selected in the \code{"RAMP"} object.
	}
\item{\dots}{Not used. Other arguments to predict. }
}

\value{
The object returned depends on type.
}

\author{
Yang Feng, Ning Hao and Hao Helen Zhang.
}


\seealso{
\code{\link{RAMP}}
}
\examples{
set.seed(0)
n = 500
p = 100 #Can be changed to a much larger number say 100000
x = matrix(rnorm(2*n*p),2*n,p)
eta = 1 * x[,1] + 2 * x[,3]  + 3*x[,6]  + 4*x[,1]*x[,3] + 5*x[,1]*x[,6]
y =  eta + rnorm(2*n)
xtr = x[1:n, ]
ytr = y[1:n]

xte = x[(n+1):(2*n), ]
yte = y[(n+1):(2*n)]

				
fit1 = RAMP(xtr, ytr)
ypred = predict(fit1,xte)


####binary prediction 
y = rbinom(2*n, 1, 1/(1+exp(-eta)))
xtr = x[1:n, ]
ytr = y[1:n]

xte = x[(n+1):(2*n), ]
yte = y[(n+1):(2*n)]
				
fit2 = RAMP(xtr, ytr,family='binomial')
ypred = predict(fit2,xte,type='class')
mean(ypred!=yte) ##classification error on the test set

ypred = predict(fit2,type='class')
mean(ypred!=ytr) ##classification error on the training set
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{prediction}
\keyword{function}% __ONLY ONE__ keyword per line
