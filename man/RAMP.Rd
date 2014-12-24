\name{RAMP}
\alias{RAMP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Regularization Algorithm under Marginality Principle (RAMP) for high dimensional generalized quadratic regression.
%%  ~~function to do ... ~~
}
\description{
	RAMP
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RAMP(X, y, family = "gaussian", penalty = "LASSO", gamma, inter = TRUE,
 eps = 1e-15, tune = "EBIC", lam.list, lambda.min.ratio, max.iter = 100, 
 max.num, n.lambda = 100, ebic.gamma = 1, refit = TRUE, trace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
input matrix, of dimension nobs x nvars; each row is an observation vector.
}
  \item{y}{
response variable, of dimension nobs x 1. non-negative counts for 
\code{family="poisson"}, binary for \code{family="binomial"}.
}
  \item{family}{
response type. Default is "gaussian".
}
  \item{penalty}{
Currently, only \code{LASSO} is implemented.
}
  \item{gamma}{
concavity parameter. Experimental feature.
}
  \item{inter}{
whether to select interaction effects. Default is TRUE.

}
  \item{eps}{
the precision used to test the convergence. Default is 1e-15.

}
  \item{tune}{
tuning parameter selection method. "AIC", "BIC" and "EBIC" are available options. Default is EBIC.

}
  \item{lam.list}{
a user supplied \eqn{\lambda} sequence. typical usage is to have the program compute its own \code{lambda} sequence based on \code{lambda.min.ratio} and \code{n.lambda}. supplying a value of \eqn{\lambda} overrides this.
}
  \item{lambda.min.ratio}{
optional input. smallest value for \code{lambda}, as a fraction of \code{max.lam}, the (data derived) entry value. the default depends on the sample size \code{n} relative to the number of variables \code{p}. if \code{n} > \code{p}, the default is 0.0001. otherwise, the default is 0.01.
}
  \item{max.iter}{
maximum number of iteration in the computation. Default is 100.

}
  \item{max.num}{
optional input. maximum number of nonzero coefficients.
}
  \item{n.lambda}{
the number of \code{lambda} values. Default is 100.

}
  \item{ebic.gamma}{
the gamma parameter value in the EBIC criteria. Default is  1.

}
  \item{refit}{
whether to perform a MLE refit on the selected model. Default is TRUE.

}
	\item{trace}{
whether to trace the fitting process. Default is FALSE.
	}
}
\value{
\item{a0}{intercept vector of length(\code{lambda}).}
\item{beta}{nvar x length(\code{lambda}) matrix of coefficients.}
}
\author{
Yang Feng, Ning Hao and Hao Helen Zhang.
}

\seealso{
 \code{\link{predict.RAMP}}
}
\examples{
set.seed(0)
n = 500
p = 100 #Can be changed to a much larger number say 100000
x = matrix(rnorm(n*p),n,p)
eta = 1 * x[,1] + 2 * x[,3]  + 3*x[,6]  + 4*x[,1]*x[,3] + 5*x[,1]*x[,6]
y =  eta + rnorm(n)

fit1 = RAMP(x, y)
fit1    ###examine the results


y = rbinom(n, 1, 1/(1+exp(-eta)))
fit2 = RAMP(x, y, family='binomial')  ###for binary response

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{function}