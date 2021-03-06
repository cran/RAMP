% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.RAMP.R
\name{predict.RAMP}
\alias{predict.RAMP}
\title{Model prediction based on a fitted RAMP object.}
\usage{
\method{predict}{RAMP}(object, newdata = NULL, type = c("link",
  "response", "class"), allpath = FALSE, ...)
}
\arguments{
\item{object}{Fitted \code{'RAMP'} model object.}

\item{newdata}{Matrix of new values for \code{x} at which predictions 
are to be made, without the intercept term.}

\item{type}{Type of prediction required. 
Type \code{'response'} gives the fitted values for \code{'gaussian'}, 
fitted probabilities for \code{'binomial'}, 
fitted mean for \code{'poisson'}, 
and the fitted relative risk for \code{'cox'}. 
Type \code{'link'} returns the linear predictors for 
\code{'binomial'}, \code{'poisson'} and \code{'cox'} models;
 for \code{'gaussian'} models it is equivalent to type 
 \code{'response'}. Type \code{'class'} applies only to \code{'binomial'} models,
  and produces the class label corresponding to the maximum probability (0-1 labels).}

\item{allpath}{allpath = T will output all the predictions on the solution path. allpath = FALSE will only output the one the criterion selected in the \code{'RAMP'} object.}

\item{\dots}{Not used. Other arguments to predict.}
}
\value{
The object returned depends on type.
}
\description{
Similar to the usual predict methods, 
this function returns predictions from a fitted \code{'RAMP'} object.
}
\seealso{
\code{\link{RAMP}},\code{\link{print.RAMP}}
}
