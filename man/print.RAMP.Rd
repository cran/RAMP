\name{print.RAMP}
\alias{print.RAMP}
\title{Result summary of a fitted RAMP object.
}
\description{
Similar to the usual print methods, this function summarize results from a fitted \code{"RAMP"} object.
}
\usage{
\method{print}{RAMP}(x,digits = max(3, getOption("digits") - 3),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 Fitted \code{"RAMP"} model object.
}
  \item{digits}{
The number of significant digits for the coefficient estimates.
}
\item{\dots}{Not used. Other arguments to predict. }
}

\value{
No value is returned.
}

\author{
Yang Feng, Ning Hao and Hao Helen Zhang.
}


\seealso{
\code{\link{RAMP}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ function }% __ONLY ONE__ keyword per line
