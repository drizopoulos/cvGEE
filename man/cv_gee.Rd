\name{cv_gee}
\alias{cv_geee}

\title{
Proper Scoring Rules for Categorical Data
}

\description{
Calculates the logarithmic, quadratic/Brier and spherical based on a fitted mixed model
for categorical data.
}

\usage{
cv_gee(object, rule = c("all", "quadratic", "logarithmic", "spherical"), 
  max_count = 500, conf_int = FALSE, level = 0.95, 
  M = 30, fold = 5, seed = 1L, 
  return_data = FALSE)
}

\arguments{
\item{object}{an object inheriting from class \code{"geeglm"} of the \bold{geepack}.}
\item{max_count}{numeric scalar or vector denoting the maximum count up to which to calculate 
    probabilities; this is relevant for count response data.}
\item{return_newdata}{logical; if \code{TRUE} the values of the scoring rules are 
    ruturned as extra columns of the \code{data} behind \code{object}.}
}

\value{
A data.frame with (extra) columns the values of the logarithmic, quadratic and spherical
scoring rules calculated based on the fitted model and the observed responses in 
\code{newdata} or \code{newdata2}.
}

\references{
Carvalho, A. (2016). An overview of applications of proper scoring rules. 
\emph{Decision Analysis} \bold{13}, 223--242. doi:10.1287/deca.2016.0337
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}


\examples{
\donttest{
NA
NA
NA
}
}