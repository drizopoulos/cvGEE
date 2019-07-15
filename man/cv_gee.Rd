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
\item{rule}{character string indicating the type of scoring rule to be used.}
\item{max_count}{numeric scalar or vector denoting the maximum count up to which to calculate 
    probabilities; this is relevant for count response data.}
\item{conf_int}{logical; if \code{TRUE} confidence interval for the scoring rule values are provided.}
\item{level}{numeric scalar indicating the level of the confindence interval; used only when \code{cont_int = TRUE}.}
\item{M}{numeric scalar denoting how many time the split of the data in \code{fold}s will be performed. The reported scoring rules values are the average over the \code{M} replicates.}
\item{fold}{numeric scalar indicating the number of folds used in the cross-validation procedure.}
\item{seed}{numeric scalre providing the seed used in the procedure.}
\item{return_data}{logical; if \code{TRUE} the values of the scoring rules are 
    ruturned as extra columns of the \code{data} behind \code{object}.}
}

\value{
A list or a data.frame with elements or (extra) columns the values of the logarithmic, quadratic and spherical
scoring rules calculated based on the GEE object.
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