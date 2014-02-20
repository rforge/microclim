\name{centres}
\alias{centres}
\docType{data}
\title{
Centring values for function micEuroclim
}
\description{
Typical values for the predictor variables, to be used as defaults when no data are specified.
}
\usage{data(centres)}
\format{
  A data frame with 10 observations on the following 4 variables.
  \describe{
    \item{\code{mean_open}}{a numeric vector}
    \item{\code{mean_forest}}{a numeric vector}
    \item{\code{daily_open}}{a numeric vector}
    \item{\code{daily_forest}}{a numeric vector}
  }
}
\details{
These values were the means of the variables in the original calibration data set (see \strong{Source}).
}
\source{
Gunton et al. (submitted) "Predicting ground temperatures across European landscapes"
}
\examples{
data(centres)
str(centres)
}
\keyword{datasets}
