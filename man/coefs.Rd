\name{coefs}
\alias{coefs}
\docType{data}
\title{
Coefficient values for function micEuroclim
}
\description{
These coefficient values were obtained through the statistical modelling process described in the \strong{Source} reference below.
}
\usage{data(coefs)}
\format{
  A data frame with 23 observations on the following 4 variables, representing 4 predictive models.
  \describe{
    \item{\code{mean_open}}{a numeric vector}
    \item{\code{mean_forest}}{a numeric vector}
    \item{\code{daily_open}}{a numeric vector}
    \item{\code{daily_forest}}{a numeric vector}
  }
}
\details{
The function micEuroclim currently makes use of all these coefficient values. If alternative estimates are obtained from some other source, an alternative data frame following this structure could be substituted. The row names are parsed to determine the linear predictor for the predictive model, so alternative interaction terms may simply be specified here using row names of the format "X1*X2".
}
\source{
Gunton et al. (submitted) "Predicting ground temperatures across European landscapes"
}
\examples{
data(coefs)
str(coefs)
}
\keyword{datasets}
