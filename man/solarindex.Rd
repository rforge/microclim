\name{solarindex}
\alias{solarindex}
\title{
Calculate proportion of solar beam incident on a slope
}
\description{
Calculates the proportion of the intensity of a direct solar beam that will strike a surface of given slope(s) and aspect(s).
}
\usage{
solarindex(slope, aspect, azi, alt)
}
\arguments{
  \item{slope}{
Slope, in degrees from horizontal
}
  \item{aspect}{
Aspect, in degrees clockwise from north
}
  \item{azi}{
Solar azimuth, in degrees clockwise from north
}
  \item{alt}{
Solar altitude, in degrees from the horizon
}
}
\details{
}
\value{
Proportion of direct solar radiation incident on a surface. If any of the arguments comprise more than a single value (possibly rasters), they should be equal to each other in dimensions and extent and the output will be of corresponding size.
}
\references{
}
\author{
Jon Bennie
}

\seealso{
\code{\link{solalt}}, 
\code{\link{solazi}}
}
\examples{
slopes <- c(0:45)
incidence <- solarindex(slope=slopes, aspect=45, azi=180, alt=30) 
plot(incidence ~ slopes, type="l", xlab="Slope (degrees)")
}