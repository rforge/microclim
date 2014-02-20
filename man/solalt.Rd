\name{solalt}
\alias{solalt}
\alias{solazi}
\title{
Calculate the altitude and azimuth of the sun for a given time and latitude
}
\description{
Calculates the angular altitude of the sun above the horizon for a given (solar) time, day of year and latitude
}
\usage{
solalt(time, day, lat)
solazi(time, day, lat)
}
\arguments{
  \item{time}{
solar time, in hours from midnight
}
  \item{day}{
day of the year (0 to 365)
}
  \item{lat}{
latitude, in degrees
}
}

\value{
The sun's angular altitude above the horizon, or azimuth from north, in degrees.
}
\references{
}
\author{
Jon Bennie
}

\seealso{
 \code{\link{solartime}}
}
\examples{
times <- c(1:24)
alt <- solalt(time=times, day=175, lat=52)
azi <- solazi(time=times, day=175, lat=52)
par(mfrow=c(2,1))
plot(alt ~ times, type="l", xlab="Time", main="Solar altitude")
plot(azi ~ times, type="l", xlab="Time", main="Solar azimuth")
}
