\name{solartime}
\alias{solartime}
\title{
Obtain local solar time from clock time
}
\description{
Calculates solar time at a given location, given the clock time
}
\usage{
solartime(localtime, lat, lon, merid = 0, dst = 0, day)
}

\arguments{
  \item{localtime}{
Local clock time, in hours from midnight
}
  \item{lat}{
Latitude, in degrees
}
  \item{lon}{
Longitude, in degrees
}
  \item{merid}{
Longitude of local standard time meridian, in degrees east from Greenwich. The default (0) corresponds to the Greenwich mean time meridian; for Central European time (most of Western Europe except UK, Eire and Portugal), use merid = 15.
}
  \item{dst}{
Local adjustment for summertime. For countries which advance their clocks in spring, this is 1 during the summer and 0 in winter.
}
  \item{day}{
Day of the year (1 to 365)
}
}

\value{
Vector of local solar time, in hours from midnight.
}
\references{
}
\author{
Jon Bennie
}

\examples{
day <- c(1:365)
soltime <- solartime(localtime=12, lat=52, lon=1, merid=0, dst=0, day=day)
plot(soltime ~ day, type="l") 	# Plot the equation of time
}


