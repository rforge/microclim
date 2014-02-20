\name{swrad}
\alias{swrad}
\title{
Predict shortwave radiation on a surface
}
\description{
Standard prediction for shortwave radiation on a surface of given slope and aspect at a given time, given an observed value for shortwave radiation on a horizontal surface.
}
\usage{
swrad(sw, slope, aspect, lat, day, elevation = 0, soltime, localtime, lon, merid = 0, dst = 0, trans0 = 0.2, ...)
}
\arguments{
  \item{sw}{
Shortwave radiation (in Watts per square metre) incident on horizontal ground. May be obtained from regional weather station or reanalysis data.
}
  \item{slope}{
Slope, in degrees from horizontal
}
  \item{aspect}{
Aspect, in degrees clockwise from north
}
  \item{lat}{
Latitude, in degrees
}
  \item{lon}{
Longitude, in degrees east from Greenwich
}
  \item{soltime}{
Local solar time, in hours from midnight
}
  \item{localtime}{
Local clock time, in hours from midnight (used if \code{soltime} is missing)
}
  \item{day}{
Day of the year (1 - 365)
}
  \item{elevation}{
Elevation of the site, in m
}
  \item{merid}{
Longitude of local standard time meridian, in degrees east from Greenwich. Thus the default (0) corresponds to the Greenwich mean time meridian; for Central European time (most of Western Europe except UK, Eire and Portugal), use merid = 15.
}
  \item{dst}{
Correction for summer time (=1 if local time has been adjusted for summer daylight-saving time, =0 if not)
}
  \item{trans0}{
typical transmission (proportion) of short-wave radiation through overcast sky. If in doubt, leave at 0.2.
}
  \item{\dots}{
Additional arguments to pass to subsidiary functions
}
}
\details{
The time when \code{sw} is observed should be specified by either \code{localtime} or \code{soltime}; the output will pertain to the same time. The function does two things: (1) it standardises the shortwave radiation with respect to the expected shortwave radiation at the specified time, and it applies \code{\link{solarindex}} to obtain a prediction for the specified slope and aspect. The calculations involve converting the horizontally-incident radiation to direct beam radiation with reference to the Sun's height; if this is less than 20 degrees, a warning is given to indicate that results may be unreliable.  In that case, it may be better to use function \code{\link{srad.raster}} with an estimate for transmission based on cloud cover or other data.
}
\value{
Predicted short-wave radiation, in Watts per square metre. If any of the arguments comprise more than a single value (possibly rasters), they should be equal to each other in dimensions and extent and the output will be of corresponding size. For predictions over a landscape raster, see \code{\link{swrad.raster}}.
}
\author{
Jon Bennie
}

\seealso{
\code{\link{swrad.raster}}
}
\examples{
swrad(sw=300, slope=30, aspect=180, lat=54, day=172, soltime=14)
}