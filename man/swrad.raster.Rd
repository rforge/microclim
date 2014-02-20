\name{swrad.raster}
\alias{swrad.raster}
\alias{swrad2.raster}
\alias{srad.raster}
\alias{trans}

\title{
Shortwave radiation predictions
}
\description{
Functions to predict shortwave solar radiation
}
\usage{
\method{swrad}{raster}(dtm, sw, lat, day, elevation, soltime, localtime, lon, merid = 0, dst = 0, trans0 = 0.2, ...)
\method{swrad2}{raster}(dtm, sw, lat, day, elevation, soltime, soltime.sw, localtime, localtime.sw = localtime, lon, merid = 0, dst = 0, trans0 = 0.2, ...)
\method{srad}{raster}(dtm, trans, lat, day, elevation, soltime, localtime, lon, merid = 0, dst = 0, trans0 = 0.2, ...)
trans(sw, elevation = 0, lat, day, soltime, localtime, lon = 0, merid = 0, dst = 0, trans0 = 0.2, ...)
}

\arguments{
  \item{dtm}{
Digital terrain model, as a raster object
}
  \item{sw}{
Shortwave radiation (in Watts per square metre) incident on horizontal ground. May be obtained from regional weather station or reanalysis data.
}
  \item{trans}{
Transmission of short-wave radiation, as a proportion of the difference between the transmission of clear and overcast skies. Typically estimated by proportional cloud cover.
}
  \item{lat}{
Latitude, in degrees. If missing, it will be taken from dtm, assuming its y extent is calibrated in degrees
}
  \item{day}{
Day of the year (1 - 365)
}
  \item{elevation}{
Average elevation of the site (in m), passed to the function \code{beamrad}. If missing, the mean of \code{dtm} will be calculated. For \code{elevation = NA}, \code{swrad.raster}, \code{swrad2.raster} and \code{srad.raster} will use the raster \code{dtm}, which may result in slower performance.
}
  \item{soltime}{
Local solar time at which predictions are required, in hours from midnight (only used if \code{localtime} is missing)
}
  \item{soltime.sw}{
Local solar time at the time to which \code{sw} pertains
}
  \item{localtime}{
Local clock time at which predictions are required, in hours from midnight 
}
  \item{localtime.sw}{
Local clock time at the time to which \code{sw} pertains (used if \code{soltime.sw} is missing)
}
  \item{lon}{
For interpreting \code{localtime} and \code{localtime.sw}: Longitude, in degrees east from Greenwich
}
  \item{merid}{
For interpreting \code{localtime} and \code{localtime.sw}: Longitude of local standard time meridian, in degrees east from Greenwich. Thus the default (0) corresponds to the Greenwich mean time meridian; for Central European time (most of Western Europe except UK, Eire and Portugal), use merid = 15.
}
  \item{dst}{
For interpreting \code{localtime} and \code{localtime.sw}: Correction for summer time (=1 if local time has been adjusted for summer daylight-saving time, =0 if not)
}
  \item{trans0}{
typical transmission (proportion) of diffuse short-wave radiation through overcast sky. If in doubt, leave at 0.2.
}
  \item{\dots}{
Additional arguments to pass to subsidiary functions (currently \code{\link{beamrad}} and \code{\link{diffuserad}}).
}
}
\details{
\code{swrad.raster} estimates short-wave radiation on a landscape given an estimate for short-wave radiation incident on the horizontal plane at the same point in time. \code{swrad2.raster} allows the estimate for horizontal radiation to apply to a different time (\code{solartime.sw} or \code{localtime.sw}). 
\cr
These functions use \code{sw} to obtain a shortwave radiation value standardised with respect to the expected shortwave radiation at the specified time. This calculation makes reference to the Sun's expected height; if this is less than 10 degrees, a warning is given to indicate that results may be unreliable, since small errors in the solar angle may be magnified. In that case, it may be better to use function \code{srad} with an estimate for transmission based on cloud cover or other data.
\cr
\code{srad.raster} and \code{trans} together provide the same functionality as \code{swrad.raster}: \code{srad.raster} uses a given estimate, \code{trans}, for cloud cover (i.e. percentage transmission above that which pertains with overcast skies) instead of requiring \code{sw}, while \code{trans} gives such an estimate for transmission based on \code{sw} with reference to expectations derived from \code{\link{beamrad}} and \code{\link{diffuserad}}.
}
\value{
\code{swrad.raster} and \code{swrad2.raster} each return a list: 
\item{swrad}{raster of short-wave radiation estimates in Watts per square metre} 
\item{trans}{single-value estimate of transmission through the sky (see \code{trans} below).}

\code{srad.raster} returns a raster of shortwave radiation estimates in Watts per square metre.
\cr
\code{trans} returns estimated proportional transmission of shortwave radiation through the sky. If any of its arguments comprise more than a single value (possibly rasters), they should be equal to each other in dimensions and extent and the output will be of corresponding size. So long as \code{sw} is less than the predicted potential shortwave radiation, values will fall in the range 0 - 1, with 0 corresponding to overcast sky. Otherwise, values > 1 are possible.
}
\references{

}
\author{
Jon Bennie, Richard Gunton
}

\seealso{
\code{\link{beamrad}},
\code{\link{swrad}}
}

\examples{
trans(sw=600, lat=54, day=172, soltime=12)

if(require(raster)){
data(volcano)
dtm <- raster(volcano,xmn=0,xmx=610,ymn=0,ymx=870,crs="+proj=equirectangular")
dtm.srad <- srad.raster(dtm, lat=-36.9, trans=0.7, day=172, soltime=12)
dtm.swrad <- swrad.raster(dtm, sw=400, lat=-36.9, day=172, soltime=12)
dtm.swrad2 <- swrad2.raster(dtm, sw=400, lat=-36.9, day=172, soltime=15, soltime.sw=12)

par(mfrow=c(2,2))
plot(dtm,col=terrain.colors(255))
plot(dtm.srad,col=heat.colors(255))
plot(dtm.swrad$swrad,col=heat.colors(255))
plot(dtm.swrad2$swrad,col=heat.colors(255))
}
}