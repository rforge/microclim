\name{micEuroclim}
\alias{micEuroclim}
\title{
Predict soil surface temperature anywhere in Europe
}
\description{
Predicts soil surface temperature (monthly mean or daily at 13.40, in open grassland or under trees) for any location in Europe on a given date and/or under specified weather conditions.
}
\usage{
micEuroclim(dtm, coast, lat, cloud, day = 197, month, air, conv, wind, air.mean, sw, soltime.sw=12, localtime.sw, lon, merid=0, dst=0, pred = c("mean_open", "mean_trees", "daily_open", "daily_trees"), daily = 12, trans0 = 0.2, scale.threshold = 0.5, cleanup = FALSE, reminders = TRUE, ...)
}
\arguments{
  \item{dtm}{
Digital terrain model, as a \code{\link{raster}} object
}
  \item{coast}{
Distances from the nearest sea (in km), in case \code{dtm} contains no sea (blocks of \code{NA} values; see \code{\link{coastal}}). Either a raster of coastal distances of the same extent as \code{dtm}, or a single approximate value to be used across the whole map (e.g. if the nearest sea is > 100 km away). See \code{\link{coastal}} for further details.
}
  \item{lat}{
Latitude (average) of the site. Only used if \code{dtm} does not have a "LonLat" projection.
}
  \item{cloud}{
Proportional cloud cover (see Details). Either a raster of the same extent as \code{dtm}, or a single value to be used across the whole map.
}
  \item{day}{
Day of the year (1 - 365) at which a prediction is required
}
  \item{month}{
Month (1 - 12) at which a prediction is required
}
  \item{air}{
Temperature of the air at 2 - 10m above ground (in degrees Celsius). See Details.
}
  \item{conv}{
Convection variable, calculated as {air temperature anomaly} * sqrt(windspeed), where the first term is air temperature minus annual mean air temperature, and the wind speed is in in m/s. See Details.
}
  \item{wind}{
Wind speed (m/s). See Details.
}
  \item{air.mean}{
Mean annual air temperature (in Celsius degrees)
}
  \item{sw}{
Shortwave radiation (in Watts per square metre) incident on horizontal ground at \code{localtime} (to be passed to \code{\link{trans}} or \code{\link{swrad2.raster}}). May be obtained from regional weather station or reanalysis data.
}
  \item{soltime.sw}{
For \code{sw}: Local solar time, in hours from midnight, at which any value for \code{sw} is given.
}
  \item{localtime.sw}{
For \code{sw}: Local clock time, in hours from midnight, at which any value for \code{sw} is given. Overrides any value for \code{soltime.sw}.
}
  \item{lon}{
For interpreting \code{localtime.sw}: Longitude, in degrees east from Greenwich
}
  \item{merid}{
For interpreting \code{localtime.sw}: Longitude of local standard time meridian, in degrees east from Greenwich. Thus the default (0) corresponds to the Greenwich Mean Time meridian. For Central European Time (most of Western Europe except UK, Eire and Portugal), use merid = 15.
}
  \item{dst}{
For interpreting \code{localtime.sw}: Correction for summer time (=1 if local time has been adjusted for summer daylight-saving time, =0 if not)
}
  \item{pred}{
String or index specifying what microclimate predictions are required: one of (1) monthly means at open sites, (2) monthly means at tree-covered sites, (3) predictions for 13.40 at open sites, or (4) predictions for 13.40 at tree-covered sites.
}
  \item{daily}{
Number of daily intervals (temporal resolution) at which calculations should be made when mean solar radiation is calculated.
}
  \item{trans0}{
Typical transmission (proportion) of diffuse short-wave radiation through overcast sky. If in doubt, leave at 0.2.
}
  \item{scale.threshold}{
Extent (in degrees latitude) of \code{dtm} below which a single, average, value for its latitude may be used, in order to speed computations and save memory.
}
  \item{cleanup}{
Logical: should rasters for estimated short-wave radiation be deleted (to save memory), and thus not returned in the output object?
}
  \item{reminders}{
Logical: should reminders about assumptions concerning input data be printed on the screen when the function runs?
}
  \item{\dots}{
Additional arguments to pass to component functions (currently \code{coastal}, \code{trans}, \code{srad.raster}, \code{swrad2.raster}).
}
}
\details{
The minimum requirement is a raster of elevations (\code{dtm}), with specified projection, and a day; in this case the output will give relative temperature deviations for the specified time of year: these can be either monthly means or instantaneous estimates at 13.40, for either open grassland or tree-covered sites, as specified by \code{pred}.  If a value for \code{air} is additionally be provided, the output comprises actual predictions; these will be improved if an estimate for \code{conv}, or \code{wind} together with \code{air.mean}, is provided.\cr
\cr
If any values for \code{air} or \code{conv} are provided, they must have appropriate temporal support for the type of predictions required. For \code{pred=1} or \code{pred=2}, they should simply be monthly means. For \code{pred=3}, \code{air} must be an estimate for the time when the prediction is required (i.e. 13:40) and \code{conv} must be a morning mean (averaged over the previous midnight to noon). For \code{pred=4}, \code{air} must be a morning mean and \code{conv} must be estimated for the time when the prediction is required (i.e. 13:40).\cr
\cr
If \code{conv} is not specified and values for \code{wind} are instead provided, there are further requirements for \code{wind} and \code{air}: \code{wind} should be a vector of typical daily wind speeds for the month or date in question, evenly spaced throughout the day, and \code{air} should be an equal-length vector of corresponding air temperatures.\cr
\cr
A typical or particular value for \code{sw} (short-wave radiation on level ground) can be provided, along with the time (defined by \code{localtime} with \code{merid} and \code{dst}; otherwise assumed to be solar noon) to which it applies. This will be used to calculate standard levels of incident radiation and also to estimate cloud cover, if \code{cloud} is not specified.
}
\value{
\item{microclim}{A raster of predicted temperatures (degrees Celsius)}
\item{rad}{Predicted short-wave radiation values (Watts per sq. m): a raster, unless cleanup=TRUE}
\item{cloud}{Predicted cloud cover (proportion)}
\item{trans}{Value used for sky transmission (proportion); this will equal \code{1- cloud*(1- trans0)} if \code{sw} is not provided. Otherwise, it is obtained using \code{\link{trans}}.}
\item{coast}{Distance from the sea (km): a raster if this was provided or if dtm contained NA values representing sea}
\item{air}{Air temperature used for predictions (degrees Celsius)}
\item{conv}{Convection value used for predictions}
\item{day}{Day of the year to which predictions apply}
\item{call}{The function call}
}
\references{
Gunton et al. (submitted) \strong{Predicting ground temperatures across European landscapes}
}
\author{
Richard Gunton
}

\seealso{
 \code{\link{coastal}}, 
 \code{\link{swrad2.raster}}, 
 \code{\link{srad.raster}}, 
 \code{\link{trans}}
}
\examples{
if(require(raster)) {
data(volcano)
dtm <- raster(volcano,xmn=0,xmx=610,ymn=0,ymx=870,crs="+proj=equirectangular")
micro3 <- micEuroclim(dtm, lat=36.9, coast=4.5, pred=3, day=172)
# imagining this site were in southern Spain!

par(mfrow=c(1,3))
plot(dtm,col=terrain.colors(255))
plot(micro3$microclim)
}
}