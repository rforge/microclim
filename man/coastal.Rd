\name{coastal}
\alias{coastal}
\title{
Create a raster of coastal distances
}
\description{
An input raster is converted to a raster of estimated distances from the sea.
}
\usage{
coastal(dtm, nosea, notsea = 2)
}
\arguments{
  \item{dtm}{
A raster with a specified projection, where any sea is represented by NA values.
}
  \item{nosea}{
Value to be used (km) if dtm contains no NA values. See details for the default.
}
  \item{notsea}{
Approximate width (km) that NA regions should exceed in order to be counted as sea. This determines the map resolution that will be used, so only gives approximate control (bands of width >= 2*notsea will definitely be counted as sea).  If sea appears at the edge of \code{dtm} and there are no other NA values (e.g. for rivers) in the map, then \code{notsea} could be reduced to ensure recognition of even a narrow band of sea.
}
}
\details{
If the raster supplied contains regions of NA values at least 2km wide, these will be interpreted as sea and the returned raster comprises the distance of each cell that was not NA from the nearest sea. Otherwise, the value "nosea" will be used to fill the whole raster, or if it is missing, the diagonal dimension of the raster will be used as a (bad) guess.
}
\value{
A raster of distances (in kilometres) from the sea, with the same properties as dtm.
}
\references{
Gunton et al. (submitted) "Predicting ground temperatures across European landscapes"
}
\author{
Richard Gunton
}

\seealso{
\code{
 \link{distance}
}
}
\examples{
require(raster)
dtm <- raster(ncol=60,nrow=60,ext=extent(56,57,-3,-2))
dtm[] <- 1
# Add a narrow band of NAs - e.g. a river:
dtm[1,] <- NA 	
# Have a broad band of NAs - representing sea:
dtm[,55:60] <- NA
dtm.coastal <- coastal(dtm)
plot(dtm.coastal)
}