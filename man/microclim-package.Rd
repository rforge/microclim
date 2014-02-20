\name{microclim-package}
\alias{microclim-package}
\alias{microclim}
\docType{package}
\title{
Predicting microclimate
}
\description{
A suite of functions to calculate variables useful in microclimate modelling. Currently the main predictions are for solar radiation (globally) and ground temperature (in Europe).
}
\details{
\tabular{ll}{
Package: \tab microclim\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-03\cr
License: \tab GPL-2\cr
}
The principal functions at present are \code{\link{srad.raster}} and \code{\link{swrad}}, both for predicting short-wave radiation, and \code{\link{micEuroclim}}, for predicting ground temperatures. \code{micEuroclim} calls most of the other functions in this package.
}
\author{
Jon Bennie, Richard Gunton
Maintainer: Richard Gunton <fbsrmg@leeds.ac.uk>
}
\references{
Gunton et al. (submitted) "Predicting ground temperatures across European landscapes"
}

\keyword{ package }

\examples{
if(require(raster)){
data(volcano)
dtm <- raster(volcano,xmn=0,xmx=610,ymn=0,ymx=870,crs="+proj=equirectangular")

# Specify some NA cells to constitute sea:
dtm[78:87,52:61] <- NA

par(mfrow=c(2,2))
plot(dtm,col=terrain.colors(255))

coastdist <- coastal(dtm,notsea=0.1)
plot(coastdist, col=topo.colors(255))

srad <- srad.raster(dtm, lat=36.9, trans=0.7, day=172, soltime=12)
plot(srad)

micro <- micEuroclim(dtm, lat=36.9, pred=3, day=172)
# ...imagining this site were in southern Spain.
plot(micro$microclim)
}
}