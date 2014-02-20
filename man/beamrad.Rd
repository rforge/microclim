\name{beamrad}
\alias{beamrad}
\alias{diffuserad}
\title{
Calculate short-wave solar irradiances
}
\description{
Calculates clear-sky direct beam solar irradiance, and diffuserad diffuse radiation, for specified conditions.
}
\usage{
beamrad(day, alt, elevation = 0, linke = 3)
diffuserad(day, alt, linke = 3)
}
\arguments{
  \item{day}{
day of the year, from 1 to 365 (366 in leap years)
}
  \item{alt}{
angle of sun above horizon in degrees
}
  \item{elevation}{
height above sea level of site (default is sea level)
}
  \item{linke}{
the Linke turbidity factor - an index of atmospheric turbidity. If in doubt, leave it at the default of 3.
}
}
\details{
If any of the arguments comprise more than a single value, they should be equal to each other in dimensions and extent and the output will be of the corresponding size. In particular, \code{alt} could be a raster.
}
\value{
Short-wave radiation values in Watts per square metre.
}
\references{
Suri, M. & Hofierka, J. (2007) Transactions in GIS, 8: 175-190
}
\author{
Jon Bennie
}
\examples{
beamrad(day=197, alt=74, elevation=0, linke=3)
beamrad(day=197, alt=0, elevation=0, linke=3)
diffuserad(day=197, alt=74, linke=3)
diffuserad(day=197, alt=0, linke=3)
}