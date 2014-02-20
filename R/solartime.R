solartime <-
function(localtime, lat, lon, merid=0, dst=0, day)
{
Bn <- 2 * pi * (day - 81) / 364
eot <- 9.87 * sin(2 * Bn) - 7.53 * cos(Bn) - 1.5 * sin(Bn)
solartime <- localtime + (4 / 60) * ((merid - lon) * pi / 180) + (eot / 60) - dst
solartime
}
