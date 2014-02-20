solazi <-
function(time, day, lat)
{
t <- 0.261799 * (time - 12)
declin <- (pi * 23.5 / 180) * cos(2 * pi * ((day - 171) / 365.25))
Sinh <- sin(declin) * sin(lat * pi / 180) + cos(declin) * cos(lat * pi / 180) * cos(t)
h <- (atan(Sinh / sqrt(1 - Sinh * Sinh)))
Sinazi <- cos(declin) * sin(t) / cos(h)
cosazi <- (sin(lat * pi / 180) * cos(declin) * cos(t) - cos(pi * lat / 180) * sin(declin)) / sqrt((cos(declin) * sin(t)) ^ 2 + (sin(pi * lat / 180) * cos(declin) * cos(t) - cos(pi * lat / 180) * sin(declin)) ^ 2)
solazi <- 180 + (180 * atan(Sinazi / sqrt(1 - Sinazi * Sinazi))) / pi
solazi[cosazi<0 & Sinazi<0] <- 180-solazi[cosazi<0 & Sinazi<0]
solazi[cosazi<0 & Sinazi>=0] <- 540-solazi[cosazi<0 & Sinazi>=0]
solazi
}
