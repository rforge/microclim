solalt <-
function(time, day, lat)
{
t <- 0.261799 * (time - 12)
declin <- (pi * 23.5 / 180) * cos(2 * pi * ((day - 171) / 365.25))
Sinh <- sin(declin) * sin(lat * pi / 180) + cos(declin) * cos(lat * 3.141 / 180) * cos(t)
solalt <- (180 * atan(Sinh / sqrt(1 - Sinh * Sinh))) / pi
solalt
}
