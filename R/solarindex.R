solarindex <-
function(slope,aspect,azi,alt)
{
alt <- alt * (pi/180)
zen <- pi/2 - alt
azi <- azi * (pi/180)
s <- slope * (pi/180)
a <- aspect * (pi/180)
index <- cos(zen) * cos(s) + sin(zen) * sin(s) * cos(azi - a)
index[alt<=0] <- 0
index[index<=0] <- 0
index
}
