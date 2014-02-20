beamrad <-
function(day,alt,elevation=0,linke=3)
{
solarconstant <- 1367
dayangle <- (2*pi*day)/365.25
e <- 1 + 0.03344 * cos(dayangle - 0.048869)
## G is extraterrestrial irradiance
G <- e * solarconstant

deltaalt <- 0.061359 * (0.1594 + 1.123 * alt + 0.065656 * alt^2) / (1+ 28.9344 * alt + 277.3971 * alt^2)
href <- alt + deltaalt
ppo <- exp(-1*elevation/8434.5)

## m is relative optical airmass
m <- ppo/(sin(href*(pi/180))+0.50572*(href+6.07995)^-1.6364)
m[alt<=1] <- 999

## dr is Rayleigh optical thickness at air mass m
dr <- 1/(6.6296 + 1.7513 * m - 0.1202 * m^2 + 0.0065 * m^3 - 0.00013 * m^4)
dr[m<=20] <- 1/(10.4 + 0.718 * m[m<=20])

## B is beam irradiance normal to the solar beam
B <- G * exp(-0.8662 * linke * m * dr)
B[alt<=1] <- 0
B
}
