diffuserad <-
function(day,alt,linke=3)
{
solarconstant <- 1367
dayangle <- (2*pi*day)/365.25
e <- 1 + 0.03344 * cos(dayangle - 0.048869)
## G is extraterrestrial irradiance
G <- e * solarconstant
Tn <- -0.015843 + 0.030543 * linke + 0.0003797 * linke^2
Ai <- 0.26463 - 0.061581 * linke + 0.0031408 * linke^2
if (Ai < 0.022) Ai <- 0.0022/Tn
Aii <- 2.04020 + 0.018945 * linke - 0.011161 * linke^2
Aiii <- -1.3025 + 0.039231 * linke + 0.0085079 * linke^2
Fd <- Ai + Aii * sin(alt*pi/180) + Aiii * sin(alt*pi/180)^2
D <- G * Tn * Fd
D[D<=0] <- 0
D
}
