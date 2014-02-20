trans <-
function(sw,elevation=0,lat,day,soltime,localtime,lon=0,merid=0,dst=0,trans0=0.2,...)
{
    if(!missing(localtime))  soltime <- solartime(localtime,lat=lat,lon=lon,merid=merid,dst=dst,day=day)
    alt <- solalt(soltime,day,lat)
    if(alt < 10) warning("Solar altitude less than 10 degrees; transmission estimate may be unreliable.")
    beam <- beamrad(day,alt,elevation,...)
    diff <- diffuserad(day,alt,...)
    trans <- (sw - trans0 * diff) / (beam * sin(alt/180*pi) - trans0 * diff)
    trans[trans < 0] <- 0
    trans
}
