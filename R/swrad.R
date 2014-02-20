swrad <-
function(sw, slope, aspect, lat, day, elevation=0, soltime, localtime, lon, merid=0, dst=0, trans0=0.2, ...)
{
    if(missing(soltime))  soltime <- solartime(localtime,lat,lon,merid,dst,day)
    alt <- solalt(soltime,day,lat)
    if(alt < 20)  warning("Solar altitude at the time given is less than 20 degrees; results may be unreliable.")
    azi <- solazi(soltime,day,lat)      
    beam <- beamrad(day,alt,elevation,...)
    diff <- diffuserad(day,alt,...)
    swflux <- beam * solarindex(0,0,azi,alt)  + diff
    slopeflux <- beam * solarindex(slope,aspect,azi,alt) + diff
    trans <- (sw - trans0 * diff) / (swflux - trans0 * diff)
    swrad <- trans0 * diff + trans * (slopeflux - trans0 * diff)
    list(swrad=swrad, trans=trans)
}
