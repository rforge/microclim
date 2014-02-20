srad.raster <-
function(dtm,trans,lat,day,elevation,soltime,localtime,lon,merid=0,dst=0,trans0=0.2,...)
{
    if(!missing(localtime)) {
        soltime <- solartime(localtime,lat=lat,lon=lon,merid=merid,dst=dst,day=day) 
    } else  if(missing(soltime))  stop("Either localtime or soltime must be specified.")
    if(missing(elevation)) elevation <- cellStats(dtm,mean,na.rm=TRUE)
    if(is.na(elevation)) elevation <- dtm
    if(missing(lat)) {
        if (!isLonLat(dtm))  stop("Dimensions of dtm are not specified as latitude & longitude; lat must be supplied.")
        lat <- (ymax(dtm) + ymin(dtm)) /2
    }
    slope <- terrain(dtm, opt='slope') /pi*180
    aspect <- terrain(dtm, opt='aspect') /pi*180
    alt <- solalt(soltime,day=day,lat=lat)
    azi <- solazi(soltime,day=day,lat=lat)
    beam <- beamrad(day=day,alt=alt,elevation=elevation,...)
    diff <- diffuserad(day=day,alt=alt,...)
    slopeflux <- beam * solarindex(slope=slope,aspect=aspect,azi=azi,alt=alt) + diff
    srad <- trans0 * diff + trans *(slopeflux - trans0 * diff)
    srad
}
