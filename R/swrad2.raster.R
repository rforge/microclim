swrad2.raster <-
function(dtm,sw,lat,day,elevation,soltime,soltime.sw,localtime,localtime.sw=localtime,lon,merid=0,dst=0,trans0=0.2,...)
{
    if(!missing(localtime)) {
        soltime <- solartime(localtime,lat=lat,lon=lon,merid=merid,dst=dst,day=day)
    } else  if(missing(soltime))  stop("Either localtime or soltime must be specified.")
    if(!missing(localtime.sw))  soltime.sw <- solartime(localtime.sw,lat=lat,lon=lon,merid=merid,dst=dst,day=day)
    if(missing(elevation))  elevation <- cellStats(dtm,mean,na.rm=TRUE)
    if(is.na(elevation))  elevation <- dtm
    if(missing(lat)) {
        if (!isLonLat(dtm))  stop("Dimensions of dtm are not specified as latitude & longitude; lat must be supplied.")
        lat <- (ymax(dtm) + ymin(dtm)) /2
    }
    slope <- terrain(dtm, opt='slope') /pi*180
    aspect <- terrain(dtm, opt='aspect') /pi*180
    alt <- solalt(soltime,day=day,lat=lat)
    alt.sw <- solalt(soltime.sw,day=day,lat=lat)
    if(alt.sw < 10)  warning("Solar altitude at the time given is less than 10 degrees; results may be unreliable.")
    azi <- solazi(soltime,day=day,lat=lat)
    beam.sw <- beamrad(day=day,alt=alt.sw,elevation=elevation,...)
    diff.sw <- diffuserad(day=day,alt=alt.sw,...)
    flux.sw <- beam.sw * sin(alt.sw/180*pi) + diff.sw
    trans <- (sw - trans0 * diff.sw) / (flux.sw - trans0 * diff.sw)
    beam <- beamrad(day=day,alt=alt,elevation=elevation,...)
    diff <- diffuserad(day=day,alt=alt,...)
    slopeflux <- beam * solarindex(slope=slope,aspect=aspect,azi=azi,alt=alt) + diff
    swrad <- trans0 * diff + trans *(slopeflux - trans0 * diff)
    list(swrad=swrad, trans=trans)
}
