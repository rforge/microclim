coastal <-
function(dtm, nosea, notsea=2)
{
    maxDist <- pointDistance(c(xmax(dtm),ymax(dtm)), c(xmin(dtm),ymin(dtm)), longlat=TRUE)
    if(missing(nosea))  nosea <- round(maxDist)  else  nosea <- nosea*1000
    coast <- nosea
    fact <- ceiling(notsea*1000 / maxDist * sqrt(ncol(dtm)^2 + nrow(dtm)^2))
    if(is.na(cellStats(dtm,sum,na.rm=FALSE))) {
        df <- data.frame(id=NA, v=-1)
        dtm.ag <- if(fact > 1) aggregate(dtm, fun=median, fact=fact, na.rm=TRUE) else dtm
        if(is.na(cellStats(dtm.ag,sum,na.rm=FALSE))) {
            coast <- subs(dtm.ag, y=df, subsWithNA=TRUE)
            coast <- distance(coast)
            coast <- disaggregate(coast, fact=fact, method="bilinear")
            coast <- crop(coast, extent(dtm))
            if(!compareRaster(coast,dtm,stopiffalse=FALSE))  stop("Sorry - failed to obtain coastal distances. Please adjust extent of dtm and try again!")
        }
    }
    coast[coast < 0] <- 0
    coast /1000
}
