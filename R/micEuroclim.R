micEuroclim <- function(dtm, coast, lat, cloud, day = 197, month, air, conv, 
    wind, air.mean, sw, soltime.sw = 12, localtime.sw, lon, merid = 0, 
    dst = 0, pred = c("mean_open", "mean_trees", "daily_open", 
        "daily_trees"), daily = 12, trans0 = 0.2, scale.threshold = 0.5, 
    cleanup = FALSE, reminders = TRUE, ...) 
{
    dai <- daily/2
    if (is.character(pred)) 
        pred <- match(pred[1], c("mean_open", "mean_trees", "daily_open", 
            "daily_trees"))
    if (!is.numeric(pred)) 
        stop("pred should be one of c(\"mean_open\",\"mean_trees\",\"daily_open\",\"daily_trees\") or a number.")
    if (pred < 1 || pred >= 5) 
        stop("pred should be between 1 and 4.")
    if (reminders && !missing(air) && length(air) == 1) 
        cat("air input taken to be", c("monthly mean.", "monthly mean.", 
            "at 13.40.", "morning mean;")[pred], fill = TRUE)
    if (reminders && !missing(conv) && length(conv) == 1) 
        cat("conv input taken to be", c("monthly mean.", "monthly mean.", 
            "morning mean.", "at 13.40.")[pred], fill = TRUE)
    if (reminders && !missing(wind) && length(wind) != length(air)) 
        warning("wind is not the same length as air, so will not be used.")
    if (missing(conv) && !missing(wind) && !missing(air) && !missing(air.mean) && 
        length(wind) == length(air)) {
        conv <- mean((air - air.mean) * sqrt(wind))
        if (reminders && length(air) > 1) 
            cat("\nwind and air input values taken to be", c("evenly spaced throughout the day.", 
                "evenly spaced throughout the day.", "evenly spaced from midnight to midday.", 
                "around 13.40 on different days.")[pred], fill = TRUE)
    }
    main.names <- c("int", "dtm", "air", "conv", "rad^0.5", "rad^2", 
        "rad", "cloud", "coast^0.5", "lat^0.5")
    data(centres, envir = environment())
    if (any(is.na(match(row.names(centres), main.names)))) 
        stop("row.names(centres) includes some terms that are not available.")
    if (!isLonLat(dtm)) {
        if (missing(lat) || (missing(coast) && !missing(air))) 
            stop("Dimensions of dtm are not specified as latitude & longitude; lat and coast must be supplied.")
        if (lat < 35 || lat > 70) 
            warning("Results may be unreliable: the model is only calibrated for Europe (c. 35 < lat < 70 deg, -10 < lon < 30 deg).")
    }
    else {
        ymin <- ymin(dtm)
        ymax <- ymax(dtm)
        if (ymin < 35 || ymax > 70) 
            warning("Results may be unreliable: the model is only calibrated for Europe (c. 35 < lat < 70 deg, -10 < lon < 30 deg).")
        if (scale.threshold < (ymax - ymin)) 
            lat <- setValues(dtm, rep(seq(ymax, ymin, length.out = nrow(dtm)), 
                each=ncol(dtm)))
        else lat <- (ymin + ymax)/2
    }
    if (missing(coast)) 
        coast <- coastal(dtm, ...)
    else {
        if (is.raster(coast)) {
            if (projection(dtm) == "NA") {
                coast <- cellStats(coast, mean, na.rm = TRUE)
                warning("Projection of dtm is not specified; mean value of coast will be used for coastal distance.")
            }
            else compareRaster(coast, dtm)
        }
    }
    mids <- c(15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 
        258, 288.5, 319, 349.5)
    if (missing(month)) 
        month <- 1 + day%/%30
    if (missing(day)) 
        day <- mids[month]
    elevation <- cellStats(dtm, mean, na.rm = TRUE)
    cloud0 <- centres["cloud", pred]
    transmit <- 1 - cloud0 * (1 - trans0)
    if (!missing(cloud)) {
        if(max(as.vector(cloud),na.rm=TRUE) > 1)  stop("cloud should be between 0 and 1.")
        transmit <- 1 - cloud * (1 - trans0)
    }
    if (!missing(sw)) {
        transmit <- trans(sw = sw, elevation = elevation, lat = lat, 
            day = day, soltime = soltime.sw, localtime = localtime.sw, 
            lon = lon, merid = merid, dst = dst, trans0 = trans0, 
            ...)
        cloud0 <- (1 - transmit)/(1 - trans0)
    }
    if (missing(cloud)) 
        cloud <- cloud0
    time <- 13.66
    rad <- srad.raster(dtm = dtm, trans = transmit, elevation = elevation, 
        lat = lat, day = day, soltime = time, lon = 0, merid = 0, 
        dst = 0, trans0 = trans0, ...)
    if (pred <= 2) {
        time <- time + 12/dai * c(-dai:-1, 1:(dai - 1))
        for (h in 1:length(time)) rad <- addLayer(rad, srad.raster(dtm = dtm, 
            trans = transmit, elevation = elevation, lat = lat, 
            day = day, soltime = time[h], lon = 0, merid = 0, 
            dst = 0, trans0 = trans0, ...))
        rad <- mean(rad)
    }
    if (pred >= 3) {
        if (!missing(sw)) 
            rad <- swrad2.raster(dtm = dtm, sw = sw, lat = lat, 
                day = day, elevation = elevation, soltime = time, 
                soltime.sw = soltime.sw, localtime.sw = localtime.sw, 
                lon = lon, merid = merid, dst = dst, trans0 = trans0, 
                ...)$swrad
        if (pred == 4) {
            time <- time - 12/dai * c(2:dai - 1)
            if (!missing(sw)) {
                for (h in 1:length(time)) suppressWarnings(rad <- addLayer(rad, 
                  swrad2.raster(dtm = dtm, sw = sw, elevation = elevation, 
                    lat = lat, day = day, soltime = time[h], 
                    soltime.sw = soltime.sw, localtime.sw = localtime.sw, 
                    lon = lon, merid = merid, dst = dst, trans0 = trans0, 
                    ...)$swrad))
            }
            else {
                for (h in 1:length(time)) rad <- addLayer(rad, 
                  srad.raster(dtm = dtm, trans = transmit, elevation = elevation, 
                    lat = lat, day = day, soltime = time[h], 
                    lon = 0, merid = 0, dst = 0, trans0 = trans0, 
                    ...))
            }
            rad <- mean(rad)
        }
    }
    SWrad <- cellStats(rad, mean)
    relative <- missing(air)
    int <- 1
    dtm <- dtm/1000
    if (missing(conv)) 
        conv <- centres["conv", pred]
    if (missing(air)) 
        air <- centres["air", pred]
    maindata <- do.call(c, as.list(parse(text = row.names(centres))))
    mains <- gsub("^", "", row.names(centres), fixed = TRUE)
    main.names <- gsub("^", "", main.names, fixed = TRUE)
    names(maindata) <- mains
    for (i in 1:length(maindata)) maindata[[i]] <- maindata[[i]] - 
        centres[i, pred]
    data(coefs, envir = environment())
    coefs.vars <- unique(unlist(strsplit(row.names(coefs), "*", 
        fixed = TRUE)))
    if (any(is.na(match(coefs.vars, main.names)))) 
        stop("row.names(coefs) includes some terms that are not available.")
    coefs <- subset(coefs, subset = coefs[, pred] != 0, select = pred)
    attach(maindata, warn.conflicts = FALSE)
    data <- do.call(c, as.list(parse(text = row.names(coefs))), 
        envir = as.environment(2))
    detach(maindata)
    rm(maindata)
    names(data) <- row.names(coefs)
    if (cleanup) 
        rad <- SWrad
    MCM <- 0
    for (i in 1:length(data)) MCM <- MCM + coefs[i, ] * data[[i]]
    if (relative) 
        MCM <- MCM - cellStats(MCM, mean, na.rm = TRUE)
    if (reminders) 
        cat("Relative "[relative], "Predictions for ", c("mean_open", 
            "mean_trees", "daily_open", "daily_trees")[pred], 
            " for month ", month, " / day ", day, "; mean vertical SW radiation = ", 
            SWrad, " W per sq. m.\n", sep = "")
    MCM <- list(microclim = MCM, rad = rad, cloud = cloud, trans = transmit, 
        coast = coast, air = air, conv = conv, day = day, call = match.call())
    return(MCM)
    class(MCM) <- "microclim"
}
