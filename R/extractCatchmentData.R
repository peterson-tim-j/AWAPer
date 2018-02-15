#' Extracts daily catchment mean and variance from netCDF files of the Bureau of Meteorology (Australia) national gridded climate data.
#'
#' \code{extractCatchmentData} extracts catchment average climate data from netCDF files containing Australian climate data.
#'
#' This function extracts the AWAP climate data for each point or polygon, and for the latter the daily spatial mean and variance of
#' each climate metric are calculated. The calculation of the spatial mean uses the fraction of each AWAP grid cell within the polygon.
#' The variance calculation does not use the fraction of the grid cell and returns NA if there are <2 grid cells in the catchment boundary.
#' Prior to the catchment averaging and variance, Morton's areal potential ET (PET) is also calculated; after which the mean and
#' variance PET is calculated. Morton's PET is calculated using the ET.MortonCRAE() function from the Evapotranspiration package at
#' a monthly time-step and using the AWAP solar radiation. For both points and polygons, the monthly PET estimate is then interpolated
#' using a spline to a daily time step (using zoo:na.spline()) and then constrained to >=0.
#'
#' Also, when "catchments" is points (not polygons), then the netCDF grids are interpolate using bilinear interpolation of
#' the closest 4 grid cells.
#'
#' Lastly, data is extracted for all time points and no temporal infilling is undertaken if the grid cells are blank.
#'
#' @param ncdfFilename is a full file name (as string) to the netCDF file.
#' @param ncdfSolarFilename is the full file name (as string) to the netCDF file.
#' @param extractFrom is a date string specifying the start date for data extraction. The default is \code{"1900-1-1"}.
#' @param extractTo is a date string specifying the end date for the data extraction. The default is today's date as YYYY-MM-DD.
#' @param getPrecip logical variable for extracting precipitation. Default is \code{TRUE}.
#' @param getTmin logical variable for extracting Tmin. Default is \code{TRUE}.
#' @param getTmax logical variable for extracting Tmax. Default is \code{TRUE}.
#' @param getVprp logical variable for extracting vapour pressure. Default is \code{TRUE}.
#' @param getSolarrad logical variable for extracting solar radiation. Default is \code{TRUE}.
#' @param getMortonsPET logical variable for calculating Morton's potential ET. Note, to calculate set \code{getTmin=T}, \code{getTmax=T},
#' \code{getVprp=T} and \code{getSolarrad=T}. Default is \code{TRUE}.
#' @param DEM is either the full file name to a ESRI ASCII grid (as lat/long and using GDA94) or a raster class grid object. The DEM is used
#' for the calculation of Morton's PET. The Australian 9 second DEM can be loaded using \code{data(DEM_9s)}. For details see
#' \url{https://www.data.gov.au/dataset/geodata-9-second-dem-and-d8-digital-elevation-model-version-3-and-flow-direction-grid-2008}
#' @param catchments is either the full file name to an ESRI shape file of points or polygons (latter assumed to be catchment boundaries) or a shape file
#' already imported using readShapeSpatial(). Either way the shape file must be in long/lat (i.e. not projected) use the ellipsoid GRS 80.
#'
#' @return
#' When "catchments" are polygons, the returned variable is list variables containing two data.frames, one of the catchment average daily climate
#' metrics and another of the catchment variance daily climate metrics.
#'
#' When "catchments" are points, the returned variable is a data.frame containing daily climate data at each point.
#'
#' @seealso
#' \code{\link{makeNetCDF_file}} for building the NetCDF files of daily climate data.
#'
#' @examples
#' # Load required packages.
#' library(sp);library(raster);library(chron);library(ncdf4);
#' library(maptools);library(Evapotranspiration);library(AWAPer)
#'
#' # Download the 9 second Australian DEM.
#' DEM_9s = getDEM()
#'
#' # Load example cacthment boundaries.
#' data("catchments")
#'
#' # Extract all climate data and calculate Morton's PET.
#' # Note, the input "catchments" can also be a file to a ESRI shape file.
#' climateData = extractCatchmentData(catchments=catchments,DEM=DEM_9s)
#'
#' # Extract the daily catchment average data.
#' climateDataAvg = climateData$catchmentAvg
#'
#' # Extract the daily catchment variance data.
#' climateDataVar = climateData$catchmentVar
#'
#' # Export data to .csv files
#' write.csv(climateDataAvg,'warrionClimateAvg.csv')
#' write.csv(climateDataVar,'warrionClimateVar.csv')
#'
#' @export
extractCatchmentData <- function(
  ncdfFilename=file.path(getwd(),'AWAP.nc'),
  ncdfSolarFilename=file.path(getwd(),'AWAP_solar.nc'),
  extractFrom = as.Date("1900-01-01","%Y-%m-%d"),
  extractTo  = as.Date(Sys.Date(),"%Y-%m-%d"),
  getPrecip = TRUE,
  getTmin = TRUE,
  getTmax = TRUE,
  getVprp = TRUE,
  getSolarrad = TRUE,
  getMortonsPET = TRUE,
  DEM='',
  catchments='')  {

  # Build time points to update
  if (is.character(extractFrom))
    extractFrom = as.Date(extractFrom,'%Y-%m-%d');
  if (is.character(extractTo))
    extractTo = as.Date(extractTo,'%Y-%m-%d');
  if (extractFrom >= extractTo)
    stop('The update dates are invalid. extractFrom must be prior to extractTo')
  timepoints2Extract = seq( as.Date(extractFrom,'%Y-%m-%d'), by="day", to=as.Date(extractTo,'%Y-%m-%d'))
  if (length(timepoints2Extract)==0)
    stop('The dates to extract produce a zero vector of dates of zero length. Check the inputs dates are as YYYY-MM-DD')


  # If all required inputs are to be extracted for Mortions PET
  if (getMortonsPET) {
    if (!getTmin || !getTmax || !getVprp || !getSolarrad)
      stop('Calculation of Mortons PET requires extractions of Tmin, Tmax, vapor pressure and solar radiation.')
  }

  # Open NetCDF grids
  awap <- nc_open(ncdfFilename)
  if (getSolarrad)
    awap_solar <- nc_open(ncdfSolarFilename)

  # Checking the DEM
  if (getMortonsPET) {
    if (is.character(DEM)) {
      if (!file.exists(DEM))
        stop(paste('The following input file for the DEM could not be found:',DEM))

      # Read in DEM
      DEM = read.asciigrid(DEM,  proj4string = CRS("+proj=longlat +ellps=GRS80"))

      # Convert to Raster object
      DEM = raster(DEM)
    } else if (!is(DEM_9s,"RasterLayer")) {
      stop('The input for the DEM is not a raster layer object.')
    }
  }

  # Open file with polygons
  if (is.character(catchments)) {
    if (!file.exists(catchments))
      stop(paste('The following input file for the catchments could not be found:',catchments))

    # Read in polygons
    catchments = readShapeSpatial(catchments, force_ring=TRUE)

  } else if (!is(catchments,"SpatialPolygonsDataFrame") && !is(catchments,"SpatialPointsDataFrame")) {
    stop('The input for "catchments" must be a file name to a shape file or a SpatialPolygonsDataFrame or a SpatialPointsDataFrame object.')
  }

  # Check projection of catchments
  if (is.na(proj4string(catchments))) {
    message('WARNING: The projection string of the catchment boundaries is NA. Setting to +proj=longlat +ellps=GRS80.')
    proj4string(catchments) = '+proj=longlat +ellps=GRS80'
  }
  if (proj4string(catchments) != '+proj=longlat +ellps=GRS80') {
    message('WARNING: The projection string of the catchment boundaries does not appear to be +proj=longlat +ellps=GRS80. Attempting to transform coordinates...')
    catchments = spTransform(catchments,CRS('+proj=longlat +ellps=GRS80'))
  }

  # Check if catchments are points or a polygon.
  isCatchmentsPolygon=TRUE;
  if (is(catchments,"SpatialPointsDataFrame")) {
    isCatchmentsPolygon=FALSE;
  }

  # Get netCDF geometry
  timePoints <- ncvar_get(awap, "time")
  nTimePoints = length(timePoints);
  Long <- ncvar_get(awap, "Long")
  Lat <- ncvar_get(awap, "Lat")
  if (getSolarrad) {
    timePoints_solar <- ncvar_get(awap_solar, "time")
    nTimePoints_solar = length(timePoints_solar);
    Long_solar <- ncvar_get(awap_solar, "Long")
    Lat_solar <- ncvar_get(awap_solar, "Lat")
  }

  # Convert the ncdf time to R time.
  # This is achieved by spliting the time units string into fields.
  # Adapted form http://geog.uoregon.edu/bartlein/courses/geog607/Rmd/netCDF_01.htm.
  tunits <- ncatt_get(awap, "time", "units")
  tustr <- strsplit(tunits$value, " ")
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth = as.integer(unlist(tdstr)[2])
  tday = as.integer(unlist(tdstr)[3])
  tyear = as.integer(unlist(tdstr)[1])
  timePoints_R = as.Date(chron(timePoints, origin = c(month=tmonth, day=tday, year=tyear),format=c(dates = "y-m-d")));
  if (getSolarrad) {
    tunits <- ncatt_get(awap_solar, "time", "units")
    tustr <- strsplit(tunits$value, " ")
    tdstr <- strsplit(unlist(tustr)[3], "-")
    tmonth = as.integer(unlist(tdstr)[2])
    tday = as.integer(unlist(tdstr)[3])
    tyear = as.integer(unlist(tdstr)[1])
    timePoints_solar_R = as.Date(chron(timePoints_solar, origin = c(tmonth, tday, tyear),format=c(dates = "y-m-d")));
  }

  # Write summary of Net CDF data.
  message('Extraction data summary:');
  message(paste('    NetCDF non-solar radiation climate data exists from ',format.Date(min(timePoints_R),'%Y-%m-%d'),' to ', format.Date(max(timePoints_R),'%Y-%m-%d')));
  if (getSolarrad)
    message(paste('    NetCDF solar radiation data exists from ',format.Date(min(timePoints_solar_R),'%Y-%m-%d'),' to ', format.Date(max(timePoints_solar_R),'%Y-%m-%d')));

  # Limit the extraction time points to the data range
  extractFrom = max(c(extractFrom,min(timePoints_R)));
  if (getSolarrad) {
    extractTo = min(c(extractTo,max(timePoints_R),max(timePoints_solar_R)));
  } else {
    extractTo = min(c(extractTo,max(timePoints_R)));
  }

  # Recalculate the time points to extract.
  timepoints2Extract = seq( as.Date(extractFrom,'%Y-%m-%d'), by="day", to=as.Date(extractTo,'%Y-%m-%d'))

  message(paste('    Data will be extracted from ',format.Date(extractFrom,'%Y-%m-%d'),' to ', format.Date(extractTo,'%Y-%m-%d'),' at ',length(catchments),' catchments '));
  message('Starting data extraction:')

  # Get one netCDF layer.
  precipGrd = raster(ncdfFilename, band=nTimePoints, varname='precip',lvar=3)
  if (getSolarrad) {
    solarGrd = raster(ncdfSolarFilename, band=nTimePoints_solar, varname='solarrad',lvar=3)
  }

  # Build a matrix of catchment weights, lat longs, and a loopup table for each catchment.
  message('... Building catchment weights:');
  if (isCatchmentsPolygon) {

    w.all = c();
    longLat.all = matrix(NA,0,2)
    catchment.lookup = matrix(NA,length(catchments),2);

    for (i in 1:length(catchments)) {
      if (i%%10 ==0 ) {
        message(paste('   ... Building weights for catchment ', i,' of ',length(catchments)));
        removeTmpFiles(h=0)
      }
      w = rasterize(catchments[i,], precipGrd,getCover=T)

      # Extract the mask values (i.e. fraction of each grid cell within the polygon.
      w2 = getValues(w);
      filt = w2>0
      wLongLat = coordinates(w)[filt,]
      w=w[filt]

      # Normalise the weights
      w = w/sum(w);

      # Add to data set of all catchments
      if (length(w.all)==0) {
        catchment.lookup[i,] = c(1,length(w));
        w.all = w;
        longLat.all = wLongLat;
      } else {
        catchment.lookup[i,] = c(length(w.all)+1,length(w.all)+length(w));
        w.all = c(w.all, w)
        longLat.all = rbind(longLat.all, wLongLat);
      }
    }
  } else {

    # For point data, set weights to 1 and coordinates from point locations
    w.all = rep(1,length(catchments))
    longLat.all = cbind(as.numeric(coordinates(catchments)[,1]),as.numeric(coordinates(catchments)[,2]))
    catchment.lookup = cbind(rep(1,length(catchments)),rep(1,length(catchments)));
  }

  removeTmpFiles(h=0)

  if (getSolarrad && getMortonsPET) {
    message('... Extracted DEM elevations.')
    DEMpoints <- raster::extract(DEM, longLat.all)
    if (any(is.na(DEMpoints))) {
      warning('NA DEM values were derived. Please check the projections of the DEM and shape file.')
    }
  }

  # Initialise the outputs
  precip = matrix(NA,length(timepoints2Extract),length(w.all))
  tmin =  precip;
  tmax =  precip;
  vprp =  precip;
  extractYear = c();
  extractMonth = c();
  extractDay = c();

  # Initialise the outputs. NOTE: the matrix is initiliased to have the same
  # number of rows as the non-solar data.
  if (getSolarrad) {
    solarrad = matrix(NA,length(timepoints2Extract),length(w.all))
    extractYear_solar = c();
    extractMonth_solar = c();
    extractDay_solar = c();
  }

  # Set the interpolation method.
  if (isCatchmentsPolygon) {
    interpMethod='simple';
  } else {
    interpMethod='bilinear';
  }

  message('... Starting to extract data across all catchments:')
  for (j in 1:length(timepoints2Extract)){

    # Find index to the date to update within the net CDF grid
    ind = as.integer(difftime(timepoints2Extract[j], as.Date("1900-1-1",'%Y-%m-%d'),units = "days" ))+1

    if (j%%10000 ==0 ) {
      message(paste('    ... Extracting data for time point ', j,' of ',length(timepoints2Extract)));
    }
    if (getPrecip)
      precip[j,1:length(w.all)]  <- raster::extract(raster(ncdfFilename, band=ind, varname='precip',lvar=3, method=interpMethod), longLat.all)
    if (getTmin)
      tmin[j,1:length(w.all)]  <- raster::extract(raster(ncdfFilename, band=ind, varname='tmin',lvar=3, method=interpMethod), longLat.all)
    if (getTmax)
      tmax[j,1:length(w.all)]  <- raster::extract(raster(ncdfFilename, band=ind, varname='tmax',lvar=3, method=interpMethod), longLat.all)
    if (getVprp)
      vprp[j,1:length(w.all)]  <- raster::extract(raster(ncdfFilename, band=ind, varname='vprp',lvar=3, method=interpMethod), longLat.all)

    # Get date of extracted grid.
    extractDate = getZ(raster(ncdfFilename, band=ind,lvar=3));
    extractYear[j] = as.integer(format(as.Date(extractDate,'%Y-%m-%d'),'%Y'));
    extractMonth[j] = as.integer(format(as.Date(extractDate,'%Y-%m-%d'),'%m'));
    extractDay[j] = as.integer(format(as.Date(extractDate,'%Y-%m-%d'),'%d'));

    # Extract the non-solar data
    #----------------------------------------------------------------------------------------
    if (getSolarrad) {
      # Find index to the date to update within the net CDF grid
      ind = as.integer(difftime(timepoints2Extract[j], as.Date("1990-1-1",'%Y-%m-%d'),units = "days" ))+1

      if (ind>0) {
        # Get ncdf grid
        solarrad[j,1:length(w.all)]  <- raster::extract(raster(ncdfSolarFilename, band=ind, varname='solarrad',lvar=3, method=interpMethod), longLat.all)

        # Get date of extracted grid.
        extractDate = getZ(raster(ncdfSolarFilename, band=ind,lvar=3));
        extractYear_solar[j] = as.integer(format(as.Date(extractDate,'%Y-%m-%d'),'%Y'));
        extractMonth_solar[j] = as.integer(format(as.Date(extractDate,'%Y-%m-%d'),'%m'));
        extractDay_solar[j] = as.integer(format(as.Date(extractDate,'%Y-%m-%d'),'%d'));
      }
    }
  }

  # Close netCDF connection
  if (!getPrecip || !getTmin || !getTmax || !getVprp)
    nc_close(awap)
  if (getSolarrad)
    nc_close(awap_solar)


  if (getSolarrad) {
    # Set non-senible values to NA
    solarrad[solarrad <0] = NA;
    solarrad_interp = solarrad;

    # Calculate the average dailysolar radiation for each day of the year.
    message('... Calculating mean daily solar radiation <1990-1-1')
    monthdayUnique = sort(unique(extractMonth*100+extractDay_solar));
    day = as.integer(format(timepoints2Extract, "%d"));
    month = as.integer(format(timepoints2Extract, "%m"));
    monthdayAll = month*100+day;
    solarrad_avg = matrix(NA, length(monthdayUnique), length(w.all));
    for (j in 1:length(monthdayUnique)) {
      ind = monthdayAll==monthdayUnique[j];
      if (sum(ind)==1) {
        solarrad_avg[j,] = solarrad[ind,];
      } else {
        solarrad_avg[j,] = apply(na.omit(solarrad[ind,]),2,mean)
      }
    }

    # Assign the daily average solar radiation to each day prior to 1 Jan 1990
    for (j in 1:length(timepoints2Extract)) {
      if (timepoints2Extract[j]<as.Date('1990-1-1','%Y-%m-%d')) {
        ind = monthdayUnique==monthdayAll[j];
        solarrad_interp[j,] = solarrad_avg[ind,]
      }
    }

    # Linearly interpolate time points without a solar radiation value.
    message('... Linearly interpolating gaps in daily solar.')
    for (j in 1:length(w.all)) {
      filt = is.na(solarrad_interp[,j])
      x = 1:length(timepoints2Extract);
      xpred = x[filt];
      x = x[!filt];
      y = solarrad_interp[!filt,j];
      ypred=approx(x,y,xpred,method='linear', rule=2);
      solarrad_interp[filt,j] = ypred$y;
    }
  }

  # Loop though each catchment and calculate the catchment averager and variance.
  if (getMortonsPET) {
    # load ET package constants
    data("constants")
  }
  message('... Calculating catchment weighted daily data.')
  for (i in 1:length(catchments)) {

    if (j%%100 ==0 ) {
      message(paste('    ... Calculating catchment ', i,' of ',length(catchments)));
    }



    # Calculate catchment stats and add data to the data.frame for all catchments
    #----------------------------------------------------------------------------------------
    catchmentAvgTmp = data.frame(CatchmentID=catchments[i,1],year=extractYear,month=extractMonth,day=extractDay);
    catchmentVarTmp = data.frame(CatchmentID=catchments[i,1],year=extractYear,month=extractMonth,day=extractDay);

    # Get the weights for the catchment
    ind = catchment.lookup[i,1]:catchment.lookup[i,2]
    w = w.all[ind]


    # Calculate Morton's PET at each grid cell and time point. NOTE, va is divided by 10 to go from hPa to Kpa
    if (getMortonsPET) {
      message('        ... Calculating Mortons areal PET unsing Evaportranspiration package CREA function.')

      # Initialise outputa
      mortAPET = matrix(NA,length(timepoints2Extract),length(ind))

      # Loop through each grid cell
      k=0;
      for (j in ind) {

        k=k+1
        if (k%%25 ==0 ) {
          message(paste('           ... Calculating PET for grid cell', k,' of ',length(w)));
        }

        # Check lat, Elev and precip are finite scalers.
        if (any(!is.finite(precip[,j])) || !is.finite(DEMpoints[j]) || !is.finite(longLat.all[j,2])) {
          message(paste('WARNING: Non-finite input values detected for catchment',i,' at grid cell',j))
          message(paste('   Elevation value:' ,DEMpoints[j]))
          message(paste('   Latitude value:' ,longLat.all[j,2]))
          mortAPET[,k] = NA;
          next
        }



        # Build data from of daily climate data
        dataRAW = data.frame(Year =  as.integer(format.Date(timepoints2Extract,"%Y")), Month= as.integer(format.Date(timepoints2Extract,"%m")), Day= as.integer(format.Date(timepoints2Extract,"%d")), Tmin=tmin[,j], Tmax=tmax[,j], Rs=solarrad_interp[,j], va=vprp[,j]/10.0, Precip=precip[,j])

        # Convert to required format for ET package
        dataPP=ReadInputs(c("Tmin","Tmax","Rs","Precip","va"),dataRAW,stopmissing = c(20,20,20))

        # Update constants for the site
        constants$Elev = DEMpoints[j]
        constants$lat = longLat.all[j,2]
        constants$lat_rad = longLat.all[j,2]/180.0*pi

        # Call  ET package
        results <- ET.MortonCRAE(dataPP, constants,est="potential ET", ts="monthly",solar="data",Tdew=FALSE, message='no');

        # Get the last day of each month
        last.day.month = as.Date(as.yearmon(time(results$ET.Monthly)), frac = 1)

        # Calculate average ET per day of each month
        days.per.month = as.integer(format.Date(as.Date(as.yearmon(time(results$ET.Monthly)), frac = 1),'%d'))
        monthly.ET.as.daily = zoo( as.numeric(results$ET.Monthly/days.per.month), last.day.month)

        # Spline interpolate Monthly average ET
        timepoints2Extract.as.zoo = zoo(NA,timepoints2Extract);
        mortAPET.tmp = na.spline(merge(monthly.ET.as.daily, dates=timepoints2Extract.as.zoo)[, 1])
        filt = time(mortAPET.tmp)>=start(timepoints2Extract.as.zoo) & time(mortAPET.tmp)<=end(timepoints2Extract.as.zoo)
        mortAPET.tmp = pmax(0.0, as.numeric( mortAPET.tmp));
        mortAPET.tmp = mortAPET.tmp[filt]
        mortAPET[,k] = mortAPET.tmp;

      }

    }

    # Check if there are enough grid cells to calculate the variance.
    calcVariance =  F;
    if (length(ind)>1)
      calcVariance =  T;

    # Apply weights
    if (getPrecip) {
      catchmentAvgTmp$precip_mm = apply(t(t(precip[,ind]) * w),1,sum,na.rm=TRUE);
      if (calcVariance) {
        catchmentVarTmp$precip_mm = apply(precip[,ind],1,var,na.rm=TRUE);
      } else {
        catchmentVarTmp$precip_mm = NA;
      }
    }
    if (getTmin) {
      catchmentAvgTmp$Tmin = apply(t(t(tmin[,ind]) * w),1,sum,na.rm=TRUE);
      if (calcVariance) {
        catchmentVarTmp$Tmin = apply(tmin[,ind],1,var,na.rm=TRUE);
      } else {
        catchmentVarTmp$Tmin = NA;
      }
    }
    if (getTmax) {
      catchmentAvgTmp$Tmax = apply(t(t(tmax[,ind]) * w),1,sum,na.rm=TRUE);
      if (calcVariance) {
        catchmentVarTmp$Tmax = apply(tmax[,ind],1,var,na.rm=TRUE);
      } else {
        catchmentVarTmp$Tmax = NA;
      }
    }
    if (getVprp) {
      catchmentAvgTmp$vprp = apply(t(t(vprp[,ind]) * w),1,sum,na.rm=TRUE);
      if (calcVariance) {
        catchmentVarTmp$vprp = apply(vprp[,ind],1,var,na.rm=TRUE);
      } else {
        catchmentVarTmp$vprp = NA;
      }
    }
    if (getSolarrad) {
      catchmentAvgTmp$solarrad = apply(t(t(solarrad[,ind]) * w),1,sum,na.rm=TRUE);
      catchmentAvgTmp$solarrad_interp = apply(t(t(solarrad_interp[,ind]) * w),1,sum,na.rm=TRUE);
      if (calcVariance) {
        catchmentVarTmp$solarrad = apply(solarrad[,ind],1,var,na.rm=TRUE);
        catchmentVarTmp$solarrad_interp = apply(solarrad_interp[,ind],1,var,na.rm=TRUE);
      } else {
        catchmentVarTmp$solarrad = NA;
        catchmentVarTmp$solarrad_interp = NA;
      }

    }
    if (getMortonsPET) {
      catchmentAvgTmp$MortonsPET_mm = apply(t(t(mortAPET) * w),1,sum,na.rm=TRUE);
      if (calcVariance) {
        catchmentVarTmp$MortonsPET_mm = apply(mortAPET,1,var,na.rm=TRUE);
      } else {
        catchmentVarTmp$MortonsPET_mm = NA;
      }
    }


    if (i==1) {
      catchmentAvg = catchmentAvgTmp;
      catchmentVar = catchmentVarTmp;
    } else {
      catchmentAvg = rbind(catchmentAvg,catchmentAvgTmp)
      catchmentVar = rbind(catchmentVar,catchmentVarTmp)
    }
  }
  # end for-loop

  if (isCatchmentsPolygon) {
    return(list(catchmentAvg=catchmentAvg, catchmentVar=catchmentVar))
  } else {
    return(catchmentAvg)
  }


  message('Data extraction FINISHED..')
}
