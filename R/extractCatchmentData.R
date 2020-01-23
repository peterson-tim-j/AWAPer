#'
#' \code{extractCatchmentData} extracts catchment average climate data from netCDF files containing Australian climate data.
#'
#' @description
#' extractCatchmentData extracts the AWAP climate data for each point or polygon, and for the latter the daily spatial mean and variance (or user defined function) of
#' each climate metric is calculated.
#'
#' @details
#' The calculation of the spatial mean uses the fraction of each AWAP grid cell within the polygon.
#' The variance calculation (or user defined function) does not use the fraction of the grid cell and returns NA if there are <2 grid cells in the catchment boundary.
#' Prior to the catchment averaging and variance, evapotranspiration (ET) can also calculated; after which the mean and
#' variance PET is calculated.
#'
#' The data extraction will by default be undertaken from 1/1/1900 to yesterday, even if the netCDF grids were only
#' built for a subset of this time period. If the latter situation applies, it is recommended that the extraction start
#' and end dates are input by the user.
#'
#' The ET can be calculated using one of eight methods at a user defined time-step. When the time-step is monthly or annual then
#' the ET estimate is linearly interpolated to a daily time step (using zoo:na.spline()) and then constrained to >=0. In calculating ET, the input data
#' is pre-processed using Evapotranspiration::ReadInputs() such that missing days, missing enteries and abnormal values are interpolated
#' (by default) with the former two interpolated using the "DoY average", i.e. replacement with same day-of-the-year average. Additionally, when AWAP solar
#' radiation is required for the ET function, data is only available from 1/1/1990. To derive ET values <1990, the average solar radiation for each day of the year from
#' 1/1/990 to "extractTo" is derived (i.e. 365 values) and then applied to each day prior to 1990. Importantly, the estimates of ET <1990
#' are dependent upon the end date extracted. Re-running the estimation of ET with a later extractTo data will change the estimates of ET
#' prior to 1990.
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
#' @param getET logical variable for calculating Morton's potential ET. Note, to calculate set \code{getTmin=T}, \code{getTmax=T},
#' \code{getVprp=T} and \code{getSolarrad=T}. Default is \code{TRUE}.
#' @param DEM is either the full file name to a ESRI ASCII grid (as lat/long and using GDA94) or a raster class grid object. The DEM is used
#' for the calculation of Morton's PET. The Australian 9 second DEM can be loaded using \code{getDEM()}.
#' @param catchments is either the full file name to an ESRI shape file of points or polygons (latter assumed to be catchment boundaries) or a shape file
#' already imported using readShapeSpatial(). Either way the shape file must be in long/lat (i.e. not projected) use the ellipsoid GRS 80.
#' @param spatial.function.name character string for the function name applied to estimate the daily spatial spread in each variable. The default is \code{var}.
#' @param ET.function character string for the evapotranspiration function to be used. The methods that can be derived from the AWAP data are are \code{\link[Evapotranspiration]{ET.Abtew}},
#' \code{\link[Evapotranspiration]{ET.HargreavesSamani}}, \code{\link[Evapotranspiration]{ET.JensenHaise}}, \code{\link[Evapotranspiration]{ET.Makkink}}, \code{\link[Evapotranspiration]{ET.McGuinnessBordne}}, \code{\link[Evapotranspiration]{ET.MortonCRAE}} ,
#' \code{\link[Evapotranspiration]{ET.MortonCRWE}}, \code{\link[Evapotranspiration]{ET.Turc}}. Default is \code{\link[Evapotranspiration]{ET.MortonCRAE}}.
#' @param ET.Mortons.est character string for the type of Mortons Et estimate. For \code{ET.MortonCRAE}, the options are \code{potential ET},\code{wet areal ET} or \code{actual areal ET}.
#'  For \code{ET.MortonCRWE}, the options are \code{potential ET} or \code{shallow lake ET}. The default is \code{potential ET}.
#' @param ET.Turc.humid logical variable for the Turc function using the humid adjustment.See \code{\link[Evapotranspiration]{ET.Turc}}. For now this is fixed at \code{F}.
#' @param ET.timestep character string for the evapotranpiration time step. Options are \code{daily},  \code{monthly}, \code{annual} but the options are dependent upon the chosen \code{ET.function}. The default is \code{monthly}.
#' @param ET.interp_missing_days T or F, indicating if missing days should be interpolated for PET calculation. Default is \code{T}. See \code{\link[Evapotranspiration]{ReadInputs}}
#' @param ET.interp_missing_entries T or F, indicating if missing data entries should be interpolated for PET calculation. Default is \code{T}. See \code{\link[Evapotranspiration]{ReadInputs}}
#' @param ET.interp_abnormal T or F, indicating if abnormal valuses should be interpolated for PET calculation. Default is \code{T}. See \code{\link[Evapotranspiration]{ReadInputs}}
#' @param ET.constants list of constants from Evapotranspiration package required for ET calculations. To get the data use the command \code{data(constants)}. Default is \code{list()}.
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
#' # The example shows how to extract and save data.
#' # For an additional example see \url{https://github.com/peterson-tim-j/AWAPer/blob/master/README.md}
#' #---------------------------------------
#'
#' # Set dates for building netCDFs and extracting data.
#' startDate = as.Date("2000-01-01","%Y-%m-%d")
#' endDate = as.Date("2000-02-28","%Y-%m-%d")
#'
#' # Set names for netCDF files (in the system temp. directory).
#' ncdfFilename = tempfile(fileext='.nc')
#' ncdfSolarFilename = tempfile(fileext='.nc')
#'
#' # Build netCDF grids and over a defined time period.
#' \donttest{
#' file.names = makeNetCDF_file(ncdfFilename=ncdfFilename,
#'              ncdfSolarFilename=ncdfSolarFilename,
#'              updateFrom=startDate, updateTo=endDate)
#'
#' # Load example cacthment boundaries.
#' data("catchments")
#'
#' # Get the constants required for ET estimation.
#' data(constants,package='Evapotranspiration')
#'
#' # Set file name for DEM file to the system temp. directory.
#' DEM.file = tempfile(fileext='.asc')
#'
#' # Download and import the Australian 9 second DEM.
#' DEM_9s = getDEM(workingFolder=dirname(DEM.file),
#'          DEMfilename=basename(DEM.file))
#'
#' # Extract precip data.
#' # Note, the input "catchments" can also be a file to a ESRI shape file.
#' climateData = extractCatchmentData(ncdfFilename=file.names$ncdfFilename,
#'               ncdfSolarFilename=file.names$ncdfSolarFilename,
#'               extractFrom=startDate, extractTo=endDate,
#'               catchments=catchments,DEM=DEM_9s, ET.constants=constants)
#'
#' # Extract the daily catchment average data.
#' climateDataAvg = climateData$catchmentAvg
#'
#' # Extract the daily catchment variance data.
#' climateDataVar = climateData$catchmentvar
#' }
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
  getET = TRUE,
  DEM='',
  catchments='',
  spatial.function.name = 'var',
  ET.function = 'ET.MortonCRAE',
  ET.Mortons.est = 'potential ET',
  ET.Turc.humid=F,
  ET.timestep = 'monthly',
  ET.interp_missing_days=T,
  ET.interp_missing_entries=T,
  ET.interp_abnormal=T,
  ET.constants=list())  {

  # Open NetCDF grids
  awap <- ncdf4::nc_open(ncdfFilename)
  if (getSolarrad)
    awap_solar <- ncdf4::nc_open(ncdfSolarFilename)

  # Check if the required variable is within the netcdf file.
  netCDF.variables = names(awap$var)
  if (getTmin & !any(netCDF.variables=='tmin'))
    stop('getTmin is true but the netCDF file was not built with tmin data. Rebuild the netCDF file.')
  if (getTmax & !any(netCDF.variables=='tmax'))
    stop('getTmax is true but the netCDF file was not built with tmax data. Rebuild the netCDF file.')
  if (getVprp & !any(netCDF.variables=='vprp'))
    stop('getVprp is true but the netCDF file was not built with vprp data. Rebuild the netCDF file.')
  if (getPrecip & !any(netCDF.variables=='precip'))
    stop('getPrecip is true but the netCDF file was not built with precip data. Rebuild the netCDF file.')

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

  # Check ET inputs
  if (getET) {
    if (length(ET.constants)==0)
      stop('ET.constants must be input from Evapotranspiration package using the command: data(constants)')

    # Check the ET function is one of the accetable forms
    ET.function.all =c('ET.Abtew', 'ET.HargreavesSamani', 'ET.JensenHaise','ET.Makkink', 'ET.McGuinnessBordne','ET.MortonCRAE' , 'ET.MortonCRWE','ET.Turc')
    if (!any(ET.function == ET.function.all)) {
      stop(paste('The ET.function must be one of the following:',ET.function.all))
    }

    # Check the ET function is one of the accetable forms
    ET.timestep.all =c('daily', 'monthly', 'annual')
    if (!any(ET.timestep == ET.timestep.all)) {
      stop(paste('The ET.timestep must be one of the following:',ET.timestep.all))
    }

    # Check the appropriate time step is used.
    if ( (ET.function == 'ET.MortonCRAE' || ET.function == 'ET.MortonCRWE') && ET.timestep=='daily' ) {
      stop('The ET.timstep must be monthly or annual when using ET.MortonCRAE or ET.MortonCRWE')
    }

    # Build a data from of the inputs required for each ET function.
    ET.inputdata.req = data.frame(ET.function=ET.function.all, Tmin=rep(F,length(ET.function.all)), Tmax=rep(F,length(ET.function.all)), va=rep(F,length(ET.function.all)), Rs=rep(F,length(ET.function.all)), Precip=rep(F,length(ET.function.all)) )
    ET.inputdata.req[1,2:6] =  c(T,T,F,T,F)       #ET.Abtew
    ET.inputdata.req[2,2:6] =  c(T,T,F,F,F)       #ET.HargreavesSamani
    ET.inputdata.req[3,2:6] =  c(T,T,F,T,F)       #ET.JensenHaise
    ET.inputdata.req[4,2:6] =  c(T,T,F,T,F)       #ET.Makkink
    ET.inputdata.req[5,2:6] =  c(T,T,F,F,F)       #ET.McGuinnessBordne
    ET.inputdata.req[6,2:6] =  c(T,T,T,T,T)       #ET.MortonCRAE
    ET.inputdata.req[7,2:6] =  c(T,T,T,T,T)       #ET.MortonCRWE
    ET.inputdata.req[8,2:6] =  c(T,T,F,T,F)       #ET.Turc
    ind = which(ET.function == ET.function.all)
    ET.inputdata.filt = ET.inputdata.req[ind,]

    # Get list of required ET variable names
    ET.var.names = colnames(ET.inputdata.filt)[2:6]

    # If all required inputs are to be extracted for Mortions PET
    if (ET.inputdata.filt$Tmin[1] && !getTmin)
      stop('Calculation of ET for the given function requires extractions of tmin (i.e. set getTmin=T)')
    if (ET.inputdata.filt$Tmax[1] && !getTmax)
      stop('Calculation of ET for the given function requires extractions of tmax (i.e. set getTmax=T)')
    if (ET.inputdata.filt$va[1] && !getVprp)
      stop('Calculation of ET for the given function requires extractions of tmax (i.e. set getVprp=T)')
    if (ET.inputdata.filt$Precip[1] && !getPrecip)
      stop('Calculation of ET for the given function requires extractions of precip (i.e. set getPrecip=T)')
  }

  # Checking the DEM
  if (getET) {
    if (is.character(DEM)) {
      if (!file.exists(DEM))
        stop(paste('The following input file for the DEM could not be found:',DEM))

      # Read in DEM
      DEM = sp::read.asciigrid(DEM,  proj4string = sp::CRS("+proj=longlat +ellps=GRS80"))

      # Convert to Raster object
      DEM = raster::raster(DEM)
    } else if (!methods::is(DEM,"RasterLayer")) {
      stop('The input for the DEM is not a raster layer object.')
    }
  }

  # Open file with polygons
  if (is.character(catchments)) {
    if (!file.exists(catchments))
      stop(paste('The following input file for the catchments could not be found:',catchments))

    # Read in polygons
    catchments = maptools::readShapeSpatial(catchments, force_ring=TRUE)

  } else if (!methods::is(catchments,"SpatialPolygonsDataFrame") && !methods::is(catchments,"SpatialPointsDataFrame")) {
    stop('The input for "catchments" must be a file name to a shape file or a SpatialPolygonsDataFrame or a SpatialPointsDataFrame object.')
  }

  # Check projection of catchments
  if (is.na(sp::proj4string(catchments))) {
    message('WARNING: The projection string of the catchment boundaries is NA. Setting to +proj=longlat +ellps=GRS80.')
    sp::proj4string(catchments) = '+proj=longlat +ellps=GRS80'
  }
  if (sp::proj4string(catchments) != '+proj=longlat +ellps=GRS80') {
    message('WARNING: The projection string of the catchment boundaries does not appear to be +proj=longlat +ellps=GRS80. Attempting to transform coordinates...')
    catchments = sp::spTransform(catchments,sp::CRS('+proj=longlat +ellps=GRS80'))
  }

  # Check if catchments are points or a polygon.
  isCatchmentsPolygon=TRUE;
  if (methods::is(catchments,"SpatialPointsDataFrame")) {
    isCatchmentsPolygon=FALSE;
  }

  # Get netCDF geometry
  timePoints <- ncdf4::ncvar_get(awap, "time")
  nTimePoints = length(timePoints);
  Long <- ncdf4::ncvar_get(awap, "Long")
  Lat <- ncdf4::ncvar_get(awap, "Lat")
  if (getSolarrad) {
    timePoints_solar <- ncdf4::ncvar_get(awap_solar, "time")
    nTimePoints_solar = length(timePoints_solar);
    Long_solar <- ncdf4::ncvar_get(awap_solar, "Long")
    Lat_solar <- ncdf4::ncvar_get(awap_solar, "Lat")
  }

  # Convert the ncdf time to R time.
  # This is achieved by spliting the time units string into fields.
  # Adapted form http://geog.uoregon.edu/bartlein/courses/geog607/Rmd/netCDF_01.htm.
  tunits <- ncdf4::ncatt_get(awap, "time", "units")
  tustr <- strsplit(tunits$value, " ")
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth = as.integer(unlist(tdstr)[2])
  tday = as.integer(unlist(tdstr)[3])
  tyear = as.integer(unlist(tdstr)[1])
  timePoints_R = as.Date(chron::chron(timePoints, origin = c(month=tmonth, day=tday, year=tyear),format=c(dates = "y-m-d")));
  if (getSolarrad) {
    tunits <- ncdf4::ncatt_get(awap_solar, "time", "units")
    tustr <- strsplit(tunits$value, " ")
    tdstr <- strsplit(unlist(tustr)[3], "-")
    tmonth = as.integer(unlist(tdstr)[2])
    tday = as.integer(unlist(tdstr)[3])
    tyear = as.integer(unlist(tdstr)[1])
    timePoints_solar_R = as.Date(chron::chron(timePoints_solar, origin = c(tmonth, tday, tyear),format=c(dates = "y-m-d")));
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

  if (extractTo < extractFrom) {
    stop('extractTo date is less than the extractFrom date.')
  }

  # Recalculate the time points to extract.
  timepoints2Extract = seq( as.Date(extractFrom,'%Y-%m-%d'), by="day", to=as.Date(extractTo,'%Y-%m-%d'))

  message(paste('    Data will be extracted from ',format.Date(extractFrom,'%Y-%m-%d'),' to ', format.Date(extractTo,'%Y-%m-%d'),' at ',length(catchments),' catchments '));
  message('Starting data extraction:')

  # Get one netCDF layer.
  precipGrd = raster::raster(ncdfFilename, band=nTimePoints, varname='precip',lvar=3)
  if (getSolarrad) {
    solarGrd = raster::raster(ncdfSolarFilename, band=nTimePoints_solar, varname='solarrad',lvar=3)
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
        raster::removeTmpFiles(h=0)
      }
      w = raster::rasterize(catchments[i,], precipGrd,getCover=T)

      # Extract the mask values (i.e. fraction of each grid cell within the polygon.
      w2 = raster::getValues(w);
      filt = w2>0
      wLongLat = sp::coordinates(w)[filt,]
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
    longLat.all = cbind(as.numeric(sp::coordinates(catchments)[,1]),as.numeric(sp::coordinates(catchments)[,2]))
    catchment.lookup = cbind(seq(1,length(catchments),by=1),seq(1,length(catchments),by=1));
  }

  raster::removeTmpFiles(h=0)

  if (getSolarrad && getET) {
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

    if (j%%(365*5) ==0 ) {
      message(paste('    ... Extracting data for time point ', j,' of ',length(timepoints2Extract)));
    }
    if (getPrecip)
      precip[j,1:length(w.all)]  <- raster::extract(raster::raster(ncdfFilename, band=ind, varname='precip',lvar=3, method=interpMethod), longLat.all)
    if (getTmin)
      tmin[j,1:length(w.all)]  <- raster::extract(raster::raster(ncdfFilename, band=ind, varname='tmin',lvar=3, method=interpMethod), longLat.all)
    if (getTmax)
      tmax[j,1:length(w.all)]  <- raster::extract(raster::raster(ncdfFilename, band=ind, varname='tmax',lvar=3, method=interpMethod), longLat.all)
    if (getVprp)
      vprp[j,1:length(w.all)]  <- raster::extract(raster::raster(ncdfFilename, band=ind, varname='vprp',lvar=3, method=interpMethod), longLat.all)

    # Get date of extracted grid.
    extractDate = raster::getZ(raster::raster(ncdfFilename, band=ind,lvar=3));
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
        solarrad[j,1:length(w.all)]  <- raster::extract(raster::raster(ncdfSolarFilename, band=ind, varname='solarrad',lvar=3, method=interpMethod), longLat.all)

        # Get date of extracted grid.
        extractDate = raster::getZ(raster::raster(ncdfSolarFilename, band=ind,lvar=3));
        extractYear_solar[j] = as.integer(format(as.Date(extractDate,'%Y-%m-%d'),'%Y'));
        extractMonth_solar[j] = as.integer(format(as.Date(extractDate,'%Y-%m-%d'),'%m'));
        extractDay_solar[j] = as.integer(format(as.Date(extractDate,'%Y-%m-%d'),'%d'));
      }
    }
  }

  # Close netCDF connection
  if (!getPrecip || !getTmin || !getTmax || !getVprp)
    ncdf4::nc_close(awap)
  if (getSolarrad)
    ncdf4::nc_close(awap_solar)


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
      } else if (length(stats::na.omit(solarrad[ind,]))>0) {
        solarrad_avg[j,] = apply(stats::na.omit(solarrad[ind,]),2,mean)
      }
    }

    # Assign the daily average solar radiation to each day prior to 1 Jan 1990
    for (j in 1:length(timepoints2Extract)) {
      if (timepoints2Extract[j]<as.Date('1990-1-1','%Y-%m-%d')) {
        ind = monthdayUnique==monthdayAll[j]

        # Only assign averages if available eg Feb 29 may not be available in average record if the post
        # 1990 extraction period does not include a leap year.
        if (sum(ind)>0) {
          solarrad_interp[j,] = solarrad_avg[ind,]
        }
      }
    }

    # Linearly interpolate time points without a solar radiation value.
    message('... Linearly interpolating gaps in daily solar.')
    for (j in 1:length(w.all)) {
      filt = is.na(solarrad_interp[,j])
      x = 1:length(timepoints2Extract);
      xpred = x[filt];

      # Interpolate if any NAs
      if (length(xpred)>0) {
        x = x[!filt]
        y = solarrad_interp[!filt,j]

        # Interpolate if at least 2 non-NA obs.
        if (length(y)>1) {
          ypred=stats::approx(x,y,xpred,method='linear', rule=2)
          solarrad_interp[filt,j] = ypred$y
        }
      }

    }
  }

  # Loop though each catchment and calculate the catchment averager and variance.
  if (getET) {
    if (!exists('ET.constants') || !is.list(ET.constants))
    stop('ET.constants is empty or not a list.')
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
    if (getET) {
      # Initialise outputa
      mortAPET = matrix(NA,length(timepoints2Extract),length(ind))

      # Loop through each grid cell
      k=0;
      for (j in ind) {

        k=k+1
        if (k==1 || k%%25 ==0 ) {
          message(paste('           ... Calculating PET for grid cell', k,' of ',length(w)));
        }

        # Check lat, Elev and precip are finite scalers.
        if (any(!is.finite(precip[,j])) || !is.finite(DEMpoints[j]) || !is.finite(longLat.all[j,2])) {
          message(paste('WARNING: Non-finite input values detected for catchment',i,' at grid cell',j))
          message(paste('   Elevation value:' ,DEMpoints[j]))
          message(paste('   Latitude value:' ,longLat.all[j,2]))
          ind.nonfinite = which(!is.finite(precip[,j]))
          if (length(ind.nonfinite)>0)
            message(paste('   Precipiation nonfinite value (first):' ,precip[ind.nonfinite[1],j]))
          mortAPET[,k] = NA;
          next
        }

        # Build data from of daily climate data
        dataRAW = data.frame(Year =  as.integer(format.Date(timepoints2Extract,"%Y")), Month= as.integer(format.Date(timepoints2Extract,"%m")), Day= as.integer(format.Date(timepoints2Extract,"%d")),
                             Tmin=tmin[,j], Tmax=tmax[,j], Rs=solarrad_interp[,j], va=vprp[,j]/10.0, Precip=precip[,j])

        # # Filter out columns not needed.
        # if (!ET.inputdata.filt$Tmin)
        #   dataRAW$Tmin = NULL
        # if (!ET.inputdata.filt$Tmax)
        #   dataRAW$Tmax = NULL
        # if (!ET.inputdata.filt$Rs)
        #   dataRAW$Rs = NULL
        # if (!ET.inputdata.filt$va)
        #   dataRAW$va = NULL
        # if (!ET.inputdata.filt$Precip)
        #   dataRAW$Precip = NULL

        # Convert to required format for ET package
        dataPP=Evapotranspiration::ReadInputs(ET.var.names ,dataRAW,constants=NA,stopmissing = c(99,99,99),
                          interp_missing_days=ET.interp_missing_days, interp_missing_entries=ET.interp_missing_entries, interp_abnormal=ET.interp_abnormal,
                          missing_method='DoY average', abnormal_method='DoY average', message = "no")

        # Update constants for the site
        constants$Elev = DEMpoints[j]
        constants$lat = longLat.all[j,2]
        constants$lat_rad = longLat.all[j,2]/180.0*pi

        # Call  ET package
        if (ET.function=='ET.Abtew') {
          results <- Evapotranspiration::ET.Abtew(dataPP, ET.constants, ts=ET.timestep,solar="data",AdditionalStats='no', message='no');
        } else if(ET.function=='ET.HargreavesSamani') {
          results <- Evapotranspiration::ET.HargreavesSamani(dataPP, ET.constants, ts=ET.timestep,AdditionalStats='no', message='no');
        } else if(ET.function=='ET.JensenHaise') {
          results <- Evapotranspiration::ET.JensenHaise(dataPP, ET.constants, ts=ET.timestep,solar="data",AdditionalStats='no', message='no');
        } else if(ET.function=='ET.Makkink') {
          results <- Evapotranspiration::ET.Makkink(dataPP, ET.constants, ts=ET.timestep,solar="data",AdditionalStats='no', message='no');
        } else if(ET.function=='ET.McGuinnessBordne') {
          results <- Evapotranspiration::ET.McGuinnessBordne(dataPP, ET.constants, ts=ET.timestep,AdditionalStats='no', message='no');
        } else if(ET.function=='ET.MortonCRAE') {
          results <- Evapotranspiration::ET.MortonCRAE(dataPP, ET.constants,est=ET.Mortons.est, ts=ET.timestep,solar="data",Tdew=FALSE, AdditionalStats='no', message='no');
        } else if(ET.function=='ET.MortonCRWE') {
          results <- Evapotranspiration::ET.MortonCRWE(dataPP, ET.constants,est=ET.Mortons.est, ts=ET.timestep,solar="data",Tdew=FALSE, AdditionalStats='no', message='no');
        } else if(ET.function=='ET.Turc') {
          results <- Evapotranspiration::ET.Turc(dataPP, ET.constants, ts=ET.timestep,solar="data",humid=F, AdditionalStats='no', message='no');
        }


        # Interpolate monthly or annual data
        if (ET.timestep=='monthly' || ET.timestep=='annual') {
          # Get the last day of each month
          last.day.month = zoo::as.Date(zoo::as.yearmon(stats::time(results$ET.Monthly)), frac = 1)

          # Get days per month
          days.per.month = as.integer(format.Date(zoo::as.Date(zoo::as.yearmon(stats::time(results$ET.Monthly)), frac = 1),'%d'))

          # Set the first month to the start date for extraction.
          start.day.month = as.numeric(format(timepoints2Extract,"%d"))[1]
          days.per.month[1] = days.per.month[1] - start.day.month + 1

          # Set the last month to the end date for extraction.
          end.day.month = as.numeric(format(timepoints2Extract,"%d"))[length(timepoints2Extract)]
          days.per.month[length(days.per.month)] = end.day.month

          # Calculate average ET per day of each month
          monthly.ET.as.daily = zoo::zoo( as.numeric(results$ET.Monthly/days.per.month), last.day.month)

          # Spline interpolate Monthly average ET
          timepoints2Extract.as.zoo = zoo::zoo(NA,timepoints2Extract);
          mortAPET.tmp = zoo::na.approx(merge(monthly.ET.as.daily, dates=timepoints2Extract.as.zoo)[, 1], rule=2)
          filt = stats::time(mortAPET.tmp)>=stats::start(timepoints2Extract.as.zoo) & stats::time(mortAPET.tmp)<=stats::end(timepoints2Extract.as.zoo)
          mortAPET.tmp = pmax(0.0, as.numeric( mortAPET.tmp));
          mortAPET.tmp = mortAPET.tmp[filt]
          mortAPET[,k] = mortAPET.tmp;
        } else {
          mortAPET[,k] = results$ET.Daily
        }

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
        catchmentVarTmp$precip_mm = apply(precip[,ind],1,spatial.function.name,na.rm=TRUE);
      } else {
        catchmentVarTmp$precip_mm = NA;
      }
    }
    if (getTmin) {
      catchmentAvgTmp$Tmin = apply(t(t(tmin[,ind]) * w),1,sum,na.rm=TRUE);
      if (calcVariance) {
        catchmentVarTmp$Tmin = apply(tmin[,ind],1,spatial.function.name,na.rm=TRUE);
      } else {
        catchmentVarTmp$Tmin = NA;
      }
    }
    if (getTmax) {
      catchmentAvgTmp$Tmax = apply(t(t(tmax[,ind]) * w),1,sum,na.rm=TRUE);
      if (calcVariance) {
        catchmentVarTmp$Tmax = apply(tmax[,ind],1,spatial.function.name,na.rm=TRUE);
      } else {
        catchmentVarTmp$Tmax = NA;
      }
    }
    if (getVprp) {
      catchmentAvgTmp$vprp = apply(t(t(vprp[,ind]) * w),1,sum,na.rm=TRUE);
      if (calcVariance) {
        catchmentVarTmp$vprp = apply(vprp[,ind],1,spatial.function.name,na.rm=TRUE);
      } else {
        catchmentVarTmp$vprp = NA;
      }
    }
    if (getSolarrad) {
      catchmentAvgTmp$solarrad = apply(t(t(solarrad[,ind]) * w),1,sum,na.rm=TRUE);
      catchmentAvgTmp$solarrad_interp = apply(t(t(solarrad_interp[,ind]) * w),1,sum,na.rm=TRUE);
      if (calcVariance) {
        catchmentVarTmp$solarrad = apply(solarrad[,ind],1,spatial.function.name,na.rm=TRUE);
        catchmentVarTmp$solarrad_interp = apply(solarrad_interp[,ind],1,spatial.function.name,na.rm=TRUE);
      } else {
        catchmentVarTmp$solarrad = NA;
        catchmentVarTmp$solarrad_interp = NA;
      }

    }
    if (getET) {
      catchmentAvgTmp$ET_mm = apply(t(t(mortAPET) * w),1,sum,na.rm=TRUE);
      if (calcVariance) {
        catchmentVarTmp$ET_mm = apply(mortAPET,1,spatial.function.name,na.rm=TRUE);
      } else {
        catchmentVarTmp$ET_mm = NA;
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
    catchmentAvg = list(catchmentAvg, catchmentVar)
    names(catchmentAvg) = c('catchmentAvg', paste('catchment',spatial.function.name,sep=''))
  }

  message('Data extraction FINISHED..')
  return(catchmentAvg)
}
