#'
#' \code{extractCatchmentData} extracts catchment average climate data from netCDF files containing Australian climate data.
#'
#' @description
#' extractCatchmentData extracts the AWAP climate data for each point or polygon. For the latter, either the daily spatial mean and variance (or user defined function) of
#' each climate metric is calculated or the spatial data is returned.
#'
#' @details
#' Daily data is extracted and can be aggregated to a weekly, monthly, quarterly, annual or a user-defined timestep using a user-defined funcion
#' (e.g. sum, mean, min, max as defined by \code{temporal.function.name}). The temporally aggreated data at each grid cell is then used to derive the spatial
#' mean or the spatial variance (or any other function as defined by \code{spatial.function.name}).
#'
#' The calculation of the spatial mean uses the fraction of each AWAP grid cell within the catchment polygon.
#' The variance calculation (or user defined function) does not use the fraction of the grid cell and returns NA if there are <2 grid cells in the catchment boundary.
#' Prior to the spatial aggregation, evapotranspiration (ET) can also calculated; after which, say, the mean and
#' variance PET can be calculated.
#'
#' The data extraction will by default be undertaken from 1/1/1900 to yesterday, even if the netCDF grids were only
#' built for a subset of this time period. If the latter situation applies, it is recommended that the extraction start
#' and end dates are input by the user.
#'
#' The ET can be calculated using one of eight methods at a user defined calculation time-step; that is the \code{ET.timestep} defines the
#' time step at which the estimates are serived and differs from the output timestep as defined by \code{temporal.function.name}). When \code{ET.timestep} is monthly or annual then
#' the ET estimate is linearly interpolated to a daily time step (using zoo:na.spline()) and then constrained to >=0. In calculating ET, the input data
#' is pre-processed using Evapotranspiration::ReadInputs() such that missing days, missing enteries and abnormal values are interpolated
#' (by default) with the former two interpolated using the "DoY average", i.e. replacement with same day-of-the-year average. Additionally, when AWAP solar
#' radiation is required for the ET function, data is only available from 1/1/1990. To derive ET values <1990, the average solar radiation for each day of the year from
#' 1/1/990 to "extractTo" is derived (i.e. 365 values) and then applied to each day prior to 1990. Importantly, in this situation the estimates of ET <1990
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
#' already imported using readShapeSpatial(). Either way the shape file must be in long/lat (i.e. not projected), use the ellipsoid GRS 80, and the first column should be a unique ID.
#' @param temporal.timestep character string for the time step of the output data. The options are \code{daily}, \code{weekly}, \code{monthly}, \code{quarterly},
#' \code{annual}  or a user-defined index for, say, water-years (see \code{xts::period.apply}). The default is \code{daily}.
#' @param temporal.function.name character string for the function name applied to aggregate the daily data to \code{temporal.timestep}.
#' Note, NA values are not removed from the aggregation calculation. If this is required then consider writing your own function. The default is \code{mean}.
#' @param spatial.function.name character string for the function name applied to estimate the daily spatial spread in each variable. If \code{NA} or \code{""} and \code{catchments} is a polygon, then
#' the spatial data is returned. The default is \code{var}.
#' @param interpMethod character string defining the method for interpolating the gridded data (see \code{raster::extract}). The options are: \code{'simple'}, \code{'bilinear'} and \code{''}. The default
#' is \code{''}. This will set the interpolation to \code{'simple'} when \code{catchments} is a polygon(s) and to \code{'bilinear'} when \code{catchments} are points.
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
#' When \code{catchments} are polygons and \code{spatial.function.name} is not \code{NA} or \code{""}, then the returned variable is a list variable containing two data.frames. The first is the areal aggreggated climate
#' metrics named \code{catchmentTemporal.} with a suffix as defined by \code{temporal.function.name}). The second is the measure of spatial variability
#' named \code{catchmentSpatial.} with a suffix as defined by \code{spatial.function.name}).
#'
#' When \code{catchments} are polygons and \code{spatial.function.name} does equal \code{NA} or \code{""}, then the returned variable is a \code{sp::SpatialPixelsDataFrame} where the first colum is the catchment IDs
#' and the latter columns are the results for each variable at each time point as defined by \code{temporal.timestep}.
#'
#' When \code{catchments} are points, the returned variable is a data.frame containing daily climate data at each point.
#'
#' @seealso
#' \code{\link{makeNetCDF_file}} for building the NetCDF files of daily climate data.
#'
#' @examples
#' # The example shows how to extract and save data.
#' # For an additional example see \url{https://github.com/peterson-tim-j/AWAPer/blob/master/README.md}
#' #---------------------------------------
#' library(sp)
#'
#' # Set dates for building netCDFs and extracting data.
#' # Note, to reduce runtime this is done only a fortnight (14 days).
#' startDate = as.Date("2000-01-01","%Y-%m-%d")
#' endDate = as.Date("2000-01-14","%Y-%m-%d")
#'
#' # Set names for netCDF files.
#' ncdfFilename = tempfile(fileext='.nc')
#' ncdfSolarFilename = tempfile(fileext='.nc')
#'
#' # Build netCDF grids and over a defined time period.
#' # Only precip data is to be added to the netCDF files.
#' # This is because the URLs for the other variables are set to zero.
#' \donttest{
#' file.names = makeNetCDF_file(ncdfFilename=ncdfFilename,
#'              ncdfSolarFilename=ncdfSolarFilename,
#'              updateFrom=startDate, updateTo=endDate,
#'              urlTmin=NA, urlTmax=NA, urlVprp=NA, urlSolarrad=NA)
#'
#' # Load example catchment boundaries and remove all but the first.
#' # Note, this is done only to speed up the example runtime.
#' data("catchments")
#' catchments = catchments[1,]
#'
#' # Extract daily precip. data (not Tmin, Tmax, VPD, ET).
#' # Note, the input "catchments" can also be a file to a ESRI shape file.
#' climateData = extractCatchmentData(ncdfFilename=file.names$ncdfFilename,
#'               ncdfSolarFilename=file.names$ncdfSolarFilename,
#'               extractFrom=startDate, extractTo=endDate,
#'               getTmin = FALSE, getTmax = FALSE, getVprp = FALSE,
#'               getSolarrad = FALSE, getET = FALSE,
#'               catchments=catchments,
#'               temporal.timestep = 'daily')
#'
#' # Extract the daily catchment average data.
#' climateDataAvg = climateData$catchmentTemporal.mean
#'
#' # Extract the daily catchment variance data.
#' climateDataVar = climateData$catchmentSpatial.var
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
  temporal.timestep = 'daily',
  temporal.function.name = 'mean',
  spatial.function.name = 'var',
  interpMethod = '',
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

  # Check time step
  temporal.timestep.options = c('daily','weekly','monthly','quarterly','annual', 'period')
  if (!is.character(temporal.timestep)) {
    if (!is.integer(temporal.timestep))
      stop('temporal.timestep must be a character string or an integer vector.')

    if (any(temporal.timestep<0 || temporal.timestep > length(timepoints2Extract)))
      stop(paste('temporal.timestep must be a character string or an integer vector with values >0 and <= the maximum days of data extracted (ie',
                 length(timepoints2Extract),')'))

    temporal.timestep.index = temporal.timestep
    temporal.timestep = 'period'

  } else if (!any(which(temporal.timestep.options == temporal.timestep))) {
    stop('temporal.timestep must be daily, weekly, monthly, quarterly or annual or an integer vector.')
  }

  # Check temporal analysis funcion is valid.
  data.junk = t(as.matrix(stats::runif(100, 0.0, 1.0)*100))
  result = tryCatch({
    apply(data.junk, 1,FUN=temporal.function.name)
  }, warning = function(w) {
    message(paste("Warning temporal.function.name produced the following",w))
  }, error = function(e) {
    stop(paste('temporal.function.name produced an error when applied using test data:',e))
  }, finally = {
    rm(data.junk)
  }
  )

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
    if (ET.inputdata.filt$Rs[1]) {
      if (!file.exists(ncdfSolarFilename)) {
        stop('The input ncdfSolarFilename is required for the selected ET function.')
      }
      if (!getSolarrad) {
        stop('Calculation of ET for the given function requires extractions of solar radiation (i.e. set getSolarrad=T)')
      }
    }
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

  # Check the interpolation method.
  if (interpMethod!='' && interpMethod!='simple' && interpMethod!='bilinear') {
    stop('The input for "interpMethod" must either "simple" or "bilinear".')
  }
  if (interpMethod=='') {
    if (isCatchmentsPolygon) {
      interpMethod='simple';
    } else {
      interpMethod='bilinear';
    }
  }

  # Check if the spatial data should be returned or analysed.
  do.spatial.analysis=T
  if (isCatchmentsPolygon && is.na(spatial.function.name) || (is.character(spatial.function.name) && spatial.function.name=='')) {
    do.spatial.analysis=F
  }

  if (do.spatial.analysis) {
    # Check spatial analysis funcion is valid.
    data.junk = t(as.matrix(stats::runif(100, 0.0, 1.0)*100))
    result = tryCatch({
      apply(data.junk, 1,FUN=spatial.function.name)
    }, warning = function(w) {
      message(paste("Warning spatial.function.name produced the following",w))
    }, error = function(e) {
      stop(paste('spatial.function.name produced an error when applied using test data:',e))
    }, finally = {
      rm(data.junk)
    }
    )
  }

  # Open solar rad netCDF file
  if (getSolarrad)
    awap_solar <- ncdf4::nc_open(ncdfSolarFilename)

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
  message('... Building catchment weights');
  if (isCatchmentsPolygon) {

    w.all = c();
    longLat.all = matrix(NA,0,2)
    catchment.lookup = matrix(NA,length(catchments),2);

    for (i in 1:length(catchments)) {
      if (i%%10 ==0 ) {
        message(paste('   ... Building weights for catchment ', i,' of ',length(catchments)));
        raster::removeTmpFiles(h=0)
      }

      # Extract the weights for grid cells within the catchments polygon.
      # Note, the AWAP raster grid is cropped to the extent of the catchment polygon.
      # This was undertaken to improve the run time performance but more importantly to overcome an error
      # thrown by raster::rasterize when the raster is large (see https://github.com/rspatial/raster/issues/192).
      # This solution should work when the catchments polygon is not very large (e.g. not a reasonable fraction of the
      # Australian land mass)
      w = raster::rasterize(x=catchments[i,], y=raster::crop(precipGrd, catchments[i,], snap='out'), fun='last',getCover=T)

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

  message('... Starting to extract data across all catchments:')
  for (j in 1:length(timepoints2Extract)){

    # Find index to the date to update within the net CDF grid
    ind = as.integer(difftime(timepoints2Extract[j], as.Date("1900-1-1",'%Y-%m-%d'),units = "days" ))+1

    if (j%%(365*5) ==0 ) {
      message(paste('    ... Extracting data for time point ', j,' of ',length(timepoints2Extract)));
    }
    if (getPrecip)
      precip[j,1:length(w.all)]  <- raster::extract(raster::raster(ncdfFilename, band=ind, varname='precip',lvar=3), longLat.all, method=interpMethod)
    if (getTmin)
      tmin[j,1:length(w.all)]  <- raster::extract(raster::raster(ncdfFilename, band=ind, varname='tmin',lvar=3), longLat.all, method=interpMethod)
    if (getTmax)
      tmax[j,1:length(w.all)]  <- raster::extract(raster::raster(ncdfFilename, band=ind, varname='tmax',lvar=3), longLat.all, method=interpMethod)
    if (getVprp)
      vprp[j,1:length(w.all)]  <- raster::extract(raster::raster(ncdfFilename, band=ind, varname='vprp',lvar=3), longLat.all, method=interpMethod)

    if (getPrecip) {
      extractDate = raster::getZ(raster::raster(ncdfFilename, band=ind,lvar=3,varname='precip'))
    } else if (getTmin) {
      extractDate = raster::getZ(raster::raster(ncdfFilename, band=ind,lvar=3,varname='tmin'))
    } else if (getTmax) {
      extractDate = raster::getZ(raster::raster(ncdfFilename, band=ind,lvar=3,varname='tmax'))
    } else if (getVprp) {
      extractDate = raster::getZ(raster::raster(ncdfFilename, band=ind,lvar=3,varname='vprp'))
    }

    # Get date of extracted grid.
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
        solarrad[j,1:length(w.all)]  <- raster::extract(raster::raster(ncdfSolarFilename, band=ind, varname='solarrad',lvar=3), longLat.all, method=interpMethod)

        # Get date of extracted grid.
        extractDate = raster::getZ(raster::raster(ncdfSolarFilename, band=ind,lvar=3, varname='solarrad'));
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
      } else {
        if (ncol(solarrad)==1) {
          solarrad_avg[j,] = mean(solarrad[ind,],na.rm=T)
        } else {
          solarrad_avg[j,] = apply(stats::na.omit(solarrad[ind,]),2,mean)
        }

      }
    }

    # Assign the daily average solar radiation to each day prior to 1 Jan 1990
    for (j in 1:length(timepoints2Extract)) {
      if (timepoints2Extract[j]<as.Date('1990-1-1','%Y-%m-%d')) {
        ind = which(monthdayUnique==monthdayAll[j])
        if (length(ind)==1) {
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

        # Convert to required format for ET package.
        # Note, missing_method changed from NULL to "DoY average" because individual grid cells can be NA. eg
        dataPP=Evapotranspiration::ReadInputs(ET.var.names ,dataRAW,constants=NA,stopmissing = c(99,99,99),
                                              interp_missing_days=ET.interp_missing_days, interp_missing_entries=ET.interp_missing_entries, interp_abnormal=ET.interp_abnormal,
                                              missing_method="DoY average", abnormal_method='DoY average', message = "no")

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

    # Calculate catchment stats and add data to the data.frame for all catchments
    #----------------------------------------------------------------------------------------
    extractDate = ISOdate(year=extractYear,month=extractMonth,day=extractDay)


    # Do the temporal aggregation
    #-----------------------------
    precip.xts = NA
    tmin.xts = NA
    tmax.xts = NA
    vprp.xts = NA
    solarrad.xts = NA
    solarrad_interp.xts = NA
    ET.xts = NA
    if (getPrecip) {
      precip.xts = xts::as.xts(precip[,ind], order.by=extractDate)
      precip.xts <-
        switch(
          which(temporal.timestep.options == temporal.timestep),
          precip.xts, # daily timestep - do nothing
          xts::apply.weekly(precip.xts, apply, 2, temporal.function.name),  # weekly timestep
          xts::apply.monthly(precip.xts, apply, 2, temporal.function.name), # monthly timestep
          xts::apply.quarterly(precip.xts, apply, 2, temporal.function.name), # quarterly timestep
          xts::apply.yearly(precip.xts, apply, 2, temporal.function.name), # annual timestep
          xts::period.apply(precip.xts, INDEX=temporal.timestep.index, apply, 2, temporal.function.name), # user defined period
        )
    }

    if (getTmin) {
      tmin.xts = xts::as.xts(tmin[,ind], order.by=extractDate)
      tmin.xts <-
        switch(
          which(temporal.timestep.options == temporal.timestep),
          tmin.xts, # daily timestep - do nothing
          xts::apply.weekly(tmin.xts, apply, 2, temporal.function.name),  # weekly timestep
          xts::apply.monthly(tmin.xts, apply, 2, temporal.function.name), # monthly timestep
          xts::apply.quarterly(tmin.xts, apply, 2, temporal.function.name), # quarterly timestep
          xts::apply.yearly(tmin.xts, apply, 2, temporal.function.name), # annual timestep
          xts::period.apply(tmin.xts, INDEX=temporal.timestep.index, apply, 2, temporal.function.name), # user defined period
        )
    }

    if (getTmax) {
      tmax.xts = xts::as.xts(tmax[,ind], order.by=extractDate)
      tmax.xts <-
        switch(
          which(temporal.timestep.options == temporal.timestep),
          tmax.xts, # daily timestep - do nothing
          xts::apply.weekly(tmax.xts, apply, 2, temporal.function.name),  # weekly timestep
          xts::apply.monthly(tmax.xts, apply, 2, temporal.function.name), # monthly timestep
          xts::apply.quarterly(tmax.xts, apply, 2, temporal.function.name), # quarterly timestep
          xts::apply.yearly(tmax.xts, apply, 2, temporal.function.name), # annual timestep
          xts::period.apply(tmax.xts, INDEX=temporal.timestep.index, apply, 2, temporal.function.name), # user defined period
        )
    }

    if (getVprp) {
      vprp.xts = xts::as.xts(vprp[,ind], order.by=extractDate)
      vprp.xts <-
        switch(
          which(temporal.timestep.options == temporal.timestep),
          vprp.xts, # daily timestep - do nothing
          xts::apply.weekly(vprp.xts, apply, 2, temporal.function.name),  # weekly timestep
          xts::apply.monthly(vprp.xts, apply, 2, temporal.function.name), # monthly timestep
          xts::apply.quarterly(vprp.xts, apply, 2, temporal.function.name), # quarterly timestep
          xts::apply.yearly(vprp.xts, apply, 2, temporal.function.name), # annual timestep
          xts::period.apply(vprp.xts, INDEX=temporal.timestep.index, apply, 2, temporal.function.name), # user defined period
        )
    }

    if (getSolarrad) {
      solarrad.xts = xts::as.xts(solarrad[,ind], order.by=extractDate)
      solarrad.xts <-
        switch(
          which(temporal.timestep.options == temporal.timestep),
          solarrad.xts, # daily timestep - do nothing
          xts::apply.weekly(solarrad.xts, apply, 2, temporal.function.name),  # weekly timestep
          xts::apply.monthly(solarrad.xts, apply, 2, temporal.function.name), # monthly timestep
          xts::apply.quarterly(solarrad.xts, apply, 2, temporal.function.name), # quarterly timestep
          xts::apply.yearly(solarrad.xts, apply, 2, temporal.function.name), # annual timestep
          xts::period.apply(solarrad.xts, INDEX=temporal.timestep.index, apply, 2, temporal.function.name), # user defined period
        )

      solarrad_interp.xts = xts::as.xts(solarrad_interp[,ind], order.by=extractDate)
      solarrad_interp.xts <-
        switch(
          which(temporal.timestep.options == temporal.timestep),
          solarrad_interp.xts, # daily timestep - do nothing
          xts::apply.weekly(solarrad_interp.xts, apply, 2, temporal.function.name),  # weekly timestep
          xts::apply.monthly(solarrad_interp.xts, apply, 2, temporal.function.name), # monthly timestep
          xts::apply.quarterly(solarrad_interp.xts, apply, 2, temporal.function.name), # quarterly timestep
          xts::apply.yearly(solarrad_interp.xts, apply, 2, temporal.function.name), # annual timestep
          xts::period.apply(solarrad_interp.xts, INDEX=temporal.timestep.index, apply, 2, temporal.function.name), # user defined period
        )
    }

    if (getET) {
      ET.xts = xts::as.xts(mortAPET, order.by=extractDate)
      ET.xts <-
        switch(
          which(temporal.timestep.options == temporal.timestep),
          ET.xts, # daily timestep - do nothing
          xts::apply.weekly(ET.xts, apply, 2, temporal.function.name),  # weekly timestep
          xts::apply.monthly(ET.xts, apply, 2, temporal.function.name), # monthly timestep
          xts::apply.quarterly(ET.xts, apply, 2, temporal.function.name), # quarterly timestep
          xts::apply.yearly(ET.xts, apply, 2, temporal.function.name), # annual timestep
          xts::period.apply(ET.xts, INDEX=temporal.timestep.index, apply, 2, temporal.function.name), # user defined period
        )
    }
    #-----------------------------

    # Determine the number of time steps (for creation of the output)
    if (xts::is.xts(precip.xts)) {
      timesteps=zoo::index(zoo::as.zoo(precip.xts))
    } else if (xts::is.xts(tmin.xts)) {
      timesteps=zoo::index(zoo::as.zoo(tmin.xts))
    } else if (xts::is.xts(tmax.xts)) {
      timesteps=zoo::index(zoo::as.zoo(tmax.xts))
    } else if (xts::is.xts(vprp.xts)) {
      timesteps=zoo::index(zoo::as.zoo(vprp.xts))
    } else if (xts::is.xts(solarrad.xts)) {
      timesteps=zoo::index(zoo::as.zoo(solarrad.xts))
    } else if (xts::is.xts(ET.xts)) {
      timesteps=zoo::index(zoo::as.zoo(ET.xts))
    }
    nTimeSteps = length(timesteps)

    # Do spatial aggregation if a spatial function is provided. If spatial aggregation
    # is not required, then build up a data.frame (each time step and data type as one column)
    # for latter conversion to a spatial object.
    #-----------------------------
    if (do.spatial.analysis) {
      # Create data frame for final results
      catchmentAvgTmp = data.frame(CatchmentID=catchments[i,1],year=as.integer(format(timesteps,"%Y")) ,month=as.integer(format(timesteps,"%m")),day=as.integer(format(timesteps,"%d")),row.names = NULL);
      catchmentVarTmp = catchmentAvgTmp

      # Check if there are enough grid cells to calculate the variance.
      calcVariance =  F;
      if (length(ind)>1)
        calcVariance =  T;

      # Do spatial aggregation using grid cell weights
      #-----------------------------
      if (getPrecip) {
        catchmentAvgTmp$precip_mm = apply(t(t(precip.xts) * w),1,sum,na.rm=TRUE)
        if (calcVariance) {
          catchmentVarTmp$precip_mm = apply(precip.xts,1,spatial.function.name,na.rm=TRUE);
        } else {
          catchmentVarTmp$precip_mm = NA;
        }
      }
      if (getTmin) {
        catchmentAvgTmp$Tmin = apply(t(t(tmin.xts) * w),1,sum,na.rm=TRUE);
        if (calcVariance) {
          catchmentVarTmp$Tmin = apply(tmin.xts,1,spatial.function.name,na.rm=TRUE);
        } else {
          catchmentVarTmp$Tmin = NA;
        }
      }
      if (getTmax) {
        catchmentAvgTmp$Tmax = apply(t(t(tmax.xts) * w),1,sum,na.rm=TRUE);
        if (calcVariance) {
          catchmentVarTmp$Tmax = apply(tmax.xts,1,spatial.function.name,na.rm=TRUE);
        } else {
          catchmentVarTmp$Tmax = NA;
        }
      }
      if (getVprp) {
        catchmentAvgTmp$vprp = apply(t(t(vprp.xts) * w),1,sum,na.rm=TRUE);
        if (calcVariance) {
          catchmentVarTmp$vprp = apply(vprp.xts,1,spatial.function.name,na.rm=TRUE);
        } else {
          catchmentVarTmp$vprp = NA;
        }
      }
      if (getSolarrad) {
        catchmentAvgTmp$solarrad = apply(t(t(solarrad.xts) * w),1,sum,na.rm=TRUE);
        catchmentAvgTmp$solarrad_interp = apply(t(t(solarrad_interp.xts) * w),1,sum,na.rm=TRUE);
        if (calcVariance) {
          catchmentVarTmp$solarrad = apply(solarrad.xts,1,spatial.function.name,na.rm=TRUE);
          catchmentVarTmp$solarrad_interp = apply(solarrad_interp.xts,1,spatial.function.name,na.rm=TRUE);
        } else {
          catchmentVarTmp$solarrad = NA;
          catchmentVarTmp$solarrad_interp = NA;
        }

      }
      if (getET) {
        catchmentAvgTmp$ET_mm = apply(t(t(ET.xts) * w),1,sum,na.rm=TRUE);
        if (calcVariance) {
          catchmentVarTmp$ET_mm = apply(ET.xts,1,spatial.function.name,na.rm=TRUE);
        } else {
          catchmentVarTmp$ET_mm = NA;
        }
      }

    } else {
      # Do not undertake spatial aggregation. Instead collate
      # results for each grid cell.

      nGridCells = catchment.lookup[i,2] - catchment.lookup[i,1] +1
      catchmentAvgTmp = data.frame(CatchmentID=rep(catchments@data[i,1],nGridCells))
      catchmentVar = NA
      if (getPrecip) {
        newDataCol = as.data.frame(t(precip.xts))
        colnames(newDataCol) = paste('Precip',format.Date(colnames(newDataCol),'%Y_%m_%d'),sep='_')
        catchmentAvgTmp = cbind.data.frame(catchmentAvgTmp, newDataCol)
      }
      if (getTmin) {
        newDataCol = as.data.frame(t(tmin.xts))
        colnames(newDataCol) = paste('Tmin',format.Date(colnames(newDataCol),'%Y_%m_%d'),sep='_')
        catchmentAvgTmp = cbind.data.frame(catchmentAvgTmp, newDataCol)
      }
      if (getTmax) {
        newDataCol = as.data.frame(t(tmax.xts))
        colnames(newDataCol) = paste('Tmax',format.Date(colnames(newDataCol),'%Y_%m_%d'),sep='_')
        catchmentAvgTmp = cbind.data.frame(catchmentAvgTmp, newDataCol)
      }
      if (getVprp) {
        newDataCol = as.data.frame(t(vprp.xts))
        colnames(newDataCol) = paste('vprp',format.Date(colnames(newDataCol),'%Y_%m_%d'),sep='_')
        catchmentAvgTmp = cbind.data.frame(catchmentAvgTmp, newDataCol)
      }
      if (getSolarrad) {
        newDataCol = as.data.frame(t(solarrad.xts))
        colnames(newDataCol) = paste('solarrad',format.Date(colnames(newDataCol),'%Y_%m_%d'),sep='_')
        catchmentAvgTmp = cbind.data.frame(catchmentAvgTmp, newDataCol)

        newDataCol = as.data.frame(t(solarrad_interp.xts))
        colnames(newDataCol) = paste('solarrad_interp',format.Date(colnames(newDataCol),'%Y_%m_%d'),sep='_')
        catchmentAvgTmp = cbind.data.frame(catchmentAvgTmp, newDataCol)
      }
      if (getET) {
        newDataCol = as.data.frame(t(ET.xts))
        colnames(newDataCol) = paste('ET_mm',format.Date(colnames(newDataCol),'%Y_%m_%d'),sep='_')
        catchmentAvgTmp = cbind(catchmentAvgTmp, newDataCol)
      }

    }


    if (i==1) {
      catchmentAvg = catchmentAvgTmp;
      if (do.spatial.analysis) {
        catchmentVar = catchmentVarTmp;
      }
    } else {
      catchmentAvg = rbind(catchmentAvg,catchmentAvgTmp)
      if (do.spatial.analysis) {
        catchmentVar = rbind(catchmentVar,catchmentVarTmp)
      }
    }
  }
  # end for-loop

  if (isCatchmentsPolygon) {
    if (do.spatial.analysis) {
      catchmentAvg = list(catchmentAvg, catchmentVar)
      names(catchmentAvg) = c(paste('catchmentTemporal.',temporal.function.name,sep=''), paste('catchmentSpatial.',spatial.function.name,sep=''))
    } else {
      # Convert data to  a spatial grid (SpatialPixelsDataFrame)
      gridCoords = data.frame(Long=longLat.all[,1], Lat=longLat.all[,2])
      catchmentAvg = cbind.data.frame(gridCoords, catchmentAvg)
      sp::coordinates(catchmentAvg) <- ~Long+Lat
      sp::gridded(catchmentAvg) <- T
    }
  }

  message('Data extraction FINISHED..')
  return(catchmentAvg)
}
