#' Build netCDF files of the Bureau of Meteorology (Australia) national gridded climate data.
#'
#' \code{makeNetCDF_file} builds two netCDF files containing Australian climate data.
#'
#' This function creates two netCDF files of daily climate data. One file contains precipitation, minimum
#' daily temperature, maximum daily temperature and vappour pressure. It should span from 1/1/1900 to today
#' and be ~360GB. The second file contains the solar radiation and started from 1/1/1990 and be ~24GB and
#' spatial gaps are infilled using a 3x3 moving average repeated 3 times. The climate data
#' is sourced from the  Bureau of Meteorology Australian Water Availability Project
#' (\url{http://www.bom.gov.au/jsp/awap/}.  For details see Jones et al. (2009).
#'
#' The output from this function is required for all data extraction functions within this package and must
#' be ran prior.
#'
#' The function can be used to build netCDF files from stratch or to update existng netCDF files previously
#' derived from this function. To not build or update a variable, set its respective URL to \code{NA}.
#'
#' @param ncdfFilename is a file path (as string) and name to the netCDF file. The default is \code{file.path(getwd(),'AWAP.nc')}.
#' @param ncdfSolarFilename is the file path (as string) and name to the netCDF file. The default is \code{file.path(getwd(),'AWAP_solar.nc')}.
#' @param updateFrom is a date string specifying the start date for the AWAP data. If
#'  \code{ncdfFilename} and \code{ncdfSolarFilename} are specified and exist, then the netCDF grids will be
#'  updated with new data from \code{updateFrom}. To update the files from the end of the last day in the file set \code{updateFrom=NA}. The default is \code{"1900-1-1"}.
#' @param updateTo is a date string specifying the end date for the AWAP data. If
#'  \code{ncdfFilename} and \code{ncdfSolarFilename} are specified and exist, then the netCDF grids will be
#'  updated with new data to \code{updateFrom}. The default is today's date as YYYY-MM-DD.
#' @param workingFolder is the file path (as string) in which to download the AWAP grid files. The default is \code{getwd()}.
#' @param keepFiles is a logical scalar to keep the downloaded AWAP grid files. The default is \code{FALSE}.
#' @param urlPrecip URL to the folder containing the AWAP daily precipittaion grids.
#'  The default is \url{http://www.bom.gov.au/web03/ncc/www/awap/rainfall/totals/daily/grid/0.05/history/nat/}.
#' @param urlTmin URL to the folder containing the AWAP daily minimum temperature grids.
#'  The default is \url{http://www.bom.gov.au/web03/ncc/www/awap/temperature/minave/daily/grid/0.05/history/nat/}.
#' @param urlTmax URL to the folder containing the AWAP daily maximum temperature grids.
#'  The default is \url{http://www.bom.gov.au/web03/ncc/www/awap/temperature/maxave/daily/grid/0.05/history/nat/}.
#' @param urlVprp URL to the folder containing the AWAP daily vapour pressure grids.
#'  The default is \url{http://www.bom.gov.au/web03/ncc/www/awap/vprp/vprph15/daily/grid/0.05/history/nat/}.
#' @param urlSolarrad URL to the folder containing the AWAP daily solar radiation grids.
#'  The default is \url{http://www.bom.gov.au/web03/ncc/www/awap/solar/solarave/daily/grid/0.05/history/nat/}.
#'
#' @seealso \code{\link{extractCatchmentData}} for extracting catchment daily average and variance data.
#'
#' @references
#' David A. Jones, William Wang and Robert Fawcett, (2009), High-quality spatial climate data-sets for Australia,
#' Australian Meteorological and Oceanographic Journal, 58 , p233-248.
#'
#' @examples
#' # Load required packages.
#' library(sp);library(raster);library(chron);library(ncdf4);
#' library(maptools);library(Evapotranspiration);library(AWAPer)
#'
#' # Example 1. Build netCDF grids for all existing time points.
#' makeNetCDF_file()
#'
#' # Example 2. Build netCDF grids for ONLY precipitation data at all existing time points.
#' makeNetCDF_file(urlTmin = NA,urlTmax = NA, urlVprp = NA, urlSolarrad = NA)
#'
#' # Example 3. Build netCDF grids for all data but only over a defined time period.
#' # Note if the netCDFs have already been created, then the netCDF files will be updated.
#' makeNetCDF_file(updateFrom='2000-1-1',updateTo='2002-12-31')
#'
#' # Example 4. Update the netCDF grids (from example 3) from the end dates within the netCDF files to the current date.
#' makeNetCDF_file(updateFrom=NA)
#'
#' @export
makeNetCDF_file <- function(
  ncdfFilename=file.path(getwd(),'AWAP.nc'),
  ncdfSolarFilename=file.path(getwd(),'AWAP_solar.nc'),
  updateFrom = as.Date("1900-01-01","%Y-%m-%d"),
  updateTo  = as.Date(Sys.Date(),"%Y-%m-%d"),
  workingFolder=getwd(),
  keepFiles=FALSE,
  urlPrecip = 'http://www.bom.gov.au/web03/ncc/www/awap/rainfall/totals/daily/grid/0.05/history/nat/',
  urlTmin = 'http://www.bom.gov.au/web03/ncc/www/awap/temperature/minave/daily/grid/0.05/history/nat/',
  urlTmax = 'http://www.bom.gov.au/web03/ncc/www/awap/temperature/maxave/daily/grid/0.05/history/nat/',
  urlVprp = 'http://www.bom.gov.au/web03/ncc/www/awap/vprp/vprph15/daily/grid/0.05/history/nat/',
  urlSolarrad = 'http://www.bom.gov.au/web03/ncc/www/awap/solar/solarave/daily/grid/0.05/history/nat/')  {


  # Load appropriate packages that are used
  #library(R.utils)
  #library(ncdf4)
  #library(raster)


  # Check that the ncdf files
  doUpdate = F;
  doUpdate_solar = F;
  if (file.exists(ncdfFilename) && file.exists(ncdfSolarFilename)) {
    doUpdate = TRUE
    doUpdate_solar = TRUE;
    message('Starting to update both netCDF files.')
  } else if (file.exists(ncdfFilename) && !file.exists(ncdfSolarFilename)) {
    message('Starting to update the main climate data netCDF file and will build the solar radiation netcdf file.')
  } else if (!file.exists(ncdfFilename) && file.exists(ncdfSolarFilename)) {
    message('Starting to update the solar radiation netcdf file and will build the main climate data netcdf file.')
    doUpdate_solar = TRUE;
  } else {
    message('Starting to build both netCDF files.')
  }

  # Assess if the main climate data netCDF grid need building/updating
  doMainCDF = TRUE
  if (is.na(urlPrecip) && is.na(urlTmin) && is.na(urlTmax) && is.na(urlVprp)) {
    doMainCDF = FALSE;
    warning('The main climate data netCDF file will not be built or updated.')
  }

  # Test the AWAP downloading
  filedate_str = '20000101'
  haveGridGeometry = FALSE;
  haveGridGeometry_solar = FALSE;
  if (doMainCDF) {
    if (!is.na(urlPrecip)) {
      message('... Testing downloading of AWAP precip. grid')
      url = paste(urlPrecip,filedate_str, filedate_str,'.grid.Z',sep='')
      destFile = file.path(workingFolder,paste('precip.',filedate_str,'.grid.Z',sep=''))
      didFail = download.file(url,destFile, quiet = FALSE)
      if (didFail)
        error('Downloading the grid failed. Please check the URL and the internet connection.')

      # Get the grid geometry of the non solar data
      message('... Getting grid gemoetry from file.')
      system(paste('znew -f ',destFile));
      destFile = file.path(workingFolder,paste('precip.',filedate_str,'.grid.gz',sep=''))
      raw<-textConnection(readLines(a<-gzfile(destFile)));
      headerData = readLines(raw,n=6)
      nCols =as.integer(sub('ncols', '', headerData[1]))
      nRows = as.integer(sub('nrows', '', headerData[2]));
      SWLong = as.numeric(sub('xllcenter', '', headerData[3]));
      SWLat = as.numeric(sub('yllcenter', '', headerData[4]));
      DPixel = as.numeric(sub('cellsize', '', headerData[5]));
      nodata = as.numeric(sub('nodata_value', '', headerData[6]));
      close(raw)
      file.remove(destFile)
      haveGridGeometry = TRUE;
    }

    # Test the AWAP downloading
    if (!is.na(urlTmin)) {
      message('... Testing downloading of AWAP tmin grid')
      url = paste(urlTmin,filedate_str, filedate_str,'.grid.Z',sep='')
      destFile = file.path(workingFolder,paste('tmin.',filedate_str,'.grid.Z',sep=''))
      didFail = download.file(url,destFile, quiet = FALSE)
      if (didFail)
        error('Downloading the grid failed. Please check the URL and the internet connection.')

      # Get grid geometry if not available from precip
      if (is.na(urlPrecip)) {
        # Get the grid geometry of the non solar data
        message('... Getting grid gemoetry from file.')
        system(paste('znew -f ',destFile));
        destFile = file.path(workingFolder,paste('tmin.',filedate_str,'.grid.gz',sep=''))
        raw<-textConnection(readLines(a<-gzfile(destFile)));
        headerData = readLines(raw,n=6)
        nCols =as.integer(sub('ncols', '', headerData[1]))
        nRows = as.integer(sub('nrows', '', headerData[2]));
        SWLong = as.numeric(sub('xllcenter', '', headerData[3]));
        SWLat = as.numeric(sub('yllcenter', '', headerData[4]));
        DPixel = as.numeric(sub('cellsize', '', headerData[5]));
        nodata = as.numeric(sub('nodata_value', '', headerData[6]));
        close(raw)
      }
      file.remove(destFile)
      haveGridGeometry = TRUE;
    }

    # Test the AWAP downloading
    if (!is.na(urlTmax)) {
      message('... Testing downloading of AWAP tmax grid')
      url = paste(urlTmax,filedate_str, filedate_str,'.grid.Z',sep='')
      destFile = file.path(workingFolder,paste('tmax.',filedate_str,'.grid.Z',sep=''))
      didFail = download.file(url,destFile, quiet = FALSE)
      if (didFail)
        error('Downloading the grid failed. Please check the URL and the internet connection.')

      # Get grid geometry if not available from precip
      if (is.na(urlPrecip)) {
        # Get the grid geometry of the non solar data
        message('... Getting grid gemoetry from file.')
        system(paste('znew -f ',destFile));
        destFile = file.path(workingFolder,paste('tmax.',filedate_str,'.grid.gz',sep=''))
        raw<-textConnection(readLines(a<-gzfile(destFile)));
        headerData = readLines(raw,n=6)
        nCols =as.integer(sub('ncols', '', headerData[1]))
        nRows = as.integer(sub('nrows', '', headerData[2]));
        SWLong = as.numeric(sub('xllcenter', '', headerData[3]));
        SWLat = as.numeric(sub('yllcenter', '', headerData[4]));
        DPixel = as.numeric(sub('cellsize', '', headerData[5]));
        nodata = as.numeric(sub('nodata_value', '', headerData[6]));
        close(raw)
        haveGridGeometry = TRUE;
      }
      file.remove(destFile)
    }

    # Test the AWAP downloading
    if (!is.na(urlVprp)) {
      message('... Testing downloading of AWAP vapour pressure grid')
      url = paste(urlVprp,filedate_str, filedate_str,'.grid.Z',sep='')
      destFile = file.path(workingFolder,paste('vprp.',filedate_str,'.grid.Z',sep=''))
      didFail = download.file(url,destFile, quiet = FALSE)
      if (didFail)
        error('Downloading the grid failed. Please check the URL and the internet connection.')

      # Get grid geometry if not available from precip
      if (is.na(urlPrecip)) {
        # Get the grid geometry of the non solar data
        message('... Getting grid gemoetry from file.')
        system(paste('znew -f ',destFile));
        destFile = file.path(workingFolder,paste('vprp.',filedate_str,'.grid.gz',sep=''))
        raw<-textConnection(readLines(a<-gzfile(destFile)));
        headerData = readLines(raw,n=6)
        nCols =as.integer(sub('ncols', '', headerData[1]))
        nRows = as.integer(sub('nrows', '', headerData[2]));
        SWLong = as.numeric(sub('xllcenter', '', headerData[3]));
        SWLat = as.numeric(sub('yllcenter', '', headerData[4]));
        DPixel = as.numeric(sub('cellsize', '', headerData[5]));
        nodata = as.numeric(sub('nodata_value', '', headerData[6]));
        close(raw)
        haveGridGeometry = TRUE;
      }
      file.remove(destFile)

    }
  }


  # Test the AWAP solar downloading
  if (!is.na(urlSolarrad)) {
    message('... Testing downloading of AWAP solar grid')
    url = paste(urlSolarrad,filedate_str, filedate_str,'.grid.Z',sep='')
    destFile = file.path(workingFolder,paste('solarrad.',filedate_str,'.grid.Z',sep=''))
    didFail = download.file(url,destFile, quiet = FALSE)
    if (didFail)
      error('Downloading the grid failed. Please check the URL and the internet connection.')

    # Get the grid geometry of the non solar data
    message('... Getting grid gemoetry from file')
    system(paste('znew -f ',destFile));
    destFile = file.path(workingFolder,paste('solarrad.',filedate_str,'.grid.gz',sep=''))
    raw<-textConnection(readLines(a<-gzfile(destFile)));
    headerData = readLines(raw,n=6)
    nCols_solar =as.integer(sub('ncols', '', headerData[1]))
    nRows_solar = as.integer(sub('nrows', '', headerData[2]));
    SWLong_solar = as.numeric(sub('xllcenter', '', headerData[3]));
    SWLat_solar = as.numeric(sub('yllcenter', '', headerData[4]));
    DPixel_solar = as.numeric(sub('cellsize', '', headerData[5]));
    nodata_solar = as.numeric(sub('nodata_value', '', headerData[6]));
    close(raw)
    file.remove(destFile)
    haveGridGeometry_solar = TRUE;
  }

  # Create net CDF files
  if (doMainCDF) {

    if (!doUpdate) {

      message('... Building AWAP netcdf file.')

      # Convert and check input timedates
      if (is.na(updateFrom) || nchar(updateFrom)==0) {
        updateFrom = as.Date("1900-01-01","%Y-%m-%d")
      } else if (is.character(updateFrom))
        updateFrom = as.Date(updateFrom,'%Y-%m-%d');
      if (is.character(updateTo))
        updateTo = as.Date(updateTo,'%Y-%m-%d');
      if (updateFrom >= updateTo)
        stop('The update dates are invalid. updateFrom must be prior to updateTo')

      # Set data time points
      timepoints = seq( as.Date("1900-01-01","%Y-%m-%d"), by="day", to=as.Date(Sys.Date(),"%Y-%m-%d"))

      # Get x and y vectors (dimensions)
      Longvector = seq(SWLong, by=DPixel,length.out = nCols)
      Latvector = seq(SWLat, by=DPixel,length.out = nRows)

      # define dimensions
      londim <- ncdim_def("Long","degrees",vals=Longvector)
      latdim <- ncdim_def("Lat","degrees",vals=Latvector)
      timedim <- ncdim_def("time","days since 1900-01-01 00:00:00.0 -0:00",unlim=T, vals=0:(length(timepoints)-1), calendar='standard')

      # define variables
      fillvalue <- NA
      dlname <- "min daily temperature"
      tmin_def <- ncvar_def("tmin","deg_C",list(londim,latdim,timedim),fillvalue,dlname,prec="single")
      dlname <- "min daily temperature"
      tmax_def <- ncvar_def("tmax","deg_C",list(londim,latdim,timedim),fillvalue,dlname,prec="single")
      dlname <- "vapour pressure"
      vprp_def <- ncvar_def("vprp","hPa",list(londim,latdim,timedim),fillvalue,dlname,prec="single")
      dlname <- "precipitation"
      precip_def <- ncvar_def("precip","mm",list(londim,latdim,timedim),fillvalue,dlname,prec="single")

      # create netCDF file and put arrays
      ncout <- nc_create(ncdfFilename,list(tmin_def,tmax_def,vprp_def,precip_def),force_v4=T)

      # put additional attributes into dimension and data variables
      ncatt_put(ncout,"Long","axis","X")
      ncatt_put(ncout,"Lat","axis","Y")
      ncatt_put(ncout,"time","axis","T")

      # add global attributes
      ncatt_put(ncout,0,"title","BoM AWAP daily climate data")
      ncatt_put(ncout,0,"institution","Data: BoM, R Code: Uni. of Melbourne")
      history <- paste("R Code: T.J. Peterson", date(), sep=", ")

      # Get netcdf time points
      timePoints_netCDF <- ncvar_get(ncout, "time")

      # Convert netCDF time points to date
      tunits <- ncatt_get(ncout, "time", "units")
      tustr <- strsplit(tunits$value, " ")
      tdstr <- strsplit(unlist(tustr)[3], "-")
      tmonth = as.integer(unlist(tdstr)[2])
      tday = as.integer(unlist(tdstr)[3])
      tyear = as.integer(unlist(tdstr)[1])
      timePoints_R = chron(timePoints_netCDF, origin = c(tmonth, tday, tyear));

    } else {
      # open netcdf file
      ncout <- nc_open(ncdfFilename, write=T)

      # Get netcdf time points
      timePoints_netCDF <- ncvar_get(ncout, "time")

      # Convert netCDF time points to date
      tunits <- ncatt_get(ncout, "time", "units")
      tustr <- strsplit(tunits$value, " ")
      tdstr <- strsplit(unlist(tustr)[3], "-")
      tmonth = as.integer(unlist(tdstr)[2])
      tday = as.integer(unlist(tdstr)[3])
      tyear = as.integer(unlist(tdstr)[1])
      timePoints_R = chron(timePoints_netCDF, origin = c(tmonth, tday, tyear));

      # Set updateFrom to the end of the netCDF file if updateFrom is NA or ''.
      if (is.character(updateFrom)) {
        if (nchar(updateFrom)>0) {
          updateFrom = min(c(max(as.Date(timePoints_R)),as.Date(updateFrom,'%Y-%m-%d')));
        } else {
          updateFrom = max(as.Date(timePoints_R));
        }
      } else {
        if (is.na(updateFrom)) {
          updateFrom = max(as.Date(timePoints_R));
        }
      }

      # Set updateTo to the min of the input data and now.
      if (is.character(updateTo)) {
        if (nchar(updateTo)>0) {
          updateTo = min(c(as.Date(Sys.Date(),"%Y-%m-%d"),as.Date(updateTo,'%Y-%m-%d')));
        } else {
          updateTo = as.Date(Sys.Date(),"%Y-%m-%d")
        }
      } else {
        if (is.na(updateTo)) {
          updateTo = as.Date(Sys.Date(),"%Y-%m-%d")
        }
      }
    }

    # Check input dates
    if (updateFrom >= updateTo)
      stop('The update dates are invalid. updateFrom must be prior to updateTo')

    timepoints2Update = seq( as.Date(updateFrom,'%Y-%m-%d'), by="day", to=as.Date(updateTo,'%Y-%m-%d'))

    if (length(timepoints2Update)==0)
        stop('The dates to update produce a zero vector of dates of zero length. Check the inputs dates are as YYYY-MM-DD')

    message(paste('... NetCDF data will be  extracted from ',format.Date(updateFrom,'%Y-%m-%d'),' to ', format.Date(updateTo,'%Y-%m-%d')));

    # Start Filling the netCDF grid.
    message('... Starting to add data AWAP netcdf file.')
    for (date in 1:length(timepoints2Update)){

      # Get datestring for input filenames
      datestring<-format(timepoints2Update[date], "%Y%m%d")
      message(paste('Working on grid time point:', format(timepoints2Update[date], "%Y-%m-%d")))

      # Find index to the date to update within the net CDF grid
      ind = as.integer(difftime(timepoints2Update[date], as.Date("1900-1-1",'%Y-%m-%d'),units = "days" ))+1

      # Update timePoints_netCDF time vector
      timePoints_netCDF[ind] = ind-1;

      # Download data
      #----------------
      didFail_precip=1
      if (!is.na(urlPrecip)) {
        url = paste(urlPrecip,datestring, datestring,'.grid.Z',sep='')
        destFile_precip = file.path(workingFolder,paste('precip.',datestring,'.grid.Z',sep=''))
        didFail_precip = download.file(url,destFile_precip, quiet = T)
        if (didFail_precip==0) {
          system(paste('znew -f ',destFile_precip));
          destFile_precip = gsub('.Z', '.gz', destFile_precip)
        }
      }

      didFail_tmin=1
      if (!is.na(urlTmin)) {
        url = paste(urlTmin,datestring, datestring,'.grid.Z',sep='')
        destFile_tmin = file.path(workingFolder,paste('tmin.',datestring,'.grid.Z',sep=''))
        didFail_tmin = download.file(url,destFile_tmin, quiet = T)
        if (didFail_tmin==0) {
          system(paste('znew -f ',destFile_tmin));
          destFile_tmin = gsub('.Z', '.gz', destFile_tmin)
        }
      }

      didFail_tmax=1
      if (!is.na(urlTmax)) {
        url = paste(urlTmax,datestring, datestring,'.grid.Z',sep='')
        destFile_tmax = file.path(workingFolder,paste('tmax.',datestring,'.grid.Z',sep=''))
        didFail_tmax = download.file(url,destFile_tmax, quiet = T)
        if (didFail_tmax==0) {
          system(paste('znew -f ',destFile_tmax));
          destFile_tmax = gsub('.Z', '.gz', destFile_tmax)
        }
      }

      didFail_vprp=1
      if (!is.na(urlVprp)) {
        url = paste(urlVprp,datestring, datestring,'.grid.Z',sep='')
        destFile_vprp = file.path(workingFolder,paste('vprp.',datestring,'.grid.Z',sep=''))
        didFail_vprp = download.file(url,destFile_vprp, quiet = T)
        if (didFail_vprp==0) {
          system(paste('znew -f ',destFile_vprp));
          destFile_vprp = gsub('.Z', '.gz', destFile_vprp)
        }
      }
      #----------------

      # Get precip grid and add to Net CDF grid
      if (!is.na(urlPrecip) && file.exists(destFile_precip) && didFail_precip==0) {
        raw<-textConnection(readLines(a<-gzfile(destFile_precip)));
        AWAPgrid<- as.matrix(t(read.table(raw, skip=6, nrow=nRows)))
        close(raw);
        AWAPgrid <- AWAPgrid[,ncol(AWAPgrid):1]
        ncvar_put( ncout, "precip", AWAPgrid, start=c(1, 1, ind), count=c(nCols, nRows, 1), verbose=F )
      }
      if (!is.na(urlPrecip) && file.exists(destFile_precip) && !keepFiles)
        file.remove(destFile_precip)

      # Get tmin grid and add to Net CDF grid
      if (!is.na(urlTmin) && file.exists(destFile_tmin) && didFail_tmin==0) {
        raw<-textConnection(readLines(a<-gzfile(destFile_tmin)));
        AWAPgrid<- as.matrix(t(read.table(raw, skip=6, nrow=nRows)))
        close(raw);
        AWAPgrid <- AWAPgrid[,ncol(AWAPgrid):1]
        ncvar_put( ncout, "tmin", AWAPgrid, start=c(1, 1, ind), count=c(nCols, nRows, 1), verbose=F )
      }
      if (!is.na(urlTmin) && file.exists(destFile_tmin) && !keepFiles)
        file.remove(destFile_tmin)

      # Get tmax grid and add to Net CDF grid
      if (!is.na(urlTmax) && file.exists(destFile_tmax) && didFail_tmax==0) {
        raw<-textConnection(readLines(a<-gzfile(destFile_tmax)));
        AWAPgrid<- as.matrix(t(read.table(raw, skip=6, nrow=nRows)))
        close(raw);
        AWAPgrid <- AWAPgrid[,ncol(AWAPgrid):1]
        ncvar_put( ncout, "tmax", AWAPgrid, start=c(1, 1, ind), count=c(nCols, nRows, 1), verbose=F )
      }
      if (!is.na(urlTmax) && file.exists(destFile_tmax) && !keepFiles)
        file.remove(destFile_tmax)

      # Get vapour pr grid and add to Net CDF grid
      if (!is.na(urlVprp) && file.exists(destFile_vprp) && didFail_vprp==0) {
        raw<-textConnection(readLines(a<-gzfile(destFile_vprp)));
        AWAPgrid<- as.matrix(t(read.table(raw, skip=6, nrow=nRows)))
        close(raw);
        AWAPgrid <- AWAPgrid[,ncol(AWAPgrid):1]
        ncvar_put( ncout, "vprp", AWAPgrid, start=c(1, 1, ind), count=c(nCols, nRows, 1), verbose=F )
      }
      if (!is.na(urlVprp) && file.exists(destFile_vprp) && !keepFiles)
        file.remove(destFile_vprp)


    }

    # Updtate netCDF time variable
    ncvar_put(ncout, "time",timePoints_netCDF)

    # close the file, writing data to disk
    nc_close(ncout)
  }


  # Example extractions
  #tmp_array <- ncvar_get(ncout,"tmax",start=c(1,1,date),count=c(nCols,nRowsPrec,1))
  #long <- ncvar_get(ncout,"Long")
  #lat <- ncvar_get(ncout,"Lat")
  #image(long,lat,tmp_array,  col=rev(brewer.pal(10,"RdBu")))

  # BUILD SOLAR DATA NETCDF
  #----------------------------------------------

  if (!haveGridGeometry_solar)
    warning('Thesolar radiation data netCDF file will not be built or updated.')

  # Create net CDF files
  if (haveGridGeometry_solar) {

    # Set data time points
    timepoints = seq( as.Date("1990-01-01","%Y-%m-%d"), by="day", to=as.Date(Sys.Date(),"%Y-%m-%d"))

    if (!doUpdate_solar) {

      Longvector = seq(SWLong_solar, by=DPixel_solar,length.out = nCols_solar)
      Latvector = seq(SWLat_solar, by=DPixel_solar,length.out = nRows_solar)

      # define dimensions
      londim <- ncdim_def("Long","degrees",vals=Longvector)
      latdim <- ncdim_def("Lat","degrees",vals=Latvector)
      timedim <- ncdim_def("time","days since 1990-01-01 00:00:00.0 -0:00",unlim=T, vals=0:(length(timepoints)-1), calendar='standard')

      # define variables
      fillvalue <- NA
      dlname <- "Solar radiation"
      solar_def <- ncvar_def("solarrad","MJ/m^2",list(londim,latdim,timedim),fillvalue,dlname,prec="single")

      # create netCDF file and put arrays
      ncout <- nc_create(ncdfSolarFilename,list(solar_def),force_v4=T)

      # put additional attributes into dimension and data variables
      ncatt_put(ncout,"Long","axis","X")
      ncatt_put(ncout,"Lat","axis","Y")
      ncatt_put(ncout,"time","axis","T")

      # put additional attributes into dimension and data variables
      ncatt_put(ncout,"Long","axis","X")
      ncatt_put(ncout,"Lat","axis","Y")
      ncatt_put(ncout,"time","axis","T")

      # add global attributes
      ncatt_put(ncout,0,"title","BoM AWAP daily solar data")
      ncatt_put(ncout,0,"institution","Data: BoM, R Code: Uni of. Melbourne")
      history <- paste("R Code: T.J. Peterson", date(), sep=", ")

      # Get netcdf time points
      timePoints_netCDF <- ncvar_get(ncout, "time")

      # Convert netCDF time points to date
      tunits <- ncatt_get(ncout, "time", "units")
      tustr <- strsplit(tunits$value, " ")
      tdstr <- strsplit(unlist(tustr)[3], "-")
      tmonth = as.integer(unlist(tdstr)[2])
      tday = as.integer(unlist(tdstr)[3])
      tyear = as.integer(unlist(tdstr)[1])
      timePoints_R = chron(timePoints_netCDF, origin = c(tmonth, tday, tyear));

      timepoints2Update = seq( as.Date(updateFrom,'%Y-%m-%d'), by="day", to=as.Date(updateTo,'%Y-%m-%d'))

    } else {

      # open netcdf file
      ncout <- nc_open(ncdfSolarFilename, write=T)

      # Get netcdf time points
      timePoints_netCDF <- ncvar_get(ncout, "time")

      # Convert netCDF time points to date
      tunits <- ncatt_get(ncout, "time", "units")
      tustr <- strsplit(tunits$value, " ")
      tdstr <- strsplit(unlist(tustr)[3], "-")
      tmonth = as.integer(unlist(tdstr)[2])
      tday = as.integer(unlist(tdstr)[3])
      tyear = as.integer(unlist(tdstr)[1])
      timePoints_R = chron(timePoints_netCDF, origin = c(tmonth, tday, tyear));

      # Set updateFrom to the min of the end of the netCDF file and updateFrom
      if (is.character(updateFrom)) {
        if (nchar(updateFrom)>0) {
          updateFrom = min(c(max(as.Date(timePoints_R)), as.Date(updateFrom,'%Y-%m-%d')));
        } else {
          updateFrom = max(as.Date(timePoints_R));
        }
      } else {
        if (is.na(updateFrom)) {
          updateFrom = max(as.Date(timePoints_R));
        }
      }

      # Set updateTo to the end of the netCDF file if updateTo is NA or ''.
      if (is.character(updateTo)) {
        if (nchar(updateTo)>0) {
          updateTo = min(c(as.Date(Sys.Date(),"%Y-%m-%d"), as.Date(updateTo,'%Y-%m-%d')));
        } else {
          updateTo = as.Date(Sys.Date(),"%Y-%m-%d")
        }
      } else {
        if (is.na(updateTo)) {
          updateTo = as.Date(Sys.Date(),"%Y-%m-%d")
        }
      }

    }

    # Check input dates
    if (updateFrom >= updateTo)
      stop('The update dates are invalid. updateFrom must be prior to updateTo')

    timepoints2Update = seq( as.Date(updateFrom,'%Y-%m-%d'), by="day", to=as.Date(updateTo,'%Y-%m-%d'))

    if (length(timepoints2Update)==0)
      stop('The dates to update produce a zero vector of dates of zero length. Check the inputs dates are as YYYY-MM-DD')

    message(paste('... NetCDF Solar data will be  extracted from ',format.Date(updateFrom,'%Y-%m-%d'),' to ', format.Date(updateTo,'%Y-%m-%d')));

    # Start Filling the netCDF grid.
    message('... Starting to add data AWAP Solar netcdf file.')

    # Start Filling the netCDF grid.
    for (date in 1:length(timepoints2Update)){

      # Get datestring for input filenames
      datestring<-format(timepoints2Update[date], "%Y%m%d")
      message(paste('Working on solar grid time point:', format(timepoints2Update[date], "%Y-%m-%d")))

      # Find index to the date to update within the net CDF grid
      ind = as.integer(difftime(timepoints2Update[date], as.Date("1990-1-1",'%Y-%m-%d'),units = "days" ))+1

      # Update timePoints_netCDF time vector
      timePoints_netCDF[ind] = ind-1;

      # Download the file
      didFail=1
      if (!is.na(urlSolarrad)) {
        url = paste(urlSolarrad,datestring, datestring,'.grid.Z',sep='')
        destFile_solarrad = file.path(workingFolder,paste('solarrad',datestring,'.grid.Z',sep=''))
        didFail = download.file(url,destFile_solarrad, quiet = T)
        if (didFail==0) {
          system(paste('znew -f ',destFile_solarrad));
          destFile_solarrad = gsub('.Z', '.gz', destFile_solarrad)
        }
      }

      # Get vapour pr grid and add to Net CDF grid
      if (file.exists(destFile_solarrad) && didFail==0) {
        # Import file
        raw<-textConnection(readLines(a<-gzfile(destFile_solarrad)));
        AWAPgrid<- raster(as.matrix(t(read.table(raw, skip=6, nrow=nRows_solar))))
        close(raw);

        # Infill NA values of grid by taking the local average and convert back to matrix.
        AWAPgrid = focal(AWAPgrid, w=matrix(1,3,3), fun=mean, na.rm=TRUE, NAonly=TRUE)
        AWAPgrid = focal(AWAPgrid, w=matrix(1,3,3), fun=mean, na.rm=TRUE, NAonly=TRUE)
        AWAPgrid = focal(AWAPgrid, w=matrix(1,3,3), fun=mean, na.rm=TRUE, NAonly=TRUE)
        #AWAPgrid = as.matrix(AWAPgrid);

        # Reorder grid for ncdf
        AWAPgrid <- AWAPgrid[,ncol(AWAPgrid):1]

        # Add to ncdf
        ncvar_put( ncout, "solarrad", AWAPgrid, start=c(1, 1, ind), count=c(nCols_solar, nRows_solar, 1), verbose=F )
      }

      if (file.exists(destFile_solarrad)  && !keepFiles)
        file.remove(destFile_solarrad)
    }

    # Updtate netCDF time variable
    ncvar_put(ncout, "time",timePoints_netCDF)

    # close the file, writing data to disk
    nc_close(ncout)
  }

  #nc_open(ncfname)
  #tmp_array <- ncvar_get(ncout,"solarrad",start=c(1,1,date),count=c(nCols,nRowsPrec,1))
  #long <- ncvar_get(ncout,"Long")
  #lat <- ncvar_get(ncout,"Lat")
  #image(long,lat,tmp_array,  col=rev(brewer.pal(10,"RdBu")))
}
