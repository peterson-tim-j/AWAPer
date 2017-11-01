#' Downloads and imports Geoscience Australia 9s DEM.
#'
#' \code{getDEM} get Australian 9s DEM.
#'
#' This function downlaod the Geoscience Australia 9 second DEM and then imports the grid.
#' The DEM is required for the calculation of evaportranspiration within \code{extractCatchmentData}. For details of the DEM see
#' \url{https://www.data.gov.au/dataset/geodata-9-second-dem-and-d8-digital-elevation-model-version-3-and-flow-direction-grid-2008}
#'
#' @param workingFolder is the file path (as string) in which to download the zip file. The default is \code{getwd()}.
#' @param keepFiles is a logical scalar to keep the downloaded zip file and extracted DEM ASCII file. The default is \code{FALSE}.
#' @param urlDEM URL to the folder containing the Geoscience Australia 9s DEM.
#'  The default is \url{https://s3-ap-southeast-2.amazonaws.com/elvis.ga.gov.au/elevation/9secPackagedData/DEM-9S_ESRI.zip}.
#'
#' @return
#' A RasterLayer DEM for Asutralia.
#'
#' @seealso \code{\link{extractCatchmentData}} for extracting catchment daily average and variance data.
#'
#' @examples
#' # Load required packages.
#' library(sp);library(raster);library(chron);library(ncdf4);
#' library(maptools);library(Evapotranspiration);library(AWAPer)
#'
#' # Download the DEM.
#' DEM_9s = getDEM()
#'
#' # Plot the DEM.
#' image(DEM_9s, xlab='Long.',ylab='Lat.')
#'
#' # Save DEM for next time it is needed.
#' save(DEM_9s,file="DEM.RData" )
#'
#' @export
getDEM <- function(workingFolder=getwd(),
                   urlDEM = 'https://s3-ap-southeast-2.amazonaws.com/elvis.ga.gov.au/elevation/9secPackagedData/DEM-9S_ESRI.zip',
                   keepFiles=F) {

  # Download .zip
  message('... Downloading DEM .zip file')
  destFile = file.path(workingFolder,paste('DEM-9S_ESRI.zip',sep=''))
  didFail = download.file(urlDEM,destFile, quiet = FALSE)

  if (didFail) {
    message('Downloading the DEM failed. Please check the URL and your internet connection.')
    stop()
  }

  # Unzip file
  message('... Extracting dem-9s.asc from the zip file')
  unzip(destFile, files="Data_9secDEM_D8/dem-9s.asc",junkpaths=T)

  message('... Reading in dem-9s.asc')
  DEM <- read.asciigrid('dem-9s.asc',colname='DEM',  proj4string = CRS("+proj=longlat +ellps=GRS80"))

  message('... Converting grid to a raster data type.')
  DEM = raster(DEM)

  if (!keepFiles) {

    if (file.exists(destFile))
      file.remove(destFile)

    if (file.exists('dem-9s.asc'))
      file.remove('dem-9s.asc')

  }

  return(DEM)

}

