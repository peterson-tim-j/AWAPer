#' Downloads and imports Geoscience Australia 9s DEM.
#'
#' \code{getDEM} get Australian 9s DEM.
#'
#' getDEM downloads the Geoscience Australia 9 second DEM and then imports the grid.
#'
#' @details
#' The DEM is required for the calculation of evaportranspiration within \code{extractCatchmentData}. For details of the DEM see
#' \url{https://www.data.gov.au/dataset/geodata-9-second-dem-and-d8-digital-elevation-model-version-3-and-flow-direction-grid-2008}
#'
#' @param workingFolder is the file path (as string) in which to download the zip file. The default is \code{getwd()}.
#' @param urlDEM URL to the folder containing the Geoscience Australia 9s DEM.
#'  The default is taken from \code{getURLs()$DEM}.
#' @param DEMfilename is the file name for the DEM (as string). The default is \code{'dem-9s.asc'}.
#' @param keepFiles is a logical scalar to keep the downloaded zip file and extracted DEM ASCII file. The default is \code{FALSE}.
#'
#' @return
#' A RasterLayer DEM for Asutralia.
#'
#' @seealso \code{\link{extractCatchmentData}} for extracting catchment daily average and variance data.
#'
#' @examples
#' # Download the DEM.
#' \donttest{
#' # Set file name for DEM file to the system temp. directory.
#' DEM.file = tempfile(fileext='.asc')
#'
#' # Download and import the Australian 9 second DEM.
#' DEM_9s = getDEM(workingFolder=dirname(DEM.file),
#'          DEMfilename=basename(DEM.file))
#' }
#' @export
getDEM <- function(workingFolder=getwd(),
                   urlDEM = getURLs()$DEM,
                   DEMfilename = 'dem-9s.asc',
                   keepFiles=F) {

  # Download .zip
  message('... Downloading DEM .zip file')
  destFile = file.path(workingFolder,paste('DEM-9S_ESRI.zip',sep=''))
  didFail = utils::download.file(urlDEM,destFile, quiet = FALSE)

  if (didFail) {
    message('Downloading the DEM failed. Please check the URL and your internet connection.')
    stop()
  }

  # get list of files within zip
  zip.names=utils::unzip(destFile, list=T,junkpaths=T)
  ind = which(regexpr(DEMfilename,zip.names$Name) != -1)

  if (length(DEMfilename)==0)
    stop(paste('The following DEM file name could not be found in the zip file:',DEMfilename))

  DEMfilename.2unzip = zip.names$Name[ind]

  # Unzip file
  message('... Extracting DEM from the zip file')
  utils::unzip(destFile, files=DEMfilename.2unzip,junkpaths=T)

  message('... Reading in DEM file')
  DEM <- sp::read.asciigrid(DEMfilename.2unzip,colname='DEM',  proj4string = sp::CRS("+proj=longlat +ellps=GRS80"))

  message('... Converting grid to a raster data type.')
  DEM = raster::raster(DEM)

  if (!keepFiles) {

    if (file.exists(destFile))
      file.remove(destFile)

    if (file.exists(DEMfilename.2unzip))
      file.remove(DEMfilename.2unzip)

  }

  return(DEM)

}

