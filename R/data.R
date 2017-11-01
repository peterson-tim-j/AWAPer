#' Australian 9 second DEM
#'
#' \code{data("DEM_9s")} load Australian 9s DEM.
#'
#' Loads the Geoscience Australia 9 second DEM. For details see
#' \url{https://www.data.gov.au/dataset/geodata-9-second-dem-and-d8-digital-elevation-model-version-3-and-flow-direction-grid-2008}
#'
#' @seealso
#' \code{\link{extractCatchmentData}} for extracting catchment average climate data.
#'
#' @examples
#' # Load required packages.
#' library(sp);library(raster);library(chron);library(AWAPer);library(ncdf4);
#' library(maptools);library(Evapotranspiration)
#'
#' # Load the 9 second Australian DEM.
#' data("DEM_9s")
#'
#' # Plot the DEM.
#' image(DEM_9s, xlab='Long.',ylab='Lat.')
#'
"DEM_9s"

#' Example catchment boundary polygon.
#'
#' \code{data("catchmentBoundaries")} load example catchment boundary polygons.
#'
#' Loads two example catchment boundaries as a SpatialPolygonsDataFrame. The catchments are Creswick Creek (ID 407214, Vic., Australia, see http://www.bom.gov.au/water/hrs/#id=407214) and
#' Bet Bet Creek (ID 407220, Vic., Australia, see http://www.bom.gov.au/water/hrs/#id=407220). The catchments can be used to extract catchment average climate data usng \code{extractCatchmentData}
#'
#' @seealso
#' \code{\link{extractCatchmentData}} for extracting catchment average climate data.
#'
#' @examples
#' # Load required packages.
#' library(sp);library(raster);library(chron);library(AWAPer);library(ncdf4);
#' library(maptools);library(Evapotranspiration)
#'
#' # Load example cacthment boundaries.
#' data("catchmentBoundaries")
#'
#' # Load the 9 second Australian DEM.
#' data("DEM_9s")
#'
#' # Plot the catchment boundaries.
#' image(DEM_9s, xlab='Long.',ylab='Lat.')
#' plot(catchments,add=T)
#'
"catchmentBoundaries"
