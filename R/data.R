#' @title Example catchment boundary polygons.
#' @description Two example catchment boundaries as a SpatialPolygonsDataFrame. The catchments are Creswick Creek (ID 407214, Vic., Australia, see http://www.bom.gov.au/water/hrs/#id=407214) and
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
#' data("catchments")
#'
#' # Load the 9 second Australian DEM.
#' data("DEM_9s")
#'
#' # Plot the catchment boundaries.
#' image(DEM_9s, xlab='Long.',ylab='Lat.')
#' plot(catchments,add=T)
#'
#' @export
"catchments"
